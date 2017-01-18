import dplkit.role.librarian
import dplkit.role.narrator
import dplkit.role.filter
from datetime import datetime,timedelta
import os
import re
import time
import numpy
import lg_dpl_toolbox.filters.substruct as substruct
import lg_base.core.array_utils as hau
import dplkit.role.decorator
import numpy as np
import copy
from bottleneck import nanmax

attrrep=dict(nyquist_velocity=None,
             gate_spacing=None,
             range_gate_spacing='gate_spacing',
             pulse_repetition_frequency=None,
             radar_operating_frequency='operating_frequency')


def addCommonAttributes(x,attrs):
    for attr,dest in attrrep.items():
        if attr not in attrs:
            continue
        v=attrs[attr][:].split()
        val=hau.T_Array(np.ones(x.times.shape)*float(v[0]))
        if len(v)>1:
            if v[1] in ('(GHz)','GHz'):
                val*=1e9
        setattr(x,dest or attr,val)

def make_radar_constant(lambda_radar):

        
        #dielectric constants squared --these are the values that were used
        #to compute the reflectivity, ARM assumes 0-deg C values
        if lambda_radar < 5.0e-3: #for 3mm radar
             k_water_sq = 0.7 # from ARM SWACR handbook
        elif lambda_radar <  1.0e-2: #for 8.6 mm radar     
             k_water_sq = 0.88 #from ARM KAZR handbook    
        else: #x band
             k_water_sq = 0.93 #ARM value for X-band
             
        k_sq_ice = 0.176
       
        #k_water_sq=0.92;  
        #k_sq_ice=0.197;

        #radar_constant = 2 * numpy.pi**5 * k_water_sq / (3 * lambda_radar**4)
        radar_constant = np.pi**4 * k_water_sq /( 4.0 * lambda_radar**4)

        radar_constant=radar_constant*1e-18; #convert from mm^6/m^3 to m^6/m^3
        return radar_constant

class RadarPrefilter(dplkit.role.filter.aFilter):
    """ RadarPrefilter for MMCR, KAZR, and WACR
        converts reflectivity to backscatter, and rescales content appropriately
    """
    def __init__(self,source):
        super(self.__class__,self).__init__(source)
        self.framestream=source
        if source.provides!=None:
            self.provides=source.provides.copy()
        if self.provides!=None and 'Reflectivity' in self.provides:
            self.provides['Backscatter']=self.provides['Reflectivity'].copy()
            self.provides['Backscatter']['shortname']='Backscatter'
            self.provides['Backscatter']['units']="1/(m sr)"
        #also copied elsewhere
        lambda_radar=self.framestream.radarLambda# 8.6e-3;

        self.radar_constant=make_radar_constant(lambda_radar)



    def reflectivityToBackscatterCrossSection(self, ref): 
        """ ret=reflectivityToBackscatterCrossSection(ref)
convert dBZ*100 to radar backscatter cross section.
note that by definition the 'radar cross section' is the backscatter
cross section times 4pi.
radar constant from Donovan--J. Geophys. Res. vol 106 p 27426 Nov
2001
see definition on p 284 of van de Hulst and equations on page 432
and page 70 van de Hulst
this code generates radar_backscat in units of 1/(m str) which
conforms to lidar usage conventions.

dbz in terms of reflectivity Z
dBz=10*log10(Z mm^6/m^3) 
converting units on reflectivity from  mm^6/m^3 to m^6/m^3
dBz=10*log10(Z*10^18 m^6/m^3)
Z=1e-18*10^(dBz/10)
note that radar merge files contain dBz values mulitplied by 100
        """
        #Radar calibration kludge for magic
        magic_start_time = datetime(2013,4,1)
        magic_end_time   = datetime(2013,10,30)
        if self.framestream.starttime > magic_start_time \
               and self.framestream.starttime < magic_end_time:
            #kludge&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
            #Reflectivity calibration adjustment for magic mwacr or kazer data
        
            if self.framestream.radarLambda < 5e-3 :
                reflectivity_adj = 0.0 #-3.5 # based on 16 jul 13, 8-9UTC  #-5.4 aug 5
                print
                print 'WACR reflectivity adjusted by ',reflectivity_adj, ' dBz in RadarFilters.py'
                print
                ref = ref + reflectivity_adj
            elif self.framestream.radarLambda < 1e-2:   #kzar wavelength
                reflectivity_adj = 0.0 # +1.0 aug 5
                print
                print 'KAZR reflectivity adjusted by ',reflectivity_adj, ' dBz in RadarFilters.py'
                print
                ref = ref + reflectivity_adj
            else:
                raise RuntimeError, 'No adjustment radar adjustment factor for this wavelength'                 
            
        ret=self.radar_constant*pow(10.0,(ref/10)); #backscat cross section
        return ret

    def process(self):
        isMmcr=(self.framestream.radarType=='MMCR') #FIXME better way to discover the division
        isKazr=(self.framestream.radarType.startswith('KAZR')) #FIXME better way to discover the division
        isMwacr=('WACR' in self.framestream.radarType) #FIXME better way to discover the division
        for f in self.framestream:
            f=copy.copy(f)
            if hasattr(f,'heights'):
            #    f._altitudevarname='heights'
                f.heights=hau.Z_Array(f.heights,dtype=numpy.float64)
            if hasattr(f,'Reflectivity'):
                f.Reflectivity=hau.TZ_Array(f.Reflectivity,dtype=numpy.float64)
                if isMmcr:
                    f.Reflectivity/=100.0
                setattr(f,'Backscatter',hau.TZ_Array(self.reflectivityToBackscatterCrossSection(f.Reflectivity)))
            if hasattr(f,'SpectralWidth'):
                f.SpectralWidth=hau.TZ_Array(f.SpectralWidth,dtype=numpy.float64)
                if isMmcr:
                    f.SpectralWidth/=1000.0
            if hasattr(f,'MeanDopplerVelocity'):
                f.MeanDopplerVelocity=hau.TZ_Array(f.MeanDopplerVelocity,dtype=numpy.float64)
                if isMmcr:
                    f.MeanDopplerVelocity/=1000.0
                if isKazr or isMwacr:
                    f.MeanDopplerVelocity*=-1.0
            yield f

class RadarBackscatterToReflectivity(dplkit.role.filter.aFilter):
    """ Radar Backatter to Reflectivity
        converts backscatter back to reflectivity, which is preferred after interpolation since reflectivity is a log-scale value
    """
    def __init__(self,source):
        super(self.__class__,self).__init__(source)
        self.framestream=source
        #self.provides=source.provides
  
        #also copied elsewhere
        lambda_radar=self.framestream.radarLambda#8.6e-3;

        self.radar_constant=make_radar_constant(lambda_radar)

    def backscatterCrossSectionToReflectivity(self, ref): 
        """ ret=reflectivityToBackscatterCrossSection(ref)
convert dBZ*100 to radar backscatter cross section.
note that by definition the 'radar cross section' is the backscatter
cross section times 4pi.
radar constant from Donovan--J. Geophys. Res. vol 106 p 27426 Nov
2001
see definition on p 284 of van de Hulst and equations on page 432
and page 70 van de Hulst
this code generates radar_backscat in units of 1/(m str) which
conforms to lidar usage conventions.

dbz in terms of reflectivity Z
dBz=10*log10(Z mm^6/m^3) 
converting units on reflectivity from  mm^6/m^3 to m^6/m^3
dBz=10*log10(Z*10^18 m^6/m^3)
Z=1e-18*10^(dBz/10)
note that radar merge files contain dBz values mulitplied by 100
        """
        ret=numpy.log10(ref/self.radar_constant)*10.0; #backscat cross section
        return ret

    def process(self):
        for f in self.framestream:
            if hasattr(f,'Backscatter'):
                f=copy.copy(f)
                setattr(f,'Reflectivity',hau.TZ_Array(self.backscatterCrossSectionToReflectivity(f.Backscatter)))
            yield f
