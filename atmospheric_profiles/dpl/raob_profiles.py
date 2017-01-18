import dplkit.role.librarian
import dplkit.role.narrator
import dplkit.role.filter
from datetime import datetime,timedelta
import numpy
import dplkit.role.decorator
import lg_base.core.array_utils as hau
import os,calendar

from scipy.optimize import brentq
from netCDF4 import Dataset
import numpy as np
import copy

def cal_frost_point(dew_point):
    """Hyland-Wexler calculation of frost point temp (K) given
       the dew point temperature (K).
       Brents method is used to find the point where the vapor
       pressure over water is equal to the vapor pressure over
       ice"""

    nalts = dew_point.shape[0]
    frost_point = np.NaN * np.zeros_like(dew_point)
    log_es = -0.58002206e4 / dew_point + 1.3914993 - 0.48640239e-1 \
        * dew_point + 0.41764768e-4 * dew_point ** 2 - 0.14452093e-7 \
        * dew_point ** 3 + 0.65459673e1 * np.log(dew_point)

    for j in range(nalts):
        if dew_point[j] < 273.15 and np.isfinite(dew_point[ j]):
            try:
                frost_point[ j] = brentq(delta_es, 100, 273, args=log_es[j])
            except ValueError: 
                frost_point[ j] = np.NaN
        else:
            frost_point[ j] = np.NaN
    return frost_point


def delta_es(x, log_es):
    g = log_es - (-0.56745359e4 / x + 6.39225247 - 0.96778430e-2 * x
                  + 0.62215701e-6 * x ** 2 + 0.2074825e-8 * x ** 3
                  - 0.94840240e-12 * x ** 4 + 0.41635019e1 * np.log(x))
    return g


def getAll(nc,varname):
    iarr=[]
    var=nc.variables[varname]
    for d in var.dimensions:
        iarr.append(slice(None,None))#slice(nc_dims[d][0],nc_dims[d][1]))
    #print 'reading variable',alias,'start',starta,'count',counta
    #get requested variable from netcdf
    try:
        if len(iarr)==0:
            tmp=var.getValue()
        elif len(iarr)==1:
            tmp=var[iarr[0]]
        else:
            tmp=var[tuple(iarr)]
    except IndexError:
        print 'Index error reading variable ',varname#,' in file ',filename
        tmp=None
    return tmp

@dplkit.role.decorator.autoprovide(frameclass=hau.Time_Z_Group,reuseGenerator=False)
class dpl_raob_narr(dplkit.role.narrator.aNarrator):
    def __init__(self,hoststream,interval_start_time,interval_end_time=None,expire_duration=None):
        self.hoststream=hoststream
        self.interval_start_time=interval_start_time
        self.interval_end_time=interval_end_time
        self.expire_duration=expire_duration

    def badSounding(self,sounding,comparesounding=None):
        if comparesounding is not None:
            if sounding.times==comparesounding.times:
                return True
        if sum(sounding.pressures>1.0)<2:
            return True
        if not (sounding.top>6000.0):
            return True
        #print vars(sounding)
        return False

    def allsoundings(self,nc,rec,fromidx=0):

        rs =  hau.Time_Z_Group()
        
        times = getAll(nc,'synTime')

        # combine mandatory and sig height measurements
       
        heights = np.hstack((getAll(nc,'htMan'), getAll(nc,'htSigT')))

        epoch = datetime(1970,1,1,0,0,0)
        t_mask = times<1e36
        for i in range(len(times)):
           t_mask[i] = times[i]<1e36 and any(heights[i,:] < 1e36)

        times=times[t_mask]
         
        times = [ epoch + timedelta(seconds=soff) for soff in times[:] ]
        
     # select times, one prior to start time --> last profile in file

        #indices = np.arange(len(times))
        #start_index = max(indices[times <= start_time])
        #rs.times = zeros(len(times) - start_index)

        wmosta = getAll(nc,'wmoStat')  # wmo station number
        stalong = getAll(nc,'staLon')
        stalat = getAll(nc,'staLat')
      
     # combine mandatory and sig height measurements
        
        temps = np.hstack((getAll(nc,'tpMan'), getAll(nc,'tpSigT')))
        pressures = np.hstack((getAll(nc,'prMan'), getAll(nc,'prSigT')))
        dew_points = np.hstack((getAll(nc,'tdMan'), getAll(nc,'tdSigT')))
        wind_dir = np.hstack((getAll(nc,'wdMan'), getAll(nc,'wdSigT')))
        wind_spd = np.hstack((getAll(nc,'wsMan'), getAll(nc,'wsSigT')))
        heights = heights[t_mask,:]        
        temps = temps[t_mask,:]
        pressures = pressures[t_mask,:]
        dew_points = dew_points[t_mask,:]
        wind_dir = wind_dir[t_mask,:]
        wind_spd = wind_spd[t_mask,:]
        frost_points=dew_points.copy()

        pressures[heights > 1e5] = np.NaN
        temps[heights > 1e5] = np.NaN
        dew_points[heights > 1e5] = np.NaN

        dew_points = temps - dew_points

        for i in range(heights.shape[0]):
            indices = np.argsort(heights[i,:])
            heights[i,:] = heights[i,indices]
            temps[i,:] = temps[i,indices]
            pressures[i,:] = pressures[i,indices]
            dew_points[i,:] = dew_points[i,indices]
            wind_dir[i,:] = wind_dir[i,indices]
            wind_spd[i,:] = wind_spd[i,indices]
            frost_points[i,:] = cal_frost_point(dew_points[i,:])

        for tidx in range(fromidx,len(times)):
            keepidx=tidx#np.arange(tidx,tidx+1)
            temp_sounding=hau.Time_Z_Group()
            temp_sounding.type = rec['type']
            temp_sounding.sounding_type = rec['type']
            temp_sounding.times= times[tidx]
            temp_sounding.sample_time= times[tidx]
            temp_sounding.latitude=stalat[keepidx]
            temp_sounding.longitude=stalong[keepidx]
            temp_sounding.stalat=stalat[keepidx]
            temp_sounding.stalong=stalong[keepidx]
            temp_sounding.altitudes=hau.Z_Array(heights[tidx,:])
            amask=temp_sounding.altitudes<90000
            if not np.any(amask):
                continue
            temp_sounding.altitudes=temp_sounding.altitudes[amask]
            temp_sounding.temps= hau.Z_Array(temps[keepidx,:])[amask]
            temp_sounding.pressures= hau.Z_Array(pressures[keepidx,:])[amask]
            temp_sounding.dew_points = hau.Z_Array(dew_points[keepidx,:])[amask]
            temp_sounding.frost_points = hau.Z_Array(frost_points[keepidx,:])[amask]
            temp_sounding.wind_dir = hau.Z_Array(wind_dir[keepidx,:])[amask]
            temp_sounding.wind_spd = hau.Z_Array(wind_spd[keepidx,:])[amask]
            temp_sounding.station_id=rec['id']
            temp_sounding.sounding_id=rec['id']
            temp_sounding.top = max(temp_sounding.altitudes)
            temp_sounding.bot = min(temp_sounding.altitudes)
            yield temp_sounding



    def read(self):
        last_sounding=None
        priorpath=None
        priorlen=0
        for rec in self.hoststream:
            nc=Dataset(rec['path'],'r')
            for si,sounding in enumerate(self.allsoundings(nc,rec,0 if priorpath!=rec['path'] else priorlen)):
                priorlen=si+1
                if self.badSounding(sounding,last_sounding):
                    continue
                if last_sounding is not None and sounding.times>self.interval_start_time:
                    print 'yielding sounding for ',last_sounding.times
                    yield last_sounding
                    print 'NEXT sounding'
                last_sounding=sounding
            priorpath=rec['path']
        if last_sounding is None:
            raise RuntimeError('ran out of radiosondes')
        yield last_sounding

@dplkit.role.decorator.autoprovide(frameclass=dict,reuseGenerator=False)
class dpl_raobfile_narr(dplkit.role.narrator.aNarrator):
    def __init__(self,instrument,sounding_id,interval_start_time,interval_end_time=None,expire_duration=None):
        import lg_dpl_toolbox.core.archival as arc
        self.station=sounding_id
        self.type='NOAA raob'
        self.basedir=arc.get_path_to_data(instrument)
        self.filename='sondes.'+sounding_id+'.nc'
        self.start=datetime(interval_start_time.year,interval_start_time.month,1,0,0,0)
        if interval_end_time is not None:
            self.end=datetime(interval_end_time.year,interval_end_time.month,1,0,0,0)
        else:
            self.end=None
        self.expire_duration=expire_duration
     
    def read(self):
        starttime=copy.copy(self.start)
        while True:
            dur=timedelta(days=calendar.monthrange(starttime.year,starttime.month)[1])
            ret={'path':os.path.join(self.basedir,'%04i' % starttime.year,'%02i' % starttime.month, self.filename),
                    'filename':self.filename,'start':starttime,'width':dur,'type':self.type,'id':self.station}
            if os.path.exists(ret['path']):
                yield ret
            starttime=starttime+dur
            if self.end is not None and starttime>self.end:
                break
            if self.end is None:
                if starttime>datetime.utcnow():
                    starttime=starttime-dur


class dpl_raob(dplkit.role.librarian.aLibrarian):
    def __init__(self,instrument,sounding_id,requested_altitudes=None):
        self.instrument=instrument
        self.sounding_id=sounding_id
        self.requested_altitudes=requested_altitudes

    def search(self,*args,**kwargs):#none may just mean end is determined by caller
        ret = dpl_raob_narr(dpl_raobfile_narr(self.instrument,self.sounding_id,*args,**kwargs),*args,**kwargs)
        if self.requested_altitudes is not None:
            import atmospheric_profiles.dpl.dpl_temperature_profiles as dtp
            ret=dtp.dpl_radiosonderesample(ret,'altitudes',self.requested_altitudes,{'temps':[50,500],'pressures':[1,1000],
                'dew_points':[150,500],'frost_points':[50,500],'altitudes':None})
        return ret


def main():
    lib=dpl_raob('bagohsrl','GRB')
    for x in lib(datetime(2015,9,3,0,0,0),datetime(2015,9,15,0,0,0)):
        print vars(x).keys()

if __name__ == '__main__':
    main()
