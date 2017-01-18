from collections import namedtuple
import lg_base.core.array_utils as hau
import numpy as np
from datetime import datetime
from lg_base.core.locate_file import locate_file
import json
import multiple_scattering_utilities as msu
try:
    from bottleneck import nanmean,nansum,anynan
except ImportError:
    print
    print 'No bottleneck.nanmean available! Falling back to SLOW scipy.stats.nanmean'
    print
    from scipy.stats import nanmean
    from numpy import nansum
    def anynan(x):
        return np.any(np.isnan(x))  

class multiple_scattering(object):
    def __init__(self,parameters,particle_parameters,*args,**kwargs):#update calling parameters to anything needed to be kept here
        #dpl object will expose the parameters object. if the name 'parameters' changes, dpl objects 'multiple_scattetering_parameters' function must be updated to return it
        #  this is so the artist will have access to these values
        
        self.load_mscatter_parameters(parameters,particle_parameters)
        self.parameters = parameters
        self.ms_obj = None #used as accumulator for computation of time averaged ms profiles
        


    def load_mscatter_parameters(self,parms,particle_parameters):
       
        systemOnly=(parms==None)
      
        if parms==None:
            parms='multiple_scatter_parameters_default.json'
        if isinstance(parms,basestring):
            from lg_base.core.locate_file import locate_file
            parms=json.load(open(locate_file(parms,systemOnly=systemOnly),'r'))
        if not isinstance(parms,dict):
            raise RuntimeError('multiple_scatter_parameters need to be a json filename or a dictionary')
        self.multiple_scatter_parameters=parms

        # if using lidar_radar derived particle info get particle size info from "particle_parameters.json"
        print
        if not particle_parameters == None \
                  and self.multiple_scatter_parameters['particle_info_source'] == 'measured' :
            self.multiple_scatter_parameters['p180_water'] \
                        = particle_parameters['p180_water']['value']
            self.multiple_scatter_parameters['p180_ice'] \
                        = particle_parameters['p180_ice'] 
            self.multiple_scatter_parameters['h2o_depol_threshold'] \
                        = particle_parameters['h2o_depol_threshold']
            self.multiple_scatter_parameters['alpha_water'] \
                        = particle_parameters['alpha_water']
            self.multiple_scatter_parameters['alpha_ice'] \
                        = particle_parameters['alpha_ice']
            self.multiple_scatter_parameters['g_water'] \
                        = particle_parameters['g_water']
            self.multiple_scatter_parameters['g_ice'] \
                        = particle_parameters['g_ice']
            print 'multiple scatttering particle mode diameter from lidar-radar measured profile'
            print 'particle parameters for multiple scattering taken from "particle_parameters.json"'
        elif particle_parameters == None \
               and not self.multiple_scatter_parameters['particle_info_source'] == 'constant' :
            raise RuntimeError, '"particle_parameters.json" not present and multiple_scatter["particle_info_source"] !="constant"'
            if self.multiple_scatter_parameters['processing_mode'] == '1d':
                print '    processing single 1-d profile as mean over time interval requested'
            else:
                print '    processing 2-d (time,altitude) multiple scatter correction'
        else:
            print 'multiple scattering parmeters taken from "multiple_scatter_parameters.json"'
            if self.multiple_scatter_parameters['processing_mode'] == '1d':
                print '    processing single 1-d mean profile over requested time interval'
            else:
                print '    processing 2-d (time,altitude) multiple scattering'
            print '    mode_diameter_water  = ',self.multiple_scatter_parameters['mode_diameter_water'], '(m)'
            print '    mode_diameter_ice    = ',self.multiple_scatter_parameters['mode_diameter_ice'], '(m)' 
        print  '    starting altitude   = ', self.multiple_scatter_parameters['lowest_altitude'], '(m)'
        print  '    p180_water          = ', self.multiple_scatter_parameters['p180_water'] 
        print  '    p180_ice            = ', self.multiple_scatter_parameters['p180_ice']
        print  '    h2o_depol_threshold = ', self.multiple_scatter_parameters['h2o_depol_threshold']*100.0, '(%)'
        print  '    alpha_water         = ', self.multiple_scatter_parameters['alpha_water']
        print  '    alpha_ice           = ', self.multiple_scatter_parameters['alpha_ice']
        print  '    g_water             = ', self.multiple_scatter_parameters['g_water']
        print  '    g_ice               = ', self.multiple_scatter_parameters['g_ice']
        print
            
    def __call__(self,rs_inv,rs_particle,calvals):#update to process,and return a completed frame
        rs_multiple_scattering=hau.Time_Z_Group(rs_inv.times.copy(),timevarname='times',altname='msl_altitudes')

        N1 = 2  
        N2 = self.multiple_scatter_parameters['highest_order']
        
        start_alt = self.multiple_scatter_parameters['lowest_altitude']
        
        p180_water = self.multiple_scatter_parameters['p180_water']
        p180_ice   = self.multiple_scatter_parameters['p180_ice']
        h2o_depol_threshold = self.multiple_scatter_parameters['h2o_depol_threshold']
       

        p180_ice = self.multiple_scatter_parameters['p180_ice']
        p180_water = self.multiple_scatter_parameters['p180_water']
        second_wavelength =self.multiple_scatter_parameters['second_wavelength']

        wavelength = calvals['wavelength']*1e-9

        #assert(rs_particle!=None or self.multiple_scatter_parameters['particle_info_source'] == 'constant')

       
        
        if 1:   #self.multiple_scatter_parameters['processing_mode'] == '1d':
            if self.multiple_scatter_parameters['particle_info_source'] == 'constant':
                #in this case the mode diameter will be reset to ice or water values from multiple_scattering.json file
                #depending on h2o_depol_threshold
                mode_diameter = None
            else:
                #use lidar-radar retrieved particle sizes
                #print 'particle'
                #print dir(rs_particle)
                mode_diameter = rs_particle.mode_diameter.copy()

            #accumulates sums for 1-d average multiple scatter profile
            self.ms_obj = msu.ms_sums(mode_diameter,rs_inv.beta_a_backscat
                         ,rs_inv.linear_depol,self.multiple_scatter_parameters,self.ms_obj)

            self.ms_obj.beta_water[self.ms_obj.beta_water < 0] = 0.0
            self.ms_obj.beta_ice[self.ms_obj.beta_ice < 0] = 0.0
            
            beta_total = self.ms_obj.beta_water +self.ms_obj.beta_ice
            
            total_samples = self.ms_obj.n_samples_water +self.ms_obj.n_samples_ice
            
            #when no ice or water data points are present averages must be zero
            #compute weighted averages of beta, diameter when ice and water are present
            #ms_obj.beta_water and ms_obj.beta_ice have the sum of beta_backscatter for ice and water
            #self.ms_obj.n_samples_water[self.ms_obj.n_samples_water == 0] = np.infty
            #self.ms_obj.n_samples_ice[self.ms_obj.n_samples_ice == 0] = np.infty
            
            
            #compute ave beta_extinction profile from sum beta_backscat profiles
            extinction_profile = (self.ms_obj.beta_water/p180_water + self.ms_obj.beta_ice/p180_ice)/total_samples
            diameter = (self.ms_obj.diameter_ice + self.ms_obj.diameter_water)/beta_total
            
          
            #convert altitudes into range
            ranges = rs_inv.msl_altitudes.copy()
            zenith_angle = np.abs(calvals['telescope_roll_angle_offset'])*np.pi/180.0
            ranges = ranges/np.cos(zenith_angle)
            start_range = start_alt/np.cos(zenith_angle)
            end_range = rs_inv.msl_altitudes[-1]/np.cos(zenith_angle)
            if start_range >= end_range:
               raise RuntimeError(' start altitude'+str(np.int(start_range))+ ' is above highest data point')
            
            ms_ratios_profile = msu.msinteg(N1,N2,start_range \
                     ,end_range,self.multiple_scatter_parameters['step_size'], extinction_profile, diameter 
                     ,ranges,wavelength,self.multiple_scatter_parameters,calvals) 

            rs_multiple_scattering.ranges = ms_ratios_profile[:,0]
            rs_multiple_scattering.ms_ratios_profile = ms_ratios_profile
            rs_multiple_scattering.extinction_profile = hau.Z_Array(extinction_profile[np.newaxis,:])   
            rs_multiple_scattering.weighted_diameter = hau.Z_Array(diameter[np.newaxis,:])
            rs_multiple_scattering.msl_altitudes = rs_inv.msl_altitudes.copy()
            rs_multiple_scattering.wavelength = wavelength

            #compute multiple scattering for a second wavelength if it is provided
            #assume no change is extinction cross section
            if second_wavelength:
                 ms_ratios_profile_2 = msu.msinteg(N1,N2,start_range \
                     ,end_range,self.multiple_scatter_parameters['step_size'], extinction_profile, diameter 
                     ,ranges,second_wavelength,self.multiple_scatter_parameters,calvals) 

                 rs_multiple_scattering.ms_ratios_profile_2 = ms_ratios_profile_2
                 rs_multiple_scattering.second_wavelength = second_wavelength
           
        if self.multiple_scatter_parameters['processing_mode'] == '2d': #do multiple scatter calculation for all profiles in frame
            print 'begining 2d multiple scatter processing'
            #estimate extinction based on backscatter phase function 
            beta = rs_inv.beta_a_backscat.copy()
            beta = beta/p180_water
            beta[rs_inv.linear_depol>self.multiple_scatter_parameters['h2o_depol_threshold']] \
                 = beta[rs_inv.linear_depol>self.multiple_scatter_parameters['h2o_depol_threshold']]*p180_water/p180_ice
            beta[beta < 0]=0.0
            if self.multiple_scatter_parameters['particle_info_source'] == 'constant':
                mode_diameter = np.ones_like(rs_inv.beta_a_backscat) \
                                * self.multiple_scatter_parameters['mode_diameter_water']
                mode_diameter[rs_inv.linear_depol > self.multiple_scatter_parameters['h2o_depol_threshold']] \
                           = self.multiple_scatter_parameters['mode_diameter_ice']

        
            #convert altitudes into ranges
            ranges = rs_inv.msl_altitudes.copy()
            zenith_angle = np.abs(calvals['telescope_roll_angle_offset'])*np.pi/180.0
            ranges = ranges/np.cos(zenith_angle)
            start_range = start_alt/np.cos(zenith_angle)
            end_range = rs_inv.msl_altitudes[-1]/np.cos(zenith_angle)
            
            for i in range(rs_inv.beta_a_backscat.shape[0]):
                print 'Computing multiple scattering for ' ,rs_inv.times[i]
                if self.multiple_scatter_parameters['particle_info_source'] == 'constant':
                    ratios = msu.msinteg(N1,N2,start_range \
                         ,end_range,self.multiple_scatter_parameters['step_size'],beta[i,:],mode_diameter[i,:] 
                         ,ranges,wavelength,self.multiple_scatter_parameters,calvals)
                else: #get mode diameter from lidar_radar measured values
                   ratios = msu.msinteg(N1,N2,start_range \
                         ,end_range,wavelength,self.multiple_scatter_parameters['step_size']
                         ,beta[i,:],rs_particle.mode_diameter[i,:] 
                         ,ranges,self.multiple_scatter_parameters,calvals)
                   
                #load values into output array 
                if not hasattr(rs_multiple_scattering,'ms_ratio_total'):
                    rs_multiple_scattering.ms_ratio_total = np.zeros((beta.shape[0],ratios.shape[0]))
                
                rs_multiple_scattering.ms_ratio_total[i,:] =nansum(ratios[:,2:],1)
               
                rs_multiple_scattering.msl_altitudes = rs_inv.msl_altitudes.copy()
            rs_multiple_scattering.ms_ratio_total =hau.TZ_Array(rs_multiple_scattering.ms_ratio_total)
                
       
        
        #rs_multiple_scattering.start_altitude = start_alt
       
        return rs_multiple_scattering

   
  
