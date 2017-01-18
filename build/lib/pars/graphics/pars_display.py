import lg_base.graphics.graphics_toolkit as gt
import numpy as np
import lg_base.core.array_utils as hau
try:
    from bottleneck import nanmean,nansum
except:
    print
    print 'No bottleneck.nanmean available! Falling back to SLOW scipy.stats.nanmean'
    print
    from scipy.stats import nanmean
    from numpy import nansum

def show_pars(instrument,display_defaults,rs,usetimes,usealts,figs):
            #toplevel
            if usealts==None and hasattr(rs,'heights'):
                usealts=rs.heights
            if usetimes==None and hasattr(rs,'times'):
                usetimes=rs.times
            #delta_t=np.zeros(rs.times.shape)
            #make vector of time differences between data points
            delta_t = hau.T_Array(np.array([(rs.times[x+1]-rs.times[x]).total_seconds() for x in range(rs.times.size-1)]+[0.0]))
            if delta_t.size>1:
                delta_t[-1]=nanmean(delta_t[:-2])
    
            if figs==None:
                figs=gt.figurelist()
          
                
            #start here
           
            if instrument.find('pars2S1') >= 0:
                color ='r'
                instrment = 'Pars2S1'
            else:
                color = 'b'
                instrument = 'Pars2S2'

            if display_defaults.enabled('parsivel_precip_rate'):                
                if hasattr(rs,'precip_rate'):
                    gt.plot_vs_time('parsivel_precip_rate'
                       ,instrument
                       ,usetimes
                       ,[rs.precip_rate]  
                       ,[color]
                       ,[2]
                       ,[]
                       ,''    
                       ,'Parsivel precip rate (mm/hr)'
                       ,[]
                       ,'Parsivel precip rate vs time'   
                       ,False #FIXME
                       ,display_defaults
                       ,figs)
                    print 'rendering Parsivel precip rate'               
                else:
                    print
                    print 'No Parsivel precip rate image--precip rate field not found'
                    
            if display_defaults.enabled('parsivel_accumulated_precip'):
                if hasattr(rs,'precip_rate'):
                    #precip rate in units of mm/hr--convert seconds to hours
                    temp = (rs.precip_rate * delta_t/3600.0).copy()
                    temp[np.isnan(temp)]=0.0
                    accumulated_precip = np.cumsum(temp,0)   
                    gt.plot_vs_time('parsivel_accumulated_precip'
                       ,instrument
                       ,usetimes
                       ,[accumulated_precip]  
                       ,[color]
                       ,[2]
                       ,[]
                       ,''    
                       ,'Parsivel precip (mm)'
                       ,[]
                       ,'Parsivel precip vs time'   
                       ,False #FIXME
                       ,display_defaults
                       ,figs)
                    print 'rendering Parsivel precip'               
                else:
                    print
                    print 'No Parsivel precip image--precip_rate field not found'

            if display_defaults.enabled('parsivel_median_volume_diameter'):
                if hasattr(rs,'median_volume_diameter'):
                    temp = rs.median_volume_diameter.copy()
                    temp[temp <0] = np.NaN
                    gt.plot_vs_time('parsivel_median_volume_diameter'
                       ,instrument
                       ,usetimes
                       ,[temp]  
                       ,[color]
                       ,[2]
                       ,None
                       ,''    
                       ,'median volume diameter'
                       ,'mm'
                       ,'parsivel median vol. dia. vs time'   
                       ,False #FIXME
                       ,display_defaults
                       ,figs)
                    print 'rendering Parsivel median volume diameter'               
                else:
                    print
                    print 'No Parsivel median volume diameter plot--field not found'
                    
            if display_defaults.enabled('parsivel_liquid_water_content'):
                if hasattr(rs,'liquid_water_content'):
            
                    gt.plot_vs_time('parsivel_liquid_water_content'
                       ,instrument
                       ,usetimes
                       ,[rs.liquid_water_content * 1e-3]  #convert from mm^3/m^3 to gr/m^3  
                       ,[color]
                       ,[2]
                       ,None
                       ,''    
                       ,'Liquid water content'
                       ,'mm'
                       ,'Parsivel LWC vs time'   
                       ,False #FIXME
                       ,display_defaults
                       ,figs)
                    print 'rendering Parsivel liquid water content'               
                else:
                    print
                    print 'No Parsivel liquid water content plot--field not found'
                    
            if display_defaults.enabled('parsivel_radar_reflectivity'):
                if hasattr(rs,'equivalent_radar_reflectivity'):
                    temp = rs.equivalent_radar_reflectivity.copy()
                    temp[temp<-50] = np.NaN
                    gt.plot_vs_time('parsivel_radar_reflectiviy'
                       ,instrument
                       ,usetimes
                       ,[temp]   
                       ,[color]
                       ,[2]
                       ,None
                       ,''    
                       ,'Radar reflectivity'
                       ,'dBZ'
                       ,'Parsivel radar reflectivity vs time'   
                       ,False #FIXME
                       ,display_defaults
                       ,figs)
                    print 'rendering Parsivel radar reflectivity'               
                else:
                    print
                    print 'No Parsivel equivalent radar reflectivity--field not found'
            if 0:
              print 'particle size'
              print rs.particle_size.shape
              print rs.particle_size
              print 'particle width'
              print rs.class_size_width.shape
              print rs.class_size_width
              print 'fall velocity'
              print rs.raw_fall_velocity.shape
              print rs.raw_fall_velocity
            
            if 0: #display_defaults.enabled('parsivel_fall_velocity'):
                if hasattr(rs,'raw_fall_velocity'):
            
                    gt.plot_vs_time('parsivel_fall_velocity'
                       ,instrument
                       ,usetimes
                       ,[rs.raw_fall_velocity]   
                       ,[color]
                       ,[2]
                       ,None
                       ,''    
                       ,'Raw fall velocity'
                       ,'m/s'
                       ,'Parsivel fall velocity vs time'   
                       ,False #FIXME
                       ,display_defaults
                       ,figs)
                    print 'rendering Parsivel raw fall velocity'               
                else:
                    print
                    print 'No Parsivel raw fall velocity plot--field not found'
                    
            if 0: # display_defaults.enabled('parsivelfall_velocity_vs_time'):
                if hasattr(rs,'raw_spectrum'):

                    fall_velocity = rs.raw_spectrum.copy()
                    
                    fall_velocity[fall_velocity < -9000.0] = 0.0
                    #sum over particle size
                    fall_velocity = np.sum
                
                    gt.plot_vs_time('parsivel_fall_velocity_vs_time'
                        ,instrument
                        ,usetimes
                        ,[fall_velocity]  
                        ,[color]
                        ,[2]
                        ,[]
                        ,''    
                        ,'Parsivel fall velocity (m/s)'
                        ,[]
                        ,'Parsivel fall velocity vs time'   
                        ,False #FIXME
                        ,display_defaults
                        ,figs)
                    print 'rendering Parsivel raw fall velocity'               
                else:
                     print
                     print 'No Parsivel raw fall velocity plot--field not found'








            if display_defaults.enabled('parsivel_size_spectrum'):
                if hasattr(rs,'raw_spectrum'):
                    raw_spectrum = rs.raw_spectrum.copy()
                    raw_spectrum[raw_spectrum < -9000] = 0.0
                    #sum over time
                    raw_spectrum = np.sum(raw_spectrum,0)
                    #sum over velocities
                    raw_spectrum = np.sum(raw_spectrum,1)
                   
                    np.set_printoptions(threshold=np.NaN)

                   

                    gt.plot_xy('parsivel_size_spectrum'  #plot name
                        ,'Parsivel'
                        ,rs.times
                        ,[rs.particle_size[:19]]
                        ,[raw_spectrum[:19]]
                        ,[color]
                        ,[]
                        ,[]
                        ,['-']
                        ,[2]
                        ,[]
                        ,'upper right'
                        ,'particle size '
                        ,'mm'
                        ,'number'
                        ,None
                        ,'Raw spectrum'
                        ,'' #text_str
                        ,None #'text_position_x
                        ,None #text_position_y
                        ,None #text_angle
                        ,display_defaults
                        ,figs)
                else:
                    print
                    print 'Parsivel raw spectum plot not plotted--variable not found'

            if display_defaults.enabled('parsivel_fall_velocity_spectrum'):
                if hasattr(rs,'raw_spectrum'):
                    raw_spectrum = rs.raw_spectrum.copy()
                    raw_spectrum[raw_spectrum < -9000] = 0.0
                    #sum over time
                    raw_spectrum = np.sum(raw_spectrum,0)
                    #sum over sizes
                    raw_spectrum = np.sum(raw_spectrum,1)
                    
                    np.set_printoptions(threshold=np.NaN)

                    gt.plot_xy('parsivel_fall_velocity_spectrum'  #plot name
                        ,'Parsivel'
                        ,rs.times
                        ,[rs.raw_fall_velocity]
                        ,[raw_spectrum]
                        ,[color]
                        ,[]
                        ,[]
                        ,['-']
                        ,[2]
                        ,[]
                        ,'upper right'
                        ,'fall velocity'
                        ,'m/s'
                        ,'number'
                        ,None
                        ,'fall spectrum'
                        ,'' #text_str
                        ,None #'text_position_x
                        ,None #text_position_y
                        ,None #text_angle
                        ,display_defaults
                        ,figs)
                else:
                    print
                    print 'Parsivel fall velocity spectum  not plotted--variable not found'

            if display_defaults.enabled('parsivel_number_detected_particles_vs_time'):
                if hasattr(rs,'number_detected_particles'):
                    print 'rendering parsivel number_detected_particles_vs_time'
                    gt.plot_vs_time('parsivel_number_detected_particles_vs_time'
                        ,instrument
                        ,usetimes
                        ,[rs.number_detected_particles]  
                        ,[color]
                        ,[2]
                        ,[]
                        ,''    
                        ,'Number of particles'
                        ,None
                        ,'Parsivel number of detected particles vs time'   
                        ,False #FIXME
                        ,display_defaults
                        ,figs)
                            
                else:
                   print
                   print 'Pareivel number_detected_particles not plotted--variable not found'
                   
            if display_defaults.enabled('parsivel_number_density'):
                if hasattr(rs,'number_density_drops'):
                    number_density = rs.number_density_drops.copy()
                    number_density[number_density < -9000] = 0.0
                    #sum over time
                    number_density = np.mean(number_density,0)
                        
                    #np.set_printoptions(threshold=np.NaN)
            
                    gt.plot_xy('parsivel_number_density'
                        ,instrument
                        ,rs.times
                        ,[rs.particle_size[:19]]
                        ,[number_density[:19]]
                        ,[color]
                        ,[]
                        ,[]
                        ,['-']
                        ,[2]
                        ,[]
                        ,'upper right'
                        ,'particle size '
                        ,'mm'
                        ,'number density'
                        ,'1/(m^3 *mm)'
                        ,'Number density'
                        ,'' #text_str
                        ,None #'text_position_x
                        ,None #text_position_y
                        ,None #text_angle
                        ,display_defaults
                        ,figs)
                else:
                    print
                    print 'Parsivel number density not plotted--variable not found'
