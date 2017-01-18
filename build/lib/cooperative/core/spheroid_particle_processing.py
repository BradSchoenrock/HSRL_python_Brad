
from collections import namedtuple
import effective_diameter_prime as edp
import lg_base.core.array_utils as hau
import spheroid_utilities as su
import numpy as np
from datetime import datetime
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

def process_spheroid_particle(rs_inv,rs_radar,particle_parameters,lambda_radar,entire_frame,
                              sounding=None,size_dist=None):
            """
            process_spheroid_particle(rs_inv,rs_radar,particle_parameters,lambda_radar,entire_frame,
                              sounding=None,p180_water=None,size_dist=None):
            generate and return the particle measurements based on a given hsrl inverted data,
            radar (and its lambda), and particle parameters dictionary
            """

            #create timez group and add heights            
            rs_particle=hau.Time_Z_Group(rs_inv.times.copy(),timevarname='times',altname='heights')
            setattr(rs_particle,'heights',rs_inv.msl_altitudes.copy())
            setattr(rs_particle,'delta_t',rs_inv.delta_t.copy())

            #remove points where lidar signal is noise dominated by setting to
            #very small value.
            #clipped_beta_a_back=rs_inv.beta_a_backscat.copy()
            #if 0: #if hasattr(rs_inv,'std_beta_a_backscat'):
            #    clipped_beta_a_back[clipped_beta_a_back<(2*rs_inv.std_beta_a_backscat)]=-np.inf
            #else:
            #    print 'No std_beta_a_backscat statistics to filter particle measurements'
            #clipped_beta_a_back[np.logical_not(np.isfinite(rs_inv.beta_a_backscat))]=-np.inf;
            
            rs_particle.q_backscatter = np.NaN * np.zeros_like(rs_inv.beta_a_backscat)
            rs_particle.phase = np.zeros_like(rs_inv.beta_a_backscat)
            rs_particle.phase[rs_inv.linear_depol > particle_parameters['h2o_depol_threshold']] = 1
            rs_particle.phase[np.isnan(rs_inv.beta_a_backscat)] = np.NaN

            #set aspect ratio parameter for ice filled bins
            rs_particle.zeta = np.ones(rs_inv.beta_a_backscat.shape)
            rs_particle.zeta[rs_inv.linear_depol > particle_parameters['h2o_depol_threshold']] \
                          = particle_parameters['zeta']                       

            print 'Extinction cross section for particle size calculations derived from ' \
                        ,particle_parameters['ext_source']
            print 'Extinction due nonprecipitating aerosols = '\
                        ,particle_parameters['background_aerosol_bs'],'1/(m sr)'

            
            
            #store the mask field with the particle info
            rs_particle.qc_mask = rs_inv.qc_mask.copy()

            clipped_beta_a_backscat = rs_inv.beta_a_backscat.copy()
            copy_radar_backscatter = rs_radar.Backscatter.copy()
            #clipped_beta_a_backscat = copy_beta_a.copy()
            clipped_beta_a_backscat = clipped_beta_a_backscat \
                      - particle_parameters['background_aerosol_bs']
            clipped_beta_a_backscat[clipped_beta_a_backscat < 0] = np.NaN

            #create an empty mode_diameter array     
            rs_particle.mode_diameter = np.zeros_like(rs_inv.beta_a_backscat)

            #create an empty array for extinction--used only for particle calculations
            #bs_ratio_to_dmode will return extinction cross section in clipped_beta_a
            clipped_beta_a = np.NaN * np.zeros_like(rs_inv.beta_a_backscat)

            #water
            #compute mode diameter, extinction cross section, and backscatter efficeincy
            #from radar and lidar backscatter cross sections using mie theory and assumed
            #size distribution to predict mode diameter and q_backscatter for points
            #identified as water.
            if particle_parameters['radar_model'] == "Mie":
                rs_particle.mode_diameter, clipped_beta_a, rs_particle.q_backscatter\
                          ,rs_particle.dstar \
                      = size_dist.dmode_from_radar_lidar_mie(copy_radar_backscatter\
                          ,clipped_beta_a_backscat)
             
            else:
                #use only Rayliegh approx solution--particle_parameter['radar_model']=="Rayleigh"
                #mode diameter is computed for all points assuming everything is water
                #subsequent calculation will replace ice phase points.
                if particle_parameters['ext_source'] == 'ext':
                    clipped_beta_a = rs_inv.extinction_aerosol.copy()
                elif particle_parameters['ext_source'] == 'bs/p180':
                    clipped_beta_a = clipped_beta_a_backscat/particle_parameters['p180_water']
                else:
                    print 'particle_parameters=',particle_parameters['ext_source'],' not supported'
                    print 'in spheroid_particle_processing'
                    print j
                clipped_beta_a[np.isnan(clipped_beta_a_backscat)]=np.NaN
                phase = np.zeros_like(rs_inv.beta_a_backscat)
                zeta  = np.ones_like(rs_inv.beta_a_backscat)
                
                rs_particle.mode_diameter = size_dist.dmode_from_lidar_radar_rayleigh(
                           rs_particle.mode_diameter
                           ,clipped_beta_a
                           ,copy_radar_backscatter
                           ,zeta
                           ,phase)
            
            #ice
            #compute extinction cross section for ice points using backscatter phase function
            clipped_beta_a[rs_particle.phase==1] = \
                clipped_beta_a_backscat[rs_particle.phase==1]/particle_parameters['p180_ice']
            zeta = np.zeros_like(clipped_beta_a)
            zeta[rs_particle.phase==1] = particle_parameters['zeta']

            
            #derive mode_diameter directly from radar backscatter and lidar extinction
            #cross sections for parts of image populated by ice
            rs_particle.mode_diameter[rs_particle.phase==1] = size_dist.dmode_from_lidar_radar_rayleigh(\
                rs_particle.mode_diameter[rs_particle.phase==1] \
                ,clipped_beta_a[rs_particle.phase==1],copy_radar_backscatter[rs_particle.phase==1]\
                ,zeta[rs_particle.phase==1],rs_particle.phase[rs_particle.phase==1])


            #creates effective_diameter_prime array from mode diameter
            rs_particle.effective_diameter_prime = \
              size_dist.deff_prime(rs_particle.mode_diameter,rs_particle.phase,zeta)    

            rs_particle.effective_diameter = size_dist.eff_diameter(\
                                  rs_particle.mode_diameter,rs_particle.phase)
            
            rs_particle.mean_diameter = size_dist.mean_diameter(\
                                  rs_particle.mode_diameter,rs_particle.phase)

            #compute liquid water content for bins with phase == 0
            #bins with phase > 0 will return with NaN's
            if particle_parameters['radar_model'] == "Mie":
                 rs_particle.LWC = su.liquid_water_content_mie(rs_particle.effective_diameter
                       ,clipped_beta_a,rs_particle.q_backscatter)
                 rs_particle.p180_extinction = rs_inv.beta_a_backscat / rs_particle.q_backscatter
                
            else:
                if particle_parameters['ext_source'] == 'bs/p180':
                    rs_particle.extinction_aerosol= rs_inv.beta_a_backscat / particle_parameters['p180_water']
                    clipped_beta_a = rs_particle.extinction_aerosol.copy()
                else:    
                    clipped_beta_a = rs_inv.extinction_aerosol.copy()        
                clipped_beta_a[np.isnan(clipped_beta_a_backscat)]=np.NaN
                rs_particle.LWC = np.NaN * np.zeros_like(rs_particle.effective_diameter)
                su.liquid_water_content_ext_approx(rs_particle.LWC,rs_particle.effective_diameter
                       ,clipped_beta_a,rs_particle.phase)
                rs_particle.p180_extinction = rs_inv.beta_a_backscat / particle_parameters['p180_water']
            
            rs_particle.extinction_aerosol = rs_inv.extinction_aerosol.copy()
            #compute ice water water content for bins with phase > 0 (kg/m^3)
            #return in LWC array bins with phase > 0
            su.ice_water_content(rs_particle.LWC,rs_particle.effective_diameter
                     ,clipped_beta_a,rs_particle.phase)

            

            if hasattr(rs_radar,'vertically_averaged_doppler'):
                rs_radar.raw_MeanDopplerVelocity = rs_radar.MeanDopplerVelocity.copy()
                motion_correction = np.transpose(rs_radar.vertically_averaged_doppler\
                                    *np.transpose(np.ones_like(rs_radar.MeanDopplerVelocity))) 
                rs_radar.MeanDopplerVelocity -= motion_correction                              
                              
            if sounding!=None:
                s_time = datetime.utcnow()
                
                rs_particle.rw_fall_velocity,rs_particle.mw_fall_velocity \
                     ,rs_particle.model_spectral_width,rs_particle.nw_fall_velocity\
                     = su.weighted_fall_velocity( 
                     rs_particle.mode_diameter
                    ,particle_parameters
                    ,rs_particle.zeta
                    ,sounding.temps
                    ,sounding.pressures
                    ,rs_particle.phase,size_dist)                                                    
                print 'time for fall_velocity = ',datetime.utcnow() - s_time

            # compute precip rate (m/s) #rain_rate = 1/density
            # (m^3/kg) * LWC (kg/m^3) * fall_velocity (m/s) #using
            # Doppler velocity rs_particle.hsrl_radar_dv_precip_rate =
            # 0.001 * rs_particle.LWC * rs_radar.MeanDopplerVelocity

            #using raw Doppler to give precip rate in m/s
            rs_particle.hsrl_radar_dv_precip_rate = 0.001 * rs_particle.LWC * rs_radar.MeanDopplerVelocity
            #remove points with Doppler folding
            rs_particle.hsrl_radar_dv_precip_rate[rs_radar.MeanDopplerVelocity < -2.0] = np.NaN

            if sounding!=None:
                #using modeled mass weighted velocity and dividing by the density of water,
                #                              1000 kg/m^3, to give precip_rate in m/s
                rs_particle.hsrl_radar_precip_rate = 0.001 * rs_particle.LWC * rs_particle.mw_fall_velocity
            
           
            #retype all these fields to a proper TZ_Array
            #for f in ['effective_diameter_prime']:
            for f in vars(rs_particle).keys():
                v=getattr(rs_particle,f)
                if isinstance(v,hau.Z_Array):
                    continue #properly typed. continue
                elif isinstance(v,np.ndarray):
                    if len(v.shape)==2:
                        setattr(rs_particle,f,hau.TZ_Array(v))
                    elif len(v.shape)==1:
                        print '1 Dimensional Variable '+f+' will be changed to a T_Array. FIXME!!!!'
                        setattr(rs_particle,f,hau.T_Array(v))
                    else:
                        raise RuntimeError("I don't know what to type particle array "+f+" with dimensions "+repr(v.shape))
                else:
                    pass #not an array type. should be safe to ignore
            """
            #compute fitting error of computed radar weighted fall velocity
            #to measured Doppler veleocity.
            temp = rs_radar.Backscatter.copy()
            temp[np.isnan(rs_radar.Backscatter)]=0.0
            temp[rs_inv.msl_altitudes>400]=0.0
            fitting_error = np.sqrt(nanmean((rs_particle.rw_fall_velocity[temp >1e-9] \
                                  -  rs_radar.MeanDopplerVelocity[temp >1e-9])**2))
            print
            print rs_radar.times[0],'  --->  ' ,rs_radar.times[-1]
            print 'fitting_error (m/s)= ',fitting_error
            print
            """

            'rs_particle--spher'
            print dir(rs_particle)
            return rs_particle
