
from collections import namedtuple
import lidar_radar_eff_diameter as lred
import lg_base.core.array_utils as hau
import numpy as np

def process_mass_dimension_particle(rs_inv,rs_radar,particle_parameters,lambda_radar,entire_frame):
            """
            generate and return the particle measurements based on a given hsrl inverted data, radar (and its lambda), and particle parameters dictionary
            """

            ParticleParameters=namedtuple('ParticleParameters',','.join(particle_parameters.keys()))
            pparm=ParticleParameters(**particle_parameters)#pparm is a structure of the particle parameters, instead of a dictionary 'particle_parameters'

            #create timez group and add heights            
            rs_particle=hau.Time_Z_Group(rs_inv.times.copy(),timevarname='times',altname='heights')
            setattr(rs_particle,'heights',rs_inv.msl_altitudes.copy())

            #remove points where lidar signal is noise dominated by setting to
            #very small value.
            clipped_beta_a_back=rs_inv.beta_a_backscat.copy()
            if hasattr(rs_inv,'std_beta_a_backscat'):
                clipped_beta_a_back[clipped_beta_a_back<(2*rs_inv.std_beta_a_backscat)]=-np.inf
            else:
                print 'No std_beta_a_backscat statistics to filter particle measurements'
            clipped_beta_a_back[np.logical_not(np.isfinite(rs_inv.beta_a_backscat))]=-np.inf;
            used_beta_a=None
            if hasattr(rs_inv,'beta_a'):
                used_beta_a=rs_inv.beta_a.copy()
                #mask beta_a?
          
            rs_particle.effective_diameter_prime,used_beta_a = \
                lred.lidar_radar_eff_diameter_prime(
                   beta_a_backscat=clipped_beta_a_back
                  ,radar_backscat=rs_radar.Backscatter
                  ,depol=rs_inv.linear_depol
                  ,h2o_depol_threshold=particle_parameters['h2o_depol_threshold']
                  ,beta_a=used_beta_a
                  ,p180_ice=particle_parameters['p180_ice']
                  ,lambda_radar=lambda_radar)
            
            rs_particle.effective_diameter,rs_particle.num_particles,rs_particle.LWC,rs_particle.mean_diameter=\
                 lred.d_eff_from_d_eff_prime(
                     rs_particle.effective_diameter_prime
                     ,clipped_beta_a_back
                     ,rs_inv.linear_depol
                     ,used_beta_a
                     ,pparm)


              
            #rs_particle.hsrl_radar_rain_rate = 3600 * 10* .0001 * rs_particle.LWC * rs_radar.MeanDopplerVelocity
            #convert to mks units LWC kg/m^3, Doppler m/s, rain_rate m/s
            rs_particle.hsrl_radar_rain_rate = 0.001 * rs_particle.LWC * rs_radar.MeanDopplerVelocity
             
            #retype all these fields to a proper TZ_Array
            for f in ['hsrl_radar_rain_rate','effective_diameter_prime','effective_diameter','num_particles','LWC','mean_diameter']:
                if hasattr(rs_particle,f):
                    setattr(rs_particle,f,hau.TZ_Array(getattr(rs_particle,f)))
            

            return rs_particle
