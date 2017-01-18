import numpy as np
import scipy.special as gm
import scipy.interpolate as scinterp

#def lidar_radar_eff_diameter_prime(lidar_scat,radar_scat,phase,lambda_radar):
def lidar_radar_eff_diameter_prime(lidar_scat,radar_backscat,phase,lambda_radar):
    """ Calculate effective diameter prime from HSRL and Radar Backscatters
    
    lidar_scat      = lidar scattering cross section,  1/m 
    radar_backscat  = radar backscatter cross section, 1/(m str)
    phase          = particle phase, array like beta_a, 0==water, 1= ice
    lambda_radar   = radar wavelength (m)

    returns:
    eff_diameter_prime = lidar-radar effective diameter (microns)
                       

 
    Assumptions:

      - cloud at a given data point is either all water or all ice
      - particles are large compared to the lidar wavelength such that the optical scattering
        cross section is twice the projected are
      - radar scattering is in the Rayliegh regime"""
    
    k_sq_water=0.93;  #dielectric constant squared for water
    k_sq_ice=  .176;   #dielectric constant squared for ice
    
    #dielectric constant for water or ice depending on phase
    k_sq_array = k_sq_water * np.ones_like(lidar_scat)
    k_sq_array[phase > 0] = k_sq_ice

    #value used before adding 3/2 factor in radar cross section
    #constant = lambda_radar * (3.0/4.0)**0.25 /np.pi
    #radar_factor_array = constant * k_sq_array**-0.25 

    #radar_factor_array = (lambda_radar/np.pi) *(2.0* k_sq_array)**-0.25 
    radar_factor_array = lambda_radar * (2.0/(k_sq_array *np.pi**3))**0.25 
     
    #eff_diameter_prime is the fundamental size dependent quantity derived 
    #from lidar and radar scattering cross sections.
    #<volume of particle^2>/<projected area of particle> = eff_diameter_prime^4 *pi/9
    #see Donovan and van Lammeran, JGR, Vol 106 no D21, Nov 16, 2001.
    #for more on this definition which we state in terms of diameter rather
    #than raduis. In the special case of monodisperse spherical particles
    #the effective_diameter_prime is equal to the effective diameter. 
    #In general, the converstion to effective diameter is a function of size
    #distribution parameters and ice crystal shape.
     
    lidar_scat[lidar_scat <=0]=np.NaN
    #radar_scat[radar_scat<=0]=np.NaN
    radar_backscat[radar_backscat <= 0] = np.NaN
    
    #eff_diameter_prime = radar_factor_array * (radar_scat/lidar_scat)**.25 
    eff_diameter_prime = radar_factor_array * (radar_backscat/lidar_scat)**0.25
    
    eff_diameter_prime[eff_diameter_prime<0.1e-6]=-np.inf;
    return eff_diameter_prime
         
