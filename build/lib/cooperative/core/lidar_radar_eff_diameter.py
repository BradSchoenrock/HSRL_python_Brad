import numpy as np
import scipy.special as gm
import scipy.interpolate as scinterp
import lg_base.core.decoratortools as nt
import logging

LOG = logging.getLogger(__name__)
"""
def lidar_radar_eff_diameter_prime(beta_a_backscat,radar_backscat,depol
                     ,h2o_depol_threshold,beta_a=None,p180_ice=None,lambda_radar=None):
       Calculate effective diameter prime from HSRL and Radar Backscatters
    
    :param beta_a_backscat: lidar backscatter cross section, (1/(m sr))
    :param radar_backscat: radar backscatter cross section, (1/(m sr))
    :param depol: lidar measured depolarization, this can be either linear or circular depolarization as long as h2o_depol_threshold is set for proper pol type.
    :param beta_a: scattering cross section array array computed from measured od.
    :param p180_ice: p(180)/4pi = ice crystal backscatter phase function if p180_ice=None (ie no number is supplied) beta_a computed from the measured od will be used, otherwise beta_a will be set equal to beta_a_backscat/p180_ice. This is useful because the measured beta_a is sensitive to averaging lengths and times, with short averages apt to have excessive noise and long averages apt to mix regions of low and high scattering the computed array can be replaced by a user supplied constant value. It is recomended that this be based on an independent analysis of the data to determine an appropriate value for p(180)/4pi in ice regions of the cloud. 
    :param h2o_depol_threshold: depol < h2o_depol_threshold indicates liquid water cloud
    :param eff_diameter_prime: lidar-radar effective diameter (microns) =(9*<Volume^2>/(pi*<Area>))^.25 where <Volume^2> and <Area> refer to particle averages as defined by Donovan and Lammeren JGR,v106, D21, pp27425
                    
 Assumptions:

      - cloud at a given data point is either all water or all ice
      - when depolarization is < h2o_depol_threshold, cloud is water
      - backscatter phase function for water droplets, p(180)/4pi = 0.05
      - particles are large compared to the lidar wavelength such that the optical scattering cross section is twice the projected area of the particle.
      - particles are small compared to the radar wavelength
      - radar attenuation is negligible
    

    assert(lambda_radar!=None)
    k_sq_water=0.93;  #dielectric constant squared for water
    k_sq_ice=  .176;   #dielectric constant squared for ice

    # set separate constants for water and ice regions.
    #dielectric constant
    k_sq_array=k_sq_ice*np.ones_like(beta_a_backscat)
    k_sq_array[depol< h2o_depol_threshold]=k_sq_water

    radar_factor_array=3*(lambda_radar**4.0)/(64.0*(np.pi**4.0)*k_sq_array)
     
     
    #eff_diameter_prime is the fundamental size dependent quantity derived 
    #from lidar and radar scattering cross sections.
    #<volume of particle^2>/<projected area of particle> = eff_diameter_prime^4 *pi/9
    #see Donovan and van Lammeran, JGR, Vol 106 no D21, Nov 16, 2001.
    #for more on this definition which we state in terms of diameter rather
    #than raduis. In the special case of monodisperse spherical particles
    #the effective_diameter_prime is equal to the effective diameter. 
    #In general, the converstion to effective diameter is a function of size
    #distribution parameters and ice crystal shape.
     

    beta_a_backscat[beta_a_backscat<=0]=np.NaN
    radar_backscat[radar_backscat<=0]=np.NaN
     
     
    eff_diameter_prime=2.0*((8.0*np.pi/3.0)*radar_factor_array*radar_backscat/beta_a)**.25; 
    eff_diameter_prime[eff_diameter_prime<1e-6]=-np.inf;
    
    return eff_diameter_prime,beta_a
    """
def d_prime_to_eff_table(ice_distribution):
    """ diameter prime to effective diameter prime table
     Sucessive values in the output arrays correspond to increasing values of the parameter (1/b) in the gamma size distribution defined in "ice_distribution" using power law descriptions of crystal volume and area as a function of particle size. size distribution of ice crystals
           N(D) ~ D^alpha * exp(-b*D^g_ice)
     where D=diameter of smallest sphere enclosing ice crystal and where we shall assume that D and Dr are expressed in microns.

     Projected area of ice crystals
       for D < Dr;    A=sigma_a*Dr^(2-delta_a1)*D^delta_a1   
           D>= Dr;    A=sigma_a*Dr^(2-delta_a2)*D^delta_a2
     Volume of ice in crystal
           D < Dr;    V=sigma_v*Dr^(3-delta_v1)*D^delta_v1
           D>= Dr;    V=sigma_v*Dr^(3-delta_v2)*D^delta_v2
     ave_area= ave area of ice particle corresponding to d_eff_prime (microns^2).
     sample call: ::

    s=struct('alpha_ice',2,'g_ice',1,'sigma_a',1,'delta_a1',2,'delta_a2',2,'sigma_v',1,'delta_v1',3,'delta_v2',3,'h2o_depol_threshold',.1,'Dr',100);[dep,de]=d_prime_to_eff_table(s);

    :structure ice_distribution
    -   .alpha_ice
    -   .g_ice    
    -   .sigma_a
    -   .delta_a1  
    -   .delta_a2 
    -   .sigma_v
    -   .delta_v1  
    -   .delta_v2 
    -   .h2o_depol_threshold  (range [0,1])
    -   .Dr       
    """

    alpha_ice=ice_distribution.alpha_ice;
    g_ice=ice_distribution.g_ice;
    sigma_a=ice_distribution.sigma_a;
    delta_a1=ice_distribution.delta_a1;
    delta_a2=ice_distribution.delta_a2;
    sigma_v=ice_distribution.sigma_v;
    delta_v1=ice_distribution.delta_v1;
    delta_v2=ice_distribution.delta_v2;
    h2o_depol_threshold=ice_distribution.h2o_depol_threshold;
    Dr=ice_distribution.Dr;


    #find maximum value of x needed to cover particles up to 1 cm.
    xmax=100;
    deff_prime_max=0;
    while deff_prime_max<10000:
        [deff_prime_max,ave_area]=compute_deff_prime(alpha_ice,g_ice,sigma_a,delta_a1,delta_a2 
                    ,sigma_v,delta_v1,delta_v2,Dr,xmax);
        xmax=xmax*1.2;
    LOG.debug(xmax)

    # make vector of x=(1/b)*Dr^gamma
    x=np.logspace(1, np.log(xmax)/np.log(10), 200);

    #clear x
    #xmax=max(3,log10(1000^g_ice))
    # x=logspace(-2, 2*xmax, 400);


    [deff_prime,ave_area]=compute_deff_prime(alpha_ice,g_ice,sigma_a,delta_a1,delta_a2 
                    ,sigma_v,delta_v1,delta_v2,Dr,x);

    [deff,d_mean]=compute_deff(alpha_ice,g_ice,sigma_a,delta_a1,delta_a2 
                    ,sigma_v,delta_v1,delta_v2,Dr,x);
  

    if 0:
      pass
      #figure(2)
      #hold off
      #plot(deff_prime_raw,'+b')
      #hold on
      #plot(deff_raw*1.05,'+r')
      #grid
 
    #deff=interp(deff_prime_raw,deff_raw,1:500,'linear',nan);  

    if 0:
      pass
      #figure(3)
      #plot(deff_prime_raw,deff_raw,'+',1:500,deff,'r')
      #ax=axis;
      #axis([1 500 0 deff(500)])
      #grid

    if 0:
        pass
        #figure(1)
        #max(deff_prime)
        #max(deff)
        #max(ave_area)
        #plot(deff_prime,deff/deff_prime,deff_prime,ave_area/((pi/4)*deff**2))
        #ylabel('d_e_f_f / d_e_f_f prime')
        #xlabel('d_e_f_f prime')
        #ax=axis;
        #axis([1 500 0 ax(4)*1.1])
        #h=line([Dr Dr], [0 ax(4)*1.1]);
        #set(h,'color','r')
        #grid

    #b_table=deff;
    return deff_prime,deff,ave_area,d_mean


def compute_deff_prime(alpha_ice,g_ice,sigma_a,delta_a1,delta_a2,sigma_v,delta_v1,delta_v2,Dr,x):


    alpha_a1_g_ice=(alpha_ice+delta_a1+1)/g_ice;
    alpha_a2_g_ice=(alpha_ice+delta_a2+1)/g_ice;

    t_k=(1/x)*(Dr**g_ice); #coordinate transform for limit of integration

    gamma_fun_a1=gm.gamma(alpha_a1_g_ice);
    inc_gamma_fun_a1=gm.gammainc(t_k,alpha_a1_g_ice);
    gamma_fun_a2=gm.gamma(alpha_a2_g_ice);
    inc_gamma_fun_a2=gm.gammainc(t_k,alpha_a2_g_ice);

    alpha_v1_g_ice=(alpha_ice+2*delta_v1+1)/g_ice;
    alpha_v2_g_ice=(alpha_ice+2*delta_v2+1)/g_ice;



    gamma_fun_v1=gm.gamma(alpha_v1_g_ice);
    inc_gamma_fun_v1=gm.gammainc(t_k,alpha_v1_g_ice);
    gamma_fun_v2=gm.gamma(alpha_v2_g_ice);
    inc_gamma_fun_v2=gm.gammainc(t_k,alpha_v2_g_ice);


    A=Dr**(-2*delta_v1)*inc_gamma_fun_v1*gamma_fun_v1;
    B=Dr**(-2*delta_v2)*(1-inc_gamma_fun_v2)*gamma_fun_v2;
    C=Dr**(-delta_a1)*inc_gamma_fun_a1*gamma_fun_a1;
    G=Dr**(-delta_a2)*(1-inc_gamma_fun_a2)*gamma_fun_a2;

    area_factor=(C*x**(delta_a1/g_ice)+G*x**(delta_a2/g_ice));

    deff_prime=(A*x**(2*delta_v1/g_ice)+B*x**(2*delta_v2/g_ice))/area_factor;

    deff_prime=Dr*deff_prime**0.25*(sigma_v**2/sigma_a)**0.25;


    #compute average area of particle in square microns.
    ave_area=sigma_a*(np.pi/4)*Dr**2*area_factor/gm.gamma((alpha_ice+1)/g_ice);
    ave_area *=1e-12
    LOG.debug('ave_area')
    LOG.debug(ave_area)
    
    return deff_prime,ave_area


def compute_deff(alpha_ice,g_ice,sigma_a,delta_a1, delta_a2,sigma_v,delta_v1,delta_v2,Dr,x):

    t_k=(1/x)*Dr**g_ice;

    gamma_fun_a1=gm.gamma((alpha_ice+delta_a1+1)/g_ice);
    inc_gamma_fun_a1=gm.gammainc(t_k,(alpha_ice+delta_a1+1)/g_ice);
    gamma_fun_a2=gm.gamma((alpha_ice+delta_a2+1)/g_ice);
    inc_gamma_fun_a2=gm.gammainc(t_k,(alpha_ice+delta_a2+1)/g_ice);

    gamma_fun_v1=gm.gamma((alpha_ice+delta_v1+1)/g_ice);
    inc_gamma_fun_v1=gm.gammainc(t_k,(alpha_ice+delta_v1+1)/g_ice);
    gamma_fun_v2=gm.gamma((alpha_ice+delta_v2+1)/g_ice);
    inc_gamma_fun_v2=gm.gammainc(t_k,(alpha_ice+delta_v2+1)/g_ice);

    A=Dr**(-delta_v1)*inc_gamma_fun_v1*gamma_fun_v1;
    B=Dr**(-delta_v2)*(1-inc_gamma_fun_v2)*gamma_fun_v2;
    C=Dr**(-delta_a1)*inc_gamma_fun_a1*gamma_fun_a1;
    G=Dr**(-delta_a2)*(1-inc_gamma_fun_a2)*gamma_fun_a2;



    deff=(A*x**(delta_v1/g_ice)+B*x**(delta_v2/g_ice))/(C*x**(delta_a1/g_ice)+G*x**(delta_a2/g_ice));
    deff=deff*Dr*sigma_v/sigma_a;

    #for ice regions
    #dmean=x.^g_ice*gamma((alpha_ice+2)/g_ice)/gamma((alpha_ice+1)/g_ice);
    #corrected 28-Nov-07 ewe
    d_mean=x**(1/g_ice)*gm.gamma((alpha_ice+2)/g_ice)/gm.gamma((alpha_ice+1)/g_ice);


    return deff,d_mean

def d_eff_from_d_eff_prime(eff_diameter_prime,beta_a_backscat,depol,beta_a,pparm):
    """ Effective Diameter from Effective Diameter Prime

    :param beta_a_backscat: lidar backscatter cross section, 1/(m sr)
    :param beta_a:          lidar scattering cross section, 1/m
    :param depol:           lidar measured circular depolarization. 
    :param pparm:           structure containing the following elements(shown with sample values) 
              alpha_ice: 2
                  g_ice: 1
                sigma_a: 1
               delta_a1: 2
               delta_a2: 2
                sigma_v: 1
               delta_v1: 3
               delta_v2: 3
    h2o_depol_threshold: 0.1000
                     Dr: 100
            alpha_water: 2
                g_water: 1
               p180_ice: 0.0350
   where:
    h2o_depol_threshold = linear depolarization threshold, linear_depol<threshold identified as liquid. 
    alpha_water         = user supplied size distribution shape parameter alpha for liquid
    alpha_ice           = user supplied size
    g_water             = user supplled size distribution parameter g for liquid
    g_ice               = user supplied size distribution shape parameter g for ice
    p180_ice            = p(180)/4pi = backscatter phase function for ice crystals this can be left empty and beta_a will be used directly rather than beta_a_backscat/p180_ice.
    sigma_a, delta_a    = projected area of ice crystal = sigma_a*(pi/4)*D^delta_a where D is maximum dimension of ice crystal see Donovan and van Lammeren, JGR, Vol 106, pp 27433, no D21, Nov 16 2001. for spheres, sigma_a=1, delta_a=2
    sigma_v, delta_v   = volume of ice crystal = sigma_v*(pi/6)D^delta_a for spheres sigma_v=1, delta_v=3
    eff_diameter_prime = lidar-radar effective diameter (microns) = (9*<volume^2>/(pi*<projected area>))^.25
    LWC                = liquid water content gr/m^3
    mean_diameter      = mean value of particle diameter(ie mean of max dimension)
    eff_diameter       = effective diameter (microns) = 2*<volume>/(3*<area>) computed from lidar and radar backscatter cross sections assuming and a particle sizes given by a modified gamma distribution as presented by Deirmendjian, 'Electromagnetic Scattering on Spherical Polydispersions', Elsevier, NY, 1969:
                     n(D) = a * D^alpha * exp(-b*D^g)      (equation 1) 
    
                     Where:
                      D    = Maximum dimension of particle
                      a    = Num_particles*g*b^((alpha+1)/g)*gamma((alpha+1)/g)
                      n    = number of particles per unit volume
                      b    = parameter computed from lidar-radar signal ratio
                      gamma= the gamma function
     
 Assumptions:
      - Size distribution given by equation 1
      - cloud at a given data point is either all water or all ice
      - when depolarization is < h2o_depol_threshold, cloud is water
      - backscatter phase function for water droplets, p(180)/4pi = 0.05
      - particles are large compared to the lidar wavelength such that the optical scattering cross section is twice the projected area of the particle.
      - particles are small compared to the radar wavelength
      - radar attenuation is negligible
      - the volume of ice in a particle is related to its max dimension, D, by:                   
                 Volume = sigma_v * pi/6 *Dr^(3-delta_v)* D^delta_v    
                    for water, sigma_v =1, delta_v=3 thus Volume= pi/6 * D^3
      - the projected area of an ice particle is related to D by:
                 Area = sigma_a *pi/4 *Dr^(2-delta_a)* D^delta_a        
                    for water, sigma_a =1, delta_a=2, thus Area=(pi/4) * D^2
      - because the power laws for volume and area are often different for small and large particles two values can be specified for delta_a and delta_v:
           delta_a=delta_a1, delta_v=delta_v1   for D < Dr    (microns)
           delta_a=delta_a2, delta_v=delta_v2   for D>= Dr   
          
      - sigma_a=the area fill fraction at D=Dr, it is the projected area of an ice crystal of size D=Dr divided by (pi/4)*Dr^2.
      - sigma_v=volume fill fraction at D=Dr, it is the volume of ice within an ice particle of size D=Dr divided by (pi/6)*Dr^3.
    """

    
    h2o_depol_threshold=pparm.h2o_depol_threshold;

    p180_ice=pparm.p180_ice;
    alpha_water=pparm.alpha_water;
    g_water=pparm.g_water;

    #density_ice=0.92; #g/cm^3
    specific_gravity_ice = 0.92 

    #this returns vectors of d_prime and the corresponding values of d_eff
    #which can be used as a lookup table and covert from microns to meters.
    [d_prime_table,d_eff_table,ave_area_table,mean_diameter_table]=d_prime_to_eff_table(pparm); 
    #figure(1000);
    # plot(d_prime_table,mean_diameter_table);grid;xlabel('d_e_f_f prime');ylabel('d_e_f_f')
    # pause
    #convert from microns to meters
    d_prime_table=d_prime_table/1e6;
    d_eff_table=d_eff_table/1e6;
    ave_area_table=ave_area_table/1e12;
    mean_diameter_table=mean_diameter_table/1e6;

    #convert eff_diameter_prime from microns to meters
    #eff_diameter_prime=eff_diameter_prime/1e6;

    #convert d_eff_prime in those portions of the cloud dominated by ice 
    #to d_eff by interpolating lookup table
    eff_diameter=scinterp.interp1d(d_prime_table,d_eff_table,bounds_error=False)(eff_diameter_prime);


    #eff_diameter=eff_diameter1d.reshape((ntimes,nalts));
    #eff_diameter_prime=eff_diameter_prime.reshape((ntimes,nalts));


    if (p180_ice!=None and np.isfinite(p180_ice)) or beta_a==None:
        beta_a=beta_a_backscat/p180_ice;
      

    #compute beta_a in water clouds using a fixed value of p180.
    beta_a[depol< h2o_depol_threshold]=beta_a_backscat[depol<h2o_depol_threshold]/.05;
      

    #set size distribution factor in water regions
    c_factor=gm.gamma((alpha_water+4)/g_water)/gm.gamma((alpha_water+3)/g_water);
    c_factor=c_factor*(gm.gamma((alpha_water+3)/g_water)/gm.gamma((alpha_water+7)/g_water))**.25;
    eff_diameter[depol<h2o_depol_threshold]=c_factor*eff_diameter_prime[depol<h2o_depol_threshold]
     


    #assume particles large compared to lidar lidar wavelength
    #beta_a = num_particles * 2 * <area per particle>;

    #in ice regions
    area_per_particle=np.zeros_like(eff_diameter_prime);
    area_per_particle[depol>=h2o_depol_threshold]=scinterp.interp1d(d_prime_table,ave_area_table,bounds_error=False)(eff_diameter_prime[depol>=h2o_depol_threshold]);
    
    num_particles=beta_a/(2*area_per_particle); #1/m^3
    
      
    #eff_diameter = 3*<volume per particle>/(2*<area per particle>)
    #LWC   = num_particles * <volume per particle> * density of water
    #multiply by 10^6 to convert m^3 to cm^3 yielding gr/m^3
     
    #LWC = (2.0/3.0)*density_ice*num_particles * eff_diameter *area_per_particle *1e6;
    LWC = (2.0/3.0)*1000.0*specific_gravity_ice*num_particles * eff_diameter *area_per_particle
     
   
    
    #in water regions
    area_per_particle[depol<h2o_depol_threshold]=(np.pi/4.0)*(eff_diameter[depol<h2o_depol_threshold]*gm.gamma((alpha_water+3.0)/g_water)/
        (gm.gamma((alpha_water+4.0)/g_water)))**2*gm.gamma((alpha_water+3)/g_water)/gm.gamma((alpha_water+1.0)/g_water)
     
    num_particles[depol<h2o_depol_threshold]=beta_a_backscat[depol<h2o_depol_threshold]/(.05*2.0*area_per_particle[depol<h2o_depol_threshold]) #1/m^3
     
     
    #eff_diameter = 3*<volume per particle>/(2*<area per particle>)
    #LWC   = num_particles * <volume per particle> * density of water
    #multiply by 10^6 to convert m^3 to cm^3 yielding gr/m^3

    
    LWC[depol<h2o_depol_threshold]=(2.0/3.0)*1000.0*num_particles[depol<h2o_depol_threshold]* eff_diameter[depol<h2o_depol_threshold]*area_per_particle[depol<h2o_depol_threshold] 
     
    #compute mean diameter of particles at each point in image
    #ice regions of cloud
    #mean_diameter=np.zeros_like(eff_diameter_prime);
    mean_diameter=scinterp.interp1d(d_prime_table,mean_diameter_table,bounds_error=False)(eff_diameter_prime);
 
    
    #for water regions of cloud
    mean_dia_const=gm.gamma((alpha_water+3.0)/g_water)/gm.gamma((alpha_water+4.0)/g_water)*gm.gamma((alpha_water+2.0)/g_water)/gm.gamma((alpha_water+1.0)/g_water)
               
    mean_diameter[depol<h2o_depol_threshold]=mean_dia_const*eff_diameter[depol<h2o_depol_threshold];
    mean_diameter[ np.logical_not(np.isfinite(mean_diameter)) ]=np.nan;
    # mean_diameter*=1e6; #change units to microns
    
    LOG.debug(h2o_depol_threshold)
    
     
    #eff_diameter(~isreal(eff_diameter))=nan;      
    eff_diameter[ np.logical_not(np.isfinite(eff_diameter)) ]=-np.inf;
    #eff_diameter*=1e6; #change units to microns
    #num_particles/=1000.0; #change to number particles/liter.
    
    return eff_diameter,num_particles,LWC,mean_diameter
  
