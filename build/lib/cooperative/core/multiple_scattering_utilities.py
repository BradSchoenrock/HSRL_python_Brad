import numpy as np
import lg_base.core.array_utils as hau
import scipy.special as ss
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

def ms_sums(diameter,beta,depol,multiple_scatter_parameters,ms_obj):
        
        '''  ms_sums(diameter,beta,depol,multiple_scatter_parameters,ms_obj):

             diameter = particle diameter profile or profiles
             beta     = extinction or backscatter cross section profile or profiles
             depol    = depolarization profile or profiles
             returns:
             ms_obj   = summed profiles
               "       .diameter_ice
               "       .diameter_water
               "       .beta_ice
               "       .beta_water
               "       .n_samples_ice
               "       .n_samples_water  '''
        is_ice =np.zeros_like(beta)
        is_water = np.zeros_like(beta)
        is_ice[depol >= multiple_scatter_parameters['h2o_depol_threshold']]= 1.0
        is_water[depol < multiple_scatter_parameters['h2o_depol_threshold']] = 1.0

        #diameter_ice = diameter * is_ice
        #diameter_water = diameter * is_water
        beta_ice = beta * is_ice
        beta_water = beta * is_water
       
        #compute diameter * beta for ice and water components
        if not diameter == None:
            diameter_ice = diameter * beta * is_ice
            diameter_water = diameter * beta * is_water
        #if diameter is not supplied from lidar-radar retrieval get it from constants in multiple_scatter_parameters
        else:     
            diameter_ice = hau.Z_Array(np.ones_like(beta) * multiple_scatter_parameters['mode_diameter_ice'] * beta * is_ice)
            diameter_water = hau.Z_Array(np.ones_like(beta) * multiple_scatter_parameters['mode_diameter_water'] * beta * is_water)
        
        if ms_obj == None:
            ms_obj = hau.Time_Z_Group()
            setattr(ms_obj,'beta_ice',hau.Z_Array(nansum(beta_ice ,0)))
            setattr(ms_obj,'beta_water',hau.Z_Array(nansum(beta_water ,0)))
            setattr(ms_obj,'diameter_ice',hau.Z_Array(nansum(diameter_ice,0)))
            setattr(ms_obj,'diameter_water',hau.Z_Array(nansum(diameter_water,0)))
            setattr(ms_obj,'n_samples_ice',hau.Z_Array(sum(~np.isnan(beta_ice)*is_ice ,0)))
            setattr(ms_obj,'n_samples_water',hau.Z_Array(sum(~np.isnan(beta_water)*is_water,0)))
        else:
            ms_obj.beta_ice  += nansum(beta_ice,0)
            ms_obj.beta_water += nansum(beta_water,0) 
            ms_obj.diameter_ice +=  nansum(diameter_ice,0)
            ms_obj.diameter_water += nansum(diameter_water,0)
            ms_obj.n_samples_ice += sum(~np.isnan(beta) * is_ice,0)
            ms_obj.n_samples_water += sum(~np.isnan(beta) * is_water,0)
       
        return ms_obj

def msinteg(N1,N2,firstR,lastR,inc,beta,mode_dia 
             ,ranges,wavelength,multiple_scattering_parameters,calvals):
    """
       msinteg(N1,N2,firstR,lastR,inc,beta,mode_dia 
             ,ranges,multiple_scattering_parameters,calvals,wavelength)
       N1         = highest order of scattering to compute.
       N2         = lowest order of scattering to compute.
       firstR     = begining of the range interval to compute (meters),
                    optical depth below firstR does not contribute to computed result.
       lastR      = end of the range interval to compute (meters).
       inc        = number of range indices between ms computations (normally inc=1).
       beta       = column vector with scattering cross sections (1/m), 
       mode_dia   = column vector with mode diameter at each range (m)
                    if supplied as single value it will be expanded to a column vector.
       ranges     = column vector containing the ranges at which beta_c, beta_R, and
                    mode_dia are specified. Answers will be returned at these points. 
                    ranges must be equally spaced. Ranges must be supplied in meters,
       wavelength = lidar wavelength (m)             
       multiple_scatter_parameters['particle_info_source'] == 'measured' | 'constant'
                                  ['p180_water'] 
                                  ['p180_ice'] 
                                  ['h2o_depol_threshold']  water when depol < threshold
                                  ['alpha_water'] 
                                  ['alpha_ice']
                                  ['g_water'] 
                                  ['g_ice'] 
                                  where: N(D) ~ D^alpha * exp(-alpha/g *(D/mode_dia)^g)
       calvals =  dictionary of calibration constants for particular hsrl
       

       #old doc header
       function ms=msinteg(id,N1,N2,firstR,lastR,inc,beta_c,r_particle...
             ,ranges,wavelength,FOVT,d_beam,FOVR,d_rec,P180,P18_2,P18_n);


       Python only version of mutliple scatter code-uses ms_int.py a recursive
       integration routine. The mode diameter of the gamma distribution is used to
       compute the square of the Gaussian diffraction peak width.
      
       
       See 'http://lidar.ssec.wisc.edu/mscat/derivations' for details. 
 
       id         = integer used in output file name,usually the starting shot number.
       N1         = highest order of scattering to compute.
       N2         = lowest order of scattering to compute.
       firstR     = begining of the range interval to compute (meters),
                    optical depth below firstR does not contribute to computed result.
       lastR      = end of the range interval to compute (meters).
       inc        = number of range indices between ms computations (normally inc=1).
       ranges     = column vector containing the ranges at which beta_c, beta_R, and
                    mode_dia are specified. Answers will be returned at these points. 
                    ranges must be equally spaced. Ranges must be supplied in meters, 
       beta       = column vector with scattering cross sections (1/m), 
       mode_dia   = column vector with mode diameter at each range (m)
                    if supplied as single value it will be expanded to a column vector.
       alpha,gam  = gamma size parameters, N(D)=D^alpha *exp(-b*D^gam)      
       lambda     = wavelength (m)
       FOVT       = laser divergence in radians.
       d_beam     = 1/e full with of the laser beam at telecope exit
       FOVR       = reciever fov in radians.
       d_rec      = diameter of receiving mirror (m)
       P180       = backscatter phase function for single scattering.
       P18_2      = backscatter phase function for second order scattering
       P18_n      = backscatter phase function for all higher orders of scattering.
       beta_source is a string which will be appended to field name of results.
       it is expected to indicate the algorithm used to compute the extinction
       cross section profile used in the calculations.

    """



    #convert angles to half angles--inputs in form of full angle divergence
    #and program in terms of half angles.
  
    FOVR = calvals['telescope_fov'] / 2.0
    FOVT = calvals['laser_beam_divergence'] / 2.0
    r_beam = calvals['1/e_laser_beam_width'] / 2.0
    r_rec  = calvals['telescope_diameter'] / 2.0
    #wavelength = calvals['wavelength']*1e-9
    alpha = multiple_scattering_parameters['alpha_water']
    gam = multiple_scattering_parameters['g_water'] 
   
   
    nranges=ranges.shape[0]
    
    #find index of first and last range involved in this computation
    first_pt = 1
    for i in range(nranges):
        if ranges[i] <= firstR:
             first_pt = i
        if ranges[i] <= lastR:
	     last_pt = i

    #convert mode diameter to diffraction peak width squared
    #area weighted diameter-squared from mode diameter of gamma distribution
    ave_D_squared = mode_dia**2 * (gam/alpha)**2 * ss.gamma((alpha+3)/gam) /ss.gamma((alpha+1)/gam)
    #theta_s = (4*wavelength^2/(pi^2 *<D^2>)
    theta_sq = (2 * wavelength/np.pi)**2 / ave_D_squared

    #if mode_dia is provided as a single value expand into array
    #if not isinstance(mode_dia,hau.Z_Array):
    #    if not isinstance(mode_dia,np.ndarray):
    #        theta_sq = theta_sq * np.ones_like(ranges)
   
    nranges = len(ranges)
    beta_c = beta.copy()
    dr = ranges[1] - ranges[0]

    # cloud optical depth
    tau_c = np.zeros_like(ranges)
    beta_c[np.isnan(beta_c)]=0.0
    tau_c[first_pt:last_pt] = np.cumsum(beta_c[first_pt:last_pt])-beta_c[first_pt]/2.0 \
                               -beta_c[first_pt:last_pt]/2.0
    tau_c = tau_c * dr


    #sum of Rayleigh and cloud optical depth
    tau_t = np.cumsum(beta_c[first_pt:last_pt])-beta_c[first_pt]/2.0 \
            -beta_c[first_pt:last_pt]/2.0
    tau_t = tau_t * dr



    
    tbeta = np.transpose(beta_c)
    #time0=cputime;
    #ms = np.zeros((ranges.shape[0],N2-N1+3))
    #compute for subset of ranges if inc > 1
    #inc =3
    indices = range(first_pt,last_pt,inc)
    ms=np.zeros((len(indices),N2+1))
    ms[:,0]=ranges[indices]
    s_time =datetime.utcnow()
    for n_th in range(N1,N2+1):
        j=0
        for i in indices:
            divlR2=(FOVT * ranges[i]+r_beam)**2
            fovR2 = (FOVR * ranges[i]+r_rec)**2
            ms[j,n_th] = ms_int(n_th,ranges \
                           ,first_pt,i, tbeta,theta_sq,divlR2,fovR2)
            ms[j,n_th]=tau_c[i]**(n_th-1) / ss.gamma(n_th) - ms[j,n_th]
            #print 'Range=%g N=%g tau=%g,ms=%g' %(ranges[i],n_th,tau_c[i]**(n_th-1)/ss.gamma(n_th),ms[j,n_th])
            
            j=j+1
    #interpolate back to input resolution
    #ms_ratios_profile[number_of_range_bins,number_of_orders_of_scattering + 1]
    ms_ratios_profile = np.zeros((beta.shape[0],N2+1))
    ms_ratios_profile[first_pt:last_pt,0]=ranges[first_pt:last_pt]
    for i in range(2,N2+1):
        ms_ratios_profile[first_pt:last_pt,i] = np.interp(ranges[first_pt:last_pt],ranges[indices],ms[:,i])
    ms_ratios_profile = hau.Z_Array(ms_ratios_profile)    
    print 'time for multiple scatter cal = ',datetime.utcnow() - s_time        
    return ms_ratios_profile


def ms_int(n_th,ranges,first_pt,last_pt,beta_c,theta_sq,divlR2,fovR2):
    """ms_size_int(n_th,ranges,first_pt,last_pt,beta_c,divlR2,fovR2)

       Compute the multiple integral for nth order scattering with a recursive
       python routine--note this does not contain the tau^(n-1)/tau(n-1)! term 
       or factors for P(180,n)/P(180) and overlap function.  
       See 'http:\\lidar.ssec.wisc.edu/mscat/derivation' for details.
  
       n_th    = order of scattering.
       ranges  = range vector (column vector, in meters).
       first_pt= range index at which to start computation.
                 optical depth below firstR does not contribute to computed result.
       last_pt = range index of the turn-around slab (ie range to single scatter).
       beta_c  = scattering cross section at each range (column vector).
       
       divergence * range to single scatter/(2*lambda))^2.  q=
       (pi*divl*ranges(last_pt)/(lambda*2))^2; diveregence in radains
       and range in meters, lambda in microns.  fovR2 = (pi*full
       field-of-view * range to single scatter/(2*lambda))^2.  fovR2 =
       (pi*fov*ranges(last_pt)/2)^2; fov in radians and range in
       meters, lambda in microns.
    """
    
    
    dr = ranges[1]-ranges[0]
    npts = last_pt-first_pt+1
    kern = np.zeros_like(ranges)
    if n_th == 2:        #do innermost range and particle size integrals
        for i in range(first_pt,last_pt):     # integrate over range                               
            kern[i-first_pt] = np.exp(-fovR2 \
                    /(theta_sq[i] * (ranges[last_pt]-ranges[i])**2  + divlR2))
          
          
        #no contribution at first_pt
        ms_integral = dr * (nansum(beta_c[first_pt:last_pt] * kern[0:npts-1]) \
            -beta_c[last_pt] * kern[npts-1]/2)        
        
    else:               #do one of outer integrals
        n_th = n_th-1
        sint = np.zeros_like(ranges)        
        for i in range(first_pt,last_pt):    
             s = divlR2 + theta_sq[i] * (ranges[last_pt]-ranges[i])**2 
             sint[i] = ms_int(n_th,ranges,i,last_pt,beta_c,theta_sq,s,fovR2)

        first_term = beta_c[first_pt] * sint[first_pt]/2.0
        if np.isnan(first_term):
            first_term = 0.0
        last_term  = beta_c[last_pt]*sint[last_pt]/2.0
        if np.isnan(last_term):
            last_term = 0.0
        ms_integral = dr * (nansum(beta_c[first_pt:last_pt] * sint[first_pt:last_pt])\
            -first_term - last_term)
            #-beta_c[first_pt] * sint[first_pt]/2.0-beta_c[last_pt]*sint[last_pt]/2.0) 
    
    return ms_integral
  
if __name__ == '__main__':
    try:
         ranges = np.arange(0,1000.0,7.5)
         beta_R = np.ones_like(ranges)*1e-6
         beta_c = np.zeros_like(ranges)
         beta_c[10:50] = 1e-3
         mode_dia  = np.ones_like(ranges) * 50.0e-6
         wavelength =532e-9
         FOVT = 50e-6
         FOVR = 100.0e-6
         P180 = 0.05
         P18_2 = 0.05
         P18_n = 0.05
         beta_source = 'test'
         [z,ms] = msinteg('test',2,4,0.0,600.0,1,beta_c,beta_R,mode_dia,alpha,gam 
             ,ranges,wavelength,FOVT,0.2,FOVR,0.4,P180,P18_2,P18_n,beta_source)
         if 1:
             print
             print "range   2'nd     3'rd    4'th     5'th"
             for i in range(z.shape[0]):
                 print '%4.1f   %6.3f ' %(z[i],ms[i,2]),
                 j=3
                 while j < ms.shape[1]:
                    print '%6.3f' %(ms[i,j]),
                    j = j +1
                 print 
             import matplotlib.pylab as plt
             n_order = ms.shape[1]
             colors = ['b','g','r','k','c']
             legend_num =['2','3','4','5','6']
             plt.figure(1)
             for order in range(ms.shape[1]):
                if order >1: 
                  plt.plot(ms[:,order],z,color=colors[order-2])
             plt.legend(legend_num[:order-1])
             plt.grid(True)
             plt.ylabel('Altitude')
             plt.xlabel('multiple/single ratio')
             plt.show()
       
    except RuntimeError, msg:
        print 'WARNING: ', msg
    except OSError, msg:
        print 'WARNING: ', msg


