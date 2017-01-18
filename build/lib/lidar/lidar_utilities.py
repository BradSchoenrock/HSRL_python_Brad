
import numpy as np
import scipy
from scipy.optimize import leastsq, curve_fit
import lg_base.core.array_utils as hau
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

def geometry_correction(apply_to,rs,rs_cal,nbins,s_bin,apply_geo_corr=None):
    """geometry_correction(selected_vars,rs,rs_cal,nalts,s_bin,wfov_mol_ratio=Noneapply_geo_corr=None)
       apply_to       = apply geo correction to these variables
       rs             = structure containing channel counts
       rs_cal         = structure containing geo_corr data
       nbins          = number of altitude bins to process
       s_bin          = range bin containing laser pulse
       wfov_mol_ratio = ratio of molecular_wfov_counts to molecular counts (optional correction)
       apply_geo_corr = 
       """
   

    for field in apply_to:
        if hasattr(rs,field):
            ones_array = np.ones(getattr(rs,field).shape)
            break    
    if apply_geo_corr is None or apply_geo_corr > 0:
        geocorr = rs_cal.geo.geo_correction[:nbins - s_bin] * ones_array[:, s_bin:nbins] 
    elif apply_geo_corr == 0:
        #r-squared correction
        geocorr = (rs_cal.geo.data[:(nbins-s_bin),0]**2
                     * ones_array[:,s_bin:nbins])/1e6
    else:
        #no correction applied
        return rs
    
    if hasattr(rs,'telescope_pointing') and hasattr(rs_cal,'ngeo') \
               and (rs.telescope_pointing<0.1).any():
        mask = (rs.telescope_pointing < 0.1)
        indices = arange(ones_array.shape[0])
        indices= indices[mask]
        geocorr[indices,:] = rs_cal.ngeo.data[:nalts - s_bin, 1]\
                    * ones_array[indices, s_bin:nalts]
    for field in apply_to:    
       if hasattr(rs,field):
          temp = getattr(rs, field)
          #print 'field_________________________',field
          if 0:  #field == 'nitrogen_counts_low':
             import matplotlib.pylab as plt
             bin_vec = np.arange(temp.shape[1])
             plt.figure(77777)
             plt.plot(bin_vec,temp[0,:],'b', s_bin + np.arange(nbins-s_bin),geocorr[0,:],'g')
          temp[:,s_bin:nbins]*=  geocorr
          if 0:  #field == 'nitrogen_counts_low':
             plt.plot(bin_vec,np.nanmean(temp,0),'r')
                 
          setattr(rs,field,temp)
       
      
    #add geo_corr to rs and offset by s_bin to align with range scale of data buffers
    geo_corr = np.zeros(ones_array.shape[1])
    geo_corr[s_bin:nbins] = rs_cal.geo.data[:nbins-s_bin,1].copy()
    rs.geo_corr = hau.Z_Array(geo_corr)

    
    #make full resolution geo_corr array with sbin offset to match data arrays
    rs.geo_corr_array = hau.TZ_Array(np.zeros(ones_array.shape))  
    rs.geo_corr_array[:,s_bin:nbins] = geocorr
    return rs



def clear_first_bins(select,rs,process_control,consts):
    first_bin = process_control.get_value('first_bin_to_process','bin_number')
    first_bin = first_bin + consts['data_bin_containing_laser_pulse']
    for name in select:
         var = getattr(rs,name)
         var[:,:first_bin] = np.NaN
         setattr(rs,name,var)
    return


def dark_count_correction_prefire(chan_sel_dict,rs,corr_adjusts=None):
    """
       dark_count_correction_prefire(chan_sel_dict,rs,corr_adjusts=None)
       uses dark count extracted from range bins prior to laser pulse 
       for each profile and subtracts this value for each profile

       chan_sel_dict = dict(chanel_names = "dark_count_names" ......) 
                       correct channel_names with corresponding dark_counts
       rs       = structure containing channel counts to be corrected
       corr_adjusts = optional dark count scaling used for testing
       
    """
 
    for channel_name,dark_count_name in chan_sel_dict.iteritems():
       if hasattr(rs,channel_name): 
          channel_counts = getattr(rs, channel_name)
          if not corr_adjusts == None and corr_adjusts.has_key(dark_count_name):
              dark_counts = getattr(rs,dark_count_name)*corr_adjusts[dark_count_name]
          else:
               dark_counts = getattr(rs,dark_count_name)
          ones_array = np.transpose(np.ones_like(channel_counts))
          dark_corr_array = np.transpose(ones_array * dark_counts)
          channel_counts -= dark_corr_array 
          setattr(rs,channel_name,channel_counts)
    return rs

def dark_count_correction_from_signal(chan_sel_dict,rs,start_index,end_index):
    """
       dark_count_correction_from_signal(apply_to,rs,start_index,end_index)
       finds nanmean of signal in altitude bins between start and end indices 
       for each profile and subtracts this value for each profile
       Also add a vector T_Array vector of dark counts for each profile.
       
       chan_sel_dict = dict(chanel_names = "dark_count_names" ......) 
                       correct channel_names with corresponding dark_counts
       rs       = structure containing variables to be corrected
       start_index = array index for start of dark average
       end_index = array index for end of dark count average
    """
 
    #for field in apply_to:
    for channel_name,dark_count_name in chan_sel_dict.iteritems():    
       #if hasattr(rs,field):
       if hasattr(rs,channel_name):
          temp_d = getattr(rs, channel_name).copy()
          ones_array = np.transpose(np.ones_like(temp_d))
          dark_corr = hau.T_Array(np.nanmean(temp_d[:,start_index:end_index],1))
          dark_corr_array = np.transpose(ones_array * dark_corr)
          temp_d -= dark_corr_array
          #write corrected counts back to rs
          setattr(rs,channel_name,temp_d)
          #print 'field entering',channel_name
                
          #add dark count field to rs
          setattr(rs,dark_count_name,dark_corr)
          
    return rs

def Rayleigh_cross_section(wavelength,pressures,temps,altitudes):
    """
       Rayleigh_lidar_return(wavelength,,pressures,temperatures,altitudes):
       compute Rayleigh scattering cross section
       see R Holz thesis for this equation giving the Rayleigh scattering cross section.
       beta=3.78e-6*press/temp at a wavelength of 532 nm
       then rescale to actual wavelength
    """
    
    nalts = altitudes.shape[0]
    beta_r = hau.Z_Array(np.zeros( nalts ))
    
    #Rayleigh scatteromg cross section at 532 nm
    beta_r[:nalts] = 3.78e-6 * pressures[:nalts]  /temps[:nalts]

    #Rayleigh scattering  cross section
    beta_r = hau.Z_Array(beta_r * (532.0 / wavelength)**4)
    
    return beta_r

def integrated_backscatter(beta_a_backscat,altitudes):
    #compute integrated aerosol backscatter  
    temp = beta_a_backscat.copy()
    temp[np.isnan(temp)]=0.0
    dz = altitudes[1]-altitudes[0]
    return np.cumsum(temp,1)*dz   
    


def lidar_return(beta_r,beta_e_trans,beta_e_rec,altitudes):
    """
       lidar_return(beta_r,beta_e_trans,beta_e_rec,altitudes
       beta_r       = backscatter cross section profile(1/(m sr)
       beta_e_trans = extinction cross section at transmit wavelength (1/m)
       beta_e_rec   = extinction cross section at receive wavelength (1/m)
       altitudes    = vector of altitudes at which to make calculations (m)
       atten_bs     = attenuated backscatter cross_section_profile (1/m)
       
    """
    
    delta_r = np.zeros_like(altitudes)
    delta_r[:-1] = altitudes[1:]-altitudes[:-1]
    
    #optical depths on tranmit and receive paths
    od_trans = np.cumsum(beta_e_trans*delta_r)
    od_rec   = np.cumsum(beta_e_rec*delta_r)
    atten_bs = hau.Z_Array(beta_r * exp(- od_trans -od_rec))

    return atten_bs

def Raman_n2_backscatter_cross_section(beta_r,wavelength_elastic,wavelength_n2_raman,temps):
    """
       Raman_nitrogen_cross_section(wavelength,presures,temperatures,altitudes)
       computed from the ratios of Raman to Rayleigh scattering using
       equations 9.5 and 9.6 of Wandinger page 247 'lidar--Range Resolved Optical...
    """
    #wavenumber (cm^-1)elastic scattered light, transmitted light
    k_elastic = 1.0 / (wavelength_elastic * 1e-9)

    #wavenumber (cm^-1) for raman shifted light
    k_n2   = 1.0 / (wavelength_n2_raman * 1e-9)

    c         = 2.997e8              #speed of light m/s
    h         = 6.626e-34          #planck's constant m^2 kg / s
    kb        = 1.38e-23           #Boltzman constant m^2  kg /(s^2 K-deg)
    e0        = 8.854e-12          #permitivity of free space F/m

    e_const_sq = (4 * np.pi * e0)**2

    #kv = np.pi**2 / e0**2
    kv = (2 * np.pi)**4
    #Rayleigh n2  
    a_n2_sq     = 3.17e-60 #* e_const_sq
    gamma_n2_sq = 0.52e-60 #* e_const_sq
    
 
    #Rayleigh o2
    a_o2_sq     = 2.66e-60 #* e_const_sq
    gamma_o2_sq = 1.26e-60 #* e_const_sq

    #Raman_n2
    a_prime_sq     =  2.62e-14 #* e_const_sq
    gamma_prime_sq =  4.23e-14 #* e_const_sq


    #differential cross sections divided by k_nu
    Rayleigh_n2 = kv * k_elastic**4 * (a_n2_sq + 7.0 / 180.0 * gamma_n2_sq)
    Rayleigh_o2 = kv * k_elastic**4 * (a_o2_sq + 7.0 / 180.0 * gamma_o2_sq)
    Rayleigh_air = 0.8 * Rayleigh_n2 +0.2 * Rayleigh_o2   
   
    b_nu_sq = h /(8 * np.pi**2 * c * k_n2)
    Raman_n2 = kv * (k_elastic - k_n2)**4 * b_nu_sq /(1- np.exp(-h * c * k_n2 / (kb*temps)))\
                     *(a_prime_sq +(7.0/45.0) * gamma_prime_sq)                             
   

    ratios = Raman_n2 / Rayleigh_air

    beta_raman = ratios * beta_r
    return beta_raman,Rayleigh_air


def compute_optical_depth( Nm, beta_r_backscat, msl_altitudes, processing_defaults
                           ,constants,telescope_pointing = None):
    """uses Nm and beta_r_backscat to compute the optical
       depth profile with the optical depth profile set to zero at
       an alititude given by 'mol_norm_alt' which must provided
       in meters.
       returns:
           od = total optical depth
           od_aerosol = aerosol optical depth
           mol_norm_index = bin number at which optical depth is normalized to zero
           mol_ref_aod = estimated optical depth at altitudes[mol_norm_index]"""

  
    #mol_od_between_lidar_and_norm_alt = 0.0 
    mol_norm_alt=processing_defaults.get_value('mol_norm_alt','meters')
    mol_norm_index = len(msl_altitudes[msl_altitudes <= mol_norm_alt])
    if processing_defaults.enabled('molecular_smooth'):
       window_offset = np.float(processing_defaults.get_value('molecular_smooth'
                    ,'window_length'))/2.0
    else:
       window_offset = 0.0
    if ('installation' not in constants or constants['installation'] == 'ground' 
               or constants['installation'] == 'shipborne')\
               and mol_norm_alt < (constants['lidar_altitude']+150.0 + window_offset):
        mol_norm_alt = constants['lidar_altitude'] + 150. +window_offset
        mol_norm_index = len(msl_altitudes[msl_altitudes <= mol_norm_alt])
        lidar_alt_index = len(msl_altitudes[msl_altitudes <= constants['lidar_altitude']])
        
        print 'requested molecular norm altitude too low--  reset at ', mol_norm_alt,' m'
        processing_defaults.set_value('mol_norm_alt','meters',mol_norm_alt)              
    od = hau.TZ_Array(np.NaN * np.zeros_like(Nm))
    indices = np.arange(len(msl_altitudes))
   
    if processing_defaults.enabled('molecular_smooth'):
        window = processing_defaults.get_value('molecular_smooth','window_length')
        if mol_norm_alt < constants['lidar_altitude'] + 150 +window/2.0 :
           mol_norm_alt = constants['lidar_altitude'] + 150 +window/2.0
           mol_norm_index =len(msl_altitudes[msl_altitudes <= mol_norm_alt])
           print ' '
           print 'Warning:******requested mol_normalization_altitude not within acquired data'
           print '      setting normalization altitude = ',msl_altitudes[mol_norm_index]
           print ' '
    if mol_norm_index >= len(msl_altitudes):
           mol_norm_index= len(msl_altitudes) -1
           print ' '
           print 'Warning:******requested normalization altitude higher than requested range'
           print '     setting normalization altitude =',msl_altitudes[mol_norm_index]
       
    if mol_norm_index <= Nm.shape[1]:
        if 0:
            print
            print
            print 'alts',msl_altitudes.shape
            print 'beta_r',beta_r_backscat.shape
            print 'Nm',Nm.shape
            print 'norm_alt',msl_altitudes[mol_norm_index]
            print
            import matplotlib.pylab as plt
            plt.figure(2000)
            plt.plot(msl_altitudes,np.nanmean(Nm,0)*beta_r_backscat[mol_norm_index]/np.nanmean(Nm[:,mol_norm_index],0),msl_altitudes,beta_r_backscat)
            plt.ylabel('altitude')
            ax = plt.gca()
            ax.set_yscale('log')
            plt.grid(True)
      
        time_vec = np.ones_like(Nm[:,0])
        bin_vec =np.ones_like(beta_r_backscat)
        beta_r_array = (time_vec[:,np.newaxis] * beta_r_backscat[np.newaxis,:])     
        Nm_norm_array = Nm[:,mol_norm_index]
        Nm_norm_array = Nm_norm_array[:,np.newaxis] * bin_vec[np.newaxis,:]
        prelog_od=Nm[:, :] / (beta_r_array * Nm_norm_array) * beta_r_backscat[mol_norm_index]
        prelog_od[np.isnan(prelog_od)] = 0
        prelog_od[prelog_od<=0.0]=np.NaN
        od = -0.5 * np.log(prelog_od)
        if telescope_pointing is not None:
            tmpod=od
            od=od.copy()
            od[:,:]=np.NAN
            od[telescope_pointing<0.1,:] = -1.0 * tmpod[telescope_pointing<0.1,:]
            od[telescope_pointing>0.9,:] =        tmpod[telescope_pointing>0.9,:]
        dz = msl_altitudes[1]-msl_altitudes[0]    
        mol_od = 8.0 * np.pi * np.cumsum(beta_r_backscat)*dz/3.0
        mol_od = mol_od -mol_od[mol_norm_index]
        
        mol_od_array = time_vec[:,np.newaxis] * mol_od[np.newaxis,:] 
        od_aerosol = od - mol_od_array
    else:
        od[:,:]=np.NaN
        od_aerosol[:,:]=od.copy()
        print ' '
        print '*******requested od_normalization altitude not within acquired data'
        print ' '

    #compute optical depth below mol_norm_index
    #extrapolate from just above mol_norm_index to estimate unmeasured optical depth
    if dz >= 100:
        n_bins = 1    
    else:
        n_bins = int(100.0/dz +1)
    norm_alt_range = msl_altitudes[mol_norm_index]-constants['lidar_altitude']     
    mol_ref_aod = (od[:,mol_norm_index + n_bins] - od[:,mol_norm_index])* (norm_alt_range/(dz*n_bins))
    mol_ref_aod -= (mol_od_array[:,mol_norm_index + n_bins] - mol_od_array[:,mol_norm_index])\
                      *norm_alt_range/(dz*n_bins)
    mol_ref_aod = hau.T_Array(mol_ref_aod)
    od = hau.TZ_Array(od)
    od_aerosol = hau.TZ_Array(od_aerosol)
    return (od, od_aerosol,  mol_norm_index,mol_ref_aod)

def compute_klett_backscatter(signal,beta_r_backscat,altitudes,lidar_ratio,ref_altitude):
   """ compute_klett_backscatter(signal,beta_r_backscat,altitudes,lidar_ratio,ref_altitude)
       klett solution for backscatter cross-section
       signal         = range corrected counts profile
       beta_r_backscat= molecular backscatter cross-section at altitudes in altitudes vector
       altitudes      = bin altitudes at which signal,beta_r_backscat are defined
       ref_altitude   = altitude where aerosol scattering is absent (m) 
       LR             = assumed lidar ratio for aerosol
       
   """
   
   delta_r = altitudes[2]-altitudes[1]
   #ref_bin = np.int(processing_defaults.get_value('klett','ref_altitude') * 1000.0 / delta_r )
   #LR_1064 = processing_defaults.get_value('klett','LR_ratio_1064')
   ref_bin = np.int((ref_altitude-altitudes[0]) /delta_r)
   if ref_bin>=signal.shape[1]:
       ref_bin=-1
       print
       print 'klett ref altitude set too high, reset to ' \
             +str(altitudes[ref_bin])+ ' m'
   sig = signal.copy()
   sig[np.isnan(sig)] = 0.0
   beta_r_b = beta_r_backscat * np.ones_like(sig)

   #one dimensional integral in range
   TR = np.cumsum(beta_r_backscat * (lidar_ratio -8 * np.pi /3.0))*delta_r
   #integral between bin positions and the reference bin
   TR = TR - TR[ref_bin] 
   #make TR 2-D by adding time axis
   TR = np.ones_like(sig) *TR
   TR = np.exp(-2 * TR)
 
   sig_0 = sig[:,ref_bin] /beta_r_backscat[ref_bin]
   sig_0 = np.ones_like(altitudes) * sig_0[:,np.newaxis]
   
   int_at_R = np.ones_like(sig) * np.NaN
   int_at_R = delta_r * lidar_ratio * np.cumsum(sig * TR ,1)

   #get integral values at ref_bin as function of time 
   int_at_ref = int_at_R[:,ref_bin]
   #add an altitude axis
   int_at_ref = int_at_ref[:,np.newaxis]
   #fill all altitudes with value of integral at the ref bin
   int_at_ref = np.repeat(int_at_ref,sig.shape[1],1)
   #find value of integral evaluated between R and ref bin
   integral = int_at_R - int_at_ref
   #evaluate klett backscatter  
   total_backscat = sig * TR / (sig_0 - 2 * integral)
   backscat_klett = total_backscat - beta_r_b
      
   return backscat_klett


if __name__ == '__main__':
    Raman_nitrogen_cross_section()
