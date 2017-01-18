#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
from datetime import datetime,timedelta
import atmospheric_profiles.soundings.sounding_utilities as su
import lg_base.formats.calvals as cru
import lg_base.core.array_utils as hau
import lg_base.core.read_utilities as ru
from math import cos,acos

#try:
    #from bottleneck import nanmean,nansum
#except ImportError:
    #print
    #print 'No bottleneck.nanmean available! Falling back to SLOW scipy.stats.nanmean'
    #print
    #from scipy.stats import nanmean
    #from numpy import nansum

class Cxx(object):
   pass

def calval_info(inst):
    return cru.calvals_class(instrument=inst)

def update_calibrations(
    instrument,
    interval_start_time,
    n_range_ave,
    rs_raw,
    rs_constants,
    rs_soundings,
    rs_cal,
    need_new_quick_cal,
    sounding,
    process_defaults,
    corr_adjusts, requested_altitudes=None
    ):
    """check to see if any of the calibrations need to be updated.
       If they do, read the new calibration data. The the arrays,
       'rs_constants','sounding', and 'rs_cal' are returned in a
       dictionary"""


   
    # make empty dictonary for return of variables

    R = {}

    #flag to determine if a new quick_cal
    #need_new_quick_cal = False
   
    
    start_time_str = interval_start_time.strftime('%d-%b-%y %H:%M')

    # is it time to update system constants?

    #if interval_start_time >= rs_constants['next_cal_time']:
    #    rs_constants = cru.select_sys_constants(cal_info.calvals,
    #            interval_start_time)
    #    need_new_quick_cal = True
    #    print 'system constants reloaded at ', start_time_str
    print rs_constants['sounding_type']
    # is it time to update NOAA raob sounding?
    if rs_soundings==None:
      pass #caller gave us the right sounding... right?
    elif rs_constants['sounding_type'] == 'NOAA raob' \
             or rs_constants['sounding_type'] == 'time curtain':
        if interval_start_time >= sounding.expire_time:
            #sounding = hau.selectByTime(rs_soundings, interval_start_time)
            oldsounding=sounding
            sounding = rs_soundings.profile(interval_start_time,[],[])
            
         # if interval_start_time is later than the last sounding
         # check to see if more soundings are available
            if sounding.times==oldsounding.times:
                print 'loading new'
                rs_tail_soundings = su.sounding_archive(
                    instrument,
                    rs_constants['sounding_type'],
                    rs_constants['sounding_id'],
                    interval_start_time,requested_altitudes if requested_altitudes!=None else
                    np.arange(0,50000,7.5 * n_range_ave))
 
                # check to see if file has new soundings
                tail_soundings = rs_tail_soundings.profile(interval_start_time,[],[])
                if tail_soundings.times > sounding.times:
                    # trim any soundings repeated in tail file
                    #tail_soundings.trimTimeInterval( rs_soundings.times[-1, 0] + .01, np.inf )
                    rs_soundings.append( rs_tail_soundings )
                    sounding=tail_soundings

                    # remove soundings more than one day old
                    #rs_soundings.trimTimeInterval( interval_start_time - 1, np.inf )
                else:
                    # if no new sounding was found wait three hours before trying again
                    next_check_interval = timedelta(hours=3) 
                    sounding.expire_time = sounding.expire_time \
                               + next_check_interval
            need_new_quick_cal = True        
            print sounding.times.strftime('%d-%b-%y %H:%M'
                    ), sounding.station_id, 'sounding loaded at ', start_time_str, ' expiration at ',sounding.expire_time
    elif rs_constants['sounding_type'].find('GRIB file') >= 0:

        angle_delta = acos(cos(rs_raw.aircraft_latitude
                             - sounding.latitude)
                             * cos(rs_raw.aircraft_longitude
                             - sounding.longitude))
        distance_to_last_sounding = 6378.1 * angle_delta  # radius of earth *angle
        if interval_start_time > rs_constants['GRIB_radius'][0] / (24
                * 3600) + last_sounding_time \
            or distance_to_last_sounding > rs_constants['GRIB_radius'
                ][1]:
            sounding = su.read_grib_file(interval_start_time, 7.5
                    * n_range_ave, 50000, rs_raw.aircraft_latitude,
                    rs_raw.aircraft_longitude)
            need_new_quick_cal = True
    else:
        raise RuntimeError('Error in calvals definition of sounding_type')

    # is it time to update cal vectors? (i.e. geo_corr, diff_geo_cor,baseline,i2_scan_file)

    if interval_start_time >= rs_cal.expire_time:
        #max_r = 4000
        rs_cal.read(interval_start_time)
        print 'New cal vectors read at ', start_time_str
        need_new_quick_cal = True

    if need_new_quick_cal:
        R['rs_Cxx'] = update_Cxx(rs_constants,rs_cal,corr_adjusts,sounding,process_defaults)
        
    # this is a bit of a kludge
    R['sounding'] = sounding
    #R['rs_constants'] = rs_constants
    R['rs_cal'] = rs_cal
    #if need_new_quick_cal:
    #  R['rs_Cxx'] = rs_Cxx   
    return R


def update_Cxx(rs_constants,rs_cal,sounding,process_defaults,corr_adjusts):
    i2a_scan_corr = 1.0

    if rs_constants.has_key('i2a_scan_adjustment'):
            i2a_scan_corr *= rs_constants['i2a_scan_adjustment']
            if 'i2a_corr' in corr_adjusts:
              i2a_scan_corr *= corr_adjusts['i2a_corr']
   
    return quick_cal( rs_cal.i2scan.data
                       ,rs_cal.i2scan.Cam
                       ,rs_cal.i2scan.Cam_i2a     
                       ,sounding
                       ,rs_constants['wavelength']
                       ,process_defaults.get_value('molecular_spectrum'
                                  ,'model', key = 'process_defaults')
                       ,rs_constants['i2_scan_adjustment']
                            *corr_adjusts['i2_corr']
                       ,i2a_scan_corr)
   


def witschas_spectrum(
    temperature,
    pressure,
    wavelength,
    frequency,
    ):

 # compute Brillouin spectrum via witschas approximation
 # Applied optics 20-jan-2011, vol 50, #3, pp267
 # temperature = air temperature (deg K)
 # pressure    = air pressure (mb)
 # wavelength = in m

 # k=4*pi/(wavelength*sin(theta/2)
 # theta=180 deg,

    #print 'witschas pressures',pressure
    #print 'witschas temp',temperature
    k = 4 * np.pi / wavelength  # wavelength (m))

    Kb = 1.38044e-23  # Boltzman constant
    m_bar = 4.789e-26  # mean molecular weight of air (kg)
    eta = 1.846e-5  # shear viscosity of air (Pa/(m sec)

    nu_0 = np.sqrt(Kb * temperature / m_bar)
    x = 2 * np.pi * frequency / (np.sqrt(2) * k * nu_0)
    y = 100*pressure / (1* np.sqrt(2) * k * nu_0 * eta)  # multiplication by 100 to convert mb to Pa
    #print 'witschas y= ', y ,'P= ',pressure,'T= ',temperature
    A = 0.18526 * np.exp(-1.31255 * y) + 0.07103 * np.exp(-18.26117
            * y) + 0.74421
    
    sigma_r = 0.70813 + 0.16366 * y ** 2 + 0.19132 * y ** 3 - 0.07217 \
        * y ** 4
    
    sigma_b = 0.07845 * np.exp(-4.88663 * y) + 0.80400 \
        * np.exp(-0.15003 * y) - 0.45142
    xb = 0.80893 - 0.30208 * 0.10898 ** y

    spectrum = A / (np.sqrt(2 * np.pi) * sigma_r) * np.exp(-0.5 * (x
            / sigma_r) ** 2) + (1 - A) / (2 * np.sqrt(2 * np.pi)
            * sigma_b) * np.exp(-0.5 * ((x + xb) / sigma_b) ** 2) + (1
            - A) / (2 * np.sqrt(2 * np.pi) * sigma_b) * np.exp(-0.5
            * ((x - xb) / sigma_b) ** 2)
    return spectrum


def quick_cal( i2_scan, Cam, Cam_i2a,sounding, wavelength, method_string
               , i2_scan_corr, i2a_scan_corr):
    """quick_cal_stream(i2_scan, Cam, sounding, wavelength, method_string, i2_scan_corr)
       A function which computes hsrl calibration coefficients.
       at the altitudes specified (in meters) within 'alt_vector'
       using a precomputed iodine scan file(i2_scan_file).
       i2_scan(:,:) = input containing i2 scan info
       i2_scan(:,0) = freq (GHz)
       -------(:,1) = combined channel scan
       -------(:,2) = molecular channel scan
       -------(:,3) = theoretical i2 transmission
       -------(:,4) = measured i2/combined
       if bagohsrl with argon buffered i2 cell
       -------(:,5) = molecular i2a/combined
       -------(:,6) = molecular i2a channel scan

       Cam           = aerosol in molecular coefficent
       Cam_i2a       = aerosol in argon-buffered molecular channel (will = None when not present)
       sounding      = return structure from read_sounding_file.py contains temp profile
       wavelength    = laser wavelength (nm)
       method_string = molecular line shape ''maxwellian','tenti_s6','wirtschas'
       i2_scan_corr  = adjustment factor for i2 scan that adjusts
                       i2 molecular channel gain relative to combined channel gain
       i2a_scan_corr = adjustment factor for i2a scan that adjusts
                       i2a molecular channel gain relative to combined channel gain
       rs            = return structure containing calibration coefficents
                       calibration values are returned at altitudes
                       rs_sounding.alititudes[:]
                       
                       1            = particulate in combined channel
                       rs.Cmc[i]    = molecular in combined channel
                       rs.Cmm[i]    = molecular in molecular channel
                       rs.Cam       = particulate in molecular channel
                       rs.beta_r[i] = Rayleigh scattering cross section (1/m)"""

       
    rs = hau.Time_Z_Group()  # calibration structure

    # i2_scan=cal_vec.i2_scan
    i2_scan=i2_scan.copy()

    # if selected use theory*combined as synthetic mol scan
    print method_string
    if method_string.find('i2 theory') >= 0:
        i2_scan[:, 2] = i2_scan[:, 4] * i2_scan[:, 1]
   
    # trim i2 scan to +-4 GHz about line center

    #i2_scan = i2_scan[abs(i2_scan[:, 0]) <= 4, :]

    # rescale i2 molecular component of i2 scan if required
    if i2_scan_corr != 1.0:
        i2_scan[:, 2] = i2_scan[:, 2] * i2_scan_corr
        
    # rescale i2a molecular component of i2 scan if required    
    if i2a_scan_corr != 1.0 and i2_scan.shape[1]>6:
        i2_scan[:, 6] = i2_scan[:, 6] * i2a_scan_corr

    if 0:
       import matplotlib.pylab as plt
       plt.figure(444443)
       plt.plot(i2_scan[:,0],i2_scan[:,1:3],i2_scan[:,0],i2_scan[:,2]/i2_scan[:,1],'k')
       plt.xlabel('freq GHz')
       plt.grid(True)

    # trim i2 scan to +-4 GHz about line center
    i2_scan = i2_scan[abs(i2_scan[:, 0]) <= 4, :]

    if 0:
       
       import matplotlib.pylab as plt
       plt.figure(444444)
       plt.plot(i2_scan[:,0],i2_scan[:,1:3])
       plt.xlabel('freq GHz')
       plt.grid(True)
    
    # compute Rayleigh scattering cross section
    # see R Holz thesis for this equation giving the Rayleigh scatter cross section.
    # beta=3.78e-6*press/temp at a wavelength of 532 nm
    # then rescale to actual wavelength

    nalts = sounding.altitudes.shape[0]
    if not (nalts==sounding.pressures.shape[0] and nalts==sounding.temps.shape[0]):
      print "ERROR: SOMETHIGN BAD ABOUT SOUNDING AT TIME ",sounding.times
      return None
    rs.beta_r = hau.Z_Array(np.zeros( nalts ))
    rs.beta_r[:nalts] = 3.78e-6 * sounding.pressures[:nalts]  / sounding.temps[:nalts]
    rs.beta_r = rs.beta_r * (532.0 / wavelength) ** 4
   
    # spectral width of molecular scattering

    m_bar = 28.97 * 1.65978e-27  # average mass of an air molecule
    sigma_0 = 1 / (wavelength * 1e-9)  # number in 1/meters
    kb = 1.38044e-23  # Boltzmans constant J/(K deg)
    c = 3e8  # speed of light in m/s
    sigma = i2_scan[:, 0] * 1e9 / c  # wavenumber vector

    rs.Cmm     =  hau.Z_Array(np.zeros( nalts ))
    rs.Cmc     =  hau.Z_Array(np.zeros( nalts ))
    sample_altitudes = np.zeros(nalts)   
    #for bagohsrl with argon buffered i2 cell
    if len(i2_scan[0,:]) > 6:   
        rs.Cmm_i2a =  hau.Z_Array(np.zeros( nalts ))
        rs.Cam_i2a =  Cam_i2a
    rs.Cam = Cam 
    
    print 'RBS computed with '+method_string+  ' spectrum'
  
    spectrum_time = datetime.utcnow()        
    
    dz = sounding.altitudes[2]-sounding.altitudes[1]
    delta_i = np.int(np.ceil(300.0/dz))
    nk=int(nalts/delta_i)
    if delta_i>1 and nk<2:# if interpolation is to happen, but not enough to interpolate, force it to the edge
        nk=2
        delta_i=nalts-1
    sample_altitudes = np.zeros(nk)
    rs.msl_altitudes = sounding.altitudes.copy()

    i=0
    k=0
    while k < len(sample_altitudes):
        if not np.isfinite(sounding.temps[i]) or not np.isfinite(sounding.pressures[i]):
          i=i+delta_i
          k=k+1
          continue
        if method_string.find('maxwellian') >= 0:
            norm = m_bar * c ** 2 / (8 * sigma_0 ** 2 * kb
                    * sounding.temps[ i])
            spectrum = np.exp(-norm * sigma ** 2)
        elif method_string.find('tenti_s6') >= 0:
            from tenti_s6 import tenti_s6
            spectrum = tenti_s6(wavelength * 1e-9,sounding.temps[i],
                    sounding.pressures[ i],
                    i2_scan[:, 0] * 1e9)
            
        elif method_string.find('witschas') >= 0:
            spectrum = witschas_spectrum(sounding.temps[i],
                    sounding.pressures[ i], wavelength * 1e-9,
                    i2_scan[:, 0] * 1e9)

        spectrum = spectrum / sum(spectrum)
       
    
        sample_altitudes[k] = sounding.altitudes[i]
        rs.Cmc[ k] = sum(spectrum * i2_scan[:, 1])
        rs.Cmm[ k] = sum(spectrum * i2_scan[:, 2])
        if i2_scan.shape[1]>6:
            rs.Cmm_i2a[ k] = sum(spectrum * i2_scan[:, 6])
        i = i + delta_i
        k = k + 1

    # if Cxx computed at less than full altitude resolution 
    if delta_i >1:
       rs.Cmc = np.interp(sounding.altitudes,sample_altitudes[0:k-1]
                        ,rs.Cmc[0:k-1])
       rs.Cmm = np.interp(sounding.altitudes,sample_altitudes[0:k-1]
                        ,rs.Cmm[0:k-1])
       if hasattr(rs,'Cmm_i2a'):
           rs.Cmm_i2a = np.interp(sounding.altitudes,sample_altitudes[0:k-1]
                        ,rs.Cmm_i2a[0:k-1])
    print method_string, 'computed for ',k-1,' altitudes in '\
          , (datetime.utcnow() - spectrum_time).total_seconds(),' seconds'
    plots = 0
    if plots:
        import matplotlib.pyplot as plt
        plt.figure(600)
        plt.plot(i2_scan[:, 0], spectrum[:])
        fig = plt.grid(True)
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('Intensity')
        # ax=gca()
        # ax.set_yscale('log')

        plt.figure(601)
        plt.plot(rs.Cmm,sounding.altitudes/1000.0,'b'
                 ,rs.Cmm_i2a,sounding.altitudes/1000.0,'r')
        plt.grid(True)
        plt.xlabel('Cmm, Cmm_i2a')
        plt.ylabel('Altitude')
        plt.show()
        
    # add sounding identification to return structure
    rs.sounding_id = sounding.station_id
    rs.sounding_time = sounding.times
  
    return rs


    
# separate molecular and particulate sigals, compute calibrated profiles
def hsrl_inversion(r_msl, rs_Cxx, rs_constants,corr_adjusts,process_defaults):
    """hsrl_inversion(range_processed_returns,calibration_structure
    ,system_constants)
    
    Invert hsrl raw count data into separate molecular and particulate profiles
    and compute backcatter cross-section and depolarization from the profiles
    returned structure rs always returns:
       times               = times of records
       delta_t             = time seperation between records
       msl_altitudes       = bin altitudes in meters   
       seeded shots        = total number of seeded laser shots in profile
       beta_r_backscat     = Rayleigh scattering cross section from sounding
       Na                  = number of aerosol photon counts
       Nm                  = total number of  molecular counts(including i2a if present)
       Nm_i2               = number of molecular counts in i2 channel
       Na                  = number of particulate counts
       Ncp                 = Ncp photon counts
       linear_depol        = fractional particulate depolarization
       beta_a_backscat     = particulate aerosol backscatter cross-section
       beta_a_backscat_par = par polarization component of backscat cross-section
       beta_a_backscat_perp= perp polarization component of backscat cross-section
       integrated_backscat = cumsum of backscatter cross section in altitude
       
    if i2a channel exists the following are added to rs:
       Nm_i2a = number of molecular photons derived from i2a channel    
       Nm     =  rs.Nm_2i + rs.Nm_i2a
       Ncp_i2a
       beta_a_backscat_par_i2a/Nm
       beta_a_backscat_perp_i2a
       
    if these exist in input file they are added to rs:
       telescope_pointing
       GPS_MSL_Alt
       circular_depol
    """

    rs = hau.Time_Z_Group(like=r_msl)
    # r_msl.molecular_counts is  hau.TZ_Array (2D)
    nalts = r_msl.molecular_counts.shape[1]
    if  rs_Cxx.beta_r.shape[0] < nalts:
        print 'hsrl_inversion(): size too small on calibration arrays : rs_Cxx.beta_r = %d vs nalts = %d. padding cal with last value' % \
            (rs_Cxx.beta_r.shape[0], nalts)
        #assert(rs_Cxx.beta_r.shape[0] != nalts)
        os=rs_Cxx.beta_r.shape[0]
        for k,v in vars(rs_Cxx).items():
          if hasattr(v,'size') and v.size==os:
            ns=list(v.shape)
            ns[0]=nalts
            tmp=np.zeros(ns,dtype=v.dtype)
            tmp[:]=v[-1]
            tmp[:os]=v
            setattr(rs_Cxx,k,tmp)
    elif rs_Cxx.beta_r.shape[0] > nalts:
        print 'WARNING hsrl_inversion(): size larger on calibration arrays. may be an error : rs_Cxx.beta_r = %d vs nalts = %d' % \
            (rs_Cxx.beta_r.shape[0], nalts)
    rs.times = r_msl.times.copy()
    rs.delta_t = r_msl.delta_t.copy()
    rs.msl_altitudes = r_msl.msl_altitudes.copy()
    rs.seeded_shots=r_msl.seeded_shots.copy()
    if hasattr(r_msl,'telescope_pointing'):
        rs.telescope_pointing=r_msl.telescope_pointing.copy()

    # Rayleigh backscatter cross section profile
    rs.beta_r_backscat =  hau.Z_Array(np.zeros(nalts))
    rs.beta_r_backscat[:nalts] = rs_Cxx.beta_r[:nalts] * 3.0 / (8.0 * np.pi)
    #for normal i2 channel
    kk = 1.0 / (rs_Cxx.Cmm[:nalts] - rs_Cxx.Cmc[:nalts] * rs_Cxx.Cam)
    rs.Na =  hau.TZ_Array( kk[:nalts] * (rs_Cxx.Cmm[:nalts] * r_msl.combined_counts\
               - rs_Cxx.Cmc[ : nalts] * r_msl.molecular_counts) )
    rs.Nm =  hau.TZ_Array( kk[:nalts] * (r_msl.molecular_counts \
               - rs_Cxx.Cam * r_msl.combined_counts) )
    rs.Nm_i2 = rs.Nm.copy()
    
    #if data includes an i2a channel we generation a seperate inversion--bagohsrl only
    #systems with i2a channel record linear depolarization
    if hasattr(r_msl,'molecular_i2a_counts') and hasattr(rs_Cxx,'Cmm_i2a'):
       
            kk = 1.0 / (rs_Cxx.Cmm_i2a[:nalts] - rs_Cxx.Cmc[:nalts] * rs_Cxx.Cam_i2a)
            rs.Na_i2a = hau.TZ_Array( kk[:nalts] * (rs_Cxx.Cmm_i2a[:nalts]
                 * r_msl.combined_counts - rs_Cxx.Cmc[ : nalts] * r_msl.molecular_i2a_counts) )
            rs.Nm_i2a =  hau.TZ_Array( kk[:nalts] * (r_msl.molecular_i2a_counts \
                    - rs_Cxx.Cam_i2a * r_msl.combined_counts) )
            rs.Nm = (rs.Nm_i2 + rs.Nm_i2a)/2.0

            #bagohsrl is only hsrl with i2a channel--it measures linear depolarization
            rs.Ncp = rs_constants['combined_to_cross_pol_gain_ratio'] \
                    * (r_msl.cross_pol_counts
                    - rs_constants['polarization_cross_talk']*corr_adjusts['pol_x_talk']
                    * r_msl.combined_counts) - rs.Nm_i2 * 0.0035 / (1.0 - 0.0035)
            rs.Ncp_i2a = rs_constants['combined_to_cross_pol_gain_ratio'] \
                    * (r_msl.cross_pol_counts
                    - rs_constants['polarization_cross_talk']*corr_adjusts['pol_x_talk']
                    * r_msl.combined_counts) - rs.Nm_i2a * 0.0035 / (1.0 - 0.0035)
           
            if not rs_constants.has_key('no_depol_channel') or rs_constants['no_depol_channel']==0 :   
                rs.linear_depol_i2a = rs.Ncp_i2a / rs.Na_i2a
                #compute the linear depolarization as the average of normal and i2a values
                #note these two componets should be almost identical
                rs.linear_depol = (rs.Ncp + rs.Ncp_i2a)/(rs.Na_i2a + rs.Na)
            
                #when backscatter is small linear_depol can become indeterminate--bound values
                rs.linear_depol[rs.linear_depol < 0.0] = 0.0
                rs.linear_depol[rs.linear_depol > 0.6] = 0.6
            else:
                rs.linear_depol = np.zeros_like(rs.Nm)
                rs.linear_depol_i2a = np.zeros_like(rs.Nm)
                
            rs.beta_a_backscat_perp_i2a=rs.Na_i2a/rs.Nm_i2a \
                            *rs.linear_depol_i2a*rs.beta_r_backscat
            rs.beta_a_backscat_perp = rs.Na/rs.Nm_i2 \
                            *rs.linear_depol * rs.beta_r_backscat
            
            rs.beta_a_backscat_par_i2a = rs.Na_i2a / rs.Nm_i2a * rs.beta_r_backscat
            rs.beta_a_backscat_par = rs.Na / rs.Nm_i2 * rs.beta_r_backscat

            rs.beta_a_backscat = 0.5 *(rs.beta_a_backscat_par +rs.beta_a_backscat_perp \
                       + rs.beta_a_backscat_par_i2a + rs.beta_a_backscat_perp_i2a)

            if rs_constants.has_key('no_i2_channel') and rs_constants['no_i2_channel']==1:
                print
                print 'WARNING--I2 channel is not being used in calculations'
                print "calvals has 'no_i2_channel'== 1"
                print
                rs.linear_depol = rs.Ncp_i2a / rs.Na_i2a 
                rs.beta_a_backscat = rs.beta_a_backscat_par_i2a + rs.beta_a_backscat_perp_i2a 
            
            if not process_defaults.enabled('depolarization_is_aerosol_only'):    
                #user wants bulk depolarization aerosol combined with molecular
                #recompute depolarization
                print 'computing bulk depolarization--aerosol and molecular combined'
                rs.linear_depol = rs_constants['combined_to_cross_pol_gain_ratio'] \
                    * (r_msl.cross_pol_counts\
                    - rs_constants['polarization_cross_talk']*corr_adjusts['pol_x_talk']\
                    * r_msl.combined_counts) / r_msl.combined_counts
            #compute circular polarization from linear--good only for vertical pointing systems.    
            rs.circular_depol = 2 * rs.linear_depol / (1 - rs.linear_depol)
            
    elif hasattr(rs_Cxx,'Cmm_i2a') or rs_constants.has_key('polarization_is_linear') \
                and rs_constants['polarization_is_linear']==1:
       if hasattr(rs_Cxx,'Cmm_i2a'):
            print
            print 'hsrl_inversion(): WARNING  i2a counts found, but no calibration'
            print 'computing without i2a channel'
            print
       rs.Ncp = rs_constants['combined_to_cross_pol_gain_ratio'] \
               * (r_msl.cross_pol_counts
               - rs_constants['polarization_cross_talk']*corr_adjusts['pol_x_talk']
               * r_msl.combined_counts) - rs.Nm_i2 * 0.0035 / (1.0 - 0.0035)
       rs.linear_depol = rs.Ncp/rs.Na
            
            #when backscatter is small linear_depol can become indeterminate--bound values
       rs.linear_depol[rs.linear_depol < 0.0] = 0.0
       rs.linear_depol[rs.linear_depol > 0.6] = 0.6

       if not process_defaults.enabled('depolarization_is_aerosol_only'):    
                #user wants bulk depolarization aerosol combined with molecular
                #recompute depolarization 
                rs.linear_depol = rs_constants['combined_to_cross_pol_gain_ratio'] \
                    * (r_msl.cross_pol_counts\
                    - rs_constants['polarization_cross_talk']*corr_adjusts['pol_x_talk']\
                    * r_msl.combined_counts) / r_msl.combined_counts


       #compute circular from linear--good only for vertical pointing system
       rs.circular_depol = 2*rs.linear_depol /(1 - rs.linear_depol)
       
       rs.beta_a_backscat_perp = rs.Na/rs.Nm_i2 \
                            *rs.linear_depol * rs.beta_r_backscat
       rs.beta_a_backscat_par = rs.Na / rs.Nm * rs.beta_r_backscat
       rs.beta_a_backscat = rs.beta_a_backscat_par +rs.beta_a_backscat_perp
    else: #instrument with no i2a channel and measures circular polarization
        rs.Ncp = rs_constants['combined_to_cross_pol_gain_ratio'] \
                * (r_msl.cross_pol_counts\
                - rs_constants['polarization_cross_talk']*corr_adjusts['pol_x_talk']\
                * r_msl.combined_counts) - rs.Nm_i2 * 0.0074 / (1.0 - 0.0074)
        rs.circular_depol = rs.Ncp / rs.Na

        #when Na becomes small, circular_depol may become indeterminate
        # (and don't fault on Nans)
        #rs.circular_depol = np.nan_to_num(rs.circular_depol)
        rs.circular_depol[rs.circular_depol < 0.0] = 0.0
        rs.circular_depol[rs.circular_depol > 3.0] = 3.0

        rs.beta_a_backscat_perp=rs.Na/rs.Nm_i2*rs.circular_depol*rs.beta_r_backscat
        rs.beta_a_backscat_par = rs.Na / rs.Nm_i2 * rs.beta_r_backscat
        rs.beta_a_backscat = rs.beta_a_backscat_par + rs.beta_a_backscat_perp
    
        if not process_defaults.enabled('depolarization_is_aerosol_only'): 
           #user wants bulk depolarization containing both aerosol and molecular
           print 'computing bulk depolarization--air and aerosol together'
           rs.circular_depol = rs_constants['combined_to_cross_pol_gain_ratio'] \
                    * (r_msl.cross_pol_counts\
                    - rs_constants['polarization_cross_talk']*corr_adjusts['pol_x_talk']\
                    * r_msl.combined_counts) / r_msl.combined_counts
        rs.linear_depol = rs.circular_depol / (2
                    + rs.circular_depol)
    #compute integrated backscatter cross section
    rs.integrated_backscat = rs.beta_a_backscat.copy()
    rs.integrated_backscat[np.isnan(rs.integrated_backscat)] = 0.0
    da=rs.msl_altitudes.copy()
    da[:-1]=da[1:]-da[:-1]
    da[-1]=da[-2]
    tda=hau.TZ_Array(np.transpose(da[:,np.newaxis]*np.ones((1,rs.times.size))))
    rs.integrated_backscat = np.cumsum(rs.integrated_backscat,1) \
                 *tda               

    if hasattr(r_msl,'GPS_MSL_Alt'):
        rs.GPS_MSL_Alt = r_msl.GPS_MSL_Alt

    if hasattr(r_msl,'combined_1064_counts'):
        rs.combined_1064_counts = r_msl.combined_1064_counts\
                              * rs_constants['IR_combined_hi_gain_ratio']

    if hasattr(r_msl,'wfov_extinction_corr'):
        #transf wfov_extinction correction to output structure
        rs.wfov_extinction_corr = r_msl.wfov_extinction_corr
  
    return rs


def compare_calfiles(instrument,calfile,rs_cal):    
    """Reads calibration file, 'calfile' and plots it with the default
       calibration in effect at a given time, 'time'. The calibration file
       header and data are returned.
       """
   
    import matplotlib.pylab as plt
    plt.figure(1000)
    plt.figure(1000).canvas.set_window_title(calfile)
    plt.rcParams['figure.figsize'] = [8,6]
    plt.subplots_adjust(
            top=.9,
            bottom=.1,
            left=.1,
            right=.9,
            hspace=0,
            wspace=0,
            )

    [header,data]=ru.readascii(calfile);
    print ' '
    print 'compare "',calfile, '" with active default calfile'
    print header
    print ' '
     
    # is requested file a baseline file
    if calfile.find('.blc')>0:
        #read new baseline file and make comparison plots
        n_blc_header=header
        n_blc=data
        energy_index=n_blc_header.find('#ave energy per shot=')
        print 'energy per shot = ',n_blc_header[energy_index+ 22:energy_index + 30]
        n_blc[:, 1:] = n_blc[:, 1:] / float(n_blc_header[energy_index
                + 22:energy_index + 30])

        #don't plot beyound the end of the max altitude in the current request
        n=np.sum(n_blc[:,0]<rs_cal.baseline.data[-1,0])
        
        lines=plt.plot(n_blc[0:n,0],n_blc[0:n,1],'r',n_blc[0:n,0],n_blc[0:n,2],'c'\
                      ,n_blc[0:n,0],n_blc[0:n,3],'b',n_blc[0:n,0],n_blc[0:n,4],'g')
        plt.plot(rs_cal.baseline.data[:,0],rs_cal.baseline.data[:,1],'r' \
                 ,rs_cal.baseline.data[:,0],rs_cal.baseline.data[:,2],'c'\
                 ,rs_cal.baseline.data[:,0],rs_cal.baseline.data[:,3],'b' \
                 ,rs_cal.baseline.data[:,0],rs_cal.baseline.data[:,4],'g')
        
        plt.grid(True)
        ax = plt.gca()
        plt.setp(lines[0:4], linewidth=3)
        #plt.setp(lines[2], linewidth=3)

        ax.set_yscale('log')
        plt.axis([0, n , 1e-5, 1])
        plt.xlabel('bin number')
        plt.ylabel('Baseline correction')
        plt.legend(('comb hi','comb lo','mol','c pol'), 'upper right')
        plt.title('Energy normalized baseline correction (bold=newfile)')


        
        #plot fractional change    
        plt.figure(1001)
        plt.figure(1001).canvas.set_window_title(calfile)
        plt.rcParams['figure.figsize'] = [8,6]
        plt.subplots_adjust(
            top=.9,
            bottom=.1,
            left=.1,
            right=.9,
            hspace=0,
            wspace=0,
            )
        def_comb_hi=np.interp(n_blc[:,0],rs_cal.baseline.data[:,0],rs_cal.baseline.data[:,1])
        def_comb_lo=np.interp(n_blc[:,0],rs_cal.baseline.data[:,0],rs_cal.baseline.data[:,2])
        def_mol=np.interp(n_blc[:,0],rs_cal.baseline.data[:,0],rs_cal.baseline.data[:,3])
        def_cpol=np.interp(n_blc[:,0],rs_cal.baseline.data[:,0],rs_cal.baseline.data[:,4])
        
       
    
        lines=plt.plot(n_blc[:n,0],n_blc[:n,1]/def_comb_hi[:n],'r' \
                       ,n_blc[:n,0],n_blc[:n,2]/def_comb_lo[:n],'c'\
                       ,n_blc[:n,0],n_blc[:n,3]/def_mol[:n],'b' \
                       ,n_blc[:n,0],n_blc[:n,4]/def_cpol[:n],'g')
        
        plt.grid(True)
        plt.setp(lines[:], linewidth=3)
        ax=plt.gca()
        plt.axis([0, n, 0, 5])
        plt.xlabel('Bin number')
        plt.ylabel('new baseline / default')
        plt.legend(('comb hi','comb lo','mol','c pol'),'upper right')
        plt.title('Ratio new baseline to default baseline')

        #return the new baseline data
        return header, n_blc
        
    # is requested file an i2 scan file
    elif calfile.find('-scan-')>0 and calfile.find('.cal')>0:
        #n_i2_header=header
        n_i2_scan=data
        [n,k]=rs_cal.i2scan.data.shape
        print ' '
        print 'i2_scan header '
        print header
        print ' '
        plt.figure(1000)
        lines=plt.plot(rs_cal.i2scan.data[:,0],rs_cal.i2scan.data[:,1],'k' \
            ,n_i2_scan[:,0],n_i2_scan[:,1],'r' \
            ,rs_cal.i2scan.data[:,0],rs_cal.i2scan.data[:,2],'k' \
            ,n_i2_scan[:,0],n_i2_scan[:,2],'r')
                       
        plt.grid(True)
        ax = plt.gca()
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('Intensity')
        plt.legend(('chi/mol','default'), 'upper right')
        plt.title('I2 scan comparison')              

        #rato of channels plt.figure(1001)
        lines=plt.plot(n_i2_scan[:,0],n_i2_scan[:,2]/n_i2_scan[:,1],'b'
        ,rs_cal.i2scan.data[:,0],rs_cal.i2scan.data[:,2]/rs_cal.i2scan.data[:,1],'r'
        ,rs_cal.i2scan.data[:,0],rs_cal.i2scan.data[:,4],'k')
        plt.grid(True)
        ax = plt.gca()
        plt.setp(lines[0], linewidth=2)
        plt.axis([-15, 15 , 0,5])
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('I2/combined ratio')
        plt.legend(('new','default','theory'), 'upper right')
        plt.title('I2/combined ratio vs frequency')

        #ratio new to default 
        plt.figure(1002)
        new_ratio=np.interp(rs_cal.i2scan.data[:,0],n_i2_scan[:,0],n_i2_scan[:,2]/ \
                            n_i2_scan[:,1])
        lines=plt.plot(rs_cal.i2scan.data[:,0],new_ratio/  \
                     (rs_cal.i2scan.data[:,2]/rs_cal.i2scan.data[:,1]),'b')

        plt.grid(True)
        ax = plt.gca()
        #plt.setp(lines[0], linewidth=2)
        plt.axis([-15, 15 , 0 ,5])
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('new/default transmission')
        #plt.legend(('new/default'), 'upper right')
        plt.title('new/default transmission vs frequency')

        #return new i2 scan
        return header, n_i2_scan
        
    #is requested file a diff geofile
    elif calfile.find('.geo') > 0 and calfile.find('diff') >0 \
             and calfile.find('i2a') < 0:
           #n_d_geo_header=header
           n_d_geo=data
           n=len(data[:,0])
           plt.figure(2000)
           #note add one to default in order to match raw files
           lines=plt.plot(n_d_geo[0:n,0]/1000.0,n_d_geo[0:n,1],'b'\
                          ,rs_cal.diff_geo.data[:,0]/1000.0,rs_cal.diff_geo.data[:,1]+1,'r') 
           plt.grid(True)
           ax = plt.gca()
           plt.setp(lines[0:1], linewidth=3)
           plt.axis([0, rs_cal.diff_geo.data[-1,0]/1000.0 , 0.8,1.2])
           plt.xlabel('Altitude (km)')
           plt.ylabel('Comb hi diff geo correction')
           plt.legend(('chi/mol','default'), 'upper right')
           plt.title('Combined hi differential geo correction')

           plt.figure(2001)
           #note add one to default in order to match raw files
           lines=plt.plot(n_d_geo[0:n,0]/1000.0,n_d_geo[0:n,2],'b'\
                   ,rs_cal.diff_geo.data[:,0]/1000.0,rs_cal.diff_geo.data[:,2]+1,'r') 
           plt.grid(True)
           ax = plt.gca()
           plt.setp(lines[0:1], linewidth=3)
           plt.axis([0, rs_cal.diff_geo.data[-1,0]/1000.0 , 0.8,1.2])
           plt.xlabel('Altitude (km)')
           plt.ylabel('Diff geo correction')
           plt.legend(('clo/mol','default'), 'upper right')
           plt.title('Combined lo differential geo correction')
           
           #return new diff geo file
           data=n_d_geo
           data[0:n,1]=data[0:n,1]-1
           return header,data
       
    # if it is either a zenith pointing or nadir pointing geofile    
    elif calfile.find('.geo') > 0  and calfile.find('diff') <0:
           #new_geo_header=header
           new_geo=data
           #remove r-squared correction
           
           if not instrument == 'gvhsrl' or rs_cal.n_geo.filename == None :
               f_type = 'zenith'
               old_geo=np.ones_like(rs_cal.geo.data)
               old_geo[:,1]=1e6*rs_cal.geo.data[:,1]/(rs_cal.geo.data[:,0]*rs_cal.geo.data[:,0])
               old_geo[:,0]=rs_cal.geo.data[:,0]
           else: #for gvhsrl with nadir geofile present ask telescope pointing direction
               print '    gvhsrl may have separate geofiles for nadir and zenith pointing'
               f_type=raw_input('    pointing direction (nadir, zenith) ?  ')
               if f_type == 'nadir':
                  old_geo=np.ones_like(rs_cal.n_geo.data)
                  old_geo[:,1]=1e6*rs_cal.n_geo.data[:,1]/(rs_cal.n_geo.data[:,0] \
                            *rs_cal.n_geo.data[:,0])
                  old_geo[:,0]=rs_cal.n_geo.data[:,0]
               else:
                  old_geo=np.ones_like(rs_cal.geo.data)
                  old_geo[:,1]=1e6*rs_cal.geo.data[:,1]/(rs_cal.geo.data[:,0]*rs_cal.geo.data[:,0])
                  old_geo[:,0]=rs_cal.geo.data[:,0]
           
           
          

           plt.figure(1000)
           lines=plt.plot(new_geo[:,0]/1000.0 ,new_geo[:,1],'b'\
                    ,old_geo[:,0]/1000.0 ,old_geo[:,1],'r')       
           plt.grid(True)
           ax = plt.gca()
           plt.setp(lines[0], linewidth=3)
           ax.set_yscale('log')
           plt.axis([0, rs_cal.geo.data[-1,0]/1000 , 0.5, 100])
           plt.xlabel('Alitude (km)')
           plt.ylabel('Geo correction')
           plt.legend(('new','default'), 'upper right')

           n=calfile.rfind('/')+1
           plt.title(calfile[n:]+' vs. '+f_type+' default')

           #plot fractional change
           plt.figure(1001)
           temp_geo=np.interp(old_geo[:,0],new_geo[:,0],new_geo[:,1])
           
           
           lines=plt.plot(old_geo[:,0]/1000.0,temp_geo/old_geo[:,1],'b')       
           plt.grid(True)
           ax = plt.gca()
           plt.setp(lines[0], linewidth=3)
           plt.xlabel('Altitude (km)')
           plt.ylabel('new geo corr / default')
           plt.title(calfile[n:] + ' to '+f_type+' default geo corr')

           #add the r_squared correction to the new data
           data[:,1]=1e-6*new_geo[:,1]*(new_geo[:,0]*new_geo[:,0]) 
           return header,data
        
    elif calfile.find('i2a_mol_diff_geofile')>=0:
        print 'comparing '+ calfile + ' with default i2a_mol_diff_geofile'
        plt.figure(1001)
        plt.plot(data[:,1],data[:,0],'b'
                ,rs_cal.i2a_diff_geo.data[:,1]+1
                ,rs_cal.i2a_diff_geo.data[:,0],'r')
        plt.grid(True)
        plt.legend(['new','old'])
        ax=plt.gca()
        plt.axis([.95, 1.05, 0, 4000])
        plt.title('i2a geofile comparison')
        data[:,1] = data[:,1]-1
        return header,data

    elif calfile.find('i2a_temp_table_')>=0:
       print 'comparing '+ calfile + ' with default i2a_temp_table'
       return header,data
    
    #invalid file name
    else:
           print ' '
           print '**** "', calfile, '" was not found ********'
    return

def compute_raw_color_ratio(rs,Cxx,rs_constants,corr_adjusts):
    """compute_color_ratio(rs,Cxx,rs_constants):
       compute the ratio of 1064/532 counts after range_process
       corrections have been applied. This includes both aerosol
       and molecular contributions. Ratio is corrected for the
       difference in molecular attenuation at the two wavelengths"""
   
    mol_od = np.zeros_like(Cxx.msl_altitudes)
    delta_r = rs_constants['binwidth']*1.5e8 #binwidth in seconds times 1/2 speed of light
    mol_od[Cxx.msl_altitudes >0]  =  np.cumsum(Cxx.beta_r[Cxx.msl_altitudes >0])*delta_r
    shot_energy_532 = rs.transmitted_energy[:,np.newaxis].copy()
    if hasattr(rs,'transmitted_1064_energy'):
        shot_energy_1064 = rs.transmitted_1064_energy[:,np.newaxis].copy()
    else:
        shot_energy_1064 = shot_energy_532
    #add cross pol to 532 combined for comparison with unpolarized 1064    
    total_combined_counts = rs.combined_counts \
           + rs_constants['combined_to_cross_pol_gain_ratio'] * rs.cross_pol_counts
    raw_color_ratio = rs.combined_1064_counts  \
           /(total_combined_counts * rs_constants['IR_combined_hi_gain_ratio'] * (corr_adjusts['1064_gain'] if '1064_gain' in corr_adjusts else 1.0))\
           * shot_energy_532 / shot_energy_1064 
    raw_color_ratio = raw_color_ratio * np.exp(-30.0 * mol_od / 16.0)
    if 0:
       import matplotlib.pylab as plt
       plt.figure(222)
       plt.plot(rs.combined_counts[0,:],Cxx.msl_altitudes,'r'
                ,rs.cross_pol_counts[0,:],Cxx.msl_altitudes,'g'
                ,rs.combined_1064_counts[0,:],Cxx.msl_altitudes,'k')
       plt.grid(True)
       ax=plt.gca()
       ax.set_xscale('log')

    return raw_color_ratio

def compute_color_ratio(inv):
    """compute_color_ratio(inv):
       compute the ratio of 1064/532 backscatter"""
 
    inv.color_ratio = inv.beta_a_1064_backscat / inv.beta_a_backscat

    return

   
def compute_1064_aerosol_backscatter(rs,inv,processing_defaults,constants,corr_adjusts):
    """compute_1064_aerosol_backscatter(rs,inv,processing_defaults,constants)
    computes 1064 aerosol backscatter from 1064 channel and inverted backscatter
    cross sections using an assumed angstrom coefficent to estimate the 1064 aerosol
    extinction from the measured 532 extinction."""
    if hasattr(rs,'transmitted_1064_energy'):
        shot_energy_1064 = rs.transmitted_1064_energy[:,np.newaxis].copy()
    else:
        print
        print "WARNING: 'transmitted_1064_energy' not found--assuming == 532 energy"
        print
        shot_energy_1064 = rs.transmitted_energy[:,np.newaxis].copy()
        
    angstrom_coef = processing_defaults.get_value('color_ratio','angstrom_coef')
    ones_array = np.ones_like(rs.combined_1064_counts)
    shot_energy_532 = rs.transmitted_energy[:,np.newaxis].copy()
    rayleigh_backscat = ones_array *inv.beta_r_backscat
    
   # inverse_scat_ratio = rayleigh_backscat / inv.beta_a_backscat
    count_ratio = rs.combined_1064_counts / rs.combined_counts 
    total_532_backscat = inv.beta_a_backscat + rayleigh_backscat
    mol_od = inv.optical_depth - inv.optical_depth_aerosol
    inv.beta_a_1064_backscat = total_532_backscat \
        *(shot_energy_532 / shot_energy_1064)\
        *count_ratio / (constants['IR_combined_hi_gain_ratio'] * (corr_adjusts['1064_gain'] if '1064_gain' in corr_adjusts else 1.0))\
        *np.exp(-30.0 * mol_od / 16)\
        *np.exp(-2 * (1 - 0.5**angstrom_coef) * inv.optical_depth_aerosol)\
        - rayleigh_backscat/16.0
    inv.beta_r_1064 = inv.beta_r_backscat * np.pi / 6.0 # 2-d vector of 1064 Rayleight extinction
    return 
