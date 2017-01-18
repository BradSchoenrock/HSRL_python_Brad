import sys
import numpy as np
import copy
from lg_base.core.fmt_mpl_date import fmt_mpl_datetime, fmt_mpl_time
import hsrl.calibration.calibration_utilities as cu
from datetime import datetime,timedelta
from scipy.interpolate import UnivariateSpline
import lg_base.core.array_utils as hau
import hsrl.filters.savitzky_golay as sg
from hsrl.filters.filtered_extinction import filtered_extinction
from hsrl.filters.cubic_fit import cubic_fit
from hsrl.filters.quadratic_fit import quadratic_fit
from hsrl.filters.fifth_order_fit import fifth_order_fit
from hsrl.filters.seventh_order_fit import seventh_order_fit
from hsrl.filters.ninth_order_fit import ninth_order_fit
import hsrl.filters.constrained_lsq as clsq
import hsrl.polarization.HSRL_PolCorrection as hp
import lg_base.core.bit_ops as bo
import lidar.lidar_utilities as lu
from time import sleep

import lg_base.core.decoratortools as nt


# Try to use the much faster nanmean from bottleneck, otherwise fall back
# to the scipy.stats version

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

def deltaz(alts):
    dz = alts.copy()
    dz[:-1]=alts[1:]-alts[:-1]
    dz[-1]=dz[-2]
    return dz

def process_data( instrument, start, end, min_alt, max_alt, requested_times
               ,rs_mem, rs_Cxx, rs_constants, rs_cal, sounding, corr_adjusts, processing_defaults
               ,compute_stats ,remainder):
             
    """def process_data(instrument,start,end,min_alt,max_alt,time_res\
               ,rs_raw,rs_mem,rs_Cxx,rs_constants,rs_cal,corr_adjusts\
               ,process_defaults,compute_stats):
        instrument    = 'gvhsrl','ahsrl','nshsrl','mf2hsrl'
        start         = start time of data interval as a python datetime
        end           = end time of data interval as a python datetime
        min_alt       = process data above this altitude (m)
        max_alt       = process data up to this altitude (m)
        requested_times= time vector requested for return values
        rs_mem        = raw data buffer that must contian requested data segment
        rs_Cxx        = inversion coeficients for the requested time period
        rs_constants  = system constants for requested time period
        rs_cal        = correction vectors (geo corr etc.) for requested time period
        compute_stats = if True compute photon counting statistics
      
        calling routine must insure that rs_mem, rs_Cxx, rs_constants, and rs_cal
        contain current data and calibrations for the requested segment.
        Outputs and inputs are packaged in a single structure, 'rs', containing---
        rs_Cxx, rs_inv, rp=rp, rs_stats, and rs_raw"""

    print 'process_data: from %s to %s' % \
            (fmt_mpl_datetime(start), fmt_mpl_datetime(end))

    rs = hau.Time_Z_Group(like=rs_mem)
    max_alt = float(max_alt)
    min_alt = float(min_alt)

  

    # select just those shots in current time interval for processing
    # note that raw is a new variable not just a pointer to rs_mem
  

    raw = copy.deepcopy(rs_mem)#hau.trimTimeInterval(  rs_mem, start, end )
    if 0:
        import matplotlib.pylab as plt
        [jnk,bins] = raw.molecular_counts.shape
        plt.figure(1)
        #plt.plot(np.nanmean(raw.molecular_counts,0),np.arange(bins))
        plt.plot(np.nanmean(raw.combined_hi_counts,0)/np.nanmean(raw.molecular_counts,0)\
                 *np.nanmean(raw.molecular_counts[0,500])/np.nanmean(raw.combined_hi_counts[:,500])\
                 ,np.arange(bins))                                                                          
        plt.title('rs_raw.molecular_counts')
        plt.ylabel('bin_number')
        plt.xlabel('counts/shot/bin')
        ax=plt.gca()
        #ax.set_xscale('log')
        plt.grid(True)
        
    # range dependent processing takes place at native range resolution of data

    binwidth = rs_constants['binwidth'] * 1.5e8
    if raw is None or not hasattr(raw,'molecular_counts') or len(raw.molecular_counts) == 0:
        rs_mean=remainder
        remainder=None
        print '*************** remainder swap ',rs_mean
    else:
        max_range_bin = raw.molecular_counts.shape[1]
        max_range = max_range_bin * binwidth   
        
        # for ground based operation we can shorten max_range without signal_in_dark_cor
        if 'installation' not in rs_constants \
               or not rs_constants['installation'] == 'airborne': 
           # if signal in dark count correction is turned off
            if not processing_defaults.enabled('signal_in_dark') or corr_adjusts['signal_in_dark'] == 0:
                max_range = max_alt + rs_constants['apd_pulse_timing'][1] \
                    * 1.5e8  # speed of light/2
            else:
                # if signal_in_dark_count is enabled
                max_range = 3601 * rs_constants['binwidth'] * 1.5e8
        
        # range_process accepts rs_raw as function of range from the lidar. It does all
        # of the processing that is range dependent, r_mean is returned as function of range
        rs_mean = range_process( instrument, raw, max_range / 1000, rs_constants
            ,rs_cal, rs_Cxx, corr_adjusts ,processing_defaults)
        
        if 0:
            import matplotlib.pylab as plt
            [jnk,bins] = rs_mean.molecular_counts.shape
            plt.figure(2)
            #plt.plot(np.nanmean(raw.molecular_counts,0),np.arange(bins))
            plt.plot(np.nanmean(rs_mean.combined_hi_counts,0)/np.nanmean(rs_mean.molecular_counts,0)\
                 ,np.arange(bins))                                                                          
            plt.title('rs_raw.molecular_counts')
            plt.ylabel('bin_number')
            plt.xlabel('counts/shot/bin')
            ax=plt.gca()
            #ax.set_xscale('log')
            plt.grid(True)
        
        s_time = datetime.utcnow()
        #up to this point lidar profiles include pre-trigger bins that were used for dark count
        rs_mean = convert_range_to_msl( rs_mean, rs_constants, raw, min_alt, max_alt, compute_stats )
        print 'time for convert_range_to_msl', datetime.utcnow() - s_time
        if 0:
            import matplotlib.pylab as plt
            [jnk,bins] = rs_mean.molecular_counts.shape
            plt.figure(3)
            #plt.plot(np.nanmean(raw.molecular_counts,0),np.arange(bins))
            plt.plot(np.nanmean(rs_mean.combined_hi_counts,0)/np.nanmean(rs_mean.molecular_counts,0)\
                 ,np.arange(bins))                                                                          
            plt.title('rs_raw.molecular_counts')
            plt.ylabel('bin_number')
            plt.xlabel('counts/shot/bin')
            ax=plt.gca()
            #ax.set_xscale('log')
            plt.grid(True)
            
        # add time vectors to the r_ms structure
        rs_mean.mergeTimeVectors( raw )
        #if fewer than n_time_ave shots are available, r_ms.times
        #will return [], no need to continue range processing.
        if rs_mean.times.shape[0] == 0:
            print 'attempt to process empty chunk'
            return None,remainder
       
        # trim inversion coefficents to requested altitude range
        # This gives us an array of altitudes to compute meanAltitudes
        # 
        # NOTE: the trimmed rs_Cxx is passed back to the caller
        rs_Cxx.trimByAltitude(min_alt, max_alt)
           #the vertical resolution is still at raw value 
        rs_mean.meanAltitudes( rs_Cxx.msl_altitudes )
           #vertical resolution now averaged to requested output resolution
      
        #correct variances for range averaging 
        sqrt_n_ave_alt = np.sqrt(rs_mean.n_ave_altitude)     
        for field in ['var_raw_combined_hi_counts'
                ,'var_raw_combined_lo_counts','var_raw_molecular_i2a_counts'
                ,'var_raw_molecular_counts','var_raw_cross_pol_counts'
                ,'var_raw_combined_1064_counts']:    
            if hasattr(rs_mean,field):
                setattr(rs_mean,field, getattr(rs_mean, field)/sqrt_n_ave_alt)
                 
    if rs_mean is None:
        rs_mean=remainder
        remainder=None
    else:
        rs_mean,remainder=rs_mean.includingRemainder(remainder)

    if rs_mean is None or rs_mean.times.shape[0] == 0:
            print 'attempt to process empty chunk'
            return None,remainder
    
    #  optional for polarization processing
    if hasattr(rs_mean,'qw_rotation_angle'):
        print 'Mean QWP Rotation Rate in deg/s:  ', np.mean(np.diff(rs_mean.qw_rotation_angle))*2
    if rs_constants['quarter_wave_plate_rotation'] == 'rotating':
        pc = hp.HSRL_PolCorrection(rs_mean, rs_cal, requested_times )
        pc.computeCorrection()
   
    mask_dataonly(rs_mean,match=match_invertable)

    if requested_times is not None:
        remainder=rs_mean.hereGoneBinTimes(requested_times,remainder=remainder,withRemainder=True)
        #print 'post trim ',requested_times,rs_mean,remainder
        #sleep(5)

    if rs_constants['quarter_wave_plate_rotation'] == 'rotating':
        pc.applyCorrection(rs_mean)

    if processing_defaults.enabled('m_smooth'):
        dz=rs_mean.msl_altitudes[2]-rs_mean.msl_altitudes[1]
        window =int(np.ceil(processing_defaults.get_value('m_smooth','window')/dz))
        rs_mean = molecular_smoothing(rs_mean,rs_constants,processing_defaults)

    if sounding is not None:
        # perform hsrl inversion 
        inv=process_inv_step(rs_mean,rs_Cxx,sounding, rs_constants,corr_adjusts,processing_defaults)
       
        # Package output in one structure
        rs.rs_inv = inv
    rs.rs_mean = rs_mean
    rs.rs_raw = raw
    return rs,remainder

def expand_atmospheric(rs_inv,sounding,rs=None):
    if rs is None:
        rs = hau.Time_Z_Group(like=sounding)
    for f,v in vars(sounding).items():#('temps','pressures'):
        if f.startswith('_') or hasattr(rs,f):
            continue
        if isinstance(v,hau.TZ_Array):
            v=hau.TZ_Array(np.ones((rs_inv.times.size,1))*v[:,:])
        elif isinstance(v,hau.T_Array):
            continue
        elif isinstance(v,hau.Z_Array):
            tmp=v
            v=hau.TZ_Array(np.transpose(tmp[:,np.newaxis]*np.ones((1,rs_inv.times.size))))
            v=np.require(v,requirements=['C'])
        else:
            continue
        setattr(rs,f,v)
    return rs

def maybetuple(ax,val,count):
    ret=[]
    for x in range(count):
        if x==ax:
            ret.append(val)
        else:
            ret.append(slice(None))
    if count==1:
        return ret[0]
    return tuple(ret)

def match_invertable(name):
    for k in ('_counts','_shots','shot_count','delta_t','telescope','roll_angle','pitch_angle','GPS_MSL'):
        if k in name:
            return True
    return False

def datamask(rs):
    indices=np.arange(rs.times.shape[0])
    indices=(indices>=0)#boolean now
    if not hasattr(rs,'op_mode'):
        return indices
    indices[np.bitwise_and(rs.op_mode,16)!=0]=False
    return indices

def mask_dataonly(rs,names=None,match=None):
    indices = datamask(rs)
    if indices.all():
        return indices
    for k in (names or vars(rs).keys()):
        if not hasattr(rs,k):
            continue
        if match is not None and not match(k):
            continue
        v=getattr(rs,k)
        if not isinstance(v,hau.T_Array):
            continue
        if 'O' in str(v.dtype):
            fill = None
        elif 'uint' in str(v.dtype):
            fill=0
        elif 'int' in str(v.dtype):
            fill=-1
        else:
            fill=np.NAN
        v[maybetuple(0,~indices,len(v.shape))]=fill
    return indices

def select_dataonly(rs):
    indices=datamask(rs)
    if indices.all():
        return rs,indices
    rs_masked = hau.Time_Z_Group(like=rs)
    for k,v in vars(rs).items():
        if isinstance(v,hau.T_Array):
            v=v[maybetuple(0,indices,len(v.shape))]
        setattr(rs_masked,k,v)
    return rs_masked,indices

def process_inv_step(rs_mean,rs_Cxx,sounding,rs_constants,corr_adjusts,processing_defaults):
    #rs_mean,_=select_dataonly(rs_mean)
    inv = cu.hsrl_inversion(rs_mean, rs_Cxx, rs_constants,corr_adjusts,processing_defaults)

    #add attenuated backscatter to inv structure
    compute_atten_backscat(rs_mean,inv,rs_Cxx,rs_constants,processing_defaults)
    
    if hasattr(rs_mean,'combined_1064_counts'):
    
        #compute 1064/532 color ratio if 1064 channel data is available
        inv.raw_color_ratio=cu.compute_raw_color_ratio(rs_mean,rs_Cxx,rs_constants,corr_adjusts)
        compute_atten_1064_backscat(rs_mean,inv,rs_Cxx,processing_defaults,rs_constants)
        
    if processing_defaults.enabled('klett'):
        ref_altitude = processing_defaults.get_value('klett','ref_altitude') * 1000.0
        if ref_altitude < rs_mean.msl_altitudes[0] or ref_altitude > rs_mean.msl_altitudes[-1] :
            print
            print
            print 'WARNING---klett ref altitutde=',ref_altitude, ' is not in requested altitudes'
            print 'no klett retrieval attempted'
            print
        else:    
          #for 1064 nm
          if hasattr(rs_mean,'combined_1064_counts'):
            LR_ratio_1064 = processing_defaults.get_value('klett','lidar_ratio_1064') 
            inv.beta_a_1064_backscat_klett = lu.compute_klett_backscatter(rs_mean.combined_1064_counts \
                ,inv.beta_r_backscat/16.0,rs_mean.msl_altitudes,LR_ratio_1064,ref_altitude)
          #for 532 nm    
          if hasattr(rs_mean,'combined_counts'):
            LR_ratio_532 = processing_defaults.get_value('klett','lidar_ratio_532') 
            inv.beta_a_backscat_klett = lu.compute_klett_backscatter(rs_mean.combined_counts \
                ,inv.beta_r_backscat,rs_mean.msl_altitudes,LR_ratio_532,ref_altitude)

    
    #compute integrated backscatter and add  inv structure 
    temp = inv.beta_a_backscat.copy()
    temp[np.isnan(temp)]=0.0
    dz=deltaz(inv.msl_altitudes)
    inv.integrated_backscatter = np.cumsum(temp,1)*dz   
    
    #make qc_mask--- inv.qc_mask and rs_mean.qc_mask are the same thing
    inv.qc_mask = hau.TZ_Array((np.ones_like(rs_mean.molecular_counts)*65535).astype('uint16'),summode='or')
    #inv.qc_mask = rs_mean.qc_mask.copy()
    s_time =datetime.utcnow()
    
    if hasattr(rs_mean,'combined_1064_counts'):
        make_1064_mask(rs_mean,inv,rs_constants)
        
    if processing_defaults.enabled('signal_lost_mask'):
        mol_lost_level = processing_defaults.get_value('signal_lost_mask','lost_level')
        if not rs_constants['installation'] == 'airborne':
           if not hasattr(rs_mean,'molecular_i2a_counts'):
               #ground based lidar
               make_hsrl_mask_simple(
                  inv.qc_mask,rs_mean.molecular_counts, mol_lost_level)
           else:
               make_hsrl_mask_simple(
                  inv.qc_mask,rs_mean.molecular_counts, mol_lost_level
                  ,rs_mean.molecular_i2a_counts)
        else:
           inv.qc_mask &= make_hsrl_mask(rs_mean,rs_constants,mol_lost_level)

        #inv.qc_mask = inv.qc_mask.copy()   
    print 'time for qc_mask = ',datetime.utcnow() - s_time
    if processing_defaults.enabled('cloud_mask'):
        s_time = datetime.utcnow()
        print 'cloud mask is enabled'
                                             
        make_cloud_mask(inv,rs_mean,processing_defaults,rs_constants)
        print 'time for cloud mask = ',datetime.utcnow() - s_time
        
    #generate mask indicating when I2 lock is lost
    if processing_defaults.enabled('I2_lock_mask'):
        make_i2_lock_mask(rs_mean,inv,processing_defaults,rs_constants)
      
    # compute photon counting statistics if requested.
    if processing_defaults.enabled('compute_stats'):
        s_time = datetime.utcnow()    
        compute_photon_statistics(rs_mean,inv, rs_Cxx,
                rs_constants)
        print 'time for photon statistics = ',datetime.utcnow() - s_time
        if processing_defaults.enabled('mol_signal_to_noise_mask'):
            make_mol_signal_to_noise_mask(inv,processing_defaults,rs_constants)
        if processing_defaults.enabled('particulate_backscatter_signal_to_noise_mask'):
            make_particulate_backscatter_signal_to_noise_mask(inv,processing_defaults)
            
    # note telsecope pointing from first profile--it does not make sense
    # to mix upward and downward pointing optical depth profiles
    if hasattr(rs_mean,'telescope_pointing'):
        [inv.optical_depth, inv.optical_depth_aerosol, inv.mol_norm_index,inv.mol_ref_aod] = \
                 lu.compute_optical_depth(inv.Nm, inv.beta_r_backscat
                 ,rs_Cxx.msl_altitudes, processing_defaults,rs_constants,rs_mean.telescope_pointing)
    else:
        [inv.optical_depth, inv.optical_depth_aerosol, inv.mol_norm_index,inv.mol_ref_aod] = \
                 lu.compute_optical_depth(inv.Nm, inv.beta_r_backscat
                 ,rs_Cxx.msl_altitudes, processing_defaults,rs_constants)
    
    
    #add aerosol 1064 aerosol backscatter and color ratio
    #to inv if system has 1064nm channel
    if hasattr(rs_mean,'combined_1064_counts'):
         cu.compute_1064_aerosol_backscatter(rs_mean
                     ,inv,processing_defaults,rs_constants,corr_adjusts)
         cu.compute_color_ratio(inv)
    if processing_defaults.enabled('include_expanded_profiles',return_if_missing=True): #omitted or explicitly not enabled
        expand_atmospheric(inv,sounding,inv)
    return inv

def compute_atten_1064_backscat(rs,inv,rs_Cxx,processing_defaults,constants):
    ref_alt = constants['lidar_altitude'] \
            + processing_defaults.get_value('atten_backscat_norm_range','range')
    indices = np.arange(rs.msl_altitudes.shape[0])
    if ref_alt < rs.msl_altitudes[1]:
        ref_alt = rs.msl_altitudes[2]
    ref_index = indices[rs.msl_altitudes <= ref_alt][-1]
    norm_factor = 0.5 * inv.beta_a_backscat[:,ref_index]\
                   / rs.combined_1064_counts[:,ref_index]       
    inv.attenuated_1064_backscatter = rs.combined_1064_counts[:,:] \
              * norm_factor[:,np.newaxis]
   
    #normalize assuming angstrom ratio = 1 at rs.msl_altitutes[ref_index]
    #seeded_shots_array = rs.seeded_shots[:, np.newaxis]
    #norm_factor = 0.5 * nanmean(inv.beta_a_backscat[:,ref_index])\
    #               /nanmean(rs.combined_1064_counts[:,ref_index])             
   
    #inv.attenuated_1064_backscatter = rs.combined_1064_counts[:,:] \
    #          * norm_factor *nanmean(seeded_shots_array[:,0])/ seeded_shots_array


   
    return

def compute_atten_backscat(rs,inv,rs_Cxx,constants,processing_defaults):
 # compute atten backscat normalized to beta_a_backscat at ref alt
    # sum polarization components to get total unpolarized return
    total_counts = rs.combined_counts \
        + constants['combined_to_cross_pol_gain_ratio'] \
        * rs.cross_pol_counts
    
    norm_factor = np.ones( (rs.times.shape[0],1) )
    if not('installation' in constants) \
              or (constants['installation']== 'ground'
                 or constants['installation'] == 'shipborne') \
              or  rs.GPS_MSL_Alt[0] < 0:  # ground based lidar
        ref_alt = constants['lidar_altitude'] \
            + processing_defaults.get_value('atten_backscat_norm_range','range')
        indices = np.arange(rs.msl_altitudes.shape[0])
        if ref_alt < rs.msl_altitudes[ 1]:
            ref_alt = rs.msl_altitudes[ 2]
        ref_index = indices[rs.msl_altitudes <= ref_alt][-1]
       
        # normalize to hsrl beta, add in molecular scattering for normalization
        if constants['polarization_is_linear']:
            norm_factor = np.transpose((inv.beta_a_backscat_par[:,
                ref_index] * (1 + inv.linear_depol[:,ref_index])
                + rs_Cxx.beta_r[ref_index] * 3 / (np.pi * 8))
                / total_counts[:,ref_index])
        else:    
            norm_factor = np.transpose((inv.beta_a_backscat_par[:,
                ref_index] * (1 + inv.circular_depol[:,ref_index])
                + rs_Cxx.beta_r[ref_index] * 3 / (np.pi * 8))
                / total_counts[:,ref_index])
        inv.atten_beta_a_backscat = np.transpose(np.transpose(total_counts)\
               *norm_factor) 
    else:
        # airborne lidar
        norm_range = processing_defaults.get_value('atten_backscat_norm_range','range')
        indices = np.arange(rs.msl_altitudes.shape[0])
        ref_index = np.zeros(rs.times.shape)
        norm_factor = np.zeros( (rs.times.shape[0], 1) )
        ref_alt=np.zeros(rs.times.shape)

        # airplane altitude may change--new norm altitude value for each profile
        for i in range(len(rs.times)):
            if rs.telescope_pointing[i] > 0.9:  # telescope pointing upward
                ref_alt[i] = rs.GPS_MSL_Alt[i] + norm_range \
                    * np.cos((rs.roll_angle[i]-4.0) * np.pi / 180.0)\
                       * np.cos(rs.pitch_angle[i] * np.pi / 180.0)
            elif rs.telescope_pointing[i] < 0.1:
                # telescope is pointing downward
                ref_alt[i] = rs.GPS_MSL_Alt[i] - norm_range \
                    * np.cos((rs.roll_angle[i]+4.0) * np.pi / 180.0)\
                       * np.cos(rs.pitch_angle[i] * np.pi / 180.0)
            else:
                norm_factor[i]=np.NaN #NOOP for jumbled orientation
                continue
            # don't try to use data below the ground
            ref_index = 0
            if ref_alt[i] > 0:
                ix = indices[rs.msl_altitudes <= ref_alt[i]]
                if ix.size:
                    ref_index = ix[-1]
            
            norm_factor[i] = (inv.beta_a_backscat_par[i, ref_index]
                              * (1 + inv.circular_depol[i, ref_index])
                              + rs_Cxx.beta_r[ref_index] * 3
                              / (np.pi * 8)) / total_counts[i, ref_index]
        inv.atten_beta_a_backscat = total_counts * norm_factor
    return   
        


def molecular_smoothing(rs,constants,processing_defaults):            
        dz=rs.msl_altitudes[2]-rs.msl_altitudes[1]
        window =int(np.ceil(processing_defaults.get_value('molecular_smooth'
             ,'window_length')/dz))
        n_average = dz /np.float(constants['binwidth']*1.5e8)
        
        #window must be odd
        if window % 2 != 1:
            window = window + 1
        print 'mol smoothing window = ',window,' points, ',dz*window, '(m)'
    
        
        #if window is long enough to support this polynomial
        if window > processing_defaults.enabled('m_smooth'): 
            [ntimes,nbins]=rs.molecular_counts.shape
            
            #if ground based don't let window contain lowest bins
            first_bin_to_smooth = int(window/2)\
               + int(np.float(processing_defaults.get_value(
               'first_bin_to_process','bin_number')\
               /n_average)) 

            #if ground based
            if not constants.has_key('installation') \
                        or  (constants['installation'] == 'ground'
                             or constants['installation'] == 'shipborne'):
                first_bin_to_smooth=first_bin_to_smooth+int(constants['lidar_altitude']/dz)
                for i in range(ntimes):
                    temp=rs.molecular_counts[i,:first_bin_to_smooth].copy()
                    rs.molecular_counts[i,:]\
                       =sg.savitzky_golay(rs.molecular_counts[i,:]
                       ,window,processing_defaults.get_value('molecular_smooth'
                       ,'polynomial_order'),deriv = 0)
                    rs.molecular_counts[i,:first_bin_to_smooth]=temp
            else:   #if airborne
        
                for i in range(ntimes):
                
                    if np.isnan(rs.GPS_MSL_Alt[i]):
                        continue
                    aircraft_alt_index = int(rs.GPS_MSL_Alt[i]/dz)
                    #if downlooking
                    max_index = rs.molecular_counts.shape[1]
                    if rs.telescope_pointing[i] < 0.1:
                        indices=range(aircraft_alt_index-first_bin_to_smooth-2,
                                      min( max_index, aircraft_alt_index) )
                    elif rs.telescope_pointing[i] > 0.9:  #if uplooking
                        # limit indices to size of array
                        indices=range( min(max_index,
                                       first_bin_to_smooth+aircraft_alt_index+2) )
                    else:
                        rs.molecular_counts[i,:]=np.NAN
                        continue#NOOP for scattered pointing
                
                    temp=rs.molecular_counts[i,indices].copy()
                    rs.molecular_counts[i,:]\
                       =sg.savitzky_golay(rs.molecular_counts[i,:]
                       ,window,processing_defaults.get_value('molecular_smooth'
                       ,'polynomial_order'),deriv = 0)
                    rs.molecular_counts[i,indices]=temp
        else:
            print ' '
            print '******No optical depth smoothing applied--window is too short'
            print ' '
        
        return rs
    

def alt_kludge(times, filename):
    import lg_base.core.read_utilities as ru
    from pycdf import CDF
    nc = CDF(filename)
    delta_times = nc.var('Time').get()
    alt = nc.var('ALT_G').get()
    roll = nc.var('ROLL').get()
    pitch = nc.var('PITCH').get()
    base_time_dict = ru.convert_date_str('6-jan-2012 18:35:19')
    base_time = base_time_dict['datetime']
    times_aircraft = base_time + timedelta(seconds=delta_times + 5.0)
    GPS_MSL_Alt = np.interp(times, times_aircraft, alt)
    GPS_MSL_Alt = hau.T_Array(GPS_MSL_Alt[:, np.newaxis])
    roll_angle = np.interp(times, times_aircraft, roll)
    roll_angle = hau.T_Array(roll_angle[:, np.newaxis])
    pitch_angle = np.interp(times, times_aircraft, pitch)
    pitch_angle = hau.T_Array(pitch_angle[:, np.newaxis])

    return (GPS_MSL_Alt, roll_angle, pitch_angle)


# convert measurement from range dependence to msl alititude dependence
def getAirborneAltitudes(telescope_pointing,roll_angle,pitch_angle,telescope_roll_angle_offset,range_vec,platformaltitude): 
  telescopemult=np.ones_like(telescope_pointing,dtype='float64')
  telescopemult[:]=np.NaN #this nan's profiles that are half up half down
  telescopemult[telescope_pointing>0.9]=1.0
  telescopemult[telescope_pointing<0.1]=-1.0
  full_roll_angle=(roll_angle + telescopemult * telescope_roll_angle_offset)*np.pi/180.0
  full_pitch_angle = pitch_angle * np.pi / 180.0
  coef=np.cos(full_roll_angle) * np.cos(full_pitch_angle) * telescopemult
  ntimes=telescope_pointing.shape[0]
  altitudes=np.zeros([ntimes,range_vec.shape[0]],dtype='float64')
  for i in range(ntimes):
    altitudes[i,:] = range_vec * coef[i] + platformaltitude[i]
  return altitudes

def interp(dest,src,dat):
    return np.interp(dest,src,dat,left=np.NaN,right=np.NaN)

def resampleZ(dest,src,dat):
  ret=np.zeros([dest.shape[0]],dtype='float64')
  ret[:]=interp(dest,src,dat[:])
  return ret
    
def resampleAllTimes(dest,src,data):
  ntimes=data.shape[0]
  ret=np.zeros([ntimes,dest.shape[0]],dtype='float64')
  for i in range(ntimes):
    ret[i,:]=interp(dest,src,data[i,:])
  return ret

def resampleAllTimes2D(dest,srcs,data):
  ntimes=data.shape[0]
  ret=np.zeros([ntimes,dest.shape[0]],dtype='float64')
  uporder=np.arange(0,data.shape[1])
  downorder=np.arange(data.shape[1]-1,-1,-1)
  for i in range(ntimes):
    if np.isnan(srcs[i,0]):
        ret[i,:]=np.NaN
    else:
        if srcs[i,1]<srcs[i,0]:
            order=downorder
        else:
            order=uporder
        ret[i,:] = interp(dest,srcs[i,order],data[i,order])
  return ret

def resampleZ2D(dest,srcs,data):
  ntimes=srcs.shape[0]
  ret=np.zeros([ntimes,dest.shape[0]],dtype='float64')
  uporder=np.arange(0,data.shape[0])
  downorder=np.arange(data.shape[0]-1,-1,-1)
  for i in range(ntimes):
    if np.isnan(srcs[i,0]):
        ret[i,:]=np.NaN
    else:
        if srcs[i,1]<srcs[i,0]:
            order=downorder
        else:
            order=uporder
        ret[i,:] = interp(dest,srcs[i,order],data[order])
  return ret

def convert_range_to_msl( rp, rs_constants, rs_raw, min_alt, max_alt, compute_stats):
    """convert_range_to_msl(rp,rs_constants,rs_raw,min_alt,max_alt,compute_stats)

    transform from range based bins to MSL alitutude based bins
    rp           = structure returned by range_process.py
    rs_constants = structure returned by select_system_constants.py
    rs_raw       = raw data returned by read_raw.py
    min_alt      = lowest altitude to process (m)
    min_alt      = highest altitude to process (m)
    compute_stats= compute photon counting stats if True"""


    #bin width in meters
    binwidth = 1.5e8 * float(rs_constants['binwidth'])  

    [dark_interval_end_time, laser_pulse_time, cal_pulse_end_time] = \
        rs_constants['apd_pulse_timing']

    # laser pulse bin number
    s_bin = int(laser_pulse_time / rs_constants['binwidth'])  
    [ntimes, nalts] = rp.molecular_counts.shape

    #resampleVars=['molecular_counts','combined_hi_counts','combined_counts','combined_lo_counts',\
    #        'cross_pol_counts','molecular_i2a_counts','combined_i2a_counts','combined_wfov_counts'\
    #        ,'molecular_wfov_counts','combined_1064_counts','var_raw_molecular_counts','geo_corr'\
    #        ,'var_raw_combined_hi_counts','var_raw_cross_pol_counts']
   
    if ('installation' not in rs_constants or rs_constants['installation'] == 'ground' \
       or rs_constants['installation'] == 'shipborne') : 

        #get telescope zenith angle, degrees converted to radians    
        telescope_zenith_angle = \
               np.abs(rs_constants['telescope_roll_angle_offset'])*np.pi/180.0

        print 'ground based lidar, zenith_angle = '\
              ,telescope_zenith_angle*180/np.pi

        rs_out = rp
        #input altitude vector
        altitudes = np.arange(rp.molecular_counts.shape[1]) * binwidth \
             * np.cos(telescope_zenith_angle) + rs_constants['lidar_altitude']\
             - s_bin * np.cos(telescope_zenith_angle)*binwidth

        #output altitude vector
        nalts_out = int(max_alt / binwidth)
        
        newalts=np.arange(len(altitudes)) * binwidth
        for varn in vars(rs_out).keys():
            varv=getattr(rs_out,varn)
            if isinstance(varv,hau.TZ_Array):
                v=resampleAllTimes(newalts,altitudes,varv)
                setattr(rs_out,varn,hau.TZ_Array(v,summode=varv.summode))
            elif isinstance(varv,hau.Z_Array) and not isinstance(varv,hau.T_Array):
                v=resampleZ(newalts,altitudes,varv)
                setattr(rs_out,varn,hau.Z_Array(v,summode=varv.summode))                
        rs_out.msl_altitudes = hau.Z_Array(newalts)
        rs_out.binwidth = binwidth
        rs_out.lidar_altitude = rs_constants['lidar_altitude']
        rs_out.trimByAltitude(min_alt, max_alt)
                
        return rs_out

    elif rs_constants['installation'] == 'airborne' :

        print 'airborne lidar'

        rs_out = rp
              
        range_vec = np.arange(rp.molecular_counts.shape[1])\
                   *binwidth -  s_bin*binwidth
       
        
        #output altitude vector
        nalts_out = int(max_alt / binwidth) 
        #newalts = np.arange(len(alt_vec)) * binwidth
        newalts = np.arange(nalts_out) * binwidth
        
        altitudes=getAirborneAltitudes(rs_raw.telescope_pointing,rs_raw.roll_angle
               ,rs_raw.pitch_angle,rs_constants['telescope_roll_angle_offset']
               ,range_vec,rs_raw.GPS_MSL_Alt)
        for varn in vars(rs_out).keys():
            varv=getattr(rs_out,varn)
            if isinstance(varv,hau.TZ_Array):
                v=resampleAllTimes2D(newalts,altitudes,varv)
                setattr(rs_out,varn,hau.TZ_Array(v,summode=varv.summode))
            elif isinstance(varv,hau.Z_Array) and not isinstance(varv,hau.T_Array):
                v=resampleZ2D(newalts,altitudes,varv)
                setattr(rs_out,varn,hau.TZ_Array(v,summode=varv.summode))                
        rs_out.msl_altitudes = hau.Z_Array(newalts)
        rs_out.binwidth = binwidth

        rs_out.trimByAltitude(min_alt, max_alt)
      
        return rs_out
    
def gps_quality_check(rs_raw,rs_constants):
    """checks for gaps is aircraft GPS and INS data stream.
       If no data available, assume calvals_****.txt values.
       If data gaps use calvals_****.txt values to fill"""
    #This is called in input_translators because thats where raw can be modified
    # if gps data is available, fill any gaps with ground values:  
    if hasattr(rs_raw,'GPS_MSL_Alt'):    
        rs_raw.GPS_MSL_Alt[np.isnan(rs_raw.GPS_MSL_Alt)] = rs_constants['lidar_altitude']            
        rs_raw.latitude[np.isnan(rs_raw.GPS_MSL_Alt)] = rs_constants['latitude']
        rs_raw.longitude[np.isnan(rs_raw.GPS_MSL_Alt)] = rs_constants['longitude']       

        maximumAlt = 16000.0
        if np.any(np.greater(rs_raw.GPS_MSL_Alt, maximumAlt)):
                # replace huge altitude values with a default altitude
                print 'hsrl_read_utilities: Warning - Spikes in GPS altitude!'
                rs_raw.GPS_MSL_Alt[np.greater(rs_raw.GPS_MSL_Alt, maximumAlt)] = \
                    rs_constants['lidar_altitude']
    else:   #GPS is missing, assume aircraft is on the ground
            #get missing values from calvals_*****.txt file
        print 'mobile installation with missing fields'       
        rs_raw.GPS_MSL_Alt = hau.T_Array(rs_constants['lidar_altitude']*np.ones(rs_raw.times.shape))
        if not hasattr(rs_raw,'latitude'):
            rs_raw.latitude = hau.T_Array(rs_constants['latitude']*np.ones(rs_raw.times.shape))
        if not hasattr(rs_raw,'longitude'):
            rs_raw.longitude = hau.T_Array(rs_constants['longitude']*np.ones(rs_raw.times.shape))
        if not hasattr(rs_raw,'pitch_angle'):
            rs_raw.pitch_angle = hau.T_Array(np.zeros(rs_raw.times.shape))
        if not hasattr(rs_raw,'roll_angle'):
            rs_raw.roll_angle = hau.T_Array(np.zeros(rs_raw.times.shape))
        
    if hasattr(rs_raw,'pitch_angle'):
        try:
            if anynan(rs_raw.pitch_angle):
                print 'Pitch angle missing----replacing with zeros'
                rs_raw.pitch_angle[np.isnan(rs_raw.pitch_angle)] = \
                    np.zeros_like(rs_raw.pitch_angle[np.isnan(rs_raw.pitch_angle)])
                rs_raw.roll_angle[np.isnan(rs_raw.roll_angle)] =  \
                    np.zeros_like(rs_raw.roll_angle[np.isnan(rs_raw.roll_angle)])
        except TypeError:#not a float
            pass
    return rs_raw
    


def range_process( instrument, raw, max_range, constants
    ,rs_cal, rs_Cxx, corr_adjusts ,processing_defaults):
    """range_process(instrument, rs_in,max_range, n_range_ave, constants\
         ,rs_cal, signal_in_dark_cor, corr_adjusts,display_defaults, processing_defaults):
       completes signal correction and all processes that are dependent on range from lidar.
       The returned profiles including bins prior to the laser pulse. 
       instrument = instrument selection. eg. 'ahsrl','gvhsrl','mf2hsrl','nshsrl'
       raw                = time_z_object containing raw counts---all contents of raw must
                            be preserved unchanged for possible reprocessing of data. Data
                            will be returned in new structure 'rs'.
       max_range         = farthest range to process (km)
       constants         = dictionary containing system constants
                           created by cal_file_reader.py and select_system_contants.py
       rs_cal            = altitude dependent calibration coeficients
       signal_in_dark_cor= coeficients used to separate signal from dark count
                           , see signal_in_dark_count.py               
       corr_adjusts      = dictionary contianing scaling constants for signal corrections,
                           see 'get_scale_corrections.py'
       process_defaults  = dictionary containing processing directives
                           see 'process_defaults.json'
    """


   

    assert(rs_Cxx is not None)
    rs = hau.Time_Z_Group(like=raw)

    if 0:
        import matplotlib.pylab as plt
        
        mol = np.nanmean(raw.molecular_counts,0)
        wfov  = np.nanmean(raw.molecular_wfov_counts,0)
        bin_vec = 7.5 * np.arange(len(wfov))
        mol = mol - np.nanmean(mol[0:40])
        wfov = wfov - np.nanmean(wfov[0:40])
        mol *= (bin_vec-45*7.5)**2
        wfov *= (bin_vec-45*7.5)**2
        wfov *= np.exp(-2*bin_vec *1e-5)
        #wfov = wfov - bin_vec*wfov[900]/(900 *0.0001)
        wfov *= mol[900]/wfov[900]
        plt.figure(99999)
        plt.plot(bin_vec,wfov,'c',bin_vec,mol,'r')
        ax=plt.gca()
        ax.set_yscale('log')
        plt.grid(True)
        plt.show()
        print j
    #copy corrected raw into rs 
    for field in ['transmitted_1064_energy','transmitted_energy','seeded_shots','molecular_counts'
           ,'combined_lo_counts','combined_hi_counts','cross_pol_counts',
           'combined_wfov_counts','molecular_wfov_counts',
           'molecular_i2a_counts','combined_1064_counts','telescope_pointing']:
        if hasattr(raw,field):
            setattr(rs,field,getattr(raw,field).copy())
            setattr(rs,'raw_'+field,getattr(raw,field).copy())
   
    # compute bin number of laser pulse
    [dark_interval_end_time, laser_pulse_time, cal_pulse_end_time] = \
        constants['apd_pulse_timing']
    bin_duration = constants['binwidth']
    s_bin = int(laser_pulse_time / bin_duration)  # laser pulse bin number
    #dark_interval_end_bin = int(dark_interval_end_time / bin_duration)- 1

    nalts = raw.molecular_counts.shape[1]

    #save geo corrected raw counts as 'var_xxx' in rs so that they get averaged without
    #other range processing for use in compute_photon_statistics. We also multiply
    #by the square of the geocorrection to account for the geocorrection in the
    #signal used compute_phothon_statistics()
    if processing_defaults.enabled('compute_stats'):
        ones_array = np.ones(raw.molecular_counts.shape)
        # bin 0 of geo_correction is defined as occurring at the laser pulse
        geocorr = ones_array.copy()
        geocorr[:,s_bin:] = rs_cal.geo.data[:nalts-s_bin, 1] * ones_array[:,s_bin:]
    
        for field in ('molecular_counts','combined_lo_counts'
           ,'combined_hi_counts','cross_pol_counts','combined_wfov_counts'
           ,'molecular_wfov_counts','molecular_i2a_counts','combined_1064_counts'):
           if hasattr(raw,field):
               setattr(rs,'var_raw_'+field,getattr(raw,field)*geocorr*geocorr)       
            
    #counts arrays are the average number of counts in a data raw acquistion interval
    #of raw.times[2]-raw.times[1] while seeded_shots is the total number of laser pulses
    #the acquisition interval prior to preaveraging in preprocess_raw.py

    #note: this does not compensate for the pileup correction--in very high count areas this
    #will under estimate the varience because actual counts are multipled by a pileup correction
    #in the preprocess_raw.py routine  

    #counts have been pileup corrected in preprocessing
    #do dark correction for all channels
   
    s_time =datetime.utcnow()
    dark_count_correction(instrument,raw,rs,rs_Cxx,corr_adjusts,processing_defaults,constants)
    print 'time for dark correction = ',datetime.utcnow() - s_time

    # gain correction for nadir pointing in airborne operation
    # this is expected to be a very small correction with little
    # impact on signal statitics
    if 'installation' in constants and constants['installation'] == 'airborne' \
           and constants['nadir_comb_gain_adjustment'] != 1.0:
        print 'Apply nadir gain adjustment'
        print 'nadir gain adj= ', constants['nadir_comb_gain_adjustment']
        ix = np.arange(rs.telescope_pointing.shape[0])
        indices = ix[rs.telescope_pointing[:] < 0.1]
        nadir_gain_adj = constants['nadir_comb_gain_adjustment']
        rs.combined_lo_counts[indices, :] *= nadir_gain_adj
        rs.combined_hi_counts[indices, :] *= nadir_gain_adj
    
    #np.set_printoptions(threshold='nan')
   
    #do baseline correction
    rs = baseline_correction(rs,rs_cal,nalts,corr_adjusts,constants)
                 
    # correct for differential geometry between 1064 and 532 nm channels
    rs = diff_1064_532_geometry_correction(rs,rs_cal,nalts,processing_defaults
                                      ,corr_adjusts)
    if 0:
        import matplotlib.pylab as plt
        plt.figure(67)
        plt.plot(np.nanmean(rs.combined_hi_counts,0),np.arange(len(rs.combined_hi_counts[0,:])),'r'
                    ,np.nanmean(rs.molecular_counts,0),np.arange(len(rs.molecular_counts[0,:])),'b')
        ax=plt.gca()
        ax.set_xscale('log') 
    #do combined-molecular differential geo correction if available
    rs = diff_geometry_correction(rs,rs_cal,nalts,processing_defaults
                                      ,corr_adjusts)
    if 0:
        import matplotlib.pylab as plt
        plt.figure(68)
        plt.plot(np.nanmean(rs.combined_hi_counts,0),np.arange(len(rs.combined_hi_counts[0,:])),'r'
                    ,np.nanmean(rs.molecular_counts,0),np.arange(len(rs.molecular_counts[0,:])),'b')
        ax=plt.gca()
        ax.set_xscale('log') 
   
    # Matt Add:  do cross polarization differential geometry correction
    rs = diff_cp_geometry_correction(rs,rs_cal,nalts,processing_defaults
                                      ,corr_adjusts)
      
    # do i2a differential geo correction if present and relavent to instrument
    if hasattr(rs,'molecular_i2a_counts') and corr_adjusts['i2a_dif_geo_corr'] > 0:
        rs = i2a_diff_geo_correction(rs,rs_cal,corr_adjusts)

    #create combined_counts from combined_hi and combined_lo profiles
    rs = merge_combined_hi_and_lo(rs,constants)
    if 0:
        import matplotlib.pylab as plt
        plt.figure(69)
        plt.plot(np.nanmean(rs.combined_hi_counts,0),np.arange(len(rs.combined_hi_counts[0,:])),'r'
                 ,np.nanmean(rs.molecular_counts,0),np.arange(len(rs.molecular_counts[0,:])),'b'
                 ,np.nanmean(rs.combined_lo_counts,0),np.arange(len(rs.combined_lo_counts[0,:])),'c'
                 ,np.nanmean(rs.cross_pol_counts,0),np.arange(len(rs.cross_pol_counts[0,:])),'g'
                 ,np.nanmean(rs.combined_counts,0),np.arange(len(rs.combined_counts[0,:])),'k')
        ax=plt.gca()
        ax.set_xscale('log')
        #plt.show()

        print 'cp/mol'
    """
    if processing_defaults.enabled('wfov_geo_corr') and hasattr(rs,'molecular_wfov_counts'):
        #do geometry correction after adjusting geo_corr with wide-field-of-view data.
        geo_corr = rs_cal.geo.data[:4000,1]
        s_bin = np.int(constants['apd_pulse_timing'][1]/constants['binwidth'])
        wfov_ratios = np.zeros(rs.molecular_wfov_counts.shape[1])
        wfov_ratios[:-s_bin] = nanmean(rs.molecular_wfov_counts[:,s_bin:],0)\
                  / nanmean(rs.molecular_counts[:,s_bin:],0) 
        wfov_geometry_correction(rs,wfov_ratios,geo_corr,processing_defaults,constants,corr_adjusts)
    """
    #does wfov corr exist?
    if processing_defaults.enabled('wfov_corr') and hasattr(rs,'molecular_wfov_counts')\
                              and hasattr(rs_cal,'geo')\
                              and hasattr(rs_cal.geo,'wfov_mol_ratio'):
    
       
        #add pre-trigger bins to wfov_mol_ratio array provided in geofile_default_file
        #and add to structure for use in extinction processing
        calibration_wfov_mol_ratio = np.zeros(rs.molecular_counts.shape[1])
        calibration_wfov_mol_ratio[s_bin:] = \
                       rs_cal.geo.wfov_mol_ratio[:(rs.molecular_counts.shape[1]-s_bin)]
        rs.calibration_wfov_mol_ratio = hau.Z_Array(calibration_wfov_mol_ratio)
     
    # do the normal geometric  correction on the following variables
    select = ['molecular_counts','combined_lo_counts','combined_hi_counts'
                  ,'molecular_i2a_counts','combined_1064_counts','molecular_wfov_counts'
                  ,'combined_counts','cross_pol_counts']
    rs = lu.geometry_correction(select,rs,rs_cal,nalts,s_bin,corr_adjusts['geo_corr'])
        
    #mask close range bin counts
    first_bin_to_process = processing_defaults.get_value('first_bin_to_process','bin_number')
    for field in ['combined_hi_counts','combined_lo_counts','combined_wfov_counts','molecular_wfov_counts'
        'molecular_i2a_counts','molecular_counts','cross_pol_counts','combined_counts'\
        'combined_1064_counts']:
        if hasattr(rs,field):
            getattr(rs,field)[:, :(s_bin+first_bin_to_process)] = np.NaN
 
    return rs
 
def merge_combined_hi_and_lo(rs,constants):
    """merge_comb_hi_and_lo(rs,constants):
        merge low channel counts into high when high is overloaded
        combined=low channel*gain  when combined_hi > threshold
        """

    if hasattr(rs,'combined_lo_counts'):
        seeded_shots_array2=rs.seeded_shots[:, np.newaxis]
        merge_mask = rs.combined_hi_counts / seeded_shots_array2 \
               > constants['combined_channel_merge_threshhold']
        rs.combined_counts = rs.combined_hi_counts.copy()
        rs.combined_counts[merge_mask] = \
               rs.combined_lo_counts[merge_mask] \
               * constants['hi_to_low_combined_channel_gain_ratio']
    else:
        rs.combined_counts = rs.combined_hi_counts.copy()
    return rs
                         

def make_i2_lock_mask(rs_mean,inv,processing_defaults,constants):
   """clear bits 0 and bit 2 of qc_mask if seed voltage deviates more than
      the tolerance provided in process_defaults['I2_lock_mask','lock_lost_offset'] 
      The shortcell ratio is scaled by values of the unlocked and locked ratios provided
      for the current instrument in constants['shortcell_locked_ratio'][0] and [1]
      action is taken if tolerance value is empty or process_defaults['I2_lock_mask'] is
      not found"""
   if constants.has_key('shortcell_locked_ratio')  and not constants['shortcell_locked_ratio'][0]==-9999 and hasattr(rs_mean,'filtered_energy'):
       ratio_lost = processing_defaults.get_value('I2_lock_mask','lock_lost')
       ratio_warning = processing_defaults.get_value('I2_lock_mask','lock_warning')
       [ntimes,nalts]=inv.qc_mask.shape
       ones_array_t=np.ones((nalts,ntimes)).astype('uint16')
       shortcell_locked_ratio = np.float(constants['shortcell_locked_ratio'][0])
       
       #note: shortcell_ratio is normalized by constants['shortcell_locked_ratio[0] and [1]
       #this is the ratio of the mean energies
       shortcell_ratio = (rs_mean.filtered_energy[:,0] / rs_mean.nonfiltered_energy[:,0] \
            -constants['shortcell_locked_ratio'][0])\
           / (constants['shortcell_locked_ratio'][1]-constants['shortcell_locked_ratio'][0]) 

       #this is the ratio of the mean of the max filtered energies over the mean nonfiltered
       if len(rs_mean.filtered_energy.shape)>1 and rs_mean.filtered_energy.shape[1]>2:
           max_shortcell_ratio = (rs_mean.filtered_energy[:,2] / rs_mean.nonfiltered_energy[:,0] \
            -constants['shortcell_locked_ratio'][0])\
           / (constants['shortcell_locked_ratio'][1]-constants['shortcell_locked_ratio'][0]) 
       else:
           max_shortcell_ratio=None
     
       #generate mask indicating loss of I2 lock or missing nonfiltered_energy
       mask = (shortcell_ratio > ratio_lost).astype('uint16')     
       mask[:] |= (np.isnan(rs_mean.nonfiltered_energy[:,0])) 
       mask = ~((2**2+1)*mask)
       mask =  np.transpose(mask * ones_array_t)
       #generate warning mask for deviations in shortcell ratio
       #when averge of max is greater than lost
       #or average is greater than warning
       if max_shortcell_ratio is not None:
           warning = ((max_shortcell_ratio > ratio_lost)
                  | (shortcell_ratio > ratio_warning)).astype('uint16')
       else:
           warning = (shortcell_ratio > ratio_warning).astype('uint16')           
      
       warning =~(2**13 * warning)
       warning = np.transpose(warning * ones_array_t)
       inv.qc_mask[:,:] &= mask.astype('uint16')
      

       inv.qc_mask[:,:] &= warning.astype('uint16')
     
       
   return
    
def make_cloud_mask( inv, rs_mean,processing_defaults, rs_constants):
    """make_cloud_mask(rs,rs_constants, cloud_threshhold,buffer_bins
                     max_cloud_alt,full_mask)
       generates a logical mask that becomes
       false after the particulate backscatter cross section exceeds the cloud
       cloud_threshhold in units of 1/(m sr).
       buffer_bins sets mask to False prior to cloud contact
       max_cloud_alt--ignore clouds above this altitude
       full_mask---mask entire profile if any cloud is found below max_cloud_alt"""
       
    np.set_printoptions(threshold=np.NaN)
    cloud_threshhold = processing_defaults.get_value('cloud_mask','backscat_threshhold')
    full_mask = processing_defaults.get_value('cloud_mask','mask_entire_profile')
    buffer_bins = np.int(processing_defaults.get_value('cloud_mask','cloud_buffer_zone')    
                      /(rs_mean.msl_altitudes[2]-rs_mean.msl_altitudes[1]))
    [ntimes, nalts] = inv.beta_a_backscat_par.shape
    temp = inv.beta_a_backscat_par.copy()
    temp = temp +inv.beta_a_backscat_perp
    temp[:, 0] = 0.0
    temp[temp > cloud_threshhold] = np.NaN
   

    
    # does not allow for shift from ground-based (i.e. no GPS)
    # to airborne within one record
    if ('installation' not in rs_constants
        or rs_constants['installation'] == 'ground' 
        or rs_constants['installation'] == 'shipborne'):  # lidar is on the ground looking up
        start_alt = rs_constants['lidar_altitude'] + 300
        start_alt = np.max([start_alt, rs_mean.msl_altitudes[0]])
        temp[:, rs_mean.msl_altitudes < start_alt] = 0.0
        mask = np.isfinite(np.cumsum(temp, 1)).astype('uint16')
        #apply a pre-trigger buffer on nbuf data points to mask
        mask[:,:(nalts-buffer_bins)] = np.bitwise_and(mask[:,:(nalts-buffer_bins)]
                        , mask[:,buffer_bins:])
        if full_mask == 1:
           max_cloud_alt = np.float(processing_defaults.get_value('cloud_mask','max_cloud_alt')) 
           index = len(rs_mean.msl_altitudes[rs_mean.msl_altitudes < max_cloud_alt*1000.0])
           for i in range(ntimes):
               if any(mask[i,:index] == False):
                  mask[i,:] = False
                          
    else:
        # lidar is airborne
        indices = np.arange(nalts)
        mask = np.zeros_like(temp).astype('uint16')
        for i in range(ntimes):
            if rs_mean.telescope_pointing[i] > 0.9:  # telescope pointing up
                ix = indices[rs_mean.msl_altitudes <= rs_mean.GPS_MSL_Alt[i] + 250]
                if len(ix) > 0:
                    start_index = np.max(ix)
                    mask[i, start_index:] = \
                        np.isfinite(np.cumsum(temp[i, start_index:]))\
                        .astype('uint16')
                #apply a pre-trigger buffer on nbuf data points to mask
                mask[i,:(nalts-buffer_bins)] = np.bitwise_and(mask[i,:(nalts-buffer_bins)]
                                      , mask[i,buffer_bins:])   
            elif rs_mean.telescope_pointing[i] < 0.1:
                # telescope is pointing down
                ix = indices[rs_mean.msl_altitudes <= rs_mean.GPS_MSL_Alt[i] - 250]
                if len(ix) > 0:
                    start_index = np.max(ix)
                    mask[i,start_index:0:-1] = \
                        np.isfinite(np.cumsum(temp[i, start_index:0: -1]))\
                        .astype('uint16')
                    #apply a pre-trigger buffer on nbuf data points to mask
                    if buffer_bins:
                       print 'buffer bins not implemented for nadir viewing'
                       
                       
                       #mask[i,buffer_bins:nalts] = \
                       #        np.bitwise_and(mask[i,:(nalts-buffer_bins)]\
                       #        , mask[i,buffer_bins:nalts])
               
                   
    #mask for bits 0 and 7   
    mask =   ~(129*(mask==0)).astype('uint16')
    inv.qc_mask &= mask
   
    
    #rs_mean.qc_mask =inv.qc_mask.copy()

    return

def make_mol_signal_to_noise_mask(inv,processing_defaults,rs_constants):
    mol_sn_threshhold =processing_defaults.get_value('mol_signal_to_noise_mask','threshhold')
    if rs_constants.has_key('no_i2_channel') and rs_constants['no_i2_channel'] == 1:
        temp = inv.SN_i2a_mol.copy()
    else:
        temp = inv.SN_mol.copy()

    if rs_constants['installation'] != 'airborne' :  # lidar is on the ground looking up
        start_alt = rs_constants['lidar_altitude'] +200.0
        start_alt = np.max([start_alt, inv.msl_altitudes[1]])
        temp[:, inv.msl_altitudes <= start_alt] = mol_sn_threshhold
        temp[temp < mol_sn_threshhold] = np.NaN
        #generate mask clearing bits 0 and 6 if cumsum(temp,1) == NaN
        mask = np.isfinite(np.cumsum(temp,1))== 0
        #multiply by 2**0+2**4 and complement to generate mask
        mask = ~(17*mask.astype('uint16'))
        inv.qc_mask &= mask
    elif rs_constants['installation'] == 'airborne':
        indices = range(len(inv.msl_altitudes))
        temp[temp < mol_sn_threshhold] = np.NaN
        mask=np.zeros_like(temp)
        for i in range(len(inv.GPS_MSL_Alt)):
             #uplooking telescope
             if inv.telescope_pointing[i] > .9:
                  start_index = sum(inv.msl_altitudes < (inv.GPS_MSL_Alt[i] + 300))
                  #generate mask clearing bits 0 and 6 if cumsum(temp,1) == NaN
                  mask[i,start_index:] = np.isfinite(np.cumsum(temp[i,start_index:]))== 0
             #downlooking telescope
             elif inv.telescope_pointing[i] < 0.1:
                 if inv.msl_altitudes[-1] > inv.GPS_MSL_Alt[i]+200 :
                     # consider only those altitudes more than 200 m below the plane
                     start_index = sum(inv.msl_altitudes < (inv.GPS_MSL_Alt[i] - 200))
                 else:
                     #if airplane is higher that 200m above the requested altitude 
                     start_index = np.int(sum(inv.msl_altitudes < inv.msl_altitudes[-1] - 200))

                 indices = range(start_index,-1,-1)
                 mask[i,indices] = np.isnan(np.cumsum(temp[i,indices]))
        #multiply by 2**0+2**4 and complement to generate mask
        mask = ~(17*mask.astype('uint16'))
        inv.qc_mask &= mask
    return

def make_particulate_backscatter_signal_to_noise_mask(inv,processing_defaults):
    SN_threshold =processing_defaults.get_value('particulate_backscatter_signal_to_noise_mask','threshhold')
    mask = np.zeroslike(inv.beta_a_backscat)
    mask[inv.SN_beta_a_backscat < SN_threshold] = 1
    #multiply by 2**0+2**5 and complement to generate mask
    mask = ~(33*mask.astype('uint16'))
    inv.qc_mask &= mask
    return

def make_1064_mask(rs_mean,rs_inv,rs_constants):
    temp = np.ones_like(rs_inv.qc_mask).astype('uint16')
    if hasattr(rs_mean,'IRDetectorShutterIn'):
        val=rs_mean.IRDetectorShutterIn[:,np.newaxis].astype('uint16') #1 where bad, 0 where good
        if 'IRDetectorShutterInverse' in rs_constants and rs_constants['IRDetectorShutterInverse']!=0: 
            val=1-val #inverted, so its good when shutter is not in. 1-0 = 1, 1-1 = 0
        mask = val * temp #make a 2-D matrix, true where bad
        mask = ~((2**15+2**0) * mask.astype('uint16'))
        rs_inv.qc_mask &= mask
    return

    
def make_hsrl_mask(rs,rs_constants,mol_lost_level):
    """make_hsrl_mask(rs, rs_constants, mol_lost_level)
       generates a logical mask that becomes false after the molecular signal first falls
       below the mol_lost_level. Mask values are stored in a bit field as follows:
       qc_mask, bit[0] = logical and of mask bits 0-->15
       qc_mask, bit[1] =
       qc_mask, bit[2] =
       qc_mask, bit[3] =
       qc_mask, bit[4] = cleared for altitude bins after mol S/N ratio fall below threshold
       qc_mask, bit[5] =
       qc_mask, bit[6] = cleared for altitude bins after molecular counts in bin <= mol_lost threshhold
       qc_mask, bit[12]= cleared for altitude bins after backscatter cross section > cloud_threshhold
       qc_mask, bit[13]= cleared for i2_lost_warning
       qc_mask, bit[15]= cleared for 1064 detector shutter closed"""
    
    # np.set_printoptions(threshold=np.NaN)

    [ntimes, nalts] = rs.molecular_counts.shape
    temp = rs.molecular_counts.copy()
    temp[:, 0] = mol_lost_level
    temp[temp < mol_lost_level] = np.NaN

    # does not allow for shift from ground-based (i.e. no GPS)
    # to airborne within one record
   
    if rs_constants['installation'] != 'airborne' :  # lidar is on the ground looking up
        start_alt = rs_constants['lidar_altitude'] +100.0
        start_alt = np.max([start_alt, rs.msl_altitudes[0]])
        temp[:, rs.msl_altitudes < start_alt] = mol_lost_level

        #generate mask clearing bits 0 and 6 if cumsum(temp,1) == NaN
        mask = np.isfinite(np.cumsum(temp,1))== 0
        #multiply by 2**0+2**6 and complement to generate mask
        mask = ~(65*mask.astype('uint16'))
            
    else:
        # lidar is airborne
        indices = np.arange(nalts)
        mask = np.zeros_like(temp).astype('uint16')
        for i in range(ntimes):
            if rs.telescope_pointing[i] > 0.9:  # telescope pointing up
                ix = indices[rs.msl_altitudes <= rs.GPS_MSL_Alt[i] + 250]
                if len(ix) > 0:
                    start_index = np.max(ix)
                    mask[i, start_index:] = \
                        ~(65*(np.isfinite(np.cumsum(temp[i, start_index:]))==0))\
                              .astype('uint16')
            elif rs.telescope_pointing[i] < 0.1:
                # telescope is pointing down
                ix = indices[rs.msl_altitudes <= rs.GPS_MSL_Alt[i] - 250]
                if len(ix) > 0:
                    start_index = np.max(ix)
                    
                    mask[i, start_index:0:-1] = \
                        ~(65*(np.isfinite(np.cumsum(temp[i, start_index:0:-1]))==0))\
                              .astype('uint16')
   
    return mask.astype('uint16')
 






def make_hsrl_mask_simple(qc_mask,molecular_counts,mol_lost_level,i2a_molecular_counts=None):
    """make_hsrl_mask_numba(rs, rs_constants, mol_lost_level)
       generates a logical mask that becomes false after the molecular signal first falls
       below the mol_lost_level. Mask values are stored in a bit field as follows:
       qc_mask, bit[0] = logical and of all other mask bits
       qc_mask, bit[1] =
       qc_mask, bit[2] =
       qc_mask, bit[3] =
       qc_mask, bit[4] =
       qc_mask, bit[5] =
       qc_mask, bit[6] = cleared for altitude bins after molecular counts in bin <= mol_lost threshhold
       qc_mask, bit[12]= cleared for altitude bins after backscatter cross section > cloud_threshhold"""

    # np.set_printoptions(threshold=np.NaN)

    if i2a_molecular_counts is None:
       m_counts = molecular_counts
    else:
       m_counts = molecular_counts +i2a_molecular_counts

    mask = np.uint16(65470)
    [ntimes, nalts] = molecular_counts.shape
    for i in range(ntimes):
       flag =1
       for j in range(nalts):
          if flag == 0:
              #qc_mask[i,j:] = np.logical_and(qc_mask[i,j:] , mask )
              qc_mask[i,j:] = np.bitwise_and(qc_mask[i,j:] , mask) 
              break
          #elif not true for NaN's --- thus ignores NaN's before start of data   
          elif molecular_counts[i,j] <= mol_lost_level:
             flag = 0
    return  









    
def  spline_smoothing(array,std_array,smoothing):
    """Uses a smoothing spline to fit array values,(applies to 2nd dimension only)
    array     = input and return array, array(ntimes,nbins)
    std_array = standard deviation estimate, std_array(ntimes,nbins)
                std_array must not contain NaN's,infinities or values <=0.
    smoothing = controls amount of smoothing
              = 0.0 , no smoothing
              = 1.0 ,default value"""

    if smoothing > 0:
        [ntimes,nbins]=array.shape
        bins=range(0,nbins)
        for i in range(0,ntimes):
            weights=1.0/std_array[i,:]
            weights=weights/np.mean(weights)
            sfactor= (nbins-np.sqrt(nbins))*smoothing
            temp =array[i,:].copy()
            temp[np.isnan(temp)]=0.0
            #print array[i,:]
            if 0:
                plt.figure(6000)
                plt.plot(bins,array[i,:],'r')
                plt.grid(True)
            s = UnivariateSpline(bins, array[i,:],w=weights,k=3, s=sfactor)     
            array[i,:]=s(bins)
            if 0:
                plt.figure(6001)
                plt.plot(bins,temp,bins,array[i,:],'r')
                plt.grid(True)
                plt.show()
    return array




def compute_photon_statistics(rs,inv, rs_Cxx, rs_constants ):
    """ rs=compute_photon_statistics(inv,rs_Cxx,rs_constants)
        compute errors due to photon counting statistics.
        rs_in  = structure containing var_raw_xxxx_counts arrays
              these are geo corrected counts with no background or
              baseline correction.
        rs_out
        """
  
    [ntimes, nalts] = rs.var_raw_molecular_counts.shape
    ones_array = np.ones((ntimes,nalts))
    Cmm = ones_array * np.transpose(rs_Cxx.Cmm[:nalts])
    Cmc = ones_array * np.transpose(rs_Cxx.Cmc[:nalts])
    Cam = ones_array * rs_Cxx.Cam
    if hasattr(rs_Cxx,'Cam_i2a'):
       Cam_i2a = ones_array * rs_Cxx.Cam_i2a
       Cmm_i2a = ones_array * np.transpose(rs_Cxx.Cmm_i2a[:nalts])

    
    cpol_gain = rs_constants['combined_to_cross_pol_gain_ratio']
   
    var_mol = rs.var_raw_molecular_counts
    var_comb= rs.var_raw_combined_hi_counts
    var_cp  = rs.var_raw_cross_pol_counts
    if hasattr(rs,'var_raw_molecular_i2a_counts'):
        var_mol_i2a = rs.var_raw_molecular_i2a_counts
    if hasattr(rs,'var_raw_combined_1064_counts'):
        var_combined_1064_counts = rs.var_raw_combined_1064_counts

    Scp = rs.cross_pol_counts
    Sc = rs.combined_hi_counts
    Sm = rs.molecular_counts
    
    if hasattr(rs,'var_raw_molecular_i2a_counts'):
        var_mol_i2a = rs.var_raw_molecular_i2a_counts
        Sm_i2a = rs.molecular_i2a_counts

    # variance of the scatttering ratio
    mol_f = Sm - Cam * Sc
    mol_f_sqrd = mol_f**2
    aero_f = Cmm * Sc - Cmc * Sm
   
    
    SR_i2_std = np.sqrt(
           var_comb*((Cmm/mol_f-Cam*aero_f/mol_f_sqrd)**2
                     + (cpol_gain * Scp * Cmm * Cam / mol_f_sqrd)**2) 
           +var_mol * ((aero_f/mol_f_sqrd -Cam/mol_f)**2
                    + (cpol_gain * Scp * Cmm / mol_f_sqrd)**2)
           +var_cp *(cpol_gain * Cmm/mol_f)**2)
    if hasattr(rs_Cxx,'Cam_i2a') and hasattr(rs,'molecular_i2a_counts'):
        mol_i2a_f = Sm - Cam * Sc
        mol_i2a_f_sqrd = mol_f**2
        aero_i2a_f = Cmm * Sc - Cmc * Sm

        SR_i2a_std = np.sqrt(
           var_comb*((Cmm_i2a/mol_f-Cam_i2a*aero_i2a_f/mol_i2a_f_sqrd)**2
                     + (cpol_gain * Scp * Cmm_i2a * Cam_i2a / mol_i2a_f_sqrd)**2) 
           +var_mol_i2a * ((aero_i2a_f/mol_i2a_f_sqrd -Cam_i2a/mol_i2a_f)**2
                    + (cpol_gain * Scp * Cmm_i2a / mol_i2a_f_sqrd)**2)
           +var_cp *(cpol_gain * Cmm_i2a/mol_i2a_f)**2)
        #note SR is computed from average of i2 and i2a determinations
        SR_std =0.5 * np.sqrt(SR_i2_std**2 + SR_i2a_std**2)

    else: #if no i2a channel total scattering ratio is computed from i2 channel
        SR_std = SR_i2_std
        
    #Signal to noise ratio of molecular count
    SN_mol = Sm/np.sqrt((Cam/Cmm)**2 * var_comb +1/Cmm**2 * var_mol)
   
    if hasattr(rs_Cxx,'Cam_i2a') and hasattr(rs,'molecular_i2a_counts'):
         SN_i2a_mol = Sm_i2a/np.sqrt((Cam_i2a/Cmm_i2a)**2 * var_comb +1/Cmm_i2a**2 * var_mol_i2a)
         setattr(inv,'SN_i2a_mol',SN_i2a_mol)#FIXME add to another layer somewhere
    #standard deviation of backscatter cross section
    setattr(inv,'std_beta_a_backscat',SR_std * inv.beta_r_backscat)
    #standard deviation of the backscatter ratio
    setattr(inv,'SR_std',SR_std)#FIXME add to mean layer
    #signal-to-noise ratio for beta_a_backscatter
    setattr(inv,'SN_beta_a_backscat', inv.beta_a_backscat / (SR_std * inv.beta_r_backscat))
    #signal-to-noise ratio for the molecular count profile
    setattr(inv,'SN_mol',SN_mol)#FIXME add to mean layer
    return 
        


def baseline_correction(rs,rs_cal,nalts,corr_adjusts,constants):
    """baseline_correction(rs,rs_cal,nalts,corr_adjusts)
       rs                 = structure containing channel count data
       rs_cal             = structure containing calibration baseline data
       nalts              = number of altitude bins to process
       corr_adjusts       = contains scale factors for baseline correction
       constants          = contains baseline adj factors from calvals_xxxx.txt
       baseline.data[:,0] = bin number
                    [:,1] = combined_hi baseline
                    [:,2] = combined_lo baseline
                    [:,3] = molecular baseline
                    [:,4] = cross pol baseline
                    [:,5] = molecular_i2a baseline
                    [:,6] = combined_1064 baseline
       """
    
    if rs_cal.qw_baseline.data is None:   
        if rs_cal.baseline.data.shape[0]<nalts:
            print 'baseline is too short! extending with 0'
            tmp=rs_cal.baseline.data
            rs_cal.baseline.data=np.zeros([nalts,tmp.shape[1]])
            rs_cal.baseline.data[:tmp.shape[0],:]=tmp[:,:]
        rs.cross_pol_counts -= rs_cal.baseline.data[:nalts, 4][np.newaxis,:] \
            * rs.transmitted_energy[:, np.newaxis] \
            * corr_adjusts['c_pol_baseline'] * constants['baseline_adjust'][3] 
        rs.combined_hi_counts -= rs_cal.baseline.data[:nalts, 1][np.newaxis,:] \
            * rs.transmitted_energy[:, np.newaxis] \
            * corr_adjusts['comb_hi_baseline'] * constants['baseline_adjust'][0]
        if hasattr(rs,'combined_lo_counts'):
            rs.combined_lo_counts -= rs_cal.baseline.data[:nalts, 2][np.newaxis,:] \
                 * rs.transmitted_energy[:, np.newaxis] \
                 * corr_adjusts['comb_lo_baseline'] * constants['baseline_adjust'][1]
            
        if hasattr(rs,'combined_1064_counts') and hasattr(rs,'transmitted_1064_energy'):
            if rs_cal.baseline.data.shape[1] == 7:
              rs.combined_1064_counts -= rs_cal.baseline.data[:nalts, 6][np.newaxis,:] \
                 * rs.transmitted_1064_energy[:, np.newaxis] \
                 * corr_adjusts['comb_1064_baseline'] * constants['baseline_adjust'][5]
            else:
              print "WARNING Baseline can't be applied to 1064 channel"
              
            
        if hasattr(rs,'molecular_i2a_counts'):
            if rs_cal.baseline.data.shape[1] == 6 or rs_cal.baseline.data.shape[1] == 7:
              rs.molecular_i2a_counts -= rs_cal.baseline.data[:nalts, 5][np.newaxis,:] \
                * rs.transmitted_energy[:, np.newaxis] \
                * corr_adjusts['mol_i2a_baseline'] * constants['baseline_adjust'][4]
            else:
              print "WARNING Baseline can't be applied to i2a channel"
        rs.molecular_counts -= rs_cal.baseline.data[:nalts, 3][np.newaxis,:] \
            * rs.transmitted_energy[:, np.newaxis] \
            * corr_adjusts['mol_baseline'] * constants['baseline_adjust'][2]
       
    #if both qw_baseline and baseline corrections exist along with the rotation
    #angle of the quarter waveplate, this is gvhsrl data taken with Matt Hayman's
    #modification   
    elif not rs_cal.qw_baseline.data is None and hasattr(rs,'qw_rotation_angle'):
        
        
        qwp_theta = np.mod(rs.qw_rotation_angle.astype(int),360)
        
        #add quarter wave plate rotation angle to the rs structure
        rs.qwp_theta = hau.T_Array(qwp_theta)
        
        rs.cross_pol_counts -= rs_cal.baseline.data[:nalts, 4][np.newaxis,:] \
            * rs.transmitted_energy[:, np.newaxis] * corr_adjusts['c_pol_baseline']\
            * constants['baseline_adjust'][3]\
            * rs_cal.qw_baseline.data[qwp_theta,4][:,np.newaxis]
        rs.combined_hi_counts -= rs_cal.baseline.data[:nalts, 1][np.newaxis,:] \
            * rs.transmitted_energy[:, np.newaxis] * corr_adjusts['comb_hi_baseline']\
            * constants['baseline_adjust'][0]\
            * rs_cal.qw_baseline.data[qwp_theta,1][:,np.newaxis]
        if hasattr(rs,'combined_lo_counts'):
            rs.combined_lo_counts -= rs_cal.baseline.data[:nalts, 2][np.newaxis,:] \
                 * rs.transmitted_energy[:, np.newaxis] * corr_adjusts['comb_lo_baseline']\
                 * constants['baseline_adjust'][1]\
                 * rs_cal.qw_baseline.data[qwp_theta,2][:,np.newaxis]
        rs.molecular_counts -= rs_cal.baseline.data[:nalts, 3][np.newaxis,:] \
            * rs.transmitted_energy[:, np.newaxis] * corr_adjusts['mol_baseline']\
            * constants['baseline_adjust'][2]\
            * rs_cal.qw_baseline.data[qwp_theta,3][:,np.newaxis]
    else:
        print 
        print "baseline_correction:********missing qw_baseline file**************"
        print
   
    return rs   

def diff_geometry_correction(rs,rs_cal,nalts
                         ,process_defaults,corr_adjusts):
    """diff_geometry_correction(rs,rs_cal,nalts
                             ,process_defaults,corr_adjusts):
       rs               structure containing channel counts
       rs_cal           = structure containing dif_geo data
       nalts            = number of altitudes to process
       process_defaults = processing request from process_control.json

       profiles include pre-trigger dark count bins.
       """
   
    if corr_adjusts['dif_geo_corr'] > 0 :
       # is the correction really available?
       if rs_cal.diff_geo.data  is not None:
           ones_array = np.ones_like(rs.molecular_counts)
           # combined hi
           rs.combined_hi_counts[:, :nalts] /= (ones_array[:, :nalts] \
           + rs_cal.diff_geo.data[:nalts, 1] * ones_array[:, :nalts] \
               * corr_adjusts['dif_geo_corr'])
        
           # combined lo
           if hasattr(rs,'combined_lo_counts'):
              rs.combined_lo_counts[:, :nalts] /= (ones_array[:, :nalts] \
                     + rs_cal.diff_geo.data[:nalts, 2] * ones_array[:, :nalts] \
                     * corr_adjusts['dif_geo_corr'])

           if hasattr(rs,'molecular_i2a_counts') and rs_cal.diff_geo.data.shape[1]==4:
               rs.molecular_i2a_counts[:, :nalts] /= (ones_array[:, :nalts] \
                     + rs_cal.diff_geo.data[:nalts, 5] * ones_array[:, :nalts] \
                     * corr_adjusts['dif_geo_corr'])
       else :
           print 'diff_geometry_correction: NOTE: skipping diff geo correction - none available'
    return rs

def  diff_1064_532_geometry_correction(rs,rs_cal,nalts,processing_defaults
                                      ,corr_adjusts):
    """ diff_1064_532_geometry_correction(rs,rs_cal,nalts,processing_defaults
                                      ,corr_adjusts)
        correct combined_1064_counts for differential geometry vs 532nm combined
    """
   
    if corr_adjusts['diff_1064_532_geo_corr'] > 0 and hasattr(rs,'combined_1064_counts') :
       # is the correction really available?
       if rs_cal.diff_1064_532_geo.data  is not None:
           ones_array = np.ones_like(rs.molecular_counts)
           # combined 1064 counts
           rs.combined_1064_counts[:, :nalts] /= (ones_array[:, :nalts] \
           + rs_cal.diff_1064_532_geo.data[:nalts, 1] * ones_array[:, :nalts] \
               * corr_adjusts['diff_1064_532_geo_corr'])
       else:
          print
          print 'combined_1064_counts in record but no rs_cal.diff_1064_532_geo data provided in month file'
          print     
    return rs

def diff_cp_geometry_correction(rs,rs_cal,nalts
                             ,process_defaults,corr_adjusts):
    """diff_cp_geometry_correction(rs,rs_cal,nalts
                             ,process_defaults,corr_adjusts):
       rs               structure containing channel counts
       rs_cal           = structure containing dif_geo data
       nalts            = number of altitudes to process
       process_defaults = processing request from process_control.json
       """

    # is the correction really available?
    if corr_adjusts is None:
        corr_adjusts={}

    if rs_cal.cpol_diff_geo.data is not None:
           ones_array = np.ones_like(rs.molecular_counts)
           # cross pol correction
        
           if rs_cal.diff_geo.data is not None and not corr_adjusts.get('dif_geo_corr',1.0)==0:
               rs.cross_pol_counts[:, :nalts] /= \
                   ((ones_array[:, :nalts] + rs_cal.cpol_diff_geo.data[:nalts,1] * ones_array[:, :nalts] \
                               * corr_adjusts.get('cpol_dif_geo',1.0))* \
                   (ones_array[:, :nalts] + rs_cal.diff_geo.data[:nalts, 1]  * ones_array[:, :nalts] \
                               * corr_adjusts.get('dif_geo_corr',1.0)))
           else:
               rs.cross_pol_counts[:, :nalts] /= \
                   ones_array[:, :nalts] + rs_cal.cpol_diff_geo.data[:nalts,1] * ones_array[:, :nalts] \
                               * corr_adjusts.get('cpol_dif_geo',1.0)
                    
    else:
           print 'cpol_diff_geometry_correction: NOTE: skipping cross polar geo correction - none available'

    return rs

def i2a_diff_geo_correction(rs,rs_cal,corr_adjusts):
    """diff_geometry_correction(rs,rs_cal,process_defaults,corr_adjusts):
       rs               =structure containing channel counts
       rs_cal           = structure containing dif_geo data
    
       
       """
    
    if corr_adjusts['i2a_dif_geo_corr'] > 0 and hasattr(rs,'molecular_i2a_counts'):
       # is the correction really available?
       if rs_cal.i2a_diff_geo.data is not None:
           ones_array = np.ones_like(rs.molecular_i2a_counts)
           na=rs.molecular_i2a_counts.shape[1]
           rs.molecular_i2a_counts[:, :na] /= (ones_array[:, :] \
           + rs_cal.i2a_diff_geo.data[:na, 1] * ones_array[:, :] \
               * corr_adjusts['i2a_dif_geo_corr'])
   
       else :
           print 'i2a_diff_geo_correction: NOTE: skipping i2a diff geo correction - none available'
    return rs           

def extract_dark_count(rs,rs_constants):
    """extract_dark_count(instrument,rs,rs_constants)
       extract dark count from bins prior to the laser pulse this
       does not apply corr_adjusts or correct for signal in dark count"""
    
    [ns,na] = rs.molecular_counts.shape 
    [dark_interval_end_time, laser_pulse_time, cal_pulse_end_time] = \
        rs_constants['apd_pulse_timing']
    bin_duration = rs_constants['binwidth']
    dark_interval_end_bin = int(dark_interval_end_time / bin_duration) - 1

    if rs_constants['dark_count_timing'] == 'first_bins':
        dark_interval = np.arange(dark_interval_end_bin)
    else : # last bin dark correction
        print ' '
        print '******************last bin dark correction******************'
        print ' '
        dark_interval = np.arange(na - 60, na - 20, 1)   

    #compute dark corrections for all channels 
    mol_dark_counts = \
        hau.T_Array(nanmean(rs.molecular_counts[:,
                               dark_interval], 1),summode=rs.molecular_counts.summode)
    rs.mol_dark_counts = mol_dark_counts[:, np.newaxis] 

    c_pol_dark_counts = \
        hau.T_Array(nanmean(rs.cross_pol_counts[:,
                               dark_interval], 1),summode=rs.molecular_counts.summode)
    rs.c_pol_dark_counts = c_pol_dark_counts[:, np.newaxis]

    c_hi_dark_counts = \
        hau.T_Array(nanmean(rs.combined_hi_counts[:
               ,dark_interval], 1),summode=rs.molecular_counts.summode)
    rs.c_hi_dark_counts = c_hi_dark_counts[:, np.newaxis]

    if hasattr(rs,'combined_lo_counts'):
        c_lo_dark_counts = \
            hau.T_Array(nanmean(rs.combined_lo_counts[:,
            dark_interval], 1),summode=rs.molecular_counts.summode)
        rs.c_lo_dark_counts = c_lo_dark_counts[:, np.newaxis]
    if hasattr(rs,'molecular_i2a_counts'):
        mol_i2a_dark_counts = \
            hau.T_Array(nanmean(rs.molecular_i2a_counts[:,
               dark_interval], 1),summode=rs.molecular_counts.summode)
        rs.mol_i2a_dark_counts=mol_i2a_dark_counts[:,np.newaxis]
    if hasattr(rs,'combined_wfov_counts'):
        c_wfov_dark_counts = \
            hau.T_Array(nanmean(rs.combined_wfov_counts[:,
            dark_interval], 1),summode=rs.molecular_counts.summode)
        rs.c_wfov_dark_counts = c_wfov_dark_counts[:, np.newaxis]
    if hasattr(rs,'molecular_wfov_counts'):
        m_wfov_dark_counts = \
            hau.T_Array(nanmean(rs.molecular_wfov_counts[:,
            dark_interval], 1),summode=rs.molecular_counts.summode)
        rs.m_wfov_dark_counts = m_wfov_dark_counts[:, np.newaxis]

    if hasattr(rs,'combined_1064_counts'):
        #combined_1064_dark_counts are obtained from the end of the data interval
        #because seed laser light is present in the normal dark interval
        #dark interval = [start_bin,end_bin]
        bin_limits = rs_constants['IR_dark_interval']
        if rs.combined_1064_counts.shape[1] >= bin_limits[1]:
            IR_dark_interval =np.arange(bin_limits[0],bin_limits[1])
            combined_1064_dark_counts = \
               hau.T_Array(nanmean(rs.combined_1064_counts[:,IR_dark_interval], 1)\
               ,summode=rs.molecular_counts.summode)
            rs.combined_1064_dark_counts = combined_1064_dark_counts[:, np.newaxis]
        else:
            print
            print 'WARNING---altitude range too short for end of record dark count'
            print '          1064 dark count set to zero'
            print
            combined_1064_dark_counts = np.zeros(rs.combined_1064_counts.shape[0])
            rs.combined_1064_dark_counts = hau.T_Array(combined_1064_dark_counts[:,np.newaxis])
    return

def dark_count_correction(instrument,raw,rs,rs_Cxx,corr_adjusts,processing_defaults,rs_constants):
    """dark_count_correction(rs,corr_adjusts)
    adds dark count and raw counts fields to rs and subtracts
    backgroud from all channels
       instrument   = instrument name (e.g. mf2hsrl)
       rs           = structure containing channel counts
       rs_Cxx       = structure contianing calibration coef
       corr_adjusts = scaling factors for calibration constants
       rs_constants = constants derived from calval_xxxx.txt file
    """
  
  
    if processing_defaults is None:
        pass  
    # remove signal from previous laser pulse in dark_counts
    elif processing_defaults.enabled('signal_in_dark') and corr_adjusts['signal_in_dark'] !=0 \
           and rs_constants['installation'] == 'airborne':
        print
        print 'signal in dark correction not implemented for airborne system'
        print
    elif processing_defaults.enabled('signal_in_dark') \
             and corr_adjusts['signal_in_dark'] != 0 and rs_Cxx is not None:
        signal_in_dark_count(instrument,rs,rs_Cxx,corr_adjusts,rs_constants)
    
    #apply dark counts to all channels
    [ns,na] = rs.molecular_counts.shape
    ones_array = np.ones((ns, na))
   
    rs.combined_hi_counts = rs.combined_hi_counts \
             - (raw.c_hi_dark_counts * ones_array) * corr_adjusts['comb_hi_dark_count']
    rs.molecular_counts = rs.molecular_counts \
              - (raw.mol_dark_counts * ones_array) * corr_adjusts['mol_dark_count']
    if hasattr(rs,'combined_lo_counts'):
        rs.combined_lo_counts = rs.combined_lo_counts \
             - (raw.c_lo_dark_counts * ones_array) * corr_adjusts['comb_lo_dark_count']
    if hasattr(rs,'combined_wfov_counts'):
        rs.combined_wfov_counts = rs.combined_wfov_counts \
             - (raw.c_wfov_dark_counts * ones_array) * corr_adjusts['comb_wfov_dark_count']
    if hasattr(rs,'molecular_wfov_counts'):
        rs.molecular_wfov_counts = rs.molecular_wfov_counts \
             - (raw.m_wfov_dark_counts * ones_array) * corr_adjusts['mol_wfov_dark_count']
    if hasattr(rs,'molecular_i2a_counts'):
        rs.molecular_i2a_counts = rs.molecular_i2a_counts \
        - (raw.mol_i2a_dark_counts * ones_array) * corr_adjusts['mol_i2a_dark_count']
    if hasattr(rs,'cross_pol_counts'):
        rs.cross_pol_counts = rs.cross_pol_counts \
            - (raw.c_pol_dark_counts * ones_array) * corr_adjusts['c_pol_dark_count']
    if hasattr(rs,'combined_1064_counts'):
        rs.combined_1064_counts = rs.combined_1064_counts \
            - (raw.combined_1064_dark_counts * ones_array) * corr_adjusts['comb_1064_dark_count']      
    return rs

def  signal_in_dark_count(instrument,rs,rs_Cxx,corr_adjusts,rs_constants):
        """signal_in_dark_count(instrument,rs,rs_Cxx,corr_adjusts,rs_constants)
           subtracts a correction for the molecular signal resulting from the
           previous laser pulse in the current dark count"""
       
        [nshots, nalts] = rs.molecular_counts.shape
        signal_in_dark_cor = np.zeros((nshots, 6))
        range_at_dark_interval = 1.5e8 \
            / float(rs_constants['laser_rep_rate'])
        binwidth_meters = 1.5e8 * float(rs_constants['binwidth'])
        pulse_delay_meters = float(rs_constants['apd_pulse_timing'][1]) * 1.5e8
        dark_interval_meters = float(rs_constants['apd_pulse_timing'][0]) * 1.5e8
        projection = np.cos(rs_constants['telescope_roll_angle_offset'] * np.pi
                            / 180.0)

        alt_3000 = 3000 * binwidth_meters * projection \
                + rs_constants['lidar_altitude'] - pulse_delay_meters

        # altitude when return appears in middle of dark count
        alt_dark = projection * range_at_dark_interval \
                - dark_interval_meters + rs_constants['lidar_altitude'] \
                - dark_interval_meters / 2


        
        # cal coeficient index at bin 3000
        cal_index_3000 = ((rs_Cxx.msl_altitudes - alt_3000)
                              ** 2).argmin(axis=0)

        # cal coeficient index at altitude where dark count is measured
        cal_index_dark = ((rs_Cxx.msl_altitudes - alt_dark)
                              ** 2).argmin(axis=0)

        # components [0] and [1] apply to molecular signal
        signal_in_dark_cor[:, 0] = rs_Cxx.Cmm[cal_index_dark] \
                * rs_Cxx.beta_r[cal_index_dark] / (rs_Cxx.Cmm[cal_index_3000]\
                * rs_Cxx.beta_r[cal_index_3000]) \
                * ((3000 * binwidth_meters - dark_interval_meters)
                   / range_at_dark_interval) ** 2
        signal_in_dark_cor[:, 1] = np.ones(nshots)- signal_in_dark_cor[:, 0]

        # components [2] and [3] apply to combined_hi signal
        signal_in_dark_cor[:, 2] = rs_Cxx.Cmm[cal_index_dark] \
                * rs_Cxx.beta_r[cal_index_dark] / (rs_Cxx.Cmm[cal_index_3000]\
                * rs_Cxx.beta_r[cal_index_3000]) \
                * ((3000 * binwidth_meters - dark_interval_meters)
                   / range_at_dark_interval) ** 2
        signal_in_dark_cor[:, 3] = np.ones(nshots) - signal_in_dark_cor[:, 2]
        if hasattr(rs,'molecular_i2a_counts')\
               and hasattr(rs_Cxx,'Cmm_i2a'):
            # components [4] and [1] apply to molecular signal
            signal_in_dark_cor[:, 4] = rs_Cxx.Cmm_i2a[cal_index_dark] \
                * rs_Cxx.beta_r[cal_index_dark] / (rs_Cxx.Cmm_i2a[cal_index_3000]\
                * rs_Cxx.beta_r[cal_index_3000]) \
                * ((3000 * binwidth_meters - dark_interval_meters)
                   / range_at_dark_interval) ** 2
            signal_in_dark_cor[:, 5] = np.ones(nshots) - signal_in_dark_cor[:, 4]
        
        # compute mean signal between bins 2900 and 3100
        combined_hi_3000 = \
            hau.T_Array(nanmean(rs.combined_hi_counts[:,2900:3100], axis=1))
        #combined_hi_3000 = combined_hi_3000[:,np.newaxis] 
        molecular_3000 = \
            hau.T_Array(nanmean(rs.molecular_counts[:,2900:3100], axis=1))
        

        # apply correction to molecular and combined_hi, assume no effect
        # on cross_pol or combined_lo. 
        rs.mol_dark_counts = (rs.mol_dark_counts[:,0] - molecular_3000
             * signal_in_dark_cor[:, 0]) / signal_in_dark_cor[:, 1]
        
        rs.c_hi_dark_counts = (rs.c_hi_dark_counts[:,0] - combined_hi_3000
             * signal_in_dark_cor[:, 2]) / signal_in_dark_cor[:, 3]
        rs.mol_dark_counts = rs.mol_dark_counts[:,np.newaxis]
        rs.c_hi_dark_counts = rs.c_hi_dark_counts[:,np.newaxis]
        if hasattr(rs,'molecular_i2a_counts'):
            molecular_i2a_3000 = \
                 hau.T_Array(nanmean(rs.molecular_i2a_counts[:,2900:3100], axis=1))
            rs.mol_i2a_dark_counts = (rs.mol_i2a_dark_counts[:,0] - molecular_i2a_3000
                 * signal_in_dark_cor[:, 4]) / signal_in_dark_cor[:, 5]
            rs.mol_i2a_dark_counts = rs.mol_i2a_dark_counts[:,np.newaxis]
        return     

def extract_cal_pulse(rs,constants):
    """extract_cal_pulse(rs,constants,extracopy=None)
       sums scattered light signal pulse that occurs as laser pulse
       exits the system. This is used for calibration"""
    

    [dark_interval_end_time, laser_pulse_time, cal_pulse_end_time] = \
        constants['apd_pulse_timing']
    bin_duration = constants['binwidth']
    s_bin = int(laser_pulse_time / bin_duration)  # laser pulse bin number

    cw_i2scan='enable_cw_i2scan' in constants and constants['enable_cw_i2scan']
    #if a i2_scans generated with cw seedlaser
    #use entire record after end of normal cal pulse
    if cw_i2scan:
        dark_interval_end_bin = int(cal_pulse_end_time / bin_duration)- 1
        cal_pulse_end_bin = len(rs.molecular_counts)
    else:    
        dark_interval_end_bin = int(dark_interval_end_time / bin_duration)- 1
        cal_pulse_end_bin = int(np.ceil(cal_pulse_end_time / bin_duration))
    
    if hasattr(rs,'molecular_counts'):
        rs.molecular_cal_pulse = np.sum(rs.molecular_counts[:
            ,dark_interval_end_bin:cal_pulse_end_bin], 1)
        if not cw_i2scan:
            rs.molecular_cal_pulse -= (cal_pulse_end_bin-dark_interval_end_bin) * rs.mol_dark_counts[:,0]
        rs.molecular_cal_pulse = \
            hau.T_Array(rs.molecular_cal_pulse \
                    / (1.0* rs.seeded_shots))
    if hasattr(rs,'combined_hi_counts'):    
        rs.combined_hi_cal_pulse = np.sum(rs.combined_hi_counts[:
                       , dark_interval_end_bin:cal_pulse_end_bin], 1)
        if not cw_i2scan:
            rs.combined_hi_cal_pulse -= (cal_pulse_end_bin-dark_interval_end_bin) * rs.c_hi_dark_counts[:,0]            
        rs.combined_hi_cal_pulse = \
                      hau.T_Array(rs.combined_hi_cal_pulse \
                    / (1.0* rs.seeded_shots))
    if hasattr(rs,'combined_lo_counts'):
        rs.combined_lo_cal_pulse = np.sum(rs.combined_lo_counts[:
              , dark_interval_end_bin:cal_pulse_end_bin], 1)
        if not cw_i2scan:
            rs.combined_lo_cal_pulse -= (cal_pulse_end_bin-dark_interval_end_bin) * rs.c_lo_dark_counts[:,0]           
        rs.combined_lo_cal_pulse = \
                    hau.T_Array(rs.combined_lo_cal_pulse \
                    / (1.0* rs.seeded_shots)) 
    if hasattr(rs,'molecular_i2a_counts'):
        rs.molecular_i2a_cal_pulse = np.sum(rs.molecular_i2a_counts[:
            , dark_interval_end_bin:cal_pulse_end_bin], 1)
        if not cw_i2scan:
            rs.molecular_i2a_cal_pulse -= (cal_pulse_end_bin-dark_interval_end_bin) * rs.mol_i2a_dark_counts[:,0]            
        rs.molecular_i2a_cal_pulse = \
                    hau.T_Array(rs.molecular_i2a_cal_pulse \
                    / (1.0* rs.seeded_shots))
    if hasattr(rs,'combined_1064_counts'):
        rs.combined_1064_cal_pulse = np.sum(rs.combined_1064_counts[:
            , dark_interval_end_bin:cal_pulse_end_bin], 1) # \
            #-(cal_pulse_end_bin-dark_interval_end_bin) * rs.combined_1064_dark_counts[:,0]               
        rs.combined_1064_cal_pulse = \
                    hau.T_Array(rs.combined_1064_cal_pulse \
                    / (1.0* rs.seeded_shots))
        
    return

   

def pileup_correction(rs,constants,corr_adjusts):
    """ pileup_correction(rs,constants,corr_adjusts)
    Use detector dead time and signal count rates with non-paralizing
    tp apply pileup corrections.
    rs           = structure containing raw count profiles
    constants    = system contants from calvals_xxxx.json file
    corr_adjusts = adjustment factors used to test sensitivity to
                   calibration constants
    """


    if not hasattr(rs,'molecular_counts') or rs.molecular_counts.shape[0]==0:
        return
    ones_array = np.ones(rs.molecular_counts.shape)

    #find total accumulation time for each bin  
    seeded_shots_array = rs.seeded_shots 
    bin_duration = constants['binwidth']
    seeded_shots_array2 = seeded_shots_array[:, np.newaxis]
    total_bin_time = hau.Z_Array(bin_duration) * seeded_shots_array2
    if hasattr(rs,'cross_pol_counts'):
        p_corr = constants['cross_pol_dead_time'] * corr_adjusts['c_pol_pileup'] \
            * rs.cross_pol_counts / total_bin_time
        p_corr[p_corr > .99] = .95
        rs.cross_pol_counts /= (ones_array - p_corr)
   
        
    if hasattr(rs,'molecular_counts'):
        p_corr = constants['molecular_dead_time'] * corr_adjusts['mol_pileup'] \
           * rs.molecular_counts / total_bin_time
        p_corr[p_corr > .99] = .95
        rs.molecular_counts /= (ones_array - p_corr)

    if hasattr(rs,'molecular_i2a_counts'):
        p_corr = constants['molecular_i2a_dead_time'] * corr_adjusts['mol_i2a_pileup'] \
           * rs.molecular_i2a_counts / total_bin_time
        p_corr[p_corr > .99] = .95
        rs.molecular_i2a_counts /= (ones_array - p_corr)

    
    if hasattr(rs,'combined_lo_counts'):
        p_corr = constants['combined_lo_dead_time'] * corr_adjusts['comb_lo_pileup'] \
            * rs.combined_lo_counts / total_bin_time
        p_corr[p_corr > .99] = .95
        rs.combined_lo_counts /= (ones_array - p_corr)
   
    if hasattr(rs,'combined_1064_counts'):
        p_corr = constants['combined_1064_dead_time'] * corr_adjusts['comb_1064_pileup'] \
            * rs.combined_1064_counts / total_bin_time
        p_corr[p_corr > .99] = .95
        rs.combined_1064_counts /= (ones_array - p_corr)

    if hasattr(rs,'molecular_wfov_counts'):
        p_corr = constants['molecular_wfov_dead_time'] * corr_adjusts['mol_wfov_pileup'] \
            * rs.molecular_wfov_counts / total_bin_time
        p_corr[p_corr > .99] = .95
        rs.molecular_wfov_counts /= (ones_array - p_corr)    

    if hasattr(rs,'combined_hi_counts'):    
        p_corr = constants['combined_hi_dead_time'] * corr_adjusts['comb_hi_pileup'] \
             * rs.combined_hi_counts / total_bin_time
        p_corr[p_corr > .99] = .95
        rs.combined_hi_counts /= (ones_array - p_corr)
   
    return


