import numpy as np
import copy
import os
import lidar.lidar_utilities as lu
from datetime import datetime, timedelta
import hsrl.filters.savitzky_golay as sg
import matplotlib.pylab as plt

def calibration_path_for(instrument,date,process_defaults,alternate_cal_dir=None):
    import lg_dpl_toolbox.core.archival as hru
    basedir=hru.get_path_to_data(instrument,date)
    if alternate_cal_dir==None:
        if process_defaults!=None:
            alternate_cal_dir = process_defaults.get_value('alternate_cal_dir','full_dir_path') 
            if alternate_cal_dir == 'None':
                alternate_cal_dir = None
        if alternate_cal_dir==None:
            alternate_cal_dir=os.getenv('HSRL_ALTERNATE_CAL_DIR',None)
    if alternate_cal_dir!=None:
        basedir=os.path.join(alternate_cal_dir,instrument)
    ret=os.path.join(basedir,'%04i' % date.year,'%02i' % date.month,'%02i' % date.day,"calibration")
    print 'path ',ret
    try:
      os.makedirs(ret)
      os.chmod(ret,0775)
    except:
      pass
    return ret

def overlap_correction(altitudes,lidar_ratio,fit_range,dc_counts,beta_r,beta_a_backscat,consts):
    """
       overlap_correction(delta_r,lidar_ratio,fit_range,counts,beta_r,beta_a_backscat,consts)
       dc_counts = dark and pileup corrected counts profile, sbin pre-laser bins removed
                   dc_counts must be supplied with raw bin spacing--1D vector
       beta_r = molecular extinction cross section profile,  corresponding to altitudes--1D vector
       beta_a_backscat = aerosol backscatter cross section profile, corresponding to altitudes--1D vector
       fit_range = range at which to normalize expected signal
    """
 
    delta_r = consts['binwidth'] * 1.5e8 #dt * 1/2 speed of light
    lidar_altitude = consts['lidar_altitude']
    fit_index = np.int(fit_range/delta_r)
    
    #put variables on range scale with bin spacing = delta_r
    ranges = delta_r * np.arange(dc_counts.shape[0])
    beta_a = np.interp(ranges,altitudes - lidar_altitude,beta_a_backscat)
    beta_m = np.interp(ranges,altitudes - lidar_altitude,beta_r)
    counts = dc_counts.copy()
  
    beta_a[np.isnan(beta_a)] = 0.0  
    optical_depth_aerosol = delta_r * np.cumsum(beta_a * lidar_ratio)
    optical_depth_molecular = delta_r * np.cumsum(beta_m)
    expected_signal =beta_m * np.exp(-2.0 *(optical_depth_aerosol\
                            +optical_depth_molecular))/ranges**2
    
    overlap_corr = np.ones_like(dc_counts)
    temp_counts = expected_signal  *  counts[fit_index] \
             / expected_signal[fit_index]
    overlap_corr = temp_counts / counts
    return overlap_corr, fit_index   #first bin at laser pulse




def user_input(channel_name):
           print
           print 'input estimated lidar ratio for ' + channel_name
           print 'c to clear figure'
           print 'CR if done'
           lidar_ratio = raw_input('lidar ratio = ?')
           if len(lidar_ratio)==0:
               #lidar_ratio = None
               #fit_range = None
               return None , None
               #temp_counts = counts.copy()
               #temp_counts[:index_z_low] *= geo_corr[:index_z_low]
               #return temp_counts,geo_corr

           elif lidar_ratio == 'c':
               f=plt.figure(6000)
               f.clear()
               lidar_ratio = raw_input('lidar ratio = ?')
               lidar_ratio = np.float(lidar_ratio)
           else:
               lidar_ratio = np.float(lidar_ratio) 
           fit_range = np.float(raw_input('replace with expected signal below (km) = ?')) * 1000.0
           return lidar_ratio,fit_range
       
def corr_wfov_overlap(counts,wfov_name,rs_p,ranges,consts):
       """correct overlap of wfov signal after dark correction and removal of pre-lase bins"""       
       while 1:
          
           lidar_ratio,fit_range = user_input(wfov_name)
           #user input returns none when finished
           if lidar_ratio == None:
               temp_counts = counts.copy()
               temp_counts[:index_z_low] *= geo_corr[:index_z_low]
               return temp_counts,geo_corr
       
               
           beta_r = (rs_p.inv.beta_r_355 + rs_p.inv.beta_r_387)/2.0
          
           geo_corr,index_z_low = overlap_correction(rs_p.inv.altitudes,lidar_ratio,fit_range,counts,beta_r,rs_p.inv.beta_a_backscat[0,:],consts)
              
           if 1:
               plt.ion()
               plt.figure(6000)
               plt.plot(counts,ranges,'b',counts[:index_z_low]*geo_corr[:index_z_low],ranges[:index_z_low],'r')
               plt.grid(True)
               plt.title(wfov_name)
               plt.ylabel('Range (m)')
               plt.xlabel('counts')
               ax = plt.gca()
               ax.set_xscale('log')

               plt.ion()
               plt.figure(6001)
               plt.plot(ranges,geo_corr,'b',ranges[:index_z_low],geo_corr[:index_z_low],'r')
               plt.grid(True)
               plt.title(wfov_name)
               plt.ylabel('Range (m)')
               plt.legend('all','fit_range')
               plt.xlabel('geo_corr--no rsqd')
               ax = plt.gca()
               ax.set_xscale('log')
              
         
           
      

def make_wfov_geofile(chan_sel_dict,rs_in,consts,process_control,corr_adjusts=None):
    """
       make_wfov_geofile(chan_sel_dict,rs,sbin,consts)
       chan_sel_dict = dict(chanel_names = "dark_count_names" ......) 
                       correct channel_names with corresponding dark_counts
       rs       = structure containing channel counts to be corrected
    """
   
    if not hasattr(rs_in, 'raman_rawprofile') \
          or not hasattr(rs_in,'raman_inv'):
        print
        print 'Missing structure'
        print 'Raman wfov geofile requires both "raman_rawprofile"'
        print 'and "raman_inv"'
        print
        return
  
    
    rs = copy.deepcopy(rs_in.raman_rawprofile)
    rs_p = copy.deepcopy(rs_in.raman_profile)

    select = dict(elastic_counts_low='elastic_counts_low'
                ,nitrogen_counts_low ='nitrogen_counts_low')
    lu.clear_first_bins(select,rs,process_control,consts)
    print 'rs',dir(rs)
    print
    print 'rs_p',dir(rs_p)


    #chanels and associated dark counts in rs
    chan_sel_dict = dict(sum_nitrogen_counts_high = 'nitrogen_high_dark_counts'
                  ,sum_nitrogen_counts_low = 'nitrogen_low_dark_counts')
    start_index=0
    sbin = consts['data_bin_containing_laser_pulse']
    end_index = sbin - 2
    lu.dark_count_correction_from_signal(chan_sel_dict,rs,start_index,end_index)
    
    channel_name = 'sum_nitrogen_counts_low'
    wfov_counts = getattr(rs,channel_name).astype('float')
    wfov_counts = wfov_counts[0,:]
    
    delta_r = np.float(consts['binwidth']) * 1.5e8  #binwidth in seconds * 1/2 * speed of light
    bin_vec = np.arange(wfov_counts.shape[0])
    ranges = delta_r * bin_vec
    print 'ranges',ranges.shape   

    #strip pre-lase bins
    wfov_counts[:-sbin] = wfov_counts[sbin:]

    wfov_counts,wfovgeo = corr_wfov_overlap(wfov_counts,channel_name,rs_p,ranges,consts)


 
    wfovgeo *=  ranges**2 / 1.0e6
    gcorr = wfovgeo[:,np.newaxis]
    channel_name_list =['nitrogen_low']
    if 1:
        import matplotlib.pylab as plt
        print ranges
        print 'ranges',ranges[sbin:].shape
        print 'n2_counts',rs.sum_nitrogen_counts_low[0,sbin:].shape
        plt.figure(88888)
        plt.plot(ranges[:-sbin],rs.sum_nitrogen_counts_low[0,sbin:] * wfovgeo[:-sbin])

    write_geofile('raman',rs.times,ranges,gcorr,'test',channel_names=channel_name_list)
        
    return

def write_geofile(instrument,times,bin_ranges,geo_corr,geo_type,channel_names = None
                  ,assumptions = None, process_defaults=None):
    """write_geofile(instrument,times,bin_ranges,geo_corr,type,lidar_ratio = None)
        write geofile in "..../year/mm/dd/calibration/" directory
        instrument = instrument name (eg. bagohsrl)
        times = times used in generating geo correction data
        bin_ranges = vector of bin ranges (m)
        geo_corr = vector of geo corrections at specified ranges
        geo_type = 'std'--computed from ratio of mol return and Rayleigh profile
                 = 'wfov'--computed from ratio of wfov and molecular return 
                 = 'composite'--uses combination of 'std' and 'wfov'
        assumptions = dict['lidar_ratio'] = lidar_ratio
                          ['lidar_ratio_alt_lmt'] = apply lidar_ratio below this alt
                          ['lidar_ratio_bin_lmt'] = bin number for above alt limit
                          ['linear_extrap_bin_lmt'] = do linear extrap above this bin
                          ['linear_extrap_rate']    = slope of linear extrapolation
 
      """


    start_time_str = times[0].strftime('%d-%b-%y %H:%M')
    end_time_str   = times[-1].strftime('%H:%M')

    #get path to directory


    #dir_path=calibration_path_for(instrument,times[0],process_defaults)
    dir_path = '/home/eloranta/typedarray_hsrl_python/raman/config'

    #define start time for use in filename
    str=times[0].strftime("%Y%m%dT%H%M")

    filename=os.path.join(dir_path,'geofile_'+str+'.geo')
    print filename
    fileid=open(filename,'w')
    os.chmod(filename,0664)
    print >>fileid, '#'+geo_type+' geo_corr, data from %s -->%s UTC' %(start_time_str,end_time_str)
    print >>fileid, '# file created  %s' %datetime.now()
    if assumptions:
        for item in assumptions:
           print item
        print >>fileid, '#lidar ratio assumed = %2.1f below %2.1f km (bin_number=%i)'\
            %(assumptions['lidar_ratio'],assumptions['lidar_ratio_alt_lmt']\
               ,assumptions['lidar_ratio_bin_lmt'])
        if assumptions.has_key('linear_extrap_above_bin'):
            print >>fileid, '#linear extrapolation at bins greater than %i '\
             %(assumptions['linear_extrap_above_bin'])
            print >>fileid, '#linear_extrapolation_rate = %6.4f'%(assumptions['linear_extrap_rate'])    
    print '#'        

    nbins = geo_corr.shape[0]
    nchan = geo_corr.shape[1]
    print geo_corr.shape
    if not channel_names == None:
        print >>fileid, '#Range  ',
        for i in range(nchan):
            print i, channel_names[i]
            print >>fileid,  ' %s ' %(channel_names[i]),
        print >>fileid, ' '
    else:
        print >>fileid, '#range   geometry correction'

        
    for i in range(nbins):
        print >>fileid, '%11.4e '  %(bin_ranges[i]),
        for k in range(nchan):
             print >>fileid, '%12.6e ' %(geo_corr[i,k]), 
        print >>fileid, ' '
    
    
    fileid.close()

    print "\nnew geofile created with '" +geo_type+"'"
    print filename
      
    return
