#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import os.path
import traceback
from time import sleep
import numpy as np
import copy
from datetime import datetime,timedelta
#from matplotlib.dates import date2num
#from matplotlib.dates import num2date
import lg_base.core.array_utils as hau
# formatting utils for matplotlib datetimes
from lg_base.core.fmt_mpl_date import fmt_mpl_datetime, fmt_mpl_time
from hsrl.filters.filtered_extinction import filtered_extinction
import hsrl.data_stream.input_translators as it
import hsrl.data_stream.preprocess_raw as ppr
import hsrl.data_stream.preprocess_level2 as ppl2
import lg_base.core.decoratortools as nt
import lidar.lidar_utilities as lu
# Try to use the much faster nanmean from bottleneck, otherwise fall back
# to the scipy.stats version
try:
    from bottleneck import nanmean,nansum
except ImportError:
    print
    print 'No bottleneck.nanmean available! Falling back to SLOW scipy.stats.nanmean'
    print
    from scipy.stats import nanmean
    from numpy import nansum

import hsrl.calibration.calibration_utilities as cu
import processing_utilities as pu


import json

from lg_base.core.accumulation import accumulate

def generate_ave_profiles(rs_mean,qc_mask,rs_Cxx,rs_constants,processing_defaults
                             ,sel_telescope_dir,corr_adjusts,old_profiles=None):
    """create average of profiles dependent on telescope pointing direction
                     , telescope pointing may be 'all', 'zenith', or 'nadir' """
   
    try:
        [ntimes, nalts] = rs_mean.molecular_counts.shape
    except AttributeError:
        #raise RuntimeError, \
        #    "molecular counts missing"
        ntimes=0
        nalts=500
    #if rs_mean.times.shape[0] == 0:
        #raise RuntimeError, \
        #  "times missing"
    #    return old_profiles
    
    if qc_mask is not None and processing_defaults.get_value('averaged_profiles','apply_mask'):    
        #make mask array with NaN's for array elements where bit[0] of qc_mask==0 
        #all other elements of mask = 1
        mask = (np.bitwise_and(qc_mask,1)).astype('float')
        #mask is float to allow use of NaN values
        mask[mask == 0] = np.NaN
       
        print 'qc_mask applied to time averaged profiles vs altitude'
    else:
       #set mask == 1
       mask=None
       print 'qc_mask has not been applied to time averaged profiles'
    # most of the time we need to generate some/all of these profiles, 
    # because users are plotting them
    # if this is too time-consuming, we'll generate just the ones necessary for the 
    # user selected plots.

  
    indices=np.arange(rs_mean.times.shape[0])
    indices=(indices>=0)#boolean now
    if ntimes==0 and len(indices)==0:
      return old_profiles

    #average only those profiles occuring when locked to the i2 line
    #and not in cal_scan mode
    #op_mode = rs.rs_raw.op_mode.astype(int)
    #not_cal_scan = (~op_mode[:] & 32)/32
    #locked_to_i2 = ((op_mode[:] & 4)/4) 
   
    
    #bits = np.zeros((len(rs.rs_raw.times),len(bit_value)+1))
    #opmode=rs.rs_raw.op_mode.astype('uint32')
    #for i in range(len(bit_value)):
    #   bits[:,i]=(i+1)*(opmode[:] & bit_value[i])/bit_value[i]


    #this dictionary maps the long name to the shorthand used with dark_counts variables
    # the keys are used here as a list of all possible channels too
    channel_shorthand=dict(molecular_counts='mol',combined_hi_counts='c_hi',combined_lo_counts='c_lo'
                           ,combined_wfov_counts='c_wfov',combined_1064_counts='combined_1064'  
                           ,molecular_wfov_counts='m_wfov',molecular_i2a_counts='mol_i2a'
                           ,cross_pol_counts='c_pol',combined_counts='comb')

    #,combined_1064_counts='combined_1064' 
    #select desired telescope pointing direction for profiles
   
    if  (rs_constants['installation'] == 'ground' 
            or rs_constants['installation'] == 'shipborne'
            or sel_telescope_dir=='all'):
        print 'Selecting all telescope pointing directions for profiles'
        

        if processing_defaults is not None and processing_defaults.get_value('averaged_profiles','apply_mask'):
            indices[rs_mean.i2_locked  <= 0.99]=False
      
    elif sel_telescope_dir == 'zenith':
        #if telescope pointing exists, limit to selected pointing direction
        if rs_mean.telescope_pointing.shape[0]:
            print 'Selecting only zenith pointing data for profiles'
            indices[rs_mean.telescope_pointing !=1]=False
            rs_mean.telescope_pointing=rs_mean.telescope_pointing[indices]
        else:
            print 'Warning--using all shots--no telescope pointing direction in data file'
    elif sel_telescope_dir == 'nadir':
        #if telescope pointing exists, limit to selected pointing direction
        if rs_mean.telescope_pointing.shape[0]:
            print 'Selecting only nadir pointing data for profiles'
            indices[rs_mean.telescope_pointing != 0]=False
        else:
            print 'Warning--using all shots--no telescope pointing direction in data file'    
    else:         
        raise RuntimeError, \
       "Unrecognized value '%s' for telescope pointing dir--valid(all, zenith,nadir)" \
       % (sel_telescope_dir)
    #this is meant to filter out chunks from rs_mean that are calibration intervals
    if hasattr(rs_mean,'op_mode'):
        indices[np.bitwise_and(rs_mean.op_mode,16)!=0]=False
    indices=np.arange(indices.shape[0])[indices] #back to indexes
    profiles = hau.Time_Z_Group(can_append=False)
    profiles.hist=hau.Time_Z_Group()
    ft=None
    tc=0
    #tv=0
    lt=None
    if old_profiles is not None:
      #total_seeded_shots=total_seeded_shots+profiles.hist.total_seeded_shots
      ft=old_profiles.hist.ft
      tc=old_profiles.hist.tc
      #tv=old_profiles.hist.tv
      lt=old_profiles.hist.lt
    if len(indices)>0:
      if ft is None:
        ft=rs_mean.times[indices][0]          
      lt=rs_mean.times[indices][-1]+timedelta(seconds=rs_mean.delta_t[indices][-1] if not np.isnan(rs_mean.delta_t[indices][-1]) else 0)                  
      for x in indices:
          if rs_mean.times[x] is not None:
              tc=tc+1;
              #tv=tv+(rs_mean.times[x]-ft).total_seconds()
    if tc>0:
      profiles.times = hau.T_Array([ft])#+timedelta(seconds=tv/tc), ])
      profiles.start = ft
      profiles.width = lt-ft
      profiles.delta_t = hau.T_Array([profiles.width.total_seconds()])
    else:
      profiles.times=hau.T_Array([])
      profiles.start = ft
      profiles.width = timedelta(seconds=0)
      profiles.delta_t = hau.T_Array([])
    profiles.start_time = ft
    profiles.end_time = lt
    profiles.hist.ft=ft
    profiles.hist.tc=tc
    #profiles.hist.tv=tv
    profiles.hist.lt=lt
    #profiles.hist.total_seeded_shots=total_seeded_shots
    if rs_mean is not None and hasattr(rs_mean,'msl_altitudes'):
      profiles.msl_altitudes = rs_mean.msl_altitudes
    elif hasattr(old_profiles,'msl_altitudes'):
      profiles.msl_altitudes=old_profiles.msl_altitudes

    if rs_mean is not None and hasattr(rs_mean,'geo_corr'):
      profiles.geo_corr = rs_mean.geo_corr
    elif hasattr(old_profiles,'geo_corr'):
      profiles.geo_corr = old_profiles.geo_corr
      
    if rs_constants['installation'] == 'airborne' and len(indices)>0:
      if hasattr(rs_mean,'GPS_MSL_Alt'):
        profiles.min_GPS_MSL_Alt=np.min(rs_mean.GPS_MSL_Alt[indices])
        profiles.mean_GPS_MSL_Alt = nanmean(rs_mean.GPS_MSL_Alt[indices])
        profiles.max_GPS_MSL_Alt=np.max(rs_mean.GPS_MSL_Alt[indices])
        profiles.telescope_pointing =np.zeros(1)
      if hasattr(rs_mean,'telescope_pointing'):
        if (rs_mean.telescope_pointing>.95).all():
            profiles.telescope_pointing[0] = 1
        elif (rs_mean.telescope_pointing<.05).all():
            profiles.telescope_pointing[0] = 0
        else:
            profiles.telescope_pointing[0] = np.NaN
   
   
    accumulate(profiles,old_profiles,rs_mean,indices,'seeded_shots',pref='mean_',filler=hau.T_Array([0]))
    total_seeded_shots=profiles.mean_seeded_shots
    profiles.seeded_shots=total_seeded_shots.copy()
    print 'Total seeded shots for profile =',total_seeded_shots
    accumulate(profiles,old_profiles,rs_mean,indices,'transmitted_1064_energy',total_seeded_shots)
    accumulate(profiles,old_profiles,rs_mean,indices,'transmitted_energy',total_seeded_shots)
    # create TZ_Array with time dimension of '1', so hsrl_inversion doesn't choke
   
    for chan in channel_shorthand.keys():
        accumulate(profiles,old_profiles,rs_mean,indices,chan,total_seeded_shots,mask)

    #compute inverted profiles from mean count profiles
    if rs_Cxx is not None and hasattr(profiles,'molecular_counts'):
       profiles.inv = cu.hsrl_inversion(profiles
               , rs_Cxx, rs_constants,corr_adjusts,processing_defaults)
      
    elif hasattr(old_profiles,'inv'):
      profiles.inv=old_profiles.inv
    
    #adds raw_color_ratio to profiles.inv
    if hasattr(profiles,'inv') and hasattr(profiles,'combined_counts') \
               and hasattr(profiles,'combined_1064_counts'):
        if 0:
          print 'profiles'
          print dir(profiles)
          print 'inv'
          print dir(profiles.inv)
        
          import matplotlib.pylab as plt
          plt.figure(3333)
          plt.plot(profiles.combined_counts[0,:],profiles.inv.msl_altitudes,'r'
            ,profiles.cross_pol_counts[0,:],profiles.inv.msl_altitudes,'g'
            ,profiles.combined_counts[0,:]
            +rs_constants['combined_to_cross_pol_gain_ratio']*profiles.cross_pol_counts[0,:],'c'   
            ,profiles.combined_1064_counts[0,:],profiles.inv.msl_altitudes,'k')
          plt.grid(True)
          ax=plt.gca()
          ax.set_xscale('log')

          
        profiles.inv.raw_color_ratio=cu.compute_raw_color_ratio(profiles,rs_Cxx,rs_constants,corr_adjusts)

        if 0:
          plt.figure(3334)
          plt.plot(profiles.inv.raw_color_ratio[0,:],profiles.inv.msl_altitudes,'c')
          plt.grid(True)
          ax=plt.gca()
          ax.set_xscale('log')

    #generate klett profiles if requested      
    if processing_defaults is not None and processing_defaults.get_value('klett','enable'):
        ref_altitude = processing_defaults.get_value('klett','ref_altitude') * 1000.0
        if ref_altitude < profiles.inv.msl_altitudes[0] \
                        or ref_altitude > profiles.inv.msl_altitudes[-1] :
            print
            print
            print 'WARNING---klett ref altitutde=',ref_altitude, ' is not in requested altitudes'
            print 'no klett profile retrieval attempted '
            print
        else:    
            if hasattr(profiles,'combined_1064_counts'): 
               profiles.inv.beta_a_1064_backscat_klett = lu.compute_klett_backscatter(
                   profiles.combined_1064_counts
                   ,profiles.inv.beta_r_backscat/16.0
                   ,profiles.inv.msl_altitudes
                   ,processing_defaults.get_value('klett','lidar_ratio_532')
                   ,ref_altitude)
            if hasattr(profiles,'combined_counts'):
                profiles.inv.beta_a_532_backscat_klett = lu.compute_klett_backscatter(
                  profiles.combined_counts
                  ,profiles.inv.beta_r_backscat
                  ,profiles.msl_altitudes
                  ,processing_defaults.get_value('klett','lidar_ratio_532')
                  ,ref_altitude)
           
    for chan in channel_shorthand.keys():
          accumulate(profiles,old_profiles,rs_mean,indices,'var_raw_'+chan
               ,total_seeded_shots,mask)
          accumulate(profiles,old_profiles,rs_mean,indices,'raw_'+chan
               ,total_seeded_shots,mask)


    if processing_defaults is not None and processing_defaults.get_value('compute_stats','enable'):
     
        #kludge
        pf = hau.Time_Z_Group()
        for chan in channel_shorthand.keys():
          if hasattr(profiles,'sum_var_raw_'+chan):
            setattr(pf,chan,getattr(profiles,'sum_'+chan))
            setattr(pf,'var_raw_'+chan,getattr(profiles,'sum_var_raw_'+chan))
        
        if rs_Cxx is not None and hasattr(profiles,'inv'):
            pu.compute_photon_statistics(pf,profiles.inv,rs_Cxx,rs_constants)
                     
    
    if rs_Cxx is not None and hasattr(profiles,'inv'):    
      [profiles.inv.optical_depth, profiles.inv.optical_depth_aerosol
             , profiles.inv.mol_norm_index,profiles.inv.mol_ref_aod] = \
                 lu.compute_optical_depth(profiles.inv.Nm \
                 ,profiles.inv.beta_r_backscat\
                 ,profiles.msl_altitudes\
                 ,processing_defaults\
                 ,rs_constants
                 ,profiles.inv.telescope_pointing if hasattr(profiles.inv,'telescope_pointing') else None)
      
      #add 1064 aerosol backscatter and color ratio to profiles
      if hasattr(profiles,'combined_counts') and hasattr(profiles,'combined_1064_counts'): 
          cu.compute_1064_aerosol_backscatter(profiles,profiles.inv,processing_defaults
               ,rs_constants,corr_adjusts) 
          cu.compute_color_ratio(profiles.inv)
          

      if processing_defaults.get_value('extinction_processing','filter_type') == 'savitzky_golay':
        od_threshhold = processing_defaults.get_value('extinction_processing','od_threshhold')
        z_window_width = processing_defaults.get_value('extinction_processing','alt_window_length')
        order = processing_defaults.get_value('extinction_processing','polynomial_order')
        
        min_filter_alt = processing_defaults.get_value('extinction_processing','min_alt')
        if min_filter_alt < profiles.msl_altitudes[0] :
            min_filter_alt = profiles.msl_altitudes[0]
        
        adaptive = processing_defaults.get_value('extinction_processing','adaptive')
      
        t_window_width =0.0
        if profiles.inv.times.size==0:
          pass
        elif hasattr(rs_mean,'telescope_pointing'):
            profiles.inv = filtered_extinction(profiles.inv
             ,profiles.msl_altitudes,min_filter_alt,od_threshhold
             ,t_window_width, z_window_width,order,adaptive,rs_mean.telescope_pointing)       
        else:
           profiles.inv = filtered_extinction(profiles.inv
             ,profiles.msl_altitudes,min_filter_alt,od_threshhold
             ,t_window_width,z_window_width,order,adaptive)
        
      else:
        bin_delta = processing_defaults.get_value('extinction_processing','bin_delta')
        bin_delta=int(bin_delta)
        pts_to_ave= processing_defaults.get_value('extinction_processing','ext_pts_to_ave')

        if hasattr(rs_mean,'telescope_pointing'):
            profiles.inv = pu.compute_extinction(profiles.inv
              ,profiles.msl_altitudes, bin_delta,pts_to_ave
              ,rs_mean.telescope_pointing)
        else:
            profiles.inv = pu.compute_extinction(profiles.inv
              ,profiles.msl_altitudes, bin_delta,pts_to_ave)
    
       
    # raw profiles--ave photons/bin/laser_pulse
    #for chan in channel_shorthand.keys():
    #    accumulate(profiles,old_profiles,rs_raw,raw_indices,chan,raw_total_seeded_shots,pref='raw_')

    if 0:
        import matplotlib.pylab as plt
        plt.figure(898989)
        plt.plot(np.nanmean(rs_mean.combined_1064_counts,0),rs_mean.msl_altitudes,'k'
                 ,np.nanmean(rs_mean.combined_counts,0),rs_mean.msl_altitudes,'r')
        ax = plt.gca()
        ax.set_xscale('log')
   
    # dark corrected raw profiles
    for chan,shorthand in channel_shorthand.items():  
        dcchan=shorthand+'_dark_counts'
        correc='dc_'+chan
        source='raw_'+chan
        if accumulate(profiles,old_profiles,rs_mean,indices,dcchan,total_seeded_shots,extravars=[correc]):
            if hasattr(profiles,source):
                #print 'Applying dark count from raw frame to raw counts*******'
                setattr(profiles,correc,getattr(profiles,source)-getattr(profiles,dcchan))                    
            elif hasattr(old_profiles,correc):
                setattr(profiles,correc,getattr(old_profiles,correc))
                print 'Copying corrected counts because source doesn\'t exist???? WARNING ***'
                if hasattr(old_profiles,source):
                    setattr(profiles,source,getattr(old_profiles,source))
                else:
                    raise RuntimeError('Something is wrong with input channels. BUG!')
   
    if processing_defaults is not None:
      profiles.mol_norm_alt = processing_defaults.get_value('mol_norm_alt','meters')     
    return profiles


