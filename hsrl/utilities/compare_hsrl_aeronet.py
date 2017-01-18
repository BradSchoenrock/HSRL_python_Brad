from matplotlib.dates import date2num, num2date
from datetime import datetime
from pycdf import CDF
import numpy as np
from read_aeronet import  read_aeronet
import hsrl.graphics.hsrl_display as hd
import matplotlib.pylab as plt
import lg_base.graphics.graphics_toolkit as gt

import matplotlib.pyplot as plt

def compare_hsrl_aeronet(rti,aeronet_file,cmp_alt=None,plot_limit=None):
    """compare_hsrl_aeronet(rti,aeronet_file,cmp_alt=None,plot_limit=None):
       plot aeronet optical depth on existing od_vs_time plot
       aeronet_file = name of file downloaded from aeronet web site for requested time
       cmp_alt = plot hsrl vs aeronet using hsrl data from this altitude 
       plot_limit  = limit hsrl vs aeronet plot to this OD 
    """
 
    if cmp_alt == None or plot_limit == None:
       print "INPUT ERROR--usage cmp_aeronet(rti, 'aeronet_file_name',comparison_alt(km),od_plot_upper_limit)"
       return                        
    times = rti.rs.rs_inv.times
    print 'read aeronet file = ' + aeronet_file
    [AOT500,AOT532,AOT667,AOT370,aeronet_times] = read_aeronet(aeronet_file,times)
    AOT532_trimed =AOT532[aeronet_times <=rti.rs.rs_inv.times[-1]]
    AOT370_trimed = AOT370[aeronet_times <= rti.rs.rs_inv.times[-1]]
    aeronet_times_trimed = aeronet_times[aeronet_times <= rti.rs.rs_inv.times[-1]]
    """ 
    #estimate hsrl unmeasured aerosol optical depth below molecular normalization altitude
    #by extrapolating optical depth slope measured at ref_alt to ground level
    norm_index = rti.rs.rs_inv.mol_norm_index
    #calculate the number bins between ground and ref altitude--allow fractional bins
    bin_range_to_ref_alt = (rti.rs.rs_inv.msl_altitudes[norm_index] - rti.rs.rs_mean.lidar_altitude)\
                           /(rti.rs.rs_inv.msl_altitudes[1]-rti.rs.rs_inv.msl_altitudes[0])
    print 'bin_range_to_ref_alt ',bin_range_to_ref_alt
    print rti.rs.rs_inv.msl_altitudes[norm_index]
    print rti.rs.rs_mean.lidar_altitude
    print norm_index
    #overlap_aod = #_bins_grd_to_ref_alt * OD_of_one_bin_at_ref_alt
    overlap_aod = bin_range_to_ref_alt\
           * (rti.rs.rs_inv.optical_depth_aerosol[:,norm_index+1] \
             -rti.rs.rs_inv.optical_depth_aerosol[:,norm_index])
    """
    norm_index = rti.rs.rs_inv.mol_norm_index
    overlap_aod = rti.rs.rs_inv.mol_ref_aod
    #interpolate overlap_aod computed at hsrl times to aeronet times
    overlap_aod_at_aeronet_times = np.zeros(len(aeronet_times_trimed))
    for i in range(len(aeronet_times_trimed)):
        #index into times less <= aeronet_times
        index = len(times[times<= aeronet_times_trimed[i]])
        if index == 0 or index == len(times) or times[index] == aeronet_times_trimed[i]:
            overlap_aod_at_aeronet_times[i] = overlap_aod[index]
        else:
            overlap_aod_at_aeronet_times[i] = (overlap_aod[index-1] + overlap_aod[index])/2.0

    print 'aeronet start, end times = ', aeronet_times[0],aeronet_times[-1]
    print 'overlap_aods'
    print overlap_aod_at_aeronet_times
    #plot on existing od vs time plot
    gt.plot_vs_time('od_vs_time'
                 ,'ahsrl'
                 ,aeronet_times_trimed
                 ,[AOT532_trimed - overlap_aod_at_aeronet_times]
                 ,['r']
                 ,[0]
                 ,None
                 ,'upper left'    
                 ,'Optical depth above ' + str(int(rti.rs.rs_inv.msl_altitudes[norm_index])/1000.0)+'(km)'
                 ,None
                 ,'Optical depth vs time'
                 ,None   
                 ,rti.display_defaults #grabs the display defaults stored from the artist
                 ,rti.figs)
   
    #plot hsrl od vis aeronet od
    bs_threshold = 1e-4
    i=0
    #    alt_index = len(rti.rs.rs_inv.msl_altitudes[self.rs.rs_inv.msl_altitudes<= cmp_alt*1000.0])
    #    altitude = rti.rs.rs_inv.msl_altitudes[alt_index]

    j=0
    alt_index = len(rti.rs.rs_inv.msl_altitudes[rti.rs.rs_inv.msl_altitudes<= cmp_alt*1000.0])
    altitude = rti.rs.rs_inv.msl_altitudes[alt_index]

    hsrl_ods = np.NaN * np.ones(aeronet_times_trimed.shape)
    aeronet_ods = np.NaN * np.ones(aeronet_times_trimed.shape)
    hsrl_ods_cloud = np.NaN * np.ones(aeronet_times_trimed.shape)
    aeronet_ods_cloud = np.NaN * np.ones(aeronet_times_trimed.shape)
    for a_time in aeronet_times_trimed:
          time_index = len(times[times <= a_time])
          h_od = rti.rs.rs_inv.optical_depth_aerosol[time_index,alt_index]
          a_od = AOT532_trimed[i]
          if h_od < plot_limit and a_od < plot_limit \
                and rti.rs.rs_inv.qc_mask[time_index,alt_index]:
              hsrl_ods[i] = h_od
              aeronet_ods[i] = a_od
              print i,h_od, a_od
              i = i + 1
              #find clouds in hsrl data
              if any(rti.rs.rs_inv.beta_a_backscat[time_index,2:] > bs_threshold):
                 hsrl_ods_cloud[j] = h_od
                 aeronet_ods_cloud[j] = a_od
                 j = j + 1       
    one_to_one_line = np.array([0,plot_limit])
    colors=['k','c','r','g']
    linetype=['*','-','*','*']
    linewidth=[2,4,4,4]
    
    gt.plot_xy('od_hsrl_vs_aeronet_raw'  #plot name
                ,' '
                ,times
                ,[one_to_one_line,aeronet_ods,aeronet_ods_cloud]
                ,[one_to_one_line,hsrl_ods,hsrl_ods_cloud]
                ,['k','r','c']
                ,['*','o','o']
                ,[4,4,4,4]
                ,['-','o','o','o']
                ,linewidth
                ,None
                ,'upper right'
                ,'Aeronet Optical depth'
                ,None
                ,'HSRL 532 OD ' +str(int(altitude)) +'(m)'
                ,None
                ,'HSRL vs aeronet OD'
                ,['no near range dead zone correction']
                ,[plot_limit *0.5]
                ,[plot_limit *0.95]
                ,[0]
                ,rti.display_defaults
                ,rti.figs)
    
    gt.plot_xy('od_hsrl_vs_aeronet_corrected'  #plot name
                ,' '
                ,times
                ,[one_to_one_line,aeronet_ods,aeronet_ods_cloud]
                ,[one_to_one_line,overlap_aod_at_aeronet_times + hsrl_ods,hsrl_ods_cloud]
                ,colors
                ,['*','o','o']
                ,[4,4,4,4]
                ,['-','o','o','o']
                ,linewidth
                ,None
                ,'upper right'
                ,'Aeronet Optical depth'
                ,None
                ,'HSRL 532 OD 0-->'+ str(int(altitude)/1000.0) +'(km)'
                ,None
                ,'HSRL vs AERONET optical depth'
                ,['includes near range dead zone correction']
                ,[plot_limit *0.5]
                ,[plot_limit *0.95]
                ,[0]
                ,rti.display_defaults
                ,rti.figs)
    #compute standard deviation of aeronet and hsrl ods less than plot limit
    temp_x = aeronet_ods[aeronet_ods <= plot_limit]
    temp_y = overlap_aod_at_aeronet_times[aeronet_ods <= plot_limit] + hsrl_ods[aeronet_ods <=plot_limit]
    temp_x = temp_x[temp_y <= plot_limit]
    temp_y = temp_y[temp_y <= plot_limit]
    std    =  np.nanmean(np.sqrt((temp_x -temp_y)**2))
    print
    print  'std deviation ODs = ', std
    print    

    
    temp_x = aeronet_ods[np.isnan(overlap_aod_at_aeronet_times+hsrl_ods) == 0]
    temp_y = overlap_aod_at_aeronet_times[np.isnan(overlap_aod_at_aeronet_times+hsrl_ods)==0]\
               + hsrl_ods[np.isnan(overlap_aod_at_aeronet_times+hsrl_ods) == 0]
    temp_y = temp_y[np.isnan(temp_x)==0]
   
    pc = np.polyfit(temp_x, temp_y , 1)
    x_od = np.array([0,plot_limit])
    fit_od = np.polyval(pc,x_od)
    gt.plot_xy('od_hsrl_vs_aeronet_with_fit'  #plot name
                ,' '
                ,times
                ,[one_to_one_line,aeronet_ods,x_od]
                ,[one_to_one_line,overlap_aod_at_aeronet_times + hsrl_ods,fit_od]
                ,colors
                ,['*','o','o']
                ,[4,4,4,4]
                ,['-','o','-']
                ,linewidth
                ,None
                ,'upper right'
                ,'Aeronet Optical depth'
                ,None
                ,'HSRL 532 OD 0-->'+ str(int(altitude)/1000.0) +'(km)'
                ,None
                ,'HSRL vs AERONET optical depth'
                ,['std = ' + str(std)] #text_str
                ,[plot_limit * 0.5] #text_position_x
                ,[plot_limit * 0.95]  #text_position_y
                ,[0] #text_angle
                ,rti.display_defaults
                ,rti.figs)

    #create a new od vs time plot with hsrl values at comparison altitudes
    alt_index = len(rti.rs.rs_inv.msl_altitudes[rti.rs.rs_inv.msl_altitudes<= cmp_alt*1000.0])
    gt.plot_vs_time('hsrl_and_aeronet_od_vs_time'
                 ,'ahsrl'
                 ,times
                 ,[rti.rs.rs_inv.optical_depth_aerosol[:,alt_index] + overlap_aod]
                 ,['b']
                 ,[2]
                 ,None
                 ,'upper left'    
                 ,'Optical depth'
                 ,None
                 ,'Optical depth vs time'
                 ,None   
                 ,rti.display_defaults #grabs the display defaults stored from the artist
                 ,rti.figs)
    #overplot with aeronet values
    gt.plot_vs_time('hsrl_and_aeronet_od_vs_time'
                 ,'ahsrl'
                 ,aeronet_times_trimed
                 ,[AOT532_trimed]
                 ,['r']
                 ,[0]
                 ,None
                 ,'upper left'    
                 ,'Optical depth'
                 ,None
                 ,'Optical depth vs time'
                 ,None   
                 ,rti.display_defaults #grabs the display defaults stored from the artist
                 ,rti.figs)
    
   
    
    if 'rl_od_vs_time' in rti.figs:
        #don't bother with this plot if no raman optical depth plot
        try:
           gt.plot_vs_time('rl_od_vs_time'
                 ,'raman'
                 ,aeronet_times_trimed
                 ,[AOT370_trimed]
                 ,['m']
                 ,[0]
                 ,['370nm']
                 ,'upper left'    
                 ,'Optical depth'
                 ,None
                 ,'Optical depth vs time'
                 ,None   
                 ,rti.display_defaults #grabs the display defaults stored from the artist
                 ,rti.figs)
        except:
            print 'no raman optical depth vs time plot'
         
    rti.figs.showall()
    """
    #plt.show()
    
    
if __name__ == '__main__':
    hsrl_file = '/home/eloranta/processed_lidar_data/'
    hsrl_file = hsrl_file + 'bagohsrl_20140722T1200_20140722T2359_22.nc'
    aeronet_file ='/home/eloranta/frappe/'
    aeronet_file = aeronet_file + '140722_140722_BSRN_BAO_Boulder.lev10'
    compare_hsrl_aeronet(hsrl_file,aeronet_file,'22-jul-14 12:00','22-jul-14 23:59'
                         ,1.8,7.0)
   """
