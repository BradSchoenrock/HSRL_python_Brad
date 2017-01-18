import numpy as np
from datetime import timedelta

def set_z_resolution(min_alt=None, max_alt=None, binwidth=None
                     , number_y_pixels = None,manual=None,n_dr=None):
    """set_z_resolution(min_alt=None, max_alt=None, binwidth=None
                     , number_y_pixels = None,manual=None,n_dr=None):
       sets the altitude resolution of the output data, delta_z 
       min_alt   = min altitude requested 
       max_alt   = max altitude requested 
       number_y_pixels = number y pixels in image
       binwidth        = instrument altitude binwidth (meters)
       for manual mode (meters):
       manual          = dz    altitude res = dz meters 
       for bincount mode (bincount):
       n_dr            = nbins     altitude res = nbins * binwidth 
       for auto mode (best fit. requires alts, binwidth, and number_y_pixels):
       don't specifiy manual and n_dr => altitude res = (max_alt-min_alt)/number_y_pixels
       """
    if n_dr!=None:
        delta_z = np.int(n_dr) *binwidth
        n_range_ave = np.int(n_dr)
        if n_range_ave < 1 :
            n_range_ave =1
    elif manual!=None:
        delta_z = np.float(manual)
        n_range_ave = np.int(np.round(delta_z /binwidth))
        if n_range_ave <1:
            n_range_ave =1
    else:
         delta_z = ((np.float(max_alt)-np.float(min_alt))\
            / number_y_pixels)
         n_range_ave = np.int(delta_z / binwidth)
         if delta_z < binwidth :
             delta_z = binwidth
             n_range_ave = 1        
    return delta_z,n_range_ave


def set_time_resolution(delta_t=None, number_x_pixels=None
                        ,frame_dt=0,bin_dt=None,manual=None,n_dt=None):
    """set_time_resolution(start_time=None,end_time=None, number_x_pixels=None
                        ,frame_dt=0,bin_dt=None,manual=None,n_dt=None)
       sets time resolution of output data
       delta_t    = duration of data segment (timedelta)
       x_image_size_pixels = number time pixels in image
       frame_dt   = min time for data frame, (e.g. polarization processing interval)
       bin_dt     = instrument integration time for one  bin

       for manual mode:
       manual    = dt (in floating seconds)
       for bincount mode:
       n_dt = n  =  n times instrument bin time
       for automatic mode:
       don't specify n_dt or manual => (end_time - start_time) / x_image_size_pixels
                   values rounded to even number of range bin widths
       time_res = time resolution output (timedelta)
       ntime_ave= number of records to pre-average
 """

 
    if n_dt :
        time_res = np.float(n_dt)*bin_dt
        time_res = timedelta(seconds = max(time_res, bin_dt, np.float(frame_dt))) 
        ntime_ave = np.int(time_res.total_seconds()/bin_dt)
        
    elif manual :
        time_res = np.float(manual)
        time_res = timedelta(seconds = max(time_res, bin_dt, frame_dt))
        ntime_ave = np.int(time_res.total_seconds() / bin_dt)
        
    else:
        #time resolution in seconds
        time_res = timedelta(seconds=max(delta_t.total_seconds() \
                   /np.float(number_x_pixels), bin_dt, frame_dt))
        if time_res.total_seconds() < bin_dt :
            time_res = timedelta(seconds = bin_dt)
            print 'WARNING---requested time resolution less than data resolution'
           
        ntime_ave = max(np.int(0.95*time_res.total_seconds() / bin_dt),1)
    return time_res, ntime_ave
     
