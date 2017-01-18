#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import os.path
import traceback
from time import sleep
import numpy as np
from datetime import datetime,timedelta
#from matplotlib.dates import date2num
import hsrl.data_stream.hsrl_read_utilities as hru
# formatting utils for matplotlib datetimes
import lg_base.core.decoratortools as nt
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
#import hsrl.calibration.cal_read_utilities as cru


class cal_vectors(object):

    def __init__( self, instrument, time, max_range_bin, alternate_cal_dir=None ):
        """initialize calibration vector structure"""

        self.instrument = instrument
        self.max_range_bin = max_range_bin
        self.alternate_cal_dir=alternate_cal_dir
        self.read(time)
        
        
        
    def read(self, time):
        """read baseline,geo,diff_geo, and i2_scan files
           n_geo will return empty if instrument <> 'gvhsrl' """

        self.baseline = hru.read_baseline(self.instrument, time
                ,self.alternate_cal_dir,max_range_bin=self.max_range_bin
                ,old_calvec=self.baseline if hasattr(self,'baseline') else None)

        self.geo = hru.read_geo_corr(self.instrument, time
                ,self.alternate_cal_dir,max_range_bin=self.max_range_bin
                ,old_calvec=self.geo if hasattr(self,'geo') else None)
      
        self.n_geo = hru.read_nadir_geo_corr(self.instrument,time
                ,self.alternate_cal_dir,max_range_bin=self.max_range_bin
                ,old_calvec=self.n_geo if hasattr(self,'n_geo') else None)

        self.qw_baseline = hru.read_qw_baseline(self.instrument,time
                 ,self.alternate_cal_dir,old_calvec=self.qw_baseline
                 if hasattr(self,'qw_baseline') else None)
       
        self.diff_geo = hru.read_diff_geo(self.instrument, time
              ,self.alternate_cal_dir,max_range_bin=self.max_range_bin
              ,old_calvec=self.diff_geo if hasattr(self,'diff_geo') else None)
            
        self.cpol_diff_geo = hru.read_cross_poll_diff_geo(self.instrument, time
              , self.alternate_cal_dir, max_range_bin=self.max_range_bin
              ,old_calvec=self.cpol_diff_geo if hasattr(self,'cpol_diff_geo') else None)
            
        self.i2a_diff_geo = hru.read_i2a_diff_geo(self.instrument,time
               , self.alternate_cal_dir,max_range_bin=self.max_range_bin
               ,old_calvec=self.i2a_diff_geo if hasattr(self,'i2a_diff_geo') else None)
        
        self.diff_1064_532_geo = hru.read_diff_1064_532_geo(self.instrument,time
               , self.alternate_cal_dir,max_range_bin=self.max_range_bin
               ,old_calvec=self.diff_1064_532_geo if hasattr(self,'diff_1064_532_geo') else None)
        if 0:
            print self.diff_1064_532_geo.data.shape
            print len(self.diff_1064_532_geo.data[:,1])
            import matplotlib.pylab as plt
            plt.figure(999)
            plt.plot(1+self.diff_1064_532_geo.data[:,1],self.diff_1064_532_geo.data[:,0]/1000.0,'r')
            plt.xlabel('1064/532 differential geo correction')
            plt.ylabel('Altitude (km)')
            plt.grid(True)
 
        self.i2scan = hru.read_i2_scan(self.instrument, time
              , self.alternate_cal_dir,old_calvec=self.i2scan if hasattr(self,'i2scan') else None)

        self.i2a_temp_table = hru.read_i2a_temp_table(self.instrument,time
              ,self.alternate_cal_dir,old_calvec=self.i2a_temp_table if hasattr(self,'i2a_temp_table') else None)
      
        
        self.pol_cal = hru.read_pol_cal(self.instrument, time
              ,self.alternate_cal_dir,old_calvec=self.pol_cal if hasattr(self,'pol_cal') else None)

        times=[self.baseline.expire_time,
                               self.geo.expire_time,
                               self.diff_geo.expire_time,
                               self.i2scan.expire_time,
                               self.qw_baseline.expire_time,
                               self.i2a_diff_geo.expire_time,
                               self.pol_cal.expire_time,
                               self.n_geo.expire_time]

        i2offset=hru.read_i2_offset(self.instrument,time,self.alternate_cal_dir)
        if i2offset is not None:
          self.i2offset=i2offset.i2offset
          times.append(i2offset.expire_time)
          print 'I got i2offset of %i ending %s' %(i2offset.i2offset,i2offset.expire_time.strftime('%Y%m%d-%H%M%S'))


       # find first expiration time
        self.expire_time = min(*times)
       
