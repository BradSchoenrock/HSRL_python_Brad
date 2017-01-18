#!/usr/bin/python
# -*- coding: utf-8 -*-

import bz2
from os import system
import os.path
import numpy as np
from lg_base.core.fmt_mpl_date import fmt_mpl_datetime,fmt_mpl_time
import time
import os
from datetime import datetime, timedelta
import string 
#import hsrl.data_stream.main_routines as mr
import lg_base.core.array_utils as hau
#import hsrl.data_stream.input_translators as it
#from hsrl.data_stream.preprocess_and_time_ave import preprocess_and_time_ave
import lg_base.core.walkdir as walkdir

import json
import struct

from lg_base.formats.vector_table import calibration_vector_table

def pause():
    try:
        print ' '
        print 'press enter to proceed'
        input()
    except KeyboardInterrupt:
        print 'proceed'



def fetch_data( instrument, interval_start_time, interval_end_time,
    max_range_bin, requested_vars,cdf_to_hsrl,pre_process=None,dpl_librarian=None,dpl_zookeeper=None):
    """appends data to rs_mem until it contains requested time interval
       or extends to current time.
       instrument          = which data set----'ahsrl','gvhsrl','mf2hsrl','nshsrl'
       rs_mem              = append new data to this structure
       interval_start_time = start_time of requested data as a python datetime
       interval_end_time   = end_time of requested data as a python datetime
       data_end_time       = time of last measurement in rs_mem as a python datetime
       max_range_bin       = largest range bin to read
       data_request        = read this subset of data from raw netcdf
                             e.g. 'images','housekeeping'
       small_memory        =[enable,number_of_shots_to_pre_average]      """

    
    from lg_dpl_toolbox.dpl.NetCDFZookeeper import GenericTemplateRemapNetCDFZookeeper 
    from hsrl.dpl.HSRLLibrarian import HSRLLibrarian 

    # print num2date(c_time)
    if dpl_librarian:
        lib=dpl_librarian
    else:
        lib=HSRLLibrarian(instrument=instrument)
    if dpl_zookeeper:
        zoo=dpl_zookeeper
    else:
        zoo=GenericTemplateRemapNetCDFZookeeper(instrument,requested_vars, max_range_bin,pre_process)

    #rs_mem=None

    for uri in lib(start=interval_start_time,end=interval_end_time):
        # reading first record with no prior data in rs_mem?
      try:
        filename=zoo(uri)
        print 'datafile   ',filename
        rs_tail=zoo.open(filename,firsttime=interval_start_time,lasttime=interval_end_time)
        found_vars=zoo.getFoundVars(filename)
        if rs_tail is None or rs_tail.times.size==0:
            continue
    
        missing_data_check(rs_tail)
        cdf_to_hsrl(rs_tail,found_vars)
        if rs_tail.times.size==0:
            continue
        yield rs_tail
      except MemoryError as e:
        import traceback
        traceback.print_exc()
        print e
        continue

       
    #assert(len(rs_mem.times.shape) == 1)
   
    return #rs_mem

def last_fetch(cdf_to_hsrl,pre_process,zoo):#FIXME this is soooo hacky
        rs_tail=pre_process.flush()
        if rs_tail is None or rs_tail.times.size==0:
            return None
        found_vars=zoo.lastfoundvars #IT BURNS
    
        missing_data_check(rs_tail)
        cdf_to_hsrl(rs_tail,found_vars)
        return rs_tail


from lg_dpl_toolbox.core.archival import get_path_to_data

import time

def read_i2_offset(instrument, theTime,alternate_cal_dir=None,max_range_bin=None,filename=None,expire_time=None,old_calvec=None):
    import lg_base.formats.calvals as cru
    [filename, expire_time] = find_cal_file(instrument, 'i2offset', theTime
            ,alternate_cal_dir,filename=filename,expire_time=expire_time)

    if filename is None:
        return None
    if old_calvec is not None and filename == old_calvec.filename:
        return old_calvec
    if theTime is not None:
        print 'read_i2_offset: theTime = ', fmt_mpl_datetime(theTime)
    i2_offsetfile=cru.calvals_class(filename)
    if theTime is None: #only a parse test
        return None
    i2offset=i2_offsetfile(theTime)
    if i2offset is None or i2offset==0:
        return None
    from collections import namedtuple
    return namedtuple("I2Offset","i2offset expire_time")(i2offset['i2offset'],i2offset['next_cal_time'])


def read_baseline(instrument, theTime,alternate_cal_dir=None,max_range_bin=None,filename=None,expire_time=None,old_calvec=None):

   # reads fixed ca

   # get baseline correction
    if theTime is not None:
        print 'read_baseline: theTime = ', fmt_mpl_datetime(theTime)
    [filename, expire_time] = find_cal_file(instrument, 'baseline', theTime
            ,alternate_cal_dir,filename=filename,expire_time=expire_time)

    if filename is None:
        raise RuntimeError, \
            "no baseline correction file - install 'baseline_xxxxx.blc' in month directory for date %s" % (theTime)
    if old_calvec is not None and filename == old_calvec.filename:
        return old_calvec

    cal_vec = calibration_vector_table(filename,expire_time,max_range_bin=max_range_bin,requireColumns=True
        ,fallbackColumns=["bin_num","combined_hi","combined_lo","molecular","crosspol","mol_I2A","comb_1064"])
    bl_header=cal_vec.header
    #bl_corr=cal_vec.data

    #get ave 532 nm energy per shot when baseline file was created 
    energy_index = bl_header.find('#ave energy per shot=')
    if energy_index <= 0:
        print ' '
        print '************ERROR---energy not found in baseline file header ************'
        print '*****using 0.075 mJ/shot***** '
        baseline_532_energy = 0.075
    else:

     # raise RuntimeError, "no energy found in baseline correction file"
         # energy found in baseline file

        print 'index= ', energy_index, ' baseline energy (mJ/shot)= ', \
            float(bl_header[energy_index + 22:energy_index + 30])
        baseline_532_energy =float(bl_header[energy_index + 22:energy_index + 30])
        print 'baseline energy = ',baseline_532_energy
    have1064=False
    for f in cal_vec.fields:
        if '1064' in f:
            have1064=True
            continue
        v=getattr(cal_vec,f)
        v /= baseline_532_energy
    print bl_header
    if have1064:
        #get average 1064 nm energy per shot when baseline file was created
        energy_index = bl_header.find('#ave 1064 energy per shot=')
        if energy_index <= 0:
            print ' '
            print '************WARNING---1064 energy not found in baseline file header ************'
            print '*****using 532nm value ***** '
            baseline_1064_energy = baseline_532_energy 
        else:
            # 1064 energy found in baseline file
            
            print 'index= ', energy_index, ' baseline energy (mJ/shot)= ', \
                float(bl_header[energy_index + 27:energy_index + 35])
            baseline_1064_energy =float(bl_header[energy_index + 27:energy_index + 35])
            print 'baseline 1064 energy = ',baseline_1064_energy
        for f in cal_vec.fields:
            if '1064' not in f:
                continue
            v=getattr(cal_vec,f)
            v /= baseline_1064_energy
        
    print 
    print 'baseline corr-', filename
    print
    return cal_vec

def read_i2a_temp_table(instrument, time,alternate_cal_dir=None,max_range_bin=None,filename=None,expire_time=None,old_calvec=None):
   
    if not instrument == 'bagohsrl':
       print 'no i2a_temp_table for ',instrument
       return
   
    # get i2a_temp_table 
   
    [filename, geo_expire] = find_cal_file(instrument, 'i2a_temp_table', time
                , alternate_cal_dir,filename=filename,expire_time=expire_time)
    if filename is not None and old_calvec is not None and filename == old_calvec.filename:
        return old_calvec

    cal_vec = calibration_vector_table(filename,geo_expire,max_range_bin=max_range_bin)
   
    if filename is None:
        print "WARNING----no  i2a_temp_table - install 'i2a_temp_table.cal' in month directory for date %s" % (time)
    print     
    print 'i2a_temp_table-------', filename
    print
    return cal_vec

def read_geo_corr(instrument, time,alternate_cal_dir=None,max_range_bin=None,filename=None,expire_time=None,old_calvec=None):

   # get normal geo correction
   
    [filename, geo_expire] = find_cal_file(instrument, 'geo', time
                , alternate_cal_dir,filename=filename,expire_time=expire_time)
   
    if filename is None:
        raise RuntimeError, \
            "no geo correction file - install 'geofile_default.geo' in month directory for date %s" % (time)
    if old_calvec is not None and filename == old_calvec.filename:
        return old_calvec
    cal_vec = calibration_vector_table(filename,geo_expire,max_range_bin=max_range_bin,requireColumns=True,fallbackColumns=('Range','geo_correction'))
    geo_header=cal_vec.header
    geo_corr=cal_vec.data[:4000,:]
         
   # print 'geo_header = ', geo_header, 'geo_corr = ', geo_corr
   # multiply by r-squared in km^2
   
    geo_corr[:, 1] = geo_corr[:, 1] * geo_corr[:, 0] * geo_corr[:, 0]/1.0e6
    print
    print 'geofile-------', filename
    print

    #find wfov exponential correction term if it exists
    wfov_corr_index = geo_header.find('#wfov_corrected by exp(-range * ')
    if wfov_corr_index <= 0:
        wfov_exp_corr = None
        print 
        print 'no exponetial wfov correction applied'
        print 
    else:
        #wfov exponential correction is written into zeroth element of the wfov correction vector
        wfov_exp_corr =float(geo_header[wfov_corr_index +31:wfov_corr_index + 42])
        print
        print 'wfov exponential correction applied '
        print '         wfov * exp(-range *',wfov_exp_corr,')'
        geo_corr[0,2]=wfov_exp_corr
    return cal_vec

def read_pol_cal(instrument, time, alternate_cal_dir=None,max_range_bin=None,filename=None,expire_time=None,old_calvec=None):

   # get polarization calibration 
   
    [filename, expire_time] = find_cal_file(instrument, 'pol', time, alternate_cal_dir,filename=filename,expire_time=expire_time)

    if filename is not None and old_calvec is not None and filename == old_calvec.filename:
        return old_calvec

    cal_vec = calibration_vector_table(filename,expire_time,isJson=True,max_range_bin=max_range_bin)

    if filename is None:
       print '******pol_cal file not found--no quarter wave plate angle corr applied******'
         

    print 'pol_cal file-------', filename

    return cal_vec

def read_nadir_geo_corr(instrument, time,alternate_cal_dir=None,max_range_bin=None,filename=None,expire_time=None,old_calvec=None):

   # get normal geo correction
   
    [filename, n_geo_expire] = find_cal_file(instrument, 'n_geo', time
         ,alternate_cal_dir,filename=filename,expire_time=expire_time)
   
    if filename is None:
        #raise RuntimeError, \
        #   "no nadir geo correction file"
        cal_vec = calibration_vector_table(filename,datetime(2100,1,1,0,0,0),max_range_bin=max_range_bin)

    elif old_calvec is not None and filename == old_calvec.filename:
        return old_calvec
    else:
        cal_vec = calibration_vector_table(filename,n_geo_expire,max_range_bin=max_range_bin)
        n_geo_header=cal_vec.header
        n_geo_corr=cal_vec.data
         
        # print 'n_geo_header = ', n_geo_header, 'n_geo_corr = ', n_geo_corr
        # multiply by r-squared in km^2

        n_geo_corr[:, 1] = n_geo_corr[:, 1] *n_geo_corr[:, 0] * n_geo_corr[:, 0] \
              / 1e6
        print
        print 'nadir geofile-------', filename
        print
    return cal_vec

def read_diff_geo(instrument, time,alternate_cal_dir=None,max_range_bin=None,filename=None,expire_time=None,old_calvec=None):
    
   # get differential geo correction

    [filename, dgeo_expire] = find_cal_file(instrument, 'd_geo', time
                ,alternate_cal_dir,filename=filename,expire_time=expire_time)

    if filename is not None and old_calvec is not None and filename == old_calvec.filename:
        return old_calvec

    cal_vec = calibration_vector_table(filename,dgeo_expire,max_range_bin=max_range_bin)
    dgeo_header=cal_vec.header
    dgeo_corr=cal_vec.data
    if filename is None:
       print '******differential geo file was not found--no diff geo corr applied******'
    else:
        # difference from unity
        # this alows adjustment by cal_scale factor in process_data
        dgeo_corr[:, 1:4] = dgeo_corr[:, 1:4] - 1
    print
    print 'diff geofile--', filename
    print
       
    return cal_vec

def read_cross_poll_diff_geo(instrument, time,alternate_cal_dir=None,max_range_bin=None,filename=None,expire_time=None,old_calvec=None):
    
   # get differential geo correction

    [filename, expire] = find_cal_file(instrument, 'cpol_diff_geo', time
                ,alternate_cal_dir,filename=filename,expire_time=expire_time)

    if filename is not None and old_calvec is not None and filename == old_calvec.filename:
        return old_calvec

    cal_vec = calibration_vector_table(filename,expire,max_range_bin=max_range_bin)
    dgeo_header=cal_vec.header
    dgeo_corr=cal_vec.data

    if filename is None:
       print '******cross_pol differential geo file was not found--no cross pol diff geo corr applied******'
    else:
        dgeo_corr[:, 1] = dgeo_corr[:, 1] - 1
    print
    print 'cross pol diff geofile--', filename
    print    
    return cal_vec

def read_qw_baseline(instrument, time,alternate_cal_dir=None,max_range_bin=None,filename=None,expire_time=None,old_calvec=None):
    """reads baseline dependence on quarterwave plate rotation angle for
       Matt Hymman's modification to the gvhsrl
       multiplier to the baseline correction is supplied as a function of 
       angle in degrees for each instrument channel"""

    [filename, qw_baseline_expire] = find_cal_file(instrument, 'qw_baseline'
               , time,alternate_cal_dir,filename=filename,expire_time=expire_time)
    if filename is not None and old_calvec is not None and filename == old_calvec.filename:
        return old_calvec

    cal_vec = calibration_vector_table(filename,qw_baseline_expire,max_range_bin=max_range_bin)
    if filename is None:
       print '******qw_baseline file not found--no quarter wave plate angle corr applied******'

    print 'qw_baseline_file--', filename

    return cal_vec


def read_i2a_diff_geo(instrument, time,alternate_cal_dir=None,max_range_bin=None,filename=None,expire_time=None,old_calvec=None):

   # This channel is only implemented on the bagohsrl
    if not instrument == 'bagohsrl':
        cal_vec = calibration_vector_table(None,datetime(2100,1,1,0,0,0),max_range_bin=max_range_bin)
    else:

        #get i2a differential geo file
        [filename, i2a_dgeo_expire] = find_cal_file(instrument, 'i2a_d_geo', time,alternate_cal_dir,filename=filename,expire_time=expire_time)

        if filename is not None and old_calvec is not None and filename == old_calvec.filename:
            return old_calvec

        cal_vec = calibration_vector_table(filename,i2a_dgeo_expire,max_range_bin=max_range_bin)
        i2a_dgeo_header=cal_vec.header
        i2a_dgeo_corr=cal_vec.data
        if filename is None:
           print '******i2a differential geo file was not found--no i2a diff geo corr applied******'
        else:
            # difference from unity
            # this alows adjustment by cal_scale factor in process_data
            i2a_dgeo_corr[:, 1] = i2a_dgeo_corr[:, 1] - 1
            print
            print 'i2a_diff geofile--', filename
            print
    return cal_vec

def read_diff_1064_532_geo(instrument, time,alternate_cal_dir=None,max_range_bin=None,filename=None,expire_time=None,old_calvec=None):

   # This channel is implemented only on the bagohsrl system and on the ahsrl after the 2015 rebuild
    if not (instrument == 'bagohsrl' or instrument == 'ahsrl'):
        cal_vec = calibration_vector_table(None,datetime(2100,1,1,0,0,0),max_range_bin=max_range_bin)
    else:
        #get 1064_532 differential geo file
        #contains the ratio of the 1064nm/532nm channel gains as a function of range
        [filename, dgeo_1064_532_expire] = find_cal_file(instrument, 'd_geo_1064_532'
                     , time,alternate_cal_dir,filename=filename,expire_time=expire_time)

        if filename is not None and old_calvec is not None and filename == old_calvec.filename:
            return old_calvec
        
        cal_vec = calibration_vector_table(filename,dgeo_1064_532_expire,max_range_bin=max_range_bin)
        dgeo_1064_532_header=cal_vec.header
        dgeo_1064_532_corr=cal_vec.data
        if filename is None:
           print '******1064nm/532nm differential geo file was not found--no 1064_532 diff geo corr applied******'
        else:
            # difference from unity
            # this alows adjustment by cal_scale factor in process_data
            dgeo_1064_532_corr[:, 1] = dgeo_1064_532_corr[:, 1] - 1
            print
            print '1064nm/532nm differential geofile--', filename
            print
    return cal_vec

def read_i2_scan(instrument, time,alternate_cal_dir=None,max_range_bin=None,filename=None,expire_time=None,old_calvec=None):

    # get i2 scan file

    if not instrument == 'rbhsrl':
       [filename, i2_expire] = find_cal_file(instrument, 'i2', time,alternate_cal_dir,filename=filename,expire_time=expire_time)
    else: 
        [filename,i2_expire] = find_cal_file(instrument,'rb',time,alternate_cal_dir,filename=filename,expire_time=expire_time)
    if filename is None:
        raise RuntimeError, \
            "no i2-scan file - install 'i2-default-scan-XXX.cal' in month directory for date %s" % (time)
    if filename is not None and old_calvec is not None and filename == old_calvec.filename:
        return old_calvec
    cal_vec = calibration_vector_table(filename,i2_expire,max_range_bin=max_range_bin)
    i2_header=cal_vec.header
    i2_scan=cal_vec.data
    print
    print 'absorption filter------------', filename
    print
   # print i2_header
   
    cal_vec.Cam = None
    cal_vec.Cam_i2a = None
    
    lines=cal_vec.header.strip().split('\n')
    for l in lines:
        if '# Cam =' in l: 
            cal_vec.Cam = float(l.split()[3])
        elif '# Cam_i2a =' in l: 
            cal_vec.Cam_i2a = float(l.split()[3])
   
    return cal_vec


def find_cal_file(instrument, file_type, start_time,alternate_cal_dir=None,filename=None,expire_time=None):
    """returns None if cal file was not found
       this is used to find active calibration directory for this instrument
       and time"""

 # raise RuntimeError =[] for case when no file is found

 # routines requires both name and extension for match
    if filename is not None:
        return (filename,expire_time or datetime(2200,1,1,0,0,0))

    import hsrl.dpl.calibration.TableLibrarian as tl
    return tl.findCalFile(instrument,file_type,start_time,alternativePath=alternate_cal_dir,expire_time=expire_time)


def read_processed_hsrl_data(filename):
    """read a netcdf file created by write_netcdf. Data will be
       restored to python structures"""
    from pycdf import CDF, CDFError
    nc=CDF(filename)
    nc_vars = nc.variables().keys()
    for name in nc_vars:
       print name,'            \t\t  ', nc.var(name).dimensions()
       tmp = nc.var(name).get()

def missing_data_check(raw):       
    """check for miss match between size of time dimensioned variable and time array--
       set variable to np.NaN's if not matched"""
    length_time_vec = raw.times.shape[0]
    for (name,value) in vars(raw).items():
        if isinstance(value, hau.T_Array) :
            #print '    s-check ',     name, value.shape[0],raw.times.shape[0]
            if not value.shape[0] == length_time_vec:
                print
                print 'hru--size error ',name,'has ',value.size,', raw.times has ',raw.times.size,' points.'
                print
                newshape=list(value.shape)
                newshape[0]=length_time_vec
                try:
                    vars(raw)[name]= type(value)(np.NaN * np.zeros(newshape),dtype=value.dtype,summode=value.summode)
                except TypeError:
                    vars(raw)[name]= type(value)(np.zeros(newshape),dtype=value.dtype,summode=value.summode)
