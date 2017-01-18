#!/usr/bin/python
# -*- coding: utf-8 -*-
# utilities to import sounding from radiosondes or other sources.

from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta
import string
from scipy.optimize import brentq
from scipy import interpolate
import os,os.path

import lg_base.core.array_utils as hau

# atmospheric sounding object
# class temp_sounding(object):
#  pass

def time_offset(time,offset):
    time=time.replace(day=1,hour=0,minute=0,second=0,microsecond=0)
    if offset==0:
        return time
    import calendar
    if offset>0:
        _,ndays=calendar.monthrange(time.year,time.month)
        return time_offset(time+timedelta(days=ndays),offset-1)
    return time_offset(time-timedelta(days=2),offset+1)

def read_grib_file( time, bin_width, max_alt, latitude, longitude ):
    """ read_grib_file(time,bin_width,max_alt,latitude,longitude)
         create sounding from nearest grid point in GRIB file
         with GRIB filenames:
         'grib_yyyymmddhh_yyyymmddhh_L1_L2_G1_G2.grib' yyyymmddhh are
         start and end date of data in the file and L1_l2 are start and
         end latitudes (deg) and G1_G2 are start and end
         longitudes (deg)"""

    rs = []
    return rs

def getAll(nc,varname):
    iarr=[]
    var=nc.variables[varname]
    for d in var.dimensions:
        iarr.append(slice(None,None))#slice(nc_dims[d][0],nc_dims[d][1]))
    #print 'reading variable',alias,'start',starta,'count',counta
    #get requested variable from netcdf
    try:
        if len(iarr)==0:
            tmp=var.getValue()
        elif len(iarr)==1:
            tmp=var[iarr[0]]
        else:
            tmp=var[tuple(iarr)]
    except IndexError:
        print 'Index error reading variable ',varname#,' in file ',filename
        tmp=None
    return tmp

#def read_sounding_file(arg):
def read_sounding_file(instrument,sounding_type,id,start_time,requested_altitudes):
    """ read_sounding_file([instrument,sounding_type,id,start_time,alt_res,max_alt)   
     returns arrays rs.temps(sounding#, rs.dew_points(sounding#,alt_index),
     rs.wdir(sounding#,alt_index), rs.wspd(sounding#,alt_index) along with several
     scalers with sounding info instrument (e.g. 'ahsrl','gvhsrl','mf2hsrl','nshsrl')
     sounding_type may be radiosonde station id, model identifier, or other instrument
     sounding_type (e.g. 'NOAA','ARM',.......
     sounding id (for sounding_type=NOAA, this a 3-letter e.g. 'MSN')
     start_time first time for which the sounding is needed, provided as matplot time
     requested_altitudes is a vector of altitudes at which sounding values are requested (m) 
     returns temp_sounding(object) with all soundings from the current file after the starting time
     returns the last time at which this sounding can be used as rs.expire_time"""
    
    import lg_dpl_toolbox.core.archival as hru
   
    
    if sounding_type[:].find('NOAA raob') >= 0:
        rs =  hau.Time_Z_Group()
        
        time_struct = start_time
        dir_path = hru.get_path_to_data(instrument, start_time)
        filename = dir_path + '/' + '%4i' % time_struct.year + '/' \
            + '%02i' % time_struct.month + '/sondes.' + id[:] + '.nc'
        print 'sounding file--', filename
        if not os.path.exists(filename):
            return None
            #raise RuntimeError, 'sounding file %s does not exist' \
            #    % filename
        nc = Dataset(filename,'r')
        times = getAll(nc,'synTime')

        # combine mandatory and sig height measurements
       
        heights = np.hstack((getAll(nc,'htMan'), getAll(nc,'htSigT')))

        epoch = datetime(1970,1,1,0,0,0)
        t_mask = times<1e36
        for i in range(len(times)):
           t_mask[i] = times[i]<1e36 and any(heights[i,:] < 1e36)

        times=times[t_mask]
      
      
    
        times = [ epoch + timedelta(seconds=soff) for soff in times[:] ]
        
     # select times, one prior to start time --> last profile in file

        #indices = np.arange(len(times))
        #start_index = max(indices[times <= start_time])
        #rs.times = zeros(len(times) - start_index)
        rs.times =np.zeros(len(times))
        #rs.times = times[start_index:]
        #         rs.times = hau.T_Array( rs.times )
        rs.times = times[:]
        
        wmosta = getAll(nc,'wmoStat')  # wmo station number
        stalong = getAll(nc,'staLon')
        stalat = getAll(nc,'staLat')
      
     # combine mandatory and sig height measurements
        
        temps = np.hstack((getAll(nc,'tpMan'), getAll(nc,'tpSigT')))
        pressures = np.hstack((getAll(nc,'prMan'), getAll(nc,'prSigT')))
        dew_points = np.hstack((getAll(nc,'tdMan'), getAll(nc,'tdSigT')))
        wind_dir = np.hstack((getAll(nc,'wdMan'), getAll(nc,'wdSigT')))
        wind_spd = np.hstack((getAll(nc,'wsMan'), getAll(nc,'wsSigT')))
        heights = heights[t_mask,:]        
        temps = temps[t_mask,:]
        pressures = pressures[t_mask,:]
        dew_points = dew_points[t_mask,:]
        wind_dir = wind_dir[t_mask,:]
        wind_spd = wind_spd[t_mask,:]
    
 
        
        [n_soundings, n_heights] = temps.shape
    
     # defined standard atmosphere climatology for use above highest reported level
     # climate=temp_sounding()

        climate =  hau.Time_Z_Group()
        climate.altitudes = np.zeros((n_soundings, 9))
        climate.temps = np.zeros((n_soundings, 9))
        climate.pressures = np.zeros((n_soundings, 9))
        climate.dew_pt = np.zeros((n_soundings, 9))
        climate.wind_spd = np.zeros((n_soundings, 9))
        climate.wind_dir = np.zeros((n_soundings, 9))

     # find the highest valid point in each sounding

        rs.top = np.zeros((n_soundings,))
        rs.bot = np.zeros((n_soundings,))

     # climate.altitudes[0,:]=array([10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000])
             

        for i in range(n_soundings):
            mask = heights[i,:] <= 50000
            if np.any(mask == True):
                rs.top[i] = max(heights[i,mask])
                rs.bot[i] = min(heights[i, temps[i, :] != 99999])
            else:
                rs.top[i] = 0.0
                rs.bot[i] = 0.0
                
            rs.top =hau.T_Array(rs.top)
            rs.bot =hau.T_Array(rs.bot)   
            climate.altitudes[i, :] = np.array([
                10000,
                15000,
                20000,
                25000,
                30000,
                35000,
                40000,
                45000,
                50000,
                ])
            climate.temps[i, :] = np.array([
                223.1,
                216,
                216,
                221,
                226,
                237,
                251,
                265,
                270,
                ])
            climate.pressures[i, :] = np.array([
                264.3,
                120.45,
                54.75,
                25.11,
                11.71,
                5.58,
                2.77,
                1.43,
                0.759,
                ])
            climate.dew_pt[i, :] = np.NaN

       # don't use climatology lower than 2km above highest valid measurement

            climate.altitudes[climate.altitudes <= rs.top[i] + 2000] = \
                9e36
            climate.temps[climate.altitudes <= rs.top[i] + 2000] = 9e36
            climate.pressures[climate.altitudes <= rs.top[i] + 2000] = \
                9e36

     # stack the climatology on top of the observations

        heights = np.hstack((heights, climate.altitudes))
        temps = np.hstack((temps, climate.temps))
        pressures = np.hstack((pressures, climate.pressures))
        dew_points = np.hstack((dew_points, climate.dew_pt))
        wind_dir = np.hstack((wind_dir, climate.wind_dir))
        wind_spd = np.hstack((wind_spd, climate.wind_spd))
        #print heights.shape
        heights_unsorted = heights.copy()
        temps_unsorted = temps.copy()
        pressures_unsorted = pressures.copy()
        dew_points_unsorted = dew_points.copy()
        wind_dir_unsorted = wind_dir.copy()
        wind_spd_unsorted = wind_spd.copy()
       
        for i in range(heights_unsorted.shape[0]):
            indices = np.argsort(heights_unsorted[i,:])
            heights[i,:] = heights_unsorted[i,indices]
            temps[i,:] = temps_unsorted[i,indices]
            pressures[i,:] = pressures_unsorted[i,indices]
            dew_points[i,:] = dew_points_unsorted[i,indices]
            wind_dir[i,:] = wind_dir_unsorted[i,indices]
            wind_spd[i,:] = wind_spd_unsorted[i,indices]
       
     # sort combined file by height and select times of interest
        if 0:
          indices = heights.argsort(axis=1)
          index_a = np.transpose(np.transpose(np.ones(heights.shape,
                               dtype=int)) * np.arange(heights.shape[0]))
          heights = heights[index_a,indices]
          temps = temps[index_a, indices]
          pressures = pressures[index_a, indices]
          dew_points = dew_points[index_a, indices]
          wind_dir = wind_dir[index_a, indices]
          wind_spd = wind_spd[index_a, indices]

        pressures[heights > 1e5] = np.NaN
        temps[heights > 1e5] = np.NaN
        dew_points[heights > 1e5] = np.NaN


     # interpolate to altitude resolution requested
     # and remove missing data points
      
        max_alt = requested_altitudes[-1]
        max_bin = requested_altitudes.shape[0]
        alts = requested_altitudes
        
        n_soundings = len(rs.times)
        #max_bin = round(max_alt / float(alt_res)) + 1
        
     # create sounding arrays as hau class items


    
        rs.altitudes = hau.TZ_Array(np.zeros((n_soundings, max_bin)))
        rs.temps = hau.TZ_Array(np.zeros((n_soundings, max_bin)))
        rs.dew_points = hau.TZ_Array(np.zeros((n_soundings, max_bin)))
        rs.pressures = hau.TZ_Array(np.zeros((n_soundings, max_bin)))
        rs.wind_spd = hau.TZ_Array(np.zeros((n_soundings, max_bin)))
        rs.wind_dir = hau.TZ_Array(np.zeros((n_soundings, max_bin)))

        #rs.times =hau.T_Array(np.zeros(n_soundings))
        rs.times=hau.T_Array(times)
        #rs.expire_time =hau.T_Array(np.zeros(n_soundings))
        rs.stalat =hau.T_Array(np.zeros(n_soundings))
        rs.stalong =hau.T_Array(np.zeros(n_soundings))
        rs.wmosta =hau.T_Array(np.zeros(n_soundings))

        rs.stalat[:] = stalat[0]
        rs.stalong[:] = stalong[0]
        rs.wmosta[:] = wmosta[0]
        rs.station_id = id
        
     # interpolate to lidar altitude scale
       
        for i in range(n_soundings):
            
            rs.altitudes[i, :] = alts
            k = i

            rs.temps[i, :] = np.interp(alts, heights[k, temps[k, :]
                                    != 99999], temps[k, temps[k, :]
                                    != 99999])
                                    
            rs.dew_points[i, :] = np.interp(alts, heights[k, dew_points[k,
                    :] != 99999], dew_points[k, dew_points[k, :]
                    != 99999])
          
            #now using spline fit to pressures
            #should use hydrostatic equation to interpolate
            
            press1=pressures[k,pressures[k,:] !=99999]
            h1=heights[k,pressures[k,:] !=99999]
            temp=interpolate.splrep(h1[~np.isnan(press1)],press1[~np.isnan(press1)])
            rs.pressures[i,:]=interpolate.splev(alts,temp,der=0)
           
        
            rs.wind_spd[i, :] = np.interp(alts, heights[k, wind_spd[k, :]
                    != 99999], wind_spd[k, wind_spd[k, :] != 99999])
            rs.wind_dir[i, :] = np.interp(alts, heights[k, wind_dir[k, :]
                    != 99999], wind_dir[k, wind_dir[k, :] != 99999])
            #if i < n_soundings - 1:
            #    #rs.expire_time[i] = (rs.times[i] + rs.times[i + 1])
            #    rs.expire_time = rs.times[i] + timedelta(seconds=(rs.times[i+1] - rs.times[i]).total_seconds() / 2.0) + timedelta(seconds=60*30)  # add 1/2 hour to point to next sounding
            #else:
            #    #rs.expire_time[i] = rs.times[i] + 0.25 + 1 / 48.0  # new sonde profile expected in 6 1/2 hrs
            #    rs.expire_time = rs.times[i] + timedelta(days=0.25, seconds= 30*60)  # new sonde profile expected in 6 1/2 hrs

     # convert dew point depression to dew point temp

        rs.dew_points = rs.temps - rs.dew_points
        
        plots = 0  # FIXME
        if plots == 1:
            import matplotlib.pylab as plt
            plt.figure(801)

            plt.plot(temps[0, :], heights[0, :] / 1000, 'r',
                     dew_points[0, :], heights[0, :] / 1000)
            fig = plt.grid(True)
            plt.xlabel('deg-K ')
            plt.ylabel('Altitude MSL (km)')
            plt.title('Temperature, dew point')

            plt.figure(802)
            #set_printoptions(threshold=np.NaN)
            
            plt.plot(rs.temps[0, :], rs.altitudes[0, :] / 1000, 'r',
                     rs.dew_points[0, :], rs.altitudes[0, :] / 1000)
            fig = plt.grid(True)
            plt.xlabel('deg-K ')
            plt.ylabel('Altitude MSL (km)')
            plt.title('Temperature, dew point')

            plt.figure(803)
            plt.plot(pressures[0, :], heights[0, :] / 1000, 'r')
            fig = plt.grid(True)
            ax=plt.gca()
            #ax.set_xscale('log')
            plt.xlabel('deg-K ')
            plt.ylabel('Altitude MSL (km)')
            plt.title('Pressures')

            plt.figure(804)
            bin_vec =range(len(heights[0,:]))
            heights[heights>1000]=np.NaN
            plt.plot(heights[0, :] / 1000,bin_vec, 'r')
            fig = plt.grid(True)
            ax=plt.gca()
            #ax.set_xscale('log')
            plt.xlabel('altitudes')
            plt.ylabel('bins')
            plt.title('Heights')

            plt.show()
            raise RuntimeError('deliberate abort')
    else:
        print ' '
        print 'ERROR**************unknown sounding source************'
        print ' '
        rs = []
    return rs



def cal_frost_point(dew_point):
    """Hyland-Wexler calculation of frost point temp (K) given
       the dew point temperature (K).
       Brents method is used to find the point where the vapor
       pressure over water is equal to the vapor pressure over
       ice"""

    nalts = dew_point.shape[0]
    frost_point = np.NaN * np.zeros_like(dew_point)
    log_es = -0.58002206e4 / dew_point + 1.3914993 - 0.48640239e-1 \
        * dew_point + 0.41764768e-4 * dew_point ** 2 - 0.14452093e-7 \
        * dew_point ** 3 + 0.65459673e1 * np.log(dew_point)

    for j in range(nalts):
        if dew_point[j] < 273.15 and np.isfinite(dew_point[ j]):
            try:
                frost_point[ j] = brentq(delta_es, 100, 273, args=log_es[j])
            except ValueError: 
                frost_point[ j] = np.NaN
        else:
            frost_point[ j] = np.NaN
    return frost_point


def delta_es(x, log_es):
    g = log_es - (-0.56745359e4 / x + 6.39225247 - 0.96778430e-2 * x
                  + 0.62215701e-6 * x ** 2 + 0.2074825e-8 * x ** 3
                  - 0.94840240e-12 * x ** 4 + 0.41635019e1 * np.log(x))
    return g


# do not execute on import--only execute when run from cmd line

if __name__ == '__name__':

 # main

    import sys
    print 'start subroutine read_sounding file '
    args = sys.argv[1:]
    print args[1:]
    filename = read_sounding_file(args[0], args[1], args[2])
