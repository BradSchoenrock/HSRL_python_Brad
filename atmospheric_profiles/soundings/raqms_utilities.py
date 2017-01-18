#!/usr/bin/python
# -*- coding: utf-8 -*-
# utilities to import sounding from radiosondes or other sources.

import os,os.path
import numpy as np
from datetime import datetime,timedelta
from lg_base.core.fmt_mpl_date import fmt_mpl_datetime
import string 
from netCDF4 import Dataset
from scipy.optimize import brentq
from scipy import interpolate
import lg_base.core.array_utils as hau

from raob_utilities import getAll

#class raqms_object(object):
#    pass

def time_offset(time,offset):
    return time.replace(hour=0,minute=0,second=0,microsecond=0)+timedelta(days=offset)

def find_raqms_filename(instrument,start_time):
    """return file name of raqms profile curtain"""
    import lg_dpl_toolbox.core.archival as hru

    time_struct = start_time
    dir_path = '%s/%4i/%02i/raqms_soundings'  % \
        (hru.get_path_to_data(instrument, start_time),
         time_struct.year, time_struct.month)

    print "fetching raqms sounding curtain from '"+dir_path+"'",start_time
    if not os.path.exists(dir_path):
        print ("doesn't exist")
        return None
    filenames=os.listdir(dir_path)
    filenames.sort()
    filename=None
    for file in filenames:
        index=file.find('.raqms.curtain.')
        filetime=datetime(int(file[index-8:index-4]),int(file[index-4:index-2]),int(file[index-2:index]),0,0,0)
        if filetime<=time_struct:
            filename=dir_path + '/' + file
    print "using sounding data from raqms file = '",filename,"'"
            #break

    return filename


def read_raqms_file(instrument,start_time):
    """read raqms file between start and end time

    instrument - e.g. 'gvhsrl'
    start_time -  datetime object    
    
    """

    raqms=hau.Time_Z_Group(altname='model_level_alts')

    filename=find_raqms_filename(instrument,start_time)
    if not filename:
        return None
        
    nc = Dataset(filename,'r')
    times = getAll(nc,'time')
    aircraft_alts = getAll(nc,'alt')
    pressures = getAll(nc,'pressure')
    temperatures = getAll(nc,'temperature')
    model_level_alts = getAll(nc,'altitude')

    relative_humidity = getAll(nc,'rh')
    latitude=getAll(nc,'lat')
    longitude=getAll(nc,'lon')
    u_vel=getAll(nc,'uvel')
    v_vel=getAll(nc,'vvel')
    ext_total = getAll(nc,'ext_tot')
    ext_dust  = getAll(nc,'ext_dust')
    ext_salt  = getAll(nc,'ext_salt')

    base_time=datetime(start_time.year,start_time.month,start_time.day,0,0,0)
    #np.fix(start_time)
    #time=times.astype('float64')

    #convert raqms seconds from start of day to python datetimes
    #times=base_time + time/(3600.0*24.0)
    times=hau.T_Array([ base_time + timedelta(seconds=float(x)) for x in times ])

    assert(times.size>0)
  
    selectedMask=(times > start_time)
    for i,x in enumerate(selectedMask):
        fi=i
        if x:
            if i>0:
                selectedMask[i-1]=True
            break
    selectedMask[-1]=True

    selectedTimes = np.arange(times.size)[selectedMask]

    raqms.latitude=hau.T_Array(latitude[selectedTimes])
    raqms.longitude=hau.T_Array(longitude[selectedTimes])
    raqms.pressures=hau.TZ_Array(pressures[selectedTimes,:])
    raqms.temperatures=hau.TZ_Array(temperatures[selectedTimes,:])
    raqms.ext_total=hau.TZ_Array(ext_total[selectedTimes,:])
    raqms.ext_dust=hau.TZ_Array(ext_dust[selectedTimes,:])
    raqms.ext_salt=hau.TZ_Array(ext_salt[selectedTimes,:])
    raqms.relative_humidity = hau.TZ_Array(relative_humidity[selectedTimes,:])
    raqms.u_vel=hau.TZ_Array(u_vel[selectedTimes,:])
    raqms.v_vel=hau.TZ_Array(v_vel[selectedTimes,:])
    raqms.model_level_alts=hau.TZ_Array(model_level_alts[selectedTimes,:]*1000.0)
    raqms.times=times[selectedTimes]
   
    return raqms


def select_raqms_profile(soundings,request_time,requested_altitudes,offset=0):
    """selects sounding prior to request_time from soundings -- the sounding
       is returned in a Time_Z_Group as Z_arrays"""

    if soundings is None or soundings.times.size == 0:
        raise RuntimeError, "select_faqms_profile: No soundings for %s " %\
              request_time

    import atmospheric_profiles.soundings.sounding_utilities as su
    sounding = hau.Time_Z_Group()

    sounding.altitudes = hau.Z_Array(requested_altitudes)
    max_alt = requested_altitudes[-1]
    max_bin = len(requested_altitudes) 
    index = sum(soundings.times <= request_time) -1+offset

    if index<0 or index>=len(soundings.times):
        return None   
    #initialize variables for inclusion in T_Z_Group
    sounding.temps = hau.Z_Array(np.zeros((max_bin)))
    sounding.dew_points = hau.Z_Array(np.zeros(max_bin))
    sounding.frost_points = hau.Z_Array(np.zeros(max_bin))
    sounding.pressures = hau.TZ_Array(np.zeros(max_bin))
    sounding.ext_total = hau.TZ_Array(np.zeros(max_bin))
    sounding.ext_salt = hau.TZ_Array(np.zeros(max_bin))
    sounding.wind_spd = hau.TZ_Array(np.zeros(max_bin))
    sounding.wind_dir = hau.TZ_Array(np.zeros(max_bin))



    #sounding.times is a single time at this point, however it will later be included
    #in a list of all the soundings used in this processing request. In order that it
    #be treated properly it must be defined as a T_Array

    sounding.times=hau.T_Array([soundings.times[index]])
    sounding.latitude=hau.T_Array([soundings.latitude[index]])
    sounding.longitude=hau.T_Array([soundings.longitude[index]])
  
    #sounding.times=hau.T_Array([soundings.times[index]])
  
    #interpolate model levels to lidar bin altitudes

    #temp=interpolate.splrep(soundings.model_level_alts[index,-1::-1] \
    #     ,soundings.temperatures[index,-1::-1])
    #sounding.temps=interpolate.splev(sounding.altitudes,temp,der=0)

    temp=interpolate.splrep(soundings.model_level_alts[index,-1::-1] \
                            ,soundings.pressures[index,-1::-1])
    sounding.pressures=interpolate.splev(sounding.altitudes,temp,der=0)



    sounding.temps=np.interp(sounding.altitudes \
                             ,soundings.model_level_alts[index,-1::-1] \
                             ,soundings.temperatures[index,-1::-1])
    
    #calculate dew point at model levels for selected profile
    dew_pts=su.cal_dew_point(soundings.relative_humidity[index,:] \
                             ,soundings.temperatures[index,:])
    frost_pts=su.cal_frost_point(dew_pts)

    #calculate wind speed and direction from u and v
    u_vel = soundings.u_vel[index,-1::-1] 
    v_vel = soundings.v_vel[index,-1::-1]



    wind_spd =np.sqrt(u_vel**2 +v_vel**2)
    wind_dir =np.arctan(v_vel/u_vel)*180.0/np.pi



    for i in range(len(u_vel)):
        if (u_vel[i] < 0 and v_vel[i]) < 0: 
            wind_dir[i] = 180.0 - wind_dir[i]
        elif  (u_vel[i] > 0 and v_vel[i]) > 0:
            wind_dir[i]= 180.0 +wind_dir[i]
        elif u_vel[i]<0:
            wind_dir[i]= 270.0 - wind_dir[i]
        else:   
            wind_dir[i]=np.nan



    #interpolate to lidar bin altitudes
    sounding.frost_points=np.interp(sounding.altitudes \
                                    ,soundings.model_level_alts[index,-1::-1],frost_pts[-1::-1])
    sounding.dew_points=np.interp(sounding.altitudes \
                                  ,soundings.model_level_alts[index,-1::-1],dew_pts[-1::-1])
    sounding.ext_total=np.interp(sounding.altitudes\
                                 ,soundings.model_level_alts[index,-1::-1]\
                                 ,soundings.ext_total[index,-1::-1])
    sounding.ext_salt=np.interp(sounding.altitudes\
                                ,soundings.model_level_alts[index,-1::-1]\
                                ,soundings.ext_salt[index,-1::-1])
    sounding.ext_dust=np.interp(sounding.altitudes\
                                ,soundings.model_level_alts[index,-1::-1]\
                                ,soundings.ext_dust[index,-1::-1])



    sounding.wind_dir = np.interp(sounding.altitudes \
                                  ,soundings.model_level_alts[index,-1::-1],wind_dir)                           
    sounding.wind_spd = np.interp(sounding.altitudes \
                                  ,soundings.model_level_alts[index,-1::-1],wind_spd)                             

    sounding.top=sounding.altitudes[-1]
    sounding.bot=sounding.altitudes[0]

    #plt.figure(1)
    #plt.plot(temperatures,altitudes,dew_points,altitudes)

    #plt.figure(2)
    #plt.plot(ext_total,altitudes)
    #plt.show()
    return sounding
