import numpy as np
from datetime import datetime, timedelta
from scipy.optimize import brentq
import raob_utilities as raob
import raqms_utilities as raqms
import lg_base.core.array_utils as hau
import logging
LOG = logging.getLogger(__name__)
DBG=logging.debug

        
class sounding_archive(object):

    def __init__(self, instrument,sounding_type,sounding_id, start_time \
                 ,requested_altitudes,offset=0,allow_fallback=True):
                 #,alt_res,max_alt):
        """finds a file of sounding_type == 'type' containing soundings between
        start_time and end_time and reads these soundings--the following
        variables are retained:

        # single values self.instrument self.sounding_type
        # self.sounding_id self.alt_res self.max_alt
        
        #the following are 1-d arrays
        self.time      
        self.latitude  
        self.longitude
        
        #the following are 2 or 3-d arrays (time,lat,long)
        self.geo_height  
        self.temperature 
        self.pressure    
        self.dew_point   
        """

        DBG('sounding_archive: sounding_type = {0}, sounding_id = {1}, start_time={2}'
            .format(sounding_type, sounding_id, start_time))   
        self.instrument = instrument
        self.sounding_type = sounding_type
        self.sounding_id = sounding_id
        self.start_time = start_time
        #self.alt_res = alt_res
        #self.max_alt = max_alt
        self.requested_altitudes = requested_altitudes
        self.latitude = []
        self.longitude = []
        fallback_time=None
               
        if sounding_type == 'NOAA raob':
           #noaa raob soundings are returned in a Time_Z_Group
           
           self.soundings = raob.read_sounding_file(self.instrument\
                ,self.sounding_type,self.sounding_id,self.start_time\
                ,self.requested_altitudes)                                    
           fallback_time=raob.time_offset(start_time,1 if offset>0 else -1)
           #     ,self.alt_res,self.max_alt)
        elif self.sounding_type == 'time curtain':
           if self.sounding_id == 'raqms':
               #raqms soundings are provided at 1 second intervals at 35 model levels
               #with--profiles[time_index,model_levels_top_to_bottom] with time indexes
               #separated by one second--- sounding.times contains python datetimes and
               #soundings.model_level_alts contains the altitudes of the model layers.

               self.soundings = raqms.read_raqms_file(self.instrument,self.start_time)
               fallback_time=raqms.time_offset(start_time,1 if offset>0 else -1)
           else:
               raise RuntimeError, \
                        'request for unknown time curtain sounding id == %s  ' \
                        % self.sounding_id               
        else:
            raise RuntimeError, \
                        'request for unknown sounding type == %s  ' \
                        % self.sounding_type
        if self.profile(start_time,[],[],offset) is None and allow_fallback:
            tmp=sounding_archive(instrument,sounding_type,sounding_id,fallback_time,requested_altitudes,allow_fallback=False)
            if fallback_time>start_time:
                self.append(tmp)
            else:
                tmp.append(self)
                self.soundings=tmp.soundings

    def append(self,sounding_archive_tail):
        if self.soundings==None:
            self.soundings=sounding_archive_tail.soundings
        elif sounding_archive_tail.soundings==None:
            pass
        else:
            self.soundings.append(sounding_archive_tail.soundings)
    
    def profile(self,request_time,request_lat,request_long,offset=0):
        """returns a profile of temperature,pressure, dew_point, frost point at 
        time, lat, and long extracted from the sounding_archive. If request_lat and
        request_long are empty they are ignored
        request_time = python datetime for reqested sounding
        request_lat  = requested latitude for sounding--ignored if []
        request_lon  = requested longitude for sounding--ignored if []"""

        if self.soundings==None:
            return None

        print 'sounding_type= ',self.sounding_type,request_time
        if self.sounding_type == 'NOAA raob':
            temp_sounding=hau.Time_Z_Group()
             
            #find correct sounding profile out of archive file
            sounding=hau.selectByTime(self.soundings,request_time,offset=offset)
            if sounding is None:
                return None

            sounding.sounding_type=self.sounding_type
            sounding.sounding_id=sounding.station_id
            sounding.latitude=sounding.stalat
            sounding.longitude=sounding.stalong
            sounding.frost_points=cal_frost_point(sounding.dew_points)

            temp_sounding.type = self.sounding_type
            temp_sounding.times= sounding.times
            temp_sounding.altitudes=hau.Z_Array(sounding.altitudes)
            temp_sounding.temps= hau.TZ_Array(sounding.temps[np.newaxis,:])
            temp_sounding.pressures= hau.TZ_Array(sounding.pressures[np.newaxis,:])
            temp_sounding.dew_pts = hau.TZ_Array(sounding.dew_points[np.newaxis,:])
            temp_sounding.frost_pts = hau.TZ_Array(sounding.frost_points[np.newaxis,:])
            temp_sounding.wind_dir = hau.TZ_Array(sounding.wind_dir[np.newaxis,:])
            temp_sounding.wind_spd = hau.TZ_Array(sounding.wind_spd[np.newaxis,:])
            temp_sounding.wind_spd = hau.TZ_Array(sounding.wind_spd[np.newaxis,:])
            temp_sounding.station_id=sounding.station_id
            temp_sounding.top =sounding.top
        

            sounding.times=sounding.times[0]
            #sounding.expire_time=sounding.expire_time[0]
           
        elif self.sounding_type == "time curtain" \
                  and self.sounding_id == "raqms":

            temp_sounding=hau.Time_Z_Group()
            sounding = raqms.select_raqms_profile(self.soundings,request_time \
                       ,self.requested_altitudes,offset=offset)
            #                      ,self.max_alt,self.alt_res)
            if sounding==None:
                return None
            sounding.station_id =self.sounding_id

            temp_sounding.type = self.sounding_type
            temp_sounding.times= sounding.times
            temp_sounding.latitude=hau.T_Array(sounding.latitude)
            temp_sounding.longitude=hau.T_Array(sounding.longitude)
            temp_sounding.altitudes=hau.Z_Array(sounding.altitudes)
            temp_sounding.temps= hau.TZ_Array(sounding.temps[np.newaxis,:])
            temp_sounding.pressures= hau.TZ_Array(sounding.pressures[np.newaxis,:])
            temp_sounding.dew_pts = hau.TZ_Array(sounding.dew_points[np.newaxis,:])
            temp_sounding.frost_pts = hau.TZ_Array(sounding.frost_points[np.newaxis,:])
            temp_sounding.wind_dir = hau.TZ_Array(sounding.wind_dir[np.newaxis,:])
            temp_sounding.wind_spd = hau.TZ_Array(sounding.wind_spd[np.newaxis,:])
            temp_sounding.ext_total = hau.TZ_Array(sounding.ext_total[np.newaxis,:])
            temp_sounding.ext_salt = hau.TZ_Array(sounding.ext_salt[np.newaxis,:])
            temp_sounding.ext_dust = hau.TZ_Array(sounding.ext_dust[np.newaxis,:])
            temp_sounding.wind_spd = hau.TZ_Array(sounding.wind_spd[np.newaxis,:])
            temp_sounding.station_id=sounding.station_id
            temp_sounding.top =sounding.top

            sounding.times=sounding.times[0]

            #set up time to read in new sounding (in datetime)
            #sounding.expire_time = sounding.times+timedelta(seconds=5*60)
            #expire time can not be less then request time--raqms file must be missing soundings
            #set expire_time 5 minutes after request time
            #if sounding.expire_time <= request_time:
            #    print "****Warning----missing raqms sounding at ",request_time
            #    print "               using sounding from ",sounding.times
            #    sounding.expire_time = request_time + timedelta(seconds=5*60)
            
        #remove any negative pressure values
        temp_sounding.pressures[temp_sounding.pressures <0] = 0.0
        #remove NaN's from pressure and temperature values
        temp_sounding.pressures[np.isnan(temp_sounding.pressures)] = 0.0
        temp_sounding.temps[np.isnan(temp_sounding.temps)] = 1.0    
        return sounding


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
        if dew_point[j] < 273.15 and np.isfinite(dew_point[j]):
            try:
                frost_point[j] = brentq(delta_es, 100, 273, args=log_es[j])
            except ValueError:
                frost_point[j] = np.NaN
        else:
            frost_point[j] = np.NaN
    return frost_point

def delta_es(x, log_es):
    g = log_es - (-0.56745359e4 / x + 6.39225247 - 0.96778430e-2 * x
                  + 0.62215701e-6 * x ** 2 + 0.2074825e-8 * x ** 3
                  - 0.94840240e-12 * x ** 4 + 0.41635019e1 * np.log(x))
    return g

   
def es(temperature):
    """es =saturated water vapor pressure at the dry_bulb temp (mb)
       temperature (K-deg), es from Buck equation"""

    T=temperature-273.15
    #es = 6.112*np.exp((17.67*T)/(T+243.5))
    es = 6.1121*np.exp((18.678 - T/234.5)*(T/(257.14 + T)))
    return es

def cal_dew_point(relative_humidity,temperature):
    """dew point temperature K-deg"""
   
    vapor_pressure=relative_humidity*es(temperature)/100.0
    
    dew_point=243.5*np.log(vapor_pressure/6.112) \
               /(17.67 -np.log(vapor_pressure/6.112))
    dew_point=dew_point+273.15
    return dew_point

def hydrostatic_interp(s_pressures,heights,lidar_alts):
    """Use hydrostatic equation to interpolate sounding_pressures
       this version uses s_temps does not require interpolated temps
       s_temps,s_pressures and heights are the raw sounding values
       returns pressures at requested lidar_alts as l_pressures
       g = accel of gravity, 9.8 m/s/s
       R = gas constant"""
    g = 9.8 #m/s/s    accel of gravity
    R = 287.04 # J/Kg/K  Specific gas constant , dry air
    A = g/R
    i = 0  #lidar bin index
    k = 0  #sounding index
    l_pressures=np.zeros_like(lidar_alts)
    #step through lidar altitude bins
  
  
    k = 0
    for i in range(len(lidar_alts)):
        #insure that pressure defining heights bound lidar_alts[i]
        while k < len(s_pressures)-2 and lidar_alts[i] >= heights[k+1] :
            k = k + 1    
        CC = -np.log(s_pressures[k+1]/s_pressures[k])/(heights[k+1]-heights[k])
        #compute pressure at this alt assuming hydrostatic balance
        l_pressures[i]=s_pressures[k] * np.exp(-CC * (lidar_alts[i]-heights[k]))
    return l_pressures
