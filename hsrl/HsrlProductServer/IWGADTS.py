from datetime import datetime
from PyQt4.QtCore import *
from PyQt4.QtNetwork import *

import dateutil.parser

class IWGADTS:
    # Minimum list of variables for IWGADTS packet
    __MIN_VAR_LIST__ = [
        'date_time',
        'lat',
        'lon',
        'gps_msl_alt',
        'wgs_84_alt',
        'press_alt',
        'radar_alt',
        'grnd_spd',
        'true_airspeed',
        'indicated_airspeed',
        'mach_number',
        'vert_velocity',
        'true_hdg',
        'track',
        'drift',
        'pitch',
        'roll',
        'side_slip',
        'angle_of_attack',
        'ambient_temp',
        'dew_point',
        'total_temp',
        'static_press',
        'dynamic_press',
        'cabin_press',
        'wind_speed',
        'wind_dir',
        'vert_wind_speed',
        'solar_zenith',
        'sun_elev_ac',
        'sun_az_ground',
        'sun_az_ac'
    ]

    def __init__(self):
        '''
        Instantiate an IWGADTS sample.
        '''
