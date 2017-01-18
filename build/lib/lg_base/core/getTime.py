#!/usr/bin/env python
import os
import re
import datetime
import time

def get_time(year, month, day, working_dir):
    """This module finds the earliest time of applicable data available.
    
    This module works by accessing the directory listed in HSRL_CONFIG
    and then traversing that directory to find the earliest time that
    contains applicable data for the given date.
    """
    
    startTime = '235959'
    d = (working_dir +'/' + year + '/'+ month + '/' + day + '/raw')
    for root, dirs, files in os.walk(d):
        for name in files:
            f = str(os.path.join(root, name))
            if (f.endswith('.nc')):
                initTime = re.search('T\d{1,6}', f).group()
                initTime = initTime[1:]
                if (initTime < startTime):
                    startTime = initTime
    startTime = year + ' ' + month + ' ' + day + ' ' + startTime
    startTime = datetime.datetime.strptime(startTime, '%Y %m %d %H%M%S')
    return startTime
