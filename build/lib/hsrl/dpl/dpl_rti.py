#!/usr/bin/python
# -*- coding: utf-8 -*-

from hsrl.dpl.dpl_hsrl import dpl_hsrl

class dpl_rti(object):
    """ HSRL data stream inspired by DPL (see dpl_hsrl for actual DPL-based interfaces)

    Example: ::
    
        r = dpl_rti('gvhsrl',datetime.datetime(2011,8,11,0,0),datetime.datetime(2011,8,15,0,0), datetime.timedelta(seconds=5), datetime.timedelta(seconds=60*60*2),0,15000,50)
        for data in r:
            (data is the rs structure from the processing functions, maximum amount of data per loop is 'maxtimeslice')

    :param instrument:             hsrl id string (eg. 'ahsrl','gvhsrl','nshsrl','mf2hsrl').
    :param start_time_datetime:    datetime.datetime object first time to retrieve.
    :type start_time_datetime: datetime
    :param end_time_datetime:      datetime.datetime object last time to retrieve (None for open-ended, continue to get mini-slices between then-now, and now-now)
    :type end_time_datetime: datetime
    :param timeres_timedelta:      time resolution as a datetime.timedelta object (native would be timedelta(seconds=2.5)), or None to figure optimum image res from processing defaults
    :type timeres_timedelta: timedelta
    :param maxtimeslice_timedelta: datetime.timedelta object for amount of data retrieved at one time. window size (safe is 1 or 2 hours), or None for entire window at once
    :type maxtimeslice_timedelta: timedelta
    :param min_alt_m:              minimum altitude in meters to display
    :param max_alt_m:              maximum altitude in meters to display.
    :param altres_m:               altitude resolution in meters, or None to figure optimum image resfrom processing defaults

    """

    def __init__(self, instrument, start_time_datetime, end_time_datetime,timeres_timedelta=None,maxtimeslice_timedelta=None,min_alt_m=None,max_alt_m=None,altres_m=None,process_defaults=None,data_request=None):
        print '**** dpl_rti is deprecated. use dpl_hsrl'
        self.dpl=dpl_hsrl(instrument=instrument,maxtimeslice_timedelta=maxtimeslice_timedelta,process_defaults=process_defaults,data_request=data_request)
        self.timeslice=maxtimeslice_timedelta
        self.parms={'start_time_datetime':start_time_datetime,
                    'end_time_datetime':end_time_datetime,
                    'timeres_timedelta':timeres_timedelta,
                    'min_alt_m':min_alt_m,'max_alt_m':max_alt_m,
                    'altres_m':altres_m}
        import lg_base.core.array_utils as hau #import T_Array,Z_Array,TZ_Array,Time_Z_Group
        self.rs_static = hau.rs_xfer()
        #self.rs_static.processing_defaults = self.dpl.get_processing_defaults(process_defaults)[0]
        return

    def __iter__(self):
        return self.__gen()

    def __gen(self):
        catframe=None
        for frame in self.dpl(**self.parms):
            setattr(self,'rs_init',frame.rs_init)
            setattr(self,'rs_static',frame.rs_static)
            if self.timeslice!=None:
                yield frame
            else:
                if catframe==None:
                    catframe=frame
                else:
                    catframe.append(frame)
        if catframe!=None:
            yield catframe

if __name__ == '__main__':
    dplhsrl=dpl_rti('bagohsrl',datetime.datetime(2012,4,10,2,0,0),datetime.datetime(2012,4,10,12,0,0),datetime.timedelta(seconds=60),datetime.timedelta(seconds=60*30),0,20000,150)
    for i in dplhsrl:
        print i