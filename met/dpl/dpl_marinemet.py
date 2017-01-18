#!/usr/bin/python
# -*- coding: utf-8 -*-

import datetime
import numpy as np
import dplkit.role.narrator
import dplkit.role.librarian
import json

import traceback
from time import sleep

import lg_base.core.array_utils as hau
import dplkit.role.decorator

class dpl_marinemet(dplkit.role.librarian.aLibrarian):
    """ Marinemet data DPL framestream

    Example: ::

      r = dpl_hsrl(instrument='ahsrl')
      for data in r(datetime.datetime(2011,8,11,0,0),datetime.datetime(2011,8,15,0,0), timeres_timedelta=datetime.timedelta(seconds=5), min_alt_m=0,max_alt_m=15000,altres_m=50):
            (data is the rs structure from the processing functions, maximum amount of data per loop is 'maxtimeslice')

    :param instrument: hsrl id string (eg. 'ahsrl','gvhsrl','nshsrl','mf2hsrl').
    :param maxtimeslice_timedelta: datetime.timedelta object for amount of data retrieved at one time. window size (safe is 1 or 2 hours), or None for entire window at once
    :type maxtimeslice_timedelta: datetime.timedelta
    :param process_control:        process control structure or json filename (contains corrections and process defaults)
    :param data_request:           data request string   
    :param filetype: HSRL raw filetype to use (None for all (default), 'data' for only data files)

    """

    def __init__(self, instrument,*args, **kwargs):
        super(self.__class__,self).__init__(None)
        self.instrument=instrument
        import lg_base.core.open_config as oc
        self.oc=oc
        import lg_base.core.json_config as jc
        self.jc=jc
        from lg_base.core.locate_file import locate_file
        self.locate_file=locate_file
        import lg_base.core.array_utils as hau #import T_Array,Z_Array,TZ_Array,Time_Z_Group
        self.hau=hau
        from lg_dpl_toolbox.dpl.NetCDFZookeeper import GenericTemplateRemapNetCDFZookeeper 
        import MarinemetLibrarian as mm
        if instrument in ('marinemet','mf2hsrl','mf2marinemet','magmarinemet'):            
            self.instrument='marinemet'
            self.zoo=GenericTemplateRemapNetCDFZookeeper('marinemet')
            self.lib=mm.MarinemetLibrarian('mf2hsrl',['magmarinemetM1'],zoo=self.zoo)
            self.instrumentbase='mf2hsrl'
        elif instrument in ('met','mf2met','tmpmet'):       
            self.instrument='met'
            self.zoo=GenericTemplateRemapNetCDFZookeeper('met')
            self.lib=mm.MarinemetLibrarian('mf2hsrl',['tmpmetM1'],zoo=self.zoo)
            self.instrumentbase='mf2hsrl'
        else:
            raise RuntimeError('Unknown met source '+instrument)


    def __repr__(self):
        return 'DPL Marinemet Librarian (instrument="%s")' % (self.instrument)

    def search(self, start_time_datetime=None, end_time_datetime=None,reverse_padding=None,timeres_timedelta=None,window_width_timedelta=None,forimage=None,inclusive=False,timesource=None,allow_nans=False,*args, **kwargs):
        """
        :param start_time_datetime: start time (optional)
        :type start_time_datetime: datetime.datetime
        :param end_time_datetime: end time (optional) if unspecified, will continue to return frames thru now, unending
        :type end_time_datetime: datetime.datetime
        :param reverse_padding: (optional)in the case of reading up to the current time, this timedelta is always subtracted from the current time to get the most recent moment to read
        :type reverse_padding: datetime.timedelta
        :param timeres_timedelta: (optional) time resolution, or None to optimized for images
        :type timeres_timedelta: datetime.timedelta
        :param window_width_timedelta:   used with start or end time (not both) to determine preferred return width. if unspecified, will imply the other unspecified, or from now if neither
        :type window_width_timedelta: datetime.timedelta
        :param corr_adjusts: correction adjustments. if unspecified, will use default   
        :param block_when_out_of_data: (optional) if unspecified or True, will block for 'timeres_timedelta' until trying for more frames. if False, yield None until more data is available
        :param forimage: (optional) if provided, is a dict(x=##,y=##) containing preferred x and y pixel count for an image. if no resolutions are provided, these are used to create an optimal one. if not provided, lacking resolutions are native
        :param inclusive: if true, will expand window to ensure including any data that intersects the requested times
        """
        if len(args):
            print 'Unused dpl_marinemet.search args = ',args
        if len(kwargs):
            print "Unused dpl_marinemet.search kwargs = ",kwargs

        # altitude configuration notes: min and max alt are required. resolution is optional, in which case the process_control will determine pixel count, and thus resolution

        # time configuration notes: there are many modes to this operation, and two layers
        #   calibration layer: possible parameters are start, end, and window width. (implemented by dpl_calibration.parse_timewindow(), which will return in all cases start and width, and end if will terminate)
        #      specification groups possible: start+end, start+window,start,end+window,window (in order of priority)
        #      start and end specify the range of calibrations to stream
        #      if only one is specified, window is used to roll forward or back the one specified to make the other
        #         if window is not specified, it goes from start to now (defining the window)
        #      if neither are specified, start is set to now-window
        #      if end is not specified, it keeps rolling forward without termination. if it is, it WILL go to that point (even in future) and terminate
        #         if start and window are specified, and start+window is past now, start will be ignored
        #   processing layer:
        #      extra parameter of maxtimeslice specifies the largest volume of actual data to be returned (may be violated if no data is available, and filler piles up)
        #         if it is not specified, natural flow is not interrupted for the sake of data volume
        #      here, it only cares if timeres is not specified.  If it is not specified, it will follows the same steps as calibration to find the preferred window size, and use process_control to determine pixel count, and resolution

        prms={}
        from lg_dpl_toolbox.filters import time_frame,resample_altitude
        if reverse_padding!=None:
            prms['now']=datetime.datetime.utcnow()-reverse_padding
        d=time_frame.parse_timewindow(start_time_datetime,end_time_datetime,window_width_timedelta,**prms)

        if forimage!=None:
            if not isinstance(forimage,dict):
                import lg_base.core.canvas_info as ci
                forimage=ci.load_canvas_info()['canvas_pixels']
            if timeres_timedelta==None:
                timeres_timedelta=datetime.timedelta(seconds=d['windowwidth'].total_seconds()/forimage['x'])

        nar=self.lib(start=d['starttime'],end=d['endtime'])
        #FIXME too constant, and hardcoded
        if timeres_timedelta!=None and nar.metNativeTimeResolution>timeres_timedelta:
            #dropping time resolution to 'pure native'
            timeres_timedelta=None

        if inclusive:
            #remake with padding
            padAmount=(2*timeres_timedelta) if timeres_timedelta!=None else (5*nar.metNativeTimeResolution)
            if d['starttime']:
                d['starttime']-=padAmount
            if d['endtime']:
                d['endtime']+=padAmount
            nar=self.lib(start=d['starttime'],end=d['endtime'])
        elif timeres_timedelta:
            nar=self.lib(start=d['starttime']-timeres_timedelta,end=d['endtime']+timeres_timedelta)
            
        if timesource:
            nar=time_frame.MeanNarrator(nar,timesource=timesource,allow_nans=allow_nans)
        elif timeres_timedelta:
            nar=time_frame.MeanNarrator(nar,basetime=d['starttime'],timeres=timeres_timedelta,endtime=d['endtime'],allow_nans=allow_nans)
        import lg_dpl_toolbox.dpl.TimeSource as TimeSource
        nar=TimeSource.AddPseudoDeltaT(nar,'times','delta_t')

        return nar

if __name__ == '__main__':
    dplhsrl=dpl_hsrl('bagohsrl',datetime.timedelta(seconds=60*60*2))
    for i in dplhsrl(start_time_datetime=datetime.datetime(2013,1,10,0,0,0),end_time_datetime=datetime.datetime(2013,1,10,12,0,0),min_alt_m=0,max_alt_m=15000):
        print 'step'
