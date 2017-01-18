import dplkit.role.librarian
import dplkit.role.narrator
import dplkit.role.filter
from datetime import datetime,timedelta
import numpy
import dplkit.role.decorator
import lg_base.core.array_utils as hau
import os,calendar
import lg_dpl_toolbox.dpl.ARMNarrator as ARMNarrator
import re

def expand(t,td):
    if t is None:
        return t
    return t+td

@dplkit.role.decorator.exposes_attrs_in_chain(['ramanNativeTimeResolution','ramanType'])
class RamanNarrator(ARMNarrator.ARMNarrator):
    """ Narrator Initialzation for Raman.  Shouldn't be created, should only be used by the RamanLibrarian

            :param host: source librarian object
            :param basedir: base directory to find Raman data
            :param dataprefix: filename prefixes for Raman data
            :param start: start time
            :param end: end time
            :param zoo: zookeeper, if narrator should return data, not file info
            :param kwargs: addtional parameters are passed to thezookeeper on read, if provided

    """

    @property
    def ramanNativeTimeResolution(self):
        if self.timeres==None:
            if False:#self.zoo!=None:
                for fr in self:
                    if not hasattr(fr,'times') or fr.times.size<2:
                        continue
                    self.timeres=timedelta(seconds=(fr.times[-1]-fr.times[0]).total_seconds()/(fr.times.size-1))
                    break
            if self.timeres==None:
                self.timeres=timedelta(seconds=3600)
        return self.timeres

    @property
    def ramanType(self):
        return self._ramanType

    def __init__(self,host,basedir,datatype,start,end,platform=None):
        super(RamanNarrator,self).__init__(basedir,"Raman Lidar",(platform or '.*')+'rlprof'+datatype.split('_')[0],'raman',start,end,platform=platform)
        self.host=host
        self.timeres=None
        self._ramanType='raman_'+datatype

    def fixValue(self,var,attrs,globalattrs,attrname,replaceval):
        nv=None
        if attrname in attrs:
            nv=attrs[attrname]
        elif attrname in globalattrs:
            nv=globalattrs[attrname]
        if nv is None:
            return
        elif isinstance(nv,basestring):
            nv=float(nv)
        try:
            var[var==nv]=replaceval
        except TypeError:
            pass

    def preYield(self,x,attrs,found):
        if not hasattr(x,'heights') or not hasattr(x,'times') or x.times.size==0:
            return False
        
        if 'merge' in self.ramanType:
            basealt=0 #don't correct yet for merge
        elif hasattr(x,'alt'):#platform height from sealevel
            basealt=x.alt
        elif "altitude [m AMSL]" in attrs:
            basealt=float(attrs["altitude [m AMSL]"])
            setattr(x,'alt',basealt)
            print 'Got an altitude from attributes', basealt,x.alt,self.ramanType
        else:
            raise RuntimeError("Can't find MSL altitude")
        dualmerge=hasattr(x,'heights') and hasattr(x,'heights_low')
        if dualmerge:
            print 'merging low into high'
            #print 'high =',x.heights
            #print 'low  =',x.heights_low
            mergemap=[]
            merged=[0,0]
            hlowcount=x.heights_low.size
            for lidx,alt in enumerate(x.heights_low):
                idxs=np.where(x.heights==alt)[0]
                if idxs.size==0:
                    #print 'height_low',alt,'at offset',lidx,'not in high'
                    #raise RuntimeError('missing alt')
                    mergemap.append(None)
                    merged[1]=merged[1]+1
                else:
                    if idxs.size>1:
                        print 'heigh_low',alt,'has',len(idxs),'matches'
                    #print 'merging in alt',alt,'from low',lidx,'to high',idxs,'(len',idxs.size,')'
                    mergemap.append(int(idxs.ravel()[0]))
                    merged[0]=merged[0]+1
            if merged[1]>0:
                print 'Successfully merged =',merged[0],'   unmerged and dropped =',merged[1]
            delattr(x,'heights_low')

        for n,v in vars(x).items():
            #print n,v
            if n.startswith('_'):
                continue
            #v=v[altmask]
            fd = found[n] if n in found else {}
            self.fixValue(v,fd,attrs,'missing_value',np.NAN)
            self.fixValue(v,fd,attrs,'missing-data',np.NAN)
            self.fixValue(v,fd,attrs,'nan_value',np.NAN)
            units = '' if 'units' not in fd else fd['units']
            if '/' in units:
                sp=units.split('/')
                numer=sp[0]
                denom=sp[1]
            else:
                numer=units
                denom='1'
            if 'km' in denom:# THIS DOES ASSUME ANY KM/M conversions don't stick, and have to be redone here
                v/=1000.0
            if 'km' in numer:
                v*=1000.0
            if 'depol' in n and 'counts' not in n and 'uncert' not in n and (len(units)<2 or units=="unitless") :
                v*=100.0
            if n=='heights' and 'merge' not in self.ramanType:
                setattr(x,'_altitudevarname','altitudes')
                setattr(x,"altitudes",hau.Z_Array(v+basealt))
                continue
            if (n.endswith('_low') or '_low_' in n) and len(v.shape)>1 and v.shape[1]==hlowcount and dualmerge:
                tmp=v
                hname=n.replace('_low','_high')
                if not hasattr(x,hname):
                    hname=n.replace('_low','')
                if not hasattr(x,hname):
                    print "Can't find high equivalent to "+n
                v=getattr(x,hname).copy()
                def subslice(chunk=None):
                    ret=[slice(None),slice(chunk)]
                    for x in range(2,len(v.shape)):
                        ret.append(slice(None))
                    return tuple(ret)
                v[subslice()]=np.NAN
                for si,i in enumerate(mergemap):
                    if i is not None:
                        v[subslice(i)]=tmp[subslice(si)]
                setattr(x,n,v)
        dt=x.times.copy()
        dt.summode="sum"
        setattr(x,'width',dt)
        setattr(x,'start',x.times.copy())
        if dt.size>1:
            for i in range(dt.size-1):
                dt[i]=(x.times[i+1]-x.times[i]).total_seconds()
            dt[-1]=dt[-2]
            for i,dtv in enumerate(dt):
                x.start[i]-=timedelta(seconds=dtv/2.0)
        return True


class RamanLibrarian(dplkit.role.librarian.aLibrarian):
    """ Librarian Initialzation for Raman
  
            :param siteid: source site id for the data source. typically an hsrl instrument or base directory
            :param dataprefix: filename prefixes for Raman data
            :param zoo: zookeeper, if narrator should return data, not file info
    """
    def __init__(self, siteid,datatype,zoo=None):
        super(self.__class__,self).__init__()
        self.basedir=ARMNarrator.getBasedir(siteid)
        if datatype.startswith('rlprof'):
            datatype=datatype.replace('rlprof','')
        self.datatype=datatype
        self.zoo=zoo

    def search(self,start,end,*args,**kwargs):
        """ Librarian Generator function
        extra parameters given here will be passed to the returned narrator's init
        """
        ret=RamanNarrator(self,self.basedir,self.datatype,start,end)
        zoo=self.zoo
        if zoo is None and not kwargs.pop('filenames',False):
            from lg_dpl_toolbox.dpl.NetCDFZookeeper import GenericTemplateRemapNetCDFZookeeper 
            zoo=GenericTemplateRemapNetCDFZookeeper('raman_'+self.datatype)
        if zoo is not None:
            ret=ARMNarrator.ARMFileNarrator(ret,zoo,preYield=ret.preYield,*args,**kwargs)
        return ret

#!/usr/bin/python
# -*- coding: utf-8 -*-

#import datetime
import numpy as np
import dplkit.role.narrator
import dplkit.role.librarian
import json

import traceback
from time import sleep

import lg_base.core.array_utils as hau
import dplkit.role.decorator
import radar.core.radar_masking as rm
import radar.core.radar_time_filtering as rtf

from lg_dpl_toolbox.filters.time_frame import parse_timewindow
from lg_dpl_toolbox.dpl.dpl_configured_callable import dpl_configured_callable,dpl_expose_attributes_filter
from lg_dpl_toolbox.filters.dpl_rolling_window_filter import dpl_rolling_window_filter

class dpl_raman(dplkit.role.librarian.aLibrarian):
    """ Raman data DPL framestream

    Example: ::

      r = dpl_radar(instrument='mmcr')
      for data in r(datetime.datetime(2011,8,11,0,0),datetime.datetime(2011,8,15,0,0), timeres_timedelta=datetime.timedelta(seconds=5), min_alt_m=0,max_alt_m=15000,altres_m=50):
            (data is the rs structure from the processing functions, maximum amount of data per loop is 'maxtimeslice')

    :param instrument: instrument id string (eg. 'mmcr','magkazrge','magkazrmd','nsakazrge','nsakazrmd'). also supports naming the co-located HSRL to get a default radar source
    :param process_control:        process control structure or json filename (contains corrections and process defaults)

    """

    def __init__(self, site,datatype,process_control=None):#,*args, **kwargs):
        super(self.__class__,self).__init__(None)
        if datatype.startswith('rlprof'):
            datatype=datatype.replace('rlprof','')
        self.instrument='rlprof'+datatype
        from lg_base.core.locate_file import locate_file
        #self.process_control_file=locate_file(process_control or 'raman_processing_defaults.json',systemOnly=process_control==None)
        import lg_base.core.open_config as oc
        self.oc=oc
        import lg_base.core.json_config as jc
        self.jc=jc
        #import hsrl.utils.hsrl_array_utils as hau #import T_Array,Z_Array,TZ_Array,Time_Z_Group
        #self.hau=hau
        #from lg_dpl_toolbox.dpl.NetCDFZookeeper import GenericTemplateRemapNetCDFZookeeper 
        #import RamanFilters as rf
        #self.rf=rf
        #self.callableargs=kwargs
        #self.instrument=instrument
        self.lib=RamanLibrarian(site,datatype)#,zoo=self.zoo)
        self.instrumentbase=site


    def __repr__(self):
        return 'DPL Radar Librarian (instrument="%s")' % (self.instrument)

    def search(self, start_time_datetime=None, end_time_datetime=None,reverse_padding=None,timeres_timedelta=None,min_alt_m=None,max_alt_m=None,\
        altres_m=None,window_width_timedelta=None,forimage=None,inclusive=False,timesource=None,allow_nans=False,*args, **kwargs):
        """
        :param start_time_datetime: start time (optional)
        :type start_time_datetime: datetime.datetime
        :param end_time_datetime: end time (optional) if unspecified, will continue to return frames thru now, unending
        :type end_time_datetime: datetime.datetime
        :param reverse_padding: (optional)in the case of reading up to the current time, this timedelta is always subtracted from the current time to get the most recent moment to read
        :type reverse_padding: datetime.timedelta
        :param timeres_timedelta: (optional) time resolution, or None to optimized for images
        :type timeres_timedelta: datetime.timedelta
        :param min_alt_m: minimum altitude in meters
        :param max_alt_m: maximum altitude in meters
        :param altres_m: (optional) altitude resolution
        :param window_width_timedelta:   used with start or end time (not both) to determine preferred return width. if unspecified, will imply the other unspecified, or from now if neither
        :type window_width_timedelta: datetime.timedelta
        :param corr_adjusts: correction adjustments. if unspecified, will use default   
        :param block_when_out_of_data: (optional) if unspecified or True, will block for 'timeres_timedelta' until trying for more frames. if False, yield None until more data is available
        :param forimage: (optional) if provided, is a dict(x=##,y=##) containing preferred x and y pixel count for an image. if no resolutions are provided, these are used to create an optimal one. if not provided, lacking resolutions are native
        :param inclusive: if true, will expand window to ensure including any data that intersects the requested times
        """
        if len(args):
            print 'Unused dpl_radar.search args = ',args
        if len(kwargs):
            print "Unused dpl_radar.search kwargs = ",kwargs

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
        if reverse_padding!=None:
            prms['now']=datetime.utcnow()-reverse_padding
        d=parse_timewindow(start_time_datetime,end_time_datetime,window_width_timedelta,**prms)

        if forimage!=None:
            if not isinstance(forimage,dict):
                import lg_base.core.canvas_info as ci
                forimage=ci.load_canvas_info()['canvas_pixels']
            if altres_m==None:
                if min_alt_m==None:
                    min_alt_m=0
                if max_alt_m==None:
                    max_alt_m=30
                altres_m=(max_alt_m-min_alt_m)/forimage['y']
            if timeres_timedelta==None:
                timeres_timedelta=timedelta(seconds=d['windowwidth'].total_seconds()/forimage['x'])
        from lg_dpl_toolbox.filters import time_frame,resample_altitude

        def extralibparms():#scope abuse
            if 'merge' in self.instrument:
                return dict()
            return dict(firstaltitude=min_alt_m-altpad,lastaltitude=max_alt_m+altpad)

        altpad=(2*altres_m) if altres_m!=None else 0
        ramannar=self.lib(start=d['starttime'],end=d['endtime'],**extralibparms())
        #FIXME too constant, and hardcoded
        if timeres_timedelta!=None and ramannar.ramanNativeTimeResolution>timeres_timedelta:
            #dropping time resolution to 'pure native'
            timeres_timedelta=None

        if inclusive:
            #remake with padding
            padAmount=(2*timeres_timedelta) if timeres_timedelta!=None else None
            if padAmount is None or padAmount<(2*ramannar.ramanNativeTimeResolution):
                padAmount=2*ramannar.ramanNativeTimeResolution
            if d['starttime']:
                d['starttime']-=padAmount
            if d['endtime']:
                d['endtime']+=padAmount
            altpad=(2*altres_m) if altres_m!=None else 50
            ramannar=self.lib(start=d['starttime'],end=d['endtime'],**extralibparms())
        elif timeres_timedelta:#just pad the input data. don't mess with the mean narrator below
            ramannar=self.lib(start=d['starttime']-timeres_timedelta,end=d['endtime']+timeres_timedelta,**extralibparms())


        #ramannar=self.rf.RadarPrefilter(ramannar)
        if timesource is not None:
            ramannar=time_frame.MeanNarrator(ramannar,timesource=timesource,allow_nans=allow_nans)
        elif timeres_timedelta is not None:
            ramannar=time_frame.MeanNarrator(ramannar,basetime=d['starttime'],timeres=timeres_timedelta,endtime=d['endtime'],allow_nans=allow_nans)
        if 'merge' not in self.instrument and altres_m is not None:
            requested_altitudes=hau.Z_Array(np.arange(min_alt_m,max_alt_m+altres_m*0.1,altres_m))
            ramannar=resample_altitude.ResampleXd(ramannar,'altitudes',requested_altitudes,left=np.NaN,right=np.NaN)

        return ramannar

def main():
    import sys
    sp=RamanLibrarian(sys.argv[1],sys.argv[2])
    for x in sp(datetime(2015,7,10,1,0,0),None):#,filenames=True):
        print x
        #print x.altitudes

if __name__ == '__main__':
    main()
