
from datetime import datetime,timedelta
import dplkit.role.narrator
import dplkit.role.decorator
import dplkit.role.filter
import lg_base.core.array_utils as hau
import numpy
import copy
import time
import logging
import threading
import warnings
LOG = logging.getLogger(__name__)

def parse_timewindow(starttime,endtime,windowwidth=None,now=None):
    if now is None:
        now=datetime.utcnow()
    if starttime!=None:
        if endtime!=None:
            if windowwidth!=None:
                raise dplkit.role.librarian.AmbiguousQueryError,'Do not specify all start, end, and window'
            return {'starttime':starttime,'endtime':endtime,'windowwidth':endtime-starttime}
        if windowwidth:
            et=starttime+windowwidth
            if et>now:
                return{'starttime':et-windowwidth,'endtime':None,'windowwidth':windowwidth}
            return {'starttime':starttime,'endtime':None,'windowwidth':windowwidth}
        return {'starttime':starttime,'endtime':None,'windowwidth':now-starttime}
    if windowwidth is None:
        raise dplkit.role.librarian.AmbiguousQueryError,'Need at least starttime or windowwidth'
    if endtime:
        return {'starttime':endtime-windowwidth,'endtime':endtime,'windowwidth':windowwidth}
    return {'starttime':now-windowwidth,'endtime':None,'windowwidth':windowwidth}

@dplkit.role.decorator.exposes_attrs_of_field('framestream')
class TimeTrickle(object):
    """ wrapper to a frame stream that uses step() to get the current frame, or iterate depending on the time given. Typically used with sparse framestreams, like calibations and temperature profiles

    :param framestream: source frame stream
    :param timename: name of the start time in the frame
    :param endtimename: optional name of field containing end time. if is not provided, will use the timewidthname to find the end time
    :param timewidthname: optional name of field containing the frame's time width. if not provided, will iterate forward and use the next frame's start time as the end time
    """
    def __init__(self,framestream,timename='start',endtimename=None,timewidthname=None,getClosest=False):
        self.framestream=framestream
        self.iterator=iter(self.framestream)
        self.timename=timename
        self.endtimename=endtimename
        self.getClosest=getClosest
        #if endtimename==None and timewidthname==None:
        #    raise RuntimeError('Need an endtimename or timewidthname')
        self.timewidthname=timewidthname
        self.priorframe=self.__getone()
        self.nextframe=self.__getone()

    def __getone(self):
        try:
            return self.iterator.next()
        except StopIteration:
            return None

    def __getFieldNamed(self,frame,name):
        if isinstance(frame,dict):
            return frame[name]
        if hasattr(frame,name):
            return getattr(frame,name)
        print "couldn't find "+name,'in',frame
        return None

    def __timeForFrame(self,frame):
        return self.__getFieldNamed(frame,self.timename)

    def __endTimeForFrame(self,frame,nextframe=None):
        if self.getClosest:
            if nextframe is None:
                return datetime(2100,1,1,0,0,0)
            c=self.__timeForFrame(frame)
            n=self.__timeForFrame(nextframe)
            return c+timedelta(seconds=(n-c).total_seconds()/2.0)
        if self.endtimename==None:
            if self.timewidthname==None:
                if nextframe==None:
                    return datetime(2100,1,1,0,0,0)
                return self.__timeForFrame(nextframe)
            return self.__timeForFrame(frame)+self.__getFieldNamed(frame,self.timewidthname)
        return self.__getFieldNamed(frame,self.endtimename)

    @property
    def atEnd(self):
        """ :returns: true if there is no next frame
        """
        return self.nextframe==None

    @property
    def nextTime(self):
        """ :returns: time of start of next frame, or None if there is no next frame
        """
        if self.nextframe==None:
            return None
        return self.__timeForFrame(self.nextframe)

    def expireTime(self):
        if self.nextframe==None:
            return datetime(2100,1,1,0,0,0)
        return self.__endTimeForFrame(self.priorframe,self.nextframe)

    def step(self,atime):
        if self.atEnd:
            print '*** Warning at end of',str(self.framestream)
        count=0
        while (self.nextframe is not None and self.__timeForFrame(self.nextframe)<atime and \
                 self.__timeForFrame(self.nextframe)!=self.__timeForFrame(self.priorframe)) \
              or atime>=self.__endTimeForFrame(self.priorframe,self.nextframe):
            if self.nextframe is None:
                #print 'no next frame'
                return self.priorframe
            count=count+1
            #print count,'skipping ',self.__timeForFrame(self.priorframe),'current is',self.__timeForFrame(self.nextframe)
            self.priorframe=self.nextframe
            self.nextframe=self.__getone()
            #print 'read in',self.__timeForFrame(self.nextframe)
            if self.__timeForFrame(self.nextframe)==self.__timeForFrame(self.priorframe):
                break
        #print 'requested',atime,'got profile',self.__timeForFrame(self.priorframe),self.__endTimeForFrame(self.priorframe,self.nextframe),'after',count,'frames iterated. next to be',self.__timeForFrame(self.nextframe) if self.nextframe!=None else 'END'
        return self.priorframe

    def trickleGenerator(self,start,end):
        nextTime=start
        while nextTime<end:
            yield self.step(nextTime),min(end,self.expireTime())
            nextTime=self.expireTime()

    def __call__(self,atime):
        return self.step(atime)

#autprovides is only here for size of the frame.
@dplkit.role.decorator.autoprovidenested(nestedclasses=[hau.Time_Z_Group,hau.rs_xfer,dict])
class FrameCachedConcatenate(dplkit.role.filter.aFilter):
    """ DPL Filter object that accumulates frames, appending them to the first, and only yielding at the end. repeated uses will yield the cached final value

    :param stream: source stream
    """
    def __init__(self,stream,keylist=None,withoutConcatenate=None):#if asClass is none, will return a dictionary, otherwise will make a new class, and use setattr. if classdictparameter is
        super(FrameCachedConcatenate,self).__init__(stream)    
        self.stream=stream#this has provides and any exposed attributes
        self.keylist=keylist
        self.result=None
        self.isList=False if withoutConcatenate is None else withoutConcatenate
        self.lock=threading.Lock()

    @property
    def framecount(self):
        if not self.isList:
            return 1
        return len(self.result)

    def process(self):
        if self.result is None:
            with self.lock:
                if self.result is None:
                    isList=False
                    result=None
                    noFrames=True
                    for f in self.stream:
                        if self.keylist is not None and f is not None:
                            f=copy.deepcopy(f)
                            _f=f
                            if not isinstance(f,dict):
                                _f=vars(_f)
                            if True:
                                for k in _f.keys():
                                    if k.startswith("_"):
                                        continue
                                    if k in self.keylist:
                                        continue
                                    print 'dropping key',k
                                    del _f[k]
                        if result is None:
                            if hasattr(f,'append') and not self.isList:
                                result=copy.deepcopy(f)
                            else:
                                result=[f]
                                isList=True
                        else:
                            result.append(f)
                        noFrames=False
                    if noFrames:
                        result=[]
                        isList=True
                    self.result=result
                    self.isList=isList
        if self.isList:
            for f in self.result:
                yield f
        else:
            yield self.result


@dplkit.role.decorator.exposes_attrs_of_field('framestream')
class MeanNarrator(dplkit.role.narrator.aNarrator):
    """ Resample a Time_Z_Group compound framestream object based on time information

        :param framestream: source 
        :param basetime: start time of the resampling window
        :param endtime: end of the resampling window
        :param timeres: time resolution stepsize
        :param timesource: optional time source framestream, instead of using above parameters to generate one
        :param allow_nans: default false. if true, intervals that would have no data are filled appropriately. if false, these intervals are omitted entirely
    """
    def __init__(self,framestream,basetime=None,endtime=None,timeres=None,timesource=None,allow_nans=False):
        self.framestream=framestream
        self.base=basetime
        self.endtime=endtime
        self.step=timeres
        self.timesource=timesource
        self.allow_nans=allow_nans
        import lg_dpl_toolbox.dpl.TimeSource as TimeSource
        if self.timesource==None:
            self.timesource=TimeSource.TimeGenerator(start_time=basetime,end_time =endtime,time_resolution=timeres)

    def __repr__(self):
        return "MeanNarrator <"+repr(self.framestream)+">"

    def read(self):
        remainder=None
        base=None
        import lg_dpl_toolbox.dpl.TimeSource as TimeSource
        timesource=TimeSource.CompoundTimeGenerator(self.timesource)
        for f in self.framestream:
            if timesource.isDone:
                break
            if remainder==None:
                remainder=copy.deepcopy(f)
            else:
                remainder.append(f)
            t=getattr(remainder,remainder._timevarname)
            #print t.shape
            if t.shape[0]==0:
                continue
            requestedtimes=hau.T_Array(timesource.getBinsFor(starttime=base,endtime=t[-1]))
            if requestedtimes.size<2:
                continue
            lastTime=requestedtimes[-1]
            retarr=remainder
            remainder=hau.trimTimeInterval(retarr,lastTime,datetime(2200,1,1,0,0,0))
            retarr.trimTimeInterval(base,lastTime)
            retarr.hereGoneBinTimes(requestedtimes,allow_nans=self.allow_nans)
            print 'range',base,lastTime,'returning:',retarr,'remainder',remainder
            yield retarr
            base=lastTime
        if remainder!=None and timesource.end_time!=None:
            requestedtimes=hau.T_Array(timesource.getBinsFor(starttime=base,endtime=timesource.end_time,inclusive=True))
            remainder.hereGoneBinTimes(requestedtimes,allow_nans=self.allow_nans)
            if getattr(remainder,remainder._timevarname).shape[0]>0:
                print 'range',base,timesource.end_time,'returning:',remainder
                yield remainder

def getTimeStartTime(prior,current,next,frame):
    return dict(start=current,width=next-current)

def getTimeMiddleTime(prior,current,next,frame):
     ret=dict(start= prior+timedelta(seconds=(current-prior).total_seconds()/2))
     ret['width']=(current+timedelta(seconds=(next-current).total_seconds()/2))-ret['start']
     return ret

def getTimeEndTime(prior,current,next,frame):
    return dict(start=prior,width=current-prior)

defaultGetTimes=getTimeStartTime

@dplkit.role.decorator.exposes_attrs_of_field('framestream')
class TimeGinsu(dplkit.role.narrator.aNarrator):
    """ Slice Compound Frames by time into Simple Frames

        :param framestream: source 
        :param timefield: field name for time
        :param include: array of fields from the source to include in the resulting frame (default all)
        :param omit: array of fields from the source to omit from the resulting frame (default only timefield). applied to include
        :param onlyTime: if true, will only stream the start/width of the frames
        :param omitTime: if true (default), will omit the timefield if the omit paramter is not provided. this is typically done because interpolation doesn't work with datetime objects

    """
    def __init__(self,framestream,timefield,dtfield,include=None,omit=None,onlyTime=False,omitTime=True,getTimeFunction=None):
        self.framestream=framestream
        self.timefield=timefield
        self.dtfield=dtfield
        self.provides={}
        self.provides['start']={'shortname':'start','type':datetime}
        self.provides['width']={'shortname':'width','type':timedelta}
        self.getTimesFor=getTimeFunction or defaultGetTimes
        toinclude=include
        toomit=omit
        if toinclude==None:
            if onlyTime:
                toinclude=[]
            elif framestream.provides!=None:
                toinclude=framestream.provides.keys()
            else:
                toinclude=[]
        if toomit==None:
            if omitTime:
                toomit=[timefield]
            else:
                toomit=[]
        toomit.append('start')
        toomit.append('width')
        for x in toomit:
            if x in toinclude:
                toinclude.remove(x)
        if 'start' in toinclude or 'width' in toinclude:
            raise ValueError('collision of time info and frame')
        if framestream.provides!=None:
            for k in toinclude:#framestream.provides.items():
                #if issubclass(v['type'],hau.T_Array):
                if not k in framestream.provides:
                    raise ValueError('Framestream '+repr(framestream)+" doesn't provide the field "+k)
                self.provides[k]=framestream.provides[k]
        #print 'TimeGinsu',self,'has',self.provides

    def __repr__(self):
        return "TimeGinsu <"+repr(self.framestream)+">"

    def _getTimeAfter(self,times,time,prestart=-1,timesorder=None):
        tidx=prestart
        if timesorder is None:
            timesorder=range(len(times))
        rt=times[timesorder[prestart]] if (prestart>=0 and prestart<len(times)) else time
        for t in range(prestart+1,len(times)):
            if times[timesorder[t]] > time and (rt<=time or times[timesorder[t]]<rt):
                tidx=t
                rt=times[timesorder[t]]
            elif times[timesorder[t]] <= time:
                continue
            else:
                break
        if rt<=time:
            #print 'after',time,'not found'
            return (None,-1)
        #print 'after',time,'found at',tidx,'=',rt
        return (rt,tidx)

    def _frameAtIndex(self,f,times,tidx,timesorder=None):
        frame=hau.Time_Z_Group(like=f)
        if hasattr(frame,self.timefield):
            delattr(frame,self.timefield)
        #setattr(frame,'end',time)

        idxm=numpy.array([timesorder[tidx] if timesorder is not None else tidx])

        if False and hasattr(f,'start') and hasattr(f,'width'):
            pass
        elif self.dtfield==None or not hasattr(f,self.dtfield):
            #raise RuntimeError('NO DT')
            if self.dtfield!=None and not hasattr(f,self.dtfield):# and hasattr(frame.self.:
                print 'Missing dt field. tidx is ',tidx,' this is an error, since one was configured. field is',self.dtfield,'available is',repr(vars(f).keys())
                t=dict(start=times[1],width=timedelta(seconds=0))
                #raise RuntimeError('Missing dt field. tidx is ',tidx,' this is an error, since one was configured. field is',self.dtfield,'available is',repr(vars(f).keys()))
            else:
                t=self.getTimesFor(times[0],times[1],times[2],f)

            setattr(frame,'width',t['width'])#times[1]-times[0])
            setattr(frame,'start',t['start'])#times[0])
        else:
            w=getattr(f,self.dtfield)
            if isinstance(w,numpy.ndarray):
                w=w[idxm[0]]
            if not isinstance(w,timedelta):
                if numpy.isfinite(w):
                    w=timedelta(seconds=w)
                else:
                    w=timedelta(seconds=0)
            setattr(frame,'width',w)#times[1]-times[0])
            setattr(frame,'start',times[1])#times[0])


        assert(isinstance(frame.start,datetime))
        assert(isinstance(frame.width,timedelta))
        #print '%s %s' % (frame.start,frame.width)
        #lastvar=None
        #print idxm
        #print 'TimeGinsu',self,'using',self.provides
        for k in self.provides.keys():#,v in vars(f).items():
            #lastvar=k
            if k in ('start','width'):
                continue
            if not hasattr(f,k):
                LOG.warning("Frame is missing " +k)
                continue
            v=getattr(f,k)
            provides=self.provides
            if not 'type' in provides[k] and not isinstance(v,hau.Time_Z_Group):
                print 'Type not provided for '+k#+' with dictionary ',provides
            if isinstance(v,hau.Time_Z_Group):
                setattr(frame,k,self._frameAtIndex(v,times,tidx,timesorder,useProvides=provides[k]))
            elif 'type' in provides[k] and issubclass(provides[k]['type'],hau.T_Array):
                #print 'keeping',k
                if len(v.shape)==1:
                    setattr(frame,k,v[idxm])
                elif len(v.shape)>1:
                    tmp=[idxm]
                    for x in range(1,len(v.shape)):
                        tmp.append(slice(None,None))
                    setattr(frame,k,v[tuple(tmp)])
                else:
                    raise RuntimeError('Bad Shape in variable. '+k+'.shape=='+repr(v.shape)+' (dimension count = '+repr(len(v.shape))+')')
            else:
#                if k in self.provides:
                    #print 'missing',k
                setattr(frame,k,v)
                #print 'dropping non-time var',k,'type is',self.provides[k]['type'],type(v)
                pass
        return frame

    def read(self):
        priortime=None
        nexttime=None
        oldframe=None
        for f in self.framestream:
            if not hasattr(f,self.timefield) or getattr(f,self.timefield).size==0:
                continue
            times=getattr(f,self.timefield).copy()
            timesorder=times.argsort()
            if oldframe is None:#first run
                (time,tidx)=self._getTimeAfter(times,datetime(2000,1,1,0,0,0),timesorder=timesorder)
                nextidx=tidx
            else:
                nextidx=-1
            if self.dtfield!=None:
                while time is not None:
                    tf=self._frameAtIndex(f,(None,time,None),tidx,timesorder=timesorder)
                    #print 'yielding frame',tidx,time,tf#,vars(tf).keys()
                    yield tf
                    (time,tidx)=self._getTimeAfter(times,time,tidx,timesorder=timesorder)
            else:
                (nexttime,nextidx)=self._getTimeAfter(times,time,nextidx,timesorder=timesorder)
                if priortime==None and nexttime!=None:
                    priortime=time-(nexttime-time)
                while nexttime!=None:
                    if oldframe:#time is in the prior frame
                        if tidx>=0:
                            tf=self._frameAtIndex(oldframe,(priortime,time,nexttime),tidx,timesorder=oldorder)
                            #print 'yielding oldframe',tidx,time,tf
                            yield tf
                        oldframe=None
                    else:
                        tf=self._frameAtIndex(f,(priortime,time,nexttime),tidx,timesorder=timesorder)
                        #print 'yielding frame',tidx,time,tf#,vars(tf).keys()
                        yield tf
                    priortime=time
                    pidx=tidx
                    tidx=nextidx
                    time=nexttime
                    (nexttime,nextidx)=self._getTimeAfter(times,nexttime,nextidx,timesorder=timesorder)
                    if nexttime==None:
                        oldframe=copy.deepcopy(f)
                        oldorder=timesorder
                        break;
        if self.dtfield==None and oldframe and tidx>=0:
            tf=self._frameAtIndex(oldframe,(priortime,time,time+(time-priortime)),tidx,timesorder=oldorder)
            #print 'yielding oldframeEND',tidx,time,tf
            yield tf#self._frameAtIndex(oldframe,(priortime,time,time+(time-priortime)),tidx)

@dplkit.role.decorator.exposes_attrs_of_field('framestream')
@dplkit.role.decorator.exposes_attrs_in_chain(['start_time','end_time'])
class JustTime(dplkit.role.narrator.aNarrator):
    """ Extract just the Time part of the frame

        :param framestream: source 

    """
    def __init__(self,framestream,start_time=None,end_time=None):
        self.framestream=framestream
        self.provides=dict(start=dict(shortname='start',type=datetime),
                           width=dict(shortname='width',type=timedelta))
        self.start_time=start_time
        self.end_time=end_time

    def __repr__(self):
        return "JustTime <"+repr(self.framestream)+">"

    def read(self):
        ret=None
        for f in self.framestream:
            if not isinstance(f,dict):
                try:
                    f=vars(f)
                except:
                    continue
            ret=dict(start=f['start'],width=f['width'])
            if self.start_time is None:
                self.start_time=f['start']
            yield ret
        if self.end_time is None and ret is not None:
            self.end_time=ret['start']+ret['width']

@dplkit.role.decorator.exposes_attrs_of_field('framestream')
class Nearest(dplkit.role.filter.aFilter):
    """ Duplicate frames from a source to match a destination time rate

        :param framestream: source simple stream
        :param timestream: time and width simple stream
        :param maxgap: timedelta describing maximum tolerable gap in source data

    """
    def __init__(self,framestream,timestream,maxgap=None):
        super(Nearest,self).__init__(framestream)
        self.framestream=framestream
        self.timestream=timestream
        self.maxgap=maxgap

    def __repr__(self):
        return "Nearest <"+repr(self.framestream)+">"

    def diff(self,a,b):
        if a is None or b is None:
            return None
        astart=a['start']
        if 'width' in a:
            aend=astart+a['width']
        else:
            aend=astart
        bstart=b['start']
        if 'width' in b:
            bend=bstart+b['width']
        else:
            bend=bstart
        if (aend>=bstart and astart<=bstart) or (aend>=bend and astart<=bend):
            return 0.0
        if aend<bstart:
            return (bstart-aend).total_seconds()
        return (astart-bend).total_seconds()

    def process(self):
        prior=None
        priord=None
        next=None
        nextd=None
        source=iter(self.framestream)
        for _f in self.timestream:
            f=_f
            if not isinstance(f,dict):
                f=vars(f)
            while source is not None and next is None or nextd['start'] < f['start']:
                prior=next
                priord=nextd;
                try:
                    next=source.next()
                except StopIteration:
                    next=None
                    nextd=None
                    source=None
                else:
                    if isinstance(next,dict):
                        nextd=next
                    else:
                        nextd=vars(next)
            diffp=self.diff(priord,f)
            diffn=self.diff(nextd,f)
            print 'Diffs are',diffp,diffn
            if next is None:
                ret=prior
                retdiff=diffp
            elif prior is None or diffn<diffp:
                ret=next
                retdiff=diffn
            else:
                ret=prior
                retdiff=diffp
            if ret is None:
                ret=copy.deepcopy(_f)
            else:
                ret=copy.deepcopy(ret)
                if isinstance(ret,dict):
                    retd=ret
                else:
                    retd=vars(ret)
                retd['start']=f['start']
                retd['width']=f['width']
                if 'times' in retd:
                    retd['times'][:]=f['start']
                if 'delta_t' in retd:
                    retd['delta_t'][:]=f['width'].total_seconds()
                if self.maxgap is not None and retdiff>self.maxgap.total_seconds():
                    for k in [x for x in retd.keys()]:
                        if k.startswith((".",'_')):
                            continue
                        if k in ('start','width','times','delta_t'):
                            continue
                        del retd[k]
            yield ret


def main():
    from hsrl.dpl.dpl_hsrl import dpl_hsrl
    from lg_dpl_toolbox.dpl.NetCDFZookeeper import GenericTemplateRemapNetCDFZookeeper 
    import radar.dpl.MMCRMergeLibrarian as mmcr
    import radar.dpl.RadarFilters as rf
    import hsrl.dpl.dpl_artists as artists
    import functools

    import substruct
    import resample_altitude

    starttime=datetime(2006,12,24,22,0,0)
    endtime=datetime(2006,12,25,0,0,0)

    #stitcher=TimeStitch()
    #restr=SubstructRestractor('rs_mmcr')

    hsrllib=dpl_hsrl(instrument='ahsrl')
    hsrlnar=hsrllib(start_time_datetime=starttime,end_time_datetime=endtime,min_alt_m=0.0,max_alt_m=20000.0,timeres_timedelta=timedelta(seconds=30),altres_m=15)

    mmcrzoo=GenericTemplateRemapNetCDFZookeeper('eurmmcrmerge')
    mmcrlib=mmcr.MMCRMergeLibrarian('ahsrl','eurmmcrmerge.C1.c1.',zoo=mmcrzoo)
    mmcrnar=rf.RadarPrefilter(mmcrlib(start=starttime,end=endtime))
    #mmcrnar=mmcr.MMCRMergeBackscatterToReflectivity(resample_altitude.ResampleXd(TimeGinsu(substruct.SubstructExtractor(mmcrnar,None),'times',stitcherbase=stitcher),'heights',hsrlnar.getAltitudeAxis))
    mmcrnar=rf.RadarBackscatterToReflectivity(resample_altitude.ResampleXd(TimeGinsu(mmcrnar,'times'),'heights',hsrlnar.altitudeAxis.copy()))

    hsrlnarsplitter=substruct.SubstructBrancher(hsrlnar)
    #hsrlnar=TimeGinsu(substruct.SubstructExtractor(hsrlnar,'rs_inv',restractor=restr),'times',isEnd=False,stitchersync=stitcher)
    hsrlnar=TimeGinsu(hsrlnarsplitter.narrateSubstruct('rs_inv'),'times',isEnd=False)

    from dplkit.simple.blender import TimeInterpolatedMerge

    merge=TimeInterpolatedMerge(hsrlnar,[mmcrnar],allow_nans=True,channels=['heights','Reflectivity','MeanDopplerVelocity','Backscatter','SpectralWidth'])
    merge=substruct.Retyper(merge,functools.partial(hau.Time_Z_Group,timevarname='times',altname='heights'))
    #stitcher.setFramestream(merge)
    #restr.setFramestream(stitcher)
    #curs=restr
    curs=substruct.SubstructMerger('rs_inv',{
           'rs_mean':hsrlnarsplitter.narrateSubstruct('rs_mean'),
           'rs_raw':hsrlnarsplitter.narrateSubstruct('rs_raw'),
           'rs_inv':hsrlnarsplitter.narrateSubstruct('rs_inv'),
           'rs_mmcr':merge,
           'rs_init':hsrlnarsplitter.narrateSubstruct('rs_init'),
           'rs_static':hsrlnarsplitter.narrateSubstruct('rs_static'),
           'rs_Cxx':hsrlnarsplitter.narrateSubstruct('rs_Cxx'),
           'profiles':hsrlnarsplitter.narrateSubstruct('profiles',sparse=True),
        }
        ,hau.Time_Z_Group,{'timevarname':'times','altname':'heights'})
    artist=artists.dpl_images_artist(framestream=curs,
        instrument='ahsrl',max_alt=30*1000.0,processing_defaults=curs.hsrl_process_control,
        display_defaults='all_plots.json')
    curs=artist
    for frame in curs:
        print 'frame',frame
        print 'frame keys',vars(frame).keys()
        print 'rs_inv',frame.rs_inv
        print 'rs_mmcr',frame.rs_mmcr
        print 'RefShape',frame.rs_mmcr.Reflectivity.shape
        print 'MMCRTimes',frame.rs_mmcr.times.shape
        print 'invTimes',frame.rs_inv.times.shape
        print 'heights',frame.rs_mmcr.heights.shape
        for x in range(frame.rs_mmcr.times.shape[0]-1):
            if frame.rs_mmcr.times[x]==frame.rs_mmcr.times[x+1]:
                print 'dupe at',x,frame.rs_mmcr.times[x]

        #print type(frame),type(frame['Reflectivity']) if 'Reflectivity' in frame else 'no ref',type(frame['beta_a_backscat_par']) if 'beta_a_backscat_par' in frame else 'no backscat'
    time.sleep(5)

if __name__ == '__main__':
    main()
