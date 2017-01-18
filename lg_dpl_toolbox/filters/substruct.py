
import dplkit.role.narrator
import dplkit.role.blender
import dplkit.role.artist
import dplkit.role.filter
import copy
import numpy
import lg_base.core.array_utils as hau
import dplkit.role.decorator
from datetime import timedelta,datetime
import time
import os
import traceback

import threading,Queue
import weakref
import multiprocessing
import logging
LOG = logging.getLogger(__name__)

VERBOSE = LOG.debug
import functools
from collections import OrderedDict
if not hasattr(dplkit.role.decorator,'is_stateful'):
    def nulldecorator(cls):
        print 'WARNING class ',cls,' uses is_stateful from dplkit. upgrade to dplkit 0.3.0 soon'
        return cls
    setattr(dplkit.role.decorator,'is_stateful',nulldecorator)

#@dplkit.role.decorator.autoprovidenested(nestedclasses=[hau.Time_Z_Group,hau.rs_xfer,dict])
class TupleNarrator(dplkit.role.narrator.aNarrator):
    """ Narrator for a list of objects

    :param content: framestream content as a tuple or list
    """
    def __init__(self,content,provides,*args,**kwargs):
        super(TupleNarrator,self).__init__()    
        self.content=content
        self.provides=provides or dplkit.role.decorator.describe(self.content[0],*args,**kwargs)

    def read(self):
        for f in self.content:
            yield f

class NestedFrame(dplkit.role.filter.aFilter):
    """ Filter To a nested frame

    :param stream: framestream source
    :param field: fieldname to keep
    """
    def __init__(self,stream,field,provides=None):#if asClass is none, will return a dictionary, otherwise will make a new class, and use setattr. if classdictparameter is
        super(NestedFrame,self).__init__(stream)    
        self.stream=stream#this has provides and any exposed attributes
        self.fieldlist=field.split('.')
        pv=provides or stream.provides
        for f in self.fieldlist:
            if f not in pv:
                pv={}
                err= 'WARNING: Nested frame '+f+" doesn't exist in "+repr(stream)+'. in path '+('.'.join(self.fieldlist))+' have ('+(','.join(pv.keys()))+')'
                print err
                print 
                raise KeyError(err)
            pv=pv[f]
        self.provides=pv

    def __repr__(self):
        return 'NestedFrame<'+repr(self.stream)+'['+('.'.join(self.fieldlist))+']'+'>'

    def process(self):
        for f in self.stream:
            for fn in self.fieldlist:
                if not isinstance(f,dict):
                    try:
                        f=vars(f)
                    except:
                        f=None
                        break
                if not fn in f:
                    print 'Missing ',fn,' in frame'
                    #raise RuntimeError('Missing Frame '+fn)
                    f=None
                    break
                else:
                    f=f[fn]
            yield f

class NestingArtist(dplkit.role.artist.aArtist):
    #artist will deadend
    def __init__(self,framestream,nesting,initfunc):
        super(NestingArtist,self).__init__(framestream)
        self.fr=framestream
        self.nesting=nesting
        self.initfunc=initfunc

    def render(self):
        q=Queue.Queue(1)
        subframestream=iq.QueueAsIterable(q)
        subframestream=NestedFrame(subframestream,self.nesting,self.framestream.provides)
        subframestream=self.initfunc(subframestream)
        myiter=None
        for x in self.framestream:
            q.put(x)
            if myiter is None:
                myiter=iter(subframestream)
            myiter.next()
            yield x

class DropFrameContent(dplkit.role.filter.aFilter):
    """ Filter to drop unwanted frame content

    :param stream: framestream source
    :param fieldlist: list of fieldnames to drop from the frame and the provides
    """
    def __init__(self,stream,fieldlist):#if asClass is none, will return a dictionary, otherwise will make a new class, and use setattr. if classdictparameter is
        super(DropFrameContent,self).__init__(stream)    
        self.stream=stream#this has provides and any exposed attributes
        self.fieldlist=fieldlist
        self.provides=copy.copy(stream.provides)
        if isinstance(fieldlist,basestring):
            self.fieldlist=[fieldlist]
        for f in self.fieldlist:
            if f in self.provides:
                del self.provides[f]

    def process(self):
        for f in self.stream:
            f=copy.copy(f)
            for field in self.fieldlist:
                if hasattr(f,field):
                    delattr(f,field)
            yield f


class Nester(dplkit.role.filter.aFilter):
    def __init__(self,fs,name,f=None):
        super(Nester,self).__init__(fs)
        self.fs=fs
        self.f=f or {}
        self.provides={name:fs.provides}
        self.name=name

    def process(self):
        for f in self.fs:
            dupe=copy.deepcopy(self.f)
            if isinstance(dupe,dict):
                dupe[self.name]=f
            else:
                setattr(dupe,self.name,f)
            yield dupe


@dplkit.role.decorator.is_stateful
class NestingCompositer(dplkit.role.blender.aBlender):
    """ blender object to attach source streams' frames onto a primary stream's frame

    for every frame in a primary framestream, retrieve one frame from every secondary, and merge it naively into the primary's frame

    :param primary: primary source.
    :param addonstreams: dictionary of streams, keyed to where the name of the field the stream's frame will attach to the primary
    """
    def deepset(self,frame,content,field):
        if isinstance(field,basestring):
            self.deepset(frame,content,field.split('.'))
        else:
            v=frame
            if not isinstance(frame,dict):
                v=vars(frame)
            if len(field)==1:
                #if field[0] in v:
                #    print 'Nesting is replacing a subframe',v[field[0]],'with new subframe',content
                v[field[0]]=content
            else:
                if field[0] not in v:
                    v[field[0]]=OrderedDict()
                self.deepset(v[field[0]],content,field[1:])


    def __init__(self,primary,addonstreams):
        super(NestingCompositer,self).__init__([primary] + [addonstreams[x] for x in addonstreams.keys()])        
        self.primary=primary
        self.sechost=addonstreams
        self.noFramePlaceholder='NOFRAME'
        self.provides=copy.deepcopy(self.primary.provides)
        for k in self.sechost.keys():
            self.deepset(self.provides,self.sechost[k].provides,k)

    def __repr__(self):
        return 'NestingCompositer<Primary('+repr(self.primary)+'),'+repr(self.sechost)+'>'

    def __eatFrame(self,iterat):
        try:
            return iterat.next()
        except StopIteration:
            return self.noFramePlaceholder

    def combine(self):
        iterati=OrderedDict()
        for k in self.sechost.keys():
            iterati[k]=iter(self.sechost[k])
        try:
            for frame in self.primary:
                frame=copy.deepcopy(frame)
                for k in iterati.keys():
                    sf=self.__eatFrame(iterati[k])
                    if not (sf is self.noFramePlaceholder):
                        self.deepset(frame,sf,k)
                yield frame
        finally:
            del iterati


class FrameWidth(dplkit.role.filter.aFilter):
    """ Simple filter to convert a framestream into a stream of timedelta time widths of those frames

    :param framestream: source framestream
    :param fieldname: field of the source framestream's frame to use for the timewidth, if time/width not available
    """
    def __init__(self,framestream,fieldname=None):
        super(FrameWidth,self).__init__(framestream)        
        self.framestream=framestream
        self.fieldname=fieldname
        self.provides=dict(start={'type':datetime},width={'type':timedelta})

    def __repr__(self):
        return 'FrameWidth<'+repr(self.framestream)+'['+self.fieldname+']>'

    def process(self):
        start_time=None
        total_count=0.0
        for frame in self.framestream:
            if not isinstance(frame,dict):
                frame=vars(frame)
            if 'start' in frame and 'width' in frame:
                yield dict(start=frame['start'],width=frame['width'])
            elif self.fieldname not in frame:
                yield None
            else:
                v=frame[self.fieldname]
                if v.size==0:
                    yield None
                    continue
                width=v[-1]-v[0]
                if start_time==None:
                    start_time=v[0]
                total_count+=v.size
                if total_count==1:
                    dt=1.0
                else:
                    dt=(v[-1]-start_time).total_seconds()/(total_count-1)
                yield dict(start=v[0],width=(v[-1]+timedelta(seconds=dt))-v[0])#FIXME better way

class TimeDeGinsu(dplkit.role.filter.aFilter):
    """ Filter converter of a simple framestream to a compound framestream based on time

    :param countstream: FrameWidth filtered stream used to describe how wide in time each resulting compund frame should be
    :param source: source simple framestream
    :param emptyframe: if no frame can satisfy the countstream, this will be returned. default is to ignore such frames. if set to non-None, will not ignore
    :param sendEmpty: if this is set True, any value put to emptyframe will be returned for null frames, even None
    """
    def __init__(self,countstream,source,emptyframe=None,sendEmpty=None):
        super(TimeDeGinsu,self).__init__(source)        
        self.source=source
        self.countstream=countstream
        self.noFramePlaceholder="NOFRAME"
        self.emptyframe=emptyframe
        if sendEmpty==None:
            self.sendEmpty=(emptyframe!=None)
        else:
            self.sendEmpty=sendEmpty

    def __repr__(self):
        return 'TimeDeGinsu<'+repr(self.source)+'>'

    def nextFrame(self,iterat):
        try:
            return iterat.next()
        except StopIteration:
            return self.noFramePlaceholder

    def endTime(self,fr):
        if not isinstance(fr,dict):
            fr=vars(fr)
        if 'start' in fr and 'width' in fr:
            return fr['start']+fr['width']
        if 'times' in fr:
            ret=fr['times'][-1]
            if 'delta_t' in fr:
                ret+=timedelta(seconds=fr['delta_t'][-1])
        raise RuntimeError('No Times')

    def startTime(self,fr):
        if not isinstance(fr,dict):
            fr=vars(fr)
        if 'start' in fr and 'width' in fr:
            return fr['start']
        if 'times' in fr:
            ret=fr['times'][0]
        raise RuntimeError('No Times')

    def updateWidth(self,mainframe,newframe,start=None):
        if start!=None:
            mainframe.start=start
        et=self.endTime(newframe)
        mainframe.width=et-mainframe.start

    def process(self):
        iterat=iter(self.source)
        for framewidth in self.countstream:
            if framewidth is None:
                if self.sendEmpty:
                    yield self.emptyframe
            else:
                ret=self.nextFrame(iterat)
                st=None
                et=None
                try:
                    st=self.startTime(framewidth)
                    et=self.endTime(framewidth)
                except RuntimeError:
                    if self.sendEmpty:
                        yield self.emptyframe
                    continue
                while (not(ret is self.noFramePlaceholder)) and self.endTime(ret)<st:
                    ret=self.nextFrame(iterat)
                if ret is self.noFramePlaceholder:
                    if self.sendEmpty:
                        yield self.emptyframe
                else:
                    ret=copy.deepcopy(ret)
                    if self.endTime(ret)<et:
                      while True:#(not(ret is self.noFramePlaceholder)) and startTime(ret)<st:
                        f=self.nextFrame(iterat)
                        if not f is self.noFramePlaceholder:
                            st=self.startTime(ret)
                            ret.append(f)
                            self.updateWidth(ret,f,start=st)
                        if (f is self.noFramePlaceholder) or self.endTime(ret)>et:
                            break
                    yield ret


class FrameLength(dplkit.role.filter.aFilter):
    """ Simple filter to convert a framestream into a stream of int lengths of those frames

    :param framestream: source framestream
    :param fieldname: field of the source framestream's frame to use for the length
    """
    def __init__(self,framestream,fieldname=None):
        super(FrameLength,self).__init__(framestream)        
        self.framestream=framestream
        self.fieldname=fieldname
        self.provides={'type':int}

    def __repr__(self):
        return 'FrameLength<'+repr(self.framestream)+ (('['+self.fieldname+']') if self.fieldname is not None else '')+'>'

    def process(self):
        for frame in self.framestream:
          if self.fieldname is None:
            if frame is None:
                yield 0
            else:
                try:
                    yield len(frame)
                except TypeError:
                    yield 1
          else:
            if not isinstance(frame,dict):
                try:
                    frame=vars(frame)
                except TypeError:
                    yield 1
                    continue
            if self.fieldname not in frame:
                yield 0
            else:
                try:
                    yield len(frame[self.fieldname])
                except TypeError:
                    yield 1

class CountDeGinsu(dplkit.role.filter.aFilter):
    """ Filter converter of a simple framestream to a compound framestream based on count

    :param countstream: FrameLength filtered stream used to describe length of each resulting compund frame
    :param source: source simple framestream
    :param emptyframe: frame to return if countstream yields a 0-length value. default None will ignore such frames
    :param sendEmpty: if this is set True, any value put to emptyframe will be returned for 0-length frames, even None
    """
    def __init__(self,countstream,source,emptyframe=None,sendEmpty=None):
        super(CountDeGinsu,self).__init__(source)        
        self.source=source
        self.countstream=countstream
        self.noFramePlaceholder="NOFRAME"
        self.emptyframe=emptyframe
        if sendEmpty==None:
            self.sendEmpty=(emptyframe!=None)
        else:
            self.sendEmpty=sendEmpty

    def __repr__(self):
        return 'CountDeGinsu<'+repr(self.source)+'>'

    def nextFrame(self,iterat):
        try:
            if iterat['iterat'] is not None:
                return iterat['iterat'].next()
        except StopIteration:
            iterat['iterat']=None
        return self.noFramePlaceholder

    def process(self):
        iterat=dict(iterat=iter(self.source))
        for frcount in self.countstream:
            if frcount==0:
                if self.sendEmpty:
                    yield self.emptyframe
            else:
                ret=self.nextFrame(iterat)
                if ret is self.noFramePlaceholder:
                    if self.sendEmpty and not self.fixedWidth:
                        yield self.emptyframe
                else:
                    ret=copy.deepcopy(ret)
                    for x in range(1,frcount):
                        f=self.nextFrame(iterat)
                        if not f is self.noFramePlaceholder:
                            ret.append(f)
                        elif self.sendEmpty:
                            if hasattr(ret,'extend'):
                                ret.extend(1)
                            else:
                                ret.append(self.emptyframe)
                    yield ret

@dplkit.role.decorator.is_stateful
class SubstructMerger(dplkit.role.blender.aBlender):
    """ merge multiple frame streams into a new frame DEPRECATED

    :param primaryname: name in streams to set the rhythm of the resulting framestream
    :param streams: dictionary of framestreams, keyed to name of subframe within the returned frame
    :param asClass: class of the topmost frame to return. will use a dictionary as default
    :param classParms: parmeters to class initializer
    :param classDictParameter: name of parameter to pass dictionary to start initializer. if not specified, will loop using setattr
    :param countMatch: names of streams that need to match the counted length of the primary frame stream
    """
    def __init__(self,primaryname,streams,asClass=None,classParms={},classDictParameter=None,countMatch=[]):#if asClass is none, will return a dictionary, otherwise will make a new class, and use setattr. if classdictparameter is
        super(SubstructMerger,self).__init__([streams[primaryname]] + [streams[x] for x in streams.keys()])        
        self.primaryname=primaryname
        self.sechost=streams
        self.noFramePlaceholder='NOFRAME'
        self.prim=streams[primaryname]
        del streams[primaryname]
        self.provides={primaryname:self.prim.provides}
        for k in self.sechost.keys():
            self.provides[k]=self.sechost[k].provides
        self.retclass=asClass
        self.classparms=classParms
        self.classdictparm=classDictParameter
        self.countMatch=set(countMatch)

    def __repr__(self):
        return 'SubstructMerger<Primary('+self.primaryname+','+repr(self.prim)+'),'+repr(self.sechost)+'>'

    def __eatFrameByCount(self,n,track,count):
        needcount=count
        if count<=0:
            return self.noFramePlaceholder,None
        ret=track['store']
        track['store']=[]
        needcount-=len(ret)
        try:
            while track['iter']!=None:
                f=track['iter'].next()
                needcount-=1
                if needcount<0:
                    track['store'].append(f)
                    needcount=0
                    break
                ret.append(f)
                if False and needcount==0:
                    break
        except StopIteration:
            track['iter']=None
        if len(ret)>0:
            tw=dict(start=self.__getFirstTime(ret[0]))
            if tw['start']!=None:
                tw['width']=self.__getLastTime(ret[len(ret)-1])-tw['start']
            else:
                tw=None
        else:
            return self.noFramePlaceholder,None
        if needcount>0:
            VERBOSE ('eat ended short')
        VERBOSE ('Eat by count got',len(ret),'needs',needcount,'more')
        for x in range(needcount):
            ret.append(self.noFramePlaceholder)
        return ret,tw

    def __eatFrameByTime(self,n,track,lastTime):
        ret=track['store']
        track['store']=[]
        tw=None
        try:
            if len(ret)>0:
                l=len(ret)
                if self.__getLastTime(ret[l-1])>=lastTime:
                    tw={}
                    tw['start']=self.__getFirstTime(ret[0])
                    tw['width']=self.__getLastTime(ret[l-1])-tw['start']
                    if tw['start']==None:
                        return ret[-1],None
                    return ret,tw
            while track['iter']!=None:
                f=track['iter'].next()
                if self.__getFirstTime(f)==None and len(ret)==0:
                    return f,None
                if self.__getFirstTime(f)>=lastTime:
                    track['store'].append(f)
                    VERBOSE ('out of range',n,lastTime,f)
                    if len(ret)==0:
                        return self.noFramePlaceholder,None
                    return ret,tw
                if tw == None:
                    tw={}
                    tw['start']=self.__getFirstTime(f)
                tw['width']=self.__getLastTime(f)-tw['start']
                ret.append(f)
                #if self.__getLastTime(f)>=lastTime:
                #   print 'end of range',n,f
                #   return ret,tw
        except StopIteration:
            del track['iter']
            track['iter']=None
        VERBOSE ('exhuasted generator',n)
        if tw==None:
            return self.noFramePlaceholder,None
        return ret,tw

    def __hasTimes(self,fr):
        return hasattr(fr,'_timevarname') and self.__getCount(fr)>0

    def __getTime(self,fr,idx,asArray=False):
        if asArray:
            idx=numpy.array([idx])
        return getattr(fr,fr._timevarname)[idx]


    def __getFirstTime(self,fr):
        if hasattr(fr,'start') and hasattr(fr,'width'):
            return fr.start
        if self.__hasTimes(fr):
            return self.__getTime(fr,0)
        if isinstance(fr,dict) and 'start' in fr and 'width' in fr:
            return fr['start']
        #raise 'NoTime'
        return None

    def __getLastTime(self,fr,dt=None):
        if hasattr(fr,'start') and hasattr(fr,'width'):
            return fr.start+fr.width
        if self.__hasTimes(fr):
            t=self.__getTime(fr,-1)
            if dt!=None:
                return t+dt
            if self.__getCount(fr)>1:
                return t+timedelta(seconds=(self.__getTime(fr,-1)-self.__getTime(fr,0)).total_seconds()/(self.__getCount(fr)-1))
            LOG.warning ("Warning: Can't find width of frame. faking 0 width. might cause slicing issues down the road")
            return t#+timedelta(seconds=1)
        if isinstance(fr,dict) and 'start' in fr and 'width' in fr:
            return fr['start']+fr['width']
        #raise 'NoTime'
        return None

    def __getCount(self,fr):
        if hasattr(fr,'_timevarname'):
            return getattr(fr,fr._timevarname).size
        return 1

    def combine(self):
        tracks=OrderedDict()
        for k in self.sechost.keys():
            tracks[k]={}
            tracks[k]['iter']=iter(self.sechost[k])
            tracks[k]['store']=[]
        if True:#try:
          for pf in self.prim:
            ret=OrderedDict()
            ret[self.primaryname]=pf
            lastTime=self.__getLastTime(pf)
            frameCount=self.__getCount(pf)
            for k in tracks.keys():
                if k in self.countMatch:
                    frames,timerange=self.__eatFrameByCount(k,tracks[k],frameCount)
                else:
                    frames,timerange=self.__eatFrameByTime(k,tracks[k],lastTime)
                if timerange!=None and len(frames)>0 and hasattr(frames[0],'append'):
                    fms=frames
                    frames=copy.deepcopy(fms[0])
                    for x in range(1,len(fms)):
                        if fms[x] is self.noFramePlaceholder:
                            if k in self.countMatch:
                                VERBOSE( 'extending',k,'FIXME this should have been done by the resampler')
                                frames.append(hau.Time_Z_Group(like=frames,times=self.__getTime(pf,self.__getCount(frames),asArray=True)))
                            else:
                                pass#by time, or no extend
                        else:
                            frames.append(fms[x])
                    setattr(frames,'start',timerange['start'])
                    setattr(frames,'width',timerange['width'])
                if not frames is self.noFramePlaceholder:
                    ret[k]=frames
            if self.retclass==None:
                retv=ret
            elif self.classdictparm:
                d=self.classparms.copy()
                d[self.classdictparm]=ret
                retv=self.retclass(**d)
            else:
                retv=self.retclass(**self.classparms)
                for k,i in ret.items():
                    setattr(retv,k,i)
                #print 'merged is',retv
            yield retv
        #finally:
        #  for k,tr in tracks.items():
        #    if tr['iter']!=None:
        #        del tr['iter']
        #        tr['iter']=None

@dplkit.role.decorator.exposes_attrs_of_field('framestream')
class Retyper(dplkit.role.narrator.aNarrator):
    """ Filter to re-type a dictionary framestream to somethign else. used with interpolation
    """
    def __init__(self,framestream,frametype=None,frameinit={}):
        self.framestream=framestream
        #self.provides=framestream.provides
        self.frametype=frametype
        self.frameinit=frameinit

    def _newframe(self):
        return self.frametype(**self.frameinit)

    def __repr__(self):
        return 'Retyper<'+repr(self.framestream)+' to '+repr(type(self._newframe()))

    def read(self):
        for f in self.framestream:
            #print 'was:',f
            #print 'contains:',f.keys()
            if self.provides!=None:
                _f=f
                if not isinstance(_f,dict):
                    _f=vars(f)
                for k,v in _f.items():
                    if k in self.provides and 'type' in self.provides[k] and type(v)!=self.provides[k]['type']:
                        #print k,self.provides[k],type(v)
                        try:
                            nv=copy.deepcopy(self.provides[k]['type'](v))
                            _f[k]=nv#setattr(f,k,nv)
                        except ValueError:
                            print 'Value error occurred on field '+k+'. ignored'
            if self.frametype!=None:
                rf=self._newframe()
                if type(rf)!=type(f):
                    _rf=rf
                    if not isinstance(_rf,dict):
                        _rf=vars(_rf)
                    for k,v in _f.items():
                        _rf[k]=v
                    #if self.frametype==hau.Time_Z_Group and (not hasattr(rf,rf._timevarname) or getattr(rf,rf._timevarname).shape[0]==0) and hasattr(rf,'start'):
                    #    setattr(rf,rf._timevarname,hau.T_Array([rf.start]))
                    #    #print 'set time'
                    f=rf
            #print 'is now:',f
            #print 'which now contains',vars(f).keys()
            yield f

import lg_base.core.IterableQueue as iq
import threading_filter as tf
import cPickle as pickle

class LocalOrRemoteQueue(object):
    """ Queue type that uses a threaded Queue for same-process communication and multiprocessing Queue for inter-process communication
    BROKEN

    """
    def __init__(self,*args,**kwargs):
        self.local=Queue.Queue(*args,**kwargs)
        self.remote=multiprocessing.JoinableQueue(*args,**kwargs)
        self.switcher=threading.Event()#not set means use remote
        self.putlock=multiprocessing.Condition()
        self.switchset=multiprocessing.Event()
        self.putpid=None
        self.nullreturn='notaframe'

    def __repr__(self):
        if not self.switchset.is_set():
            return 'LocalOrRemoteQueue('+repr(self.local)+','+repr(self.remote)+')'
        return 'LocalOrRemoteQueue('+repr(self._realqueue())+')'

    def __del__(self):
        if self.remote!=None:
            #self.remote.join()
            #self.remote.close()
            if hasattr(self.remote,'cancel_join_thread'):
                self.remote.cancel_join_thread()

    def put(self,*args,**kwargs):
        with self.putlock:
            if not self.switchset.is_set():
                self.putpid=os.getpid()
            self._realqueue().put(*args,**kwargs)            
            self.putlock.notify_all()

    def put_nowait(self,*args,**kwargs):
        with self.putlock:
            if not self.switchset.is_set():
                self.putpid=os.getpid()
            self._realqueue().put_nowait(*args,**kwargs)            
            self.putlock.notify_all()

    def get(self,*args,**kwargs):
        if not self.switchset.is_set():
            try:
                with self.putlock:
                    ret=self.nullreturn
                    while ret is self.nullreturn:
                        try:
                            ret=self._realqueue().get_nowait(*args,**kwargs)            
                        except Queue.Empty:
                            self.putlock.wait()
                    ret=self.tryswitch(firstPut=ret)
                    self.switchset.set()
                    if not ret is self.nullreturn:
                        return ret   
            except (EOFError,IOError):
                if not self.switcher.is_set():
                    raise
                assert(self.switchset.is_set())
                pass #got switched. continue below
        return self._realqueue().get(*args,**kwargs)            

    def get_nowait(self,*args,**kwargs):
        if not self.switchset.is_set():
            try:
                ret=self._realqueue().get_nowait(*args,**kwargs)            
                ret=self.tryswitch(firstPut=ret)
                self.switchset.set()
                if not ret is self.nullreturn:
                    return ret   
            except (EOFError,IOError):
                if not self.switcher.is_set():
                    raise
                assert(self.switchset.is_set())
                pass #got switched. continue below
        return self._realqueue().get_nowait(*args,**kwargs)            

    def join(self,*args,**kwargs):
        q=self.remote
        if q!=None:
            q.join(*args,**kwargs)
        q=self.local
        if q!=None:
            q.join(*args,**kwargs)
        #self._realqueue().join(*args,**kwargs)            

    def task_done(self,*args,**kwargs):
        assert(self.switchset.is_set())
        self._realqueue().task_done(*args,**kwargs)            

    def close(self,*args,**kwargs):
        q=self._realqueue()
        if hasattr(q,'close'):
            q.close(*args,**kwargs)

    def join_thread(self,*args,**kwargs):
        q=self._realqueue()
        if hasattr(q,'join_thread'):
            q.join_thread(*args,**kwargs)

    def cancel_join_thread(self,*args,**kwargs):
        q=self._realqueue()
        if hasattr(q,'cancel_join_thread'):
            q.cancel_join_thread(*args,**kwargs)

    def tryswitch(self,firstPut):
        if self.putpid!=None and self.putpid==os.getpid():
            return self.switch(firstPut=firstPut,doFirstPut=True)#is on same pid as put. use threaded version
        return firstPut

    def switch(self,doFirstPut=False,firstPut=None):#switching to the local queue
        if self.switcher.is_set():
            return firstPut
        rem=None
        with self.putlock:
            try:
                if self.switcher.is_set():
                    return firstPut
                self.switcher.set()
                self.switchset.set()
                rem=self.remote
                #self.remote=None
                if doFirstPut:
                    self.local.put(firstPut)
                    rem.task_done()
                while True:
                    x=rem.get_nowait()
                    self.local.put(x)
                    rem.task_done()
            except Queue.Empty:
                self.remote=None
                rem.join()
                rem.close()
                rem.join_thread()
                del rem
                #self.putlock.release()
                VERBOSE ('Remote deleted')
        return self.nullreturn

    def _realqueue(self):
        if self.remote==None:
            return self.local
        if self.local==None or self.switchset.is_set():
            if self.local!=None:
                if self.putlock.acquire(False):
                    if self.remote!=None:
                        self.local=None
                    self.putlock.release()
                return self._realqueue()
            return self.remote
        with self.putlock:
            if self.remote==None or self.switcher.is_set():
                return self.local
            else:
                return self.remote         

@dplkit.role.decorator.is_stateful
@dplkit.role.decorator.exposes_attrs_of_field('host')
class QueueAsNarrator(dplkit.role.narrator.aNarrator):
    """Convert a joinable Queue object to a DPL narrator. Used in concert with a IterableToQueue object

    :param queue: queue used as a datasource
    :param attribute_source: source of exposed attributes for framestream
    :param name: descriptive name of this object
    :param provides: provides dictionary
    :param preStep: optional callable to be called before every call to read from the queue. can be used to trigger data to be queued
    :param postStep: optional callable to be called after every call to read from the queue, but before the data is yielded
    :param atStart: called just before iteration begins
    :param atExit: called at the end of the queued data
    :param atAbort: called if this narrator's generation is terminated prematurely
    """
    def __init__(self,queue,attribute_source=None,name=None,provides=None,preStep=None,postStep=None,atStart=None,atExit=None,atAbort=None):
        self.frames=queue#multiprocessing.JoinableQueue()
        self.host=attribute_source
        self.provides=provides
        self.name=name
        self.atStart=atStart
        self.atExit=atExit
        self.atAbort=atAbort
        self.preStep=preStep
        self.postStep=postStep
        self.wrapper=None
        self.itwrapper=None

    def __repr__(self):
        return 'QueueAsNarrator<'+repr(self.host)+'['+self.name+']<'+repr(self.frames)+'>>'

    def read(self):
        if self.frames is None:
            return
        if self.wrapper is None:
            self.wrapper=iq.QueueAsIterable(self.frames,atAbort=self.atAbort,atExit=self.atExit,atStart=self.atStart)
            self.itwrapper=iter(self.wrapper)
        it=self.itwrapper
        VERBOSE ('QueueAsIterable in substruct ' + self.name+' for queue '+repr(self.frames)+' is '+repr(it))
        try:
            while True:
                if self.preStep is not None:
                    self.preStep()
                f=it.next()
                if self.postStep is not None:
                    self.postStep()
                yield f
        except StopIteration:
            if self.postStep is not None:
                self.postStep()
            return
        except GeneratorExit as e:
            del it
            self.itwrapper=None
            self.wrapper=None

#this is a single-use iterable object.
@dplkit.role.decorator.exposes_attrs_of_field('framestream')
class Brancher(object):#(dplkit.role.narrator.aNarrator):
    """ Complex Framestream to Simple Framestream filter source.

    :param framestream: source framestream
    """ 
    multiprocessable=True
    def __init__(self,framestream):
        self.framestream=framestream
        self.objs=[]
        #self.i=None#iter(self.framestream)
        self.thread=None
        self.savedFrames=[]
        self.q=Queue.Queue(1)
        if self.multiprocessable:
            self.nextlock=multiprocessing.Condition()
            self.doSave=multiprocessing.Event()
            self.loadEvent=multiprocessing.Event()
            self.terminated=multiprocessing.Event()
            self.finished=multiprocessing.Event()
            self.threadstartlock=multiprocessing.Lock()
            self.threadstarted=multiprocessing.Event()
            self.queueclass=multiprocessing.JoinableQueue#LocalOrRemoteQueue
        else:
            self.nextlock=threading.Condition()
            self.doSave=threading.Event()
            self.loadEvent=threading.Event()
            self.terminated=threading.Event()
            self.finished=threading.Event()
            self.threadstartlock=threading.Lock()
            self.threadstarted=threading.Event()
            self.queueclass=Queue.Queue
        self.thread=None# iq.IterableToQueueThread(framestream,self.q,inSync=True)#self.terminator)#,#threading.Thread(target=brancher_nextloop,args=[weakref.ref(self)]) 
        self.pid=None
        #self.startThread()
        self.provides=framestream.provides
        if self.provides==None:
            raise RuntimeError('No provides provided')

    def startThread(self):
        with self.threadstartlock:
            if not self.threadstarted.is_set():
                self.threadstarted.set()
                self.thread=threading.Thread(target=self.trytask,args=(weakref.ref(self),))
                self.thread.start()
                self.pid=os.getpid()
            elif self.terminated.is_set():
                return False
            return True

    @staticmethod
    def trytask(*args,**kwargs):
        try:
            SubstructBrancher.task(*args,**kwargs)
        except:
            import traceback
            traceback.print_exc()
            raise

    @staticmethod
    def task(selfref):
        s=selfref()
        if s==None:
            return
        desc=('Thread %i - ' % (threading.currentThread().ident)) +repr(s.framestream)
        loadEvent=s.loadEvent
        terminated=s.terminated
        finished=s.finished
        objs=s.objs
        doSave=s.doSave
        savedFrames=s.savedFrames
        framestream=iter(s.framestream)
        del s
        finished.clear()
        terminated.clear()
        try:
            while True:
                if terminated.is_set() or finished.is_set():
                    break
                VERBOSE('WAITING '+desc)
                loadEvent.wait()
                VERBOSE('DONE WAITING '+desc)
                if terminated.is_set() or finished.is_set():
                    break
                VERBOSE('Getting '+desc)
                try:
                    f = framestream.next()
                except StopIteration:
                    framestream.close()
                    finished.set()
                    break
                VERBOSE('GOT '+desc)
                loadEvent.clear()
                if terminated.is_set():
                    break
                VERBOSE('PUTTING '+desc)
                if doSave.is_set():
                    savedFrames.append(f)
                elif len(savedFrames)>0:
                    while len(savedFrames)>0:
                        savedFrames.pop()
                    #savedFrames=[] #this is probably not right. it will run an iteration after the ifrst, and lose any new narration. sigqueue might fix it
                for c in objs:
                    c.put(f)
                VERBOSE('PUT '+desc)
        except Exception as e:
            framestream.close()
            VERBOSE ('Exception on substruct')
            print 'substruct raised'
            traceback.print_exc()
            #terminated.set()
            if doSave.is_set():
                VERBOSE ('Substruct is saving')
                savedFrames.append(iq.RaisedExceptionWrapper(e))
            elif len(savedFrames)>0:
                    while len(savedFrames)>0:
                        savedFrames.pop()
            else:
                print 'Brancher',selfref(),'happy with',len(objs),'streams'
            for c in objs:
                c.put(iq.RaisedExceptionWrapper(e))
            for c in objs:
                if hasattr(c,'cancel_join_thread'):
                    c.cancel_join_thread()
                if hasattr(c,'close'):
                    c.close()
            while len(objs)>0:
                objs.pop()
            return
        VERBOSE ('Substruct exiting '+desc)
        if terminated.is_set():
            for c in objs:
                if hasattr(c,'cancel_join_thread'):
                    c.cancel_join_thread()
                if hasattr(c,'close'):
                    c.close()
            while len(objs)>0:
                objs.pop()
        else:
            if doSave.is_set():
                VERBOSE ('Substruct is saving')
                savedFrames.append(iq.myStopIteration)
            else:
                savedFrames=[] #this is probably not right. it will run an iteration after the ifrst, and lose any new narration. sigqueue might fix it
            consumed=[]
            for c in objs:
                VERBOSE('putting stop '+desc+' '+repr(c))
                c.put(iq.myStopIteration)
                VERBOSE('done putting stop '+desc+' '+repr(c))
                consumed.append(c)
            for c in consumed:
                if c in objs:
                    objs.remove(c)
            for c in consumed:
                VERBOSE('JOINING '+desc+' '+repr(c))
                c.join()
                VERBOSE('Done JOINING '+desc+' '+repr(c))
        finished.set()
        VERBOSE ('substruct exited '+desc)

    def __del__(self):
        if self.pid is not None and self.pid==os.getpid() and self.thread is not None and self.thread.is_alive():
            self.terminate()
            for c in self.objs:
                try:
                    while True:
                        c.get_nowait()
                        c.task_done()
                except Queue.Empty:
                    pass
                if hasattr(c,'cancel_join_thread'):
                    c.cancel_join_thread()
            while len(self.objs)>0:
                self.objs.pop()
            if True:
                VERBOSE('JOining in '+ ('Thread %i - ' % self.thread.ident) +repr(self))
                self.thread.join() #FIXME 20140109 JPG this causes a hang if there was an exception in the running
                VERBOSE('Done joining in '+ ('Thread %i - ' % self.thread.ident)+repr(self))
       

    def queueStep(self):
        self.loadEvent.set()

    def finish(self):
        VERBOSE('Finishing substruct '+repr(threading.currentThread().ident))
        if self.finished.is_set():
            VERBOSE("ALREADY FINISHED")
        self.finished.set()
        self.queueStep()
        VERBOSE("FINISH COMPLETE")

    def terminate(self):
        if self.finished.is_set():
            VERBOSE ('No terminate. already finished')
            return
        VERBOSE('terminating substruct '+("finished already" if self.finished.is_set() else ""))
        self.terminated.set()
        self.queueStep()

    def __repr__(self):
        return (super(Brancher,self).__repr__())+'<'+repr(self.framestream)+'>'

    @staticmethod
    def callAndGet(q,callb,removefunc=None):
        class pullT(object):
            def __init__(self,q,c,r):
                self.q=q
                self.c=c
                self.r=r

            def __call__(self):
                if self.r is not None:
                    self.r(self.q)
                self.c()
                try:
                    while True:
                        self.q.get_nowait()
                        self.q.task_done()
                except Queue.Empty:
                    pass
        return pullT(q,callb,removefunc)

    def removequeue(self,q):
        if q in self.objs:
            self.objs.remove(q)

    def narrateSubstruct(self,realsubstr=None):
        """ retrieve a narrator for a subframe of the source. each returned value can only be used in one isinstance

        :param substrname: subframe to retrieve, or None for complete frame instead (used with branch forking)
        :returns: narrator object of the requested subframe stream
        """
        #print self.provides
        with self.nextlock:
            while self.doSave.is_set():
                self.nextlock.wait()
            self.doSave.set()
            q=iq.QueueCount(self.queueclass())
            self.objs.append(q)
            for f in self.savedFrames:
                q.put(f)
        args=[q]
        kwargs=dict(attribute_source=self,name="WHOLE",preStep=self.queueStep,
            atStart=self.startThread,atExit=self.callAndGet(q,self.finish),atAbort=self.callAndGet(q,self.terminate,self.removequeue))
        try:
            ret=QueueAsNarrator(provides=self.provides,*args,**kwargs)
            if realsubstr is not None:
                ret=NestedFrame(ret,realsubstr)
        finally:
            with self.nextlock:
                self.doSave.clear()
                self.nextlock.notify_all()
        #ret=tf.KeepUpstreamOnThisProcess(ret,1)
        return ret

    def __call__(self):
        return self.narrateSubstruct()

SubstructBrancher=Brancher

if __name__ == '__main__':
    main()
