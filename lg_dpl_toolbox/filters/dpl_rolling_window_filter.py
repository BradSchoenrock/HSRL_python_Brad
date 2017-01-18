import dplkit.role.filter
import dplkit.role.decorator
import copy
from datetime import datetime,timedelta
import lg_base.core.array_utils as hau
import numpy as np
import os
import json
import traceback
import Queue

from collections import namedtuple
#WindowedFilterDescription=namedtuple('WindowedFilterDescription','filter width edgemode varlist time_width modrequirement cargs kwcargs')
class WindowedFilterDescription(namedtuple('WindowedFilterDescription','filter width edgemode varlist time_width modrequirement windowinfo_perframe middleframeparametername cargs kwcargs')):
    """ Configuration object for dpl_rolling_window_filter object

     :param filter: callable object, takes in list of frames in the window, or called with a 2D array for each array in varlist
     :param width: record count. if time_width is also specified, this is the minimum record count (even if time represented exceeds time_width) or None
     :param edgemode: how to handle edges. this can be 'short' or 'fullduplicate'. with short, the start and end will allow for a shorter window without duplicating data. with fullduplicate, the window will be met, duplicating the first or last records as neccessary to increase the number of data records
     :param varlist: a list of variable names to accumulate and call the filter for each separately. or set to None for a list of frames to be passed instead
     :param time_width: timedelta object describing the window width as a factor of time. set to None if not used
     :param modrequirement: 2-element tuple of the required window width. typically set to None, or (2,1) to mean len(window)%2 must equal 1.
     :param middleframeparametername: name of parameter to pass the central frame to (default None, not passed)
     :param cargs: additional unnamed args to the filter
     :param kwcargs: additional named args to the filter
    """
    def __new__(cls,filter,width,edgemode,varlist=None,time_width=None,modrequirement=None,windowinfo_perframe="windowdesc",middleframeparametername=None,cargs=None,kwcargs=None):
        return super(WindowedFilterDescription,cls).__new__(cls,filter=filter,width=width,edgemode=edgemode,varlist=varlist,time_width=time_width,\
                 modrequirement=modrequirement,windowinfo_perframe=windowinfo_perframe,middleframeparametername=middleframeparametername,cargs=cargs,kwcargs=kwcargs)

def deepattribute(host,fieldname):
    try:
        if isinstance(fieldname,basestring):
            return deepattribute(host,fieldname.split('.'))
        if len(fieldname)==0:
            return host
        if not isinstance(host,dict):
            return deepattribute(vars(host),fieldname)
        #if fieldname[0] not in host:
        #    raise RuntimeError(fieldname,' field not found')
        return deepattribute(host[fieldname[0]],fieldname[1:])
    except KeyError as e:
        raise KeyError(fieldname)
def setdeepattribute(host,fieldname,val):
    try:
        if isinstance(fieldname,basestring):
            return setdeepattribute(host,fieldname.split('.'),val)
        if not isinstance(host,dict):
            return setdeepattribute(vars(host),fieldname,val)
        if len(fieldname)==1:
            host[fieldname[0]]=val
            return
        #if fieldname[0] not in host:
        #    raise RuntimeError(fieldname,' field not found')
        return setdeepattribute(host[fieldname[0]],fieldname[1:],val)
    except KeyError as e:
        raise KeyError(fieldname)

class dpl_rolling_window_filter(dplkit.role.filter.aFilter):
    """ Generic configurable Dpl filter object for gathering a subslice of data from a simple stream

    :param source: input iterable framestream
    :param parms: parameter instance of WindowedFilterDescription

    What this does:
    depending upon parameters in the WindowedFilterDescription, roll thru time, calling a filter for each step, including additional surrounding data as configured.
    this makes it simple to create a time-based filter that operates on a simple framestream, leaving the multiple record tracking to this module, and the filter to do its work with what it receives.
    """

    def __new__(klass,source,parms):
        if parms.time_width!=None:
            return dpl_rollingtime_window_filter(source,parms) #this does time. dpl_rolling_window_filter only does count. simpler
        return

    def __init__(self,source,parms):
        super(dpl_rolling_window_filter,self).__init__(source)
        #inherits provides from source
        self.source=source
        assert(isinstance(parms,WindowedFilterDescription))
        self.width=parms.width
        self.callable=parms.filter
        self.edgemode=parms.edgemode or ('fullduplicate' if parms.time_width==None else 'short') # can be 'fullduplicate' or 'short'
        self.varlist=parms.varlist
        self.cargs=parms.cargs or []
        self.kwcargs=parms.kwcargs or {}
        self.windowinfo_perframe=parms.windowinfo_perframe
        self.middleframeparametername=parms.middleframeparametername
        assert(self.edgemode in ('fullduplicate','short'))
        assert(self.width>0)
        self.halfwidth=int(self.width/2)
        print ('Filter window is ',self.width)
        assert(self.callable!=None)

    def docall(self,*args,**kwargs):
        try:
           return self.realdocall(*args,**kwargs)
        except Exception as e:
           import traceback
           traceback.print_stack()
           traceback.print_exc()
           raise

    def realdocall(self,myslice,templaterec):
        """ call the filter with the given slices as input

        :param myslice: list of frames in the window
        :param templaterec: record in the window to base the return off of

        if this filter was given a varlist, this will call the filter for each individual variable, with a 2D array as the content. the return of each call will replace the content of each var in templaterec and returned
        if this filter was not given a varlist, the list of frames passed to this function will instead be passed to a single call of the filter. will return what the filter returns

        :returns: simple frame result
        """
        if self.varlist==None:
            kwcargs=self.kwcargs
            if self.middleframeparametername is not None:
                kwcargs=kwcargs.copy()
                kwcargs[self.middleframeparametername]=copy.deepcopy(templaterec)
            ret=self.callable(myslice,*self.cargs,**kwcargs)
        else:
            ret=copy.deepcopy(templaterec)
            for i,v in enumerate(self.varlist):
                try:
                    av=deepattribute(myslice[0],v)
                except KeyError:
                    print 'Warning: configured key',v,'not found in frame'
                    raise
                    continue
                av=copy.deepcopy(av)
                for f in myslice[1:]:
                    av.append(deepattribute(f,v))
                setdeepattribute(ret,v,self.callable(av,*self.cargs,**self.kwcargs))
        return ret

    def process(self):
        myslice=[]
        slicequeue=Queue.Queue()
        extraloop=0
        haveWindow=False

        for f in self.source:
            myslice.append(f)
            myslen=len(myslice)
            slicequeue.put(f)
            if haveWindow==False and myslen<self.halfwidth:
                extraloop+=1
                continue
            elif haveWindow==False:
                if self.edgemode=='fullduplicate':
                    while len(myslice)<self.width:
                        myslice.insert(0,myslice[0])
                haveWindow=True
            else:
                del myslice[0]
            requiredslice=slicequeue.get_nowait()
            #print ('filter window', len(myslice),self.underHalf(myslice),self.half(myslice),self.overFull(myslice))
            yield self.docall(myslice,requiredslice)
        if extraloop>0:
            try:
                for x in range(extraloop):
                    if self.edgemode=='fullduplicate':
                        myslice.append(myslice[-1])
                    requiredslice=slicequeue.get_nowait()
                    del myslice[0]
                    #print ('filter window wrapping up',len(myslice),self.underHalf(myslice),self.half(myslice),self.overFull(myslice))
                    yield self.docall(myslice,requiredslice)
            except Queue.Empty:
                raise RuntimeError('Premature rolling window termination!')
        try:
            slicequeue.get_nowait()
            raise RuntimeError('Incomplete rolling window? BUG')
        except Queue.Empty:
            pass#good


class dpl_rollingtime_window_filter(dplkit.role.filter.aFilter):
    """ Generic configurable Dpl filter object for gathering a subslice of data from a simple stream

    :param source: input iterable framestream
    :param parms: parameter instance of WindowedFilterDescription

    What this does:
    depending upon parameters in the WindowedFilterDescription, roll thru time, calling a filter for each step, including additional surrounding data as configured.
    this makes it simple to create a time-based filter that operates on a simple framestream, leaving the multiple record tracking to this module, and the filter to do its work with what it receives.
    """
    def __init__(self,source,parms):
        super(dpl_rollingtime_window_filter,self).__init__(source)
        #inherits provides from source
        self.source=source
        assert(isinstance(parms,WindowedFilterDescription))
        self.width=parms.width
        self.callable=parms.filter
        self.modrequirement=parms.modrequirement
        self.edgemode=parms.edgemode or 'fullduplicate' # can be 'fullduplicate' or 'short'
        self.varlist=parms.varlist
        self.time_width=parms.time_width
        self.cargs=parms.cargs or []
        self.windowinfo_perframe=parms.windowinfo_perframe
        self.kwcargs=parms.kwcargs or {}
        self.middleframeparametername=parms.middleframeparametername
        assert(self.edgemode in ('fullduplicate','short'))
        assert(self.time_width.total_seconds()>0)
        self.halfwidth=timedelta(seconds=self.time_width.total_seconds()/2.0)
        print ('Filter window is ',self.time_width,' minimum of ',self.width)
        assert(self.callable!=None)

    def docall(self,*args,**kwargs):
        try:
           return self.realdocall(*args,**kwargs)
        except Exception as e:
           import traceback
           traceback.print_stack()
           traceback.print_exc()
           raise

    def realdocall(self,myslice,templaterec):
        """ call the filter with the given slices as input

        :param myslice: list of frames in the window
        :param templaterec: record in the window to base the return off of

        if this filter was given a varlist, this will call the filter for each individual variable, with a 2D array as the content. the return of each call will replace the content of each var in templaterec and returned
        if this filter was not given a varlist, the list of frames passed to this function will instead be passed to a single call of the filter. will return what the filter returns

        :returns: simple frame result
        """
        if not self.inMiddle(myslice,templaterec):
            print 'FAILED TO ALIGN WINDOW: len ',len(myslice),' index ',self.indexOf(myslice,templaterec)
        b=hau.T_Array(np.ones((1,6)))
        c=len(myslice)
        x=self.indexOf(myslice,templaterec)
        b[0,0]=c
        b[0,1]=x
        b[0,2]=self.delta_t(myslice).total_seconds()
        b[0,3]=self.delta_t(myslice[:x] if x>0 else []).total_seconds()
        b[0,4]=self.delta_t(myslice[x:(x+1)]).total_seconds()
        b[0,5]=self.delta_t(myslice[(x+1):] if x<(c-1) else []).total_seconds()
       # print 'rolling window stats',b
        if self.varlist==None:
            kwcargs=self.kwcargs
            if self.middleframeparametername is not None:
                kwcargs=kwcargs.copy()
                kwcargs[self.middleframeparametername]=copy.deepcopy(templaterec)
            ret=self.callable(myslice,*self.cargs,**kwcargs)
        else:
            ret=copy.deepcopy(templaterec)
            for i,v in enumerate(self.varlist):
                try:
                    av=deepattribute(myslice[0],v)
                except KeyError:
                    print 'Warning: configured key',v,'not found in frame'
                    raise
                    continue
                av=copy.deepcopy(av)
                for f in myslice[1:]:
                    av.append(deepattribute(f,v))
                setdeepattribute(ret,v,self.callable(av,*self.cargs,**self.kwcargs))
        _ret=ret
        if not isinstance(ret,dict):
            _ret=vars(ret)
        if b is not None and self.windowinfo_perframe is not None:
            if self.windowinfo_perframe not in _ret:
                _ret[self.windowinfo_perframe]=b
            else:
                _ret[self.windowinfo_perframe]=np.concatenate((_ret[self.windowinfo_perframe],b),1)
        return ret

    def delta_t(self,frames):
        """ get the time span of the list of frames

        :param frames: list of frames
        :returns: timedelta object
        """
        delt=0.0
        for x in frames:
            if not isinstance(x,dict):
                x=vars(x)
            if x['width'] is not None and x['width'].total_seconds()>0:
                delt+=x['width'].total_seconds()
        #print 'framewidth is ',delt
        return timedelta(seconds=delt)

    def underHalf(self,frames):
        """ figure out if the framelist is under half complete for the window

        :param frames: list of frames
        :returns: true if the count isn't half satisfied, or the timewidth isn't half satisfied
        """
        if len(frames)<1:
            return True
        if self.width!=None and len(frames)<=int(self.width/2):
            return True
        return self.delta_t(frames)<self.halfwidth

    def half(self,frames):
        """ figure out if the framelist is exactly half complete for the window as it can be

        :param frames: list of frames
        :returns: true if dropping the first frame would make this underHalf, but as given is not underHalf
        """
        if len(frames)<2:
            return not self.underHalf(frames)
        return self.underHalf(frames[1:]) and not self.underHalf(frames)

    def overFull(self,frames,requiredslice):
        """ figure out if the frame list needs to be trimmed

        :param frames: list of frames
        :param requiredslice: frame required to be near the middle
        """
        if len(frames)<1:
            return False
        if self.modrequirement is not None:
            if len(frames)%self.modrequirement[0]!=self.modrequirement[1]:
                return True
        if self.width!=None and len(frames)<=int(self.width):
            return False
        if self.modrequirement is not None and self.delta_t(frames)>self.time_width:
            return self.overFull(frames[self.modrequirement[0]:],requiredslice)
        return self.delta_t(frames)>self.time_width

    def indexOf(self,frames,frame):
        """ find a given frame in the frame list

        :param frames: list of frames
        :param frame: frame to check
        :returns: index of given frame, or None if not present
        """
        wasZero=False
        rv=None
        for i,x in enumerate(frames):
            if x is frame:
                if i==0:
                    wasZero=True
                rv=i
                if not wasZero:#if it was the first, it may be copied, so return the last. if it wasn't, return the first
                    return rv
        if wasZero and rv==len(frames)-1:
            rv=int(len(frames)/2)
        return rv

    def inMiddle(self,frames,frame,pad=0):
        x=self.indexOf(frames,frame)
        lf=len(frames)+abs(pad)
        if pad<0:
            x-=pad
        hl=lf/2
        if lf%2==1:
            return x==(hl)
        return x in (hl-1,hl)

    def rangeForIdx(self,idx,subi,maxi):
        if subi==0:
            sidx=idx
            eidx=idx+1
        else:
            i=subi-1
            hw=(i/3)
            m=i%3
            sidx=idx-hw-(1 if m in (0,2) else 0)
            eidx=idx+hw+(1 if m in (1,2) else 0)+1
        if self.edgemode == 'short':
            if sidx<0:
                sidx=0
            if eidx>maxi:
                eidx=maxi
        return sidx,eidx

    def evalval(self,frames,idx,subi):
        sidx,eidx=self.rangeForIdx(idx,subi,len(frames))
        ret=0.0
        if self.width and self.width>0 and (eidx-sidx)<self.width:
            return ret
        if sidx<0 or eidx>len(frames):
            return ret
        if self.modrequirement!=None and (eidx-sidx)%self.modrequirement[0]!=self.modrequirement[1]:
            return ret
        ret+=.1
        dt=self.delta_t(frames[sidx:eidx])
        if dt<self.time_width:
            return ret
        if self.time_width==dt:
            return 1.0
        ret+=.9/(dt.total_seconds()-self.time_width.total_seconds())
        return ret

    def matrixFor(self,frames,frameidx):
        return [self.evalval(frames,frameidx,x) for x in range(2*len(frames))]

    def maxval(self,arr):
        x=arr[0]
        idx=0
        for a in range(1,len(arr)):
            if arr[a]>=x:
                idx=a
                x=arr[a]
        return idx,x

    def getWindow(self,frames,frame,canFail=True):
        idx=self.indexOf(frames,frame)
        arr=self.matrixFor(frames,idx)
        aidx,v=self.maxval(arr)
        if v==0.0:
            if canFail:
                return None
            sidx,eidx=self.rangeForIdx(idx,0,len(frames))
        else:
            sidx,eidx=self.rangeForIdx(idx,aidx,len(frames)) #this is a brute force method to check all permutations with the current buffer, assuming the required frame is near the middle
        #print 'rolling window range ',sidx,eidx,self.delta_t(frames[sidx:eidx])
        for x in range(sidx):
            del frames[0]
        return frames[:(eidx-sidx)]

    def process(self):
        myslice=[]
        slicequeue=Queue.Queue()
        inqueue=0
        haveWindow=False
        ranLast=True
        framesSeen=0
        for f in self.source:
            framesSeen+=1
            myslice.append(f)
            slicequeue.put(f)
            inqueue+=1
            if haveWindow==False and self.underHalf(myslice):
                continue
            elif haveWindow==False:
                if self.edgemode=='fullduplicate':
                    for x in range(len(myslice)):#((self.width or 0) + inqueue)*2):
                        myslice.insert(0,myslice[0])
                haveWindow=True
            if ranLast:
                ranLast=False
                requiredslice=slicequeue.get_nowait()
                inqueue-=1
            useslice=self.getWindow(myslice,requiredslice)
            #print ('filter window', len(myslice),self.underHalf(myslice),self.half(myslice),self.overFull(myslice,requiredslice),inqueue,framesSeen)
            if useslice!=None:
                ranLast=True
                yield self.docall(useslice,requiredslice)

        if not haveWindow:
                if self.edgemode=='fullduplicate':
                    for x in range(len(myslice)):#((self.width or 0) + inqueue)*2):
                        myslice.insert(0,myslice[0])
                haveWindow=True
        if inqueue>0 or not ranLast:
            if ranLast:
                requiredslice=slicequeue.get_nowait()
                inqueue-=1
                ranLast=False
            if self.edgemode=='fullduplicate' and len(myslice)<=(self.width or inqueue):#not even enough data to fill a slice
                while len(myslice)<(self.width or (2.1*inqueue)):
                    myslice.append(myslice[-1])
            try:
                while inqueue>0 or not ranLast:
                    myslice.append(myslice[-1])
                    if ranLast:
                        requiredslice=slicequeue.get_nowait()
                        inqueue-=1
                    useslice=self.getWindow(myslice,requiredslice,canFail=False)
                    #print ('filter window wrapping up',len(myslice),self.underHalf(myslice),self.half(myslice),self.overFull(myslice,requiredslice),inqueue,framesSeen)
                    yield self.docall(useslice,requiredslice)
                    ranLast=True
            except Queue.Empty:
                raise RuntimeError('Premature rolling window termination!')
        try:
            slicequeue.get_nowait()
            raise RuntimeError('Incomplete rolling window? BUG')
        except Queue.Empty:
            pass#good

