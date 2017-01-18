
import dplkit.role.filter
import dplkit.role.decorator
import copy
from datetime import datetime,timedelta
import lg_base.core.array_utils as hau
import numpy as np
import os
import json
import traceback

from collections import namedtuple

def returnFieldOfSelfFunc(fieldname):
    def myfunc(self):
        return getattr(self,fieldname)
    return myfunc

class dpl_expose_attributes_filter(dplkit.role.filter.aFilter):
    def __new__(cls,framestream,**kwargs):
        cls=dplkit.role.decorator.exposes_attrs_in_chain(kwargs.keys())(cls)
        for f in kwargs.keys():
            setattr(cls,f,property(returnFieldOfSelfFunc('_'+f)))
        return dplkit.role.filter.aFilter.__new__(cls,framestream,**kwargs)
    
    def __init__(self,framestream,**kwargs):
        super(dpl_expose_attributes_filter,self).__init__(framestream)
        self.framestream=framestream
        for k,v in kwargs.items():
            setattr(self,'_'+k,v)

    def process(self):
        for f in self.framestream:
            yield f


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
    except KeyError:
        raise KeyError(fieldname)

class dpl_configured_callable_withprovides(dplkit.role.filter.aFilter):
    """ Generic configurable mappable Dpl filter object that calls an iterable

    :param framestream: input iterable framestream
    :param callableprocessor: callable that returns the result. this callable should take parameters given below, and any additional kwargs given here (initialization)
    :param functionReturnSetTo: if None, this DPL object will yield what the callable return.  If set to a string, will attach the return to the original frame at that name
    :param framestream_attributes: dictionary for mapping parameter names to attribute names.  i.e. a value {'parameter_a':'attribute_a'} will have callable called with parameter_a set to framestream.attribute_a
    :param frame_contents: dictionary for mapping parameter names to frame contents.  if this is None, the first parameter to the callable will contain the entire frame
    :param fillMissing: if true, will set missing parts to None

    Any other args are passed to the callableprocessor in each iteration

    example:
    dpl_configured_callable(stream,myfunction,functionReturnSetTo='rs_particle',framestream_attributes=dict(sounding='hsrl_sounding'),frame_contents=dict(rs_radar='rs_mmcr'),myconfig={'source':'orange'})

    What this does:
    for each frame in 'stream', it will make the call myfunction(sounding=stream.hsrl_sounding,rs_radar=frame.rs_mmcr,myconfig={'source':'orange'})
    because of the following:
        framestream_attributes has a dictionary {'sounding':'hsrl_sounding'}, so myfunction is called with sounding set to the framestream's 'hsrl_sounding' attribute
        frame_contents has a dictionary {'rs_radar':'rs_mmcr'}, so myfunction is called with rs_radar set to the 'rs_mmcr' piece from the current frames
        myconfig isn't a normal parameter for initializing dpl_configured_callable, so it is passed as-is to myfunction
    the result returned by myfunction will be put into the frame as 'rs_particle' (i.e. frame['rs_particle'] or frame.rs_particle), and yielded becuase functionReturnSetTo is set to 'rs_particle'

    """
    def __init__(self,framestream,callableprocessor,functionReturnSetTo=None,framestream_attributes={},frame_contents=None,fillMissing=False,provides=None,**kwargs):
        super(dpl_configured_callable_withprovides,self).__init__(framestream)
        self.framestream=framestream
        self.funcret=functionReturnSetTo
        self.callable=callableprocessor
        self.framestream_attributes=framestream_attributes
        self.frame_contents=frame_contents
        self.kwargs=copy.copy(kwargs)
        self.fillMissing=fillMissing
        if provides is not None:
            self.provides=provides
        if frame_contents!=None:
            for k,field in frame_contents.items():
                if provides is None:
                    provides=framestream.provides
                try:
                    v=deepattribute(provides,field)
                except RuntimeError:
                    if not fillMissing:
                        raise RuntimeError('Framestream '+repr(self)+' from '+repr(framestream)+" doesn't provide "+field)

    def process(self):
        """ main dpl generator function
        """
        for frame in self.framestream:
            thearrayargs=[]
            theargs=copy.copy(self.kwargs)
            for k,field in self.framestream_attributes.items():
                try:
                    theargs[k]=deepattribute(self,field)
                except KeyError:
                    if self.fillMissing:
                        theargs[k]=None
                    else:
                        raise RuntimeError('Failed to get attribute ',field,' from framestream ',repr(self.framestream))
            if self.frame_contents==None:
                thearrayargs.append(frame)
            else:
                for k,field in self.frame_contents.items():
                    try:
                        theargs[k]=deepattribute(frame,field)
                    except KeyError:
                        if self.fillMissing:
                            theargs[k]=None
                        else:
                            raise

            ret=None
            try:
                ret=self.callable(*thearrayargs,**theargs)
            except:
                traceback.print_exc()
                raise

            if ret==None:
                if self.funcret!=None:
                    yield frame
                else:
                    print 'not yielding a none value in configured callable. this should be an error'
            elif self.funcret!=None:
                frame=copy.copy(frame)
                if isinstance(frame,dict):
                    frame[self.funcret]=ret
                else:
                    setattr(frame,self.funcret,ret)
                yield frame
            else:
                yield ret


class dpl_configured_callable(dpl_configured_callable_withprovides):

    def __new__(cls,*args,**kwargs):
        autoargs={}
        for k in ('reuseGenerator','nestedclasses','frameclass'):
            if k in kwargs:
                autoargs[k]=kwargs.pop(k)
        if 'nestedclasses' not in autoargs and 'frameclass' not in autoargs:
            autoargs['nestedclasses']=[hau.Time_Z_Group,hau.rs_xfer,dict]
        if 'nestedclasses' in autoargs:
            cls=dplkit.role.decorator.autoprovidenested(**autoargs)(cls)
        elif 'frameclass' in autoargs:
            cls=dplkit.role.decorator.autoprovide(**autoargs)(cls)
        else:
            raise RuntimeError('Missing classes for autoprovides!')
        return super(dpl_configured_callable,cls).__new__(cls, *args,**kwargs)
