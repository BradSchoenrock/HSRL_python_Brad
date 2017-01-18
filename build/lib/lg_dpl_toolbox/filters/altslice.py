
import dplkit.role.filter
import numpy
import lg_base.core.array_utils as hau
import dplkit.role.decorator
from datetime import timedelta,datetime
import time
import os
import traceback
import copy

import logging
LOG = logging.getLogger(__name__)

def tupleslice(idx,val,length):
    ret=[]
    for x in range(length):
        if x==idx:
            ret.append(val)
        else:
            ret.append(slice(None))
    if len(ret)==1:
        return ret[0]
    return tuple(ret)

class AltitudeSlicing(dplkit.role.filter.aFilter):
    """ Slice altitude axis to a subset To a nested frame
    """
    def __init__(self,stream,masker=None,altvar=None):#if asClass is none, will return a dictionary, otherwise will make a new class, and use setattr. if classdictparameter is
        super(AltitudeSlicing,self).__init__(stream)    
        self.stream=stream#this has provides and any exposed attributes
        self.altvar=altvar
        self.masker=masker

    def makeMask(self,arr):
        if self.masker is None:
            return None
        return self.masker(arr)
 
    def reslice(self,dic,altvar="----"):
        if altvar=="----":
            altvar=self.altvar
        if altvar is None or altvar not in dic:
            return
        mask=self.makeMask(dic[altvar])
        if mask is None:
            return
        for k,v in dic.items():
            if isinstance(v,hau.Time_Z_Group):
                self.reslice_tzg(v)
            elif isinstance(v,dict):
                self.reslice(v)
            elif isinstance(v,hau.TZ_Array):
                dic[k]=v[tupleslice(1,mask,len(v.shape))]
            elif isinstance(v,hau.T_Array):
                pass
            elif isinstance(v,hau.Z_Array):
                dic[k]=v[tupleslice(0,mask,len(v.shape))]

    def reslice_tzg(self,structure):
        self.reslice(vars(structure),altvar=structure._altitudevarname)

    def process(self):
        for f in self.stream:
            f=copy.deepcopy(f)
            if isinstance(f,hau.Time_Z_Group):
                self.reslice_tzg(f)
            elif isinstance(f,dict):
                self.reslice(f)
            yield f

class AltitudeMask(dplkit.role.filter.aFilter):
    """ Slice altitude axis to a subset To a nested frame
    """
    def __init__(self,stream,maskname=None,masker=None,altvar=None):#if asClass is none, will return a dictionary, otherwise will make a new class, and use setattr. if classdictparameter is
        super(AltitudeMask,self).__init__(stream)    
        self.stream=stream#this has provides and any exposed attributes
        self.altvar=altvar
        self.masker=masker
        self.maskname=maskname

    def makeMask(self,arr):
        if self.masker is None:
            return None
        return self.masker(arr)

    def subslice(self,v,keys=[]):
            if isinstance(v,hau.Time_Z_Group):
                return self.reslice_tzg(v,keys=keys)
            elif isinstance(v,dict):
                return self.reslice_dic(v,keys=keys)
            else:
                try:
                    if isinstance(v,object):
                        #print 'Maybe using type ',type(v)
                        x=repr(type(v))
                        if 'class' in x and "rray" not in x:#FIXME
                        #if not isinstance(v,(basestring,function,staticmethod,tuple,list)):
                            print 'attempting to use type ',x,'at key',('.'.join(keys))
                            return self.reslice_struc(v,keys=keys)
                except:
                    pass
                #print 'Ignoring type ',type(v)
                return v

 
    def reslice_dic(self,dic,altvar="___",keys=[]):
        if self.maskname is None:
            return dic
        for k,v in dic.items():
            dic[k]=self.subslice(v,keys=keys+[k])
        if altvar=="___":
            altvar=self.altvar
        if altvar is None or altvar not in dic:
            return dic
        mask=self.makeMask(dic[altvar])
        if mask is None:
            return dic
        #print altvar,'made mask in frame','.'.join(keys),mask
        dic[self.maskname]=mask
        return dic

    def reslice_tzg(self,structure,keys=[]):
        x=vars(structure)
        #print 'alt of ','.'.join(keys),'is',structure._altitudevarname,getattr(structure,structure._altitudevarname) if hasattr(structure,structure._altitudevarname) else None
        t=self.reslice_dic(x,altvar=structure._altitudevarname,keys=keys)
        assert(t is x)
        return structure

    def reslice_struc(self,structure,keys=[]):
        x=vars(structure)
        nx=self.reslice_dic(x,keys=keys)
        assert(x is nx)
        return structure

    def process(self):
        for f in self.stream:
            if self.maskname is not None:   
                if isinstance(f,hau.Time_Z_Group):
                    f=self.reslice_tzg(f)
                elif isinstance(f,dict):
                    f=self.reslice(f)
                elif f is None:
                    pass
                else:
                    print "Content isnt a supported type "+repr(type(f))
                    f=self.reslice_struc(f)
            yield f


class RangeMask:
    def __init__(self,minalt,maxalt):
        self.minalt=minalt
        self.maxalt=maxalt

    def __call__(self,arr):
        if self.minalt is None:
            if self.maxalt is None:
                return None
            else:
                mask=(arr<=self.maxalt)
        else:
            mask=(arr>=self.minalt)
            if self.maxalt is not None:
                mask[arr>self.maxalt]=False
        return mask

class DescreteMask:
    def __init__(self,alts):
        self.alts=alts
        try:
            x=len(alts)
        except TypeError:
            self.alts=[alts]

    def __call__(self,arr):
        if self.alts is None:
            return None
        mask=(arr==0.0)
        mask[:]=numpy.array([(x in self.alts) for x in arr])
        return mask

class NearestMask:
    def __init__(self,alts):
        self.alts=alts
        try:
            x=len(alts)
        except TypeError:
            self.alts=[alts]

    def __call__(self,arr):
        if self.alts is None:
            return None
        mask=(arr==-numpy.Infinity)
        for a in self.alts:
            delta=numpy.abs(arr-a)
            #print 'deltas to',a,'vs',arr,'is',delta
            minidx=-1
            minval=numpy.Infinity
            for i,d in enumerate(delta):
                if d<minval:
                    minidx=i
                    minval=d
            if minidx>=0:
                mask[minidx]=True
        return mask

class MaskAsIndexes:
    def __init__(self,host):
        self.host=host

    def __call__(self,arr):
        mask=self.host(arr)
        mask=numpy.arange(mask.size)[mask]
        return mask

class RangeAltitudeSlicing(AltitudeSlicing):
    """ Slice altitude axis to a subset To a nested frame
    """
    def __init__(self,stream,minalt=None,maxalt=None,*args,**kwargs):#if asClass is none, will return a dictionary, otherwise will make a new class, and use setattr. if classdictparameter is
        super(RangeAltitudeSlicing,self).__init__(stream,masker=MaskAsIndexes(RangeMask(minalt,maxalt)),*args,**kwargs)

class DescreteAltitudeSlicing(AltitudeSlicing):
    """ Slice altitude axis to a subset To a nested frame
    """
    def __init__(self,stream,alts=None,*args,**kwargs):#if asClass is none, will return a dictionary, otherwise will make a new class, and use setattr. if classdictparameter is
        super(DescreteAltitudeSlicing,self).__init__(stream,masker=MaskAsIndexes(NearestMask(alts)),*args,**kwargs)    

class RangeAltitudeMask(AltitudeMask):
    """ Slice altitude axis to a subset To a nested frame
    """
    def __init__(self,stream,minalt=None,maxalt=None,*args,**kwargs):#if asClass is none, will return a dictionary, otherwise will make a new class, and use setattr. if classdictparameter is
        super(RangeAltitudeMask,self).__init__(stream,masker=MaskAsIndexes(RangeMask(minalt,maxalt)),*args,**kwargs)

class DescreteAltitudeMask(AltitudeMask):
    """ Slice altitude axis to a subset To a nested frame
    """
    def __init__(self,stream,alts=None,*args,**kwargs):#if asClass is none, will return a dictionary, otherwise will make a new class, and use setattr. if classdictparameter is
        super(DescreteAltitudeMask,self).__init__(stream,masker=MaskAsIndexes(NearestMask(alts)),*args,**kwargs)    
