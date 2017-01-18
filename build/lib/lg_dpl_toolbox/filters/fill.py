import dplkit.role.filter
import numpy
import lg_base.core.array_utils as hau
import dplkit.role.decorator
from datetime import timedelta,datetime
import time
import os
import traceback
import copy
import warnings

import logging
LOG = logging.getLogger(__name__)

def timeToFloat(base,arr):
    return numpy.array([(x-base).total_seconds() for x in arr])

class FillIn(dplkit.role.filter.aFilter):
    def __init__(self,stream,axes,maxvals=None,subframes=None,ignoreGroups=None):#if asClass is none, will return a dictionary, otherwise will make a new class, and use setattr. if classdictparameter is
        super(FillIn,self).__init__(stream)    
        self.stream=stream#this has provides and any exposed attributes
        self.realaxes=axes
        self.axes=[x for x in axes]
        self.base=self.axes[0][0]
        self.axes[0]=timeToFloat(self.base,self.axes[0])
        self.maxvals=maxvals
        self.subframes=subframes
        self.ignoreGroups=ignoreGroups

    def dofill(self,fr,subfr=None,seenObj=None,suff=[]):
        if seenObj is None:
            seenObj=set()
        if subfr is None:
            tmp=vars(fr).keys()
            subfr=[]
            for x in tmp:
                if x.startswith('_'):
                    continue
                subfr.append(x)
        tax=None
        aax=None
        for k in subfr:
            if not hasattr(fr,k):
                continue
            #print 'attr',k
            f=getattr(fr,k)
            if isinstance(f,hau.Time_Z_Group):
                if self.ignoreGroups and self.ignoreGroups(k):
                    continue
                if f in seenObj:
                    warnings.warn("Duplicate object seen in fill!")
                    continue
                seenObj.add(f)
                self.dofill(f,seenObj=seenObj,suff=suff+[k])
            elif isinstance(f,hau.Z_Array):
                if 'O' in str(f.dtype):
                    fill = None
                elif 'uint' in str(f.dtype):
                    fill=0
                elif 'int' in str(f.dtype):
                    fill=-1
                else:
                    fill=numpy.NAN
                if isinstance(f,hau.TZ_Array) and len(self.axes)>1:
                    if tax is None:
                        tax=timeToFloat(self.base,getattr(fr,fr._timevarname))
                        setattr(fr,fr._timevarname,self.realaxes[0])
                    if aax is None:
                        aax=getattr(fr,fr._altitudevarname)
                        setattr(fr,fr._altitudevarname,self.realaxes[1])
                    try:
                        fun=hau.nearest_nd([tax,aax],f,self.maxvals,fill)
                    except AssertionError:
                        print 'Assertion error on ',k,' grid fill'
                        print f.shape,'of variable incompatible with expected',[tax.size,aax.size]
                        raise
                    setattr(fr,k,fun(*self.axes))
                elif isinstance(f,hau.T_Array):
                    if tax is None:
                        tax=timeToFloat(self.base,getattr(fr,fr._timevarname))
                        setattr(fr,fr._timevarname,self.realaxes[0])
                    if k!=fr._timevarname:
                        fun=hau.nearest_nd([tax],f,self.maxvals[0:1] if self.maxvals else None,fill)
                        setattr(fr,k,fun(self.axes[0]))
                elif isinstance(f,hau.Z_Array) and len(self.axes)>1:
                    if aax is None:
                        aax=getattr(fr,fr._altitudevarname)
                        setattr(fr,fr._altitudevarname,self.realaxes[1])
                    if k!=fr._altitudevarname:
                        fun=hau.nearest_nd([aax],f,self.maxvals[1:2] if self.maxvals else None,fill)
                        setattr(fr,k,fun(self.axes[1]))

    def process(self):
        subfr=self.subframes
        for f in self.stream:
            try:
                if f is not None:
                    f=copy.deepcopy(f)
                    self.dofill(f,subfr)
            except:
                traceback.print_exc()
                raise
            yield f
