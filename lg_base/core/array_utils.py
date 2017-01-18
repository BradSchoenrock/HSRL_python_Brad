#!/usr/bin/python
# -*- coding: utf-8 -*-

""" The time_z_data structure must contain a vector "time_z_data.times"
    containing time entries express in python datetimes with
    dimensions [ntimes,0]. It may also contian an altitude vector
    "time_z_data.altitudes" with dimensions [0,nalts], Other arrays
    must have individual entries corresponding to these times and
    altitudes, with columuns refering to time and rows to altitudes.
    Entries with time and altitude dimensions must be of subclass "tz_array".
    Entries with time but no altitude dimension must be of subclass "t_array".
    Entries with altitude but no time dimension must be of subclass "z_array"
    Items that are not members of these subclasses may be included--they
    are passed through routines without change."""


"""  editor substitutions:

%s/tzu.t_array.from_array/T_Array/g
%s/tzu.tz_array.from_array(/TZ_Array(/g
%s/tzu.time_z_data()/Time_Z_Group()/g
%s/tzu.z_array.from_array/Z_Array/

rs_tail = tzu.trim_time_interval(rs_tail, data_end_time, np.inf)
rs_tail.trimTimeInterval(data_end_time, np.inf)


rs_mem = tzu.append_time_z_data_object(rs_tail, rs_mem)
rs_mem.append(rs_tail)

r_ms = tzu.nanmean_time_z_data_object(r_ms, n_time_ave, 1)
r_ms.binMean( n_time_av, 1 ) 
"""
import logging
import warnings

LOG = logging.getLogger(__name__)

import numpy as np
from datetime import timedelta,datetime
from ExtendedArray import ExtendedArray
import copy
# Try to use the much faster nanmean from bottleneck, otherwise fall back
# to the scipy.stats version

try:
    from bottleneck import nanmean,nansum,anynan
except ImportError:
    print
    print 'No bottleneck.nanmean available! Falling back to SLOW scipy.stats.nanmean'
    print
    from scipy.stats import nanmean
    from numpy import nansum
    def anynan(x):
        return np.any(np.isnan(x))  

def allbut(ax,idx,count):
    ret=[]
    for x in range(count):
        if x==ax:
            ret.append(idx)
        else:
            ret.append(slice(None))
    if count==1:
        return ret[0]
    return tuple(ret)

class nearest_nd(object):
    def __init__(self,dims,arr,maxddims=None,fill=np.NAN,copy=False):
        if copy:
            import copy
            self.dims=copy.deepcopy(dims)
            self.arr=arr.copy()
        else:
            import copy
            self.dims=copy.copy(dims)
            self.arr=arr
        self.maxddims=maxddims
        self.fill=fill
        assert(len(self.dims)<=len(self.arr.shape))
        assert(self.maxddims is None or len(self.dims)==len(self.maxddims))
        for i,n in enumerate(self.dims):
            assert(len(n)==arr.shape[i])

    def newcontent(self,arr,olddim,newdim,axis,dmax,dfill):
        idx=0
        for n in newdim:
            while idx<(olddim.size-1) and abs(olddim[idx]-n)>abs(olddim[idx+1]-n):
                idx+=1
            if dmax is not None and abs(olddim[idx]-n)>dmax:
                yield dfill
            else:
                yield arr[allbut(axis,idx,len(arr.shape))]

    def __call__(self,*args,**kwargs):
        axes=args
        assert(len(axes)<=len(self.dims))
        if not kwargs.pop('assume_sorted',False):
            axes=[x.copy() for x in args]
            for x in axes:
                x.sort()
        preret=self.arr
        axcount=len(preret.shape)
        for ax in range(len(axes)):
            sz=list(preret.shape)
            sz[ax]=len(axes[ax])
            maxdim=self.maxddims[ax] if self.maxddims is not None else None
            if maxdim is None:
                maxdim=abs(2*np.median(np.diff(self.dims[ax])))
            ret=np.ones(sz,dtype=self.arr.dtype)
            if type(self.arr)!=np.ndarray:
                ret=type(self.arr)(ret,like=self.arr)
            for i,nx in enumerate(self.newcontent(preret,self.dims[ax],axes[ax],ax
                    ,maxdim,self.fill)):
                ret[allbut(ax,i,axcount)]=nx
            preret=ret

        return ret

def safeIsInstance(obj, klass):
    """a simple alternative to 'isinstance()' based on parsing the string
       representation of type().  E.g.
       parse  "<class 'hsrl_array_utils.Z_Array'>" into 'Z_Array'
       since

       the builtin 'isinstance(obj, class)' seems confused by "import x.y.z.class"
       vs "import z.class"
    """
    s = repr(type(obj))
    obj_type = s[s.rfind('.')+1:s.rfind("'")]
    ret = (obj_type == klass.__name__)

    ret2=type(obj)==klass
    if ret!=ret2:
        print 'CLASS INSTANCE MISMATCH ',type(obj),'vs',klass
    return ret

def keepnansum(arr,*args,**kwargs):
    mask=np.ones_like(arr)
    try:
        mask[np.isnan(arr)]=0
    except TypeError:
        print 'isnan doesn\'t like',arr.dtype
        pass#return np.sum(arr,*args,**kwargs)
    summask=nansum(mask,*args,**kwargs)
    sumarr=None
    try:
        sumarr=nansum(arr,*args,**kwargs)
        sumarr[summask==0]=np.NaN
    except TypeError:
        if sumarr is None:
            warnings.warn("variable doesn't support nansum!")
            return np.sum(arr,*args,**kwargs)
        pass
    return sumarr

def nanor(arr,axis=None):
    if len(arr.shape)>3:
        raise NotImplementedError('size = '+repr(arr.shape))
    return _nanOper(np.bitwise_or,arr,axis if axis is not None else -1,0)

def nanand(arr,axis=None):
    if len(arr.shape)>3:
        raise NotImplementedError('size = '+repr(arr.shape))
    return _nanOper(np.bitwise_and,arr,axis if axis is not None else -1,0)

#@jit
def _nanOper(oper,arr,axis,filler):
    if arr.size==0:
        return filler
    if axis<0 or (arr.shape[axis]==arr.size):
        tmp=arr.ravel()
        r=tmp[0]
        for x in tmp[1:]:
            r=oper(r,x)
        return r
    sh=arr.shape
    newsh=list(sh)
    newsh[axis]=1
    ret=np.zeros(newsh,dtype=arr.dtype)
    idx=[]
    retidx=[]
    fullidx=[]
    for x in range(len(sh)):
        idx.append(np.arange(newsh[x]))
        retidx.append(np.arange(newsh[x]))
        fullidx.append(slice(None))
    idx=tuple(idx)
    retidx=tuple(retidx)
    fullidx=tuple(fullidx)
    ret[retidx]=arr[retidx]
    for x in range(1,sh[axis]):
        idx[axis][0]=x
        ret[fullidx]=oper(ret,arr[idx])

    return type(arr)(ret,dtype=arr.dtype,summode=arr.summode)


class rs_xfer(object): # this is used to avoid appending blindly as if its a Time_Z_Group

    pass
   
def maybetuple(t):
    if len(t)==0:
        return None
    if len(t)==1:
        return t[0]
    return tuple(t)

def nanfirst(arr,axis=0,filler='notaval'):
    idx=[]
    siz=[x for x in arr.shape]
    for x in range(len(arr.shape)):
        idx.append(slice(None,None))
    idx[axis]=np.arange(1)
    siz[axis]=1
    if filler!='notaval':
        bval=filler
    elif arr.dtype=='object':
        bval=None
    elif arr.dtype in ['float','float16','float32','float64','complex','complex64','complex128']:
        bval=np.NaN
    else:
        bval=0
    if arr.shape[axis]==0:
        ret=copy.deepcopy(arr)
        ret.resize(siz)
        ret[maybetuple(idx)]=bval
    else:
        ret=arr[maybetuple(idx)].copy()
        idx[axis]+=1
        while np.any(ret==bval) and idx[axis]<arr.shape[axis]:
            ret[ret==bval]=arr[maybetuple(idx)][ret==bval]
            idx[axis]+=1
    return ret


validSumModes=('mean','sum','and','or','first','none')

class Z_Array(ExtendedArray):

    """an array which ordered by altitude """

    def __new__(
        subtype,
        data,
        info=None,
        dtype=None,
        copy=False,
        summode=None,
        like=None
        ):

        assert(not(like is not None and (summode is not None or info is not None) ))
        ret = super(Z_Array,subtype).__new__(subtype, data, info=info.copy() if info is not None else (like.info.copy() if like is not None else None),
                dtype=dtype, copy=copy)
        if summode is not None:
            if summode not in validSumModes:
                raise RuntimeError('Invalid summode '+summode)
            ret.summode=summode
        if 'summode' not in ret.info:
            if ret.dtype=='object':
                ret.summode='none'
            else:
                ret.summode='mean'
        mydtype=ret.dtype
        #if (ret.sumAnd or ret.sumOr):
        #    mydtype='uint64'
        if (ret.sumSum or ret.sumMean) and ret.dtype not in ['float','float16','float32','float64','complex','complex64','complex128']:
            mydtype='float64'
        #if ret.dtype=='object' and summode==None:
        #    mydtype='object'
        if mydtype!=ret.dtype:
            osval=ret.summode
            ret = super(Z_Array,subtype).__new__(subtype, data, info=info.copy() if info is not None else info,
                dtype=mydtype, copy=copy)
            ret.summode=osval

        #print 'obj',subtype,'instance has info',ret.info
        if (ret.sumAnd or ret.sumOr) and ret.dtype not in ['bool','int','uint','int16','uint16','int32','uint32','int64','uint64']:
            raise TypeError('Floating bitfield. type is ',ret.dtype)
        if (ret.sumSum or ret.sumMean) and ret.dtype not in ['float','float16','float32','float64','complex','complex64','complex128']:
            raise TypeError('Fixed meanable. type is ',ret.dtype)
        #if (ret.sumAnd or ret.sumOr) and ret.dtype in ['int','int16','int32','int64']:
        #    print 'Recommend using unsigned for bitfields!'
        return ret

    @property
    def sumMean(self):
        return self.summode=='mean'
    @property
    def sumSum(self):
        return self.summode=='sum'
    @property
    def sumFirst(self):
        return self.summode=='first'
    @property
    def sumAnd(self):
        return self.summode=='and'
    @property
    def sumOr(self):
        return self.summode=='or'

    @property
    def sumpower(self):
        if 'sumpower' not in self.info:
            return 0
        return self.info['sumpower']

    @sumpower.setter
    def sumpower(self,val):
        if 'sumpower' in self.info and self.info['sumpower']==val:
            return
        self.info=copy.copy(self.info)
        self.info['sumpower']=val

    @property
    def summode(self):
        if 'summode' not in self.info:
            if 'sumpower' in self.info:
                if self.sumpower==0:
                    return 'mean'
                else:
                    return 'sum'
            #print 'WARNING: SUMMODE not set on an array'
            return 'mean'
        if self.info['summode']=='sum' and self.sumpower==0:
            return 'mean'
        return self.info['summode']

    @summode.setter
    def summode(self, value):
        if value not in validSumModes:
            raise RuntimeError('Invalid summode '+value)
        if self.summode==value:
            return
        self.info=copy.copy(self.info)
        if value=='sum':
            self.sumpower=1
        elif value=='mean':
            value='sum'
            self.sumpower=0
        elif 'sumpower' in self.info:
            del self.info['sumpower']
        self.info['summode'] = value

    def doSum(self,*args,**kwargs):
        if self.sumFirst:
            return nanfirst(self,*args,**kwargs)
        if self.sumSum:
            return keepnansum(self,*args,**kwargs)
        if self.sumAnd:
            return nanand(self,*args,**kwargs)
        if self.sumOr:
            return nanor(self,*args,**kwargs)
        if self.sumMean:
            return nanmean(self,*args,**kwargs)
        raise RuntimeError('Unknown summode '+self.summode)

    def __degrademean(self,other):
        #if self.sumMean and isinstance(other,Z_Array) and other.sumSum:
        #    raise RuntimeError('Dividing meaned objects by summed objects is Bad')
        oldpower=self.sumpower
        if isinstance(other,Z_Array) and other.sumpower!=0:
            self.sumpower-=other.sumpower
            if self.sumpower==oldpower:
                raise RuntimeError
            if self.sumpower!=oldpower and self.sumpower*oldpower==0:#self.sumSum and isinstance(other,Z_Array) and other.sumSum:
                LOG.debug('Degrading from %g to %g' %(oldpower,self.sumpower))
                #if self.sumpower==-1:
                #    traceback.print_stack()
                #self.summode='mean'
        return self
    def __upgrademean(self,other):
        #if self.sumSum and isinstance(other,Z_Array) and other.sumSum:
        #    raise RuntimeError('Multiplying summed objects is Bad')
        oldpower=self.sumpower
        if isinstance(other,Z_Array) and other.sumpower!=0:
            self.sumpower+=other.sumpower
            if self.sumpower==oldpower:
                raise RuntimeError
            if self.sumpower!=oldpower and self.sumpower*oldpower==0:
                LOG.debug('Upgrading from %g to %g' %(oldpower,self.sumpower))
                #traceback.print_stack()
                #self.summode='sum'
        return self
    def __powermean(self,other):
        #if self.sumSum and isinstance(other,Z_Array) and other.sumSum:
        #    raise RuntimeError('Multiplying summed objects is Bad')
        oldpower=self.sumpower

        if self.sumpower!=0:#isinstance(other,Z_Array) and other.sumpower!=0:
            self.sumpower=float(self.sumpower)*other

            if self.sumpower==oldpower:
                raise RuntimeError
            if self.sumpower!=oldpower:# and self.sumpower*oldpower==0:
                LOG.debug( 'Powering from %g to %g' %(oldpower,self.sumpower))
                #traceback.print_stack()
                #self.summode='sum'
        return self

    def __div__(self,other):
        return super(Z_Array,self).__div__(other).__degrademean(other)
    def __floordiv__(self,other):
        return super(Z_Array,self).__floordiv__(other).__degrademean(other)
    def __truediv__(self,other):
        return super(Z_Array,self).__truediv__(other).__degrademean(other)
    def __idiv__(self,other):
        return super(Z_Array,self).__idiv__(other).__degrademean(other)
    def __ifloordiv__(self,other):
        return super(Z_Array,self).__ifloordiv__(other).__degrademean(other)
    def __itruediv__(self,other):
        return super(Z_Array,self).__itruediv__(other).__degrademean(other)
    def __mul__(self,other):
        return super(Z_Array,self).__mul__(other).__upgrademean(other)
    def __imul__(self,other):
        return super(Z_Array,self).__imul__(other).__upgrademean(other)
    def __pow__(self,other):
        return super(Z_Array,self).__pow__(other).__powermean(other)
    def __ipow__(self,other):
        return super(Z_Array,self).__ipow__(other).__powermean(other)

    def __repr__(self):

        selfEmpty = True
        nshots = 0

        # self.times.shape produces an error when self is empty,
        # this allows appending onto an empty array.
        try:
            nshots=self.shape[0]
            selfEmpty = (nshots == 0)
        except AttributeError:
            pass
        except IndexError:
            return '%s' % np.array(self) #dimensionless... cant index... just do what array would do. whatever that is
        if selfEmpty:
            ss='Empty'
        elif nshots<=1:
            ss='[%s], %i step' % (self[0],nshots)
        else:
            dif=self[-1]-self[0]
            ss='[%s, %s], %i steps, ' % (self[0],self[-1],nshots)
            try:
              if hasattr(dif,'total_seconds'):
                ss+='%.1f seconds avg per' % (dif.total_seconds()/(nshots-1))
              else:
                ss+='%.1f avg per' % (dif/(nshots-1))
            except TypeError:
                pass
        return '<%s %s>' % (type(self).__name__.split('.')[-1],ss)


class T_Array(Z_Array):

    """a 2-d array ordered by time """

    def __new__(
        subtype,
        data,
        info=None,
        dtype=None,
        copy=False,
        summode=None,
        like=None
        ):

        return super(T_Array,subtype).__new__(subtype, data, info=info,
                dtype=dtype, copy=copy, summode=summode,like=like)

    def check_type(self):
        #print 'check_type(T_Array), type = ' ,type(self)
        return safeIsInstance(self,T_Array)

    def extend(self,bins,filler_value=None):
        self.__realextend(bins,filler_value)

    def __realextend(self,bins,filler_value):
        if bins==0:
            return
        originalshape=list(self.shape)
        ax=copy.deepcopy(originalshape)
        ax[0]=ax[0]+bins
        self.resize(ax,refcheck=False)
        if filler_value==None:
            if self.dtype in ['float','float16','float32','float64','complex','complex64','complex128']:
                filler_value=np.NaN
            else:
                filler_value=0
        newind=[slice(originalshape[0],None)]
        for x in range(1,len(ax)):
            newind.append(slice(None))
        newind=tuple(newind)
        self[newind]=filler_value


    def append(self,other,filler_value=np.NaN):
        if len(other.shape)==0 or other.shape[0]==0:
            return
        originalshape=list(self.shape)
        othershape=list(other.shape)
        if len(originalshape)==0:
            originalshape=othershape*0
        if len(originalshape)>len(othershape):
            raise RuntimeError('source '+repr(originalshape)+' doesn"t match append '+repr(othershape))
        self.__realappend(other,filler_value)

    def __realappend(self,other,filler_value):
        if len(other.shape)==0 or other.shape[0]==0:
            return
        originalshape=list(self.shape)
        othershape=list(other.shape)
        if len(originalshape)==0:
            originalshape=othershape*0
        #if len(originalshape)>len(othershape):
        #    raise RuntimeError
        while len(originalshape)<len(othershape):
            originalshape.append(1)
        #print self.shape
        #print other.shape
        ax=[int(originalshape[0]+othershape[0])]
        for x in range(1,len(originalshape)):
            ax.append(int(max(originalshape[x],othershape[x])))
        self.resize(ax,refcheck=False)

        if len(ax)==1:
            self[originalshape[0]:]=other[:]
        elif len(ax)==2:
            if originalshape[1]<ax[1]:
                self[:originalshape[0],originalshape[1]:]=filler_value
            if othershape[1]<ax[1]:
                self[originalshape[0]:,othershape[1]:]=filler_value
            self[originalshape[0]:,:othershape[1]]=other[:,:]

def verifyNew(val,seenvars,trace,doerror=True):
    n='.'.join(trace)
    for op,x in seenvars.items():
        if x is val and x is not None:
            v='Duplicate at '+n+' and '+op+' = '+repr(val)
            #print v
            if doerror:
                raise RuntimeError(v)
    #print 'fresh object at '+n
    seenvars[n]=val

def verifyNewRecursive(val,seenvars,trace,doerror=True):
    if val is None or not (isinstance(val,dict) or hasattr(val,'__dict__')):
        return
    verifyNew(val,seenvars,trace,doerror)
    if isinstance(val,np.ndarray):
        return
    iterx=None
    if not isinstance(val,dict):
        val=vars(val)
    iterx=val.items()
    for k,v in iterx:
        verifyNewRecursive(v,seenvars,trace+[k],doerror)


def verifyAllNew(**kwargs):
    reuseargs=dict(seenvars={})
    for v in ('doerror','seenvars'):
        if v in kwargs:
            reuseargs[v]=kwargs.pop(v)
    basetrace=kwargs.pop('trace',[])
    if isinstance(basetrace,basestring):
        basetrace=[basetrace]
    for k,x in kwargs.items():
        verifyNewRecursive(x,trace=basetrace+[k],**reuseargs)

class TZ_Array(T_Array):

    """a 2-d array which contains time information in the
    first dimension and altitude information in the
    second dimension"""

    def __new__(
        subtype,
        data,
        info=None,
        dtype=None,
        copy=False,
        summode=None,
        like=None
        ):

        return super(TZ_Array,subtype).__new__(subtype, data, info=info,
                dtype=dtype, copy=copy,summode=summode,like=like)

    def check_type(self, obj):
        #print 'type = ', type(obj)
        #print 'safeIsInstance(obj) = ', safeIsInstance(obj,TZ_Array)
        return safeIsInstance(obj,TZ_Array)

    def check_type2(self):
        #print 'check safeIsInstance 2'
        return safeIsInstance(self,TZ_Array)

class DataCollection(object):

    def __init__(self):
        """create an empty dataCollection"""
       

    def __repr__(self):
        """return a printable representation of our contents"""

        return repr(vars(self))

    def mergeTimeVectors(self, source):
        """copy time vectors and objects other than ZArray or TZArray to self from source"""

        for (name, value) in vars(source).items():

            # skip some objects
            if isinstance(value,TZ_Array):
                continue
            elif isinstance(value,T_Array):
                pass #dont skip
            elif isinstance(value, Z_Array):
                continue#skip if its a Z array and not a T array (pure Z) or is a TZ array
            #vars(self)[name] = copy.deepcopy(value)
            setattr(self,name,copy.deepcopy(value))
 

class Time_Z_Group(DataCollection):

    """ objects of this class contain a vector "Time_Z_Group.times"
    containing time entries express in python datetimes with
    dimensions [ntimes]. It may also contain an altitude vector
    "Time_Z_Group.altitudes" with dimensions [nalts], Other arrays
    must have individual entries corresponding to these times and
    altitudes, with columuns refering to time and rows to altitudes.
    Entries with time and altitude dimensions must be of subclass "TZ_Array".
    Entries with time but no altitude dimension must be of subclass "T_Array".
    Entries with altitude but no time dimension must be of subclass "Z_Array"
    Items that are not members of these subclasses may be included--they
    are passed through routines without change."""

    def __init__(self, times=None,timevarname=None,altname=None,can_append=None,widthvarname=None,like=None,startvalues={}):
        self._can_append=can_append
        self._timevarname=timevarname
        self._altitudevarname=altname
        self._widthvarname=widthvarname
        if like is not None:
            for v in ('_can_append','_timevarname','_altitudevarname'):
                if hasattr(like,v):
                    setattr(self,v,getattr(like,v))
        for k,v in startvalues.items():
            if not hasattr(self,k) or getattr(self,k) is None:
                setattr(self,k,v)
        if self._can_append is None:
            #print '*** WARNING: no timevarname selected in init. this is bad for consistenct, and may raise in the future. setting to deprecated default "times"'
            self._can_append=True
        if self._timevarname is None:
            #print '*** WARNING: no timevarname selected in init. this is bad for consistenct, and may raise in the future. setting to deprecated default "times"'
            self._timevarname='times'
        if self._altitudevarname is None:
            #print '*** WARNING: no altname selected in init. this is bad for consistenct, and may raise in the future. setting to deprecated default "msl_altitudes"'
            self._altitudevarname='msl_altitudes'
        if times is None:
            times=T_Array([],summode='first')#this also might go away too
        if not hasattr(self,self._timevarname):
            setattr(self,self._timevarname,times)

    def showDebug(self,prefix=''):
        for (name,value) in vars(self).items():
            if isinstance(value,Z_Array):
                LOG.debug( prefix+name,'=',type(value),value.dtype,value.summode,value.sumpower )#,value
            elif isinstance(value,Time_Z_Group):
                value.showDebug(prefix+name+'.')
            else:
                LOG.debug( prefix+name,'=',type(value))
        
    def binMean( self, kt, ka, complete_only=True):
        """compute means bin number space 
....   kt= number of times to ave together
....   ka= number of altitudes to ave together
....   warning--no treatment of remainder if ka or kt does not evenly divide buffer
....   returns smaller array with average values that properly treat NaN's"""

   
        willTRem=False
        willKRem=False
        if getattr(self,self._timevarname).size%kt!=0:
            if complete_only:
                warnings.warn("**** KT in binmean has nonzero mod. remainder will be lost")
            else:
                willTRem=True
        if hasattr(self,self._altitudevarname) and getattr(self,self._altitudevarname).size%ka!=0:
            if complete_only:
                warnings.warn("**** KA in binmean has nonzero mod. remainder will be lost")
            else:
                willKRem=True

    # iterate through all the items in a dictionary

        for (name, value) in vars(self).items():
            #print name, type(value)
            if isinstance(value, T_Array): #includes TZ
                # if time vector
                if len(value.shape)==1 and value.shape[0]>0 and isinstance(value[0],datetime):
                        temp = []
                        for k in range(0,value.shape[0],kt):#while kt * (k + 1) <= value.shape[0]:
                            if False:
                                bt=value[k]
                                meansec=nanmean([(x-bt).total_seconds() if x!= \
                                                 None else np.nan for x in value[k:(kt+k)]])
                                temp.append(bt+timedelta(seconds=meansec))
                            else:
                                bt=None
                                for x in value[k:(kt+k)]:
                                    if x is not None:
                                        bt=x
                                        break
                                temp.append(bt)
                        temp=np.array(temp)
                elif kt>1:
                    desti=[slice(0,0)]
                    srci=[slice(0,0)]
                    newshape=[int(value.shape[0]/kt)]
                    if willTRem:
                        newshape[0]=newshape[0]+1
                    for x in range(1,len(value.shape)):
                        desti.append(slice(None))
                        srci.append(slice(None))
                        newshape.append(value.shape[x])
                    temp = np.zeros(newshape)
                    for k in range(0,value.shape[0],kt):#while kt * (k + 1) <= value.shape[0]:
                        desti[0]=slice(k/kt,(k/kt)+1)
                        srci[0]=slice(k,(k + kt))
                        temp[maybetuple(desti)] = value[maybetuple(srci)].doSum(axis=0)
                    if willTRem:
                        k=value.shape[0]-(value.shape[0]%kt)
                        desti[0]=slice(k/kt,(k/kt)+1)
                        srci[0]=slice(k,None)
                        temp[maybetuple(desti)] = value[maybetuple(srci)].doSum(axis=0)
                else:
                    temp=copy.deepcopy(value)
                if not isinstance(value,TZ_Array):
                    #vars(self)[name] = T_Array(temp)
                    setattr(self,name,type(value)(temp,summode=value.summode,dtype=value.dtype))
                else:
                    # proceed with altitude mean
                    def myslice(altr):
                        ret=[slice(None)]
                        ret.append(altr)
                        while len(ret)<len(value.shape):
                            ret.append(slice(None))
                        return tuple(ret)
                    if ka >1:  #if altitude averaging is requested
                        temp2 = np.zeros((int(value.shape[0] / kt)+(1 if willTRem else 0),
                                          int(value.shape[1] / ka)+(1 if willKRem else 0)))
                        for k in range(0,nalts,ka):#while ka * (k + 1) < nalts:
                            temp2[myslice(k/ka)] = nanmean(temp[myslice(slice(k,(k+ka)))], 1)
                        if willKRem:
                            k=nalts-(nalts%ka)
                            temp2[myslice(k/ka)] = nanmean(temp[myslice(slice(k,None))], 1)                            
                    else:
                        temp2=temp
                        #vars(self)[name] = TZ_Array(temp2)
                    setattr(self,name,type(value)(temp2,summode=value.summode,dtype=value.dtype))
            elif isinstance(value, Z_Array):
                if ka>1:
                    temp = np.zeros(int(value.shape[0] / ka)+(1 if willKRem else 0),)
                    #if value.shape[0]%ka!=0:
                    #    print "WARNING **** KA in binmean has nonzero mod. remainder will be lost"
                    for k in range(0,value.shape[0],ka):##while ka * (k + 1) <= value.shape[0]:
                        temp[k/ka] = nanmean(value[k:(k+ka)])
                    if willKRem:
                        temp[k/ka + 1] = nanmean(value[(k+ka):])
                else:
                    temp=copy.deepcopy(value)
                #vars(self)[name] = Z_Array(temp)
                setattr(self,name,Z_Array(temp,summode=value.summode,dtype=value.dtype))

 
    def __repr__(self):

        selfEmpty = True
        nshots = 0

        # self.times.shape produces an error when self is empty,
        # this allows appending onto an empty array.
        try:
            nshots=getattr(self,self._timevarname).shape[0]
            selfEmpty = (nshots == 0)
        except AttributeError:
            pass 
        if selfEmpty:
            if hasattr(self,'start'):
                ss="START %s" % str(self.start)
                if hasattr(self,'width'):
                    ss+=" - %s" % str(self.width)
            else:
                ss='Empty'
        else:
            ss='TIME %s' % repr(T_Array(getattr(self,self._timevarname)))

        if hasattr(self,self._altitudevarname):
            ss+=', ALT %s' % repr(Z_Array(getattr(self,self._altitudevarname)))

        if not self._can_append:
            ss+=', non-appendable'

        return '<Time_Z_Group %s>' % (ss)


    def trimByAltitude(
        D,
        min_alt,
        max_alt,
        ):
        """trim all the fields in a dictionary to [min_alt>=D.altitudes<=min_alt]"""

        nalts = getattr(D,D._altitudevarname).shape[0]
        indices = np.arange(nalts)
        alt_mask = indices[np.array([ (x>=min_alt and x<=max_alt) for x in getattr(D,D._altitudevarname)[ :] ])].copy()
        # using an actual mask is much cleaner, in that it doesn't assume content. empty sets will cause an exception, this causes an empty set. old code returned a null, uninitialized set
        # iterate through all the items in a dictionary

        for name in vars(D).keys():
            value=getattr(D,name)#vars(D)[name]
          # only trim arrays with altitude components
            if isinstance(value, Time_Z_Group):
                value.trimByAltitude(min_alt,max_alt)
            elif indices.shape[0]==alt_mask.shape[0]:
                continue
            elif indices.shape[0]==0: #empty mask makes python break...
                                      #explicit typing as boolean, despite empty?
                pass #vars(R)[name] = copy.deepcopy(value)
            elif isinstance(value, TZ_Array):
                sl=[slice(None) for x in range(len(value.shape))]
                sl[1]=alt_mask
                setattr(D,name,value[tuple(sl)])#vars(D)[name] = value[:, alt_mask ]
            elif isinstance(value,T_Array):
                pass #T array has no alt slicing
            elif isinstance(value, Z_Array):
                setattr(D,name,value[alt_mask])#vars(D)[name] = value[ alt_mask ]
            #else:
            #    vars(D)[name] = copy.deepcopy(value)

    def trimTimeInterval( self, start_time=None, end_time=None,width=None):
        """trim all the fields in a dictionary to [start_time>=getattr(D,D._timevarname)<=end_time]"""

        times=getattr(self,self._timevarname)
        if times.size==0 or not self._can_append:
            indices=np.arange(times.shape[0])
            time_mask=indices
        else:
            if width is not None:
                end_time=times[-1]
                start_time=end_time-width
                end_time=None
            indices=np.arange(times.shape[0])
            if start_time==None:
                time_mask= indices[times<end_time]
            elif end_time==None:
                time_mask= indices[times>start_time]
            else:
                time_mask = indices[np.logical_and(times>=start_time,times<end_time)]
        dowork=(time_mask.size!=indices.size) and self._can_append
        ret=dowork
        #print 'time_mask',time_mask
        # iterate through all the items in a dictionary

        for (name, value) in vars(self).items():
            
            #print 'trim ----',name,type(value),value
            if isinstance(value, Time_Z_Group):
                ret=value.trimTimeInterval(start_time,end_time,width) or ret
            elif not dowork:
                pass
            elif isinstance(value,T_Array):#includes tz_array   safeIsInstance(value, T_Array):
                #print 'hau.trim', name, value.shape,indices.shape
                if not indices.shape[0] == value.shape[0]:
                    print 'TRIM ERROR----',name,'  has times with missing values'
                    print 'ntimes = ',indices.shape[0],'  length ',name,'=',value.shape[0]
                    raise RuntimeError('Time length mismatch on field '+name)
                if len(value.shape)==1:
                    setattr(self,name,value[ time_mask ])
                else:
                    r=[time_mask]
                    for x in range(1,len(value.shape)):
                        r.append(slice(None))
                    setattr(self,name,value[tuple(r)])
            else:
                pass#vars(R)[name] = copy.deepcopy(value)
        return ret

    def trimTimeAxis( self, arange):
        """trim all the fields in a dictionary"""

        times=getattr(self,self._timevarname)
        indices=np.arange(times.shape[0])
        time_mask=np.array(arange,dtype=int)

        dowork=(time_mask.size!=indices.size) and self._can_append
        ret=dowork
        #print 'time_mask',time_mask
        # iterate through all the items in a dictionary

        for (name, value) in vars(self).items():
            
            #print 'trim ----',name,type(value),value
            if isinstance(value, Time_Z_Group):
                ret=value.trimTimeAxis(arange) or ret
            elif not dowork:
                pass
            elif isinstance(value,T_Array):#includes tz_array   safeIsInstance(value, T_Array):
                #print 'hau.trim', name, value.shape,indices.shape
                if not indices.shape[0] == value.shape[0]:
                    print 'TRIM ERROR----',name,'  has times with missing values'
                    print 'ntimes = ',indices.shape[0],'  length ',name,'=',value.shape[0]
                    raise RuntimeError('Time length mismatch on field '+name)
                if len(value.shape)==1:
                    setattr(self,name,value[ time_mask ])
                else:
                    r=[time_mask]
                    for x in range(1,len(value.shape)):
                        r.append(slice(None))
                    setattr(self,name,value[tuple(r)])
            else:
                pass#vars(R)[name] = copy.deepcopy(value)
        return ret

    def getSliceTimeAxis( self, arange):
        """trim all the fields in a dictionary"""

        ret=Time_Z_Group(like=self)

        times=getattr(self,self._timevarname)
        time_mask=np.array(arange,dtype=int)

        #print 'time_mask',time_mask
        # iterate through all the items in a dictionary

        for (name, value) in vars(self).items():
            
            #print 'trim ----',name,type(value),value
            if isinstance(value, Time_Z_Group):
                setattr(ret,name,value.getSliceTimeAxis(arange))
            elif isinstance(value,T_Array):#includes tz_array   safeIsInstance(value, T_Array):
                #print 'hau.trim', name, value.shape,indices.shape
                if not times.shape[0] == value.shape[0]:
                    print 'TRIM ERROR----',name,'  has times with missing values'
                    print 'ntimes = ',times.shape[0],'  length ',name,'=',value.shape[0]
                    raise RuntimeError('Time length mismatch on field '+name)
                if len(value.shape)==1:
                    setattr(ret,name,value[ time_mask ])
                else:
                    r=[time_mask]
                    for x in range(1,len(value.shape)):
                        r.append(slice(None))
                    setattr(ret,name,value[tuple(r)])
            else:
                setattr(ret,name,value)
        return ret

    def iterateAllTimes(self,count=1):

        times=getattr(self,self._timevarname)
        for x in range(0,times.shape[0],count):
            yield self.getSliceTimeAxis(np.arange(x,min(x+count,times.shape[0])))

    def meanAltitudes( D, new_alts):
        """
        create altitude means of time_z object centered on new altitude vector
        D= time_z_object
        new_alts = altitude vector, center nameans of D on these altitudes
        returns smaller array with average values that properly treat NaN's
        It also checks for D.nave and multiplies entries according to the
        number of alt bins averaged together at that alititude"""

        #alts = np.zeros_like(new_alts)
        #alts[1:] = (new_alts[1:] + new_alts[:-1]) / 2
        
        n_alts=new_alts.shape[0]
        bin_edges = np.zeros(n_alts+1)
        bin_edges[0] = new_alts[0] - (new_alts[1]-new_alts[0])/2.0
        bin_edges[-1]= new_alts[-1] +(new_alts[-1]-new_alts[-2])/2.0
        bin_edges[1:n_alts] = (new_alts[:-1] + new_alts[1:n_alts]) / 2.0

        alts=getattr(D,D._altitudevarname)
        ind = np.arange(alts.shape[0])

    # make a list of indices to indicate which elements to average for each entry
    # in new alts vector

        #temp = alts#getattr(D,D._altitudevarname)[ :]
        indices = []
        inxa=np.arange(1)
        setattr(D,'n_ave_altitude',np.zeros(bin_edges.shape[0]-1))
        for k in range(bin_edges.shape[0]-1):
            #inxa = ind[temp > bin_edges[ k]]
            #temp2 = getattr(D,D._altitudevarname)[ inxa]
            #inx = inxa[temp2 <= bin_edges[ k + 1]]
            inx=ind[np.logical_and(alts>bin_edges[k], alts<=bin_edges[k+1])]
            indices.append(inx)
        #store the number of bins averaged in each bin of the new alt vector    
        for i in range(len(indices)):
            D.n_ave_altitude[i] = len(indices[i])
               
    # iterate through all the items in a dictionary

        for name in vars(D).keys():
            value=getattr(D,name)#vars(D)[name]
            if name == 'var_mol':  # check if photon count stats requested
                #vars(R)['n_ave_altitude'] = np.zeros( (value.shape[0],k+1 ) )
                if hasattr(D,'n_ave_altitude'):
                    nave=getattr(D,'n_ave_altitude')
                    print 'HAD OLD AVE'
                else:
                    nave=np.ones((value.shape[0],ind.shape[0]))
                setattr(D,'n_ave_altitude',np.zeros((value.shape[0],n_alts)))
        # store number of points averaged in altitude
                for i in range(k):
                    D.n_ave_altitude[:,i] = np.sum(nave[:,indices[i]],1)#len(indices[i])
            if isinstance(value, Z_Array):

          # if this is the altitude vector replace with new alts

                if name == D._altitudevarname:
                    setattr(D,name,Z_Array(new_alts,summode=value.summode,dtype=value.dtype))# vars(D)[name] = Z_Array(new_alts)
                elif isinstance(value, TZ_Array):

                    # if time-altitude array

                    # do for each profile

                    shape = list(value.shape)
                    shape[1]=new_alts.shape[0]
                    temp = type(value)(np.zeros(shape),summode=value.summode,dtype=value.dtype)
                    #k = 0

                    # do the alt means
                    def myslice(alts):
                        ret=[slice(None)]
                        ret.append(alts)
                        while len(ret)<len(shape):
                            ret.append(slice(None))
                        return tuple(ret)

                    for k in range(new_alts.shape[0]):#while k <= new_alts.shape[0] - 1:
                        temp[myslice(k)] = nanmean(value[myslice(indices[k])], 1)
                        #k = k + 1
                    setattr(D,name,temp)
                elif isinstance(value, T_Array):
                    # if time array just pass through
                    pass
                    #vars(D)[name] = copy.deepcopy(vars(D)[name])
                elif isinstance(value, Z_Array):
                    temp = type(value)(np.zeros( new_alts.shape[0]),summode=value.summode,dtype=value.dtype)
                    for k in range(new_alts.shape[0]):#while k <= new_alts.shape[0] - 2:
                        temp[k] = nanmean(value[indices[k]])
                        #k = k + 1
                    setattr(D,name,temp)#vars(D)[name] = temp

    def meanTimes( D, new_times):
        n_times=new_times.shape[0]
        bin_edges = np.array([ None for x in range(0,n_times+1) ])
        bin_edges[0] = new_times[0] - timedelta(seconds=(new_times[1]-new_times[0]).total_seconds()/2.0)
        bin_edges[-1]= new_times[-1] + timedelta(seconds=(new_times[-1]-new_times[-2]).total_seconds()/2.0)
        bin_edges[1:n_times] = [ ( new_times[x-1] + timedelta(seconds=(new_times[x]-new_times[x-1]).total_seconds() / 2.0) ) for x in range(1,n_times)]
        D.resampleTimesByBins( bin_edges, new_times)

    def hereGoneBinTimes( D,interval_moments, new_times=None,*args,**kwargs):
        return D.resampleTimesByBins( interval_moments, interval_moments[0:-1] if new_times is None else new_times,*args,**kwargs)

    def includingRemainder( D, remainder ):
        if remainder is None:
            return D,None
        r=copy.deepcopy(remainder)
        r.append(D)
        return r,None

    def resampleTimesByBins( D, bin_edges, new_times,allow_nans=False,remainder=None,withRemainder=False):
        """
        create altitude means of time_z object centered on new time vector
        D= time_z_object
        new_times = time vector, center nan means of D on these timess
        returns smaller array with average values that properly treat NaN's
        It also checks for D.nave and multiplies entries according to the
        number of alt bins averaged together at that alititude"""
        if remainder is not None:
            D.prepend(remainder)
        if withRemainder:
            rem=copy.deepcopy(D)
        else:
            rem=None

    # create time vector centered on supplied new_new vector
        ntimes=getattr(D,D._timevarname).shape[0]
        ind = np.arange(getattr(D,D._timevarname).shape[0])

    # make a list of indices to indicate which elements to average for each entry
    # in new times vector

        temp = getattr(D,D._timevarname)[ :]
        indices = []
        keepindices = []
        for k in np.arange(bin_edges.shape[0]-1) :
            #inxa = ind[temp >= bin_edges[ k]]
            #temp2 = getattr(D,D._timevarname)[ inxa]
            #inx = inxa[temp2 < bin_edges[ k + 1]]
            inx=ind[np.logical_and(temp>=bin_edges[k], temp<bin_edges[k+1])]
            indices.append(inx.copy())
            if allow_nans or inx.size>0:
                keepindices.append(k)
        if bin_edges.shape[0]==0:
            remainders=ind.copy()
        else:
            remainders=np.array(ind[temp>=bin_edges[-1]],dtype='int')
        keepindices=np.array(keepindices,dtype='int')

    # iterate through all the items in a dictionary

        for name in vars(D).keys():
            value=getattr(D,name)#vars(D)[name]
            if name == 'var_mol':  # check if photon count stats requested
                if not hasattr(D,'n_ave'):
                    nave=np.ones_like(ind)
                else:
                    nave=getattr(D,'n_ave')
                setattr(D,'n_ave',np.zeros(len(indices)))
                if withRemainder:
                    setattr(rem,'n_ave',np.zeros(1))

        # store number of points averaged in time

                for i in range(len(indices)):  # only iterate thru actual valid indices
                    D.n_ave[i] = np.sum(nave[indices[i]])#len(indices[i])
                if len(remainders)>0:
                    rem.n_ave[0]=np.sum(nave[remainders])
                if not allow_nans:
                    if keepindices.size==0:
                        D.n_ave=np.zeros(0)
                    else:
                        D.n_ave=D.n_ave[keepindices]
            if isinstance(value,Z_Array):
                # if this is the time vector replace with new times
                if name == D._timevarname:
                    if withRemainder:
                        setattr(rem,name,getattr(D,name)[remainders])
                    setattr(D,name,T_Array(new_times[keepindices],summode=value.summode,dtype=value.dtype))
                    
                # if this the altitude vector just pass through    
                elif isinstance(value, T_Array):#includes TZ
                        # if time-altitude array
                        # do for each profile

                        shape=list(value.shape)
                        shape[0]=new_times.shape[0]
                        arr_cls=type(value)
                        temp = arr_cls(np.zeros(shape),summode=value.summode,dtype=value.dtype)
                        def myslice(*args):
                            ret=list(args)
                            while len(ret)<len(shape):
                                ret.append(slice(None))
                            return tuple(ret)

                        for k in range(shape[0]):
                            if len(indices[k])==0:
                                if value.sumOr or value.sumAnd:
                                    temp[myslice(k)]=0
                                else:
                                    temp[myslice(k)]=np.NaN
                            elif len(indices[k])==1:
                                temp[myslice(k)]=value[myslice(indices[k][0])]
                            else:
                                temp[myslice(k)] = value[myslice(indices[k])].doSum(axis=0)
                        if withRemainder:
                            setattr(rem,name,getattr(D,name)[myslice(remainders)])
                        if not allow_nans:
                            #print name,keepindices.size,value.shape
                            if keepindices.size==0:
                                shape[0]=0
                                temp=arr_cls(np.zeros(shape),summode=value.summode,dtype=value.dtype)
                            else:
                                temp=temp[myslice(keepindices)]
                        setattr(D,name,temp)#vars(D)[name] = temp
        if allow_nans:
            x=Time_Z_Group(like=D)
            if D.times.size>0:
                x.times=new_times[new_times<D.times[0]]
                if x.times.size>0:
                    D.prepend(x)
            else:
                if rem.times.size>0:
                    x.times=new_times[new_times<rem.times[0]]
                else:
                    x.times=newtimes.copy()
                if x.times.size>0:
                    rem.prepend(x)

        return rem
 
    def prepend(self, other):
        othercopy=copy.deepcopy(other)
        othercopy.append(self)
        for name in vars(othercopy).keys():
            setattr(self,name,getattr(othercopy,name))

    def append(self, D,_seenvars=None,_trace=None):
        """Append Time_Z_Data D to self
        """
   

        seenvars=_seenvars
        if seenvars==None:
            seenvars={}
        trace=_trace
        if trace==None:
            trace=[]
        # if D doesn't have a times object, or the times object is 0 length, 
        # do not append
        if not isinstance(D,Time_Z_Group) or self._timevarname!=D._timevarname:
            warnings.warn("compatibility of append with "+('.'.join(_trace)))
            if D is None:
                warnings.warn("trying to append a None object to a Time_Z_Group FIXME")
                return
            if True and type(D)==list:
                if len(D)==0:
                    return#FIXME this is happening because of frame reassembly
                if len(D)==1 and isinstance(D[0],Time_Z_Group):
                    print 'FIXME array append bug...'
                    self.append(D[0])
                    return
            print type(self),self
            print type(D),D
            try:
                raise RuntimeError("MISPLACED NONE")
            except:
                import traceback
                traceback.print_exc()
                raise
        if not self._can_append:
            for k in vars(self).keys():
                if not k.startswith('_') and k not in vars(D).keys():
                    delattr(self,k)
            for k,v in vars(D).items():
                setattr(self,k,copy.deepcopy(v))
            return
        DEmpty = True
        dshots = 0
        try:
            dshots=getattr(D,D._timevarname).shape[0]
            DEmpty = (dshots == 0)
        except AttributeError:
            pass
        
        if DEmpty:
            LOG.debug( 'append() : NOTE: appending an empty Time_Z_Group' )
            #return

        selfEmpty = True
        nshots = 0

        # self.times.shape produces an error when self is empty,
        # this allows appending onto an empty array.
        try:
            nshots=getattr(self,self._timevarname).shape[0]
            selfEmpty = (nshots == 0)
        except AttributeError:
            pass 
        
        #print 'Appending:%s + %s' % (self,D)
       
        #make a list of all pre existing variables
        #used if a variable not in structure appears and
        #it will be necessary to back fill previous times with NaN's
        pre_existing_variables = vars(self).keys()
        
        #store length of pre_existing_variables in case we need it.

        
        #print 'pre_existing_variables'
        #print pre_existing_variables
        #print ' '
      

        verifyNew(self,seenvars,['source']+trace)
        verifyNew(D,seenvars,['newdata']+trace)

    
        # "normal" case - append all items in D to self 
        # iterate through all the items in a dictionary
        Dvars=vars(D).keys()
        for name in pre_existing_variables:

            if name in Dvars:#for any variable that did exist but doesn't in the new record, extend if its a T* array
                continue
            value = getattr(self,name)
            if isinstance(value, T_Array):
                #print 'extending: %s T*_Array ' % name
                try:
                    value.extend(dshots)
                except ValueError:
                    setattr(self,name,copy.deepcopy(value))
                    getattr(self,name).extend(dshots)

#               if len(vars(self)[name].shape) != 1:
#                    print 'append - time variable %s is >1D' % name


        for name in Dvars:
            #print 'append---',name
            value = getattr(D,name)
            if name in pre_existing_variables:
                mval=getattr(self,name)
            else:
                mval=None
            #print name
            if selfEmpty and (name not in pre_existing_variables):#didn't exist before, and we're empty anyway
                #print 'append adding entry',name
                try:
                    setattr(self,name,copy.deepcopy(value))
                except:
                    raise RuntimeError('Exception trying to deepcopy '+name+' of type '+str(type(value)))
            elif isinstance(value,Time_Z_Group) or isinstance(mval,Time_Z_Group):#safeIsInstance(value, Time_Z_Group) or safeIsInstance(mval, Time_Z_Group):
                if isinstance(mval,Time_Z_Group):#safeIsInstance(mval, Time_Z_Group):
                    LOG.debug( '******** GROUP append '+name )
                    try:
                        if value is not None:
                            mval.append(value,seenvars,trace+[name])
                    except RuntimeError:
                        print name,"won't append properly"
                        raise
                    LOG.debug( '******** done group append '+name )
                else: #if safeIsInstance(vars(D)[name], Time_Z_Group):
                    setattr(self,name,copy.deepcopy(value))                    
                    LOG.debug( 'append replacing entry '+name+' for destination not being a time_z_group while source is' )
            elif isinstance(value, T_Array) or isinstance(mval, T_Array): #includes TZ
                if name in pre_existing_variables and nshots>0:
                    if not isinstance(value, T_Array): #includes TZ
                        LOG.warning( 'WARNING '+name+' in append is not a T*_Array' )
                    try:
                        #print name
                        verifyNew(mval,seenvars,['source']+trace+[name])
                        verifyNew(value,seenvars,['newdata']+trace+[name])
                        try:
                            mval.append(value)
                        except RuntimeError:
                            print 'ERROR ON '+'.'.join(trace+[name])
                            raise
                    except AttributeError:
                        LOG.warning( 'Warning: array '+name+'was not a T array' )
                        mval=copy.deepcopy(T_Array(mval,summode=value.summode))
                        setattr(self,name,mval)
                        try:
                            mval.append(value)                                                
                        except RuntimeError:
                            print 'ERROR ON '+'.'.join(trace+[name])
                            raise
                    except ValueError:
                        #print 'Warning: array ',name,'must be copied to own values'
                        mval=copy.deepcopy(mval)
                        setattr(self,name,mval)
                        try:
                            mval.append(value)                        
                        except RuntimeError:
                            print 'ERROR ON '+'.'.join(trace+[name])
                            raise
                elif nshots>0:   #fill in missing previous values with NaN's
                    LOG.debug( 'hau.append:'+ name+ '--did not exist in structure*******************' )
                    #print nan_array.shape, '--size of NaN array'
                    #print vars(D)[name].shape ,'--size of new varaible'
                    sh=list(value.shape)
                    sh[0]=nshots
                    tmp = np.nan*np.ones(sh) 
                    tmp = type(value)(tmp,summode=value.summode,dtype=value.dtype)
                    tmp = copy.deepcopy(tmp)

                    tmp.append(value)#vars(D)[name])
                    #tmp = T_Array(np.hstack((nan_array, vars(D)[name])))
                    
                    
                    
                    #vars(self)[name] = tmp
                    setattr(self,name,tmp)
                else:
                    #vars(self)[name]=copy.deepcopy(vars(D)[name])
                    setattr(self,name,copy.deepcopy(value))
#               if len(vars(self)[name].shape) != 1:
#                    print 'append - time variable %s is >1D' % name
            elif isinstance(value, Z_Array) or isinstance(mval, Z_Array):#includes T and TZ, but omitted because was tested earlier
                verifyNew(mval,seenvars,['source']+trace+[name])
                verifyNew(value,seenvars,['newdata']+trace+[name])
                if not isinstance(value, Z_Array):
                    pass
                elif not isinstance(mval, Z_Array) or (anynan(mval) and not anynan(value)):
                    setattr(self,name,copy.deepcopy(value))
            elif isinstance(value,datetime):
                if mval==None:
                    setattr(self,name,copy.deepcopy(value))#keep the first
            elif isinstance(value,timedelta):
                if mval==None:
                    setattr(self,name,copy.deepcopy(value))
                else:
                    setattr(self,name,timedelta(seconds=mval.total_seconds()+value.total_seconds()))
            else:
                # don't update start time on append
                if not name == 'start_time':
                    # ordinary array - just copy it
                    #print type(vars(D)[name])
                    try:
                        #vars(self)[name] = copy.deepcopy(vars(D)[name])
                        setattr(self,name,copy.deepcopy(value))
                    except:
                        raise RuntimeError('Exception trying to deepcopy '+name+' of type '+str(type(value)))

    def check_tarray(self, name):
        for (n, value) in vars(self).items():
            #print 'check_tarray() : name = ', n, ' value = ', value, 'type(value) ', type(value)
            if n == name:
                #print 'check_tarray() MATCH! : name = ', name, ' value = ', value, 'type(value) ', type(value)
                return safeIsInstance(value, T_Array)

        return False

    def check_tarray2(self, D, name):
        for (n, value) in vars(D).items():
            #print 'check_tarray() : name = ', n, ' value = ', value, 'type(value) ', type(value)
            if n == name:
                print 'check_tarray() MATCH! : name = ', name, ' value = ', value, 'type(value) ', type(value)
                return safeIsInstance(value, T_Array)

        return False


def selectByTime(D, time,offset=0):
    """select object members nearest to a given time
       (formerly tzu.select_time_z_data_object)
       time_z objects include 2-d numpy arrays that must include D.times
       with times in python datetimes. T_Array entries 
       correspond to times and Z_Array entries to altitude.
       Elements that are not of type numpy.ndarray are copied to output 
       without change. A field bounds is added to the object to indicate
       if requested time is between the first and last entry in D.times.
       bounds is -1 before the first entry, 0 within time period and +1
       after last entry. A field expire_time is added to indicate a point
       10% beyound the half-way point between this entry and the next entry"""
    #print 'selectByTime called!!'
 
    dt = T_Array([ (n - time).total_seconds() for n in getattr(D,D._timevarname) ])
    inx=-1
    intv=None
    for i,t in enumerate(getattr(D,D._timevarname)):
        if t<time and (intv is None or intv>(time-t)):
            inx=i
            intv=time-t
    inx+=offset
    ntimes = getattr(D,D._timevarname).shape[0]

    if inx<0 or inx>=ntimes:
        return None

    # iterate through all the items in a dictionary
    R = Time_Z_Group(like=D)
    for name in vars(D).keys():
        value=getattr(D,name)#vars(D)[name]
        #print 'name = ', name, ' value = ', value, ' type = ', type(value)
        #print 'safeIsInstance(value, T_Array) = ', safeIsInstance(value, T_Array)
        #print 'safeIsInstance(type(value), T_Array) = ', safeIsInstance(type(value), T_Array)
        if isinstance(value, TZ_Array):
            # returning a single Z element from a TZ_Array
            #vars(R)[name] = copy.deepcopy(Z_Array(value[inx]))
            setattr(R,name, copy.deepcopy(Z_Array(value[inx],summode=value.summode,dtype=value.dtype)))
            #print 'created R[%s] , R[%s].shape =' % (name, vars(R)[name].shape)
        elif isinstance(value, T_Array):
            #print 'found a T_Array'
            # return a single element as a new T_Array
            setattr(R,name,copy.deepcopy(T_Array([value[inx],],summode=value.summode,dtype=value.dtype)))#vars(R)[name] = copy.deepcopy(T_Array([value[inx],]))
            #print 't_array =', name, inx,value[inx]
        else:
            #vars(R)[name] = copy.deepcopy(value)
            setattr(R,name,copy.deepcopy(value))
    #R.bounds_flag = bounds
    #R.expire_time = np.array(expire_time)
    return R


def trimTimeInterval( D,
    start_time,
    end_time,
    ):
    """trim all the fields in a dictionary to [start_time>=getattr(D,D._timevarname)<=end_time]"""

    # find start_index and end_index
    R = Time_Z_Group(like=D)
    R.append(D)
    R.trimTimeInterval(start_time,end_time)
    return R


def check_tarray3( D, name ):
    try:
        for (n, value) in vars(D).items():
            #print 'check_tarray3() : name = ', n, ' value = ', value, 'type(value) ', type(value)
            if n == name:
                print 'check_tarray3() MATCH! : name = ', name, ' value = ', value, 'type(value) ', type(value)
                return safeIsInstance(value, T_Array)
    except TypeError,emsg:  # catch passing a non-collection/dictionary
        print 'check_tarray3 :', emsg
        pass

    return False
def test():
    r1 = np.arange(10)
    t1 = T_Array(r1)
    t1 = t1 * r1
    print 'r1 =', t1, 'type = ', type(r1)
    print 't1 =', t1, 'type = ', type(t1)
    dc1 = DataCollection()
    dc1.foo = 'bar'
    dc2 = DataCollection()
    dc2.type = 'wolf'
    dc1.mergeTimeVectors(dc2)
    print dc1


if __name__ == '__main__':
    test()
