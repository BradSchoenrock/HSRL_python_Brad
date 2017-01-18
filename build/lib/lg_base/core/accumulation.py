import numpy as np
try:
    from bottleneck import nanmean,nansum
except ImportError:
    print
    print 'No bottleneck.nanmean available! Falling back to SLOW scipy.stats.nanmean'
    print
    from scipy.stats import nanmean
    from numpy import nansum

def maybetuple(x):
  if len(x)==1:
    return x[0]
  return tuple(x)

from lg_base.core.array_utils import keepnansum

def reshapemode(v,m,s):
  if m==0:
    return v
  if m==1:
    return v.reshape(s)
  tmp=np.ones(s,dtype=np.float)
  tmpshape=list(v.shape)
  while len(tmpshape)<len(s):
    tmpshape.append(1)
  alllist=[]
  sublist=[]
  for x in range(len(tmpshape)):
    alllist.append(slice(None))
    sublist.append(slice(0,tmpshape[x]))
  tmp[maybetuple(alllist)]=np.NaN
  tmp[maybetuple(sublist)]=v
  return tmp

def matchShapes(v1,v2):
    theshape=list(v1.shape)
    willchangev1=0
    willchangev2=0
    while len(v2.shape)>len(theshape):
      theshape.append(1)
      willchangev1=1
    if len(theshape)!=len(v2.shape):
      willchangev2=1
    for i,x in enumerate(v2.shape):
      if theshape[i]<x:
        theshape[i]=x
        willchangev1=2#change size
      elif theshape[i]>x:
        willchangev2=2
    return reshapemode(v1,willchangev1,theshape),reshapemode(v2,willchangev2,theshape)



def accumulate(prof,oldprof,src,indicies,field,divval=1.0,mask=None,pref='',filler=None,extravars=[]):
  """ Accumulate variables, with oldprof as the old, prof as the new

  :param extravars: if the field doesn't exist in src, or indicies is empty, will copy the variable, and those listed here from oldprof
  """
      
  if src is None or not hasattr(src,field):
    if oldprof is not None and hasattr(oldprof,'sum_'+pref+field):
      setattr(prof,'sum_'+pref+field,getattr(oldprof,'sum_'+pref+field))
      setattr(prof,pref+field,getattr(oldprof,pref+field))
      for var in extravars:
        if hasattr(oldprof,var):
          setattr(prof,var,getattr(oldprof,var))
    elif filler is not None:
      setattr(prof,'sum_'+pref+field,filler)
      setattr(prof,pref+field,filler)      
    return False
  vp=getattr(src,field)
  klass=type(vp)
  if len(vp.shape)==2:
    v=vp[indicies,:].copy()
  elif len(vp.shape)==1:
    v=vp[indicies].copy()
  else:
    raise RuntimeError('Bad shape')
  sh=list(v.shape)
  if mask is not None:
    if mask.size>0 and indicies.size>0:
      v*=mask[indicies,:]
    else:
      v*=np.NaN
  if False:#this makes it much closer to a nanmean in the end
    tmp=np.ones_like(v,dtype=np.float64)
    tmp2=np.ones_like(v,dtype=np.float64)
    tmp[np.isnan(v)]=0.0
    tmp=np.sum(tmp,0)
    tmp2=np.sum(tmp2,0)
    v=keepnansum(v,0)
    v*=tmp2/tmp
  else:
    v=keepnansum(v,0)
  if not hasattr(v,'shape'):
    v=np.array(v)
  if len(sh)==(len(v.shape)+1):
    #print sh,v.shape
    sh=[1]+list(v.shape)
    v=v.reshape(sh)
  if oldprof is not None and hasattr(oldprof,'sum_'+pref+field):
    ov=getattr(oldprof,'sum_'+pref+field)
    v,ov=matchShapes(v,ov)
    if False:
      v+=ov
    else:
      oldfinitevals=np.isfinite(ov)
      finitevals=np.isfinite(v)
      bothfinite=np.logical_and(oldfinitevals,finitevals)
      v[bothfinite]+=ov[bothfinite]
      nonfinite=np.logical_not(finitevals)
      if nonfinite.any():
        v[nonfinite]=ov[nonfinite]# if ov is also nonfinite, don't care
  if not isinstance(v,klass):
    v=klass(v)
  setattr(prof,'sum_'+pref+field,v)
  setattr(prof,pref+field,v/divval)
 
  return True

