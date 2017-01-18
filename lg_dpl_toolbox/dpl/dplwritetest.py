#!/usr/bin/python
# -*- coding: utf-8 -*-

if __name__ == '__main__':
  
  from dpl_rti import dpl_rti
  import datetime
  import dpl_create_templatenetcdf as dpl_ctnc
  from netCDF4 import Dataset
  from lg_base.core.locate_file import locate_file
  
  
  v=None
  if 0:
    n=Dataset('out.nc','a',clobber=False)
    for i in n.variables:
        print i
        n.variables[i].dpl_py_binding='dne'
    n.close()


  if 1:
    n=Dataset('outtest2.nc','w',clobber=True)
    instrument='bagohsrl'
    timeave=10
    rangeave=30
    n.instrument=instrument
    tm=datetime.datetime.utcnow()
    tm=datetime.datetime(2012,8,8,19,0,0)#tm-datetime.timedelta(seconds=.01*60*60)
    tmp=tm-datetime.timedelta(seconds=10*60*60)
    #tm=None
    retrievalslice=6*60*60
    gen=dpl_rti(instrument,tmp,tm,
                #datetime.datetime(2012,6,1,12,0,0),datetime.datetime(2012,6,1,17,0,0),
                datetime.timedelta(seconds=timeave),datetime.timedelta(seconds=retrievalslice),
                0,12000,rangeave)

    v=None

    donetest=5

    for i in gen:
        if v==None:
            v=dpl_ctnc.dpl_create_templatenetcdf(locate_file('hsrl_nomenclature.cdl'),n,i)
        v.appendtemplatedata(i)
        if i.rs_inv.times.shape[0]==0:
            donetest-=1
        else:
            donetest=5
        if donetest==0:
            break
        #raise TypeError
        #break
        n.sync()
    n.close()
