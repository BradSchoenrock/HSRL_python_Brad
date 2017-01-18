#!/usr/bin/python
# -*- coding: utf-8 -*-

import datetime
from netCDF4 import Dataset
import dpl_create_templatenetcdf as dpl_ctnc
from matplotlib.dates import date2num
import numpy
import dplkit.role.narrator
import lg_base.core.array_utils as hau
import dplkit.role.decorator

@dplkit.role.decorator.exposes_attrs_of_field('attributesSource')
@dplkit.role.decorator.autoprovidenested(nestedclasses=[hau.Time_Z_Group,hau.rs_xfer,dict])
class dpl_read_templatenetcdf(dplkit.role.narrator.aNarrator):
    def __init__(self, dplfile, group=None, getAttributes=None, *args, **kwargs):
        super(dpl_read_templatenetcdf,self).__init__(dplfile)
        self.attributesSource=getAttributes
        if len(args):
            print 'Unused args = ',args
        if len(kwargs):
            print "Unused kwargs = ",kwargs
        self.filename=dplfile
        self.group=group
        self.data=Dataset(dplfile,'r')
        (self.xlateTable,notused)=dpl_ctnc.generateDPLNCTemplateTable(self.data,group=group,inclusive=True)
        self.val=None

    def raw_netcdf(self):
        return self.data

    def __repr__(self):
        return 'DPL Read TemplateNetCDF Narrator for file "%s"' % (self.filename)

    def __iter__(self):
        return self()
    
    def read(self, *args, **kwargs):
        if self.val!=None:
            yield self.val
            return
        if len(args):
            print 'Unused args = ',args
        if len(kwargs):
            print "Unused kwargs = ",kwargs
        ret=hau.Time_Z_Group()
        delattr(ret,'times')

        data=self.data
        if self.group is not None:
            for g in self.group.split('/'):
                data=data.groups[g]

        if 'base_time' in data.variables:
            epochtime=datetime.datetime(1970,1,1,0,0,0)
            basetime=epochtime+datetime.timedelta(seconds=long(data.variables['base_time'].getValue()))
            print 'basetime is ',basetime
        if 'time_coverage_start' in data.variables and 'time_coverage_end' in data.variables:
            basetime=datetime.datetime.strptime(''.join(data.variables['time_coverage_start'][:20]), '%Y-%m-%dT%H:%M:%SZ')
        for v in self.xlateTable:
            vari=data.variables[v]
            parts=self.xlateTable[v].split('.')
            finalpart=parts[-1]
            parts=parts[0:(len(parts)-1)]
            base=ret
            for i in parts:
                if not hasattr(base,i):
                    setattr(base,i,hau.Time_Z_Group())
                    delattr(getattr(base,i),'times')
                base=getattr(base,i)
            idx=[slice(None) for x in range(len(vari.shape))]
            if len(idx)==0:
                newval=vari.getValue()
            elif len(idx)==1:
                newval=numpy.array(vari[idx[0]])
            else:
                newval=numpy.array(vari[tuple(idx)])
            if 'dpl_py_type' in vari.ncattrs():
                dpltype=vari.dpl_py_type[:]
                if dpltype=='matplotlib_num2date':
                    newval=numpy.array([(date2num(basetime+datetime.timedelta(seconds=float(d))) if d<1e35 else float('nan')) for d in newval])
                if dpltype=='python_datetime':
                    newval=numpy.array([(basetime+datetime.timedelta(seconds=float(d)) if d<1e35 else None) for d in newval])
            setattr(base,finalpart,newval)
        addpartsbylength=dict()
        fillins=(('times',),
                 ('delta_t',),
                 ('msl_altitudes','heights'))
        for ka,f in vars(ret).items():
            for ks in fillins:
                for k in ks:
                    if hasattr(f,k):
                        t=getattr(f,k)
                        if t.size==0:
                            continue
                        for nk in ks:
                            if nk not in addpartsbylength:
                                addpartsbylength[nk]=dict()
                            if t.size not in addpartsbylength[nk]:
                                addpartsbylength[nk][t.size]=t
        #print 'available parts:',addpartsbylength
        for k,f in vars(ret).items():
            dim_preferredsizes=[[],[]]
            try:
                for ak,av in vars(f).items():
                    if hasattr(av,'shape') and len(av.shape)>=2:
                        for x in range(2):
                            dim_preferredsizes[x].append(av.shape[x])
            except TypeError:
                #print 'no parts for ',k
                continue
            #print 'sizes for',k,dim_preferredsizes
            if len(dim_preferredsizes[0])==0 or len(dim_preferredsizes[1])==0:
                continue
            dim_preferred=[]
            for x in range(2):
                dim_preferred.append(int(numpy.median(dim_preferredsizes[x])))            
            #print 'sizes are',dim_preferred,'median of ',dim_preferredsizes
            for nk,kv in addpartsbylength.items(): #TOTAL HACK
                x=None if not hasattr(f,nk) else getattr(f,nk)
                if x is None or not hasattr(x,'shape') or x.size==0:
                    try:
                        if nk in ('times','delta_t'):
                            dimpreferredidx=0
                        else:
                            dimpreferredidx=1
                        count=dim_preferred[dimpreferredidx]
                        if count in kv:
                            setattr(f,nk,kv[count])
                        #else:
                        #    print nk,'doesnt exist for length',count
                    except AttributeError:
                        pass
        self.val=ret
        yield ret #this means this operates as an iterator that runs once
