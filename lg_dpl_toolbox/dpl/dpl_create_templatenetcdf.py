#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
from netCDF4 import Dataset, stringtoarr, chartostring
import math
import datetime
from matplotlib.dates import  num2date
import tempfile
from lg_dpl_toolbox.core.CDL import CDL2dict
import numpy
from collections import OrderedDict
import lg_base.core.git_tools as git_tools

STRING_LENGTH_SHORT = 32

""" Explanation of CDL Templates, as relevant to this file:
All file attributes are copied verbatim
Variables are only created if a dpl_py_binding attribute matches a field in the framestream
-- discovery is done using the provides if provides has size information
-- otherwise, this is done using the first frame
Dimensions are included as neccessary to the variables it includes.  Size of the variables are determined as follows:
-- if a non-zero number is given, it is statically that size
-- if 0 is given, it will adapt to the size to fit the frame
-- UNLIMITED results in an unlimited dimension. if the file format doesn't support more than one, the template should avoid using more than one
--- special case: withUnlimited can be set to a static number for an instance of this object to replace all instances of unlimited with a fixed dimension (good for NC3)
---- to this end, remember any stream should have only one coincident unlimited dimension, and a file may have multiple streams feeding into it.

typing:
all types are naively dropped into the variable, with a few exceptions dependant on explicit known types, and NetCDF format
a type described in dpl_py_type can be python_datetime, so it will be mapped from a python datetime object to a double for seconds since basetime.
basetime is determined by the format:
for CFRadial, it is a string representation in time_coverage_start
for other files, it is base_time as seconds from unix epoch
"""

def doKeep(val,keeparr,remarr):
    if (keeparr is not None) and val not in keeparr:
        return False
    if (keeparr is None or val not in keeparr) and remarr is not None and val in remarr:
        return False
    return True

def copyAttrs(inputd,outputobj,keepatts=None,removeatts=None,vtype=None):
    if not isinstance(inputd['attributes'],OrderedDict):
        raise RuntimeError("Attributes aren't ordered")
    for att in inputd['attributes']:
        if not doKeep(att,keepatts,removeatts):
            #print 'dropping attr %s' % att
            continue
        #print 'copying attr',att,':',type(inputd['attributes'][att]),inputd['attributes'][att]
        if att in ('_FillValue','missing_value','range','insufficient_data') and vtype is not None:
            outputobj.setncattr(att,inputd['attributes'][att].astype(vtype))            
        else:
            outputobj.setncattr(att,inputd['attributes'][att])


typecodes={
    'byte':'b',
    'char':'c',
    'short':'i2',
    'int':'i4',
    'long':'i8',
    'double':'d',
    'float':'f'
}

def filterNC(inputd,outputnc,keepvars=None,keepdims=None,keepatts=None,removevars=None,removedims=None,removeatts=None,dimsizes=None):
    #attributes first
    copyAttrs(inputd,outputnc,keepatts,removeatts)

    #then dims
    if not isinstance(inputd['dimensions'],OrderedDict):
        raise RuntimeError("Dimensions aren't ordered")
    if not isinstance(inputd['variables'],OrderedDict):
        raise RuntimeError("Variables aren't ordered")
    for dim in inputd['dimensions']:
        if not doKeep(dim,keepdims,removedims):
            #print 'dim %s dropped' % dim
            continue
        if dim in outputnc.dimensions:
            print 'WARNING: dimension '+dim+' already in output, likely from layared netcdf artists. skipping initialization'
            continue
        l=inputd['dimensions'][dim]
        if inputd['dimensions'][dim]==-1:
            l=None
        if dimsizes and dim in dimsizes:
            l=dimsizes[dim]
        if l==0:#this should be None if it is supposed to be unlimited
            raise RuntimeError('Unexpected unlimited dimension crept in. Name:'+dim)
        outputnc.createDimension(dim,l)

    #then vars
    for var in inputd['variables']:
        if not doKeep(var,keepvars,removevars):
            #print 'var %s dropped' % var
            continue
        if var in outputnc.variables:
            print 'WARNING: variable '+var+' already in output, likely from layared netcdf artists. skipping initialization'
            continue
        #print inputnc.variables[var]
        typeval=typecodes[inputd['variables'][var]['type']]
        fillv=None if '_FillValue' not in inputd['variables'][var] else inputd['variables'][var]['_FillValue']
        if fillv==None and typeval=='c':
            fillv=u'\x00'
        outputnc.createVariable(var,typeval,inputd['variables'][var]['dimensions'],fill_value=fillv)
        copyAttrs(inputd['variables'][var],outputnc.variables[var],vtype=typecodes[inputd['variables'][var]['type']])

def getsubstruct(obj,path):
    bndclass=path.split('.')
    tmp=obj
    bndpath=''
    for ss in bndclass:
        if isinstance(tmp,dict) and ss in tmp:
            tmp=tmp[ss]
        elif hasattr(tmp,ss):
            tmp=getattr(tmp,ss)
        else:
            #print 'no %s in structure (path = %s)' % (ss,bndpath)
            return None
        if len(bndpath)==0:
            bndpath=ss
        else:
            bndpath='.'.join((bndpath,ss))
    return tmp

def generateDPLNCTemplateTable(ncdata,inclusive=False,group=None):
    ret=OrderedDict()
    dims=OrderedDict()
    if isinstance(ncdata,dict):
        if group is not None:
            raise RuntimeError("Grouping isn't supported from dictionary yet")
        for n,var in ncdata['variables'].items():
            if 'dpl_py_binding' in var['attributes']:
                ret[n]=var['attributes']['dpl_py_binding']
            elif 'dpl_py_require' in var['attributes']:
                ret[n]=None
            if n in ret:
                for d in var['dimensions']:
                    if d not in dims:
                        dims[d]=ncdata['dimensions'][d]
                        if dims[d]<0:
                            dims[d]=None
    else:
        if group is not None:
            for g in group.split('/'):
                ncdata=ncdata.groups[g]
        for i in ncdata.variables:
            if hasattr(ncdata.variables[i],'dpl_py_binding'):
                ret[i]=ncdata.variables[i].dpl_py_binding
            elif hasattr(ncdata.variables[i],'dpl_py_require'):
                ret[i]=None
            if i in ret:
                for d in ncdata.variables[i].dimensions:
                    if d not in dims:
                        dims[d]=len(ncdata.dimensions[d])
                        if ncdata.dimensions[d].isunlimited():
                            dims[d]=None
    return (ret,dims)

def fitDPLNCTemplateTable(table,n,templatestruct,bindings_to_keep=None,bindings_to_drop=None):
    ret=OrderedDict()
    dimsizes=OrderedDict()
    varorder=[]
    singledimvar=[]
    othervar=[]
    for f in table:
        bnd=table[f]
        if bnd is None:
            othervar.append(f)
            continue
        if bnd=='dne' or not doKeep(bnd,bindings_to_keep,bindings_to_drop):
                continue
        tmp=getsubstruct(templatestruct,bnd)
        if tmp is not None:
            vardims=n['variables'][f]['dimensions']
            if len(vardims)>1:
                othervar.append(f)
            else:
                singledimvar.append(f)
    varorder.extend(singledimvar)
    varorder.extend(othervar)
    for f in varorder:
        bnd=table[f]
        if bnd is None:
            ret[f]=None
            for d in n['variables'][f]['dimensions']:
                if d in dimsizes:
                    pass
                elif n['dimensions'][d]>0:
                    dimsizes[d]=n['dimensions'][d]
                else:
                    del ret[f]
                    break
            continue
        if bnd=='dne' or not doKeep(bnd,bindings_to_keep,bindings_to_drop):
            continue
        tmp=getsubstruct(templatestruct,bnd)
        if tmp is not None:
            ret[f]=bnd
            #dimid=0
            if isinstance(tmp,dict) and 'shape' in tmp:
                tmpsh=[ x for x in tmp['shape'] ]
            elif hasattr(tmp,'shape'):
                tmpsh=tmp.shape
            else:
                try:
                    tmpsh=[len(tmp)]
                except:
                    tmpsh=[]
            vardims=n['variables'][f]['dimensions']
            if len(tmpsh)>len(vardims) and tmpsh[0]==1:
                tmpsh=tmpsh[1:]
            while len(tmpsh)<len(vardims):
                if len(tmpsh)>0:
                    tmpsh=tuple([1]+list(tmpsh))
                else:
                    tmpsh=(1,)
            if len(tmpsh)!=len(vardims):
                #print 'size of variable %s is different than netcdf %s.' % (bnd,f)
                #print tmpsh
                #print vardims
                pass
            for dimid,d in enumerate(vardims):
                if d in dimsizes:
                    if dimsizes[d]!=tmpsh[dimid] and not n['dimensions'][d]==-1:
                        print 'dimension %i %s is different sizes in var %s (NC %s). %i vs %i in structure' % (dimid,d,bnd,f,dimsizes[d],tmpsh[dimid])
                else:
                    if n['dimensions'][d]==-1:
                        dimsizes[d]=None
                        #print 'dimension %i %s is unlimited %s (%s) actual size %i' % (dimid,d,bnd,f,tmpsh[dimid])
                    elif n['dimensions'][d]==0:
                        dimsizes[d]=tmpsh[dimid]
                        #print 'dimension %i %s is now %i from %s (%s)' % (dimid,d,tmpsh[dimid],bnd,f)
                    else:
                        dimsizes[d]=n['dimensions'][d]
                        print 'dimension %i %s is now %i from CDL template.' %(dimid,d,dimsizes[d])
                        if dimsizes[d]!=tmpsh[dimid]:
                            print "WARNING: variable %s doesn't match size" % (d)
    return (ret,dimsizes)

def string1D(var):
    ret=''
    for x in range(var.shape[0]):
        ret=ret+var[x]#chartostring(var[:].reshape([1]+list(var[:].shape)))[0]
    while len(ret)>0 and ret[-1]=='N':
        ret=ret[:-1]
    if len(ret)==0:
        return None
    return ret

class dpl_create_templatenetcdf(object):
    def __init__(self,templatename,outputfilename_ornc,firsttemplatesource,bindings_to_keep=None,withUnlimited=None,withBasetime=None,
        doAppend=False,addAttributes={},cfradial=None,group=None,forModule=None):
        self.templatename=templatename
        self.isCFRadial=cfradial
        self.date_fmt = '%Y-%m-%dT%H:%M:%SZ'
        self.basetime = None
        self.groupname=group
        from lg_base.core.locate_file import locate_file
        print 'templatename=',self.templatename
        if isinstance(self.templatename,(list,tuple)):
            cdl=None
            for t in self.templatename:
                tc=CDL2dict(locate_file(t,forModule=forModule))
                if cdl is None:
                    cdl=tc
                elif tc is not None:
                    for k in cdl.keys():
                        cdl[k].update(tc[k])
                if cdl is None or tc is None:
                    print 'locate failed',t,forModule
                    raise RuntimeError('Failed to get file '+t)
        else:
            cdl=CDL2dict(locate_file(templatename,forModule=forModule))
        if self.isCFRadial==None:
            self.isCFRadial='Conventions' in cdl['attributes'] and 'CF/Radial' in cdl['attributes']['Conventions']
        varstopurge=None
        dimstopurge=None
        varstokeep=[]
        dimstokeep=[]
        dimsizes=OrderedDict()
        #self.fn=outputfilename
        self.bindings=OrderedDict()
        self.appendDims=OrderedDict()
        self.willclose=True
        self.addAttributes=addAttributes

        (tbindings,tdims)=generateDPLNCTemplateTable(cdl)
        (self.bindings,dimsizes)=fitDPLNCTemplateTable(tbindings,cdl,firsttemplatesource,bindings_to_keep=bindings_to_keep)
        dimstokeep.extend([d for d in dimsizes])
        varstokeep.extend([v for v in self.bindings])

        if 'time' in dimsizes and dimsizes['time']==None:
            self.appendDims['time']=0
        if self.isCFRadial:
            if 'time' in self.bindings and 'time' in dimsizes:
                self.cfradial_writetime=True
                varstokeep.append('time_coverage_start')
                varstokeep.append('time_coverage_end')
            else:
                self.cfradial_writetime=False
            varstokeep.extend(['sweep_number', 'sweep_mode', 'fixed_angle','sweep_start_ray_index','sweep_end_ray_index'])
            dimstokeep.extend(['string_length_short','sweep'])
        #if withBasetime!=None:# 'time' in self.bindings and 'time' in dimsizes:
        varstokeep.append('base_time')
        if isinstance(outputfilename_ornc,basestring):
            self.dataset=Dataset(outputfilename_ornc,'w',clobber=not doAppend)
        else:
            self.dataset=outputfilename_ornc
            self.willclose=False
        self.rootdataset=self.dataset
        if self.groupname is not None:
            grouplist=self.groupname.split('/')
            for g in grouplist:
                if g in self.dataset.groups:
                    self.dataset=self.dataset.groups[g]
                else:
                    self.dataset=self.dataset.createGroup(g)
        for dim in dimstokeep:
            if not dim in dimsizes:#explicitly kept dimension
                pass
            elif dimsizes[dim]==None:
                if not dim in self.appendDims:
                    self.appendDims[dim]=0
                if withUnlimited!=None:
                    dimsizes[dim]=withUnlimited
        filterNC(cdl,self.dataset,keepvars=varstokeep,keepdims=dimstokeep,removevars=varstopurge,removedims=dimstopurge,dimsizes=dimsizes)
        if withBasetime is not None:
            #self.basetime=withBasetime
            assert(self.setFileBasetime(withBasetime))
        else:
            self.getBasetime()# self.basetime=None
        self.dataset.sync()
        if doAppend and withUnlimited==None:
            for adim in self.appendDims.keys():# start by loading record dimensions that should be appended
                self.appendDims[adim]=len(self.dataset.dimensions[adim])

    def close(self):
        if self.isCFRadial:
            if hasattr(self,'end_time') and self.cfradial_writetime:
                self.dataset.variables['time_coverage_end'][:] = stringtoarr(self.end_time.strftime(self.date_fmt), STRING_LENGTH_SHORT)
                self.dataset.variables['sweep_end_ray_index'][0] = self.dataset.variables['time'].shape[0] - 1
            self.dataset.variables['sweep_number'][0] = 0
            self.dataset.variables['sweep_mode'][0] = stringtoarr("pointing", STRING_LENGTH_SHORT)
            self.dataset.variables['fixed_angle'][0] = 0.0
            self.dataset.variables['sweep_start_ray_index'][0] = 0
            self.dataset.sync()


    def getBasetime(self,require=False):
        if self.basetime is None:
            if self.isCFRadial:
                    v=self.dataset
                    bt=v.variables['time_coverage_start']
                    #tmp=numpy.array(bt[:])
                    #print tmp.shape
                    btv=string1D(bt)
                    ##print type(btv),btv
                    ##sh=len(btv)
                    if btv is None:#len(btv)==0:# is unset
                        #print 'basetime is not set'
                        if require:
                            raise KeyError("Basetime not set")
                        return None
                    self.basetime=datetime.datetime.strptime(btv,self.date_fmt)
            else:
                v=self.dataset
                bt=v.variables['base_time']
                if long(bt.getValue())<=0:# is unset
                    #print 'basetime is not set'
                    if require:
                        raise KeyError("Basetime not set")
                    return None
                epochtime=datetime.datetime(1970,1,1,0,0,0)
                self.basetime=epochtime+datetime.timedelta(seconds=long(bt.getValue()))
        return self.basetime

    def setFileBasetime(self,thetime):
        #self.basetime=thetime
        v=self.dataset
        btv=self.getBasetime(require=False)
        if self.isCFRadial:
            tt=thetime.replace(microsecond=0)
            if self.cfradial_writetime:
                if btv is not None:#not len(btv)==0:# is unset
                    return tt==btv
                # write start of dataset in the form '2012-06-20T00:59:31Z'
                bt=v.variables['time_coverage_start']
                bt[:] = stringtoarr(tt.strftime(self.date_fmt), STRING_LENGTH_SHORT)
                del bt
                self.getBasetime(require=True)
        else:
            if btv is not None:
                return btv==thetime
            epochtime=datetime.datetime(1970,1,1,0,0,0)
            bt=v.variables['base_time']
            if 'dpl_py_binding' in bt.ncattrs():
                del bt.dpl_py_binding
            bt.assignValue(math.trunc((thetime-epochtime).total_seconds()))
            del bt
            btv=self.getBasetime(require=True)
            #self.basetime=epochtime+datetime.timedelta(seconds=long(bt.getValue()))
            v.variables['base_time'].string = btv.strftime(self.date_fmt)#set string value
        codever,codedate=git_tools.getCodeVersion()
        if codever is not None:
            v.setncattr('codeversion',codever)
        if codedate is not None:
            v.setncattr('codedate',codedate.strftime(self.date_fmt))
        for at,val in self.addAttributes.items():
            v.setncattr(at,val) 
        v.sync()
        return True

    def appendtemplatedata(self, dat):
        v=self.dataset#Dataset(self.fn,'a',False)
        appendDimLens=OrderedDict()
        newdimlens=self.appendDims.copy()
        #print self.bindings
        #print self.appendDims
        for adim,adimstart in self.appendDims.items():# start by loading record dimensions that should be appended
            appendDimLens[adim]=adimstart#len(v.dimensions[adim])
            if adim in v.variables and adim in self.bindings:# and getting the variables associated with them stored (reinterpreted as needed)
                dvar=getsubstruct(dat,self.bindings[adim])
                vvar=v.variables[adim]
                if dvar is None:
                    print '***WARNING variable missing from frame! BAD FORM!  VAR ',self.bindings[adim]," -> ",adim
                    continue
                    raise KeyError("Can't find record variable "+self.bindings[adim]+" for dim "+adim)
                try:
                    if len(dvar)==0:
                        continue
                except TypeError:
                    dvar=(dvar,)
                if 'dpl_py_type' in vvar.ncattrs():
                    dpltype=vvar.dpl_py_type[:];
                    if dpltype=='matplotlib_num2date':
                        dvar=[ num2date(d) for d in dvar ]
                        dpltype='python_datetime'
                    if dpltype=='python_datetime':
                        if len(dvar)>0 and self.getBasetime() is None:
                            print 'setting base time from ',adim
                            self.setFileBasetime(dvar[0])
                if self.isCFRadial and adim=='time':
                    # save end_time to nearest second
                    self.end_time = dvar[-1].replace(second=dvar[-1].second, microsecond=0)

        for f in self.bindings:
            if self.bindings[f] is None:
                continue
            fiel=getsubstruct(dat,self.bindings[f])
            if fiel is not None:
                vari=v.variables[f]
                sh=vari.shape
                basesh=[0,0,0,0,0,0]
                #fielsh=fiel.shape
                dims=vari.dimensions
                if 'dpl_py_type' in vari.ncattrs():
                    dpltype=vari.dpl_py_type[:]
                    try:
                        chk=fiel.shape
                    except AttributeError:
                        try:
                            chk=len(fiel)
                        except TypeError:
                            fiel=(fiel,)
                    if dpltype=='matplotlib_num2date':
                        fiel=[ num2date(d) for d in fiel ]
                        dpltype='python_datetime'
                    if dpltype=='python_datetime':
                        if self.getBasetime() is None and len(fiel)>0:
                            self.setFileBasetime(fiel[0])
                        if self.getBasetime() is not None:
                            fiel=[ (d-self.getBasetime()).total_seconds() for d in fiel ]
                            vari.units="seconds since " + self.getBasetime().strftime(self.date_fmt)
                #fixme this is crap. should be a more interpreted way that isn't slow or dangerous
                for didx,dimname in enumerate(dims):
                    if dimname in self.appendDims:
                        basesh[didx]=appendDimLens[dimname]
                try:
                    fsh=fiel.shape
                    if len(fsh)!=len(sh):
                        newsh=[x for x in fsh]
                        while len(newsh)<len(sh):
                            newsh=[1]+newsh
                        #newsh=[(x if x>0 else 1) for x in sh]
                        #print f,sh,fsh,newsh#,fiel
                        fiel=fiel.reshape(newsh)
                        fsh=newsh
                except AttributeError:
                    try:
                        fsh=[len(fiel)]
                        if len(fsh)<len(sh):
                            newsh=[x for x in fsh]
                            while len(newsh)<len(sh):
                                newsh=[1]+newsh
                            #print f,sh,fsh,newsh#,fiel
                            fiel=numpy.array([x for x in fiel]).reshape(newsh)
                            fsh=newsh
                    except TypeError:
                        fsh=[ 1 for x in range(len(basesh))]#'scalar'
                esh=[x for x in fsh]
                fielrange=[]
                for i,x in enumerate(basesh):
                    while len(fielrange)<=i:
                        fielrange.append(slice(None))
                    while len(esh)<=i:
                        esh.append(1)
                    if len(vari.shape)>i:
                        #print 'Shape is ',i,x,vari.shape[i]
                        if vari.shape[i]==0 or x>0:#time dimension
                            pass
                        elif x>vari.shape[i]:
                            print 'WARNING variable is smaller than netcdf!',x,vari.shape[i],self.bindings[f],f
                            x=vari.shape[i]
                        elif vari.shape[i]<fsh[i]:
                            print 'WARNING variable is larger than netcdf!',x,vari.shape[i],fsh[i],self.bindings[f],f
                            fielrange[-1]=slice(0,vari.shape[i]) #this needs verificationFIXME
                            #x+=vari.shape[i]-fsh[i]+1
                    esh[i]+=x
                if isinstance(fiel,(list,tuple)):
                    if len(fielrange)>1:
                        fielrange=[fielrange[0]]
                if hasattr(fiel,'shape'):
                    fielrange=fielrange[:len(fiel.shape)]
                #print f,' from ',self.bindings[f],'dims',sh,basesh,fsh,esh
                #if f=='raw_time' and :
                #    print fiel[0]
                for dimidx,dimname in enumerate(dims):
                    if dimname in self.appendDims:
                        newdimlens[dimname]=esh[dimidx]
                if len(sh)==0:
                    vari.assignValue(fiel)
                else:
                    varii=[]
                    for x in range(len(sh)):
                        varii.append(slice(basesh[x],esh[x]))
                    varii=tuple(varii) if len(varii)>1 else varii[0]
                    fielrange=tuple(fielrange) if len(fielrange)>1 else fielrange[0]
                    #print varii,'hasattrs',vari.ncattrs()
                    #print 'with dtype ',vari.dtype.kind
                    if vari.dtype.kind in ('i','u'):
                        fiel=fiel.copy()
                        if 'missing_value' in vari.ncattrs():
                            fiel[numpy.isnan(fiel)]=vari.getncattr('missing_value')
                        elif vari.dtype.kind in ('i',):
                            fiel[numpy.isnan(fiel)]=-1
                        elif vari.dtype.kind in ('u',):
                            fiel[numpy.isnan(fiel)]=0
                    #print 'varii type ',varii, vari[varii].dtype
                    try:
                        vari[varii]=fiel[fielrange]
                    except IndexError:
                        print 'FAILED TO SLICE',self.bindings[f],f,vari.shape,varii,fielrange,
                        if hasattr(fiel,'shape'):
                            print fiel.shape
                        else:
                            print [len(fiel)]
                        raise
                    except TypeError:
                      try:
                        vari[varii]=fiel
                      except IndexError:
                        print 'FAILED TO DumpSlice',self.bindings[f],f,vari.shape,varii,fielrange,
                        if hasattr(fiel,'shape'):
                            print fiel.shape
                        else:
                            print [len(fiel)]
                        raise

        self.appendDims=newdimlens
        v.sync()
        #v.close()

