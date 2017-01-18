#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
from netCDF4 import Dataset, stringtoarr
import datetime
import tempfile
from lg_dpl_toolbox.core.CDL2NetCDF4 import CDL2NetCDF4
from lg_dpl_toolbox.core.CDL import CDL2pseudonetCDF4

def doKeep(val,keeparr,remarr):
    if keeparr and val not in keeparr:
        return False
    if (keeparr==None or val not in keeparr) and remarr!=None and val in remarr:
        return False
    return True

def copyAttrs(inputobj,outputobj,keepatts=None,removeatts=None):
    for att in inputobj.ncattrs():
        if not doKeep(att,keepatts,removeatts):
            #print 'dropping attr %s' % att
            continue
        #print 'copying attr %s' % att
        outputobj.setncattr(att,inputobj.getncattr(att))


def filterNC(inputnc,outputnc,keepvars=None,keepdims=None,keepatts=None,removevars=None,
             removedims=None,removeatts=None,dimsizes=None):
    """
         copy attributes, dimensions, and variables from inputnc to outputnc 
    """
    copyAttrs(inputnc,outputnc,keepatts,removeatts)

    #then dims
    for dim in inputnc.dimensions:
        if not doKeep(dim,keepdims,removedims):
            print 'dim %s dropped' % dim
            continue
        l=len(inputnc.dimensions[dim])
        if inputnc.dimensions[dim].isunlimited():
            l=None
        if dimsizes and dim in dimsizes:
            l=dimsizes[dim]
        outputnc.createDimension(dim,l)

    #then vars
    for var in inputnc.variables:
        if not doKeep(var,keepvars,removevars):
            print 'var %s dropped' % var
            continue
        #print inputnc.variables[var]
        outputnc.createVariable(var,inputnc.variables[var].dtype,inputnc.variables[var].dimensions,fill_value=inputnc.variables[var].fill_value)
        copyAttrs(inputnc.variables[var],outputnc.variables[var])

def get_substruct( time_z_group, path ):
    """ return the substructure named by 'path' from a time_z_group

    time_z_group : Time_Z_Group, a collection of T_Array, Z_Array and TZ_Array objects
    path : named portion of the Time_Z_Group, e.g. profiles.inv.extinction

    """

    bndclass=path.split('.')
    bndpath=''
    tmp = time_z_group
    # descend through the structure hierarchy, e.g. 'profiles.inv.extinction'
    for ss in bndclass:
        if not hasattr(tmp, ss):
            print 'no %s in structure (path = %s)' % (ss,bndpath)
            return None
        else:
            tmp=getattr(tmp,ss)
            if len(bndpath)==0:
                bndpath=ss
            else:
                bndpath='.'.join((bndpath,ss))
    return tmp

def find_bindings(template, tz_group,bindings_to_keep=None,bindings_to_drop=None):
    """
    return bindings and a dictionary of dimensions found in the tz_group structure 

    template is the sample netCDF file built from the CDL file
    tz_group is a hierarchical structure of HSRL data

    bindings is a dictionary of netcdf variables names to tz_group structure names, e.g. 
    {u'cross_counts': u'rs_raw.cross_pol_counts', }

    """
    bindings={}
    dimsizes = {}
    for v in template.variables:
        if hasattr(template.variables[v],'dpl_py_binding'):
            bnd = template.variables[v].dpl_py_binding 
            if bnd =='dne' or not doKeep(bnd,bindings_to_keep,bindings_to_drop):
                continue
            tmp=get_substruct(tz_group,bnd)
            if tmp!=None:
                bindings[v]=bnd
                dimid=0
                if hasattr(tmp,'shape'):
                    tmpsh=tmp.shape
                else:
                    tmpsh=[len(tmp)]
                vardims=template.variables[v].dimensions
                if len(tmpsh)!=len(vardims) and tmpsh[0]==1:
                    tmpsh=tmpsh[1:]
                if len(tmpsh)!=len(vardims):
                    print 'size of variable %s is different than netcdf %s.' % (bnd,v)
                    print tmpsh
                    print vardims
                for d in vardims:
                    if d in dimsizes:
                        if dimsizes[d]!=tmpsh[dimid] and not template.dimensions[d].isunlimited():
                            print 'dimension %i %s is different sizes in var %s (NC %s). %i vs %i in structure' % \
                                  (dimid,d,bnd,v,dimsizes[d],tmpsh[dimid])
                    else:
                        if template.dimensions[d].isunlimited():
                            dimsizes[d]=None
                            print 'dimension %i %s is unlimited %s (%s) actual size %i' % (dimid,d,bnd,v,tmpsh[dimid])
                        else:
                            dimsizes[d]=tmpsh[dimid]
                            print 'dimension %i %s is now %i from %s (%s)' % (dimid,d,tmpsh[dimid],bnd,v)
                    dimid+=1
    return (bindings, dimsizes)



STRING_LENGTH_SHORT = 32
class DplCreateCfradial(object):
    """
        create a CfRadial format netCDF file HSRL data
    """
    def __init__(self,templatename,output_nc,time_z_group,bindings_to_keep=None):
        """
            write new dimensions and attributes into the specified output_nc object
        """
        if not time_z_group.__dict__.has_key('rs_raw'):
            raise RuntimeError, "no data in time_z_group"
        self.dataset = output_nc

        # create a template netCDF file containing attributes, dimensions and variables  
        # from the specified template
        if False:
            fh, tmpname = tempfile.mkstemp()
            CDL2NetCDF4(templatename, tmpname )
            template=Dataset(tmpname,'r',clobber=False)
        else:
            fh=None
            template=CDL2pseudonetCDF4(templatename)
        varstopurge=None
        dimstopurge=None
        varstokeep=['sweep_number', 'sweep_mode', 'fixed_angle','sweep_start_ray_index','sweep_end_ray_index']

        (self.bindings, dimsizes)=find_bindings(template, time_z_group,bindings_to_keep=bindings_to_keep)

        dimstokeep=[d for d in dimsizes]
        # keep this dimension, even though none of the HSRL products use it, since we need it for
        # time_coverage_start and time_coverage_end
        dimstokeep.append('string_length_short')
        dimstokeep.append('sweep')
        varstokeep=[v for v in self.bindings] + varstokeep
        self.appendDims=[]
        #if 'time' in dimsizes and dimsizes['time']==None:
        self.appendDims.append('time')
        #if 'time' in self.bindings and 'time' in dimsizes:
        varstokeep.append('time_coverage_start')
        varstokeep.append('time_coverage_end')

        filterNC(template,self.dataset,keepvars=varstokeep,keepdims=dimstokeep,
                 removevars=varstopurge,removedims=dimstopurge,dimsizes=dimsizes)
        if fh!=None:
            template.close()
            os.unlink(tmpname)
            os.close(fh)
        self.dataset.sync()
        self.date_fmt = '%Y-%m-%dT%H:%M:%SZ'

    def append_data( self, tzg ):
        """
            append the data found in 'tzg' to our netCDF file
        """

        out=self.dataset
        appendDimLens={}
        for adim in self.appendDims:
            # start by loading record dimensions that should be appended
            appendDimLens[adim]=len(out.dimensions[adim])
            if adim in out.variables:
                    # store the variables associated with this record dimension (processing as needed)
                dvar=get_substruct(tzg,self.bindings[adim])
                nc_ovar=out.variables[adim]
                if dvar==None:
                    raise KeyError("Can't find record variable %s for dim %s",self.bindings[adim],adim)
                if len(dvar)==0:
                    continue
                if adim=='time':#special case for time. this is used with HSRL, 
                    #there needs to be a way to identify this case in the template, 
                    # the what why and how, so other time axes and sources work too
                    # dvar[0] is  datetime.datetime(2012, 6, 20, 0, 59, 31, 250001, tzinfo)
                    if appendDimLens[adim]==0:
                        print 'adding first record'
                        if 'dpl_py_binding' in out.variables['time_coverage_start'].ncattrs():
                            del out.variables['time_coverage_start'].dpl_py_binding
                        # compute start time to the nearest second
                        self.start_time=dvar[0].replace(microsecond=0)
                        # write start of dataset in the form '2012-06-20T00:59:31Z'
                        out.variables['time_coverage_start'][:] = \
                            stringtoarr(self.start_time.strftime(self.date_fmt), STRING_LENGTH_SHORT)

                    # save end_time to nearest second
                    self.end_time = dvar[-1].replace(second=dvar[-1].second, microsecond=0)

        for f in self.bindings:
            field=get_substruct(tzg,self.bindings[f])
            if field!=None:
                ovar=out.variables[f]
                basesh=[0,0,0,0,0,0]
                didx=0

                if 'dpl_py_type' in ovar.ncattrs():
                    dpltype=ovar.dpl_py_type[:]
                    if dpltype=='matplotlib_num2date' or dpltype=='python_datetime': #this is actually datetime, but older form is kept around to not break things JPG 20130211
                        print 'compute relative time for %s' % self.bindings[f]
                        if not hasattr(self,'start_time'):
                            bt=out.variables['time_coverage_start']
                            btv=''
                            for x in range(bt.shape[0]):
                                btv=btv+bt[x]#chartostring(var[:].reshape([1]+list(var[:].shape)))[0]
                            while len(btv)>0 and btv[-1]=='N':
                                btv=btv[:-1]
                            if len(btv)>0:
                                self.start_time=datetime.datetime.strptime(btv,self.date_fmt)
                            else:
                                self.start_time=field[0].replace(microsecond=0)
                                # write start of dataset in the form '2012-06-20T00:59:31Z'
                                out.variables['time_coverage_start'][:] = \
                                    stringtoarr(self.start_time.strftime(self.date_fmt), STRING_LENGTH_SHORT)


                        field=[ (d-self.start_time).total_seconds() for d in field ]
                        if appendDimLens["time"]==0:
                            ovar.units="seconds since " + self.start_time.strftime(self.date_fmt)
                #fixme this is crap. should be a more interpreted way that isn't slow or dangerous
                for dimname in ovar.dimensions:
                    if dimname in self.appendDims:
                        basesh[didx]=appendDimLens[dimname]
                    didx+=1
                print 'Appending variable ',f
                if len(ovar.shape)==0:
                    ovar[:]=field
                elif len(ovar.shape)==1:
                    ovar[basesh[0]:]=field
                else:
                    topsh=[None for x in range(len(basesh))]
                    for x in range(len(field.shape)):
                        topsh[x]=basesh[x]+field.shape[x]
                    print 'appending var',f,field.shape,ovar.shape,basesh,topsh
                    ovar[tuple([slice(basesh[x],topsh[x]) for x in range(len(ovar.shape))])]=field
        out.sync()

    def close(self):
        # write end of dataset time in the form '2012-06-20T00:59:31Z'
        self.dataset.variables['time_coverage_end'][:] = \
            stringtoarr(self.end_time.strftime(self.date_fmt), STRING_LENGTH_SHORT)
        self.dataset.variables['sweep_number'][0] = 0
        self.dataset.variables['sweep_mode'][0] = stringtoarr("pointing", STRING_LENGTH_SHORT)
        self.dataset.variables['fixed_angle'][0] = 0.0
        self.dataset.variables['sweep_start_ray_index'][0] = 0
        self.dataset.variables['sweep_end_ray_index'][0] = self.dataset.variables['time'].shape[0] - 1
        self.dataset.sync()
