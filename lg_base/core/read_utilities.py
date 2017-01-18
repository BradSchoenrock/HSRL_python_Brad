#!/usr/bin/python
# -*- coding: utf-8 -*-

from netCDF4 import Dataset
import bz2
from os import system
import os.path
import numpy as np
from fmt_mpl_date import fmt_mpl_datetime,fmt_mpl_time
#import time
import os
from datetime import datetime, timedelta
import string 
#import hsrl.data_stream.main_routines as mr
from open_config import open_config
import array_utils as hau
#import hsrl.data_stream.input_translators as it
#from hsrl.data_stream.preprocess_and_time_ave import preprocess_and_time_ave
import walkdir as walkdir
import warnings

import json
import struct
from collections import OrderedDict

def pause():
    try:
        print ' '
        print 'press enter to proceed'
        input()
    except:
        print 'proceed'


def convert_date_str(date_string='',twodigityearok_fubar=False):
    from matplotlib.dates import date2num, num2date

 # assume date in format "2-feb-11 18:32"
 # or "30-jan-10 16:53"
 # or "2-jul-11"

    date_string = date_string.strip()

    if date_string == 'now':
        date_string = datetime.utcnow().strftime('%d-%b-%Y %H:%M')

 # add leading zero to day if not supplied

    if not date_string[1].isdigit():
        date_string = '0' + date_string

 # add leading zero to hour if not supplied

    if date_string.find(':') > 0:
        inx = date_string.index(':')
        if date_string[inx - 2] == ' ':
            date_string = date_string[:inx - 2] + ' 0' \
                + date_string[inx - 1:]
    date_string = date_string.replace('-', '')
    date_string = date_string.replace(':', ' ')

    try:
        if len(date_string) == 7:
            if not twodigityearok_fubar:
                warnings.warn("2-digit Years in strings are deprecated",DeprecationWarning)
            dt = datetime.strptime(date_string, '%d%b%y')
        elif len(date_string) == 9:
            dt = datetime.strptime(date_string, '%d%b%Y')
        elif len(date_string) == 13:
            if not twodigityearok_fubar:
                warnings.warn("2-digit Years in strings are deprecated",DeprecationWarning)
            dt = datetime.strptime(date_string, '%d%b%y %H %M')
        elif len(date_string) == 15:
            dt = datetime.strptime(date_string, '%d%b%Y %H %M')
        elif len(date_string) == 16:
            if not twodigityearok_fubar:
                warnings.warn("2-digit Years in strings are deprecated",DeprecationWarning)
            dt = datetime.strptime(date_string, '%d%b%y %H %M %S')
        elif len(date_string) == 18:
            dt = datetime.strptime(date_string, '%d%b%Y %H %M %S')
        else:
            raise SyntaxError, \
                """could not parse %s :
valid formats: 11-mar-11, 1-mar-2011, 1mar11 00:00, 11-mar-2011 00:00:00
 """ \
                % date_string
    except ValueError:
        raise SyntaxError, \
            """could not parse %s :
valid formats: 11-mar-11, 1-mar-2011, 1mar11 00:00, 11-mar-2011 00:00:00
 """ \
            % date_string
    #dt=datetime(*time_struct[:6])
    plot_time = date2num(dt)

    return {'plot_time': plot_time, 'datetime':dt }


def convert_to_matplot_times(hsrl_netcdf_times,ncobject=None):
    from matplotlib.dates import date2num, num2date

    plot_times_tmp=[]
   
    
    n = hsrl_netcdf_times.shape[0]
    for i in range(0,n):
        if hsrl_netcdf_times[i,0] > 0:
            plot_times_tmp.append(date2num(datetime(hsrl_netcdf_times[i,0],
                                                hsrl_netcdf_times[i,1],
                                                hsrl_netcdf_times[i,2],
                                                hsrl_netcdf_times[i,3],
                                                hsrl_netcdf_times[i,4],
                                                hsrl_netcdf_times[i,5],
                  hsrl_netcdf_times[i,6]*1000+hsrl_netcdf_times[i,7])))
            
        else:
            #print '******missing netcdf time-adding 0.5 sec, last t = ' \
            #      ,plot_times_tmp[i-1]
            #plot_times_tmp.append(plot_times_tmp[i-1]+02.5/(3600*24))
            plot_times_tmp.append(np.NaN)
            
    return np.array(plot_times_tmp)

def convert_from_km_to_m(heights_in_km,ncobject=None):
    return heights_in_km[:]*1000.0

def convert_to_python_times(hsrl_netcdf_times,ncobject=None):

    plot_times_tmp=[]
   
    
    n = hsrl_netcdf_times.shape[0]
    if n>0:
        hsrl_netcdf_times=hsrl_netcdf_times[:,:]
    last_time=None
    for i in range(0,n):
        if hsrl_netcdf_times[i,0] > 0:
            try:
                plot_times_tmp.append(datetime(hsrl_netcdf_times[i,0],
                                           hsrl_netcdf_times[i,1],
                                           hsrl_netcdf_times[i,2],
                                           hsrl_netcdf_times[i,3],
                                           hsrl_netcdf_times[i,4],
                                           hsrl_netcdf_times[i,5],
                    hsrl_netcdf_times[i,6]*1000+hsrl_netcdf_times[i,7]))
                last_time=plot_times_tmp[-1]
            except ValueError as e:#invalid date
                print 'Invalid time found step',i,'last good time',last_time,'exception received',e
                plot_times_tmp.append(None)
        else:
            #print '******missing netcdf time-adding 0.5 sec, last t = ' \
            #      ,plot_times_tmp[i-1]
            #plot_times_tmp.append(plot_times_tmp[i-1]+timedelta(seconds=2.5))#plot_times_tmp.append(plot_times_tmp[i-1]+02.5/(3600*24))
            print 'Invalid time found step',i,'last good time',last_time
            plot_times_tmp.append(None)
    #v=np.array(plot_times_tmp)
   
            
    return np.array(plot_times_tmp)

def cf_to_python_times(cf_second_offset,ncobject):

    plot_times_tmp=[]
    epoch=datetime(1970,1,1,0,0,0)
    base_time=epoch+timedelta(seconds=long(ncobject.variables['base_time'].getValue()))
    cf_second_offset=cf_second_offset[:]
    try:
        return np.array([(base_time + timedelta(seconds=x)) for x in cf_second_offset[:]])
    except OverflowError:
        r=[]
        for x in cf_second_offset[:]:
            try:
                r.append(base_time+timedelta(seconds=x))
            except OverflowError:
                r.append(None)
        return np.array(r)

def str_or_set(x):
    if isinstance(x,basestring):
        return set([x])
    return set(x)

def make_dim_mask( nc, dim, dimvar,minv,maxv,translateFunctions,requested_vars):
    if dim not in nc.dimensions: #dim not even found. fail
        print 'no such dim',dim
        return (None,None)
    nc_dims={dim:len(nc.dimensions[dim])}
    dimvarname = dimvar[0] if dimvar[0]!=None else dim
    if dimvarname!=None and dimvarname in requested_vars:
        aliases=str_or_set(requested_vars[dimvarname]['aliases'])
        nc_vars=nc.variables.keys()
        nc_var_set=set(nc_vars)
        #find which one (if any) is in the netcdf file    
        aliases=aliases.intersection(nc_var_set)
    else:
        aliases=set()
    if not aliases: # no variable found to compare to
        if dimvar[0]!=None: #explicit name given. this is unbound, technically a failure
            print 'unbound dim',dim,dimvar, aliases,minv,maxv
            return (0,nc_dims[dim])
        #explicit name not given (fell back to dimension name) but still no aliases. use indexes
        if minv==None:
            minv=0
        if maxv==None or maxv>nc_dims[dim]:
            maxv=nc_dims[dim]
        #print 'bound dim to indexes:',dim,dimvar, aliases,minv,maxv
        return (minv,maxv-minv)
    for alias in aliases:
        alldim = nc.variables[str(alias)]
    if dimvar[1]:
        alldim=translateFunctions[dimvar[1]](alldim,nc)
    mask=np.arange(alldim.shape[0])
    mask=mask[alldim[mask]!=[None for x in alldim[mask]]] #cut nones
    mask=mask[alldim[mask]==alldim[mask]] #cut nans
    if minv!=None:
        if maxv!=None:
            mask=mask[np.logical_and(alldim[mask]>=minv,alldim[mask]<maxv)]
        else:
            mask=mask[alldim[mask]>=minv]
    elif maxv!=None:
        mask=mask[alldim[mask]<maxv]
    if mask.shape[0]==0:
        mini=0
        maxi=0
    else:
        mini=int(min(mask))
        maxi=int(max(mask)-mini+1)
    #print 'dimension',dim,'is found and potentially cropped to',(mini,mini+maxi),'using min',minv,'and max',maxv
    if (alldim[slice(mini,mini+maxi)]==[None for x in alldim[slice(mini,mini+maxi)]]).any():
        ret=np.arange(mini,mini+maxi)
        ret=ret[alldim[ret]!=[None for x in alldim[ret]]]
        return ret
    return (mini,maxi)

def read_raw( filename, instrument,requested_vars, max_range_bin,pre_process,verbose=False,firsttime=None,lasttime=None,firstrange=None,lastrange=None,firstaltitude=None,lastaltitude=None,doread=True,user_read_mode=None):
    """read raw data from netcdf filename
       instrument, e.g. 'gvhsrl'
       max_range_bin  = read this number of range bins
       requested_vars  = determines how many variable to read from netcdf file
                       eg. 'images' , 'housekeeping'          
       returns TZ_Array
       """


   

    # print 'start time num= ',start_time_num, num2date(start_time_num)
    # FIXME - compressed files probably don't work on OS X because of /dev/shm

    removeTempFile=False
    TempFile=None
    config=requested_vars['config']
    requested_vars=requested_vars['selected_vars']
    found_vars=None
    rs_raw=None

    if filename.find('.bz2') >= 0:
        print 'WARNING decompression should be handled elsewhere. not in readraw. are you using the zookeeper?'
        bzf=bz2.BZ2File(filename,'r')
        if os.access(os.path.join('/dev','shm'),os.W_OK):
            randomformat='>Q'
            rnd=open('/dev/urandom','r')
            rndn=rnd.read(struct.calcsize(randomformat))
            rnd.close()
            randomtag='%08x' % struct.unpack(randomformat,rndn)
            filename=os.path.join('/dev','shm','readraw_' + instrument + '_' 
                                  + randomtag + '.nc')
        else:
            filename='tmp_readraw_' + instrument + '.nc'
        fileid=open(filename,'w')
        while True:
            tmp=bzf.read(1024*1024*64)
            if len(tmp)==0:
                break
            fileid.write(tmp)
        fileid.close()
        bzf.close()
        TempFile=filename
        removeTempFile=True

    try:
        read_mode = str_or_set(config['read_mode']) if user_read_mode==None else str_or_set(user_read_mode)
        timedim = config['timedim'] if 'timedim' in config else None
        rangedim = config['rangedim'] if 'rangedim' in config else None
        altitudedim = config['altitudedim'] if 'altitudedim' in config else None
        timerecordvariable = config['timedimvar'] if 'timedimvar' in config else (None,None)
        rangerecordvariable = config['rangedimvar'] if 'rangedimvar' in config else (None,None)
        altituderecordvariable = config['altitudedimvar'] if 'altitudedimvar' in config else (None,None)

        translateFunctions={'convert_hsrltime_to_datetime':convert_to_python_times,
                    'convert_cf_basetime_to_datetime':cf_to_python_times,
                    'convert_from_km_to_m':convert_from_km_to_m}

        tzgroup_parms={}
        if timerecordvariable[0]!=None:
            tzgroup_parms['timevarname']=timerecordvariable[0]
        if altituderecordvariable[0]!=None:
            tzgroup_parms['altname']=altituderecordvariable[0]
        elif rangerecordvariable[0]!=None:
            tzgroup_parms['altname']=rangerecordvariable[0]
            #FIXME where does msl_altitudes come from (current default for alt dim)
       
        try:
            nc = Dataset(filename,'r')
        except RuntimeError as msg:
            times=[]
            rs_raw=hau.Time_Z_Group(hau.T_Array(times),**tzgroup_parms)
            print 'read_raw: error reading netCDF file:',msg, filename
            #raise
            #return rs_raw
            raise IOError('Invalid NetCDF file '+filename)
           
        found_vars = OrderedDict()
        globalattributes = OrderedDict()
        for a in nc.ncattrs():
            globalattributes[a]=getattr(nc,a)#nc.attr(a).get()

        #get list of variables in the netcdf file
        nc_vars = nc.variables.keys()
           
        # create new rs_raw
        rs_raw = hau.Time_Z_Group(**tzgroup_parms)    

      
        nc_var_set=set(nc_vars)
        
        nc_dims=OrderedDict()
        for n,d in nc.dimensions.items():
            nc_dims[n]=(0,len(d))
        if timedim:
            nc_dims[timedim]=make_dim_mask(nc,timedim,timerecordvariable,firsttime,lasttime,translateFunctions,requested_vars)
        if rangedim:
            if rangerecordvariable[0]==None and lastrange==None:
                lastrange=max_range_bin
            nc_dims[rangedim]=make_dim_mask(nc,rangedim,rangerecordvariable,firstrange,lastrange,translateFunctions,requested_vars)
            if rangerecordvariable[0]!=None or lastrange!=None:
                if nc_dims[rangedim][0]!=None and (nc_dims[rangedim][0]+nc_dims[rangedim][1])>max_range_bin:
                    nc_dims[rangedim]=(nc_dims[rangedim][0],max_range_bin-nc_dims[rangedim][0]+1)
                    if nc_dims[rangedim][1]<=0:
                        nc_dims[rangedim]=(0,0)
        if altitudedim:
            nc_dims[altitudedim]=make_dim_mask(nc,altitudedim,altituderecordvariable,firstaltitude,lastaltitude,translateFunctions,requested_vars)
        #loop through the list of requested variable names
        #note: name may be under alias
        #print requested_vars
        for k in nc_dims.keys():
            if not isinstance(nc_dims[k],tuple):
                continue
            v=[nc_dims[k][0],nc_dims[k][1]]
            if v[0]==None:
                v[0]=0
            if v[1]!=None:
                v[1]+=v[0]
            nc_dims[k]=slice(v[0],v[1])


        for name in requested_vars:
        
            #get the aliases for this name
            aliases=str_or_set(requested_vars[name]['aliases'])
            
            #find which one (if any) is in the netcdf file    
            aliases=aliases.intersection(nc_var_set)

            #print
            #print 'variable found under this name---->',aliases
            
            #aliases is a set of one or zero items
            
            if aliases:
                for alias in aliases:
                  
                    #read sel is list specifying which read_modes triger reading
                    read_sel=str_or_set(requested_vars[name]['read_sel'])
                    read_sel.add('all')
                    if len(read_sel.intersection(read_mode)):
                        alias=str(alias)
                        found_vars[name] = OrderedDict()
                        var=nc.variables[alias]
                        for k in var.ncattrs():
                            found_vars[name][k]=getattr(var,k)
                        if doread==False:
                            continue
                        #starta=[]
                        #counta=[]
                        iarr=[]
                        for d in var.dimensions:
                            #starta.append(nc_dims[d][0])
                            #counta.append(nc_dims[d][1])
                            #if nc_dims[d][0]==None or nc_dims[d][1]==0 or nc_dims[d][1]==None:
                            #    if removeTempFile:
                            #        os.unlink(TempFile)
                            #    return None,None,None
                            iarr.append(nc_dims[d])#slice(nc_dims[d][0],nc_dims[d][1]))
                        #print 'reading variable',alias,'start',starta,'count',counta
                        #get requested variable from netcdf
                        try:
                            if len(iarr)==0:
                                tmp=var.getValue()
                            elif len(iarr)==1:
                                tmp=var[iarr[0]]
                            else:
                                tmp=var[tuple(iarr)]
                        except IndexError:
                            print 'Index error reading variable ',alias,' in file ',filename
                            continue
                        #print 'read_raw-----------------',name
                        #if name == 'filtered_energy':
                        #    print tmp
                        summode=None
                        if 'summode' in requested_vars[name]:
                            summode=requested_vars[name]['summode']
                        vdims=var.dimensions
                        if name == timerecordvariable[0]:
                            if not np.size(tmp) == 0:
                                if timerecordvariable[1]:
                                    tmp = translateFunctions[timerecordvariable[1]](tmp,nc)
                                setattr(rs_raw,timerecordvariable[0],hau.T_Array(tmp,summode=summode))
                            else:    
                                rs_raw=hau.Time_Z_Group(hau.T_Array([]),**tzgroup_parms)
                                if removeTempFile:
                                    os.unlink(TempFile)
                                return rs_raw, None,None
                        elif len(vdims) >= 2 and vdims[1] == rangedim and vdims[0] == timedim:
                            setattr(rs_raw,name,hau.TZ_Array(tmp,summode=summode))
                        elif len(vdims) >= 2 and vdims[1] == altitudedim and vdims[0] == timedim:
                            setattr(rs_raw,name,hau.TZ_Array(tmp,summode=summode))
                        elif len(vdims) > 0 and vdims[0] == timedim:
                            setattr(rs_raw,name,hau.T_Array(tmp,summode=summode))
                        elif len(vdims) > 0 and vdims[0] == rangedim:
                            setattr(rs_raw,name,hau.Z_Array(tmp,summode=summode))
                        elif len(vdims) > 0 and vdims[0] == altitudedim:
                            setattr(rs_raw,name,hau.Z_Array(tmp,summode=summode))
                        else:
                            setattr(rs_raw,name,tmp)

                        #make dictionary found variables and their attributes                
                        #found_vars[name] = nc.var(alias).attributes()
                    elif verbose:
                        print name,' not found in netcdf'


                       
        nc.close()

    finally:
        if removeTempFile:
            os.unlink(TempFile)

    if doread==False:
        return None,found_vars,globalattributes
   
                    
    if pre_process!=None:     
        pre_process(rs_raw) #do preprocessing and time average
        
    return rs_raw,found_vars,globalattributes

def readascii(filename):

 # [header,data]=readascii(filename)
 # The input file can not contain blank lines. This conforms to  xmgr usage.
 # The file must begin with one or more comment lines with a '#' in the first
 # column. An ascii data file with an arbitray number of rows and columns must
 # follow. Numbers must be separated with spaces or tabs and lines must end
 # with a linefeed. The routine returns a string array 'header' which
 # contains the comment  lines and a numeric array 'data' containing the data.
    if not filename:
        print '   requested file, ",',filename,',", is empty'
        return ['', []]

    header = ''
    data = None
    fileid = file(filename, 'r')
    for line in fileid:

        if line[0] == '#':
            header = header + line
        else:
            if line.find(',')>0:
                temp = line.split(',')
            else:    
                temp = line.split()

            if data is None:
                n = len(temp)
                data = np.zeros(n, float)
                d2 = np.zeros( n, float)
                for i in range(0, n):
                    data[i] = float(temp[i])
            else:
                for i in range(0, n):
                    d2[i] = float(temp[i])
                data = np.vstack([data, d2])
    fileid.close()
    return [header, data]

