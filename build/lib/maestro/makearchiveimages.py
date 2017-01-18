#!/usr/bin/python
# -*- coding: utf-8 -*-

import matplotlib
matplotlib.use('Agg')
#from hsrl.data_stream.rti import Rti
from datetime import datetime,timedelta
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import os
import sys
import lg_base.core.open_config as oc
import json
import multiprocessing
import threading
import time
import lg_base.core.json_config as jc
from lg_base.core.locate_file import locate_file
import copy
from netCDF4 import Dataset
import traceback
from collections import OrderedDict

from PIL import Image, ImageChops

import subprocess as sp

def ensureOnlyOne(filename,extensionsToCheck):
    myext=None
    for ext in extensionsToCheck:
        if filename.endswith(ext):
            myext=ext
            break
    if myext==None:
        return
    mostname=filename[:-(len(myext))]
    for ext in extensionsToCheck:
        if ext==myext:
            continue
        fn=mostname+ext
        if os.path.exists(fn):
            os.unlink(fn)

def convArr(arr,skip=0):
    if hasattr(arr,'shape'):
        if len(arr.shape)>0 and arr.size>0:
            if skip==0:
               return [convArr(x) for x in arr]
            ret=[]
            for i,x in enumerate(arr):
               if (i%(skip+1))==0:
                  ret.append(convArr(x))
            return ret
    try:
       return float(arr)
    except:
       try:
          return (arr-datetime(2010,1,1,0,0,0)).total_seconds()
       except:
          pass
    return arr

def makedict(struc,varlist,skip=0):
    ret=OrderedDict()
    for v in varlist:
        tmp=ret
        vl=v.split('.')
        ts=struc
        for sv in vl[:-1]:
            if ts is None or not hasattr(ts,sv):
                ts=None
                break
            ts=getattr(ts,sv)
            if not sv in tmp:
                tmp[sv]=OrderedDict()
            tmp=tmp[sv]
        if ts is None or not hasattr(ts,vl[-1]):
            continue
        tmp[vl[-1]]=convArr(getattr(ts,vl[-1]),skip=skip)
    ret['skip']=skip
    return ret

def trydecompress(strv):
    if strv[0]!='{':
        import bz2 as compressor
        import base64
        #strv=base64.b64decode(strv)
        #decompress
        strv=compressor.decompress(strv)
    return json.loads(strv, object_pairs_hook=OrderedDict)

def dumpExif(filename):
    #im=Image.open(filename)
    #print im.info['exif']
    import pexif 
    jpeg=pexif.JpegFile.fromFile(filename,'ro')
    exif = jpeg.exif
    tmp=''.join(exif.primary.ExtendedEXIF.UserComment)
    print trydecompress(tmp)

def trycompress(val):
    import bz2 as compressor
    import base64
    cval=compressor.compress(val)
    #cval=base64.b64encode(cval)
    if len(cval)<len(val):
        return cval
    return val

def addExifComment(filename,val):
    import pexif
    jpeg=pexif.JpegFile.fromFile(filename,'rw')
    exif = jpeg.exif
    if isinstance(val,dict):
        import json
        val=json.dumps(val,separators=(',',':'))
    cval=trycompress(val)
    print 'exif length is ',len(cval)
    exif.primary.ExtendedEXIF.UserComment=cval
    jpeg.writeFile(filename+'tmp')
    os.rename(filename+'tmp',filename)

def recompressImage(filename,*args,**kwargs):
    im=Image.open(filename)
    im.save(filename,*args,**kwargs)

def sendFTPFile(source,host,path,user,password):
    #print 'would send ',source
    #return False
    null = file(os.path.devnull, 'w')
    rc = sp.call(['/usr/bin/ncftpput','-u',user,'-p',password,host,path,source], stdout=null, stderr=null)
    print 'send result is ',rc
    return rc==0 

 
def sendSCPFile(source,host,path,user,keyfile):
    assert(keyfile is not None)
    null = file(os.path.devnull, 'w')
    nulli = file(os.path.devnull, 'r')
    cmd=['/usr/bin/scp','-o','PreferredAuthentications=publickey','-q','-B','-i',keyfile,os.path.abspath(source),'%s@%s:%s' % (user,host,path)]
    print cmd
    rc = sp.call(cmd, stdout=null, stderr=null,stdin=nulli)#,env={})
    print 'send result is ',rc
    return rc==0 
 

def trimborder(im):
    isfilename=isinstance(im,basestring)
    if isfilename:
        imname=im
        im=Image.open(imname)
    bg = Image.new(im.mode, im.size, im.getpixel((0,0)))
    print im.size
    diff = ImageChops.difference(im, bg)
    diff = ImageChops.add(diff, diff, 2.0, -100)
    bbox = diff.getbbox()
    if bbox:
        newim=im.crop(bbox)
        if isfilename and newim:
            newim.save(imname)
            print newim.size
    else:
        newim=None

def tinythumb(valsorig,cscale,stringheader,dpi,qcmask=None):
    plt.rc('font', size=7)
    #thmbplt=plt.figure(figureno,figsize=(1.53,2.3),dpi=usedpi)
    thmbplt=plt.figure(figsize=(1.135,1.6),dpi=dpi)
    vals=copy.copy(valsorig)
    vals[vals < cscale[0]] = cscale[0]
    vals[vals > cscale[1]] = cscale[1]
    if qcmask!=None:
        vals[np.logical_not(qcmask)]=np.nan
    #imgplot=plt.imshow(np.flipud(np.transpose(vals)),aspect=.645)
    thmbplt.clear()
    ax=thmbplt.add_axes([0.0,0.0,1.0,0.9])
    imgplot=ax.imshow(np.flipud(np.transpose(vals)),aspect=1.78)
    imgplot.set_clim(cscale)
    #thmbplt.tight_layout()
    myjetcm = cm.jet

    # myjetcm.set_over('k',alpha=.2)

    #myjetcm.set_bad('k', alpha=1)
    myjetcm.set_bad('k', alpha=1)
    #myjetcm.set_under([.3, .3, .3], alpha=.1)
    myjetcm.set_under('k',alpha=1)
    imgplot.set_cmap(myjetcm)
    ax.set_title(stringheader)
    ax.set_xticks([],'')
    ax.set_yticks([],'')
    #ax=imgplot.get_axes()
    print ax
    if not plt.isinteractive():
       plt.show(block=False)
    thmbplt.set_size_inches(thmbplt.get_size_inches(),forward=True)
    thmbplt.canvas.draw()
    #help(thmbplt)
    return thmbplt
    #help(ax)

def ignoreSomeGroups(name):
    return 'raw' in name or 'profile' in name

from lg_base.core.timeout import completion_timeout

@completion_timeout(600,killOnTimeout=True)
def makeArchiveImages(instrument,datetimestart,range_km=15,full24hour=None,filename=None,reprocess=False,attenuated=False,frappe=False,ir1064=False,ismf2ship=False,completeframe=None,*args,**kwargs):
    import lg_dpl_toolbox.filters.substruct as frame_substruct
    frame_substruct.SubstructBrancher.multiprocessable=False
    from lg_dpl_toolbox.dpl.dpl_read_templatenetcdf import dpl_read_templatenetcdf
    from lg_dpl_toolbox.filters.time_frame import FrameCachedConcatenate
    from hsrl.dpl.dpl_hsrl import dpl_hsrl
    from radar.dpl.dpl_radar import dpl_radar
    from lg_dpl_toolbox.dpl.dpl_create_templatenetcdf import dpl_create_templatenetcdf
    import hsrl.dpl.dpl_artists as hsrl_artists
    import radar.dpl.dpl_artists as radar_artists
    import raman.dpl.dpl_artists as raman_artists

    import lg_dpl_toolbox.core.archival as hru
    import hsrl.graphics.hsrl_display as du
    if filename!=None:
        useFile=True
    else:
         useFile=False #True
    if not useFile or not os.access(filename,os.F_OK):
        reprocess=True
    realend=None
    #process_control=None
    if reprocess:
        print datetime
        #instrument='ahsrl'
        if useFile:
            n=Dataset(filename,'w',clobber=True)
            n.instrument=instrument
        print datetime
        realstart=datetimestart
        if full24hour==None:
            realstart=datetimestart.replace(hour=0 if datetimestart.hour<12 else 12,minute=0,second=0,microsecond=0)
        elif frappe:
            realstart=datetimestart.replace(hour=0,minute=0,second=0,microsecond=0)
        elif ismf2ship:
            realend=realstart
            realstart=realstart-timedelta(days=2.0)
        else:
            realstart=realstart-timedelta(days=1.0)
        if 'realstart' in kwargs:
            realstart=kwargs['realstart']
        if realend is None:
            realend=realstart+timedelta(days=.5 if full24hour==None else 1.0)
        if 'realend' in kwargs:
            realend=kwargs['realend']
        isHsrl=False
        isRadar=False
        isRaman=False
        instrumentprefix=None
        if instrument.endswith('hsrl'):
            dpl=dpl_hsrl(instrument=instrument,filetype='data')
            dpl.hsrl_process_control.set_value('quality_assurance','enable',False)
            #dpl.hsrl_process_control.set_value('extinction_processing','enable',False)
            dpl.hsrl_process_control.set_value('wfov_corr','enable',False)
            gen=dpl(start_time_datetime=realstart,end_time_datetime=realend,min_alt_m=0,max_alt_m=range_km*1000,with_profiles=False)
            isHsrl=True
            #process_control=gen.hsrl_process_control
            hsrlinstrument=instrument
            if os.getenv('COMPLETEFRAME',completeframe)==None:
                import lg_dpl_toolbox.filters.substruct as frame_substruct
                dropcontent=['rs_raw','rs_mean']
                gen=frame_substruct.DropFrameContent(gen,dropcontent)

        elif ('kazr' in instrument) or ('mwacr' in instrument) or instrument=='mmcr':
            dpl=dpl_radar(instrument=instrument)
            gen=dpl(start_time_datetime=realstart,end_time_datetime=realend,min_alt_m=0,max_alt_m=range_km*1000,forimage=True,allow_nans=True)
            hsrlinstrument=dpl.instrumentbase
            isRadar=True
            if 'mwacr' in instrument:
                instrumentprefix='mwacr'
            else:
                instrumentprefix='radar'
            #merge=picnicsession.PicnicProgressNarrator(dplc,getLastOf('start'), searchparms['start_time_datetime'],searchparms['end_time_datetime'],session)
            #hasProgress=True
        elif instrument.startswith('rlprof'):
            from raman.dpl.raman_dpl import dpl_raman
            import lg_dpl_toolbox.filters.time_frame as time_slicing
            dpl=dpl_raman('bagohsrl',instrument.replace('rlprof',''))
            gen=dpl(start_time_datetime=realstart,end_time_datetime=realend,min_alt_m=0,max_alt_m=range_km*1000,forimage=True,allow_nans=True,inclusive=True)
            import functools
            import lg_base.core.array_utils as hau
            from dplkit.simple.blender import TimeInterpolatedMerge
            import lg_dpl_toolbox.filters.substruct as frame_substruct
            import lg_dpl_toolbox.dpl.TimeSource as TimeSource
            import lg_base.core.canvas_info as ci
            gen=time_slicing.TimeGinsu(gen,timefield='times',dtfield=None)
            forimage=ci.load_canvas_info()['canvas_pixels']
            timesource=TimeSource.TimeGenerator(realstart,realend,time_step_count=forimage['x'])
            gen=TimeInterpolatedMerge(timesource,[gen], allow_nans=True)
            gen=TimeSource.AddAppendableTime(gen,'times','delta_t')
            gen=frame_substruct.Retyper(gen,functools.partial(hau.Time_Z_Group,timevarname='times',altname='altitudes'))
            hsrlinstrument='bagohsrl'
            instrumentprefix=instrument
            isRaman=True


    if reprocess and useFile:
        v=None

        for i in gen:
            if v==None:
                v=dpl_create_templatenetcdf(locate_file('hsrl_nomenclature.cdl'),n,i)
            v.appendtemplatedata(i)
        #raise TypeError
        #break
            n.sync()
        n.close()

        #fn='outtest.nc'

    if useFile:
        v=dpl_read_templatenetcdf(filename)

        instrument=v.raw_netcdf().instrument[:]
    else:
        v=gen

    defaultjson='archive_plots.json'
    if ismf2ship:
        defaultjson='mf2ship_plots.json'
    (disp,conf)=du.get_display_defaults(os.getenv("PLOTSJSON",defaultjson))
  
    v=FrameCachedConcatenate(v)
    import lg_dpl_toolbox.dpl.TimeSource as TimeSource
    import lg_dpl_toolbox.filters.fill as fill
    import lg_base.core.canvas_info as ci
    forimage=ci.load_canvas_info()['canvas_pixels']
    ts=TimeSource.TimeGenerator(realstart,realend,time_step_count=forimage['x'])
    v=fill.FillIn(v,[np.array([x['start'] for x in ts])],ignoreGroups=ignoreSomeGroups)

    rs=None
    for n in v:
        if rs==None:
            rs=n
        else:
            rs.append(n)
        if isHsrl:
            if attenuated:
                bsfigname='atten_backscat_image'
                bsprefix='attbscat'
                disp.set_value(bsfigname,'enable',1)
                disp.set_value('backscat_image','enable',0)
                disp.set_value('linear_depol_image','figure',bsfigname)
                #disp.set_value('circular_depol_image','figure',bsfigname)
            else:
                bsfigname='backscat_image'
                bsprefix='bscat'
            if ir1064 or (hasattr(rs,'rs_inv') and hasattr(rs.rs_inv,'color_ratio')):
                disp.set_value('raw_color_ratio_image','enable',1)
                disp.set_value('color_ratio_image','enable',1)
            #if hasattr(rs.rs_inv,'circular_depol'):
            #    disp.set_value('linear_depol_image','enable',0) #this modification works because the artist doesn't render until after the yield
            #    field='circular_depol'
            #else:
            #disp.set_value('circular_depol_image','enable',0)
            field='linear_depol'


    if isHsrl:
        #if process_control==None:
        #    print 'loading process control from json'
        #    process_control=jc.json_config(locate_file('process_control.json'),'process_defaults')
        v=hsrl_artists.dpl_images_artist(v,instrument=instrument,max_alt=None,display_defaults=disp)
    elif isRadar:
        v=radar_artists.dpl_radar_images_artist(framestream=v,instrument=v.radarType,display_defaults=disp,subframe=None)
    elif isRaman:
        v=raman_artists.dpl_raman_images_artist(framestream=v,display_defaults=disp)

    for n in v:
        pass #run once more with cached value and actual artists

    usedpi=90

    #rs=Rti('ahsrl',stime.strftime('%d-%b-%y %H:%M'),dtime.total_seconds()/(60*60),minalt_km,maxalt_km,.5,'archive_plots.json')
    imtime=realstart
    imetime=realend
    imtimestr=imtime.strftime('%e-%b-%Y ')
    file_timetag=imtime.strftime("_%Y%m%dT%H%M_")+imetime.strftime("%H%M_")+("%i" % range_km)
    if not full24hour:
        imtimestr+= ('AM' if imtime.hour<12 else 'PM' )
        file_timetag+="_am" if imtime.hour<12 else "_pm"
    filelocation=hru.get_path_to_data(hsrlinstrument,None)
    print filelocation
    filelocation=os.path.join(filelocation,imtime.strftime("%Y/%m/%d/images/"))
    try:
        os.makedirs(filelocation)
    except OSError:
        pass
    print filelocation
    if os.getenv("DEBUG",None):
        filelocation='.'
        print 'DEBUG is on. storing in current directory'
    figs=v.figs
    extensionsList=('.png','.jpg','.gif','.jpeg')
    preferredFormat='jpg'
    preferredExtension='.'+preferredFormat
    if isHsrl:
        #figs=du.show_images(instrument,rs,None,{},process_control,disp,None,None,None)

        #for x in figs:
        #    fig = figs.figure(x)
        #    fig.canvas.draw()
        
        print figs
        f=figs.figure(bsfigname)
        #f.set_size_inches(f.get_size_inches(),forward=True)
        if frappe:
            frappetag=imtime.strftime('upperair.UW_HSRL.%Y%m%d%H%M.BAO_UWRTV_')
            hiresfile=os.path.join('.',frappetag+'backscatter_%ikm%s' %(range_km,preferredExtension))
            ensureOnlyOne(hiresfile,extensionsList)
            f.savefig(hiresfile,format=preferredFormat,bbox_inches='tight')
        else:
            hiresfile=os.path.join(filelocation,bsprefix+"_depol" + file_timetag + preferredExtension)
            ensureOnlyOne(hiresfile,extensionsList)
            f.savefig(hiresfile,format=preferredFormat,bbox_inches='tight')
        if disp.enabled('raw_color_ratio_image'):
            f=figs.figure('color_ratio_image')
            #f.set_size_inches(f.get_size_inches(),forward=True)
            hiresfileir=os.path.join(filelocation,"ratioIR" + file_timetag + preferredExtension)
            ensureOnlyOne(hiresfileir,extensionsList)
            f.savefig(hiresfileir,format=preferredFormat,bbox_inches='tight')
        figs.close()

        if not full24hour and hasattr(rs,'rs_inv'):
            scale = [float(disp.get_value(bsfigname,'lo_color_lmt')),
                     float(disp.get_value(bsfigname,'hi_color_lmt'))]
            depol=100*getattr(rs.rs_inv,field)
            depolscale = [float(disp.get_value(field+'_image', 'lo_color_lmt')),
                          float(disp.get_value(field+'_image', 'hi_color_lmt'))]

            if attenuated:
                backscat=rs.rs_inv.atten_beta_a_backscat
            elif hasattr(rs.rs_inv,'beta_a_backscat'):
                backscat=rs.rs_inv.beta_a_backscat
            else:
                backscat=rs.rs_inv.beta_a_backscat_par + rs.rs_inv.beta_a_backscat_perp 
            print rs.rs_inv.times.shape
            print backscat.shape
            if disp.get_value(bsfigname,"log_linear")=='log':
                scale=np.log10(scale)
                backscat[backscat<=0]=np.NAN;
                backscat=np.log10(backscat)
            if disp.get_value(field+'_image',"log_linear")=='log':
                depol[depol<=0]=np.NAN;
                depolscale=np.log10(depolscale)
                depol=np.log10(depol)
            print backscat.shape
            qc_mask=None
            if hasattr(rs.rs_inv,'qc_mask'):
               qc_mask=np.ones_like(rs.rs_inv.qc_mask)           
               qcbits={'mol_lost':64,'mol_sn_ratio':16,'cloud_mask':128,'I2_lock_lost':4}
               for name,maskbit in qcbits.items():
                 if disp.get_value('mask_image',name):
                   qc_mask = np.logical_and(rs.rs_inv.qc_mask & maskbit > 0,qc_mask)
            #print np.sum(backscat<=0)
            f=tinythumb(backscat,scale,imtimestr,dpi=usedpi,qcmask=qc_mask)
            thumbname=os.path.join(filelocation,bsprefix + file_timetag + '_thumb'+preferredExtension)
            print thumbname
            ensureOnlyOne(thumbname,extensionsList)
            f.savefig(thumbname,format=preferredFormat,bbox_inches='tight',dpi=usedpi)
            trimborder(thumbname)
            
            f=tinythumb(depol,depolscale,imtimestr,dpi=usedpi,qcmask=qc_mask)
            thumbname=os.path.join(filelocation,'depol' + file_timetag + '_thumb'+preferredExtension)
            ensureOnlyOne(thumbname,extensionsList)
            f.savefig(thumbname,format=preferredFormat,bbox_inches='tight',dpi=usedpi)
            trimborder(thumbname)

            if ir1064 or hasattr(rs.rs_inv,'color_ratio'):
                if hasattr(rs.rs_inv,'qc_mask'):
                   qc_mask=np.ones_like(rs.rs_inv.qc_mask)           
                   qcbits={'mol_lost':64,'mol_sn_ratio':16,'cloud_mask':128,'I2_lock_lost':4,'1064_shutter':0x8000}
                   for name,maskbit in qcbits.items():
                        if disp.get_value('mask_image',name):
                            qc_mask = np.logical_and(rs.rs_inv.qc_mask & maskbit > 0,qc_mask)
                cr=rs.rs_inv.color_ratio
                scale = [float(disp.get_value('color_ratio_image','lo_color_lmt')),
                         float(disp.get_value('color_ratio_image','hi_color_lmt'))]
                if disp.get_value('color_ratio_image',"log_linear")=='log':
                    scale=np.log10(scale)
                    cr[cr<=0]=np.NAN;
                    cr=np.log10(cr)
                f=tinythumb(cr,scale,imtimestr,dpi=usedpi,qcmask=qc_mask)
                thumbname=os.path.join(filelocation,'ratioIR' + file_timetag + '_thumb'+preferredExtension)
                print thumbname
                ensureOnlyOne(thumbname,extensionsList)
                f.savefig(thumbname,format=preferredFormat,bbox_inches='tight',dpi=usedpi)
                trimborder(thumbname)


    elif isRadar:
        f=figs.figure('radar_backscatter_image')
        #f.set_size_inches(f.get_size_inches(),forward=True)
        hiresfile=os.path.join(filelocation,instrumentprefix+"_bscat" + file_timetag + preferredExtension)
        ensureOnlyOne(hiresfile,extensionsList)
        f.savefig(hiresfile,format=preferredFormat,bbox_inches='tight')
        figs.close()
        if not full24hour:
            radarscale = [float(disp.get_value('radar_backscatter_image', 'lo_color_lmt')),
                          float(disp.get_value('radar_backscatter_image', 'hi_color_lmt'))]
            radar=rs.Backscatter
            if disp.get_value('radar_backscatter_image','log_linear')=='log':
                radar[radar<=0]=np.NAN
                radarscale=np.log10(radarscale)
                radar=np.log10(radar)
            f=tinythumb(radar,radarscale,imtimestr,dpi=usedpi)
            thumbname=os.path.join(filelocation,instrumentprefix+'_bscat' + file_timetag + '_thumb'+preferredExtension)
            ensureOnlyOne(thumbname,extensionsList)
            f.savefig(thumbname,format=preferredFormat,bbox_inches='tight',dpi=usedpi)
            trimborder(thumbname)
    elif isRaman:
        if hasattr(rs,'backscatter'):
            f=figs.figure('rl_backscatter_image')
            #f.set_size_inches(f.get_size_inches(),forward=True)
            hiresfile=os.path.join(filelocation,instrumentprefix+"_bscat" + file_timetag + preferredExtension)
            ensureOnlyOne(hiresfile,extensionsList)
            f.savefig(hiresfile,format=preferredFormat,bbox_inches='tight')
        elif hasattr(rs,'beta'):
            f=figs.figure('rl_backscatter_image')
            #f.set_size_inches(f.get_size_inches(),forward=True)
            hiresfile=os.path.join(filelocation,instrumentprefix+"_beta" + file_timetag + preferredExtension)
            ensureOnlyOne(hiresfile,extensionsList)
            f.savefig(hiresfile,format=preferredFormat,bbox_inches='tight')
        if hasattr(rs,'linear_depol'):
            f=figs.figure('rl_depol_image')
            #f.set_size_inches(f.get_size_inches(),forward=True)
            hiresfile=os.path.join(filelocation,instrumentprefix+"_dep" + file_timetag + preferredExtension)
            ensureOnlyOne(hiresfile,extensionsList)
            f.savefig(hiresfile,format=preferredFormat,bbox_inches='tight')
        figs.close()
        if not full24hour:
            if hasattr(rs,'backscatter'):
                ramanscale = [float(disp.get_value('rl_backscatter_image', 'lo_color_lmt')),
                              float(disp.get_value('rl_backscatter_image', 'hi_color_lmt'))]
                raman=rs.backscatter.copy()
                #print np.nanmax(raman),np.nanmin(raman),' is ramans actual range'
                if disp.get_value('rl_backscatter_image','log_linear')=='log':
                    raman[raman<=0]=np.NAN
                    ramanscale=np.log10(ramanscale)
                    raman=np.log10(raman)
                f=tinythumb(raman,ramanscale,imtimestr,dpi=usedpi)
                thumbname=os.path.join(filelocation,instrumentprefix+'_bscat' + file_timetag + '_thumb'+preferredExtension)
                ensureOnlyOne(thumbname,extensionsList)
                f.savefig(thumbname,format=preferredFormat,bbox_inches='tight',dpi=usedpi)
                trimborder(thumbname)
            if hasattr(rs,'beta'):
                ramanscale = [float(disp.get_value('rl_backscatter_image', 'lo_color_lmt')),
                              float(disp.get_value('rl_backscatter_image', 'hi_color_lmt'))]
                raman=rs.beta.copy()
                #print np.nanmax(raman),np.nanmin(raman),' is ramans actual range'
                if disp.get_value('rl_backscatter_image','log_linear')=='log':
                    raman[raman<=0]=np.NAN
                    ramanscale=np.log10(ramanscale)
                    raman=np.log10(raman)
                f=tinythumb(raman,ramanscale,imtimestr,dpi=usedpi)
                thumbname=os.path.join(filelocation,instrumentprefix+'_beta' + file_timetag + '_thumb'+preferredExtension)
                ensureOnlyOne(thumbname,extensionsList)
                f.savefig(thumbname,format=preferredFormat,bbox_inches='tight',dpi=usedpi)
                trimborder(thumbname)
            if hasattr(rs,'linear_depol'):
                ramanscale = [float(disp.get_value('rl_depol_image', 'lo_color_lmt')),
                              float(disp.get_value('rl_depol_image', 'hi_color_lmt'))]
                raman=rs.linear_depol.copy()
                #print np.nanmax(raman),np.nanmin(raman),' is ramans actual range'
                if disp.get_value('rl_depol_image','log_linear')=='log':
                    raman[raman<=0]=np.NAN
                    ramanscale=np.log10(ramanscale)
                    raman=np.log10(raman)
                f=tinythumb(raman,ramanscale,imtimestr,dpi=usedpi)
                thumbname=os.path.join(filelocation,instrumentprefix+'_dep' + file_timetag + '_thumb'+preferredExtension)
                ensureOnlyOne(thumbname,extensionsList)
                f.savefig(thumbname,format=preferredFormat,bbox_inches='tight',dpi=usedpi)
                trimborder(thumbname)

    today=datetime.utcnow()
    if not frappe and not ismf2ship and full24hour and imetime.year==today.year and imetime.month==today.month and imetime.day==today.day:
        destpath=os.getenv("DAYSTORE",os.path.join('/','var','ftp','public_html','hsrl'))
        destfile=os.path.join(destpath,full24hour + '_current'+preferredExtension)
        fulln=hiresfile
        if not fulln:
            return
        outf=file(destfile,'w')
        inf=file(fulln,'r')
        ensureOnlyOne(destfile,extensionsList)
        outf.write(inf.read())
        inf.close()
        outf.close()
        os.unlink(fulln)
    elif ismf2ship:
        if not hiresfile:
            return
        destpath='./'#os.getenv("DAYSTORE",os.path.join('/','var','ftp','public_html','hsrl'))
        destfile=os.path.join(destpath,realend.strftime('hsrl_%Y%m%d%H%M')+preferredExtension)
        fulln=hiresfile
        outf=file(destfile,'w')
        inf=file(fulln,'r')
        ensureOnlyOne(destfile,extensionsList)
        outf.write(inf.read())
        inf.close()
        outf.close()
        os.unlink(fulln)
        recompressImage(destfile,quality=30,optimize=True)
        varlist=('rs_mean.times','rs_mean.latitude','rs_mean.longitude','rs_mean.transmitted_energy','rs_mean.seedvoltage',\
            'rs_mean.seeded_shots','rs_mean.coolant_temperature','rs_mean.laserpowervalues','rs_mean.opticalbenchairpressure',\
            'rs_mean.humidity','rs_mean.l3locking_stats','rs_mean.l3cavityvoltage','rs_mean.nonfiltered_energy','rs_mean.filtered_energy',\
            'rs_mean.builduptime','rs_mean.one_wire_temperatures')
        skip=0
        while skip<16:
           try:
              print 'skip is ',skip
              addExifComment(destfile,makedict(n,varlist,skip=skip))
              break
           except Exception as e:
              skip+=1
              print e
              traceback.print_exc()
              pass
        sendSCPFile(destfile,host='198.129.80.15',path='/ftpdata/outgoing',user='hsrl',keyfile=os.getenv('MF2_HSRL_KEY'))
    elif frappe:
        if not hiresfile:
            return
        sendFTPFile(hiresfile,host='catalog.eol.ucar.edu',path='/pub/incoming/catalog/frappe/',user='anonymous',password='jpgarcia@lidar.ssec.wisc.edu')



def matlabMake(instrument,datetimestart,range_km=15,full24hour=None,*args,**kwargs):
    import lg_dpl_toolbox.core.archival as hru
    env=os.environ
    env['DISPLAY']=':1'
    os.chdir('/var/www/cgi-support/ahsrl_matlab/data_export')
    cmd='./hsrl_calendar.sh'
    datapath=hru.get_path_to_data(instrument,None)
    if not datapath or len(datapath)<5:
        return
    params=[cmd,datapath,datetimestart.strftime('%d-%b-%Y')]
    if not full24hour:
        params.append('am' if datetimestart.hour<12 else 'pm')
        params.append('%i' % range_km)
    bk=multiprocessing.Process(target=os.execve,args=(cmd,params,env))
    bk.start()
    bk.join()
    if full24hour:
        destpath=os.getenv("DAYSTORE",os.path.join('/','var','ftp','public_html','hsrl'))
        destfile=os.path.join(destpath,full24hour + '_current'+preferredExtension)
        imagefolder=os.path.join(datapath,datetimestart.strftime('%Y/%m/%d/images'))
        fulln=None
        for f in os.listdir(imagefolder):
            if f.endswith('%2i.jpg' % range_km):
                fulln=os.path.join(imagefolder,f)
        if not fulln:
            return
        outf=file(destfile,'w')
        inf=file(fulln,'r')
        outf.write(inf.read())
        inf.close()
        outf.close()
    #os.execv('/usr/bin/reset',('usr/bin/reset',))

def storeout():
    #print 'storing out'
    return dict(out=os.dup(sys.stdout.fileno()),err=os.dup(sys.stderr.fileno()))

def putbackout(rec):
    if 'out' not in rec:
        return
    sys.stdout.flush()
    sys.stderr.flush()
    os.dup2(rec['out'],sys.stdout.fileno())
    os.dup2(rec['err'],sys.stderr.fileno())
    #print 'got out back'
    os.close(rec['out'])
    os.close(rec['err'])
    #print 'del out and err',rec['out'],rec['err']
    del rec['out']
    del rec['err']

def wrapperfunc(wrappedfunc,*args,**kwargs):
    rec=storeout()
    wrappedfunc(*args,**kwargs)
    putbackout(rec)

def bkroundMake(func,*args,**kwargs):
    g24h=('full24hour' in kwargs and kwargs['full24hour']!=None)
    gfrap=('frappe' in kwargs and kwargs['frappe'])
    pref=''
    if gfrap:
        pref='frap_'
    elif g24h:
        pref='24h_'
    if kwargs['datetimestart']==kwargs['nowtime']:
        dt=pref+'current'
    else:
        fmt=pref+'%Y.%m.%d' + ('' if g24h else '.%H')
        dt=kwargs['datetimestart'].strftime(fmt)
    fn=kwargs['instrument'] + '_' + dt + '.log' #kwargs['datetimestart'].strftime('_%Y.%m.%d.%H.log')
    if os.getenv('NOLOG',None)!=None:
        mystdout=file('/dev/null','w')
    else:
        mystdout=file(fn,'w')
    if mystdout!=None:
        os.dup2(mystdout.fileno(),sys.stdout.fileno())
        os.dup2(mystdout.fileno(),sys.stderr.fileno())
    print kwargs
    func(*args,**kwargs)

def addContinuity(thetimes,endtime,full24hour=False):
    step=timedelta(days=.5 if not full24hour else 1.0)
    cd=thetimes[-1]+step
    while cd<endtime:
        thetimes.append(cd)
        cd+=step

def ptimestring(ds,prior=None):
    pc=ds.count('.')
    fmt=[]
    repl={}
    if ds.find('.')==2:
        if prior==None:
            prior=datetime.utcnow().replace(hour=0,minute=0,second=0,microsecond=0)
        useprior=True
    else:
        useprior=False
        if pc>=2:
            fmt.append('%Y.%m.%d')
            pc-=2
    if pc>=1 or (useprior and pc>=0):
        fmt.append('%H')
        repl['hour']=None
        pc-=1
    if pc>=1 or (useprior and pc>=0):
        fmt.append('%M')
        repl['minute']=None
        pc-=1
    if pc>=1 or (useprior and pc>=0):
        st+='%S'
        repl['second']=None
        pc-=1
    st='.'.join(fmt)
    ret = datetime.strptime(ds,st)
    if useprior:
        for f in ['hour','minute','second']:
            if f in repl:
                repl[f]=getattr(ret,f)
        ret=prior.replace(**repl)
        if ret<=prior:
            ret+=timedelta(days=1)
    return ret

def main():
    if len(sys.argv)<=1:
        print """Archive Quicklooks and Thumbnail Generation
typical uses-
Make new 12-hour images for the current 12-hour period for the BagoHSRL
    makearchiveimages.py bagohsrl
Make new 12-hour images for February 10, 2014 for the GVHSRL
    makearchiveimages.py gvhsrl 2014.02.10
Make new 12-hour images for February 10 to 12, 2014 for the GVHSRL
    makearchiveimages.py gvhsrl 2014.02.10 .. 2014.02.12
Make new 12-hour images for January 10, 2014 to now (UTC) for the TMPKAZR, and current for mf2hsrl
    makearchiveimages.py tmpkazr 2014.01.10 .. mf2hsrl

less used modes:
    times specified YYYY.MM.DD.HH-YYYY.MM.DD.HH will be interpeted as starttime-endtime, non-inclusive on the end
        to make a GVHSRL image for February 10, 2014 from 5:00 to 10:00 (useful for narrow windows of non-constant operation)
            makearchiveimages.py gvhsrl 2014.02.10.05-2014.02.10.10

    instrument tags (parameters):
        ##km - instead of 15km max, use the specified number of km as a max
            example: 'gvhsrl_30km' instead of 'gvhsrl' will create images up to 30km
        att - plot attenuated backscatter instead of backscatter (will result in a different filename)
            ahsrl_att_20km will make ahsrl images for attenuated backscatter to 20km
        IR - plot IR-1064 channels also
            bagohsrl_IR will also make bagohsrl images for IR-1024, including color ratio and backscatter ratio

Useful environment variables
TIMESHIFT: set to seconds offset to shift the 'now' time by.
   typical use: set to -3600 in a crontab set to go off at 1 minute past
   the hour to ensure the prior 12hour period has a complete image

DEBUG environment variables of use
FORCE_FOREGROUND: set to a non-zero length string to have tasks run one
    at a time without forking. If this is not set, tasks will run in the
    background, as many as 6 at a time, with output to a logfile with
    the instrument name and time, or 'current' if no time specified
"""
        sys.exit(0)
    usematlab=[]#'ahsrl','gvhsrl','nshsrl']
    allowMatlab=True#os.getenv('ALLOW_MATLAB',None)
    if os.getenv('FORCE_PYTHON',None):
        usematlab=[]
    interactive=int(os.getenv('INTERACTIVE','1'))!=0
    instorder=[]
    dates={}
    tasks=[]
    maxtasks=int(os.getenv("MAXTASKS","1"))
    nowtime=datetime.utcnow()
    if os.getenv('TIMESHIFT',None):
        nowtime+=timedelta(seconds=int(os.getenv('TIMESHIFT',0)))
    currinst=''
    continuity=False
    for ds in sys.argv[1:]:
        if '.' in ds:
            if '-' in ds:#date range:
                dss=ds.split('-')
                st=ptimestring(dss[0])
                et=ptimestring(dss[1],prior=st)
                dates[currinst].append(dict(datetimestart=st,realstart=st,realend=et))
                continue
            pc=ds.count('.')
            if ds.startswith('..'):
                continuity=True
                continue
            else:
                d=ptimestring(ds)             
            if d>nowtime:
                d=nowtime
            if continuity and len(dates[currinst])>0:
                continuity=False
                addContinuity(dates[currinst],d)
            dates[currinst].append(d)
            if pc==2:
                dates[currinst].append(d+timedelta(days=.5))
        else:
            if continuity and len(dates[currinst])>0:
                continuity=False
                addContinuity(dates[currinst],nowtime)
            currinst=ds
            if currinst not in dates:
                dates[currinst]=[]
                instorder.append(currinst)
    if continuity and len(dates[currinst])>0:
        continuity=False
        addContinuity(dates[currinst],nowtime)
    for k in dates:
        if len(dates[k])==0:
            dates[k].append(nowtime)
    if interactive:
        print dates
    didARun=True#FIXME disabled this featureFalse
    rec=storeout()
    #print instorder
    #print dates
    finalstring=''
    try:
      for inst in instorder:
        parms={'instrument':inst,'nowtime':nowtime}
        if '_' in inst:
            tmp=inst.split('_')
            parms['instrument']=tmp[0]
            for i in range(1,len(tmp)):
                if tmp[i].endswith('km'):
                    parms['range_km']=int(tmp[i][0:(len(tmp[i])-2)])
                elif tmp[i]=='att':
                    parms['attenuated']=True
                elif tmp[i]=='IR':
                    parms['ir1064']=True
                elif tmp[i]=='frappeftp':
                    parms['full24hour']='frappe'
                    parms['range_km']=20
                    parms['frappe']=True
                elif tmp[i]=='mf2ship':
                    parms['ismf2ship']=True
                    parms['full24hour']='mf2ship'
                    parms['range_km']=15
                    parms['completeframe']=True
                else:
                    parms['full24hour']=tmp[i]
        for date in dates[inst]:
            if isinstance(date,datetime):
                myparms=parms.copy()
                myparms['datetimestart']=date
            else:
                myparms=date
                myparms.update(parms)
            while len(tasks)>=maxtasks:
                for t in copy.copy(tasks):
                    if not t.is_alive():
                        if hasattr(t,'exitcode') and t.exitcode!=0 and interactive:
                            finalstring+='Task '+repr(t.pid)+' returned exit code '+repr(t.exitcode)+'\n'
                        tasks.remove(t)
                if len(tasks)>=maxtasks:
                    time.sleep(30)
            if allowMatlab and (myparms['instrument'] in usematlab):
                myparms['func']=matlabMake
            else:
                myparms['func']=makeArchiveImages
            if interactive:
                print 'starting with ', myparms
            if os.getenv('FORCE_FOREGROUND',None):
                try:
                    myparms['func'](**myparms)
                except Exception as e:
                    print 'Exception',e,'occurred'
                    traceback.print_exc()
            else:
                if (interactive or didARun) and os.getenv('NOLOG',None)==None:
                    p=multiprocessing.Process(target=bkroundMake,kwargs=myparms)
                    p.start()
                else:
                    g=myparms.copy()
                    g['wrappedfunc']=bkroundMake
                    p=threading.Thread(target=wrapperfunc,kwargs=g)
                    p.start()
                    time.sleep(5)
                    didARun=True
                if interactive:
                    print 'Started',myparms
                tasks.append(p)#,filename=inst+"_cache.nc")
      if interactive:
        print 'dispatch complete. waiting for finish'
      while len(tasks)>0:
                for t in copy.copy(tasks):
                    if not t.is_alive():
                        if hasattr(t,'exitcode') and t.exitcode!=0 and interactive:
                            finalstring+='Task '+repr(t.pid)+' returned exit code '+repr(t.exitcode)+'\n'
                        tasks.remove(t)
                if len(tasks)>0:
                    time.sleep(5)
 
                                 #raw_input("press a key")
    except:
        putbackout(rec)
        print 'exception occurred'
        traceback.print_exc()
    putbackout(rec)
    if len(finalstring)>0:
        print finalstring
    if interactive:
        print 'exiting'

if __name__ == '__main__':
    main()
    if False:#True:
        outf='test.jpg'
        file(outf,'w').write(file('hsrl_current.jpg').read())
        tx=dict(temperatures=[1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,32,3,4,453,52435,235,324,5234,532,45,2345,345,23,4523,45,2345,234,523,45,3245,234,5234,52,345,2345,234,5234,523,45,23454],pressures=[.23,3,423,43,234,.234234,234,234.234,234,234.234,234,23.423,4,23])
        import json
        txj=json.dumps(tx)
        addExifComment(outf,txj)#,'ExtendedEXIF')
        print 'added ',tx
        dumpExif(outf)
