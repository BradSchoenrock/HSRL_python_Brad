#!/usr/bin/env python
import datetime
from bottleneck import nanmin,nanmax,nanmean,nanmedian
import numpy
import os
import matplotlib
matplotlib.use('Agg')

def singleVar(content):
    return dict(min=nanmin(content),max=nanmax(content),mean=nanmean(content),median=nanmedian(content),
                valid=numpy.sum(numpy.isfinite(content))*100.0/content.size)

def expl(buf,prefixes=[]):
    ret=dict()
    for i,v in buf.items():
        p=prefixes+[i]
        if isinstance(v,dict):
            ret.update(expl(v,p))
        elif isinstance(v,basestring):
            ret['_'.join(p)]=v
        else:
            ret['_'.join(p)]='%g' % v
    return ret

def writeHierarchy(fil,names,buf):
    buf=expl(buf)
    lin=[]
    if len(names)==0:
        names.extend([k for k in buf.keys()])
        names.sort()
        idxx=0
        for v in ('start','end'):
            if v in names:
                names.remove(v)
                names.insert(idxx,v)
                idxx+=1
        fil.write('#' + (", ".join(names)) + '\n')
    for n in names:
        if n in buf:
            lin.append(buf[n])
        else:
            lin.append('0')
    fil.write((", ".join(lin)) + '\n')
    fil.flush()

def logData(fil,cols,inv):
    print vars(inv).keys()
    buf=dict(start=inv.times[0].strftime('%Y%m%dT%H%M%S'),
                 end=(inv.times[0]+datetime.timedelta(seconds=inv.delta_t[0])).strftime('%Y%m%dT%H%M%S'))
    buf['sc_ratio_par'] =singleVar( inv.Na[0, :] / inv.Nm_i2[0, :])
    buf['sc_ratio_perp']=singleVar(inv.Ncp[0, :] / inv.Nm_i2[0, :])
    if hasattr(inv,'Na_i2a'):
        buf["i2a_par"] =singleVar(inv.Na_i2a[0,:] / inv.Nm_i2a[0,:])
        buf["i2a_perp"]=singleVar(inv.Ncp[0,:]    / inv.Nm_i2a[0,:])
    print buf
    writeHierarchy(fil,cols,buf)

def skimcals(inst,start,end):
    from time import sleep
    from hsrl.dpl.dpl_hsrl import dpl_hsrl
    from hsrl.dpl.calibration.dpl_calibration import dpl_singlecalibration_narr
    import hsrl.dpl.dpl_artists as artists
    import lg_dpl_toolbox.filters.time_frame as time_slicing
    import lg_dpl_toolbox.filters.substruct as frame_substruct
    import hsrl.graphics.hsrl_display as du
    dropcontent=['rs_raw','rs_mean']
    storage='cal_verify'
    try:
        os.makedirs(storage)
    except OSError:
        pass
    getfigs=('sc_ratio_profile','raw_profiles','dark_corrected_profiles','corrected_profiles','dif_geo_profiles',\
        'wfov_geo_adjust','geometry_correction','calibration_coefficients','lapse_rate','backscat_profile','sounding')
    preferredFormat='jpg'
    preferredExtension='.'+preferredFormat
    withImage=True
    display_defaults='all_plots.json'
    print 'display_defaults=',display_defaults
    #time.sleep(5)
    dplobj=dpl_hsrl(instrument=inst)
    minalt_km=0
    maxalt_km=20
    altres_m=45.0

    cals=dplobj.cal(interval_start_time=start,interval_end_time=end,min_alt_m=minalt_km*1000.0,max_alt_m=maxalt_km*1000.0,altres_m=altres_m)

    filename='cal_log_'+inst+start.strftime('_%Y%m%dT%H%M%S')+end.strftime('_%Y%m%dT%H%M%S')+'.txt'
    outf=file(os.path.join(storage,filename),"w")
    cols=[]

    for calv in cals:
        s=calv['chunk_start_time']
        e=calv['chunk_end_time']
        if (e-s).total_seconds()<300:
            continue
        mycalgen=dpl_singlecalibration_narr(inst,calv,cals.hsrl_process_control,cals.corr_adjusts)
        dplgen=dplobj(start_time_datetime=s,end_time_datetime=e,min_alt_m=minalt_km*1000.0,max_alt_m=maxalt_km*1000.0,
            with_profiles=True,do_inversion=False,calsrc=mycalgen)
        dplgen=frame_substruct.DropFrameContent(dplgen,dropcontent)
        dplgen=time_slicing.FrameCachedConcatenate(dplgen)
        if withImage:
            (disp,conf)=du.get_display_defaults(display_defaults)
            for f in disp.get_attrs():
                disp.set_value(f,'enable',0)
            for f in getfigs:
                disp.set_value(f,'enable',1)
            dplgen=artists.dpl_images_artist(framestream=dplgen,instrument=inst,max_alt=maxalt_km*1000.0,
                processing_defaults=dplgen.hsrl_process_control,display_defaults=disp)
        haveData=False
        for rs in dplgen:
            if hasattr(rs,'profiles') and hasattr(rs.profiles,'inv'):
                logData(outf,cols,rs.profiles.inv)
                haveData=True
        if withImage and haveData:
            figs=dplgen.figs
            datetag=s.strftime('%Y%m%dT%H%M%S')+e.strftime('_%Y%m%dT%H%M%S')
            for f in getfigs:
                if f in figs:
                    print 'getting fig',f
                    fig=figs.figure(f)
                    hiresfile=os.path.join(storage,'%s_%s_%s%s' %(inst,f,datetag,preferredExtension))
                    print hiresfile
                    fig.savefig(hiresfile,format=preferredFormat,bbox_inches='tight')
        #print 'sleeping'
        #sleep(10)



def main():
    import sys
    p=os.path.realpath(os.path.join(os.path.dirname(sys.argv[0]),'..','..'))
    print(p)
    sys.path.append(p)
    inst=sys.argv[1]
    if len(sys.argv)>2:
        startdate=datetime.datetime.strptime(sys.argv[2],"%Y%m%dT%H%M%S")
    else:
        startdate=datetime.datetime.utcnow()
        startdate=startdate.replace(hour=0,minute=0,second=0,microsecond=0)
    if len(sys.argv)>3:
        enddate=datetime.datetime.strptime(sys.argv[3],"%Y%m%dT%H%M%S")
    else:
        enddate=datetime.datetime.utcnow()

    skimcals(inst,startdate,enddate)

if __name__ == '__main__':
    main()