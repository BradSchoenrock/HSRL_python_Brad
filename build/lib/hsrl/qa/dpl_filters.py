import dplkit.role.filter

import lg_dpl_toolbox.filters.time_frame as tf
#this will contain filter objects
import lg_base.core.array_utils as hau
from datetime import datetime,timedelta
import numpy as np
import copy

@dplkit.role.decorator.exposes_attrs_of_field('hostsource')
class QAFlagClonedAttachFilter(dplkit.role.filter.aFilter):
    def __init__(self,qasource,hostsource,hostsource_newframe,altname=None,constantAltitude=None,splitFields=True):
        super(QAFlagClonedAttachFilter,self).__init__(qasource)
        self.qasource=qasource
        self.qaparser=qasource.flagmanager
        self.hostsource=hostsource
        self.hostsource_newframe=hostsource_newframe.split('.') if hostsource_newframe else []
        self.timename='start'
        self.dtname='width'
        self.altname=altname
        self.provides=copy.deepcopy(hostsource.provides)
        self.provides['times']=dict(shortname='times',type=hau.T_Array)
        self.provides['delta_t']=dict(shortname='delta_t',type=hau.T_Array)
        pv=self.provides
        if len(self.hostsource_newframe)>0:
            for f in self.hostsource_newframe[:-1]:
                pv=pv[f]
            if self.hostsource_newframe[-1] not in pv:
                pv[self.hostsource_newframe[-1]]={}
            pv=pv[self.hostsource_newframe[-1]]
        if len(self.hostsource_newframe)>0:
            pv['start']=dict(shortname='start',type=datetime)
            pv['width']=dict(shortname='width',type=timedelta)
            pv['times']=dict(shortname='times',type=hau.T_Array)
            pv['delta_t']=dict(shortname='delta_t',type=hau.T_Array)
            pv['altitudes']=dict(shortname='altitudes',type=hau.Z_Array)
        self.splitFields=splitFields
        if self.splitFields:
            for f in self.qaparser.flagbits.keys():
                pv['qa_'+f]=dict(shortname='qa_'+f,type=hau.TZ_Array)
        else:
            pv['qaflags']=dict(shortname='qaflags',type=hau.TZ_Array)
        self.constantAltitude=constantAltitude
        if self.altname is not None and not (self.altname in hostsource.provides):
            self.altname=None

    def process(self):
        fr=None
        flags=None
        olda=None
        qasource=tf.TimeTrickle(self.qasource,'time')
        altitudes=hau.Z_Array(self.qaparser.altitudeAxis)
        for _f in self.hostsource:
            #FIXME include angles
            f=_f
            if not isinstance(f,dict):
                f=vars(f)
            t=f[self.timename]
            a=self.constantAltitude
            if fr is None or ((not qasource.atEnd) and t>=qasource.nextTime):#if need an update to the qa record
                #print 'Getting qa source for time',t
                fr=qasource(t)
                flags=None
            if 'range_flags' in fr and fr['range_flags'] is not None:#if there is a range dependence, and a potentially non-constant altitude
                if self.altname is not None and self.altname in f:
                    a=f[self.altname]
                if a is None:
                    raise RuntimeError('Need platform altitude to merge in range-dependant qa Flags')                    
                if olda is None or a!=olda:
                    flags=None
                    olda=a

            if flags is None:#was cleared either because new flags from the qa file, or new altitude from the stream
                if 'range_flags' in fr and fr['range_flags'] is not None:
                    flags=self.qaparser.mergeVectors(fr['flags'],fr['range_flags'],a)
                else:
                    flags=fr['flags']
                flags=self.qaparser.translateToEnumeration(flags)
                flags=hau.TZ_Array(flags.reshape([1]+list(flags.shape)),dtype='int32',summode='and')
            if len(self.hostsource_newframe)>0:
                sfn=self.hostsource_newframe[-1]
                sf=f
                for xf in self.hostsource_newframe[:-1]:
                    if xf not in sf:
                        sf=None
                        break
                    sf=sf[xf]
                    if not isinstance(sf,dict):
                        sf=vars(sf)
                if sf is None:
                    yield _f
                    continue
                ret=hau.Time_Z_Group(timevarname='times',altname='altitudes')
                setattr(ret,'start',f['start'])
                setattr(ret,'width',f['width'])
                setattr(ret,'times',hau.T_Array([f['start']]))
                setattr(ret,'delta_t',hau.T_Array([f['width'].total_seconds()]))
                setattr(ret,'altitudes',copy.copy(altitudes))
                sf[sfn]=ret
                ret=vars(ret)
            else:
                ret=f
            if self.splitFields:
                for f,idx in self.qaparser.flagbits.items():
                    ret['qa_'+f]=hau.TZ_Array((flags/(10**idx))%10,dtype='int32',summode='and')
            else:
                ret['qaflags']=flags
            yield _f


@dplkit.role.decorator.exposes_attrs_of_field('timealtsource')
class QAFlagClonedSyncFilter(dplkit.role.filter.aFilter):
    def __init__(self,qasource,timealtsource,altname=None,constantAltitude=None,splitFields=True):
        super(QAFlagClonedSyncFilter,self).__init__(qasource)
        self.qasource=qasource
        self.qaparser=qasource.flagmanager
        self.timealtsource=timealtsource
        self.timename='start'
        self.dtname='width'
        self.altname=altname
        self.provides={}
        self.provides['times']=dict(shortname='times',type=hau.T_Array)
        self.provides['delta_t']=dict(shortname='delta_t',type=hau.T_Array)
        self.provides['start']=dict(shortname='start',type=datetime)
        self.provides['width']=dict(shortname='width',type=timedelta)
        self.provides['altitudes']=dict(shortname='altitudes',type=hau.Z_Array)
        self.splitFields=splitFields
        if self.splitFields:
            for f in self.qaparser.flagbits.keys():
                self.provides['qa_'+f]=dict(shortname='qa_'+f,type=hau.TZ_Array)
        else:
            self.provides['qaflags']=dict(shortname='qaflags',type=hau.TZ_Array)
        self.constantAltitude=constantAltitude
        if self.altname is not None and not (self.altname in timealtsource.provides):
            self.altname=None

    def process(self):
        fr=None
        flags=None
        olda=None
        qasource=tf.TimeTrickle(self.qasource,'time')
        altitudes=hau.Z_Array(self.qaparser.altitudeAxis)
        for f in self.timealtsource:
            #FIXME include angles
            if not isinstance(f,dict):
                f=vars(f)
            t=f[self.timename]
            a=self.constantAltitude
            if fr is None or ((not qasource.atEnd) and t>=qasource.nextTime):#if need an update to the qa record
                #print 'Getting qa source for time',t
                fr=qasource(t)
                flags=None
            if 'range_flags' in fr and fr['range_flags'] is not None:#if there is a range dependence, and a potentially non-constant altitude
                if self.altname is not None and self.altname in f:
                    a=f[self.altname]
                if a is None:
                    raise RuntimeError('Need platform altitude to merge in range-dependant qa Flags')                    
                if olda is None or a!=olda:
                    flags=None
                    olda=a
            if flags is None:#was cleared either because new flags from the qa file, or new altitude from the stream
                if 'range_flags' in fr and fr['range_flags'] is not None:
                    flags=self.qaparser.mergeVectors(fr['flags'],fr['range_flags'],a)
                else:
                    flags=fr['flags']
                flags=self.qaparser.translateToEnumeration(flags)
                flags=hau.TZ_Array(flags.reshape([1]+list(flags.shape)),dtype='int32',summode='and')
            ret=hau.Time_Z_Group(timevarname='times',altname='altitudes')
            setattr(ret,'times',hau.T_Array([t]))
            setattr(ret,'delta_t',hau.T_Array([f['width'].total_seconds()]))
            setattr(ret,'altitudes',copy.copy(altitudes))
            setattr(ret,'start',f['start'])
            setattr(ret,'width',f['width'])
            if self.splitFields:
                for f,idx in self.qaparser.flagbits.items():
                    setattr(ret,'qa_'+f,hau.TZ_Array((flags/(10**idx))%10,dtype='int32',summode='and'))
            else:
                setattr(ret,'qaflags',flags)
            yield ret

