#!/usr/bin/python
# -*- coding: utf-8 -*-

import datetime
import numpy as np
import dplkit.role.narrator
import dplkit.role.filter
import dplkit.role.librarian
import json
import os
 
import traceback
from time import sleep
import datetime
import lg_dpl_toolbox.filters.substruct as substruct
import lg_dpl_toolbox.filters.time_frame as time_frame
import lg_dpl_toolbox.dpl.TimeSource as TimeSource
import lg_base.core.array_utils as hau
import dplkit.role.decorator
import warnings
import copy
import threading
import hsrl.filters.dpl_filtered_extinction as dfe
from lg_dpl_toolbox.filters.dpl_rolling_window_filter import dpl_rolling_window_filter
from collections import OrderedDict

def time_mod(value,base,zerotime=None):
    if zerotime is None:
        return {'zero':value}
    if base is None:
        return datetime.timedelta(seconds=0)
    z=zerotime['zero']
    while z<value:
        z+=base
    while z>value:
        z-=base
    zerotime['zero']=z
    return value-z

@dplkit.role.decorator.autoprovidenested(nestedclasses=[hau.Time_Z_Group,hau.rs_xfer,dict])
class dpl_hsrl_profile_filter(dplkit.role.filter.aFilter):
    """ DPL HSRL Profiling Filter Object. generally only be created by dpl_hsrl object
        tacks on the profile subframe, accumulated over the entire interval, yielded in each frame as accumulated to that moment
    """
    def __init__(self,hsrlframestream=None,meanhsrlframestream=None,invhsrlframestream=None,multistream=None,subscopename="profiles",useraw=None):
        super(dpl_hsrl_profile_filter,self).__init__(hsrlframestream)#,self.cal_narr) #FIXME replace url with some manner of unique path
        #import dpl_calibration
        #self.dpl_calibration=dpl_calibration
        import hsrl.data_stream.profile_accumulation as mr 
        self.mr=mr
        import hsrl.data_stream.iodine_argon_temperatures as iat
        self.iat=iat
        self.useraw=useraw
        self.subscopename=subscopename
        self.framestream=hsrlframestream
        self.meanhsrlframestream=meanhsrlframestream
        self.onlyFinal=(meanhsrlframestream is not None)
        self.multistream=(multistream or self.meanhsrlframestream!=None)
        self.invframestream=invhsrlframestream
 
    def __repr__(self):
        return 'DPL HSRL Profile Framestream Narrator'

    def updateProfiles(self,profiles,mean,hsrl_constants,calframe={},qc_mask=None):
        try:
            hsrl_Cxx = calframe.pop('rs_Cxx',None)
            hsrl_cal = calframe.pop('rs_cal',None)
            hsrl_sounding = calframe.pop('sounding',None)

            pc=self.hsrl_process_control if hsrl_Cxx is not None else None
            if not hasattr(mean,'msl_altitudes') and hasattr(mean,'molecular_counts'):
                mean=copy.copy(mean)
                setattr(mean,'msl_altitudes',hau.Z_Array(np.arange(0,mean.molecular_counts.shape[1])*hsrl_constants['binwidth']*3e8/2))
            if pc is not None:
                sel_telescope_dir = pc.get_value('averaged_profiles','telescope_pointing') 
            else:
                sel_telescope_dir='all'
            profiles=self.mr.generate_ave_profiles(
                mean,
                qc_mask,
                hsrl_Cxx,hsrl_constants,
                pc,sel_telescope_dir,self.hsrl_corr_adjusts,old_profiles=profiles)
            #compute temperatures if data is available
            if hasattr(profiles,'molecular_i2a_counts') and hasattr(hsrl_Cxx,'Cam_i2a'):
                profiles.i2a_mol_ratio = self.iat.compute_i2a_mol_ratio(
                               profiles
                              ,hsrl_Cxx
                              ,self.hsrl_corr_adjusts
                              ,pc) 
                profiles.i2a_temperatures=self.iat.compute_temperature_from_i2a(
                    self.hsrl_instrument,profiles.i2a_mol_ratio
                    ,hsrl_cal.i2a_temp_table,hsrl_sounding.pressures,self.hsrl_corr_adjusts)
                profiles.i2a_temperatures[profiles.msl_altitudes <= hsrl_constants['lidar_altitude']+200] = np.nan
            if self.multistream:#this is safe only if we are top level
                self.setProvidesUsing(profiles)
        except Exception as e:
            print 'Exception occurred in profiles ',e
            traceback.print_exc()
            raise
            #do nothing
        return profiles


    def process(self):
        """ main read generator
        """

        profiles=None
        if self.multistream:
            for m,inv in itera.multiiter(self.framestream,self.invframestream):
                if m!=None and m.times.size>0:
                    lasttime=m.times[0]
                calframe=cal_source(lasttime).copy()
                meanparms={}
                if self.invframestream is not None:
                    if hasattr(inv,'qc_mask'):
                        meanparms['qc_mask']=inv.qc_mask
                    meanparms['calframe']=calframe
                profiles=self.updateProfiles(profiles,m,hsrl_constants=calframe.pop('rs_constants',None),**meanparms)
                if not self.onlyFinal:
                    yield copy.deepcopy(profiles)
                elif self.providesRunning:
                    print 'Profiles code shortcircuiting to get provides out'
                    self.doingShortCircuit()
                    yield None
            if self.onlyFinal and profiles is not None:
                yield profiles
        else:
            useraw=self.useraw
            for rs in self.framestream:
                #print rs.rs_mean.times.shape
                #if rs!=None:
                try:
                    rsf=None
                    const=None
                    qc_mask=None
                    if rsf is None and useraw in (None,False) and hasattr(rs,'rs_mean'):
                        rsf=rs.rs_mean
                        if hasattr(rs,'calv'):
                            cal=rs.calv
                            const=rs.calv['rs_constants']
                        elif hasattr(rsf,'calv'):
                            cal=rsf.calv
                            const=rsf.calv['rs_constants']
                        else:
                            const=self.hsrl_constants
                            cal=dict(rs_Cxx=self.hsrl_Cxx,
                                rs_cal=self.hsrl_cal,
                                sounding=self.hsrl_sounding)
                        useraw=False
                    if rsf is None and useraw in (None,True) and hasattr(rs,'rs_raw'):
                        rsf=rs.rs_raw
                        if hasattr(rs,'calv'):
                            const=rs.calv['rs_constants']
                        elif hasattr(rsf,'calv'):
                            const=rsf.calv['rs_constants']
                        else:
                            const=self.hsrl_constants
                        cal=dict()
                        useraw=True
                    if not useraw and hasattr(rs,'rs_inv') and hasattr(rs.rs_inv,'qc_mask'):
                        qc_mask=rs.rs_inv.qc_mask
                    if rsf is not None:
                        print self.subscopename+' PROFILES*******************'
                        profiles=self.updateProfiles(profiles,rsf,const,cal.copy(),qc_mask)
                        print self.subscopename+' DONE PROFILES**************'
                    rs=copy.copy(rs)
                    setattr(rs,self.subscopename,copy.deepcopy(profiles))
                except:
                    print 'Exception in profile'
                    traceback.print_exc()
                yield rs

@dplkit.role.decorator.autoprovide(frameclass=hau.Time_Z_Group)
class dpl_hsrl_inv_process_filter(dplkit.role.filter.aFilter):
    def __init__(self,mean_narr):
        super(dpl_hsrl_inv_process_filter,self).__init__(mean_narr)#,self.cal_narr) #FIXME replace url with some manner of unique path
        self.mean_narr=mean_narr

    def process(self):
        import hsrl.data_stream.processing_utilities as pu
        for rs_mean in self.mean_narr: #FIXME get the proper stream of cals
            yield pu.process_inv_step(rs_mean,self.hsrl_Cxx,self.hsrl_sounding,self.hsrl_constants,self.hsrl_corr_adjusts,self.hsrl_process_control)


@dplkit.role.decorator.autoprovidenested(nestedclasses=[hau.Time_Z_Group,hau.rs_xfer,dict])
class dpl_hsrl_inv_process_complex(dplkit.role.filter.aFilter):
    def __init__(self,mean_narr):
        super(dpl_hsrl_inv_process_complex,self).__init__(mean_narr)#,self.cal_narr) #FIXME replace url with some manner of unique path
        self.mean_narr=mean_narr

    def process(self):
        import hsrl.data_stream.processing_utilities as pu
        for rs in self.mean_narr: #FIXME get the proper stream of cals
            if hasattr(rs,'rs_mean'):
                if not hasattr(rs.rs_mean,'molecular_counts'):
                    rs.rs_inv=copy.deepcopy(rs.rs_mean)
                elif hasattr(rs,'calv'):
                    calv=rs.calv
                    rs.rs_inv=pu.process_inv_step(rs.rs_mean,calv['rs_Cxx'],calv['sounding'], calv['rs_constants'],\
                                        self.hsrl_corr_adjusts, self.hsrl_process_control)
                else:
                    raise RuntimeError('linking old code. there should be a calv here')
                    rs.rs_inv=pu.process_inv_step(rs.rs_mean,self.hsrl_Cxx,self.hsrl_sounding,self.hsrl_constants,self.hsrl_corr_adjusts,self.hsrl_process_control)
            yield rs


class dpl_hsrl_strip_calv(dplkit.role.filter.aFilter):
    def __init__(self,mean_narr):
        super(dpl_hsrl_strip_calv,self).__init__(mean_narr)#,self.cal_narr) #FIXME replace url with some manner of unique path
        self.mean_narr=mean_narr
        self._provlock=threading.Lock()
        self._provides=None
        _=self.provides

    def _cleaned(self,pv):
        if 'calv' in pv:
            pv=pv.copy()
            del pv['calv']
        return pv

    @property
    def provides(self):
        if self._provides is None:
            with self._provlock:
                if self._provides is None:
                    self._provides=self._cleaned(self.mean_narr.provides)
        return self._provides

    def process(self):
        for rs in self.mean_narr: #FIXME get the proper stream of cals
            if hasattr(rs,'calv'):
                rs=copy.copy(rs)
                delattr(rs,'calv')
            yield rs


#@dplkit.role.decorator.exposes_attrs_in_chain(['hsrl_instrument','hsrl_constants','hsrl_Cxx','hsrl_process_control','hsrl_corr_adjusts'])
@dplkit.role.decorator.exposes_attrs_of_field('const_narr')
@dplkit.role.decorator.exposes_attrs_in_chain(['hsrl_cal_stream'])
@dplkit.role.decorator.autoprovidenested(nestedclasses=[hau.Time_Z_Group,hau.rs_xfer,dict])
class dpl_raw_hsrl_narr(dplkit.role.narrator.aNarrator):
    """ DPL HSRL Narrator Object. should only be created by dpl_hsrl object

        :param params: parameters dictionary
        :param const_narr: calibration narration framestream object
        :param lib: raw hsrl reading library object
        :param zoo: raw hsrl zookeeper object
        :param max_range_bin: maximum range bin count (default 5000)

        exposed attributes:
        - hsrl_cal_stream (calibration stream, can be used for parallel stream collection)

        exposed field type in chain:
        
        - hsrl.dpl.calibration.dpl_calibration_narr
    """
    #def get_process_control(self):
    #    return self.cal_narr.get_process_control()
    @property
    def hsrl_cal_stream(self):
        return self.const_narr

    def __init__(self,params,const_narr,lib,zoo,max_range_bin=5000,inclusive=None):
        super(dpl_raw_hsrl_narr,self).__init__(None)#,self.cal_narr) #FIXME replace url with some manner of unique path
        #import dpl_calibration
        #self.dpl_calibration=dpl_calibration
        import lg_base.core.open_config as oc
        self.oc=oc
        self.const_narr=const_narr
        self.netcdf_defaults=None
        self.lib=lib
        self.zoo=zoo
        self.params=params
        self.inclusive=inclusive if inclusive is not None else False
        self.max_range_bin=max_range_bin
        self._localreload()

    def __repr__(self):
        return 'DPL RAW HSRL Framestream Narrator (%s)' % (self.params)

    def _localreload(self):
        netcdf_default_file=self.instrument+'_netcdf_defaults.json'
        fd = self.oc.open_config(netcdf_default_file)
        dd = json.load(fd, object_pairs_hook=OrderedDict)            
        self.netcdf_defaults=dd
        fd.close()  

    def timegen(self):
        curr=self.params['intervalTime']
        if self.params['reverse_padding'] is not None:
            curr-=self.params['reverse_padding']
        end=self.params['finalTime']
        step=datetime.timedelta(hours=1)
        while end is None or curr<end:
            n=curr+step
            if end is not None and n>=end:
                n=end
                if self.params['reverse_padding'] is not None:
                    n+=self.params['reverse_padding']
            yield dict(chunk_start_time=curr,chunk_end_time=n)
            curr=n

    def read(self):
        """ main read generator
        """
        import hsrl.data_stream.hsrl_read_utilities as hru
        import hsrl.data_stream.input_translators as it
        import hsrl.data_stream.preprocess_raw as ppr
        import hsrl.data_stream.preprocess_level2 as ppl2
        params=self.params

        intervalTime=None
        intervalEnd=None
        zoo=self.zoo
        #if params['timeres']!=None and params['timeres']<datetime.timedelta(seconds=float(self.cal_narr.hsrl_constants['integration_time'])):
        #    params['timeres']=None #pure native
        end_time_datetime=params['finalTime']
        #timemodoffset=time_mod(params['realStartTime'],params['timeres'])
        #noframe='noframe'

        cdf_to_hsrl = None
        preprocess_ave = None
        instrument=self.hsrl_instrument
        ntime_ave=1
        streamratemult=int(os.getenv('DEBUG_RAW_FRAME_WIDTH','50'))

        for calv in self.const_narr:#self.timegen(): #self.cal_narr:
            if intervalTime is None:
                intervalTime=calv['chunk_start_time']
                if self.inclusive:
                    intervalTime-=datetime.timedelta(seconds=5)
                intervalEnd=intervalTime
            chunk_end_to_use=calv['chunk_end_time']#-time_mod(calv['chunk_end_time'],params['timeres'],timemodoffset)
            rs_constants=calv['rs_constants']
            #print 'old end',calv['chunk_end_time'],'vs new end',chunk_end_to_use,'mod base',params['timeres'],'offset',timemodoffset
            if calv['chunk_end_time']==calv['chunk_start_time'] and end_time_datetime is None:
                if params['block_when_out_of_data']:
                    if 'timeres' not in params or params['timeres'] is None:
                        sleep(rs_constants['integration_time'])
                    else:
                        sleep(params['timeres'].total_seconds())
                else:
                    yield None #this is done to get out of here, and not get stuck in a tight loop
                continue
            if cdf_to_hsrl is None:
                cdf_to_hsrl = it.raw_translator(instrument,rs_constants,self.hsrl_corr_adjusts)
            else:
                cdf_to_hsrl.update_constants(rs_constants)
            while intervalTime<chunk_end_to_use:
                # BEGIN 'init' section that couldn't start without a constants set from calibration
                #initialize the preprocess and average class that operates on raw data
                #after it is read from the netcdf and before the main processing
                if zoo is None:
                    if preprocess_ave is None:
                        #listed in reverse order
                        ntime_ave=1
                        if 'quarter_wave_plate_rotation' in rs_constants and rs_constants['quarter_wave_plate_rotation'] == 'rotating':
                            ntime_ave = 1
                        elif 'timeres' in params and not params['timeres'] is None:
                            integration_time = rs_constants['integration_time']
                            ntime_ave = max(int(0.5*params['timeres'].total_seconds()/integration_time),1)
                        if ntime_ave!=1:
                            preprocess_ave=ppl2.preprocess_level2(instrument,preprocess_ave)
                            preprocess_ave=ppr.time_ave(ntime_ave,preprocess_ave)
                            preprocess_ave=ppr.time_frame(preprocess_ave)
                        else:
                            preprocess_ave=ppr.time_frame(preprocess_ave)
                    from lg_dpl_toolbox.dpl.NetCDFZookeeper import GenericTemplateRemapNetCDFZookeeper 
                    zoo=GenericTemplateRemapNetCDFZookeeper(instrument,self.netcdf_defaults, self.max_range_bin,preprocess_ave)
                #END init section
                if intervalEnd>=chunk_end_to_use:
                    print 'Breaking calibration on endtime. raw ',intervalEnd,chunk_end_to_use,end_time_datetime
                    break
                else:
                    intervalEnd=chunk_end_to_use

                print ' new raw hsrl window is ', intervalTime, ' to ' , intervalEnd

                if True:# requested_times==None or requested_times.shape[0]>0:
                    rs=None
                    try:
                        for rs_raw in hru.fetch_data(#FIXME this should use the dpl objects
                                instrument,
                                intervalTime,
                                intervalEnd,
                                self.max_range_bin,
                                self.netcdf_defaults,
                                cdf_to_hsrl,
                                dpl_librarian=self.lib,dpl_zookeeper=zoo):
                            print 'read in raw frame ',rs_raw
                            if rs_raw is not None and rs_raw.times.size>0:
                                assert(rs_raw.times[-1] is not None)
                                if True:
                                    for rs_raw1 in rs_raw.iterateAllTimes(ntime_ave*streamratemult):
                                        rs=hau.Time_Z_Group()#can_append=False)
                                        setattr(rs,'rs_raw',rs_raw1)
                                        yield rs
                                else:
                                    rs=hau.Time_Z_Group()#can_append=False)
                                    setattr(rs,'rs_raw',rs_raw)
                                    yield rs
                                if rs_raw.times[-1]>intervalTime and rs_raw.times[-1]<intervalEnd:
                                    intervalTime=rs_raw.times[-1]
                                    if hasattr(rs_raw,'delta_t') and rs_raw.delta_t.size>0:
                                        intervalTime=intervalTime+datetime.timedelta(seconds=rs_raw.delta_t[-1])
                                    else:
                                        print 'WARNING HSRL HAS NO DELTA_T'
                                        intervalTime=intervalTime+datetime.timedelta(seconds=.01)
                        intervalTime=intervalEnd
                             #if rs!=None and hasattr(rs,'profiles'):
                            #    delattr(rs,'profiles')

                    except Exception, e:
                        print 'Exception occured in raw reading'
                        print 'Exception = ',e
                        print traceback.format_exc()
                        if isinstance(e,(MemoryError,)):
                            print 'Please Adjust Your Parameters to be more Server-friendly and try again'
                            raise
                        if not isinstance(e,(AttributeError,)):
                            raise

@dplkit.role.decorator.exposes_attrs_of_field('cal_narr')
@dplkit.role.decorator.exposes_attrs_in_chain(['hsrl_cal_stream'])
@dplkit.role.decorator.autoprovidenested(nestedclasses=[hau.Time_Z_Group,hau.rs_xfer,dict])
class dpl_hsrl_narr(dplkit.role.narrator.aNarrator):
    """ DPL HSRL Narrator Object. should only be created by dpl_hsrl object

        :param params: parameters dictionary
        :param cal_narr: calibration framestream narration object
        :param timesource: time axis generation source (could be synthetic or from another stream)
        :param rawsrc: raw data source. if not provided, will create a lot of it here
        :param lib: raw hsrl reading library object only used if rawsrc is not given
        :param zoo: raw hsrl zookeeper object only used if rawsrc is not given

        exposed attributes:
        - hsrl_cal_stream (calibration stream, can be used for parallel stream collection)

        exposed field type in chain:
        
        - hsrl.dpl.calibration.dpl_calibration_narr
    """
    #def get_process_control(self):
    #    return self.cal_narr.get_process_control()

    @property
    def hsrl_cal_stream(self):
        return self.cal_narr

    def __init__(self,params,cal_narr,timesource,rawsrc=None,compute_stats=0):
        super(dpl_hsrl_narr,self).__init__(None,cal_narr) #FIXME replace url with some manner of unique path
        #import dpl_calibration
        #self.dpl_calibration=dpl_calibration
        #self.provides=libr.provides
        self.compute_stats=compute_stats
        self.rawsrc=rawsrc
        self.params=params
        self.timesource=timesource
        self.cal_narr=cal_narr

    def __repr__(self):
        return 'DPL HSRL Framestream Narrator (%s)' % (self.params)

    def read(self):
        """ main read generator
        """
        import hsrl.data_stream.processing_utilities as pu
        params=self.params
        firsttimeever=None
        intervalTime=None
        intervalEnd=None
        rawsrc=iter(self.rawsrc)
        #if params['timeres']!=None and params['timeres']<datetime.timedelta(seconds=float(self.cal_narr.hsrl_constants['integration_time'])):
        #    params['timeres']=None #pure native
        end_time_datetime=params['finalTime']
        #timemodoffset=time_mod(params['realStartTime'],params['timeres'])
        noframe='noframe'
        fullrange=False #if this is true, it will pad the start with any missing times.

        remainder=None
        cdf_to_hsrl = None
        preprocess_ave = None
        requested_times=None
        instrument=self.hsrl_instrument
        intcount=0
        rs_mem = None
        #rs=None
        timesource=TimeSource.CompoundTimeGenerator(self.timesource) if self.timesource is not None else None
        
        for calv in self.cal_narr:
            if intervalTime is None:
                firsttimeever=calv['chunk_start_time']
                intervalTime=calv['chunk_start_time']
                intervalEnd=intervalTime
            chunk_end_to_use=calv['chunk_end_time']#-time_mod(calv['chunk_end_time'],params['timeres'],timemodoffset)
            #print 'old end',calv['chunk_end_time'],'vs new end',chunk_end_to_use,'mod base',params['timeres']
            if calv['chunk_end_time']==calv['chunk_start_time'] and end_time_datetime is None:
                if params['block_when_out_of_data']:
                    if 'timeres' not in params or params['timeres'] is None:
                        sleep(calv['rs_constants']['integration_time'])
                    else:
                        sleep(params['timeres'].total_seconds())
                else:
                    yield None #this is done to get out of here, and not get stuck in a tight loop
                continue
            while intervalTime<chunk_end_to_use:
                integration_time = calv['rs_constants']['integration_time']
                doPresample=True
                #END init section
                if intervalEnd>chunk_end_to_use:
                    print 'Breaking calibration on endtime. proc ',intervalEnd,chunk_end_to_use,end_time_datetime
                    break
                else:
                    intervalEnd=chunk_end_to_use
                #print ' Absolute window is ', actualStartTime, ' to ' , params['finalTime']
                print ' prior window was ', intervalTime, ' to ' , intervalEnd, 'terminating at ',chunk_end_to_use,rs_mem
                if True:#requested_times==None or requested_times.shape[0]>0:
                    try:
                            try:
                                while rawsrc is not None:
                                    if rs_mem is not None and rs_mem.times[0]>=chunk_end_to_use  and (end_time_datetime is None or chunk_end_to_use<end_time_datetime):
                                        break
                                    tmp=rawsrc.next()
                                    if hasattr(tmp,'rs_raw'):
                                        if rs_mem is not None:
                                            rs_mem.append(tmp.rs_raw)
                                        else:
                                            rs_mem=copy.deepcopy(tmp.rs_raw)
                                    if rs_mem is not None and rs_mem.times.shape>0:
                                        break
                                    else:
                                        rs_mem=None
                            except StopIteration:
                                print 'Raw HSRL stream is ended'
                                rawsrc=None
                            if rs_mem is None or rs_mem.times.size==0:
                                rs_mem=None
                            elif rs_mem.times[0]>=chunk_end_to_use and (end_time_datetime is None or chunk_end_to_use<end_time_datetime):
                                print 'HSRL RAW skipping to next cal because of times',intervalTime,chunk_end_to_use,end_time_datetime,rs_mem.times[0]
                                break
                            else:
                                intervalEnd=rs_mem.times[-1]
                            print 'read in raw frame to mean',rs_mem,remainder
                            if rawsrc is None:
                                intervalEnd=chunk_end_to_use

                            print 'trimmed ',rs_mem
                            if timesource is not None:
                                if timesource.isDone:
                                    break
                                useMungedTimes=False #this is in case this code will need to start shifting bins (which assumes resolutions, and implies start and end of intervales, rather than explicitly to avoid overlap or underlap
                                usePrebinnedTimes=True #this goes in the other direction of munged times to say provided times are timebin borders, and the last time is the end of the last, not included, and thus expected to be the first bin on the next window. thats the fully explicit way to describe the bins in code, but standards in describing bins to the user (a single time when the bin spans a range) is not defined yet
                                inclusive=rawsrc is None and (end_time_datetime!=None and intervalEnd>=end_time_datetime)
                                timevals=hau.T_Array(timesource.getBinsFor(starttime=intervalTime,endtime=intervalEnd,inclusive=inclusive))#,inclusive=(end_time_datetime!=None and intervalEnd>=end_time_datetime)))
                                print 'Now %i intervals %s' % (timevals.size-1, "INC" if inclusive else "NOINC"),intervalTime,intervalEnd
                            elif 'timeres' in params and params['timeres'] is not None:
                                tmp=intervalTime
                                useMungedTimes=False #this is in case this code will need to start shifting bins (which assumes resolutions, and implies start and end of intervales, rather than explicitly to avoid overlap or underlap
                                usePrebinnedTimes=True #this goes in the other direction of munged times to say provided times are timebin borders, and the last time is the end of the last, not included, and thus expected to be the first bin on the next window. thats the fully explicit way to describe the bins in code, but standards in describing bins to the user (a single time when the bin spans a range) is not defined yet

                                timevals=[]
                                timevals.append(tmp)
                                while tmp<intervalEnd:# using python datetimes for making the axis is much much more precise than matplotlib floats.
                                        #print tmp, ' = ' , du.date2num(tmp) , ' = ' , (tmp-self.actualStartTime).total_seconds()
                                        tmp+=params['timeres']
                                        timevals.append(tmp)
                                        
                                #intervalEnd=tmp
                                intcount+=len(timevals)
                                if usePrebinnedTimes:
                                    intcount-=1
                                print 'Now %i intervals' % (intcount)
                                timevals=hau.T_Array(timevals)
                            else:

                                print 'Using Native timing'
                                timevals=None

                            print ' new window is ', intervalTime, ' to ' , intervalEnd

                            requested_times=timevals
                           
                            requested_chunk_times= requested_times#requested_times[requested_times >=intervalTime]

                            if requested_chunk_times is not None and len(requested_chunk_times)<2 and rawsrc is not None:
                                #if rawsrc is not None:
                                print "not enough time to process"
                                continue
                            elif rawsrc is None and rs_mem is None and remainder is None:
                                #chunk_end_to_use=intervalTime
                                #continue
                                #print ''
                                break
                 
          
                            rs_chunk,remainder = pu.process_data( instrument, intervalTime, intervalEnd
                                ,params['min_alt'], params['max_alt'], requested_chunk_times
                                , rs_mem, calv['rs_Cxx'], calv['rs_constants'], calv['rs_cal']
                                , None , self.cal_narr.hsrl_corr_adjusts, self.cal_narr.hsrl_process_control
                                , self.compute_stats,remainder=remainder)
                            rs_mem=None
                            if rs_chunk is not None and hasattr(rs_chunk,'rs_mean') and rs_chunk.rs_mean is not None and rs_chunk.rs_mean.times.size==0:
                                rs_chunk=None
                            if rs_chunk is None and rawsrc is None:
                                break
                           #print rs_chunk
                            if rs_chunk is not None and hasattr(rs_chunk,'rs_mean') and rs_chunk.rs_mean is not None and rs_chunk.rs_mean.times.size>0:
                                if fullrange and requested_chunk_times is not None:
                                    v=hau.Time_Z_Group(like=rs_chunk.rs_mean)
                                    v.times=hau.T_Array(requested_chunk_times[requested_chunk_times<rs_chunk.rs_mean.times[0]])
                                    if v.times.size>0:
                                        rs_chunk.rs_mean.prepend(v)
                                rs_chunk.calv=calv

                                yield rs_chunk
                                intervalTime=intervalEnd
                    except Exception, e:
                        print 'Exception occured in update_cal_and_process'
                        print 'Exception = ',e
                        print traceback.format_exc()
                        if isinstance(e,(MemoryError,)):
                            print 'Please Adjust Your Parameters to be more Server-friendly and try again'
                            raise
                        if not isinstance(e,(AttributeError,)):
                            raise
        assert(remainder is None or remainder.times.size==0)
        if fullrange and end_time_datetime is not None and timesource is not None and (not timesource.isDone or requested_times is not None) and firsttimeever!=intervalTime:#either timesource indicates it wasn't run completely or requested times wasn't cleared
            requested_times=hau.T_Array(timesource.getBinsFor(starttime=intervalTime,endtime=end_time_datetime))
            if requested_times is not None and len(requested_times)>1:
                print 'NO DATA to end from ',intervalTime,' to ',end_time_datetime #FIXME IS THIS USED? JPG 20160504 
                print "times to use are ",requested_times[:-1]
                rs= hau.Time_Z_Group()
                rs.rs_mean=hau.Time_Z_Group()
                rs.rs_mean.times=hau.T_Array(requested_times[:-1]).copy() #causes the time axis to be stored, but all others may be implied MISSING
                setattr(rs.rs_mean,'delta_t',hau.T_Array(np.zeros(rs.rs_mean.times.shape)))
                yield rs

@dplkit.role.decorator.exposes_attrs_of_field('cal')
class dpl_hsrl(dplkit.role.librarian.aLibrarian):
    """ HSRL data DPL framestream

    Example: ::

      r = dpl_hsrl(instrument='gvhsrl')
      for data in r(datetime.datetime(2011,8,11,0,0),datetime.datetime(2011,8,15,0,0), timeres_timedelta=datetime.timedelta(seconds=5), min_alt_m=0,max_alt_m=15000,altres_m=50):
            (data is the rs structure from the processing functions)

    :param instrument: hsrl id string (eg. 'ahsrl','gvhsrl','nshsrl','mf2hsrl').
    :param process_control:        process control structure or json filename (contains corrections and process defaults)
    :param data_request:           data request string   
    :param filetype: HSRL raw filetype to use (None for all (default), 'data' for only data files)

    exposes field:

    - cal (dpl_calibration librarian object)

    """

    def __params__(self,params,calvals):# a lot of duplicate set up from rti.py
        import lg_base.core.canvas_info as ci
        canvas_info=ci.load_canvas_info()
        native_res=calvals['binwidth']*3e8/2

        if params['altres'] is None and params['forimage']: #unspecified, load from processing_defaults
            altrange=params['max_alt']-params['min_alt']
            params['altres']=max(native_res,float(altrange)/float(canvas_info['canvas_pixels']['y']))
            #print self.altres
 
        if params['altres'] is None:
            n_range_ave=1.0
        else:
            n_range_ave = np.ceil(params['altres']/native_res) 
        alt_res = native_res * n_range_ave#FIXME this is filling in a potentially None (pure native) value
        params['deriv_altres']=alt_res
        params['deriv_range_ave']=n_range_ave
        #print 'range ave', alt_res, n_range_ave
        return params

    def __init__(self, instrument,process_control=None,data_request=None, filetype=None,*args, **kwargs):
        super(dpl_hsrl,self).__init__(None)
        self.instrument=instrument
        self.data_request=data_request
        from hsrl.dpl.calibration import dpl_calibration
        self.dpl_calibration=dpl_calibration
        import lg_base.core.open_config as oc
        self.oc=oc
        from lg_dpl_toolbox.dpl.NetCDFZookeeper import GenericTemplateRemapNetCDFZookeeper 
        from hsrl.dpl.HSRLLibrarian import HSRLLibrarian
        self.lib=HSRLLibrarian(instrument=instrument,filetype=filetype)
        self.instrument=self.lib.dataprefix
        self.zoo=kwargs.pop('zoo',None)
        self.cal=self.dpl_calibration.dpl_calibration(instrument=self.instrument,process_control=process_control)
        self.netcdf_defaults=None
        self.qa=None
        if instrument in ['bagohsrl']:
            from hsrl.qa.dpl_narrators import dpl_hsrl_qa
            self.qa=dpl_hsrl_qa(instrument)
        self._localreload()
        if len(args):
            print 'Unused dpl_hsrl args = ',args
        if len(kwargs):
            print "Unused dpl_hsrl kwargs = ",kwargs

    def _localreload(self):
        netcdf_default_file=self.instrument+'_netcdf_defaults.json'
        fd = self.oc.open_config(netcdf_default_file)
        dd = json.load(fd, object_pairs_hook=OrderedDict)            
        self.netcdf_defaults=dd
        fd.close()

    def reload(self):
        self.cal.reload()
        self._localreload()

    def __repr__(self):
        return 'DPL HSRL Framestream (instrument="%s",calibration lib=%s)' % (self.instrument,repr(self.cal))

    def search(self, start_time_datetime=None, end_time_datetime=None,reverse_padding=None,timeres_timedelta=None,min_alt_m=None,max_alt_m=None,altres_m=None,window_width_timedelta=None,
        corr_adjusts=None,block_when_out_of_data=True,forimage=True,inclusive=None, mol_norm_alt_m=None,timesource=None,raw_only=False,cal_only=False,with_profiles=True,do_inversion=True,calibration_overrides=None,
        requested_altitudes=None,calsrc=None,constsrc=None,sounding_source=None,*args, **kwargs):
        """
        :param start_time_datetime: start time (optional)
        :type start_time_datetime: datetime.datetime
        :param end_time_datetime: end time (optional) if unspecified, will continue to return frames thru now, unending
        :type end_time_datetime: datetime.datetime
        :param reverse_padding: (optional)in the case of reading up to the current time, this timedelta is always subtracted from the current time to get the most recent moment to read
        :type reverse_padding: datetime.timedelta
        :param timeres_timedelta: (optional) time resolution, or None to optimized for images
        :type timeres_timedelta: datetime.timedelta
        :param min_alt_m: minimum altitude in meters
        :param max_alt_m: maximum altitude in meters
        :param altres_m: (optional) altitude resolution
        :param window_width_timedelta:   used with start or end time (not both) to determine preferred return width. if unspecified, will imply the other unspecified, or from now if neither
        :type window_width_timedelta: datetime.timedelta
        :param corr_adjusts: correction adjustments. if unspecified, will use default   
        :param block_when_out_of_data: (optional) if unspecified or True, will block for 'timeres_timedelta' until trying for more frames. if False, yield None until more data is available. this only effects behavior if end_time_datetime is None/unspecified
        :param forimage: (optional) True (default) will implicitly set *res if not specified to match image configuration. if false, unspecified resolution will result in native resolution
        :param inclusive: if true, will expand window to ensure including any data that intersects the requested times (NOT IMPLEMENTED)
        """
        if len(args):
            print 'Unused dpl_hsrl.search args = ',args
        if len(kwargs):
            print "Unused dpl_hsrl.search kwargs = ",kwargs

        # altitude configuration notes: min and max alt are required. resolution is optional, in which case the process_control will determine pixel count, and thus resolution

        # time configuration notes: there are many modes to this operation, and two layers
        #   calibration layer: possible parameters are start, end, and window width. (implemented by dpl_calibration.parse_timewindow(), which will return in all cases start and width, and end if will terminate)
        #      specification groups possible: start+end, start+window,start,end+window,window (in order of priority)
        #      start and end specify the range of calibrations to stream
        #      if only one is specified, window is used to roll forward or back the one specified to make the other
        #         if window is not specified, it goes from start to now (defining the window)
        #      if neither are specified, start is set to now-window
        #      if end is not specified, it keeps rolling forward without termination. if it is, it WILL go to that point (even in future) and terminate
        #         if start and window are specified, and start+window is past now, start will be ignored
        #   processing layer:
        #      here, it only cares if timeres is not specified.  If it is not specified, it will follows the same steps as calibration to find the preferred window size, and use process_control to determine pixel count, and resolution

        params={}
        #params['windowwidth']=window_width_timedelta
        #params['intervalTime']=start_time_datetime#pytz.UTC.localize(start_time_datetime);
        import lg_base.core.canvas_info as ci
        canvas_info=ci.load_canvas_info()
        if timesource is not None:
            ts=timesource
        elif timeres_timedelta is None and not forimage:#NATIVE
            ts=None
        else:
            ts=TimeSource.TimeGenerator(start_time=start_time_datetime,end_time=end_time_datetime,width=window_width_timedelta,
                time_resolution=timeres_timedelta,time_step_count=None if (not forimage) else float(canvas_info['canvas_pixels']['x']))
        if ts is None:
            tmp=TimeSource.TimeGenerator(start_time=start_time_datetime,end_time=end_time_datetime,width=window_width_timedelta,
                time_resolution=timeres_timedelta,time_step_count=None if (not forimage) else float(canvas_info['canvas_pixels']['x']))
            params['realStartTime']=tmp.start_time
            params['intervalTime']=tmp.start_time
            params['finalTime']=tmp.end_time
            params['timeres']=tmp.time_resolution
        elif hasattr(ts,'start_time'):
            params['realStartTime']=ts.start_time
            params['intervalTime']=ts.start_time
            params['finalTime']=ts.end_time
            params['timeres']=ts.time_resolution if hasattr(ts,'time_resolution') else None
            #params['windowwidth']=ts.
        else:
            params['timeres']=None

        params['reverse_padding']=reverse_padding
        #actualStartTime=params['realStartTime']
        #params['finalTime']=end_time_datetime
        #params['timeres']=timeres_timedelta
        params['forimage']=forimage

        params['min_alt']=min_alt_m
        params['max_alt']=max_alt_m
        params['altres']=altres_m
        params['block_when_out_of_data']=block_when_out_of_data
        instrument=self.instrument
        intcount=0
        #fixme this should deprecate the params structure, making only deliberate parts exposed to the narrator as needed, via the calibration (altitude) and timesource (time), and timeslice/blocking (local)
        #ts=None
        ret=None

        const_narr=calsrc or constsrc
        if const_narr is None:
            const_narr=self.cal(interval_start_time=params['intervalTime'],reverse_padding=params['reverse_padding'],interval_end_time=params['finalTime'],
                            corr_adjusts=corr_adjusts,mol_norm_alt=mol_norm_alt_m,nocal=True)
            if cal_only and raw_only:
                return const_narr
        if not cal_only:
            ret=dpl_raw_hsrl_narr(params=params,const_narr=const_narr,lib=self.lib,zoo=self.zoo,inclusive=inclusive)
        if not raw_only:

            params=self.__params__(params,const_narr.hsrl_constants_first)
            cal_narr=calsrc
            if cal_narr is None:
                cal_narr=self.cal(min_alt_m=params['min_alt'],max_alt_m=params['max_alt'],altres_m=params['deriv_altres'],useconsts=const_narr
                    ,calibration_overrides=calibration_overrides,requested_altitudes=requested_altitudes,soundingsource=sounding_source)
            elif sounding_source is not None or requested_altitudes is not None or calibration_overrides is not None:
                warnings.warn("Not using sounding source, altitudes, or calibration overrides as provided.")
            if cal_only:
                return cal_narr
            ret=dpl_hsrl_narr(params=params,cal_narr=cal_narr,rawsrc=ret,timesource=ts)

            if 'rs_mean' in ret.provides:#mean filter
                windowparms=dfe.mean_filter_setup(ret.hsrl_process_control,ret.hsrl_constants_first,ret.provides['rs_mean'],cal_narr.provides)

                if windowparms is not None:
                    import lg_dpl_toolbox.filters.time_frame as time_slicing
                    import lg_dpl_toolbox.filters.substruct as frame_substruct
                    splitter=frame_substruct.SubstructBrancher(ret)
                    sliced=time_slicing.TimeGinsu(splitter.narrateSubstruct('rs_mean'),timefield='times',dtfield='delta_t',omitTime=False) #break it up
                    subslice=splitter.narrateSubstruct('rs_mean')
                    masterslice=splitter.narrateSubstruct(None)

                    if windowparms is not None:
                        if not isinstance(windowparms,(list,tuple)):
                            windowparms=[windowparms]
                        for w in windowparms:
                            sliced=dpl_rolling_window_filter(sliced,w) #run the filters

                    sliced=frame_substruct.CountDeGinsu(frame_substruct.FrameLength(subslice,'times'),sliced) #re-assemble in same chunk size
                    ret=frame_substruct.NestingCompositer(masterslice,dict(rs_mean=sliced))


            if do_inversion:
                ret=dpl_hsrl_inv_process_complex(ret)

                if 'rs_inv' in ret.provides:
                  windowparms=dfe.inv_filter_setup(ret.hsrl_process_control,ret.hsrl_constants_first,ret.provides['rs_inv'],ret.provides['rs_mean'])
                  qa=self.qa
                  if not ret.hsrl_process_control.enabled('quality_assurance',return_if_missing=True): #omitted or explicitly not enabled
                        qa=None

                  if (windowparms is not None or qa is not None):
                    import lg_dpl_toolbox.filters.time_frame as time_slicing
                    import lg_dpl_toolbox.filters.substruct as frame_substruct
                    splitter=frame_substruct.SubstructBrancher(ret)
                    sliced=time_slicing.TimeGinsu(splitter.narrateSubstruct('rs_inv'),timefield='times',dtfield='delta_t',omitTime=False) #break it up
                    subslice=splitter.narrateSubstruct('rs_inv')
                    masterslice=splitter.narrateSubstruct(None)

                    if windowparms is not None:
                        if not isinstance(windowparms,(list,tuple)):
                            windowparms=[windowparms]
                        for w in windowparms:
                            sliced=dpl_rolling_window_filter(sliced,w) #run the filters

                    if qa is not None:
                        sliced=qa(start_time_datetime= params['realStartTime'], end_time_datetime=params['finalTime'],hostsource=sliced,hostsource_newframe='qa_flags')

                    sliced=frame_substruct.CountDeGinsu(frame_substruct.FrameLength(subslice,'times'),sliced) #re-assemble in same chunk size
                    ret=frame_substruct.NestingCompositer(masterslice,dict(rs_inv=sliced))

        if with_profiles:
            ret=dpl_hsrl_profile_filter(ret,subscopename='raw_profiles',useraw=True)
            if not raw_only:
                ret=dpl_hsrl_profile_filter(ret)
        ret=dpl_hsrl_strip_calv(ret)
        return ret

if __name__ == '__main__':
    dplhsrl=dpl_hsrl('bagohsrl',datetime.timedelta(seconds=60*60*2))
    for i in dplhsrl(start_time_datetime=datetime.datetime(2013,1,10,0,0,0),end_time_datetime=datetime.datetime(2013,1,10,12,0,0),min_alt_m=0,max_alt_m=15000):
        print 'step'
