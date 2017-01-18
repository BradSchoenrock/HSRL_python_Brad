import dplkit.role.librarian
import dplkit.role.narrator
import dplkit.role.filter
import lg_base.formats.calvals as cru
from datetime import datetime,timedelta
import os
import re
import time
import numpy
import lg_dpl_toolbox.filters.substruct as substruct
import lg_base.core.array_utils as hau
import dplkit.role.decorator
import numpy as np
import copy
from bottleneck import nanmax
import lg_dpl_toolbox.filters.time_frame as time_frame
import traceback
import raman.core.raman_processing_utilities as rif
from lg_base.formats.vector_table import calibration_vector_table
import lg_base.core.read_utilities as hru
import lg_base.core.git_tools as git_tools
import lg_dpl_toolbox.dpl.calvals_narrator as cvn

getRollingDescription=rif.getRollingDescription

@dplkit.role.decorator.exposes_attrs_in_chain(['raman_constants_first','raman_platform','raman_process_control','timeinfo'])
class dpl_constants_narr(cvn.dpl_calvals_narrator):
    """ Raman Constants Framestream Narrator. should only be created from the dpl_calibration object

        :param instrument: instrument name
        :param timeinfo: time window info dictionary
        :param calvals: calvals object
        :param process_control: processing control json structure

        exposed attributes:

        - raman_constants_first (calvals constants from start of window)
        - raman_platform (raman instrument platform name)
        - raman_process_control (static raman processing control json structure)
    """

    @property
    def raman_platform(self):
        return self.platform
    @property
    def raman_process_control(self):
        return self.rprocess_control
    @property
    def raman_corr_adjusts(self):
        return self.rcorr_adjusts

    @property
    def raman_constants_first(self):
        return self.constants_first

    def __init__(self,platform,timeinfo,process_control,corr_adjusts=None):
        super(dpl_constants_narr, self).__init__(timeinfo,cru.cal_file_reader(platform))                
        self.platform=platform
        self.rprocess_control=process_control
        if self.rprocess_control is None:
            self.rprocess_control='raman_process_control.json'
        if isinstance(self.rprocess_control,basestring):
            import lg_base.core.locate_file as lf
            import lg_base.core.json_config as jc
            self.rprocess_control=jc.json_config(lf.locate_file(self.rprocess_control,systemOnly=False),'raman_process_defaults')
        self.rcorr_adjusts=self.rprocess_control.get_dict('corr_adjusts').copy()
        if corr_adjusts!=None:
            self.rcorr_adjusts.update(corr_adjusts)

@dplkit.role.decorator.exposes_attrs_of_field('consts')
@dplkit.role.decorator.exposes_attrs_in_chain(['raman_sourcepath'])
@dplkit.role.decorator.autoprovidenested(nestedclasses=[hau.Time_Z_Group,hau.rs_xfer,dict,rif.cal_vector_object,calibration_vector_table])
class dpl_calibration_tables_narr(dplkit.role.narrator.aNarrator):
    """ HSRL Calibration Framestream Narrator. should only be created from the dpl_calibration object

        :param consts: calvals stream
        :param max_range_bin: maximum range bin
        :param requested_altitudes: requested altitudes vector
        :param alternate_cal_dir: alternate location to get calibration tables

        exposed attributes:
  
        - hsrl_cal (calibration tables structure)

        exposes field:

        - consts (dpl_constants_narr narrator object)
            - hsrl_constants (calvals constants, regularly updated)
            - hsrl_instrument (hsrl instrument name)
            - hsrl_process_control (static hsrl processing control json structure)
            - hsrl_corr_adjusts (static correction adjustments)
    """
  
    def __init__(self,basedir,consts,alternate_cal_dir=None):
        self.sourcepath=basedir
        self.consts=consts
        self.alternate_cal_dir=alternate_cal_dir

    @property
    def raman_sourcepath(self):
        return self.sourcepath

    def read(self):
        """generator function.
        """
        if True:
            alternate_cal_dir=self.alternate_cal_dir
            if alternate_cal_dir==None:
                alternate_cal_dir = self.raman_process_control.get_value('alternate_cal_dir','full_dir_path') 
                if alternate_cal_dir == 'None':
                     alternate_cal_dir = None
            if alternate_cal_dir==None:
                alternate_cal_dir=os.getenv('RAMAN_ALTERNATE_CAL_DIR',None)
            exposeBasedir=alternate_cal_dir
            #if exposeBasedir!=None:
            #    if hasattr(self,'hsrl_instrument'):
            #        exposeBasedir=os.path.join(exposeBasedir,self.hsrl_instrument)#fixme
            #    else:
            #        exposeBasedir=os.path.join(exposeBasedir,self.raman_platform)#fixme

        rs_cal = None
        #from hsrl.dpl.HSRLLibrarian import HSRLLibrarian

        #basedir = HSRLLibrarian(instrument=self.instrument).basedir
        #import lg_base.core.git_tools as git_tools
        chunk_end_time=None
        chunk_start_time=None
        #print 'beginning const iteration using ',self.consts
        for constframe in self.consts:#interval_end_time==None or chunk_start_time<interval_end_time:
          #print 'const frame is ',constframe
          rs_constants=constframe['rs_constants']
          nextconsttime=constframe['chunk_end_time']
          if chunk_start_time==None:
              chunk_start_time=constframe['chunk_start_time']
          while chunk_end_time==None or chunk_end_time<nextconsttime:
            rs_cal=rif.load_raman_tables(self.raman_sourcepath,self.raman_platform,rs_constants,chunk_start_time,alternate_cal_dir=exposeBasedir)
            chunk_end_time= min([nextconsttime, rs_cal.expire_time])
            if chunk_end_time<=chunk_start_time:
                print 'WARNING raman dpl_table_calibration trying to use 0-length window. source is behind',nextconsttime, rs_cal.expire_time,chunk_end_time,chunk_start_time
                chunk_end_time=nextconsttime
                if chunk_end_time<=chunk_start_time:
                    break
            retdict=constframe.copy()
            retdict.update(dict(chunk_start_time=chunk_start_time,
                                chunk_end_time=chunk_end_time,
                                rs_cal=rs_cal))
            try:
                    gitver,gitdate=git_tools.getCurrentHash(os.path.join(exposeBasedir or self.raman_sourcepath,'%04i' % chunk_start_time.year,'%02i' % chunk_start_time.month))
                    if gitver!=None: #none is returned if not supported by python
                        #retdict['gitversion']=gitver
                        retdict['raman_cal_gitversion']=gitver
            except ValueError:
                pass#not managed, no version
            if hasattr(rs_cal,'geo') and rs_cal.geo!=None and hasattr(rs_cal.geo,'data') and rs_cal.geo.data is not None:
                retdict['geo_corr']=rs_cal.geo.data[:,1]
            yield retdict
            chunk_start_time=chunk_end_time
  
        return

@dplkit.role.decorator.exposes_attrs_of_field('tables')
@dplkit.role.decorator.exposes_attrs_in_chain(['altitudeAxis'])
@dplkit.role.decorator.autoprovidenested(nestedclasses=[hau.Time_Z_Group,hau.rs_xfer,dict,rif.cal_vector_object,calibration_vector_table])
class dpl_calibration_narr(dplkit.role.narrator.aNarrator):
    """ HSRL Calibration Framestream Narrator. should only be created from the dpl_calibration object

        :param tables: caltables stream
        :param requested_altitudes: requested altitudes vector
        :param sounding_source: optional sounding source. will use configured from calvals if not provided at init

        exposed attributes:
  
        - altitudeAxis
        - hsrl_sounding (current sounding)
        - hsrl_Cxx (calibration structure)

        exposes field:

        - tables (dpl_calibration_tables_narr narrator object)
    """

    @property
    def altitudeAxis(self):
        return self.requested_altitudes
  
    def __init__(self,tables,requested_altitudes=None,sounding_source=None):
        self.tables=tables
        self.requested_altitudes=requested_altitudes
        self.sounding_source=sounding_source
        if sounding_source is None and requested_altitudes is None:
            raise RuntimeError('Need a sounding source or an altitude axis')
        assert(not (requested_altitudes is not None and sounding_source is not None))
        if self.sounding_source is None:
            self.sounding_source=self.make_default_sounding()


    def make_default_sounding(self):
        import hsrl.dpl.calibration.dpl_dynamic_atmospheric_profile_narrator as ddapn
        return ddapn.dpl_dynamic_hsrl_atmospheric_profile_narrator(self.raman_platform,self.tables.timeinfo,self.requested_altitudes,edgepadding=timedelta(days=5),
            calvals=cru.cal_file_reader(self.raman_platform),soundingdatapath=self.raman_sourcepath)
        import atmospheric_profiles.dpl.dpl_temperature_profiles as dtp#import dpl_soundingarchive,dpl_virtualradiosonde
        interval_start_time=self.tables.timeinfo['starttime']
        interval_end_time=self.tables.timeinfo['endtime']
        rs_constants=self.raman_constants_first
        requested_altitudes=None if not hasattr(self,'requested_altitudes') else self.requested_altitudes
        if rs_constants['sounding_type'] in ('virtual','remote_virtual','model','virtual_remote','virtual_cache'):
            raise NotImplementedError(rs_constants['sounding_type'])
            from hsrl.dpl.HSRLLibrarian import HSRLLibrarian
            from lg_dpl_toolbox.dpl.NetCDFZookeeper import GenericTemplateRemapNetCDFZookeeper
            #import lg_dpl_toolbox.filters.time_frame as time_frame# TimeGinsu
            tempzoo=GenericTemplateRemapNetCDFZookeeper(self.instrument,user_read_mode='position')#,keepfields=fields)
            templib=HSRLLibrarian(instrument=self.instrument,zoo=tempzoo)#site=16)#None,datatype)
            m=templib(start=interval_start_time,end=interval_end_time)#,filetype='data')
            delt=15#90
            if 'virtual_sounding_update_interval' in rs_constants:
                delt=rs_constants['virtual_sounding_update_interval']
            elif 'sounding_update_interval' in rs_constants:
                delt=rs_constants['sounding_update_interval']
            delt=timedelta(minutes=delt)
            gribcache=os.getenv('GRIB_CACHE','/arcueid/data/grib_cache')
            try:
                os.makedirs(gribcache)
            except:
                pass
            if not os.path.exists(gribcache):
                raise RuntimeError("Can't Access GRIB Cache directory %s. If this isn't an expected error, set environment variable \"GRIB_CACHE\" to something better and try again." % gribcache)
            parameters=dict(requested_altitudes=requested_altitudes \
                ,remote=False ,download=True)
            extrasearchparams=dict()
            if 'remote' in rs_constants['sounding_type']:
                parameters['remote']=True
            if rs_constants['sounding_type'] in ('model','virtual_cache'):
                parameters['download']=False
            if rs_constants['sounding_type'] == 'model':
                parameters['format']='model.grib2'
                parameters['predict_horizon']=24*365
                if 'model_filename' in rs_constants:
                    parameters['format']=rs_constants['model_filename']
                print 'VRS Forecast model will use grib2 file at '+os.path.join(gribcache,parameters['format'])
            if 'static_cal_time' in rs_constants and rs_constants['static_cal_time'] is not None and len(rs_constants['static_cal_time'])>2:
                extrasearchparams['fixed_time']=hru.convert_date_str(rs_constants['static_cal_time'])['datetime']
            if 'sounding_horizon' in rs_constants and rs_constants['sounding_horizon'] is not None:
                if 'predict_horizon' not in parameters or rs_constants['sounding_horizon']>parameters['predict_horizon']:
                    parameters['predict_horizon']=rs_constants['sounding_horizon']
            if 'predict_horizon' in parameters:
                    print "Predict horizon is ",parameters['predict_horizon'],'hours'
                    interval_start_time=interval_start_time-timedelta(hours=parameters['predict_horizon'])
            soundinglib=dtp.dpl_virtualradiosonde(rs_constants['sounding_type'],gribcache,delt,**parameters)
            if m.provides!=None and 'latitude' in m.provides:
                print 'Got latitude stream'
                m=time_frame.TimeGinsu(m,'times',None)
                soundingnarr=soundinglib(m,starttime=interval_start_time,**extrasearchparams)
            else:
                print 'no position stream/ using config'
                soundingnarr=soundinglib(starttime=interval_start_time,endtime=interval_end_time,fixed_position=(rs_constants['latitude'],rs_constants['longitude']),**extrasearchparams)
            ret=soundingnarr
        elif rs_constants['sounding_type'] in ('sparc','ssec'):
            import atmospheric_profiles.dpl.sparc_profiles as sparc#import dpl_soundingarchive,dpl_virtualradiosonde
            soundinglib=sparc.SPARCSondeLibrarian(self.raman_sourcepath,requested_altitudes)
            ret=soundinglib(interval_start_time,interval_end_time)
        elif rs_constants['sounding_type'] in ('arm',):
            import atmospheric_profiles.dpl.arm_profiles as armp#import dpl_soundingarchive,dpl_virtualradiosonde
            soundinglib=armp.ARMSondeLibrarian(self.raman_sourcepath,requested_altitudes)
            ret=soundinglib(interval_start_time,interval_end_time)
        else:
            expire_duration=None
            if 'sounding_update_interval' in rs_constants:
                expire_duration=timedelta(minutes=rs_constants['sounding_update_interval'])
            elif 'installation' not in rs_constants or rs_constants['installation']=='ground':
                expire_duration=timedelta(hours=1)
            elif rs_constants['installation']=='airborne': #these are the defaults based on platform
                expire_duration=timedelta(minutes=5)
            elif rs_constants['installation']=='shipborne':
                expire_duration=timedelta(minutes=60)
            else:
                raise RuntimeError('Installation of '+rs_constants['installation']+' in calvals is unknown')
            soundinglib=dtp.dpl_soundingarchive(self.raman_sourcepath,
                rs_constants['sounding_type'],
                rs_constants['sounding_id'],
                requested_altitudes,expire_duration=expire_duration)
            ret=soundinglib(interval_start_time,interval_end_time)
        if hasattr(self,'requested_altitudes'):
            import lg_dpl_toolbox.filters.resample_altitude as altitude_resampling
            ret=altitude_resampling.ResampleXd(ret,'altitudes',self.requested_altitudes)
        return ret

    def read(self):
        """generator function.
        """

        soundingsource=time_frame.TimeTrickle(self.sounding_source,timename='sample_time',getClosest=True)

        chunk_start_time=None

        chunk_end_time=None
        sounding=None
        priorsounding=None
        priorconstants=None
        priorcal=None
        #print 'beginning const iteration using ',self.consts
        for calframe in self.tables:
          #print 'const frame is ',constframe
          rs_constants=calframe['rs_constants']
          rs_cal=calframe['rs_cal']
          nextcaltime=calframe['chunk_end_time']
          if chunk_start_time==None:
            chunk_start_time=calframe['chunk_start_time']
          for sounding,chunk_end_time in soundingsource.trickleGenerator(chunk_start_time,nextcaltime):
            pos=''
            if 'installation' in rs_constants and rs_constants['installation']!='ground':
                pos=' ; %fN %fE' % (sounding.latitude,sounding.longitude)
            print 'sounding= "%s" [ %s%s ]'  % (sounding.station_id[:],sounding.times.strftime('%d-%b-%y %H:%M'),pos)
            if chunk_end_time<=chunk_start_time:
                chunk_end_time=nextcaltime
            if chunk_end_time==chunk_start_time:
                print 'Ignored. chunk is length',(chunk_end_time-chunk_start_time)
                print 'times are',chunk_start_time,chunk_end_time,nextcaltime
                break
            try:
                if priorsounding is not sounding or priorconstants is not rs_constants or priorcal is not rs_cal:
                    rs_Cxx = rif.raman_quick_cal(rs_constants,rs_cal,sounding)
                    priorsounding=sounding
                    priorconstants=rs_constants
                    priorcal=rs_cal
            except RuntimeError:
                if chunk_start_time>datetime.utcnow():
                    print 'Repeating last calibration. Requested interval extends beyond now (and known calibrations)'
                else:
                    print 'calibration error for an interval in non-future space.'
                    raise #reraise this.
            if chunk_end_time<=chunk_start_time:
                print 'WARNING dpl_calibration_narr trying to use 0-length window. source is behind',chunk_end_time,chunk_start_time,nextcaltime
                chunk_end_time=nextcaltime
                if chunk_end_time<=chunk_start_time:
                    break
                #chunk_end_time=use_end_time
                #break
            retdict= calframe.copy()
            retdict.update(dict(chunk_start_time=chunk_start_time,
                                chunk_end_time=chunk_end_time,
                                sounding=sounding,
                                Cxx=rs_Cxx))
            print 'Yielding calframe ',chunk_start_time,chunk_end_time,sounding.times
            yield retdict
            chunk_start_time=chunk_end_time
  
        return

@dplkit.role.decorator.exposes_attrs_of_field('raman_calibration_tables_stream')
@dplkit.role.decorator.exposes_attrs_in_chain(['raman_calibration_tables_stream'])
@dplkit.role.decorator.autoprovidenested(nestedclasses=[hau.Time_Z_Group,hau.rs_xfer,dict])
class RamanRangedFilter(dplkit.role.filter.aFilter):
    """ Raman Range-based processing
    """
    def __init__(self,datasource,process_control=None,calsource=None,datasourcescope=None,destinationscope=None):
        super(RamanRangedFilter,self).__init__(datasource)
        self.cals=calsource
        self.datasource=datasource
        if self.cals is None:
            if hasattr(self.datasource,'raman_constants_stream'):
                self.cals=self.datasource.raman_constants_stream
            else:
                self.cals=dpl_constants_narr(datasource.platform,datasource.timeinfo,process_control)
            self.cals=dpl_calibration_tables_narr(datasource.arm_sourcepath,self.cals)
        else:
            assert(process_control is None)
        import cPickle as pickle
        self.rif=rif
        self.pickle=pickle
        self.sourcename=datasourcescope
        self.destination=destinationscope

    @property
    def raman_calibration_tables_stream(self):
        return self.cals

    def process(self):
        #calsource=time_frame.TimeTrickle(self.cals,timename='chunk_start_time',endtimename='chunk_end_time')
        rif=self.rif
        pickle=self.pickle
        lasttime=None
        calsource=time_frame.TimeTrickle(self.cals,timename='chunk_start_time',endtimename='chunk_end_time')
        for f in self.datasource:
            orig=f
            if self.sourcename is not None:
                if isinstance(f,dict) and self.sourcename in f:
                    f=f[self.sourcename]
                elif hasattr(f,self.sourcename):
                    f=getattr(f,self.sourcename)
                else:
                    print "No Raman Ranged to process. can't find subframe "+self.sourcename
                    continue
            if f is None or f.times.size==0:
                print "No Raman Ranged to process. Raman Merge is empty"
                continue
            print "Raman Merge to Ranged running"
            lasttime=f.times[0]
            caldictionary=calsource(lasttime)
            pname=('frame','cal frame')
            p1=[]
            p1.append(pickle.dumps(f))
            p1.append(pickle.dumps(caldictionary))
            ret=rif.process_ranged_raman(caldictionary['rs_constants'],f,self.raman_process_control,caldictionary['rs_cal'],self.raman_corr_adjusts)
            p2=[]
            p2.append(pickle.dumps(f))
            p2.append(pickle.dumps(caldictionary))
            for i in range(len(p1)):
                if p1[i]!=p2[i]:
                    print 'Detected modification of '+pname[i]+'. FAIL'
                    raise RuntimeError('Modification of Raman Ranged '+pname[i])

            hau.verifyAllNew(ramancal=caldictionary,ramanmerge=f,ramanrange=ret)

            if ret is None:
                print "Raman Range is none"
                import time
                time.sleep(5)
            elif hasattr(f,'alt'):#platform height from sealevel
                setattr(ret,'_altitudevarname','altitudes')
                setattr(ret,"altitudes",hau.Z_Array(f.heights+f.alt))

            if self.destination is not None:
                if isinstance(orig,dict):
                    orig[self.destination]=ret
                else:
                    setattr(orig,self.destination,ret)
                yield orig
            else: 
                yield ret


@dplkit.role.decorator.exposes_attrs_of_field('raman_calibration_stream')
@dplkit.role.decorator.exposes_attrs_in_chain(['raman_calibration_stream'])
@dplkit.role.decorator.autoprovidenested(nestedclasses=[hau.Time_Z_Group,hau.rs_xfer,dict])
class RamanInvertFilter(dplkit.role.filter.aFilter):
    """ Raman Inverting
    """
    def __init__(self,datasource,process_control=None,calsource=None,datasourcescope=None,destinationscope=None,soundingsource=None,requested_altitudes=None):
        super(RamanInvertFilter,self).__init__(datasource)
        self.cals=calsource
        self.datasource=datasource
        if self.cals is None:
            if hasattr(self.datasource,'raman_calibration_tables_stream'):
                self.cals=self.datasource.raman_calibration_tables_stream
            else:
                if hasattr(self.datasource,'raman_constants_stream'):
                    self.cals=self.datasource.raman_constants_stream
                    assert(process_control is None)
                else:
                    self.cals=dpl_constants_narr(datasource.platform,datasource.timeinfo,process_control)
                self.cals=dpl_calibration_tables_narr(datasource.arm_sourcepath,self.cals)
            if soundingsource is None and requested_altitudes is None:
                requested_altitudes=datasource.altitudeAxis
            self.cals=dpl_calibration_narr(self.cals,requested_altitudes=requested_altitudes,sounding_source=soundingsource)
        else:
            assert(process_control is None)
            assert(requested_altitudes is None)
            assert(soundingsource is None)
        import raman.core.raman_processing_utilities as rif
        import cPickle as pickle
        self.rif=rif
        self.pickle=pickle
        self.sourcename=datasourcescope
        self.destination=destinationscope

    @property
    def raman_calibration_stream(self):
        return self.cals

    def process(self):
        #calsource=time_frame.TimeTrickle(self.cals,timename='chunk_start_time',endtimename='chunk_end_time')
        rif=self.rif
        pickle=self.pickle
        lasttime=None
        calsource=time_frame.TimeTrickle(self.cals,timename='chunk_start_time',endtimename='chunk_end_time')
        for orig in self.datasource:
          
            if self.sourcename is None:
                highf=orig
            else:
                if isinstance(orig,dict) and self.sourcename in orig:
                    highf=orig[self.sourcename]
                elif hasattr(orig,self.sourcename):
                    highf=getattr(orig,self.sourcename)
                else:
                    print "No Raman Inv to process. can't find subframe "+self.sourcename
                    highf=None
            if highf is not None and highf.times.size==0:
                highf=None
            if highf is None:
                print "No Raman Inv to process. Raman Merge is empty"
                continue
            print "Raman Merge to INV running"
            lasttime=highf.times[0]
            caldictionary=calsource(lasttime)
            pname=('frame','cal frame')
            p1=[]
            p1.append(pickle.dumps(highf))
            p1.append(pickle.dumps(caldictionary))
            ret=rif.process_raman(caldictionary['rs_constants'],highf,self.raman_process_control,caldictionary['rs_cal']
                ,caldictionary['Cxx'],self.raman_corr_adjusts)
            p2=[]
            p2.append(pickle.dumps(highf))
            p2.append(pickle.dumps(caldictionary))
            for i in range(len(p1)):
                if p1[i]!=p2[i]:
                    print 'Detected modification of '+pname[i]+'. FAIL'
                    raise RuntimeError('Modification of Raman '+pname[i])

            #hau.verifyAllNew(ramancal=caldictionary,sounding=sounding,ramanmergehigh=highf,ramanmerge=lowf,ramaninv=ret)

            if ret is None:
                print "Raman INV is none"
                import time
                time.sleep(5)

            if self.destination is not None:
                if isinstance(orig,dict):
                    orig[self.destination]=ret
                else:
                    setattr(orig,self.destination,ret)
                yield orig
            else: 
                yield ret


@dplkit.role.decorator.autoprovidenested(nestedclasses=[hau.Time_Z_Group,hau.rs_xfer,dict])
class dpl_raman_profile_filter(dplkit.role.filter.aFilter):
    """ DPL HSRL Profiling Filter Object. generally only be created by dpl_hsrl object
        tacks on the profile subframe, accumulated over the entire interval, yielded in each frame as accumulated to that moment
    """
    def __init__(self,framestream,srcscopename=None,subscopename=None,calsource=None):
        super(dpl_raman_profile_filter,self).__init__(framestream)#,self.cal_narr) #FIXME replace url with some manner of unique path
        import raman.core.raman_profiles as rp
        self.rp=rp
        self.subscopename=subscopename
        self.srcscopename=srcscopename
        self.framestream=framestream
        self.calsource=calsource
        if self.calsource is None:
            if hasattr(self,'raman_calibration_stream'):
                self.calsource=framestream.raman_calibration_stream
            elif hasattr(self,'raman_calibration_tables_stream'):
                self.calsource=framestream.raman_calibration_tables_stream
            elif hasattr(self,'raman_constants_stream'):
                self.calsource=framestream.raman_constants_stream
        self.onlyFinal=(subscopename is None)
 
    def __repr__(self):
        return 'DPL Raman Profile Framestream Narrator'

    def updateProfiles(self,profiles,mean,constants,calframe={},qc_mask=None):
        try:
            Cxx = calframe.pop('Cxx',None)
            rs_cal = calframe.pop('rs_cal',None)
            raman_sounding = calframe.pop('sounding',None)
            process_control= None if not hasattr(self,'raman_process_control') else self.raman_process_control
            corr_adjusts= None if not hasattr(self,'raman_corr_adjusts') else self.raman_corr_adjusts

            pc=None#self.raman_process_control if hsrl_Cxx is not None else None
            profiles=self.rp.accumulate_raman_profiles(constants,mean,qc_mask,Cxx=Cxx,rs_cal=rs_cal
                ,process_control=process_control,corr_adjusts=corr_adjusts,old_profiles=profiles)
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
        calsource = None if self.calsource is None else time_frame.TimeTrickle(self.calsource,timename='chunk_start_time',endtimename='chunk_end_time')
        for f in self.framestream:
            if f is None:
                continue
            origf=f
            if self.srcscopename is not None:
                print 'Updating profile for scope ',self.srcscopename
                if hasattr(f,self.srcscopename):
                    f=getattr(f,self.srcscopename)
                else:
                    f=None
            if f is None or f.times.size==0:
                continue
            lasttime=f.times[0]
            profparams={}
            if calsource is not None:
                calframe=calsource(lasttime).copy()
            else:
                calframe=dict()#rs_constants=self.raman_constants_first)
            profparams['calframe']=calframe
            profiles=self.updateProfiles(profiles,f,constants=calframe.pop('rs_constants',None),**profparams)
            if not self.onlyFinal:
                ret=copy.deepcopy(profiles)
                if self.subscopename is not None:
                    tmp=copy.copy(origf)
                    setattr(tmp,self.subscopename,ret)
                    ret=tmp
                yield ret
            elif self.providesRunning:
                if profiles is not None:#this is safe only if we are top level
                    self.setProvidesUsing(profiles)
                    print 'Profiles code shortcircuiting to get provides out'
                    self.doingShortCircuit()
                    yield None
        if self.onlyFinal and profiles is not None:
            yield profiles


@dplkit.role.decorator.autoprovidenested(nestedclasses=[hau.Time_Z_Group,hau.rs_xfer,dict])
class dpl_raman_inverted_profile_filter(dplkit.role.filter.aFilter):
    """ DPL HSRL Profiling Filter Object. generally only be created by dpl_hsrl object
        tacks on the profile subframe, accumulated over the entire interval, yielded in each frame as accumulated to that moment
    """
    def __init__(self,framestream,srcscopename=None,subscopename=None,constsource=None):
        super(dpl_raman_inverted_profile_filter,self).__init__(framestream)#,self.cal_narr) #FIXME replace url with some manner of unique path
        import raman.core.raman_profiles as rp
        self.rp=rp
        self.subscopename=subscopename
        self.srcscopename=srcscopename
        self.constsource=constsource
        if self.constsource is None:
            if hasattr(self,'raman_constants_stream'):
                self.constsource=framestream.raman_constants_stream
        self.framestream=framestream
        self.onlyFinal=(subscopename is None)
 
    def __repr__(self):
        return 'DPL Raman Profile Framestream Narrator'

    def updateProfiles(self,profiles,mean,constants,qc_mask=None):
        try:
            process_control= None if not hasattr(self,'raman_process_control') else self.raman_process_control
            corr_adjusts= None if not hasattr(self,'raman_corr_adjusts') else self.raman_corr_adjusts

            pc=None#self.raman_process_control if hsrl_Cxx is not None else None
            profiles=self.rp.accumulate_raman_inverted_profiles(constants,mean,qc_mask
                ,process_control=process_control,corr_adjusts=corr_adjusts,old_profiles=profiles)
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
        constsource = None if self.constsource is None else time_frame.TimeTrickle(self.constsource,timename='chunk_start_time',endtimename='chunk_end_time')
        for f in self.framestream:
            if f is None:
                continue
            origf=f
            if self.srcscopename is not None:
                print 'Updating profile for scope ',self.srcscopename
                if hasattr(f,self.srcscopename):
                    f=getattr(f,self.srcscopename)
                else:
                    f=None
            if f is None or f.times.size==0:
                continue
            lasttime=f.times[0]
            profparams={}
            if constsource is not None:
                calframe=constsource(lasttime).copy()
            else:
                calframe=dict()#rs_constants=self.raman_constants_first)
            profiles=self.updateProfiles(profiles,f,constants=calframe.pop('rs_constants',None),**profparams)
            if not self.onlyFinal:
                ret=copy.deepcopy(profiles)
                if self.subscopename is not None:
                    tmp=copy.copy(origf)
                    setattr(tmp,self.subscopename,ret)
                    ret=tmp
                yield ret
            elif self.providesRunning:
                if profiles is not None:#this is safe only if we are top level
                    self.setProvidesUsing(profiles)
                    print 'Profiles code shortcircuiting to get provides out'
                    self.doingShortCircuit()
                    yield None
        if self.onlyFinal and profiles is not None:
            yield profiles

