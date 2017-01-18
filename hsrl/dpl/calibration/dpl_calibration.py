
import dplkit.role.narrator
import dplkit.role.librarian
import dplkit.role.filter
import hsrl.calibration.calibration_utilities as cu
import numpy as np
import lg_base.core.open_config as oc
from lg_base.core.locate_file import locate_file
import json
from datetime import datetime,timedelta
import hsrl.calibration.cal_vectors as mr 
import lg_base.core.json_config as jc
import lg_dpl_toolbox.filters.time_frame as time_frame
import dplkit.role.decorator
import os
import lg_base.core.array_utils as hau
import hsrl.data_stream.hsrl_read_utilities as hru
from collections import OrderedDict
import warnings
import copy
import time

import lg_dpl_toolbox.dpl.calvals_narrator as cvn

@dplkit.role.decorator.exposes_attrs_in_chain(['hsrl_constants','hsrl_constants_first','hsrl_Cxx','hsrl_instrument','hsrl_process_control',\
    'hsrl_corr_adjusts','timeinfo','hsrl_cal','altitudeAxis','hsrl_sounding'])
@dplkit.role.decorator.autoprovidenested(nestedclasses=[hau.Time_Z_Group,hau.rs_xfer,dict,mr.cal_vectors,hru.calibration_vector_table],reuseGenerator=False)
class dpl_singlecalibration_narr(dplkit.role.narrator.aNarrator):

    def __init__(self,instrument,calframe,process_control,corr_adjusts):
        self.instrument=instrument
        self.calframe=calframe
        self.process_control=process_control
        self.corr_adjusts=corr_adjusts
        self._timeinfo=time_frame.parse_timewindow(self.calframe['chunk_start_time'],self.calframe['chunk_end_time'])

    @property
    def hsrl_constants(self):
        warnings.warn("retrieving HSRL_CONSTANTS attribute is deprecated. It changes, so you should be using a data stream directly or indirectly, or use HSRL_CONSTANTS_FIRST instead. Detailed plan pending")
        return self.hsrl_constants_first

    @property
    def hsrl_constants_first(self):
        return self.calframe['rs_constants']

    @property
    def hsrl_Cxx(self):
        warnings.warn("retrieving HSRL_CXX attribute is deprecated. It changes, so you should be using a data stream directly or indirectly. Detailed plan pending")
        return self.calframe['rs_Cxx'] if 'rs_Cxx' in self.calframe else None

    @property
    def hsrl_instrument(self):
        return self.instrument

    @property
    def hsrl_process_control(self):
        return self.process_control

    @property
    def hsrl_corr_adjusts(self):
        return self.corr_adjusts

    @property
    def timeinfo(self):
        return self._timeinfo

    @property
    def hsrl_cal(self):
        warnings.warn("retrieving HSRL_CAL attribute is deprecated. It changes, so you should be using a data stream directly or indirectly. Detailed plan pending")
        if 'rs_cal' not in self.calframe:
            raise AttributeError('rs_cal')
        return self.calframe['rs_cal'] if 'rs_cal' in self.calframe else None

    @property
    def altitudeAxis(self):
        if 'sounding' not in self.calframe:
            raise AttributeError('altitudeAxis')
        return self.calframe['sounding'].altitudes if 'sounding' in self.calframe else None

    @property
    def hsrl_sounding(self):
        warnings.warn("retrieving HSRL_SOUNDING attribute is deprecated. It changes, so you should be using a data stream directly or indirectly. Detailed plan pending")
        if 'sounding' not in self.calframe:
            raise AttributeError('sounding')
        return self.calframe['sounding'] if 'sounding' in self.calframe else None

    def read(self):
        yield self.calframe

@dplkit.role.decorator.exposes_attrs_in_chain(['hsrl_constants','hsrl_constants_first','hsrl_Cxx','hsrl_instrument','hsrl_process_control','hsrl_corr_adjusts','timeinfo'])
class dpl_constants_narr(cvn.dpl_calvals_narrator):
    """ HSRL Constants Framestream Narrator. should only be created from the dpl_calibration object

        :param instrument: instrument name
        :param timeinfo: time window info dictionary
        :param calvals: calvals object
        :param process_control: processing control json structure
        :param corr_adjusts: correction adjustments (optional)
        :param mol_norm_alt: molecular normalziation altitude (optional)

        exposed attributes:

        - hsrl_constants_first (first time's calvals constants)
        - hsrl_Cxx (set to None)
        - hsrl_instrument (hsrl instrument name)
        - hsrl_process_control (static hsrl processing control json structure)
        - hsrl_corr_adjusts (static correction adjustments)
    """

    @property
    def hsrl_instrument(self):
        return self.instrument
    @property
    def hsrl_process_control(self):
        return self.process_control
    @property
    def hsrl_corr_adjusts(self):
        return self.corr_adjusts
    @property
    def hsrl_constants(self):
        warnings.warn("retrieving HSRL_CONSTANTS attribute is deprecated. It changes, so you should be using a data stream directly or indirectly. Detailed plan pending")
        return self.rs_constants
    @property
    def hsrl_constants_first(self):
        return self.constants_first

    @property
    def hsrl_Cxx(self):
        warnings.warn("retrieving HSRL_CXX attribute is deprecated. It changes, so you should be using a data stream directly or indirectly. Detailed plan pending")
        return None 

    def __init__(self,instrument,timeinfo,calvals,process_control,corr_adjusts=None,mol_norm_alt=None):
        super(dpl_constants_narr, self).__init__(timeinfo,calvals)                
        self.instrument=instrument
        self.process_control=process_control
        if self.process_control==None:
            self.corr_adjusts=corr_adjusts
        else:
            self.process_control=self.process_control.copy()
            self.corr_adjusts=process_control.get_dict('corr_adjusts')
            if corr_adjusts!=None:
                self.corr_adjusts.update(corr_adjusts)
        self.rs_constants=self.calvals.select_time(self.timeinfo['starttime'])
        if mol_norm_alt!=None:
            assert(self.process_control!=None)
            self.process_control.set_value('mol_norm_alt','meters',mol_norm_alt)

@dplkit.role.decorator.exposes_attrs_in_chain(['hsrl_cal'])
class dpl_override_calibration_tables_narr(dplkit.role.filter.aFilter):
    """ HSRL Override Calibration Framestream filter. should only be created from the dpl_calibration object

        :param cal: cal stream
        :param max_range_bin: maximum range bin
        :param alternate_cal_dir: alternate location to get calibration tables

        other parmaeters are key -> filename to overload calibrations. Will error on init if the keys are unknown

        exposed attributes:
  
        - hsrl_cal (calibration tables structure)
    """
  
    def __init__(self,fs,max_range_bin,alternate_cal_dir=None,**kwargs):
        super(dpl_override_calibration_tables_narr, self).__init__(fs)        
        self.fs=fs
        self.overrides=dict()
        keymap=dict(baseline=hru.read_baseline,
            geo=hru.read_geo_corr,
            n_geo=hru.read_nadir_geo_corr,
            qw_baseline=hru.read_qw_baseline,
            diff_geo=hru.read_diff_geo,
            cp_d_geo=hru.read_cross_poll_diff_geo,
            i2a_diff_geo=hru.read_i2a_diff_geo,
            i2scan=hru.read_i2_scan,
            i2a_temp_table=hru.read_i2a_temp_table,
            pol_cal=hru.read_pol_cal)

        for k,v in kwargs.items():
            if k not in self.provides['rs_cal']:
                raise RuntimeError('Calibration table '+k+" isn't in the calibration stream already. Typo?")
            if k not in keymap:
                raise RuntimeError('Calibration table '+k+" not known to dpl_override_calibration_narr. update the code with a function reference")
            if isinstance(v,datetime): #date given
                t=keymap[k](self.hsrl_instrument,v,alternate_cal_dir,max_range_bin=max_range_bin)
            elif isinstance(v,basestring): #filename given
                t=keymap[k](self.hsrl_instrument,None,None,max_range_bin=max_range_bin,filename=v)
            elif isinstance(v,hru.calibration_vector_table):
                t=v
            else:
                raise RuntimeError('Unknown type for override on calibration table '+k)
            t._expire_time=datetime(2100,1,1,0,0,0) #make sure even time-read ones never expire
            self.overrides[k]=t
 
    @property
    def hsrl_cal(self):
        warnings.warn("retrieving HSRL_CAL attribute is deprecated. It changes, so you should be using a data stream directly or indirectly. Detailed plan pending")
        return self.cal

    def process(self):
        """generator function.
        """

        geo_corr=None
        if 'geo' in self.overrides:
            geo=self.overrides['geo']
            if hasattr(geo,'data') and geo.data is not None:
                geo_corr=geo.data[:,1]

        for calframe in self.fs:#interval_end_time==None or chunk_start_time<interval_end_time:
            #print 'const frame is ',constframe
            dupe=copy.copy(calframe)
            dupe['rs_cal']=copy.copy(dupe['rs_cal'])
            self.cal=dupe['rs_cal']
            for k,v in self.overrides.items():
                setattr(self.cal,k,v)
            if geo_corr is not None:
                dupe['geo_corr']=geo_corr
            yield dupe
  

@dplkit.role.decorator.exposes_attrs_of_field('consts')
@dplkit.role.decorator.exposes_attrs_in_chain(['hsrl_cal'])
@dplkit.role.decorator.autoprovidenested(nestedclasses=[hau.Time_Z_Group,hau.rs_xfer,dict,mr.cal_vectors,hru.calibration_vector_table])
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
  
    def __init__(self,consts,max_range_bin,alternate_cal_dir=None):
        self.consts=consts
        self.rs_constants=consts.hsrl_constants_first
        #self.timeinfo=consts.timeinfo
        self.max_range_bin=max_range_bin
        self.alternate_cal_dir=alternate_cal_dir
        self.cal=None

    @property
    def process_control(self):
        return self.hsrl_process_control
    @property
    def instrument(self):
        return self.hsrl_instrument
    @property
    def hsrl_cal(self):
        warnings.warn("Warning: retrieving HSRL_CAL attribute is deprecated. It changes, so you should be using a data stream directly or indirectly. Detailed plan pending")
        return self.cal

    def read(self):
        """generator function.
        """
        alternate_cal_dir=self.alternate_cal_dir
        if alternate_cal_dir==None:
            alternate_cal_dir = self.process_control.get_value('alternate_cal_dir','full_dir_path') 
            if alternate_cal_dir == 'None':
                 alternate_cal_dir = None
        if alternate_cal_dir==None:
            alternate_cal_dir=os.getenv('HSRL_ALTERNATE_CAL_DIR',None)
        exposeBasedir=alternate_cal_dir
        if exposeBasedir!=None:
            exposeBasedir=os.path.join(exposeBasedir,self.hsrl_instrument)

        rs_cal = None
        from hsrl.dpl.HSRLLibrarian import HSRLLibrarian

        basedir = HSRLLibrarian(instrument=self.instrument).basedir
        import lg_base.core.git_tools as git_tools
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
            try:
                if rs_cal==None:
                    rs_cal=mr.cal_vectors(self.instrument,chunk_start_time,self.max_range_bin,alternate_cal_dir)
                    print 'New cal vectors read at ', chunk_start_time
                elif rs_cal.expire_time <= chunk_start_time:
                    rs_cal.read(chunk_start_time)
                    print 'New cal vectors read at ', chunk_start_time
            except RuntimeError:
                if chunk_start_time>datetime.utcnow():
                    print 'Repeating last calibration. Requested interval extends beyond now (and known calibrations)'
                else:
                    print 'calibration error for an interval in non-future space.'
                    raise #reraise this.
            chunk_end_time= min([nextconsttime, rs_cal.expire_time])
            if chunk_end_time<=chunk_start_time:
                warnings.warn("dpl_table_calibration trying to use 0-length window. source is behind")
                print "behind times",nextconsttime, rs_cal.expire_time,chunk_end_time,chunk_start_time
                chunk_end_time=nextconsttime
                if chunk_end_time<=chunk_start_time:
                    break
                #chunk_end_time=use_end_time
            #sounding=soundingsource(chunk_start_time)
            #setattr(self,'hsrl_constants',rs_constants)#FIXME exposing these in this way, when they are dynamic parts of the framestream is Bad Bad Bad
            setattr(self,'cal',rs_cal)
            retdict=constframe.copy()
            retdict.update(dict(chunk_start_time=chunk_start_time,
                                chunk_end_time=chunk_end_time,
                                rs_cal=rs_cal))
            try:
                    gitver,gitdate=git_tools.getCurrentHash(os.path.join(exposeBasedir or basedir,'%04i' % chunk_start_time.year,'%02i' % chunk_start_time.month))
                    if gitver!=None: #none is returned if not supported by python
                        retdict['gitversion']=gitver
                        retdict['caltables_gitversion']=gitver
            except ValueError:
                pass#not managed, no version
            if hasattr(rs_cal,'geo') and rs_cal.geo!=None and hasattr(rs_cal.geo,'data') and rs_cal.geo.data is not None:
                retdict['geo_corr']=rs_cal.geo.data[:,1]
            yield retdict
            chunk_start_time=chunk_end_time
  
        return

@dplkit.role.decorator.exposes_attrs_of_field('sounding_source')
@dplkit.role.decorator.exposes_attrs_of_field('tables')
@dplkit.role.decorator.exposes_attrs_in_chain(['hsrl_sounding','hsrl_Cxx'])
@dplkit.role.decorator.autoprovidenested(nestedclasses=[hau.Time_Z_Group,hau.rs_xfer,dict,mr.cal_vectors,hru.calibration_vector_table])
class dpl_calibration_narr(dplkit.role.narrator.aNarrator):
    """ HSRL Calibration Framestream Narrator. should only be created from the dpl_calibration object

        :param tables: caltables stream
        :param sounding_source: optional sounding source. will use configured from calvals if not provided at init

        exposed attributes:
  
        - hsrl_sounding (current sounding)
        - hsrl_Cxx (calibration structure)

        exposes field:

        - tables (dpl_calibration_tables_narr narrator object)
        - sounding_source (atmospheric profile source. should expose altitudeAxis)
    """
  
    def __init__(self,tables,sounding_source):
        self.tables=tables
        self.sounding_source=sounding_source
        if sounding_source is None:
            raise RuntimeError('Need a sounding source')
        self.Cxx=None
        self.sounding=None

    @property
    def process_control(self):
        return self.hsrl_process_control
    @property
    def corr_adjusts(self):
        return self.hsrl_corr_adjusts
    @property
    def instrument(self):
        return self.hsrl_instrument
    @property
    def hsrl_sounding(self):
        warnings.warn("retrieving HSRL_SOUNDING attribute is deprecated. It changes, so you should be using a data stream directly or indirectly. Detailed plan pending")
        return self.sounding
    @property
    def hsrl_Cxx(self):
        warnings.warn("retrieving HSRL_CXX attribute is deprecated. It changes, so you should be using a data stream directly or indirectly. Detailed plan pending")
        return self.Cxx

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
                    rs_Cxx = cu.update_Cxx(rs_constants,rs_cal,sounding,self.process_control,self.corr_adjusts)
                    if rs_Cxx is None:
                        continue
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
                warnings.warn("dpl_calibration_narr trying to use 0-length window. source is behind")
                print "behind times",chunk_end_time,chunk_start_time,nextcaltime
                chunk_end_time=nextcaltime
                if chunk_end_time<=chunk_start_time:
                    break
                #chunk_end_time=use_end_time
                #break
            #sounding=soundingsource(chunk_start_time)
            #setattr(self,'hsrl_constants',rs_constants)#FIXME exposing these in this way, when they are dynamic parts of the framestream is Bad Bad Bad
            setattr(self,'sounding',sounding)# should only be init parameters that can be exposed this way
            setattr(self,'Cxx',rs_Cxx)
            retdict= calframe.copy()
            retdict.update(dict(chunk_start_time=chunk_start_time,
                                chunk_end_time=chunk_end_time,
                                sounding=sounding,
                                rs_Cxx=rs_Cxx))
            print 'Yielding calframe ',chunk_start_time,chunk_end_time,sounding.times
            yield retdict
            chunk_start_time=chunk_end_time
  
        return


@dplkit.role.decorator.exposes_attrs_in_chain(['hsrl_calvals','hsrl_instrument','hsrl_process_control','hsrl_corr_adjusts'])
class dpl_calibration(dplkit.role.librarian.aLibrarian):
    """HSRL DPL Calibration Frame Stream

     :param instrument: instrument name
     :param process_control: process parameters file or structure. if not specified, uses default
     :param max_range_bin: maximum range bin, for trimming calibrations to reduce work to do

     exposes attributes:

     - hsrl_calvals
     - hsrl_instrument
     - hsrl_process_control
     - hsrl_corr_adjusts
    """
    def __init__(self,instrument,process_control=None,max_range_bin=5000, alternate_cal_dir=None, calibration_overrides=None,*args, **kwargs):
        if len(args):
            print 'Unused dpl_calibration args = ',args
        if len(kwargs):
            print "Unused dpl_calibration kwargs = ",kwargs.keys()
        super(dpl_calibration, self).__init__()
        self.instrument=instrument
        self.system_config=(process_control is None)
        if process_control is None:
            process_control='process_control.json'
        if isinstance(process_control,basestring):
            self.process_control_file=process_control
            self.process_control=None
        elif isinstance(process_control,jc.json_config):
            self.process_control_file=None
            self.process_control=process_control
        else:
            raise RuntimeError('process_control parameter must be None, a filename, or a json_config object')
        self.max_range_bin=max_range_bin
        self.alternate_cal_dir=alternate_cal_dir
        self.calibration_overrides=calibration_overrides
        import dpl_dynamic_atmospheric_profile_narrator as ddapn
        self.default_sounding_lib=ddapn.dpl_dynamic_hsrl_atmospheric_profile_librarian(self.instrument,edgepaddingIntervals=2)
        self.reload()

    def reload(self):
        self.reload_calvals()
        self.reload_process_control()
        #self.change_sounding_calvals()

    def change_sounding_calvals(self,**kwargs):
        self.default_sounding_lib.replace_constants(**kwargs)

    def reload_calvals(self):
        self.calvals=cu.calval_info(self.instrument)
    def reload_process_control(self):
        if self.process_control_file is not None:
            self.process_control=jc.json_config(locate_file(self.process_control_file,systemOnly=self.system_config),'process_defaults')
        if self.system_config and os.getenv('OVERRIDE_SYSTEM_DEFAULTS',None) is None:
            altpath=self.process_control.get_value('alternate_cal_dir','full_dir_path')
            if not (altpath is None or altpath=="None"):
                print "System Default process_defaults.json should not specify an alternative calibration path."
                hconfig=os.getenv("HSRL_CONFIG","a custom location, set HSRL_CONFIG to that path")
                print "Strongly recommended one of the two options:"
                print "- copy process_defaults.json to %s, customize the json there, and specify this non-default parameter at the commandline" % (hconfig)
                print "- or set HSRL_ALTERNATE_CAL_DIR to %s" % (altpath)
                print "Or, alternatively, do all of the following:"
                print "set OVERRIDE_SYSTEM_DEFAULTS to yes"
                if "custom location" in hconfig:
                    hconfig="~/hsrl_config"
                    print "make directory "+hconfig
                    print "set HSRL_CONFIG to "+hconfig
                print "copy process_defaults.json to "+hconfig
                print "modify it without renaming it, and it becomes the new system default for you"
                #import time
                #time.sleep(10)


    @property
    def hsrl_calvals(self):
        return self.calvals
    @property
    def hsrl_instrument(self):
        return self.instrument
    @property
    def hsrl_process_control(self):
        return self.process_control
    @property
    def hsrl_corr_adjusts(self):
        return self.hsrl_process_control.get_dict('corr_adjusts')

    def __repr__(self):
        return 'DPL HSRL Calibration Librarian (instrument="%s",max_range_bin="%i")' % (self.instrument,self.max_range_bin)

    def search(self,interval_start_time=None,interval_end_time=None,reverse_padding=None,min_alt_m=0,max_alt_m=50000,altres_m=0,\
        corr_adjusts=None,window_width_timedelta=None,mol_norm_alt=None,nocal=False, useconsts=None,alternate_cal_dir=None,calibration_overrides=None,
        requested_altitudes=None,*args, **kwargs):
        """
        DPL Generator function for starting a frame stream
        
        :param interval_start_time: start time (optional) if not specified, is calculated from window width
        :type interval_start_time: datetime
        :param interval_end_time: end time (optional) if not specified, is a continuous livestream
        :type interval_end_time: datetime
        :param reverse_padding: period to optionally subtract from now, in the case the window may bump against current data, or live stream
        :type reverse_padding: timedelta
        :param min_alt_m: minimum altitude in meters
        :param max_alt_m: maximum altitude in meters
        :param altres_m: altitude resolution in meters
        :param corr_adjusts: adjustments dictionary
        :param window_width_timedelta: (optional) used with start time or end time or now to describe the preferred window. 
        :type window_width_timedelta: timedelta
        """
        soundingsource=kwargs.pop('soundingsource',None)
        if len(args):
            print 'Unused dpl_calibration.read args = ',args
        if len(kwargs):
            print "Unused dpl_calibration.read kwargs = ",kwargs.keys()

        if useconsts is not None: #FIXME THIS NEEDS CLEANUP
            ret=useconsts
        else:
            if reverse_padding is None:
                reverse_padding=timedelta(seconds=0)
            timeinfo=time_frame.parse_timewindow(interval_start_time,interval_end_time,window_width_timedelta,datetime.utcnow()-reverse_padding)
            timeinfo['reverse_padding']=reverse_padding
            ret=dpl_constants_narr(instrument=self.instrument,timeinfo=timeinfo,calvals=self.calvals,
                    process_control=self.process_control,corr_adjusts=corr_adjusts,mol_norm_alt=mol_norm_alt)

        if nocal:
            return ret
        if soundingsource is not None:
            if requested_altitudes is not None:
                raise RuntimeError("Can't provide sounding source and altitudes.")
        else:
            if requested_altitudes is None:
                if altres_m==0:
                    altres_m=ret.hsrl_constants_first['binwidth']*3e8/2
                requested_altitudes=np.arange(min_alt_m,max_alt_m+altres_m*.1,altres_m)
            soundingsource=self.default_sounding_lib(ret.timeinfo,requested_altitudes)

        ret= dpl_calibration_tables_narr(ret,max_range_bin=self.max_range_bin, alternate_cal_dir=alternate_cal_dir or self.alternate_cal_dir)
        if calibration_overrides is not None or self.calibration_overrides is not None:
            if self.calibration_overrides:
                overrides=self.calibration_overrides.copy()
                if calibration_overrides:
                    overrides.update(calibration_overrides)
            else:
                overrides=calibration_overrides
            ret=dpl_override_calibration_tables_narr(ret,self.max_range_bin,alternate_cal_dir or self.alternate_cal_dir,**overrides)
        ret= dpl_calibration_narr(ret,soundingsource)
        return ret

def quickparse(val):
    if val is None:
        return datetime.utcnow()
    if isinstance(val,datetime):
        return val
    if isinstance(val,timedelta):
        return datetime.utcnow()-val
    assert(isinstance(val,basestring))
    try:
        return datetime.utcnow()-datetime(minutes=float(val))
    except ValueError:
        pass
    for f in ('%Y.%m.%d','%Y.%m.%dT%H:%M:%S'):
        try:
            return datetime.strptime(val,f)
        except ValueError:
            pass
    warnings.warn('couldn\'t parse "'+val+'"')
    return None

def quickrun(inst,start=None,end=None):
    calcount=0
    start=quickparse(start)
    end=quickparse(end)
    starttime=datetime.utcnow()
    cals=[]
    for f in dpl_calibration(inst)(interval_start_time=start,interval_end_time=end):
        print '.',
        cals.append(f)
    print len(cals),'calibrations for ',inst,start,end
    print 'took',datetime.utcnow()-starttime
    for i,cal in enumerate(cals):
        print i,cal['chunk_start_time'],cal['chunk_end_time']

def main():
    import sys
    if len(sys.argv)==1:
        print 'usage: %s instrument [start [end]]'
        ex=['gvhsrl','2012.02.10T16:00:00','2012.02.10T18:00:00']
        print 'Here is ',ex
        quickrun(*ex)
        
    else:
        quickrun(*sys.argv[1:])
        return
    #c=dpl_calibration('bagohsrl')
    #for i in c(interval_start_time=datetime(2013,2,1,0,0,0),interval_end_time=datetime(2013,2,2,0,0,0)):
    #    print i
    c=dpl_calibration('gvhsrl')
    dp=c(interval_start_time=datetime(2012,2,10,16,0,0),interval_end_time=datetime(2012,2,10,18,0,0))
    count=0
    for i in dp:
        #print i
        count+=1
    print count,'records'
    count=0
    for i in dp:
        #print i
        count+=1
    print count,'records second run'

if __name__ == '__main__':
    main()
