import dplkit.role.narrator
import dplkit.role.librarian
from datetime import datetime,timedelta
import sys,os 
import warnings

#FIXME this needs to be generalized, and not hsrl-bound
all_sounding_calvals=('sounding_type','virtual_sounding_update_interval','sounding_update_interval',
            'model_filename','static_cal_time','sounding_horizon','installation','sounding_id','latitude','longitude')

class dpl_dynamic_hsrl_atmospheric_profile_librarian(dplkit.role.librarian.aLibrarian):
    def __init__(self,instrument,edgepaddingIntervals=None,calvals=None,soundingdatapath=None):
        import lg_base.formats.calvals as cru
        import lg_dpl_toolbox.dpl.calvals_narrator as cvn
        self.cvn=cvn
        self.instrument=instrument
        self.soundingdatapath=soundingdatapath
        self.paddingends=edgepaddingIntervals or 0;
        self.calvals=calvals
        if self.calvals is None:
            self.calvals=cru.calvals_class(instrument=instrument,varlist=all_sounding_calvals)
        self.newconstants=None

    def replace_constants(self,sounding_type=None,**kwargs):
        if sounding_type is None:
            self.newconstants=None
        else:
            self.newconstants=dict(sounding_type=sounding_type)
            for k,v in kwargs.items():
                if k in all_sounding_calvals:
                    self.newconstants[k]=v
                else:
                    err="Parameter "+k+" isn't a sounding calval. Ignoring!"
                    warnings.warn(err)
                    raise RuntimeError(err)

    def search(self,timeinfo,requested_altitudes):
        calnar=self.cvn.dpl_calvals_narrator(timeinfo,self.newconstants or self.calvals)
        return dpl_dynamic_hsrl_atmospheric_profile_narrator(self.instrument,calnar,requested_altitudes,self.paddingends,self.soundingdatapath)

@dplkit.role.decorator.exposes_attrs_of_field('calnar')
@dplkit.role.decorator.exposes_attrs_of_field('first_sounder')#this provides is shared
class dpl_dynamic_hsrl_atmospheric_profile_narrator(dplkit.role.narrator.aNarrator):
    def __init__(self,instrument,calnar,requested_altitudes=None,edgepaddingIntervals=None,soundingdatapath=None):
        import lg_base.formats.calvals as cru
        import lg_dpl_toolbox.dpl.calvals_narrator as cvn
        self.instrument=instrument
        self.soundingdatapath=soundingdatapath
        self.paddingends=edgepaddingIntervals or 0;
        self.requested_altitudes=requested_altitudes
        self.calnar=calnar
        self.first_sounder=self.make_sounding(self.calnar.constants_first,self.calnar.timeinfo['starttime'],self.calnar.timeinfo['endtime'])

    def read(self):
        exc=None
        for calv in self.calnar:
            exc=None
            print 'getting soundings for ',self.instrument,calv
            try:
                for sounding in self.make_sounding(calv['rs_constants'],calv['chunk_start_time'],calv['chunk_end_time'],
                    isFirst=calv['chunk_start_time']<=self.timeinfo['starttime'],isLast=calv['chunk_end_time']>=self.timeinfo['endtime']):
                    yield sounding
            except Exception as e:
                print "Exception occurred in radiosonde sourcing. jumping to next"
                exc=e
        if exc is not None:
            raise

    def sounding_from(self,sounding_type,starttime=None,endtime=None,**kwargs):
        return self.make_sounding(fakeconstants,starttime or self.calnar.timeinfo['starttime'],endtime or self.calnar.timeinfo['endtime'],isFirst=True,isLast=True)

    def make_sounding(self,rs_constants,interval_start_time,interval_end_time,isFirst=False,isLast=False):
        import atmospheric_profiles.dpl.dpl_temperature_profiles as dtp#import dpl_soundingarchive,dpl_virtualradiosonde
        if rs_constants['sounding_type'] in ('virtual','remote_virtual','model','virtual_remote','virtual_cache','virtual cache'):
            from hsrl.dpl.HSRLLibrarian import HSRLLibrarian
            from lg_dpl_toolbox.dpl.NetCDFZookeeper import GenericTemplateRemapNetCDFZookeeper
            import lg_dpl_toolbox.filters.time_frame as time_frame# TimeGinsu
            m=None
            if 'gps_stream' not in rs_constants or rs_constants['gps_stream']=='hsrl':
                tempzoo=GenericTemplateRemapNetCDFZookeeper(self.instrument,user_read_mode='position')#,keepfields=fields)
                templib=HSRLLibrarian(instrument=self.instrument,zoo=tempzoo)#site=16)#None,datatype)
                m=templib(start=interval_start_time,end=interval_end_time)#,filetype='data')
                m=time_frame.TimeGinsu(m,'times',None)
            elif rs_constants['gps_stream']=='surfacemet':
                raise NotImplemented('GPS stream from automatic surfacemet CSV not implemented')
            else:
                raise RuntimeError('VRS needs a GPS stream. '+rs_constants['gps_stream']+' is not a known type')
            delt=5#90
            #print 'default'
            if 'installation' not in rs_constants or rs_constants['installation']=='ground':
                #print 'ground default'
                delt=60
            elif rs_constants['installation']=='shipborne':
                #print 'shipdefault'
                delt=15
            if 'virtual_sounding_update_interval' in rs_constants:
                delt=rs_constants['virtual_sounding_update_interval']
            elif 'sounding_update_interval' in rs_constants:
                delt=rs_constants['sounding_update_interval']
            delt=timedelta(minutes=delt)
            if isFirst:
                interval_start_time=interval_start_time-delt*self.paddingends
            if isLast and interval_end_time is not None:
                interval_end_time=interval_end_time+delt*self.paddingends
            gribcache=os.getenv('GRIB_CACHE','/arcueid/data/grib_cache')
            try:
                os.makedirs(gribcache)
            except:
                pass
            if not os.path.exists(gribcache):
                raise RuntimeError("Can't Access GRIB Cache directory %s. If this isn't an expected error, set environment variable \"GRIB_CACHE\" to something better and try again." % gribcache)
            parameters=dict(requested_altitudes=self.requested_altitudes \
                ,remote=False ,download=True)
            extrasearchparams=dict()
            if 'remote' in rs_constants['sounding_type']:
                parameters['remote']=True
            if rs_constants['sounding_type'] in ('model','virtual_cache','virtual cache'):
                parameters['download']=False
            if rs_constants['sounding_type'] == 'model':
                parameters['format']='model.grib2'
                parameters['predict_horizon']=24*365
                if 'model_filename' in rs_constants:
                    parameters['format']=rs_constants['model_filename']
                print 'VRS Forecast model will use grib2 file at '+os.path.join(gribcache,parameters['format'])
            if 'static_cal_time' in rs_constants and rs_constants['static_cal_time'] is not None and len(rs_constants['static_cal_time'])>2:
                import lg_base.core.read_utilities as hru
                extrasearchparams['fixed_time']=hru.convert_date_str(rs_constants['static_cal_time'])['datetime']
            if 'sounding_horizon' in rs_constants and rs_constants['sounding_horizon'] is not None:
                if 'predict_horizon' not in parameters or rs_constants['sounding_horizon']>parameters['predict_horizon']:
                    parameters['predict_horizon']=rs_constants['sounding_horizon']
            if 'predict_horizon' in parameters:
                    print "Predict horizon is ",parameters['predict_horizon'],'hours'
                    interval_start_time=interval_start_time-timedelta(hours=parameters['predict_horizon'])
            print "VRS SAMPLING FREQUENCY: ",delt
            soundinglib=dtp.dpl_virtualradiosonde(rs_constants['sounding_type'],gribcache,delt,**parameters)
            if m.provides!=None and 'latitude' in m.provides:
                print 'Got latitude stream'
                soundingnarr=soundinglib(m,starttime=interval_start_time,**extrasearchparams)
            else:
                fixed_position=(rs_constants['latitude'],rs_constants['longitude'])
                print 'no position stream/ using config',fixed_position
                soundingnarr=soundinglib(starttime=interval_start_time,endtime=interval_end_time,fixed_position=fixed_position,**extrasearchparams)
            ret=soundingnarr
        elif rs_constants['sounding_type'] in ('sparc','ssec'):
            import atmospheric_profiles.dpl.sparc_profiles as sparc#import dpl_soundingarchive,dpl_virtualradiosonde
            delt=timedelta(hours=3)
            if isFirst:
                interval_start_time=interval_start_time-delt*self.paddingends
            if isLast and interval_end_time is not None:
                interval_end_time=interval_end_time+delt*self.paddingends
            soundinglib=sparc.SPARCSondeLibrarian(self.soundingdatapath or self.instrument,self.requested_altitudes)
            ret=soundinglib(interval_start_time,interval_end_time)
        elif rs_constants['sounding_type'] in ('arm',):
            import atmospheric_profiles.dpl.arm_profiles as armp#import dpl_soundingarchive,dpl_virtualradiosonde
            delt=timedelta(hours=3)
            if isFirst:
                interval_start_time=interval_start_time-delt*self.paddingends
            if isLast and interval_end_time is not None:
                interval_end_time=interval_end_time+delt*self.paddingends
            soundinglib=armp.ARMSondeLibrarian(self.soundingdatapath or self.instrument,self.requested_altitudes)
            ret=soundinglib(interval_start_time,interval_end_time)
        elif rs_constants['sounding_type'] in ('NOAA raob',):
            import atmospheric_profiles.dpl.raob_profiles as raob#import dpl_soundingarchive,dpl_virtualradiosonde
            delt=timedelta(hours=12)
            if isFirst:
                interval_start_time=interval_start_time-delt*self.paddingends
            if isLast and interval_end_time is not None:
                interval_end_time=interval_end_time+delt*self.paddingends
            soundinglib=raob.dpl_raob(self.soundingdatapath or self.instrument,rs_constants['sounding_id'],self.requested_altitudes)
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
            if expire_duration:
                if isFirst:
                    interval_start_time=interval_start_time-expire_duration*self.paddingends
                if isLast and interval_end_time is not None:
                    interval_end_time=interval_end_time+expire_duration*self.paddingends
            soundinglib=dtp.dpl_soundingarchive(self.soundingdatapath or self.instrument,
                rs_constants['sounding_type'],
                rs_constants['sounding_id'],
                self.requested_altitudes,expire_duration=expire_duration)
            ret=soundinglib(interval_start_time-timedelta(days=3),interval_end_time)
        return ret


def main():
    import sys
    if len(sys.argv)<2:
        print 'usage: %s instrument [starttime [endtime]]' % (sys.argv[0])
        print 'times take the format YYYYMMDDTHHMMSS'
    fmt='%Y%m%dT%H%M%S'
    import lg_dpl_toolbox.filters.time_frame as tf
    print sys.argv[1:]
    sttime=datetime.strptime(sys.argv[2],fmt) if len(sys.argv)>2 else (datetime.utcnow()-timedelta(days=3))
    entime=datetime.strptime(sys.argv[3],fmt) if len(sys.argv)>3 else None
    timewindow=tf.parse_timewindow(sttime,entime)
    print sys.argv[1],timewindow
    for s in dpl_dynamic_hsrl_atmospheric_profile_narrator(sys.argv[1],timewindow):
        if False:
            print vars(s).keys(),s
        else:
            keys=vars(s).keys()
            keys.sort()
            for k in ('time','latit','longi','bot','top'):
                for mk in keys:
                    if mk.startswith('_'):
                        continue
                    if k in mk:
                        print mk+':',getattr(s,mk),' ',
            print ''


if __name__ == '__main__':
    sys.path.append(os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),'..','..','..'))
    main()