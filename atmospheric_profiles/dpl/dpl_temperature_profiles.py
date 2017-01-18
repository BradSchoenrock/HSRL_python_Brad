import dplkit.role.librarian
import dplkit.role.narrator
import dplkit.role.filter
from datetime import datetime,timedelta
import numpy
import dplkit.role.decorator
import lg_base.core.array_utils as hau
import os,calendar

from atmospheric_profiles.soundings.sounding_utilities import hydrostatic_interp


from scipy.interpolate import UnivariateSpline

class fake_rs_soundings(object):
    def __init__(self,narrator):
        import lg_dpl_toolbox.filters.time_frame as dt
        self.dt=dt
        self.narr=narrator
        self.source=dt.TimeTrickle(self.narr,timename='times',getClosest=True)

    def profile(self,time,d1,d2):
        return self.source(time)

    def append(self,d1):
        raise RuntimeError("Append sounding shouldn't have to be called since it's using a DPL construct")


@dplkit.role.decorator.exposes_attrs_in_chain(['altitudeAxis'])
@dplkit.role.decorator.autoprovide(frameclass=hau.Time_Z_Group,reuseGenerator=False)
class dpl_singlesounding_narr(dplkit.role.narrator.aNarrator):

    def __init__(self,frame):
        self.frame=frame
        assert('altitudes' in self.provides)

    @property
    def altitudeAxis(self):
        return self.frame.altitudes

    def read(self):
        yield self.frame

@dplkit.role.decorator.exposes_attrs_in_chain(['altitudeAxis'])
@dplkit.role.decorator.autoprovide(frameclass=hau.Time_Z_Group,reuseGenerator=False)
class dpl_soundingarchive_narr(dplkit.role.narrator.aNarrator):
    def __init__(self,host,interval_start_time,interval_end_time,expire_duration):
        self.host=host
        self.interval_start_time=interval_start_time
        self.interval_end_time=interval_end_time
        self.expire_duration=expire_duration
        self.requested_altitudes=self.host.requested_altitudes
        import atmospheric_profiles.soundings.sounding_utilities as su
        self.su=su

    @property
    def altitudeAxis(self):
        return self.requested_altitudes

    def badSounding(self,sounding,comparesounding=None):
        if sum(sounding.pressures>1.0)<2:
            return True
        if comparesounding is not None:
            if sounding.times==comparesounding.times:
                return True
        #print vars(sounding)
        return False

    def read(self):
        rs_soundings=None
        use_time=self.interval_start_time
        req_time=self.interval_start_time#-(self.expire_duration or timedelta(days=1))
        offset=0
        last_sounding=None
        nonetries=5
        while self.interval_end_time==None or use_time<self.interval_end_time:
            #if (self.interval_end_time is None and req_time>datetime.utcnow()) or (self.interval_end_time is not None and req_time>self.interval_end_time):#+timedelta(days=5)):
            #    break
            if rs_soundings is not None:
                sounding=rs_soundings.profile(req_time,[],[],offset=offset)
                if sounding==None:
                    rs_soundings=None
                else:
                    offset=1#after first successful sounding retrieval, keep getting one after the time of that one
                    sounding.sample_time=sounding.times
            if rs_soundings is None or (last_sounding is not None and sounding.times==last_sounding.times):
                new_soundings = self.su.sounding_archive(#look for more
                                                         self.host.instrument,
                                                         self.host.sounding_type,
                                                         self.host.sounding_id,
                                                         req_time,
                                                         self.requested_altitudes,offset=offset)
                #print 'read sounding'
                #    new_soundings=None
                if last_sounding is not None and (new_soundings is None or new_soundings.profile(req_time,[],[],offset=offset) is None):#old sonde
                    if self.interval_end_time is not None:
                        break
                    print 'yielding old sounding ',last_sounding.times,' because end time is ',self.interval_end_time
                    yield last_sounding
                    continue
                #print 'setting'
                if rs_soundings is None:
                    if new_soundings is None:
                        if nonetries==0:
                            print 'No soundings found. breaking the bad loop'
                            break
                        nonetries=nonetries-1
                    rs_soundings=new_soundings
                elif new_soundings is not None:
                    rs_soundings.append(new_soundings)
                continue
            if self.badSounding(sounding,last_sounding):
                if req_time<sounding.times:
                    req_time=sounding.times
                else:
                    req_time+=timedelta(minutes=5)
                print 'bad sounding',sounding.times,'now looking for',req_time
                continue
            if last_sounding is not None and sounding.times>self.interval_start_time:
                print 'yielding sounding for ',last_sounding.times
                yield last_sounding
            print 'NEXT sounding'
            last_sounding=sounding
            #sounding.times=use_time
            if use_time<sounding.times:
                use_time=sounding.times
            req_time=use_time
            if self.expire_duration is not None:
                req_time+=self.expire_duration
        if last_sounding is None:
            raise RuntimeError('ran out of radiosondes')
        yield last_sounding

class dpl_soundingarchive(dplkit.role.librarian.aLibrarian):
    def __init__(self,instrument,sounding_type,sounding_id,requested_altitudes,expire_duration=None):
        self.instrument=instrument
        self.sounding_type=sounding_type
        self.sounding_id=sounding_id
        self.requested_altitudes=requested_altitudes
        self.expire_duration=expire_duration

    def search(self,interval_start_time,interval_end_time=None,expire_duration=None):#none may just mean end is determined by caller
        return dpl_soundingarchive_narr(self,interval_start_time,interval_end_time,expire_duration or self.expire_duration)

@dplkit.role.decorator.exposes_attrs_in_chain(['altitudeAxis'])
class dpl_radiosonderesample(dplkit.role.filter.aFilter):
    def __init__(self,framestream,altitudename,requested_altitudes,fielddescriptions,method='linear'):
        super(self.__class__,self).__init__(framestream)
        self.framestream=framestream
        self.requested_altitudes=requested_altitudes
        self.altitudename=altitudename
        self.fielddescriptions=fielddescriptions
        self.method=method
        assert(method in ('linear',))

    @property
    def altitudeAxis(self):
        return self.requested_altitudes

    def process(self):
        destlen=self.requested_altitudes.shape[0]
        for frame in self.framestream:
            if frame!=None:
                alts=getattr(frame,self.altitudename).copy()
                altmask=numpy.isfinite(alts)
                for varn in vars(frame).keys():#{'tdry':[-150,200],'pres':[1,1000],'rh':[0,100]}.items():
                    if varn.startswith('_'):
                        continue
                    if varn==self.altitudename:
                        continue
                    v=getattr(frame,varn)
                    if not isinstance(v,hau.Z_Array):
                        continue
                    goodrange=None
                    if varn in self.fielddescriptions:
                        goodrange=self.fielddescriptions[varn]
                    cl=type(v)
                    realmask=numpy.logical_and(altmask,numpy.isfinite(v))
                    if goodrange!=None:
                        if goodrange[0]!=None:
                            realmask=numpy.logical_and(realmask,v>=goodrange[0])
                        if goodrange[1]!=None:
                            realmask=numpy.logical_and(realmask,v<=goodrange[1])
                    order=numpy.argsort(alts[realmask])
                    realmask=numpy.arange(alts.size)[realmask][order]
                    alt=alts[realmask]
                    v=v[realmask]
                    sourcelen=alt.shape[0]
                    shp=list(v.shape)
                    shpi=-1
                    for i in range(len(shp)):
                        if shp[i]==sourcelen and shp[i]>0:
                            shpi=i
                            shp[i]=destlen
                    if shpi<0:
                        cantcontinue=True
                        break
                    
                    if varn in ['pressures']: #log interp provides best fit to pressure
                       #newval=numpy.exp(numpy.interp(self.requested_altitudes.ravel() \
                       #                              ,alt.ravel(),numpy.log(v.ravel())))
                       #hydrostatic fit
                       newval = hydrostatic_interp(v.ravel()
                                             ,alt.ravel(),self.requested_altitudes.ravel())
                    else:
                        newval=numpy.interp(self.requested_altitudes.ravel(),alt.ravel(),v.ravel())
                    newval=cl(newval.reshape(shp))
                    #print 'profile resample',varn,'from',v,'to',newval
                    setattr(frame,varn,newval)
                    if 0: #varn == 'pressures':
                         import matplotlib.pylab as plt
                         plt.figure(777777)
                         plt.plot(newval,self.requested_altitudes/1000.0)
                         ax = plt.gca()
                         ax.set_xscale('log')
                         plt.xlabel('pressure')
                         plt.ylabel('altitude(km)')
                setattr(frame,self.altitudename,self.requested_altitudes.copy())
            yield frame


@dplkit.role.decorator.autoprovide(frameclass=hau.Time_Z_Group,reuseGenerator=False)
class dpl_virtualradiosondearchive_narr(dplkit.role.narrator.aNarrator):
    def __init__(self,basedir,starttime,endtime,filename='NWPSondes.nc',zoo=None):
        self.basedir=basedir
        self.start=starttime
        self.end=endtime
        tmpstart=starttime-timedelta(days=1)
        self.starttime=datetime(tmpstart.year,tmpstart.month,1,0,0,0)
        self.zoo=zoo
        self.filename=filename

    def read(self):
        t=self.starttime
        while t<datetime.utcnow() and (self.end==None or self.end>=t):
            folder = os.path.join(self.basedir,'%04i' % t.year,'%02i',t.month)
            monthrange=calendar.monthrange(t.year,t.month)
            monthdur=timedelta(days=monthrange[1])
            r={'path':os.path.join(folder,self.filename),'filename':self.filename,'start':t,'width':monthdur}
            if self.zoo:
                yield self.zoo.open(self.zoo(r),firsttime=self.start,lasttime=self.end)
            else:
                yield r
            t=t+monthdur


class dpl_virtualradiosondearchive(dplkit.role.librarian.aLibrarian):
    def __init__(self,instrument,requested_altitudes,do_interpolate=True,zoo=None):
        self.instrument=instrument
        self.requested_altitudes=requested_altitudes
        self.do_interpolate=do_interpolate and (requested_altitudes is not None)
        self.zoo=zoo
        if zoo==None:
            raise RuntimeError('Need to make a defualt vrsa zookeeper')

    def search(self,starttime,endtime=None):
        #print 'Loading hru for HSRLLibrarian'
        import lg_dpl_toolbox.core.archival as hru
        basedir=hru.get_path_to_data(self.instrument,None)
        ret=dpl_virtualradiosondearchive_narr(basedir,starttime,endtime,self.zoo)
        import lg_dpl_toolbox.filters.time_frame as dt
        ret=dt.TimeGinsu(ret,'times')
        if self.do_interpolate:
            ret=dpl_radiosonderesample(ret,'altitudes',self.requested_altitudes,{'temps':[50,500],'pressures':[1,1000],'dew_points':[150,500],'frost_points':[50,500],'altitudes':None})
        return ret

@dplkit.role.decorator.autoprovide(frameclass=hau.Time_Z_Group,reuseGenerator=False)
class dpl_virtualradiosonde_narr(dplkit.role.narrator.aNarrator):
    def __init__(self,name,vrs,stream,expire_duration=None,fixed_time=None,fixed_position=None,starttime=None,endtime=None,useNarrator=False):
        assert(stream or expire_duration)
        self.stream=stream
        self.name=name
        self.vrs=vrs
        self.narrate=useNarrator
        self.expire_duration=expire_duration
        self.static_time=fixed_time
        self.fixed_position=fixed_position
        self.starttime=starttime
        self.endtime=endtime
        import atmospheric_profiles.soundings.sounding_utilities as su
        self.su=su
        if stream!=None:
            if 'latitude' not in stream.provides or 'longitude' not in stream.provides:
                raise RuntimeError('No coordinates in position stream')
        else:
            if fixed_position is None or starttime==None:
                raise RuntimeError('No coordinates specified')

    def __timeposstream(self):
        if self.stream:
            sentStart=False
            nextSend=None
            for f in self.stream:
                if not (hasattr(f,'longitude') and hasattr(f,'latitude')):
                    raise RuntimeError('source stream lied. coordinates missing from position stream')
                if f.longitude[0]>360 or f.longitude[0]<-180 or f.latitude[0]>90 or f.latitude[0]<-90:
                    continue
                if nextSend is not None and nextSend>f.start:
                    continue
                if not sentStart:
                    if self.starttime is not None:
                        yield dict(datetime=self.starttime,start=self.starttime,latitude=f.latitude[0],longitude=f.longitude[0])
                    sentStart=True
                yield dict(datetime=f.start,start=f.start,latitude=f.latitude[0],longitude=f.longitude[0])
                if self.expire_duration is not None:
                    nextSend=f.start+self.expire_duration
                if self.endtime is not None and f.start>=self.endtime:
                    return
        else:
            from collections import namedtuple
            t=self.starttime
            while self.endtime==None or t<=self.endtime:
                yield dict(datetime=t,start=t,latitude=self.fixed_position[0],longitude=self.fixed_position[1])
                t+=self.expire_duration

    def convert(self,profile,posframe):
        sounding=hau.Time_Z_Group()
        if 'cache' in self.name:
            sounding.sounding_type='virtual'
            sounding.sounding_id='Cached Forecast'
            sounding.station_id='Cached Forecast'
        elif 'virtual' in self.name:
            sounding.sounding_type='virtual'
            sounding.sounding_id='NWP Virt'
            sounding.station_id='NWP Virt'
        else:
            sounding.sounding_type='model'
            sounding.sounding_id='Static GRIB'
            sounding.station_id='Static GRIB'
        sounding.latitude=profile['lat'][0]
        sounding.longitude=profile['lon'][0]
        sounding.times=datetime(1970,1,1,0,0,0)+timedelta(seconds=profile['base_time'])+timedelta(seconds=profile['time_offset'][0])
        if posframe is None:                
            sounding.sample_latitude=sounding.latitude
            sounding.sample_longitude=sounding.longitude
            sounding.sample_time=sounding.times
        else:
            sounding.sample_latitude=posframe['latitude']
            sounding.sample_longitude=posframe['longitude']
            sounding.sample_time=posframe['start']
        minlen=profile['alt'].size
        for k,v in profile.items():
            if hasattr(v,'shape'):
                if profile['alt'].size != v.size:
                    print 'WARNING SOUNDING ATTRIBUTE '+k+' has a length',v.size,'while alts are',profile['alt'].size
                    profile[k]=v[:profile['alt'].size]
                    if minlen>v.size:
                        minlen=v.size
        if minlen<profile['alt'].size:
            print "TRIMMING TO SHORTENED SOUNDING ATTRIBUTE of length",minlen
            for k,v in profile.items():
                if hasattr(v,'shape'):
                    profile[k]=v[:minlen]

        sounding.top=numpy.max(profile['alt'])
        sounding.bot=numpy.min(profile['alt'])
        sounding.altitudes=hau.Z_Array(profile['alt'])
        if 'tkel' in profile:
            sounding.temps=hau.Z_Array(profile['tkel'])
        else:
            sounding.temps=hau.Z_Array(profile['tdry'])+273.15
        sounding.pressures=hau.Z_Array(profile['pres'])
        sounding.dew_points=hau.Z_Array(self.su.cal_dew_point(hau.Z_Array(profile['rh']),sounding.temps))
        sounding.frost_points=hau.Z_Array(self.su.cal_frost_point(sounding.dew_points))
        return sounding

    def isNewSounding(self,s1,s2):
        #if s1.times!=s2.times:
        #    return True
        #if s1.latitude!=s2.latitude:
        #    return True
        #if s1.longitude!=s2.longitude:
        #    return True
        #if s1.top!=s2.top and s1.bot!=s2.bot:
        #    return True
        if s1.temps.size!=s2.temps.size:
            return True
        if (s1.temps!=s2.temps).any():
            return True
        return False            

    def isNewProfile(self,s1,s2):
        #if s1.times!=s2.times:
        #    return True
        #if s1.latitude!=s2.latitude:
        #    return True
        #if s1.longitude!=s2.longitude:
        #    return True
        #if s1.top!=s2.top and s1.bot!=s2.bot:
        #    return True
        if s1['tdry'].size!=s2['tdry'].size:
            return True
        if (s1['tdry']!=s2['tdry']).any():
            return True
        return False            

    def read(self):
        sounding=None
        if self.narrate:
            gen=self.vrs(self.__timeposstream())
        else:
            gen=self.__timeposstream()
        priorsounding=None
        if not self.narrate:
            for f in gen:
                #print 'gen step'
                if f==None or not ('latitude' in f and 'longitude' in f):
                    print 'Sounding doesnt exist or has no position'
                    continue
                #print 'Getting',f['latitude'],f['longitude']
                try:
                    profile=self.vrs(f['start'] if self.static_time==None else self.static_time,lat=f['latitude'],lon=f['longitude']).copy()
                except ValueError as e:
                    print 'Value Error occurred in VRS:',e            
                    continue

                sounding=profile#self.convert(profile,f)
                if priorsounding is None or self.isNewProfile(priorsounding,sounding):
                    yield self.convert(sounding,f)
                    priorsounding=sounding
        else:
            for f in gen:
                #print 'gen step'
                #print 'Narrating ',f['latitude'],f['longitude']
                if f is None:
                    print 'Sounding doesnt exist'
                    continue
                profile=f
                if priorsounding is None or self.isNewProfile(priorsounding,sounding):
                    priorsounding=sounding
                sounding=profile#self.convert(profile,None)
                if priorsounding and self.isNewProfile(priorsounding,sounding):
                    yield self.convert(priorsounding,None)
                #else:
                #    print 'Sounding isnt new'
            if sounding is not None:
                print "PAST THE LAST PROFILE AVAILABLE. Last profile you'll ever get"
                yield self.convert(sounding,None)
        if sounding is None:
            raise RuntimeError('No VRS Sondes for requested time range!')


class fakevrs(object):
    def __init__(self,alwaysreturn):
        self.r=alwaysreturn

    def __call__(self,*args,**kwargs):
        return self.r.copy()

def pressfromtemp(**kwargs):
    return kwargs.pop('pressure_tdry')

def temper(**kwargs):
    return kwargs.pop('tdry')

class dpl_virtualradiosonde(dplkit.role.librarian.aLibrarian):
    def __init__(self,name,cache_dir,time_threshhold,requested_altitudes=None,download=True,do_interpolate=True,narrate=True,remote=False,**kwargs):
        initargs=kwargs
        self.name=name
        self.requested_altitudes=requested_altitudes
        self.expire_duration=time_threshhold
        self.do_interpolate=do_interpolate and (requested_altitudes is not None)
        assert('levels' not in initargs)
        assert('level_type' not in initargs)
        from lg_base.core.decoratortools import package_min_version
        if package_min_version('virtual_radiosonde_source','0.1.3'):
            initargs['channels']={'pres':pressfromtemp,'pressure':pressfromtemp,'tkel':temper,'tdry':'Temperature','rh':'Relative humidity','alt':'Geopotential Height'}
            initargs['level_type']='pressure'
            initargs['levels']=None
        else:
            pass
        if time_threshhold!=None:
            assert('time_threshold' not in initargs)
            initargs['time_threshold']=time_threshhold
        if 'download_method' not in initargs:
            initargs['download_method']='rsync'
        if remote:
            assert('format' not in initargs)
            initargs['format']='rsync://guesstimator.ssec.wisc.edu/gfs/%(yyyy)04d/%(jjj)03d/gfs.p5.%(yyyy)04d%(mm)02d%(dd)02d_%(hh)02d_%(ooo)03d.ldm.grib2'
        if cache_dir!=None:
            if False and narrate:
                self.doNarrate=True
                from virtual_radiosonde_source.vrsNarrator import VirtualRadiosondeNarrator
                self.vrs=VirtualRadiosondeNarrator(cache=cache_dir,download=download, **initargs)
            else:
                self.doNarrate=False
                from virtual_radiosonde_source.vrs import VirtualRadiosonde
                self.vrs=VirtualRadiosonde(cache=cache_dir,download=download, **initargs)
            if package_min_version('virtual_radiosonde_source','0.1.3'):
                self.vrs.levels=None
        else:
            self.vrs=fakevrs({'tdry':numpy.array([5,numpy.NaN,5,5]),'pres':numpy.array([4,4,4,4]),'alt':numpy.array([-500,2000,7000,33000])})

    def search(self,latlontime_stream=None,expire_duration=None,fixed_time=None,fixed_position=None,starttime=None,endtime=None):
        ret=dpl_virtualradiosonde_narr(self.name,self.vrs,latlontime_stream,
                                       expire_duration=expire_duration if (self.expire_duration!=None and expire_duration!=None) else self.expire_duration,fixed_time=fixed_time,
                                       fixed_position=fixed_position,starttime=starttime,endtime=endtime,useNarrator=self.doNarrate)
        if self.do_interpolate:
            ret=dpl_radiosonderesample(ret,'altitudes',self.requested_altitudes,
                {'temps':[50,500],'pressures':[1,1000],'dew_points':[150,500],'frost_points':[50,500],'altitudes':None},
                method='linear')
        return ret

def main():
    import sys
    from lg_dpl_toolbox.dpl.NetCDFZookeeper import GenericTemplateRemapNetCDFZookeeper
    from hsrl.dpl.HSRLLibrarian import HSRLLibrarian
    from lg_dpl_toolbox.filters.time_frame import TimeGinsu
    datatype='ahsrl' if len(sys.argv)<2 else sys.argv[1]
    et=datetime.utcnow() if len(sys.argv)<4 else datetime.strptime(sys.argv[3],'%Y%m%dT%H%M%S')
    st=(et-timedelta(days=.5)) if len(sys.argv)<3 else datetime.strptime(sys.argv[2],'%Y%m%dT%H%M%S')
    #fields=['times','telescope_position','telescope_rotation','telescope_rotation_measured','telescope_elevation','telescope_accelerometer_raw']#,'superseedlasercontrollog','laserpowervalues']
    zoo=GenericTemplateRemapNetCDFZookeeper(datatype,user_read_mode='position')#,keepfields=fields)
    lib=HSRLLibrarian(instrument=datatype,zoo=zoo)#site=16)#None,datatype)
    m=lib(start=st,end=et)#,filetype='data')
    m=TimeGinsu(m,'times')
    vr_lib=dpl_virtualradiosonde('/home/jpgarcia/sonde/VirtualRadiosondeFromNWP/tmp_cache',timedelta(minutes=5))
    m=vr_lib(m,numpy.arange(0,30000+.1,15),timedelta(minutes=5))
    for f in m:
        print vars(f)
        print f.temps.shape,f.temps

if __name__ == '__main__':
    main()
