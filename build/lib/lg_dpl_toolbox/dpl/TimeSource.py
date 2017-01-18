import dplkit.role.narrator
import dplkit.role.filter
import lg_dpl_toolbox.filters.time_frame as time_frame
from datetime import datetime,timedelta
import copy
import lg_base.core.array_utils as hau
import numpy

class AddPseudoDeltaT(dplkit.role.filter.aFilter):
    def __init__(self,framestream,timefieldname,dtfieldname):
        super(AddPseudoDeltaT,self).__init__(framestream)
        self.s=framestream
        self.timefieldname=timefieldname
        self.dtfieldname=dtfieldname
        self.provides=copy.copy(framestream.provides)
        self.provides[dtfieldname]=dict(shortname=dtfieldname,type=hau.T_Array)

    def process(self):
        priorFrame=None
        priorFrameD=None
        averageDiff=0.0
        for f in self.s:
            assert(f is not None)
            f=copy.copy(f)
            _f=f
            if not isinstance(_f,dict):
                _f=vars(_f)
            tf=_f[self.timefieldname]
            dtf=hau.T_Array(numpy.zeros(tf.shape),summode='sum')
            _f[self.dtfieldname]=dtf
            for i in range(tf.size-1):
                dtf[i]=(tf[i+1]-tf[i]).total_seconds()
            if priorFrame is not None:
                if priorFrameD[self.dtfieldname].size>0 and tf.size>0:
                    priorFrameD[self.dtfieldname][-1]=(tf[0]-priorFrameD[self.timefieldname][-1]).total_seconds()
                yield priorFrame
            priorFrame=f
            priorFrameD=_f
            if priorFrameD[self.dtfieldname].size>1:
                priorFrameD[self.dtfieldname][-1]=priorFrameD[self.dtfieldname][-2]
        if priorFrame is not None:
            yield priorFrame

class AddAppendableTime(dplkit.role.filter.aFilter):
    def __init__(self,framestream,timefieldname,dtfieldname):
        super(AddAppendableTime,self).__init__(framestream)
        self.s=framestream
        self.timefieldname=timefieldname
        self.dtfieldname=dtfieldname
        self.provides=copy.copy(framestream.provides)
        self.provides[timefieldname]=dict(shortname=timefieldname,type=hau.T_Array)
        self.provides[dtfieldname]=dict(shortname=dtfieldname,type=hau.T_Array)

    def process(self):
        for f in self.s:
            f=copy.copy(f)
            _f=f
            if not isinstance(f,dict):
                _f=vars(f)
            #print _f.keys()
            _f[self.timefieldname]=hau.T_Array([_f['start']],summode='first')#+timedelta(seconds=_f['width'].total_seconds()/2.0)])
            _f[self.dtfieldname]=hau.T_Array([_f['width'].total_seconds() if _f['width']!=None and _f['width'].total_seconds()>0.0 else numpy.NaN],summode='sum')
            yield f


@dplkit.role.decorator.exposes_attrs_in_chain(['start_time','end_time','time_resolution'])
class TimeGenerator(dplkit.role.narrator.aNarrator):
    def __init__(self,start_time=None,end_time=None,width=None,time_resolution=None,time_step_count=None,timeinfo=None):
        self.timeinfo=timeinfo or time_frame.parse_timewindow(start_time,end_time,width,datetime.utcnow())
        if time_resolution is not None:
            self.resolution=time_resolution
        elif time_step_count is not None:
            self.resolution=timedelta(seconds=self.timeinfo['windowwidth'].total_seconds()/float(time_step_count))
        else:
            self.resolution=None
        self.provides=dict(
                start=dict(shortname='start',type=datetime),
                width=dict(shortname='width',type=timedelta)
            )

    @property
    def start_time(self):
        return self.timeinfo['starttime']

    @property
    def end_time(self):
        return self.timeinfo['endtime']
    @property
    def time_resolution(self):
        return self.resolution

    def read(self):
        st=self.timeinfo['starttime']
        en=self.timeinfo['endtime']
        while en is None or st<en:
            yield dict(start=st,width=self.resolution)
            st=st+self.resolution

@dplkit.role.decorator.exposes_attrs_of_field('gen')
class CompoundTimeGenerator(object):
    def __init__(self,generator):
        self.gen=generator
        self.iterat=iter(generator)
        self._isDone=False
        self.prior=None
        self.history=[]
        self.getTo(datetime(1980,1,1,0,0,0))
        self.first=self.prior

    @property
    def start_time(self):
        return self.first

    @property
    def isDone(self):
        return self._isDone

    def getTo(self,time):
        ret=[]
        while self.prior is None or self.prior['start']<=time and not self._isDone:
            if self.prior is not None:
                ret.append(self.prior)
                if (self.prior['start']+self.prior['width'])>time:
                    break
            try:
                self.prior=self.iterat.next()
                if not isinstance(self.prior,dict):
                    self.prior=vars(self.prior)
                #print 'keys in prior=',self.prior.keys()
                self.prior=dict(start=self.prior['start'],width=self.prior['width'])
                #print 'prior is now',self.prior
            except StopIteration:
                self._isDone=True
                print '***********timesource finished iterating'
                return ret
        return ret

    def getBinsFor(self,endtime,starttime=None,inclusive=False):
        ret=self.getTo(endtime)
        if starttime is not None:
            tmp=self.history+ret
            ret=[]
            for x in tmp:
                if x['start']+x['width']>starttime and not x in ret:
                    ret.append(x)
            self.history=copy.copy(ret)
        if len(ret)==0:
            assert(self.prior is not None)
            if not inclusive:
                if self._isDone:
                    return [self.prior['start']+self.prior['width']]
                return [self.prior['start']]
            return [self.prior['start']] + ([self.prior['start']+self.prior['width']] if (inclusive and endtime>self.prior['start']) or self._isDone else [])
        return [x['start'] for x in ret] + ([ret[-1]['start']+ret[-1]['width']] if (inclusive and endtime>ret[-1]['start']) or self._isDone else [])

def main():
    tg=TimeGenerator(start_time=datetime(2014,1,1,0,0,0),end_time=datetime(2014,1,2,0,0,0),time_resolution=timedelta(hours=1))
    ctg=CompoundTimeGenerator(tg)
    print ctg.getBinsFor(endtime=datetime(2014,1,1,11,59,0),inclusive=False)
    print ctg.getBinsFor(endtime=datetime(2014,1,1,11,59,0),inclusive=False)
    print ctg.getBinsFor(endtime=datetime(2014,1,1,12,0,0),inclusive=False)
    print ctg.getBinsFor(endtime=datetime(2014,1,1,15,59,0),inclusive=True)
    print ctg.getBinsFor(endtime=datetime(2014,1,1,15,59,0),inclusive=True)
    print ctg.getBinsFor(endtime=datetime(2014,1,1,16,0,0),inclusive=True)


if __name__ == '__main__':
    main()
