from datetime import datetime,timedelta
import sys,os

def doFetch(cachepath,start_time=None,end_time=None,remote=True):
    import atmospheric_profiles.dpl.dpl_temperature_profiles as dtp
    moreargs=dict(remote=remote,download=True)
    if os.getenv('FORCE_MODEL',None) is not None:
        moreargs['format']=os.getenv('FORCE_MODEL',None)
        moreargs['remote']=False
        moreargs['download']=False
        moreargs['predict_horizon']=24*365
    soundinglib=dtp.dpl_virtualradiosonde("virtual_queueing",cachepath,timedelta(minutes=10),\
        do_interpolate=False,**moreargs)
    for x in soundinglib(starttime=start_time,endtime=end_time,fixed_position=(45.0,145.0)):
        #print x.temps
        pass

formats=('%Y-%m-%dT%H:%M:%S','%Y-%m-%d')

def tryparsetime(st,now=None):
    try:
        dur=float(st)
        return (now or datetime.utcnow())+timedelta(hours=dur)
    except ValueError:
        pass
    for f in formats:
        try:
            return datetime.strptime(st,f)
        except ValueError:
            pass
    print 'add format to',formats
    raise RuntimeError('Dont know how to parse time '+st)

def main():
    if len(sys.argv)==1:
        print 'Usage:\n'
        print '\tfill_cache.py cache_path hours  (will fill from now to so many hours in the future)'
        print '\tfill_cache.py cache_path starttime  (will fill from start time to now)'
        print '\tfill_cache.py cache_path endtime  (will fill from now to end time)'
        print '\tfill_cache.py cache_path starttime endtime  (will fill from start time to end time)'
        print 'available time formats:'
        for f in formats:
            print '\t%s' % (f)
        return
    cachepath=sys.argv[1]
    now=datetime.utcnow()
    if len(sys.argv)==3:
        t=tryparsetime(sys.argv[2],now)
        if t>now:
            end_time=t
            start_time=now
        else:
            start_time=t
            end_time=now
    else:
        start_time=tryparsetime(sys.argv[2])
        end_time=tryparsetime(sys.argv[3],start_time)
    doFetch(cachepath,start_time=start_time,end_time=end_time)

if __name__ == '__main__':
    print os.path.join(os.path.dirname(sys.argv[0]),'..','..')
    sys.path.append(os.path.join(os.path.dirname(sys.argv[0]),'..','..'))
    main()