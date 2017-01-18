#!/usr/bin/env python
from datetime import datetime,timedelta
#import hsrl.utils.json_config as jc
import sys,os
from time import sleep
import matplotlib.pyplot as plt
import threading
import multiprocessing
import dplkit.role.filter
from functools import partial

def dodelay(timeout=0):
    plt.ginput(timeout=timeout)

class waitforready(object):
    def __init__(self,chunksize=5.0,single_process=False):
        if not single_process:
            self.cond=multiprocessing.Event()
            self.lock=multiprocessing.Lock()
        else:
            self.cond=threading.Event()
            self.lock=threading.Lock()
        self.chunksize=chunksize
        #self.set()

    def set(self):
        with self.lock:
            print "@@@@@@@@ set"
            self.cond.set()
        sleep(self.chunksize)

    @property
    def isready(self):
        return self.cond.is_set()
    
    def clear(self):
        with self.lock:
            print "@@@@@@@@ clear"
            ret=self.cond.is_set()
            self.cond.clear()

    def waitforready(self,max_timeout=None,chunksize=None):
        if chunksize is None:
            chunksize=self.chunksize
        while max_timeout is None or max_timeout>0.0:
            if self.isready:
                print "@@@@@@@@ is ready"
                break
            plt.ginput(timeout=max_timeout if max_timeout is not None and max_timeout < chunksize else chunksize)
            print '@@@@@@@@@@@@@@@@@ ginput waited'
            if max_timeout is not None:
                max_timeout-=chunksize
        return self.clear()


def raw_input2(prompt=None,timeout=None):
    import sys
    from select import select
    if prompt!=None:
        sys.stdout.write(prompt)
        sys.stdout.flush()
    if True:
        dodelay(timeout)
        timeout=None
    rlist, _, _ = select([sys.stdin], [], [], timeout)
    if rlist:
        s = sys.stdin.readline()
        return s
    print
    return None

class signalready_filter(dplkit.role.filter.aFilter):
    def __init__(self,stream,readier):
        super(signalready_filter,self).__init__(stream)
        self.readier=readier
        self.stream=stream

    def process(self):
        for x in self.stream:
            if self.readier is not None:
                self.readier.set()
            yield x

class input_filter(dplkit.role.filter.aFilter):
    def __init__(self,stream,prompt,sleeptime,interactive=False,readier=None):
        super(input_filter,self).__init__(stream)
        self.stream=stream
        self.prompt=prompt
        self.sleeptime=sleeptime
        self.interactive=interactive
        self.readier=readier

    def process(self):
        import matplotlib.pyplot as plt
        if self.interactive:
            plt.ion()
        else:
            plt.ioff()
        for f in self.stream:
            #plt.show(block=False)
            yield f
            dodelay(0.1)
            if self.interactive:
                if raw_input2(self.prompt,self.sleeptime) is not None:
                    return
            elif self.readier is not None:
                self.readier.waitforready()
            else:
                dodelay(self.sleeptime)


def livestream(instrument,hours,starttime=None,process_control=None,display_defaults='min_plots.json',maxalt_km=15,minalt_km=0,now_delay=None):
    """
    Quick and simple livestream of HSRL data using DPL constructs.

    :param instrument: instrument name ('gvhsrl','bagohsrl', etc)
    :param hours: numeric number of hours to display
    :param starttime: optional starttime, if stream should begin in the past. it will still stream thru all data, and approach now regardless of this
    :type starttime: datetime
    :param process_control: process control file or object for processing parameters
    :param display_defaults: graphics output configuration file
    :param maxalt_km: maximum altitude in km. default is 15
    :param minalt_km: minimum altitude in km. default is 0
    :param now_delay: offset to consider the now time. "now" minus this value is the rolling end of the window.
    :type now_delay: timedelta
    """

    import lg_dpl_toolbox.filters.substruct as frame_substruct
    frame_substruct.SubstructBrancher.multiprocessable=False 
    from hsrl.dpl.dpl_hsrl import dpl_hsrl
    import hsrl.dpl.dpl_artists as artists
    import lg_dpl_toolbox.dpl.dpl_artists as tools_artists
    import lg_dpl_toolbox.filters.threading_filter as tf
    window=timedelta(seconds=60*60*hours)
    print 'display_defaults=',display_defaults
    #time.sleep(5)
    figcontainer=None
    sigg=waitforready(single_process=True)
    dplobj=dpl_hsrl(instrument=instrument,process_control=process_control)
    dplgen=dplobj(start_time_datetime=starttime,reverse_padding=timedelta(seconds=60) if now_delay==None else now_delay,
        window_width_timedelta=window,min_alt_m=minalt_km*1000.0,max_alt_m=maxalt_km*1000.0,with_profiles=False)
    #dplgen=signalready_filter(dplgen,sigg)
    dplgen=tf.forking_filter(dplgen,2)
    dplgen=tf.aThreading_filter(dplgen,2)
    #dplgen=signalready_filter(dplgen,sigg)
    maxtimeout=30.0
    dplgen=tf.wait_filter(dplgen,waitfunc=partial(plt.ginput,timeout=maxtimeout/6),maxduration=timedelta(seconds=maxtimeout))
    artist=artists.dpl_images_artist(framestream=tools_artists.dpl_window_caching_filter(dplgen,window,includeIncompleteFrames=True),
        instrument=instrument,max_alt=maxalt_km*1000.0,processing_defaults=dplgen.hsrl_process_control,
        display_defaults=display_defaults,figurecontainer=figcontainer)
    #artist=tf.aThreading_filter(artist,1)
    #sleeptime=10
    #artist=input_filter(artist,"Pausing for %i seconds. Press enter to quit: " %(sleeptime),sleeptime,readier=sigg)
    return artist

def loadandrun(*args,**kwargs):
    fr=livestream(*args,**kwargs)
    for x in fr:
        fr.figs.showall()
        pass

def onaProcess(target,*args,**kwargs):
    import multiprocessing
    t=multiprocessing.Process(target=target,args=args,kwargs=kwargs)
    t.start()
    return t

def onaThread(target,*args,**kwargs):
    import threading
    t=threading.Thread(target=target,args=args,kwargs=kwargs)
    t.start()
    return t


def main():
    addpath=os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]),os.path.pardir,os.path.pardir))
    sys.path.insert(0,addpath)
    idx=2
    parms=dict(hours=2,now_delay=timedelta(seconds=30),display_defaults='min_plots.json')
    sleeptime=10
    while idx<len(sys.argv) and sys.argv[idx][0]=='-':
        sw=sys.argv[idx]
        parm =sys.argv[idx+1]
        idx+=2
        if sw in ('-u','--update'):
            sleeptime=float(parm)
        elif sw in ('-a','--altmax','--maxalt'):
            parms['maxalt_km']=float(parm)
        elif sw in ('-d','--display','--displaydefaults'):
            parms['display_defaults']=parm
        elif sw in ('-h','--duration','--hours'):
            parms['hours']=float(parm)
        elif sw in ('-i','--delay','--ignore'):
            parms['now_delay']=timedelta(seconds=float(parm))
        elif sw in ('-s','--start'):
            try:
                parms['starttime']=datetime.utcnow()-timedelta(hours=float(parm))
            except ValueError:
                parms['starttime']=datetime.strptime(parm,'%Y%m%dT%H%M')
        elif sw ('-p','--processcontrol','--parameters','--process'):
            parms['process_control']=parm
        else:
            RuntimeError('Unknown parameter %s : %s' % (sw,parm))

    if len(sys.argv)!=idx:
        print 'unused parameters: '+repr(sys.argv[idx:])
        print 'usage: %s instrument [-a maxalt_km] [-d display_json] [-h hours_duration]'
        print '\t -s --start :\tstart time as hours before now, or YYYYMMDDThhmm'
        print '\t -a --altmax :\tspecify maximum altitude in km (default 15)'
        print '\t -d --display :\tspecify a display json to use (default min_plots.json)'
        print '\t -h --hours :\thours to display (default 2)'
        print '\t -i --ignore :\tignore x seconds prior to now (default 30)'
        print '\t -u --update :\tseconds between updates (default 10)'
        print '\t -p --process :\tprocess control json filename'
        return 0
    instrument=sys.argv[1]
    artist=livestream(instrument,**parms)

    for x in artist:
        print 'wake'

    figs=artist.figs
    del artist

    figurenames=[f for f in figs]
    figurenames.sort()
    print figurenames


if __name__ == '__main__':
    main()
