import threading, Queue
import multiprocessing
import dplkit.role.filter
from datetime import datetime,timedelta
import copy
import traceback
import logging
import weakref
import os

LOG = logging.getLogger(__name__)

import cPickle as pickle

import lg_base.core.IterableQueue as iq
from collections import namedtuple


class taggedwrapper(object):
    def __init__(self,field,name):#wrapper needs loads and dumps
        self.field=field
        self.name=name

    def tag(self,obj,tagval):
        tagval=('(PID %i Thr %i) ' % (os.getpid(),threading.currentThread().ident%1000)) +tagval
        if isinstance(obj,dict):
            if self.field not in obj:
                obj[self.field]=tagval
            else:
                obj[self.field]+=', then '+tagval
        else:
            try:
                if not hasattr(obj,self.field):
                    setattr(obj,self.field,tagval)
                else:
                    setattr(obj,self.field,getattr(obj,self.field)+', then '+tagval)
            except:
                pass
        return obj

    def dumps(self,obj):
        ret=self.tag(obj,'Dumped on '+self.name)
        return ret

    def loads(self,obj):
        ret=self.tag(obj,'Loaded on '+self.name)
        if True:
            pass
        elif isinstance(ret,dict):
            print ret[self.field]
        else:
            try:
                print getattr(ret,self.field)
            except:
                pass
        return ret

ThreadRecord=namedtuple('ThreadRecord','queue thread pid')

class aThreading_filter(dplkit.role.filter.aFilter):
    """Thread based filter. pushes the host stream to another thread from the reader

    *** WARNING This filter is extremely experimental, and can cause hangs or unexpected results because of the asynchronous behavior and process separation

    :param stream: source stream to put on another thread
    :param flow: queue flow amount. if 0, will always operate fully synced
    :param freeflow: if set True, flow will be set to infinite, allowing the source to run unhindered of consumption rate
    :param name: optional descriptive name
    :param fullsync: if specified, will force a fullsync (blocks source until after consumer is done with frame) or not (cache data as it comes in, without blocking unless flow amount is hit)
    """
    def __init__(self,stream,flow=5,freeflow=False,name=None,fullsync=None,*args,**kwargs):
        super(aThreading_filter,self).__init__(stream)
        self.stream=stream

        if fullsync==None:
            fullsync=(os.getenv('SYNC')=='YES')
        if fullsync:
            LOG.debug('FULLSYNC')
            flow=0  
            freeflow=False

        self.flow=flow
        self.name=name or "unnamed"
        self.queues=[]
        self.inSync=False
        #self.eventclass=threading.Event
        self.threadclass=threading.Thread
        self.queueclass=Queue.Queue
        #self.event=None

        if (freeflow or self.flow<0):
            self.flow=0
        elif self.flow==0:
            self.inSync=True
            self.flow=1

    def __repr__(self):
        return 'Threading Filter<'+repr(self.stream)+'>'

    def newthread(self):
        #if self.event==None:
        #    self.event=self.eventclass()
        inqueue=self.queueclass(self.flow)
        #hostwrapper=taggedwrapper('transport_tag',self.name+' host')
        #threadwrapper=taggedwrapper('transport_tag',self.name+' thread')
        #inqueue=iq.WrappedQueue(inqueue,pickle.dumps,pickle.loads)
        #inqueue=iq.WrappedQueue(inqueue,threadwrapper.dumps,threadwrapper.loads)
        #inqueue=iq.ThreadingQueue(inqueue,terminator=self.event)
        #inqueue=iq.WrappedQueue(inqueue,hostwrapper.dumps,hostwrapper.loads)
        inqueue=iq.QueueCount(inqueue)
        thread=iq.IterableToQueueThread(self.stream,inqueue,inSync=self.inSync,threadclass=self.threadclass)#,terminator=self.event)
        return ThreadRecord(inqueue,thread,os.getpid())

    def process(self,hostObject=False,**kwargs):
        args=self.newthread()

        self.queues.append(args)

        if not self.inSync:
            args.thread.start() #start now, rather than when the iterable is operated on

        ret=args.thread.getQueueAsIterable(**kwargs)#,atDone=args[1].join)
        if not hostObject:
            ret=iter(ret)
        LOG.debug('QueueAsIterable in '+self.name+' for queue '+repr(args.queue)+' is '+repr(ret))
        return ret
          
    def __del__(self):
        #if self.event!=None:
        #    self.event.set()
        queues=copy.copy(self.queues)
        del self.queues
        c=iq.cleanupQueue()
        LOG.debug (self.name+' at del stage on pid '+repr(os.getpid()))
        for record in queues:
            if record.pid!=os.getpid():
                continue
            LOG.debug ('Killing thread for '+self.name+' queue '+repr(record.queue))

            c.cleanup(record.thread,record.queue)#,forceGet=True)
        del queues
        LOG.debug (self.name+' finishing del stage')
        del c
        LOG.debug (self.name+' finished del stage')

class forking_filter(aThreading_filter):
    """ Filter that puts the source framestream's operation on another process entirely, using multiprocessing.JoinableQueue for transport

    *** WARNING This filter is extremely experimental, and can cause hangs or unexpected results because of the asynchronous behavior and process separation
    """
    def __init__(self,*args,**kwargs):
        super(forking_filter,self).__init__(*args,**kwargs)
        #self.eventclass=multiprocessing.Event
        self.threadclass=multiprocessing.Process
        self.queueclass=multiprocessing.JoinableQueue
 
    def __repr__(self):
        return 'Forking Filter<'+repr(self.stream)+'>'

class wait_filter(dplkit.role.filter.aFilter):
    def __init__(self,stream,waitfunc,maxduration=None):
        super(wait_filter,self).__init__(stream)
        assert(isinstance(stream,aThreading_filter))
        self.stream=stream
        self.waitfunc=waitfunc
        self.maxduration=maxduration

    def process(self):
        token="EMPTY"
        iv=self.stream.process(hostObject=True)
        lastupdate=datetime.utcnow()
        for f in iv:
            if f is not token:
                yield f
                if iv.blocking:
                    iv.setblocking(False,voidtoken=token)
                if self.maxduration is not None:
                    lastupdate = datetime.utcnow()
            else:
                if self.maxduration is not None and datetime.utcnow()-lastupdate>self.maxduration:
                    iv.setblocking(True)
                else:
                    print "@@@  IN WAITING",datetime.utcnow()-lastupdate
                    self.waitfunc()
                    print "@@@ OUT WAITING"
                    #iv.setblocking(True)
 
QueueEventTrio=namedtuple('QueueEventTrio','queue event finish')
ThreadQueuePair=namedtuple('ThreadQueuePair','thread queue')

def tryKeepingThread(*args,**kwargs):
        try:
            keepingThread(*args,**kwargs)
        except:
            import traceback
            traceback.print_exc()
            raise


def keepingThread(myself,iterable,cmdq,queuelist,termqueues,event):
    idx=0
    while True:
        if len(termqueues)>idx:
            termqueues[idx].thread.join()#timeout=30.0)
        #if myself()==None:
        #    break
        cmdq.put(idx)
        cmdq.join()
        if myself()==None:
            #cmdq.close()
            break
        event.acquire()
        qe=queuelist[idx]
        thread=iq.IterableToQueueThread(iterable,qe,inSync=True,threadclass=threading.Thread,terminator=qe.event,finisher=qe.finish)
        if len(termqueues)>idx:
            termqueues[idx]=ThreadQueuePair(thread,qe)
        else:
            termqueues.append(ThreadQueuePair(thread,qe))
        thread.start()
        event.notify()
        event.release()
        idx=(idx+1)%len(queuelist)


class KeepUpstreamOnThisProcess(dplkit.role.filter.aFilter):
    """ force a framestream back to the process this object was created on

    *** WARNING This filter is extremely experimental, and can cause hangs or unexpected results because of the asynchronous behavior and process separation

    :param stream: source stream to keep on this process. Downstream may be pushed to another process via a forking filter, but upstream will stay here.
    :param maxInstances: number of management objects to create, which are kept across processes and used to transport data across processes
    """
    def __init__(self,stream,maxInstances=3):
        super(KeepUpstreamOnThisProcess,self).__init__(stream)
        self.queues=[]
        self.stream=stream
        self.pid=os.getpid()
        for x in range(maxInstances):
            self.queues.append(QueueEventTrio(iq.QueueCount(multiprocessing.JoinableQueue(1)),multiprocessing.Event(),multiprocessing.Event()))
        self.cmdQueue=multiprocessing.JoinableQueue(1)
        self.cmdEvent=threading.Condition()
        self.termqueues=[]
        self.thread=threading.Thread(target=tryKeepingThread,args=(weakref.ref(self),stream,self.cmdQueue,self.queues,self.termqueues,self.cmdEvent),name='KeepUpstream on this Process Thread for '+repr(os.getpid())+' using queue '+repr(self.cmdQueue)+ ' and stream '+repr(stream))
        self.thread.start()

    def __repr__(self):
        return 'KeepUpstreamOnThisProcess '+repr(self.pid)+' Filter<'+repr(self.stream)+'>'

    def process(self):
        iterable=None
        ev=None
        if os.getpid()==self.pid:
            LOG.debug( 'KUOTP is on same process. easy')
            iterable=self.stream
        else:
            LOG.debug( 'KUOTP is on different process. loading a queue as the iterable')
            self.cmdEvent.acquire()
            qid=self.cmdQueue.get(timeout=1.0)
            self.cmdQueue.task_done()
            self.cmdEvent.wait()
            qe=self.queues[qid]#termqueues[qid]
            iterable=iq.QueueAsIterable(qe.queue,inSync=True,atAbort=qe.event.set,atExit=qe.finish.set)
            LOG.debug( 'QueueAsIterable in KUOTP for queue '+repr(qe.queue)+' is '+repr(iterable))
            self.cmdEvent.release()
        return iter(iterable)
        #for f in iterable:
        #    yield f


    def __del__(self):
        if self.pid!=os.getpid():
            return
        c=iq.cleanupQueue()
        for q in self.termqueues:#only runs on original pid
            c.cleanup(q.thread,q.queue)#,forcePut=True)
        del c
        self.cmdQueue.get(timeout=1.0)
        self.cmdQueue.task_done()
        self.thread.join()#timeout=30.0)


def safeFork(_stream,*args,**kwargs):
    """ allow a stream to fork to a separate process, but only from the constructing process. this keeps all processes as immediate children of the host process

    *** WARNING This filter is extremely experimental, and can cause hangs or unexpected results because of the asynchronous behavior and process separation
    """
    stream=forking_filter(_stream,*args,**kwargs)
    stream=KeepUpstreamOnThisProcess(stream)
    return stream

def donothing(stream,*args,**kwargs):
    return stream
threading_filter=None
import os
c=os.getenv('THREADING')
if c=='NONE':
    threading_filter=donothing
    LOG.debug ('NO THREADING')
elif c=='THREAD':
    threading_filter=aThreading_filter
    LOG.debug ('SINGLE PROCESS')
elif c=='SAFE':
    threading_filter= safeFork # aThreading_filter# donothing #forking_filter# aThreading_filter# safeFork#aThreading_filter# forking_filter
    LOG.debug ('SAFE MULTI PROCESS')
else:#if c=='SAFE':
    threading_filter= forking_filter# aThreading_filter# donothing #forking_filter# aThreading_filter# safeFork#aThreading_filter# forking_filter
    #LOG.debug 'MULTI PROCESS'
