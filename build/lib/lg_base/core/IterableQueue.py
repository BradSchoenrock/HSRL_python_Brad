import Queue
import threading
import weakref
import os
import traceback
import multiprocessing
import logging
LOG = logging.getLogger(__name__)
from collections import namedtuple
import warnings

RaisedExceptionWrapper=namedtuple('RaisedExceptionWrapper','exceptionvalue')

VERBOSE = LOG.debug

myStopIteration='StopThisIterationStuff'

def currentContext():
    return repr(os.getpid())+":"+repr(threading.currentThread().ident)    

def isStopIteration(obj):
    """
    :param obj: check this object for StopIteration indication. used to send stop iteration over queues
    :returns: True if obj is a string token identifying stop iteration, or is a StopIteration object or instance
    """
    return (isinstance(obj,basestring) and obj==myStopIteration) or obj is StopIteration or isinstance(obj,StopIteration)

def delQueue(thread,timeout,queue,putqueue,doGet,doPut,forceGet,forcePut):
    if thread==None:
        return


    if timeout!=None:
        thread.join(timeout=timeout)
    if not thread.is_alive():
        return
    #if hasattr(thread,'terminate'):
    #    print 'terminating thread for '+repr(queue)
    #    thread.terminate()
    LOG.debug('Cleanup : Thread '+repr(thread)+' queue '+repr(queue)+('' if putqueue==None else (' putqueue '+ repr(putqueue)) ))
    if putqueue!=None:
        if forceGet:
            VERBOSE ('cleanup forced get on putqueue '+repr(putqueue))
            try:
                x=putqueue.get_nowait()
                VERBOSE ('GET Stop Iteration on queue. '+repr(putqueue))
                VERBOSE ('Flushed a '+repr(x))
                putqueue.task_done()
                VERBOSE ('Done on joining queue'+repr(putqueue))
            except Queue.Empty:
                pass
        if doPut:
            VERBOSE ('cleanup put on putqueue '+repr(putqueue))
            try:
                putqueue.put_nowait(myStopIteration)
                VERBOSE ('Put Stop Iteration on queue. joining '+repr(putqueue))
                #putqueue.join()
                VERBOSE ('Done on joining queue'+repr(putqueue))
            except Queue.Full:
                VERBOSE("Couldn't put stop iteration in queue "+repr(putqueue))
                pass
        if hasattr(putqueue,'close'):
            putqueue.close()
        if hasattr(putqueue,'cancel_join_thread'):
            putqueue.cancel_join_thread()
    if doGet:
        VERBOSE ('cleanup get on queue '+repr(queue))
        try:
            while True:
                x=queue.get_nowait()#trying to wake up the pushing thread, in case its blocking on a put()
                queue.task_done()
                VERBOSE ('*************flushed queue '+repr(queue)+' content '+repr(x))
        except Queue.Empty:
            pass
    if forcePut:
        VERBOSE ('cleanup forced put on queue '+repr(queue))
        try:
            queue.put_nowait(myStopIteration)#trying to wake up the pushing thread, in case its blocking on a put()
            VERBOSE ('put on queue '+repr(queue))
            #print ('*************flushed queue '+repr(queue)+' content '+repr(x))
        except Queue.Full:
            VERBOSE("couldn't put stop iteration in queue "+repr(queue))
            pass
    VERBOSE ('joining thread '+repr(thread)+" for queue "+repr(queue))
    if False and hasattr(thread,'pid') and thread.pid!=None:
        os.system('kill -9 '+repr(thread.pid))
        VERBOSE('cleanup on thread '+repr(thread)+' is killing process '+repr(thread.pid))
    else:
        thread.join()#timeout=30.0)

class cleanupQueue(object):
    def __init__(self):
        self.mythreads=[]

    def cleanup(self,thread,queue,putqueue=None,doGet=False,doPut=False,forceGet=False,forcePut=False,timeout=5.0):
        t=threading.Thread(target=delQueue,args=(thread,timeout,queue,putqueue,doGet,doPut,forceGet,forcePut),name='Cleanup Thread for '+repr(thread))
        t.start()
        self.mythreads.append(t)

    def __del__(self):
        for t in self.mythreads:
            t.join()#timeout=30.0)

class WrappedQueue(object):
    def __init__(self,queue,wrapper,unwrapper):
        self.queue=queue
        self.wrapper=wrapper
        self.unwrapper=unwrapper

    def __repr__(self):
        return 'Wrapped Queue ['+repr(self.queue)+']'

    def put(self,_obj,*args,**kwargs):
        obj=_obj
        if self.wrapper!=None:
            obj=self.wrapper(obj)
        self.queue.put(obj,*args,**kwargs)

    def put_nowait(self,_obj,*args,**kwargs):
        obj=_obj
        if self.wrapper:
            obj=self.wrapper(obj)
        self.queue.put_nowait(obj,*args,**kwargs)

    def get(self,*args,**kwargs):
        r=self.queue.get(*args,**kwargs)
        if self.unwrapper:
            r=self.unwrapper(r)
        return r

    def get_nowait(self,*args,**kwargs):
        r=self.queue.get_nowait(*args,**kwargs)
        if self.unwrapper:
            r=self.unwrapper(r)
        return r

    def join(self,*args,**kwargs):
        self.queue.join(*args,**kwargs)        

    def task_done(self,*args,**kwargs):
        self.queue.task_done(*args,**kwargs)

    def close(self,*args,**kwargs):
        if hasattr(self.queue,'close'):
            self.queue.close(*args,**kwargs)
    def join_thread(self,*args,**kwargs):
        if hasattr(self.queue,'join_thread'):
            self.queue.join_thread(*args,**kwargs)
    def cancel_join_thread(self,*args,**kwargs):
        if hasattr(self.queue,'cancel_join_thread'):
            self.queue.cancel_join_thread(*args,**kwargs)

class SimpleValue(object):
    """ Drop in replacement for a Multiprocessing Value object without the multiprocessing
    
    :param initval: initial value
    """
    def __init__(self,initval):
        self.value=initval

def doNothing(x):
    return x

class QueueCountObj(object):
    def __init__(self,queue):
        self.queue=queue
        self.context=currentContext()
        self.localGetCount=0
        if True:#isinstance(queue,type(multiprocessing.Queue())):#,type(multiprocessing.JoinableQueue()))):
            #raise RuntimeError('Multiprocessing queue')
            self.lock=multiprocessing.Lock()
            self.inqueue=multiprocessing.Value('i',0)
            self.tasksout=multiprocessing.Value('i',0)
            self.tasksdone=multiprocessing.Value('i',0)
        else:
            self.lock=threading.Lock()
            self.inqueue=SimpleValue(0)
            self.tasksout=SimpleValue(0)
            self.tasksdone=SimpleValue(0)

    def __repr__(self):
        return 'Counting Queue ['+repr(self.queue)+']'

    def put(self,obj,*args,**kwargs):
        VERBOSE (currentContext()+' Putting queue '+repr(self))
        if isStopIteration(obj):
            VERBOSE ('This is a stop iteration')
        self.queue.put(obj,*args,**kwargs)
        VERBOSE (currentContext()+' Finished Putting queue '+repr(self))
        with self.lock:
            self.inqueue.value+=1
            self.printtasks_nolock()

    def put_nowait(self,obj,*args,**kwargs):
        VERBOSE (currentContext()+' putting (nonblock) on queue ' + repr(self))
        try:
            self.queue.put_nowait(obj,*args,**kwargs)
        except Queue.Empty:
            VERBOSE (currentContext()+' FAILED putting (nonblock) queue '+repr(self))
            raise
        VERBOSE (currentContext()+' Finished putting (nonblock) queue '+repr(self))
        with self.lock:
            self.inqueue.value+=1
            self.printtasks_nolock()

    def get(self,*args,**kwargs):
        VERBOSE (currentContext()+' Getting queue '+repr(self))
        if self.context!=currentContext():
            self.context=currentContext()
            self.localGetCount=0
        r=self.queue.get(*args,**kwargs)
        VERBOSE (currentContext()+' Finished getting queue '+repr(self))
        with self.lock:
            self.localGetCount+=1
            self.inqueue.value-=1
            self.tasksout.value+=1
            self.printtasks_nolock()
        return r

    def get_nowait(self,*args,**kwargs):
        if self.context!=currentContext():
            self.context=currentContext()
            self.localGetCount=0
        VERBOSE (currentContext()+' Getting (nonblock) queue '+repr(self))
        try:
            r=self.queue.get_nowait(*args,**kwargs)
        except Queue.Empty:
            VERBOSE (currentContext()+' FAILED Getting (nonblock) queue '+repr(self))
            raise
        VERBOSE (currentContext()+' finished Getting (nonblock) queue '+repr(self))
        with self.lock:
            self.localGetCount+=1
            self.inqueue.value-=1
            self.tasksout.value+=1
            self.printtasks_nolock()
        return r

    def join(self,*args,**kwargs):
        VERBOSE (currentContext()+' Joining Queue '+repr(self))
        self.printtasks()
        self.queue.join(*args,**kwargs)
        VERBOSE (currentContext()+' Joined queue '+repr(self))

    def task_done(self,*args,**kwargs):
        VERBOSE (currentContext()+' Task Done on queue ' + repr(self))
        #print self.context+'='+repr(self.localGetCount)+' on '+currentContext()
        #assert(self.localGetCount>0)
        #assert(self.context==currentContext())
        self.queue.task_done(*args,**kwargs)
        VERBOSE (currentContext()+' finished task_done on queue ' + repr(self))
        with self.lock:
            self.localGetCount-=1
            self.tasksout.value-=1
            self.tasksdone.value+=1
            self.printtasks_nolock()

    def close(self,*args,**kwargs):
        VERBOSE (currentContext()+' Close on queue ' + repr(self))
        if hasattr(self.queue,'close'):
            self.printtasks()
            self.queue.close(*args,**kwargs)
    def join_thread(self,*args,**kwargs):
        VERBOSE (currentContext()+' join_thread on queue ' + repr(self))
        if hasattr(self.queue,'join_thread'):
            self.printtasks()
            self.queue.join_thread(*args,**kwargs)
        VERBOSE (currentContext()+' finished join_thread on queue ' + repr(self))

    def cancel_join_thread(self,*args,**kwargs):
        VERBOSE (currentContext()+' cancel_join on queue ' + repr(self))
        if hasattr(self.queue,'cancel_join_thread'):
            self.printtasks()
            self.queue.cancel_join_thread(*args,**kwargs)

    def printtasks_nolock(self):
        if self.tasksdone<20:
            VERBOSE (repr(self)+' has '+repr(self.inqueue.value)+' inqueue, '+repr(self.tasksout.value)+' in progress, '+repr(self.tasksdone.value)+' completed')

    def printtasks(self):
        with self.lock:
            self.printtasks_nolock()        

    def __del__(self):
        VERBOSE (currentContext()+' Deleting '+repr(self))
        self.printtasks()
        VERBOSE ('Printed status. (lock is ok) '+repr(self))

QueueCount=  QueueCountObj# doNothing# Obj

class QueueAsIterable(object):
    """ Convert a Joinable Queue object into an iterable object
    
    :param queue: the queue to iterate over
    :param inSync: dead parameter
    :param atAbort: optional callable if the iterating generator is prematurely broken from the caller (downstream). typically to signal the other side of the Queue object to stop.
    :param atExit: optional callable to be called at the natural termination point
    :param atStart: optional callable to be called just as the iterator begins iterating
    """
    def __init__(self,queue,inSync=False,atAbort=None,atExit=None,atStart=None
            ,nonblocking=False,voidtoken=None):
        self.queue=queue
        self.atExit=atExit
        self.atAbort=atAbort
        self.inSync=inSync
        self.atStart=atStart
        self.lock=multiprocessing.Lock()
        self.nonblocking=nonblocking
        self.voidtoken=voidtoken
        #self.started=False
        self.pid=None

    def __repr__(self):
        return 'Iterable Queue ['+repr(self.queue)+']'

    def wouldblock(self):
        return self.queue.empty()

    @property
    def blocking(self):
        return not self.nonblocking

    def setblocking(self,bl,voidtoken='placeholder'):
        if voidtoken != 'placeholder':
            self.voidtoken=voidtoken
        self.nonblocking=not bl

    def __iter__(self):
        VERBOSE ('Using iterable generator for queue '+repr(self.queue))
        #if self.started:
        #    print "Can't use a QueueAsIterable object twice!"
        #assert(not self.started)
        doException=None
        wasError=False
        with self.lock:
            #self.started=True
            self.pid=os.getpid()
            needTaskDoneException=False
            normalExit=False
            doAbort=False
            doRun=True
            if self.atStart!=None:
                doRun=self.atStart()
                if doRun is None:
                    doRun=True
            try:
                while doRun:
                    try:
                        f=self.queue.get(block=not self.nonblocking)
                    except Queue.Empty:
                        yield self.voidtoken
                        continue
                    except:# (KeyboardInterrupt,IOError):
                        VERBOSE ('Queue '+repr(self.queue)+' was closed')
                        needTaskDoneException=False
                        raise
                    else:
                        if isStopIteration(f):
                            normalExit=True
                            return
                        elif isinstance(f,RaisedExceptionWrapper):
                            self.queue.task_done()
                            raise f.exceptionvalue
                        #print repr(f)+' from queue '+repr(queue)
                        else:
                            yield f
                            self.queue.task_done()
            except GeneratorExit as e:
                needTaskDoneException=False #broke at yield
                doAbort=True
                doException=e
                try:
                    self.queue.task_done()
                    f=self.queue.get(block=False)
                    if isStopIteration(f):
                        doAbort=False
                    self.queue.task_done()
                except Queue.Empty:
                    pass
            except Exception, e:
                VERBOSE ('Iterable generator exception on queue '+repr(self.queue))
                traceback.print_exc()
                doAbort=True
                doException=e
                wasError=True
            finally:
                if doAbort:
                    self.doAbort("Generator Exit" if not wasError else "Exception Exit",wasError,needTaskDoneException)
                if self.atExit is not None:
                    self.atExit()
                self.atAbort=None
                if normalExit:
                    self.queue.task_done()
        if doException is not None:
            raise doException

    def doAbort(self,description,isError=False,doTaskDone=False):
                if self.atAbort is not None:
                    self.atAbort()
                    self.atAbort=None
                if doTaskDone:
                    self.queue.task_done()
                #no flush needed since the exception should be from the queue anyway, and the last message from it
                try:
                    while self.queue is not None:
                        x=self.queue.get_nowait()
                        self.queue.task_done()
                        if isError:
                            warnings.warn(description+' lead to Flushed '+repr(x)+' from queue '+repr(self.queue)+'.... THIS SHOULDN"T HAPPEN EVER')
                        else:
                            VERBOSE(description+' lead to Flushed '+repr(x)+' from queue '+repr(self.queue))                            
                except Queue.Empty:
                    pass

    def __del__(self):
        if os.getpid()!=self.pid:
            return
        VERBOSE ('deleting queue as iterable '+repr(self))
        self.doAbort("Deletion")
        r=repr(self)
        #if not self.inSync:
        #    self.queue.task_done() ##FIXME why isn't this needed?
        if hasattr(self.queue,'close'):
            VERBOSE('closing '+r)
            self.queue.close()
            VERBOSE('closed '+r)
        if hasattr(self.queue,'cancel_join_thread'):
            VERBOSE('cancel_join_thread '+r)
            self.queue.cancel_join_thread()
            VERBOSE('cancel_join_threaded '+r)
        VERBOSE('clearing '+r)
        self.queue=None
        VERBOSE('cleared '+r)


class reprSet(object):
    def __init__(self,*args):
        self.args=args

    def __repr__(self):
        ret=''
        for x in self.args:
            if isinstance(x,basestring):
                ret+=x
            else:
                ret+=repr(x)
        return ret

def trythreadmain(*args,**kwargs):
        try:
            threadmain(*args,**kwargs)
        except:
            import traceback
            traceback.print_exc()
            raise


def threadmain(iterableref,queue,inSync,terminator,finisher,threadtypename):
        iterable=iterableref[0]
        if False:
            del iterableref[0]
            assert(len(iterableref)==0)
        #name=repr(os.getpid())+":"+repr(threading.currentThread().ident%1000)+' Threaded Using Iterable ['+repr(iterable)+'] to queue ['+repr(queue)+']'
        name=reprSet(os.getpid(),":",'%i' % (threading.currentThread().ident),' Threaded Using Iterable [',iterable,'] to queue [',queue,']')
        VERBOSE ( repr(name)+' queue '+threadtypename+' started THD')
        normalExit=True
        iterations=0
        if terminator.is_set():
            VERBOSE ( repr(name)+' queue '+threadtypename+' about to rapid exit THD')
            normalExit=False
        elif finisher.is_set():
            VERBOSE ( repr(name)+' queue '+threadtypename+' about to rapid exit THD')
        else:
            try:
                realiterable=iter(iterable)
                for i in realiterable:
                    iterations+=1
                    if terminator.is_set():
                        normalExit=False
                        break
                    queue.put(i)#,timeout=30.0)
                    if inSync:
                        if terminator.is_set():
                            normalExit=False
                            break
                        queue.join()
                    if terminator.is_set():
                        normalExit=False
                        break
                    if finisher.is_set():
                        break
                if False and normalExit==False:
                    try:
                        x=realiterable.next()
                    except StopIteration:
                        #print 'Almost exited abnormally by terminator'
                        normalExit=True
            except Exception as e:
                traceback.print_exc()
                queue.put(RaisedExceptionWrapper(e))#,timeout=30.0)
                normalExit=False
        if hasattr(realiterable,'close'):
            realiterable.close()
        del iterable
        if normalExit:
            VERBOSE (repr(name)+' queue '+threadtypename+' about to normal exit after '+repr(iterations)+' steps THD') 
            #if not terminator.is_set():
            queue.put(myStopIteration)
            #if hasattr(queue,'close'):
            #    queue.close()
                    #print (repr(name)+' queue put successfully. Now join thread')
            if inSync:# and not terminator.is_set():
                    #print (repr(name)+' thread joined')
                #print (repr(name)+' queue put successfully. Now join')
                queue.join()
                if False and hasattr(queue,'join_thread'):
                    queue.close()
                    queue.join_thread()
            #elif hasattr(queue,'cancel_join_thread'):
            #    queue.cancel_join_thread()
            elif hasattr(queue,'cancel_join_thread'):
                    queue.close()
                    queue.join_thread()
                    queue.cancel_join_thread()
            else:
                #print (repr(name)+' not joining. insync='+repr(inSync)+' terminated='+repr(terminator.is_set()))
                pass
        else:
            VERBOSE (repr(name)+' queue '+threadtypename+' about to abnormal exit after '+repr(iterations)+' steps THD')
            #queue.put(myStopIteration)
            if hasattr(queue,'close'):
                queue.close()
            #if hasattr(queue,'join_thread'):
            #    queue.join_thread()
            if hasattr(queue,'cancel_join_thread'):
                 queue.cancel_join_thread()
        VERBOSE (repr(name)+' queue '+threadtypename+' exited THD')

class IterableToQueueThread(object):
    """Creates a thread to push the flow of an Iterable object into a joinable queue object

    :param iterable: the iterable source
    :param queue: the queue destination
    :param inSync: if True, after every put(), the queue will be join()ed, waiting for the consumer to indicate task_done(). this means this thread will only operate when the consumer is idle, and usually waiting for content
    :param threadclass: thread class to use when creating this instance. Typically a threading.Thread or a multiprocessing.Process.
    :param terminator: Event class object (from multiprocessing or threading) that, when set, will cause the thread to rapidly exit
    :param finisher: Event class object that, when set, will cause the thread to cleanly exit.
    """
    def __init__(self,iterable,queue,inSync=False,threadclass=None,terminator=None,finisher=None):
        #mythreadclass= threadclass or threading.Thread
        self.event=terminator or multiprocessing.Event()
        self.finisher=finisher or multiprocessing.Event()
        self.iterableref=[iterable]
        self.queue=queue
        self.inSync=inSync
        self.threadclass=threadclass
        self.multiprocessing=None
        self.thread=None
        #self.pid=None

    def __repr__(self):
        return 'IterableToQueueThread('+repr(self.ident)+')['+repr(self.thread)+' with queue '+repr(self.queue)+']'

    @property
    def ident(self):
        """ identity of the thread or process """
        if not hasattr(self.thread,'ident'):
            return None
        return self.thread.ident

    @property
    def pid(self):
        """ pid of the object if multiprocessing, or None """
        if self.multiprocessing:
            return self.thread.pid
        return None

    def getQueueAsIterable(self,**kwargs):
        """ Create an Iterable that mirrors this one, fed by the queue given at initialization time

        :returns: iterable object
        """
        return QueueAsIterable(self.queue,inSync=self.inSync,atExit=self.finish,atAbort=self.terminate
            ,atStart=self.start if self.inSync else None,**kwargs)

    def is_alive(self):
        """ :returns: true if the servicing thread is running """
        return self.thread!=None and self.thread.is_alive()

    def start(self):
        assert(not self.is_alive())
        self.event.clear()
        self.finisher.clear()
        self.thread= self.threadclass(target=trythreadmain,args=(self.iterableref,self.queue,self.inSync,self.event,self.finisher,repr(self.threadclass)),name="IterableToQueue for "+repr(self.iterableref)+' to '+ repr(self.queue))
        self.multiprocessing=hasattr(self.thread,'pid')
        self.thread.start()
        self.parentpid=os.getpid()

    def join(self,*args,**kwargs):
        VERBOSE ('Joining '+repr(self.thread))
        if self.thread.is_alive():
            self.thread.join(*args,**kwargs)
        VERBOSE ('Finished joining '+repr(self.thread))

    def finish(self):
        #assert(self.parentpid==os.getpid())
        VERBOSE (repr(self)+' is being finished in context '+currentContext())
        self.finisher.set()

    def terminate(self):
        #assert(self.parentpid==os.getpid() or self.pid==os.getpid())
        self.finish()
        VERBOSE (repr(self)+' is being terminated in context '+currentContext())
        self.event.set()
        found=False
        try:
            while False:
                x=self.queue.get_nowait()
                self.queue.task_done()
                #found=True
                VERBOSE ('**********'+repr(self)+' flushed '+repr(x))
        except Queue.Empty:
            pass
        if found:
            try:
                self.queue.put_nowait(myStopIteration)
                #self.queue.join()
            except Queue.Full:
                pass
        if False and hasattr(self.thread,'terminate'):
            self.thread.terminate()
        if False and self.pid!=None:
            VERBOSE ('terminate is killing process '+repr(self.pid)+' from terminate of iterable to queue thread')
            os.system('kill -9 '+repr(self.pid))

    def __del__(self):
        if self.parentpid!=None and self.parentpid==os.getpid():
            if self.is_alive():
                VERBOSE ('Terminating on dealloc')
                self.terminate()
            self.join()#timeout=30.0)