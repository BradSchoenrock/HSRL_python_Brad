import os
from functools import wraps
import threading
import multiprocessing

def _killme():
    import signal
    os.kill(os.getpid(),signal.SIGKILL)

class wrapper:
    def __init__(self,func,seconds,error_message,doFork=False,onTimeout=None,killOnTimeout=False,noReturn=False):
        self.func=func
        self.seconds=seconds
        self.error_message=error_message
        self.onTimeout=onTimeout
        self.killOnTimeout=killOnTimeout
        self.forked=doFork
        self.noReturn=noReturn

    def handleTimeout(self,forkhost=False):
        if self.forked and not forkhost:
            _killme()
        print(self.error_message)
        try:
            if self.onTimeout is not None:
                if callable(self.onTimeout):
                    return self.onTimeout()
                else:
                    raise self.onTimeout
        finally:
            if self.killOnTimeout:
                _killme()
        raise RuntimeError(self.error_message)

    def _watchdogthread(self,ev):
        if not ev.wait(self.seconds):
            self.handleTimeout()

    def __call__(self,*args, **kwargs):
        if self.forked:
            return self.forkedcall(args,kwargs)
        return self.call(args,kwargs)

    def forkedcall(self,args,kwargs):
        ev=multiprocessing.Event()
        q=None
        if not self.noReturn:
            q=multiprocessing.Queue()
        ev.clear()
        bk=multiprocessing.Process(target=self.newprocesscall,args=(ev,q,args,kwargs))
        bk.start()
        bk.join()
        if ev.is_set():
            if q is not None:
                return q.get()
            return
        self.handleTimeout(forkhost=True);

    def newprocesscall(self,ev,q,args,kwargs):
        #ev.clear()
        res=self.call(args,kwargs,ev)
        if q is not None:
            q.put(res)
        #ev.set()

    def call(self,args,kwargs,_ev=None):
        result=None
        ev=_ev or threading.Event()
        x=threading.Thread(target=self._watchdogthread,args=[ev])
        ev.clear()
        x.start()
        try:
            result = self.func(*args, **kwargs)
        finally:
            ev.set()
        x.join()
        return result

def completion_timeout(seconds, error_message = 'Function call timed out',onTimeout=None,killOnTimeout=False):
    '''
    Decorator that provides timeout to a function
    '''
    def decorated(func):
        if os.getenv("NO_TIMEOUT",None) is not None:
            return func
        return wraps(func)(wrapper(func,seconds,error_message,doFork=False,onTimeout=onTimeout,killOnTimeout=killOnTimeout))
    return decorated

def completion_timeout_forked(seconds, error_message = 'Function call timed out',onTimeout=None,noReturn=True):
    '''
    Decorator that provides timeout to a function
    '''
    def decorated(func):
        if os.getenv("NO_TIMEOUT",None) is not None:
            return func
        return wraps(func)(wrapper(func,seconds,error_message,doFork=True,onTimeout=onTimeout,noReturn=noReturn))
    return decorated
