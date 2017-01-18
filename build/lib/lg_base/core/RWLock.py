import fcntl
import threading
import os
from datetime import datetime,timedelta
import functools
#copied from picnic.RWLock

class EnterExitWrapper:
    def __init__(self,enter,exit,needexit=None):
        self.enter=enter
        self.exit=exit
        self.needexit=needexit
        self.exitres=None

    def __enter__(self):
        self.exitres=self.enter()
        return self.exitres

    def __exit__(self, type, value, traceback):
        if self.needexit==None or self.needexit(self.exitres):
            self.exit()
        self.exitres=None
        return False

class RWLock:

    def __enter__(self):
        return self.exclusiveLock()

    def __exit__(self, type, value, traceback):
        self.unlock()
        return False

    def __init__(self):
        self.cond=threading.Condition()
        self.waitingwritelocks=set()
        self.writelockThread=None
        self.readlockThreads=set()

    @staticmethod
    def justret(x):
        return x

    @property
    def shared(self):
        return EnterExitWrapper(self.sharedLock,self.unlock)

    def withTimeout(self,timeout):
        return EnterExitWrapper(functools.partial(self.exclusiveLock,timeout),self.unlock,RWLock.justret)

    def sharedWithTimeout(self,timeout):
        return EnterExitWrapper(functools.partial(self.sharedLock,timeout),self.unlock,RWLock.justret)

    def sharedLock(self,timeout=None):
        with self.cond:
            assert(not threading.currentThread() in self.readlockThreads)
            assert(threading.currentThread()!=self.writelockThread)
            if timeout!=None:
                stime=datetime.utcnow()
                etime=stime+timedelta(seconds=timeout)
            while self.writelockThread!=None and len(self.waitingwritelocks)>0:
                if timeout!=None:
                    n=datetime.utcnow()
                    if n>=etime:
                        return False
                    timeout=(etime-n).total_seconds()
                self.cond.wait(timeout=timeout)
            self.readlockThreads.add(threading.currentThread())
        return True

    def exclusiveLock(self,timeout=None):
        with self.cond:
            assert(not threading.currentThread() in self.readlockThreads)
            assert(threading.currentThread()!=self.writelockThread)
            self.waitingwritelocks.add(threading.currentThread())
            if timeout!=None:
                stime=datetime.utcnow()
                etime=stime+timedelta(seconds=timeout)
            while self.writelockThread!=None or len(self.readlockThreads)>0:
                if timeout!=None:
                    n=datetime.utcnow()
                    if n>=etime:
                        self.waitingwritelocks.remove(threading.currentThread())
                        return False
                    timeout=(etime-n).total_seconds()
                self.cond.wait(timeout=timeout)
            self.waitingwritelocks.remove(threading.currentThread())
            self.writelockThread=threading.currentThread()
        return True

    def unlock(self):
        with self.cond:
            if self.writelockThread==threading.currentThread():
                self.writelockThread=None
            else:
                self.readlockThreads.remove(threading.currentThread())
            self.cond.notifyAll()

class SharedRWLock:

    rwlockslock=threading.Lock()
    rwlocks=dict()

    @staticmethod
    def getLockFor(fn):
        if os.path.exists(fn):
            fn=os.path.realpath(fn)
        with SharedRWLock.rwlockslock:
            if not fn in SharedRWLock.rwlocks:
                SharedRWLock.rwlocks[fn]=RWLock()
            return SharedRWLock.rwlocks[fn]

    def __enter__(self):
        return self.exclusiveLock()

    def __exit__(self, type, value, traceback):
        self.unlock()
        return False

    def __init__(self,fn):
        self.fn=fn
        self.fi=file(fn,'a')
        self.lock=self.getLockFor(fn)
        self.fd=self.fi.fileno()
        self.sharedLockCount=0
        self.sharedPending=False
        self.sharedLockCountLock=threading.Condition()

    def __getfilelock(self,shared,nonblocking=False):
        if shared:
            with self.sharedLockCountLock:
                while self.sharedPending:
                    if nonblocking:
                        return False
                    self.sharedLockCountLock.wait()
                if self.sharedLockCount>0:
                    self.sharedLockCount+=1
                    return True
                self.sharedPending=True
        try:
            lv=fcntl.LOCK_SH if shared else fcntl.LOCK_EX
            if nonblocking:
                lv|=fcntl.LOCK_NB
            fcntl.flock(self.fd,lv)
        except IOError:#failed to lock
            if shared:
                with self.sharedLockCountLock:
                    self.sharedPending=False
                    self.sharedLockCountLock.notifyAll()
            return False
        if shared:
            with self.sharedLockCountLock:
                self.sharedPending=False
                self.sharedLockCount+=1
                self.sharedLockCountLock.notifyAll()
        return True

    def __getfilelockby(self,shared,bytime):
        while bytime>datetime.utcnow():
            if self.__getfilelock(shared=shared,nonblocking=True):
                return True
        return False

    @property
    def shared(self):
        return EnterExitWrapper(self.sharedLock,self.unlock)

    def withTimeout(self,timeout):
        return EnterExitWrapper(functools.partial(self.exclusiveLock,timeout),self.unlock,RWLock.justret)

    def sharedWithTimeout(self,timeout):
        return EnterExitWrapper(functools.partial(self.sharedLock,timeout),self.unlock,RWLock.justret)

    def sharedLock(self,timeout=None):
        if timeout!=None:
            deadline=datetime.utcnow()+timedelta(seconds=timeout)
        if not self.lock.sharedLock(timeout=timeout):
            return False
        if timeout==None:
            return self.__getfilelock(shared=True)
        else:
            if not self.__getfilelockby(shared=True,bytime=deadline):
                self.lock.unlock()
                return False
            return True

    def exclusiveLock(self,timeout=None):
        if timeout!=None:
            deadline=datetime.utcnow()+timedelta(seconds=timeout)
        if not self.lock.exclusiveLock(timeout=timeout):
            return False
        if timeout==None:
            ret=self.__getfilelock(shared=False)
            return ret
        else:
            if not self.__getfilelockby(shared=False,bytime=deadline):
                self.lock.unlock()
                return False
            return True

    def unlock(self):
        with self.sharedLockCountLock:
            if self.sharedLockCount>0:#this all expects that multiple shared locking is handled by the RWLock, as well as errors related to the write lock
                self.sharedLockCount-=1
            if self.sharedLockCount==0:
                fcntl.flock(self.fd,fcntl.LOCK_UN)
        self.lock.unlock()