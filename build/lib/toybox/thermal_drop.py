from __future__ import print_function,absolute_import,division,unicode_literals
import matplotlib.pyplot as plt
from datetime import datetime,timedelta
import numpy as np
import sys,os
from bottleneck import nanmean
import threading

import sys    
import termios
import fcntl
from collections import namedtuple
import dplkit.role.filter
import dplkit.role.librarian

@dplkit.role.decorator.autoprovide(frameclass=dict)
class RawFieldSelection(dplkit.role.filter.aFilter):
    def __init__(self,source,varlist,subscope=None):
        super(RawFieldSelection,self).__init__(source)
        self.source=source
        self.subscope=subscope
        self.varlist=varlist

    def process(self):
        for f in self.source:
            print('loop')
            tempContent=dict()
            if self.subscope is not None:
                if not hasattr(f,self.subscope):
                    continue
                f=getattr(f,self.subscope)
            for v in self.varlist:
                if not hasattr(f,v):
                    print('no ',v)
                    continue
                if v not in tempContent:
                    tempContent[v]=getattr(f,v).copy()
                else:
                    tempContent[v].append(getattr(f,v))
            yield tempContent

def makelist(parts,prefix=None):
    alloptions=[]
    dropto=0
    if prefix is not None:
        alloptions.append(prefix)
        dropto=len(prefix)+1
    for p in parts:
        tmp=[]
        for x in p:
            if len(alloptions)==0:
                tmp.append(x)
            else:
                for pp in alloptions:
                    tmp.append(pp+'_'+x)
        alloptions=tmp
    ret=dict()
    for x in alloptions:
        ret[x]=x[dropto:]
    return ret

def appendContent(outdict,frame):
    for k,attr in frame.items():
        if k not in outdict:
            outdict[k]=attr.copy()
        else:
            outdict[k].append(attr)
            

def main(inst,startdate,enddate,sig=None):
    history=dict()
    parts=(('thermal1','thermal2'),
           ('records','goodrecords','errorCount'))
    alltherms=list(makelist(parts).keys())+['times']
    began=datetime.utcnow()
    import hsrl.dpl.dpl_hsrl as dpl_hsrl
    from lg_dpl_toolbox.dpl.NetCDFZookeeper import GenericTemplateRemapNetCDFZookeeper 
    zoo=GenericTemplateRemapNetCDFZookeeper(inst,forModule=dpl_hsrl,keepfields=alltherms)

    dpllib=dpl_hsrl.dpl_hsrl(inst,zoo=zoo)
    dpl=dpllib(start_time_datetime=startdate,end_time_datetime=enddate,min_alt_m=0,max_alt_m=500,
            forimage=False,raw_only=True,with_profiles=False)
    dpl=RawFieldSelection(dpl,alltherms,subscope='rs_raw')
    for fr in dpl:
        appendContent(history,fr)
    print('timed at '+repr(datetime.utcnow()-began))
    print(repr(history.keys()))
    #print(retlist)
    for x in alltherms:
        print(x)
        if 'errorCount' in x:
            history[x]=history[x].astype('float')
            history[x][history[x]<0.0]=-1
            history[x][history[x]>1e5]=-1
            print('%f %f' % (np.nanmin(history[x]),np.nanmax(history[x])))
        elif '_records' in x:
            history[x]=history[x].astype('float')
            history[x][history[x]<=0]=-1
            history[x][history[x]>1e7]=-1
        elif '_goodrecords' in x:
            history[x]=history[x].astype('float')
            history[x][history[x]<=0]=0
            history[x][history[x]>1e7]=0


    plt.figure()
    plt.title(inst+" error")
    starts=history['times']
    plt.plot(starts,history['thermal1_errorCount'],'b',
             starts,history['thermal2_errorCount'],'r')
    plt.legend(('thermal1','thermal2'))
    plt.grid()
    plt.figure()
    plt.title(inst+" badcount")
    plt.plot(starts,history['thermal1_records']-history['thermal1_goodrecords'],'b',
             starts,history['thermal2_records']-history['thermal2_goodrecords'],'r')
    plt.legend(('thermal1','thermal2'))
    plt.grid()
    print('with all graphs timed at '+repr(datetime.utcnow()-began))
    if sig is not None:
        with sig:
            sig.notifyAll()
    plt.show()

if __name__ == '__main__':
    p=os.path.realpath(os.path.join(os.path.dirname(sys.argv[0]),'..'))
    print(p)
    sys.path.append(p)
    inst=sys.argv[1]
    startdate=datetime.strptime(sys.argv[2],'%Y%m%d')#(2015,5,1,0,0,0)
    enddate=datetime.strptime(sys.argv[3],'%Y%m%d')#datetime(2015,5,18,0,0,0)
    if False:
        l=threading.Condition()
        with l:
            threading.Thread(target=main,args=(inst,startdate,enddate,l)).start()
            l.wait();
    else:
        main(inst,startdate,enddate)
    #print('press a key')
    #myGetch()

