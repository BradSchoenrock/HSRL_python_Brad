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

def myGetch():
    fd = sys.stdin.fileno()

    oldterm = termios.tcgetattr(fd)
    newattr = termios.tcgetattr(fd)
    newattr[3] = newattr[3] & ~termios.ICANON & ~termios.ECHO
    termios.tcsetattr(fd, termios.TCSANOW, newattr)

    oldflags = fcntl.fcntl(fd, fcntl.F_GETFL)
    fcntl.fcntl(fd, fcntl.F_SETFL, oldflags | os.O_NONBLOCK)

    try:        
        while 1:            
            try:
                c = sys.stdin.read(1)
                break
            except IOError: pass
    finally:
        termios.tcsetattr(fd, termios.TCSAFLUSH, oldterm)
        fcntl.fcntl(fd, fcntl.F_SETFL, oldflags)

def extremeValue(fit,isMin=None,isMax=None,fitmetric=None,fitlimit=None):
    if fitlimit is not None:
        if fitmetric is None or fitmetric>=fitlimit:
            return np.NAN
    if isMin is None and isMax is None:
        return -fit[1]/(fit[0]*2)
    if isMin is None and isMax:
        if fit[0]<0:
            return extremeValue(fit)
        return np.NAN
    if isMax is None and isMin:
        if fit[0]>0:
            return extremeValue(fit)
        return np.NAN
    assert(False)

def getExtremeIndex(value,**kwargs):
    l=value.size
    r=np.arange(l)
    allx=np.polyfit(r,value, 2,full=True)
    fit=allx[0]
    stat=allx[1][0]
    #if stat>1.0:
    #    return None,stat
    try:
        centeridx=int(extremeValue(fit,fitmetric=stat,**kwargs))#np.where(value==np.min(value))[0][0]
    except ValueError:
        centeridx=-1
    if centeridx<0 or centeridx>=l:
        return None,stat
    return centeridx,stat

def getValueAtExtremeIndex(value,arr,**kwargs):
    x,st=getExtremeIndex(value,**kwargs)
    if x is None:
        raise IndexError("oops. didn't fit")
    return arr[x],x,st

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

@dplkit.role.decorator.autoprovide(frameclass=dict)
class EachCalScan(dplkit.role.filter.aFilter):
    def __init__(self,source):
        super(EachCalScan,self).__init__(source)
        self.source=source

    def process(self):
        content=dict()
        for tempContent in self.source:
            iscal=tempContent['op_mode'].astype('int32')&(2**0)
            nextStart=0
            while True:
                calidxs=np.where(iscal[nextStart:]!=0)[0]
                if len(calidxs)==0:
                    print('no cals')
                    break
                firstcalidx=calidxs[0]+nextStart
                calchangeidxs=np.where(iscal[firstcalidx:]==0)[0]
                if len(calchangeidxs)==0:
                    calrange=slice(firstcalidx,None)
                    haveall=False
                else:
                    calrange=slice(firstcalidx,calchangeidxs[0]+firstcalidx)
                    nextStart=calchangeidxs[0]+firstcalidx
                    haveall=True
                if firstcalidx!=0 and len(content)>0:
                    yield content
                    content=dict()
                for v,val in tempContent.items():
                    if v not in content:
                        content[v]=tempContent[v][calrange].copy()
                    else:
                        content[v].append(tempContent[v][calrange])
                if 'seedlaser_temp_to_freq' in self.hsrl_constants:
                    content['seedlaser_temp_to_freq']=self.hsrl_constants['seedlaser_temp_to_freq']
                if not haveall:
                    break
                yield content
                content=dict()
        if len(content)>0:
            yield ret


@dplkit.role.decorator.autoprovide(frameclass=dict)
class DropNarrowScans(dplkit.role.filter.aFilter):
    def __init__(self,source,minimumWidthHz):
        super(DropNarrowScans,self).__init__(source)
        self.source=source
        self.minimumWidthHz=minimumWidthHz

    def process(self):
        for content in self.source:
            if abs(np.max(content['interf_freq'])-np.min(content['interf_freq']))>=self.minimumWidthHz:
                yield content

@dplkit.role.decorator.autoprovide(frameclass=dict)
class IdentifyAdjusted(dplkit.role.filter.aFilter):
    def __init__(self,source):
        super(IdentifyAdjusted,self).__init__(source)
        self.source=source

    def process(self):
        endLastCal=datetime(2000,1,1,0,0,0)
        for ret in self.source:
            ret=ret.copy()
            if 'times' in ret:
                ret['isAdj']=(ret['times'][0]-endLastCal)<timedelta(minutes=5.0)
                endLastCal=ret['times'][-1]
            else:
                ret['isAdj']=(ret['s']-endLastCal)<timedelta(minutes=5.0)
                endLastCal=ret['e']
            yield ret

def operateOnScan(content):
        title=content['times'][0].strftime('%Y%m%dT%H%M%S')+' - '+content['times'][-1].strftime('%Y%m%dT%H%M%S')
        reccount=content['molecular_cal_pulse'].size
        recmin=int(reccount*(.5-.10))
        recmax=int(reccount*(.5+.10))
        mask_molec=slice(recmin,recmax)
        freqoffset=None
        ratio=None
        m_stat=0.0
        centerfreq=0.0
        lockstat=[]
        lockstatname=[]
        lockmask=[]
        if 'filtered_energy' in content:
            ratio=content['filtered_energy'][:,0]/content['nonfiltered_energy'][:,0]
            centeridx=np.argmin(ratio[mask_molec])+recmin
            lockstat.append(ratio)
            lockstatname.append("energy ratio")
            lockmask.append(mask_molec)
        else:
            centeridx,m_stat=getExtremeIndex(content['molecular_cal_pulse'][mask_molec],\
                    isMin=True,fitlimit=1.0)
            centeridx+=recmin
            lockstat.append(content['molecular_cal_pulse'])
            lockstatname.append("molec")
            lockmask.append(mask_molec)
        if 'superseedlasercontrollog' in content and 'seedlaser_temp_to_freq' in content and not np.all(np.isnan(content['superseedlasercontrollog'][:,7])):
            freq=content['superseedlasercontrollog'][:,7]*content['seedlaser_temp_to_freq']
            freq*=1e9
            print('using seed control for frequency')
            #print sidx,ridx,l3slope[idx],abs(seedoffset[sidx]-seedoffset[ridx]),abs(interfoffset[sidx]-interfoffset[ridx])
        else:
            freq=content['interf_freq']
        centerfreq=freq[centeridx]
        if 'l3locking_stats' in content:
            mask_l3=np.where(np.logical_and(freq>=(centerfreq-3e8),freq<=(centerfreq+3e8)))[0]
            #print(mask_l3)
            if len(mask_l3)>0:
                l3slope=content['l3locking_stats'][:,0].copy()
                l3slope[l3slope>1e10]=np.NAN
                centeridx=mask_l3[np.argmin(np.abs(l3slope[mask_l3]))]
                centerfreq=freq[centeridx]
                lockstat.append(l3slope)
                lockstatname.append("l3slope")
                lockmask.append(mask_l3)
        freqoffset=freq-centerfreq
        has1064= 'combined_1064_cal_pulse' in content and np.any(np.abs(content['combined_1064_cal_pulse'])>0.0025)
        if has1064:
            freqoffset_1064=freqoffset/2.0
        earlyplot=False
        if earlyplot:
            plt.figure()
            parms=[]
            parms.extend([freqoffset,content['combined_hi_cal_pulse'],'b'])
            parms.extend([freqoffset,content['molecular_cal_pulse'],'k'])
            if has1064:
                parms.extend([freqoffset_1064,content['combined_1064_cal_pulse'],'r'])
            plt.subplot(311);plt.plot(*parms)
            plt.title(title)
            for x in range(len(lockstat)):
                plt.subplot(300+x+1+10*len(lockstat));plt.plot(freqoffset,lockstat[x],'b',freqoffset[lockmask[x]],lockstat[x][lockmask[x]],'.b')
                plt.ylabel(lockstatname[x])
        freqrange_532=3e9
        freqrange_1064=5e9
        centerfreq_532,centeridx_532,_=getValueAtExtremeIndex(content['combined_hi_cal_pulse'],freqoffset)
        centerfreq_1064=np.NAN
        centeridx_1064=-1
        if has1064:
            centerfreq_1064,centeridx_1064,_=getValueAtExtremeIndex(content['combined_1064_cal_pulse'],freqoffset_1064)
        print('centers',centerfreq_532/1e9,centerfreq_1064/1e9)
        if has1064 and np.isnan(centerfreq_1064):
            print('Failed to get lock on 1064. dropped')
            has1064=False
        print('range',freqoffset[0],freqoffset[-1])
        #if abs(freqoffset[0]-freqoffset[-1])<1e9:
        #    return (content['times'][0],content['times'][-1],np.NAN,np.NAN)
        mask_532=np.logical_and(freqoffset>=(centerfreq_532-freqrange_532),freqoffset<=(centerfreq_532+freqrange_532))
        if not np.any(mask_532):
            mask_532=np.logical_not(np.isnan(freqoffset))
            return None
            raise IndexError("Bad 532")
        c532_stat=np.polyfit(freqoffset[mask_532],content['combined_hi_cal_pulse'][mask_532], 2,full=True)
        c532f=c532_stat[0]
        c532_stat=c532_stat[1][0]
        fit532=np.polyval(c532f, freqoffset[mask_532])
        print(c532f)
        peak_532=extremeValue(c532f,isMax=True,fitmetric=c532_stat,fitlimit=1.0)
        g_lock_energy=content['transmitted_energy'][centeridx]
        g_max_energy=content['transmitted_energy'][centeridx_532]
        g_max_calpulse=content['combined_hi_cal_pulse'][centeridx_532]
        g_lock_calpulse=content['combined_hi_cal_pulse'][centeridx]
        ir_max_energy=ir_lock_energy=ir_max_calpulse=ir_lock_calpulse=peak_1064=np.NAN
        if has1064:
            mask_1064=np.logical_and(freqoffset_1064>=(centerfreq_1064-freqrange_1064),freqoffset_1064<=(centerfreq_1064+freqrange_1064))
            c1064_stat=np.polyfit(freqoffset_1064[mask_1064],content['combined_1064_cal_pulse'][mask_1064], 2,full=True)
            c1064f=c1064_stat[0]
            c1064_stat=c1064_stat[1][0]
            fit1064=np.polyval(c1064f, freqoffset_1064[mask_1064])
            ir_lock_energy=content['transmitted_1064_energy'][centeridx]
            ir_max_energy=content['transmitted_1064_energy'][centeridx_1064]
            ir_max_calpulse=content['combined_1064_cal_pulse'][centeridx_1064]
            ir_lock_calpulse=content['combined_1064_cal_pulse'][centeridx]
            if c1064_stat>1.0:
                peak_1064=np.NAN
            else:
               peak_1064=extremeValue(c1064f,isMax=True,fitmetric=c1064_stat,fitlimit=1.0)
        if not earlyplot:
            plt.figure()
            parms=[]
            parms.extend([freqoffset,content['combined_hi_cal_pulse'],'b'])
            parms.extend([freqoffset,content['molecular_cal_pulse'],'k'])
            if has1064:
                parms.extend([freqoffset_1064,content['combined_1064_cal_pulse'],'r'])
            plt.subplot(311);plt.plot(*parms)
            plt.title(title)
            for x in range(len(lockstat)):
                plt.subplot(300+10*len(lockstat)+x+3);plt.plot(freqoffset,lockstat[x],'b',freqoffset[lockmask[x]],lockstat[x][lockmask[x]],'.b')
                plt.ylabel(lockstatname[x])
        label=[]
        parms=[]
        parms.extend([freqoffset[mask_532],fit532,'.b',[peak_532,peak_532],[0,.25],'b'])
        parms.extend([freqoffset,content['combined_hi_cal_pulse'],'b'])
        label.append('comb 532 center %g dev %g' % (peak_532,c532_stat))
        parms.extend([freqoffset[mask_molec],content['molecular_cal_pulse'][mask_molec],'.k',[0.0,0.0],[0,.25],'k'])
        parms.extend([freqoffset,content['molecular_cal_pulse'],'k'])
        #label.append('mol 532 center %g dev %g' % (0.0,m_stat))
        if has1064:
            parms.extend([freqoffset_1064[mask_1064],fit1064,'.r',[peak_1064,peak_1064],[0,.25],'r'])
            parms.extend([freqoffset_1064,content['combined_1064_cal_pulse'],'r'])
            label.append('comb 1064 center %g dev %g' % (peak_1064,c1064_stat))
        plt.subplot(313);plt.plot(*parms)
        plt.xlabel(','.join(label))
        return dict(s=content['times'][0],e=content['times'][-1],g=peak_532,ir=peak_1064
            ,g_max_energy=g_max_energy,g_lock_energy=g_lock_energy,ir_max_energy=ir_max_energy,ir_lock_energy=ir_lock_energy
            ,g_max_calpulse=g_max_calpulse,g_lock_calpulse=g_lock_calpulse,ir_max_calpulse=ir_max_calpulse,ir_lock_calpulse=ir_lock_calpulse)

@dplkit.role.decorator.autoprovide(frameclass=dict)
class EtalonDriftTrack(dplkit.role.filter.aFilter):
    def __init__(self,source):
        super(EtalonDriftTrack,self).__init__(source)
        self.source=source

    def process(self):
        for x in self.source:
            r=operateOnScan(x)
            if r is not None:
                yield r
            else:
                print('IS NONE')



class etalondriftgenerator(dplkit.role.librarian.aLibrarian):
    def __init__(self,inst):
        from hsrl.dpl.dpl_hsrl import dpl_hsrl
        self.dpl=dpl_hsrl(inst,filetype='calibration')

    def search(self,start,end,makeGraph=False):
        r=self.dpl(start_time_datetime=start,end_time_datetime=end,min_alt_m=0,max_alt_m=5000,
            forimage=False,raw_only=True,with_profiles=False)
        myvars=('times','interf_freq','filtered_energy','nonfiltered_energy','l3locking_stats','superseedlasercontrollog'
            ,'combined_hi_cal_pulse','combined_1064_cal_pulse','molecular_cal_pulse','op_mode'
            ,'transmitted_1064_energy','transmitted_energy')
        r=RawFieldSelection(r,myvars,subscope='rs_raw')
        r=EachCalScan(r)
        r=DropNarrowScans(r,6e9)
        r=EtalonDriftTrack(r)
        r=IdentifyAdjusted(r)
        return r

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

def appendContent(frame,outdict,alloptions,**kwargs):
    havenull='nullcontent' in kwargs
    if havenull:
        nullcontent=kwargs.pop('nullcontent')
    for k,attr in alloptions.items():
        if k not in outdict:
            outdict[k]=[]
        if frame is not None and attr in frame:
            outdict[k].append(frame[attr])
        elif not havenull:
            raise RuntimeError('MISSING VALUE '+attr)
        else:
            outdict[k].append(nullcontent)
            

def main(inst,startdate,enddate,sig=None):
    history=dict()
    parts=(('ir','g'),
           ('lock','max'),
           ('energy','calpulse'))
    allparts=makelist(parts)
    etalonparts=makelist((parts[0],),'adj')
    alwaysparts=makelist((('s','e'),))
    alwaysparts.update(makelist((parts[0],),'center'))
    began=datetime.utcnow()
    drift=etalondriftgenerator(inst)
    for fr in drift(startdate,enddate):
        if not fr['isAdj']:
            appendContent(fr,history,alwaysparts)
            appendContent(None,history,etalonparts,nullcontent=np.NAN)
            appendContent(fr,history,allparts)
        appendContent(None,history,allparts,nullcontent=np.NAN)
        appendContent(fr,history,etalonparts)
        appendContent(fr,history,alwaysparts)
    print('timed at '+repr(datetime.utcnow()-began))

    if 's' not in history:
        print('NO HISTORY')
        return
    #print(retlist)
    plt.figure()
    plt.title("Drift")
    starts=history['s']
    plt.plot(starts,history['center_g'],'+-b',starts,history['center_ir'],'o-r',[starts[0],starts[-1]],[0.0,0.0],'k')
    plt.plot(starts,history['adj_g'],'b',starts,history['adj_ir'],'r',linewidth=3)
    plt.grid()
    plt.figure()
    plt.title("Energy")
    plt.plot(starts,history['g_max_energy'],'ob',starts,history['ir_max_energy'],'or',
        starts,history['g_lock_energy'],'+b',starts,history['ir_lock_energy'],'+r')
    plt.legend(('532 Energy at Etalon Peak','1064 Energy at Etalon Peak',
        '532 Energy at Lock Point','1064 Energy at Lock Point'))
    plt.grid()
    plt.figure()
    plt.subplot(211)
    plt.title("Counts")
    plt.plot(starts,history['g_max_calpulse'],'ob',starts,history['g_lock_calpulse'],'+b')
    plt.legend(('532 CalPulse at Etalon Peak','532 CalPulse at Lock Point'))
    plt.grid()
    plt.subplot(212)
    plt.plot(starts,history['ir_max_calpulse'],'or',starts,history['ir_lock_calpulse'],'+r')
    plt.legend(('1064 CalPulse at Etalon Peak','1064 CalPulse at Lock Point'))
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
    print('press a key')
    myGetch()

