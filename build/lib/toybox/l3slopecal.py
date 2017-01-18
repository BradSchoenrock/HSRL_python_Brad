import matplotlib.pyplot as plt
from datetime import datetime,timedelta
import numpy as np
import sys,os
from bottleneck import nanmean

def makeL3Cal(inst,start,end,graphs=None,usel3slopecenter=True):
    from hsrl.dpl.dpl_hsrl import dpl_hsrl
    dpl=dpl_hsrl(inst)
    r=dpl(start_time_datetime=start,end_time_datetime=end,min_alt_m=0,max_alt_m=5000,forimage=False,raw_only=True,with_profiles=False)
    myvars=('times','interf_freq','filtered_energy','nonfiltered_energy','l3cavityvoltage','l3locking_stats','superseedlasercontrollog','molecular_cal_pulse')
    content=dict()
    for x in r:
        for v in myvars:
            if not hasattr(x.rs_raw,v):
                print 'no ',v
            if v not in content:
                content[v]=getattr(x.rs_raw,v).copy()
            else:
                content[v].append(getattr(x.rs_raw,v))
    interfoffset=content['interf_freq']
    seedoffset=content['superseedlasercontrollog'][:,7]*r.hsrl_constants['seedlaser_temp_to_freq']
    seedoffset*=1e9
    ratio=content['filtered_energy'][:,0]/content['nonfiltered_energy'][:,0]
    l3slope=content['l3locking_stats'][:,0].copy()
    ridx=np.argmin(ratio)
    sidx=np.argmin(np.abs(l3slope))
    if usel3slopecenter:
        idx=sidx
    else:
        idx=ridx
    print sidx,ridx,l3slope[idx],abs(seedoffset[sidx]-seedoffset[ridx]),abs(interfoffset[sidx]-interfoffset[ridx])
    seedoffset-=seedoffset[idx]
    interfoffset-=interfoffset[idx]

    pc = np.polyfit(l3slope,seedoffset, 5)
    #pc[-1]=0.0
    if graphs is not None:
        #plt.subplot(411);plt.plot(content['times'],ratio)
        if 2 in graphs:
            plt.subplot(412);plt.plot(content['times'],seedoffset,'b',content['times'],interfoffset,'r')
        if 3 in graphs:
            plt.subplot(413);plt.plot(content['times'],content['l3locking_stats'][:,0],'b')
        if 1 in graphs:
            plt.subplot(411);plt.plot(seedoffset,ratio,'b')#content['molecular_cal_pulse'])
        if 4 in graphs:
            seedfit = np.polyval(pc, l3slope)
            plt.subplot(414);plt.plot(seedoffset,l3slope,'b',seedfit,l3slope,'g')
    print 'fit',pc
    return pc
    plt.show()
    #[ -7.59976307e-04   1.53739948e-01   6.58691640e+01  -1.07624030e+04 -4.34018180e+06   1.88870179e+06]   20150212JPG


def showL3Shift(inst,start,end,cal):
    from hsrl.dpl.dpl_hsrl import dpl_hsrl
    dpl=dpl_hsrl(inst,filetype='data')
    r=dpl(start_time_datetime=start,end_time_datetime=end,min_alt_m=0,max_alt_m=5000,forimage=False,raw_only=True,with_profiles=False)
    myvars=('times','interf_freq','filtered_energy','nonfiltered_energy','l3cavityvoltage','l3locking_stats','superseedlasercontrollog','molecular_cal_pulse')
    content=dict()
    for x in r:
        for v in myvars:
            if not hasattr(x.rs_raw,v):
                print 'no ',v
            if v not in content:
                content[v]=getattr(x.rs_raw,v).copy()
            else:
                content[v].append(getattr(x.rs_raw,v))
    if len(content.keys())==0:
        return None
    interfoffset=content['interf_freq']
    seedoffset=content['superseedlasercontrollog'][:,7]*r.hsrl_constants['seedlaser_temp_to_freq']
    ratio=content['filtered_energy'][:,0]/content['nonfiltered_energy'][:,0]
    l3slope=content['l3locking_stats'][:,0].copy()
    #plt.subplot(411);plt.plot(content['times'],ratio)
    plt.subplot(412);plt.plot(content['times'],seedoffset,'b',content['times'],interfoffset,'r')
    plt.subplot(413);plt.plot(content['times'],l3slope,'b')
    ret=np.polyval(cal, l3slope)
    plt.subplot(414);plt.plot(content['times'],ret,'g')
    #plt.subplot(414);plt.plot(seedoffset,l3slope,seedfit,l3slope)
    #plt.show()
    return ret

def makeTable(inst,start,end,cal,averagewindow,filename):
    v=[]
    step=start
    while step<end:
        res=showL3Shift(inst,step,step+averagewindow,cal)
        if res is not None:
            v.insert(0,(step,nanmean(res)))
        step=step+averagewindow
    v.insert(0,(step,0))
    f=file(filename,'w')

    f.write('# i2 offset correction table\n')
    coef='[' + ( ','.join([('%g' % x) for x in cal ])) + ']'
    f.write('# using calibration '+coef+'\n')
    f.write('\ni2offset(GHz)\n')
    for t,v in v:
        f.write(t.strftime('    %d-%b-%Y %H:%M , [')+('%f' % v)+']\n')
    del f

if __name__ == '__main__':
    p=os.path.realpath(os.path.join(os.path.dirname(sys.argv[0]),'..'))
    print p
    sys.path.append(p)
    makeL3Cal('mf2hsrl',datetime(2015,2,12,17,14,30),datetime(2015,2,12,17,24,30),graphs=[1],usel3slopecenter=False)
    pc=makeL3Cal('mf2hsrl',datetime(2015,2,12,17,18,30),datetime(2015,2,12,17,21,30))#,graphs=[2,3,4])
    #vals=showL3Shift('mf2hsrl',datetime(2015,1,28,0,0,0),datetime(2015,1,28,1,0,0),pc)
    startdate=datetime(2015,1,1,0,0,0)
    enddate=datetime(2015,2,1,0,0,0)
    makeTable('mf2hsrl',startdate,enddate,pc,timedelta(hours=6),
        'i2offsettable.i2off'+startdate.strftime('_s%Y%m%dT%H%M%S')+enddate.strftime('_e%Y%m%dT%H%M%S')+datetime.utcnow().strftime('_c%Y%m%dT%H%M%S'))
    #showL3Shift('mf2hsrl',datetime(2015,1,12,0,0,0),datetime(2015,1,13,0,0,0),pc)
    plt.show()


