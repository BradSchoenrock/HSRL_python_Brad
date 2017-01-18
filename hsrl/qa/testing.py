from datetime import datetime,timedelta
from collections import OrderedDict,namedtuple

def safeenc(d):
    if isinstance(d,datetime):
        return d.strftime('%Y.%m.%d %H:%M:%S')
    if isinstance(d,dict):
        r=OrderedDict()
        for i,v in d.items():
            r[i]=safeenc(v)
        return r
    if isinstance(d,list):
        r=[]
        for v in d:
            r.append(safeenc(v))
        return r
    return d

def deepgetattr(f,x):
    if not isinstance(x,list):
        return deepgetattr(f,x.split('.'))
    if len(x)==0:
        return f
    return deepgetattr(getattr(f,x[0]),x[1:])

def main():
    from dpl_filters import QAFlagClonedSyncFilter
    from dpl_narrators import QCFlagNarrator,QCLogNarrator
    import dpl_zookeeper as qazoo
    import dpl_librarian as qalib
    #from parsing import fileParser,flagParser
    from hsrl.dpl.dpl_hsrl import dpl_hsrl
    import lg_dpl_toolbox.filters.time_frame as time_slicing
    import sys
    import json

    td=timedelta(hours=2)
    et=datetime.utcnow()
    if len(sys.argv)>3:
        st=datetime.strptime(sys.argv[3],'%Y.%m.%dT%H:%M:%S')
    else:
        st=et-td
    if len(sys.argv)>4:
        et=datetime.strptime(sys.argv[4],'%Y.%m.%dT%H:%M:%S')   

    from lg_dpl_toolbox.dpl.TimeSource import TimeGenerator

    hsrllib=dpl_hsrl(instrument=sys.argv[2])

    import lg_base.core.canvas_info as ci
    canvas_info=ci.load_canvas_info()
    #process_control=hsrllib.hsrl_process_control
    number_x_pixels = canvas_info['canvas_pixels']['x']
    timesource=TimeGenerator(start_time=st,end_time=et,time_step_count=number_x_pixels)
    zookeeper=qazoo.QualityAssuranceZookeeper()#responsible for reading the files
    librarian=qalib.QualityAssuranceLibrarian(instrument=sys.argv[2])#responsible for finding the files

    hsrlnar=hsrllib(min_alt_m=0,max_alt_m=20000,timesource=timesource)

    #f=fileParser()
    #fp=flagParser(hsrlnar.altitudeAxis,binwidth=hsrlnar.hsrl_constants['binwidth'] * 1.5e8)
    #x= f.parseFile(sys.argv[1])
    x=zookeeper(librarian(start_time=st,end_time=et))#finding the files and reading them all
    print x
    json.dump(safeenc(x),file('dump.json','w'),indent=4,separators=(',', ': '))
    lognarr=QCLogNarrator(x,st,et)
    print 'Log Content for',st,'to',et,'is:'
    for e in lognarr:
        print e['header']['date'],':',e['content'].strip()
    #assert(0)
    qcnarr=QCFlagNarrator(x,altitude_axis=hsrlnar.altitudeAxis,binwidth=hsrlnar.hsrl_constants['binwidth']*1.5e8)#narrator streams the read flags as a single timestep per entry
    import lg_dpl_toolbox.filters.substruct as frame_substruct

    qcnarr=QAFlagClonedSyncFilter(qasource=qcnarr,timealtsource=timesource,timename='start')#takes the output of the narrator, and echos entries for each timestep (and potentially project range entries to altitude)
    if False:
        hsrlnarsplitter=frame_substruct.SubstructBrancher(hsrlnar)

        #timesource=time_slicing.TimeGinsu(hsrlnarsplitter.narrateSubstruct('rs_inv'),'times',onlyTime=True)#,isEnd=True)
        #qcnarr=frame_substruct.Retyper(qcnarr,hau.Time_Z_Group,dict(timevarname='times',altname='altitudes'))
        qcnarr=frame_substruct.CountDeGinsu(frame_substruct.FrameLength(hsrlnarsplitter.narrateSubstruct('rs_inv'),'times'),qcnarr)
        dplc=frame_substruct.NestingCompositer(hsrlnarsplitter.narrateSubstruct(None),dict(rs_hsrl_qc=qcnarr))
        printkeys=('rs_inv','rs_hsrl_qc','rs_hsrl_qc.qcflags.shape','rs_hsrl_qc.qcflags')
    else:
        dplc=qcnarr
        printkeys=('qcflags.shape','qcflags')

    dplc=time_slicing.FrameCachedConcatenate(dplc)

    for f in dplc:
        print vars(f).keys()
        for x in printkeys:
            print x
            print deepgetattr(f,x)

if __name__ == '__main__':
    main()