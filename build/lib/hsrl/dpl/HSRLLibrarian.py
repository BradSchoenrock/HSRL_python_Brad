import dplkit.role.librarian
import dplkit.role.narrator
import dplkit.role.decorator
from datetime import datetime,timedelta
import os
import time
import re
import sys
import plistlib
import lg_base.core.array_utils as hau
import lg_dpl_toolbox.core.archival as arc

@dplkit.role.decorator.autoprovide(frameclass=hau.Time_Z_Group)
class HSRLRawNarrator(dplkit.role.narrator.aNarrator):
    def __repr__(self):
        return 'HSRLRawNarrator(%s)' % (self.host)

    def __init__(self,host,zoo,**kwargs):
        super(HSRLRawNarrator,self).__init__(host)
        self.host=host
        self.zoo=zoo
        self.allkwargs=kwargs.copy()


    def read(self,*args,**kwargs):

        allkwargs=self.allkwargs.copy()
        allkwargs.update(kwargs)

        for r in self.host:
            x=self.zoo.open(self.zoo(r),firsttime=self.host.start,lasttime=self.host.end,*args,**allkwargs)
            if x!=None and x.times.size>0:
                yield x


@dplkit.role.decorator.autoprovide(frameclass=dict)
class HSRLFileNarrator(dplkit.role.narrator.aNarrator):
    """ HSRL Librarian Search Result object. Shouldn't be created, should only be used by the HSRLLibrarian

    :param host: source librarian object
    :param basedir: base directory to find HSRL data
    :param dataprefix: filename prefixes for HSRL data
    :param start: start time
    :param end: end time
    :param filetype: file suffix to be included
    :param zoo: zookeeper, if narrator should return data, not file info
    :param kwargs: addtional parameters are passed to thezookeeper on read, if provided

    """

    def folderForDate(self,atime):
        """ folder path for a date

        :param atime: time for the folder. only date-resolution is used
        :return: string path of the folder containing the data for the time given
        :type atime: datetime
        :rtype: str

        """

        basefolder=os.path.join(self.basedir,'%4i' % atime.year,'%02i' % atime.month,'%02i' % atime.day)
        possibleFolders=['rawnc4','zraw','raw']
        for f in possibleFolders:
            folder=os.path.join(basefolder,f)
            if os.access(folder,os.R_OK):
                break
        return folder

    def listForDate(self,atime):
        """ list datafiles for a given date window. for HSRL, the date window is a whole day

        :param atime: time for the day of data to return
        :return: list of filenames
        :type atime: datetime
        :rtype: list of str

        """

        folder=self.folderForDate(atime)
        if not os.access(folder,os.R_OK):
            return []
        try:
            flist=os.listdir(folder)
            flist.sort()
            r=[]
            for f in flist:
                if f.startswith(self.dataprefix) and (self.filetype==None or ('_'+self.filetype) in f):
                    r.append(f)
            return r
        except OSError:
            return []

    def nextDate(self, atime):
        """ Get the next relevant storage window time. for HSRL data, this is 1 day

        :param atime: time
        :return: next time interval moment
        :type atime: datetime
        :rtype: datetime

        """

        ret=atime+timedelta(days=1)
        #if ret>datetime.utcnow():
        #    time.sleep(5)
        #    return atime
        return ret

    def parseTimeFromFile(self,fname):
        """ extract file start time from the filename

        :param fname: filename
        :return: datetime from file, or None if not discovered
        :type fname: str
        :rtype: datetime

        """

        m=self.host.datematch.search(fname)
        if not m:
            print fname , ' failed'
            return None
        tmp=m.group(0)
        tmp=tmp[:5]+'.'+tmp[5:]
        return datetime.strptime(tmp,'_%Y.%m%dT%H%M%S_')

    def __repr__(self):
        return 'HSRLRawNarrator(%s,%s,%s-%s,%s)' % (self.host,self.basedir,self.start,self.end,self.filetype)

    def __init__(self,host,basedir,dataprefix,start,end,filetype):
        super(HSRLFileNarrator,self).__init__(host)
        self.host=host
        self.basedir=basedir
        self.dataprefix=dataprefix
        self.start=start
        self.end=end
        self.filetype=filetype
        self.debug = False
        if hasattr(host,'site'):
            for w in host.site['Windows']:
                if self.end!=None and self.end<w['Start']:
                    continue
                if 'End' in w and w['End']<self.start:
                    continue
                self.window=w
                if self.start<w['Start']:
                    self.start=w['Start']
                if 'End' in w and (self.end==None or self.end>w['End']):
                    self.end=w['End']

        maxIntervalsBack=14
        self.starttime=datetime(self.start.year,self.start.month,self.start.day,0,0,0)
        self.startidx=0
        idx=-1
        while idx<0 and maxIntervalsBack>0:
            for f in self:
                t=f['start']
                if t!=None and t>self.start:
                    break
                idx=idx+1
            if idx<0:
                self.starttime=self.starttime-timedelta(days=1)
                maxIntervalsBack=maxIntervalsBack-1
        if idx<0:
            idx=0
        self.startidx=idx
        #print 'search for ', self.start, ' starts at ', self.starttime, ':', idx
 

    def read(self):
        """ Narrator Generator function
        takes no parameters. initialization determines Narrator behavior
        parameters given here will be passed to the zookeeper, if narrator was given one
        """
        time=self.starttime
        fidx=self.startidx-1
        flist=None
        folder=None

        shouldsend=False
        priorfile=None
        priortime=None
        priorfolder=None

        while self.end==None or time<=self.end:
            fidx=fidx+1
            if folder==None or flist==None:
                folder=self.folderForDate(time)
                flist=self.listForDate(time)
            if self.debug:
                print "HSRLRawNarrator.read()"
                print folder
                print time
                print fidx
                print len(flist)
                print ('None' if (fidx<0 or fidx>=len(flist)) else flist[fidx]) 
                print flist
            if fidx<len(flist) and fidx>=0:
                t=self.parseTimeFromFile(flist[fidx])
                if priorfile!=None:
                    dur=t-priortime
                    if dur.total_seconds()>(60*60*1.1) and 'data' in priorfile:
                        dur=timedelta(seconds=60*60*1.1)
                    r={'path':os.path.join(priorfolder,priorfile),'filename':priorfile,'start':priortime,'width':dur}
                    if self.debug:
                        print 'HSRLRawNarrator.read() yielding ', r
                    yield r
                priorfile=flist[fidx]
                priortime=t
                priorfolder=folder
                if self.end!=None and t>self.end:
                    #print 'Done'
                    return
                shouldsend=True
            else:
                newtime=self.nextDate(time)
                #print newtime
                if self.end!=None or newtime<datetime.utcnow():
                    #print datetime.utcnow(), ' next folder'
                    fidx=fidx-len(flist)-1
                    time=newtime
                else:
                    #print datetime.utcnow(), ' now. stop'
                    if priorfile!=None:
                        dur=datetime.utcnow()-priortime
                        if dur.total_seconds()>(60*60*1.1) and 'data' in priorfile:
                            dur=timedelta(seconds=60*60*1.1)
                        r={'path':os.path.join(priorfolder,priorfile),'filename':priorfile,'start':priortime,'width':dur}
                        yield r
                    #return
                    #this code is if the last should repeatedly return
                    if len(flist)==0:  
                        fidx=-1
                        print 'start from beginning'
                    else:
                        fidx=fidx-2
                    print 'start from %i of this folder with %i files' % (fidx+1,len(flist))
                folder=None
        if priorfile!=None:
            dur=self.end-priortime
            if dur.total_seconds()>(60*60*1.1) and 'data' in priorfile:
                dur=timedelta(seconds=60*60*1.1)
            r={'path':os.path.join(priorfolder,priorfile),'filename':priorfile,'start':priortime,'width':dur}
            yield r


class HSRLLibrarian(dplkit.role.librarian.aLibrarian):
    """ HSRL Raw Data Search

    Example:
    ``lib=HSRLLibrarian(dataset='gvhsrl')``

    :param instrument: instrument name (uses ``HSRL_CONFIG``)
    :param site:       site id (index into data archive list by site deployment)
    :param dataset:    dataset id (index or short name in data archive list by instrument)
    :param dataarchive_path: specifiy location of data archive list, default is in environemnt ``HSRL_DATA_ARCHIVE_CONFIG``, or ``/etc/dataarchive.plist``

    A DPL Search Generator. Should be used with a Zookeeper to access data (FIXME should have a Narrator wrapper to actually read in the data. Zookeeper finds/makes it available on disk)
    simple usecase example (creates a csv of laser temperatures) ::

        fields=['times','superseedlasercontrollog','laserpowervalues']    
        st=datetime(2013,01,28,21,35,0)
        et=datetime(2013,01,28,22,15,0)
        lib=HSRLLibrarian(dataset='bagohsrl')
        zoo=GenericTemplateRemapNetCDFZookeeper(datatype,keepfields=fields)
        m=lib.search(start=st,end=et)#,filetype='data')
        
        outf=file('lasertemps_'+ st.strftime('%Y%m%d_%H%M%S') +'_'+ et.strftime('%Y%m%d_%H%M%S') +'.csv','w')
        outf.write('time,current,voltage,LDD_Temp,2HG_Temp,LC_ksd,LD_ksd,Amb\\n')
        for f in m:
            res=zoo(uri=f)
            r=zoo.open(res)
            for i in range(0,r.times.shape[0]):
                if r.times[i]<st or r.times[i]>et:
                    continue
                outf.write('%f,%f,%f\\n' % ((r.times[i]-st).total_seconds(),r.superseedlasercontrollog[i,4],r.superseedlasercontrollog[i,5]))
        outf.close()

    Object is reusable, as the search function (callable) doesn't modify this object

    """

    def __init__(self,filetype=None,zoo=None,**kwargs):
        super(HSRLLibrarian,self).__init__()
        self.zoo=zoo
        self.datematch=re.compile('_[0-9]{8}T[0-9]{6}_')
        self.filetype=filetype
        import lg_dpl_toolbox.core.archival as hru
        site=hru.selectSource(**kwargs)
        for i in site['Instruments']:
            p=i.lower()
            if 'hsrl' in p:
                self.dataprefix=p.split('-')[0]
                break
        self.basedir=site['Path']
        self.name=site['Name']
        if 'SiteID' in site:
            self.siteid=site['SiteID']
            self.site=site
        if 'DatasetID' in site:
            self.dataset=site['DatasetID']

    def __repr__(self):
        if hasattr(self,'site'):
            return 'HSRLLibrarian(site:%s(%s),%s,%s)' % (self.siteid,self.site['Name'],self.basedir,self.dataprefix)
        if hasattr(self,'dataset'):
            return 'HSRLLibrarian(dataset:%s(%s),%s,%s)' % (self.dataset,self.name,self.basedir,self.dataprefix)
        return 'HSRLLibrarian(instrument:%s,%s)' % (self.dataprefix,self.basedir)

    def search(self,start,end=None, *args, **kwargs):
        """ Conduct a search of this library

        :param start: start time
        :type start: datetime
        :param end: end time or None to continue to now
        :type end: datetime
        :param filetype: file suffix. default None is all, other common possibilities are 'data' or 'calibration'
        :type filetype: str
        :rtype: HSRLSearchResult (iterable)

        """

        if len(args):
            print 'Unused args = ',args
        if len(kwargs):
            print "Unused kwargs = ",kwargs
        ft=self.filetype
        if 'filetype' in kwargs:
            ft=kwargs['filetype']
            del kwargs['filetype']
        ret=HSRLFileNarrator(self,self.basedir,self.dataprefix,start,end,ft)
        if self.zoo is not None:
            ret= HSRLRawNarrator(ret,zoo=self.zoo,*args,**kwargs)
        return ret

if __name__=='__main__':
    from lg_dpl_toolbox.dpl.NetCDFZookeeper import GenericTemplateRemapNetCDFZookeeper 
    datatype='ahsrl' if len(sys.argv)<2 else sys.argv[1]
    et=datetime.utcnow() if len(sys.argv)<4 else datetime.strptime(sys.argv[3],'%Y%m%dT%H%M%S')
    st=(et-timedelta(days=.5)) if len(sys.argv)<3 else datetime.strptime(sys.argv[2],'%Y%m%dT%H%M%S')
    fields=['times','telescope_position','telescope_rotation','telescope_rotation_measured','telescope_elevation','telescope_accelerometer_raw']#,'superseedlasercontrollog','laserpowervalues']
    lib=HSRLLibrarian(instrument=datatype)#site=16)#None,datatype)
    zoo=GenericTemplateRemapNetCDFZookeeper(datatype,keepfields=fields)
    m=lib(start=st,end=et)#,filetype='data')

    outf=file('telescopeaccel_'+ st.strftime('%Y%m%d_%H%M%S') +'_'+ et.strftime('%Y%m%d_%H%M%S') +'.csv','w')
    outf.write('seconds,pos,rot,arot,aele,ax,ay,az\n')#current,voltage,LDD_Temp,2HG_Temp,LC_ksd,LD_ksd,Amb\n')
    for f in m:
        print 'from librarian:',f
        continue
        res=zoo(uri=f)
        print 'uri from zoo:',res
        r=zoo.open(res)
        print f
        print r.times
        print r
        for i in range(0,r.times.shape[0]):
            if r.times[i]<st or r.times[i]>et:
                continue
            outf.write('%f,%f,%f,%f,%f,%f,%f,%f\n' % ((r.times[i]-st).total_seconds(),
            r.telescope_position[i],r.telescope_rotation[i],
            r.telescope_rotation_measured[i],r.telescope_elevation[i],
            r.telescope_accelerometer_raw[i,0],r.telescope_accelerometer_raw[i,1],r.telescope_accelerometer_raw[i,2])
            ) 
    outf.close()
