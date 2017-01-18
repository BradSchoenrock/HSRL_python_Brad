import lg_base.core.array_utils as hau
import dplkit.role.decorator
import re,os
from datetime import datetime,timedelta
import lg_dpl_toolbox.filters.time_frame as time_frame

class basicARMTimeParser(object):
    def __init__(self):
        self.datematch=re.compile('\.[0-9]{8}\.[0-9]{6}\.')

    def __call__(self,fname):
        """ extract file start time from the filename

        :param fname: filename
        :return: datetime from file, or None if not discovered
        :type fname: str
        :rtype: datetime

        """

        m=self.datematch.search(fname)
        if not m:
            print fname , ' failed'
            return None
        tmp=m.group(0)
        tmp=tmp[:5]+'.'+tmp[5:]
        return datetime.strptime(tmp,'.%Y.%m%d.%H%M%S.')

def getBasedir(siteid):
        import lg_dpl_toolbox.core.archival as hru
        if isinstance(siteid,basestring):
            if os.access(siteid,os.F_OK):
                return siteid
            else:
                try:
                    tmp = hru.selectSource(instrument=siteid)
                except:
                    tmp = hru.selectSource(site=siteid)
        else:
            try:
                tmp = hru.selectSource(instrument=siteid)
            except:
                tmp = hru.selectSource(site=siteid)
        return tmp['Path']

class TimePath(object):
    def __init__(self,order):
        self.order=order

    def __call__(self,base,subdir,time):
        orders=('%4i' % time.year,'%02i' % time.month,'%02i' % time.day)
        parms=[base]
        parms.extend(orders[0:self.order])
        if subdir is not None:
            parms.append(subdir)
        return os.path.join(*parms)

    def startOf(self,start):
        if self.order==0:
            return start
        if self.order==1:
            return datetime(start.year,start.month,start.day,start.hour,start.minute,0)
        if self.order==2:
            return datetime(start.year,start.month,start.day,start.hour,0,0)
        if self.order==3:
            return datetime(start.year,start.month,start.day,0,0,0)
        if self.order==4:
            return datetime(start.year,start.month,1,0,0,0)
        if self.order==5:
            return datetime(start.year,1,1,0,0,0)
        raise RuntimeError("Order "+repr(self.order)+" unsupported")

    def dataDirDelta(self,time=None):
        if self.order!=0 and self.order!=3:
            raise NotImplementedError("Not implemented. need calendar")
        if self.order==0:
            return None#timedelta(days=365*100)
        if self.order==3:
            return timedelta(days=1)

    def nextTime(self,time):
        ordr=self.dataDirDelta(time)
        if ordr==None:
            return None
        return time+ordr

def basicARMPlatformParser(filename):
    if len(filename)>=3 and filename[0]!='.':
        return filename[:3]
    return None

@dplkit.role.decorator.exposes_attrs_of_field('src')
@dplkit.role.decorator.autoprovide(frameclass=hau.Time_Z_Group,reuseGenerator=False)
class ARMFileNarrator(dplkit.role.filter.aFilter):
    """ Narrator Initialzation.  Shouldn't be created

            :param src: file narrator
            :param zoo: zookeeper, if narrator should return data, not file info
            :param kwargs: addtional parameters are passed to thezookeeper on read, if provided
    """
    def __init__(self,src,zoo,preYield=None,**kwargs):
        super(ARMFileNarrator,self).__init__(src)
        self.src=src
        self.zoo=zoo
        self.init_preYield=preYield
        self.allkwargs=kwargs.copy()

    def preYield(self,x,attrs,found):#override this in subclasses
        return True

    def process(self,*args,**kwargs):
        fileOnly=kwargs.pop('fileOnly',None)
        allkwargs=self.allkwargs.copy()
        allkwargs.update(kwargs)
        for pfnr in self.src:
            x=self.zoo.open(self.zoo(pfnr),firsttime=self.src.start,lasttime=self.src.end,*args, **allkwargs)
            if x is not None:
                attrs=self.zoo.getAttributes(self.zoo(pfnr))
                found=self.zoo.getFoundVars(self.zoo(pfnr))
                if self.preYield(x,attrs,found) and (self.init_preYield is None or self.init_preYield(x,attrs,found)):
                    yield x

@dplkit.role.decorator.exposes_attrs_in_chain(['platform','arm_sourcepath','timeinfo','start','end'])
@dplkit.role.decorator.autoprovide(frameclass=dict,reuseGenerator=False)
class ARMNarrator(dplkit.role.narrator.aNarrator):
    """ Narrator Initialzation.  Shouldn't be created

            :param basedir: base directory to find data
            :param dataprefix: filename prefixes for data
            :param foldernames: folders
            :param start: start time
            :param end: end time
            :param timeparse: callable for parsing time
    """

    def __init__(self,basedir,instname,dataprefix,foldernames,start,end,timeparse=None,platform=None,platformparse=None,folderorder=3,debug=False):
        self.basedir=basedir
        self.foldernames=foldernames
        self.folderorder=TimePath(folderorder)
        self.timeinfo=time_frame.parse_timewindow(start,end)
        self._platform=platform
        if isinstance(foldernames,basestring):
            self.foldernames=[foldernames]
        if isinstance(dataprefix,basestring):
            dataprefix=[dataprefix]
        tmpdataprefix=[]
        for dp in dataprefix:
            for f in ('mf2','ns'):
                if dp.startswith(f+instname):
                    print 'USE of HSRL prefixes is deprecated. replace beginning of '+f+' with .*'
                    dp=dp.replace(f,'.*')
            for f in ('nsa','tmp','mag'):
                if dp.startswith(f+instname):
                    print 'USE of ARM prefixes is deprecated. replace beginning of '+f+' with .*'
                    dp=dp.replace(f,'.*')
            if dp not in tmpdataprefix:
                tmpdataprefix.append(dp)
        self.dataprefix=[]
        for p in tmpdataprefix:
            print 'prefix',p
            self.dataprefix.append(re.compile('^'+p))
        self.start=start
        self.end=end
        self.debug=debug or os.getenv('ARMNARR_DEBUG',None) is not None
        self.parseTimeFromFile=timeparse or basicARMTimeParser()
        self.starttime=self.folderorder.startOf(start)
        ordr=self.folderorder.dataDirDelta()
        if ordr is not None:
            self.starttime-=ordr
        if self._platform is None:
            self.platformparse=platformparse
            if self.platformparse is None:
                self.platformparse=basicARMPlatformParser
            if len(self.dataprefix)==1:
                self._platform=self.platformparse(self.dataprefix[0].pattern[1:])
 
    @property
    def arm_sourcepath(self):
        return self.basedir

    @property
    def platform(self):
        if self._platform is None:
            for f in self.read():
                self._platform=self.platformparse(f['filename'])
            print 'Platform for ARM is ',self._platform
        return self._platform

    def read(self):
        """ Narrator Generator function
        takes no parameters. initialization determines Narrator behavior
        parameters given here will be passed to the zookeeper, if narrator was given one
        """
        time=self.starttime
        fidx=-1
        pfnr=None
        while (self.end is None or time<=self.end) and time<=datetime.utcnow():
            if self.debug:
                print 'checking time ',time,self.end
            folder=None
            if self.foldernames is None:
                folder = self.folderorder(self.basedir,None,time)
            else:
                for f in self.foldernames:
                    folder=self.folderorder(self.basedir,f,time) #os.path.join(self.basedir,'%4i' % time.year,'%02i' % time.month,'%02i' % time.day,f)
                    if os.path.exists(folder):
                        break
            fidx=fidx+1
            try:
                if self.debug:
                    print 'listing path ' , folder
                flist=os.listdir(folder)
                flist.sort()
                if self.debug:
                    print flist
            except OSError:
                flist=[]
            idxoffile=-1
            fileindex=0
            while idxoffile<fidx and fileindex<len(flist):
                f=flist[fileindex]
                for pref in self.dataprefix:
                    if pref.search(f) is not None:
                        idxoffile=idxoffile+1
                fileindex=fileindex+1

            if idxoffile==fidx:
                fnr={'path':os.path.join(folder,f),'filename':f,'start':self.parseTimeFromFile(f),'width':timedelta(days=1)}
                if pfnr:
                    pfnr['width']=fnr['start']-pfnr['start']
                    if (self.end!=None and pfnr['start']>=self.end) or fnr['start']<=self.start:
                        pass
                    else:
                        yield pfnr
                pfnr=fnr
            else:
                fidx=-1
                time=self.folderorder.nextTime(time)# timedelta(days=1)
                if time is None:#no next
                    break
        if pfnr:
            if self.end:
                pfnr['width']=self.end-pfnr['start']
            else:
                pfnr['width']=datetime.utcnow()-pfnr['start']
            if (self.end!=None and pfnr['start']>=self.end) or (pfnr['start']+pfnr['width'])<=self.start:
                pass
            else:
                yield pfnr

