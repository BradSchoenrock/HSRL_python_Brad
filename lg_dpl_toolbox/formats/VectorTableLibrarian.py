import calendar
import dplkit.role.decorator
import dplkit.role.narrator
import dplkit.role.librarian
import re,os
from datetime import datetime,timedelta

class basicTimeParser(object):
    def __init__(self):
        self.datematch=re.compile('[0-9]{8}T[0-9]{4}')

    def __call__(self,fname,deftime,info):
        """ extract file start time from the filename

        :param fname: filename
        :return: datetime from file, or None if not discovered
        :type fname: str
        :rtype: datetime

        """

        if info is not None and 'hasDate' in info and info['hasDate']==False:
            return deftime
        m=self.datematch.search(fname)
        if not m:
            print fname , ' failed'
            return None
        tmp=m.group(0)
        tmp=tmp[:4]+'.'+tmp[4:]
        try:
            return datetime.strptime(tmp,'%Y.%m%dT%H%M')
        except:
            print 'Malformed date '+tmp+' vs '+'%Y.%m%dT%H%M'
            return None

def calendarDays(date):
    _,dc=calendar.monthrange(date.year,date.month)
    return timedelta(days=dc)

@dplkit.role.decorator.autoprovide(frameclass=dict)
class VectorTableNarrator(dplkit.role.narrator.aNarrator):
    """ Narrator Initialzation.  Shouldn't be created

            :param basedir: base directory to find data
            :param dataprefix: filename prefixes for data
            :param datasuffix: filename suffix for data
            :param foldernames: folders
            :param start: start time
            :param end: end time
            :param timeparse: callable for parsing time
            :param kwargs: addtional parameters are passed to thezookeeper on read, if provided
    """

    def __init__(self,basedir,instname,info,start,end,timeparse=None,completeList=False,alternativePath=None,subdir=None):
        self.basedir=basedir if alternativePath is None else os.path.join(alternativePath,instname)
        self.info=info
        self.start=start
        self.end=end
        self.subdir=subdir
        self.parseTimeFromFile=timeparse or basicTimeParser()
        self.starttime=datetime(start.year,start.month,1,0,0,0)
        self.completeList=completeList

    def read(self):
        """ Narrator Generator function
        takes no parameters. initialization determines Narrator behavior
        parameters given here will be passed to the zookeeper, if narrator was given one
        """
        time=self.starttime
        fidx=-1
        pfnr=None
        while (self.end is None or time<=self.end) and time<=datetime.utcnow():
            folder=os.path.join(self.basedir,'%4i' % time.year,'%02i' % time.month)
            if self.subdir is not None:
                folder=os.path.join(folder,self.subdir)
            if 'subpath' in self.info:
                folder=os.path.join(folder,self.info['subpath'])
            fidx=fidx+1
            try:
                #print 'listing path ' , folder
                flist=os.listdir(folder)
                flist.sort()
                #print flist
            except OSError:
                flist=[]
            idxoffile=-1
            fileindex=0
            while idxoffile<fidx and fileindex<len(flist):
                f=flist[fileindex]
                if ('prefix' not in self.info or f.startswith(self.info['prefix'])) and ('suffix' not in self.info or f.endswith(self.info['suffix'])):
                    idxoffile=idxoffile+1
                fileindex=fileindex+1

            if idxoffile==fidx:
                fnr={'path':os.path.join(folder,f),'filename':f,'start':self.parseTimeFromFile(f,time,self.info),'width':timedelta(days=1),'pathtime':time}
                if fnr['start'] is None:
                    raise RuntimeError('Bad filename date in '+fnr['path'])
                if not pfnr is None and (self.completeList or fnr['start']>pfnr['start']):
                    pfnr['width']=fnr['start']-pfnr['start']
                    if (self.end!=None and pfnr['start']>=self.end) or fnr['start']<=self.start:
                        pass
                    else:
                        yield pfnr
                if pfnr is None or pfnr['start']<=fnr['start'] or self.completeList:
                    pfnr=fnr
                else:
                    print 'WARNING FILES OUT OF ORDER ',pfnr['path'],fnr['path']
            else:
                fidx=-1
                time+=calendarDays(time)
        if pfnr:
            if self.end:
                pfnr['width']=self.end-pfnr['start']
            else:
                pfnr['width']=datetime.utcnow()-pfnr['start']
            if (self.end!=None and pfnr['start']>=self.end) or (pfnr['start']+pfnr['width'])<=self.start:
                pass
            else:
                yield pfnr


class VectorTableLibrarian(dplkit.role.librarian.aLibrarian):
    """ Librarian Initialzation for Calibration Tables
  
            :param siteid: source site id for the data source. typically an hsrl instrument or base directory
            :param datatype: file type to list
    """
    def __init__(self, instrumentname,datatype,basedir,**kwargs):
        super(VectorTableLibrarian,self).__init__()
        self.basedir=basedir
        self.datatype=datatype
        self.instrumentname=instrumentname
        self.kwargs=kwargs.copy()

    def search(self,start,end,*args,**kwargs):
        """ Librarian Generator function
        extra parameters given here will be passed to the returned narrator's init
        """
        myargs=self.kwargs.copy()
        myargs.update(kwargs)
        return VectorTableNarrator(self.basedir,self.instrumentname,self.datatype,start,end,*args,**myargs)
