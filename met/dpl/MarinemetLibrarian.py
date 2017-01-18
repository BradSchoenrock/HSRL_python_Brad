
import dplkit.role.librarian
import dplkit.role.narrator
import dplkit.role.filter
from datetime import datetime,timedelta
import os
import re
import time
import numpy
import lg_dpl_toolbox.filters.substruct as substruct
import lg_base.core.array_utils as hau
import dplkit.role.decorator
import lg_dpl_toolbox.filters.time_frame as time_frame

@dplkit.role.decorator.exposes_attrs_in_chain(['metNativeTimeResolution'])
@dplkit.role.decorator.autoprovide(frameclass=hau.Time_Z_Group)
class MarinemetNarrator(dplkit.role.narrator.aNarrator):
    """ Narrator Initialzation for Marinemet.  Shouldn't be created, should only be used by the MarinemetLibrarian

            :param host: source librarian object
            :param basedir: base directory to find Marinemet data
            :param dataprefix: filename prefixes for Marinemet data
            :param start: start time
            :param end: end time
            :param zoo: zookeeper, if narrator should return data, not file info
            :param kwargs: addtional parameters are passed to thezookeeper on read, if provided

         exposed attributes:
            - metNativeTimeResolution : (timedelta) native time resolution (single typical resolution)
    """
    @property
    def metNativeTimeResolution(self):
        if self.timeres==None:
            if self.zoo!=None:
                for fr in self:
                    if not hasattr(fr,'times') or fr.times.size<2:
                        continue
                    self.timeres=timedelta(seconds=(fr.times[-1]-fr.times[0]).total_seconds()/(fr.times.size-1))
                    break
            if self.timeres==None:
                self.timeres=timedelta(seconds=60)
        return self.timeres

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
        return datetime.strptime(tmp,'.%Y.%m%d.%H%M%S.')


    def __init__(self,host,basedir,dataprefix,start,end,zoo=None,*args,**kwargs):
        self.host=host
        self.basedir=basedir
        if isinstance(dataprefix,basestring):
            self.dataprefix=[dataprefix]
        else:
            self.dataprefix=dataprefix
        self.start=start
        self.end=end
        self.zoo=zoo
        self.allkwargs=kwargs.copy()
        self.timeres=None
        self.starttime=datetime(start.year,start.month,start.day,0,0,0)-timedelta(days=1)

    def read(self,*args,**kwargs):
        """ Narrator Generator function
        takes no parameters. initialization determines Narrator behavior
        parameters given here will be passed to the zookeeper, if narrator was given one
        """
        time=self.starttime
        fidx=-1
        pfnr=None
        allkwargs=self.allkwargs.copy()
        allkwargs.update(kwargs)
        possibleprefixes=('met','magmarinemet')
        while time<=self.end and time<=datetime.utcnow():
            for p in possibleprefixes:
                folder=os.path.join(self.basedir,'%4i' % time.year,'%02i' % time.month,'%02i' % time.day,p)
                if os.path.exists(folder):
                    break
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
                for pref in self.dataprefix:
                    if f.startswith(pref):
                        idxoffile=idxoffile+1
                fileindex=fileindex+1

            if idxoffile==fidx:
                fnr={'path':os.path.join(folder,f),'filename':f,'start':self.parseTimeFromFile(f),'width':timedelta(days=1)}
                if pfnr:
                    pfnr['width']=fnr['start']-pfnr['start']
                    if (self.end!=None and pfnr['start']>=self.end) or fnr['start']<=self.start:
                        pass
                    elif self.zoo==None:
                        yield pfnr
                    else:
                        print pfnr
                        print self.zoo(pfnr)
                        x=self.zoo.open(self.zoo(pfnr),firsttime=self.start,lasttime=self.end,*args, **allkwargs)
                        if x!=None:
                            yield x
                pfnr=fnr
            else:
                fidx=-1
                time+=timedelta(days=1)
        if pfnr:
            if self.end:
                pfnr['width']=self.end-pfnr['start']
            else:
                pfnr['width']=datetime.utcnow()-pfnr['start']
            if (self.end!=None and pfnr['start']>=self.end) or (pfnr['start']+pfnr['width'])<=self.start:
                pass
            elif self.zoo==None:
                yield pfnr
            else:
                x=self.zoo.open(self.zoo(pfnr),firsttime=self.start,lasttime=self.end)
                if x!=None:
                    yield x



class MarinemetLibrarian(dplkit.role.librarian.aLibrarian):
    """ Librarian Initialzation for Marinemet
  
            :param siteid: source site id for the data source. typically an hsrl instrument or base directory
            :param dataprefix: filename prefixes for Marinemet data
            :param zoo: zookeeper, if narrator should return data, not file info
    """
    def __init__(self, siteid,dataprefix,zoo=None):
        super(self.__class__,self).__init__()
        import lg_dpl_toolbox.core.archival as hru
        if isinstance(siteid,basestring):
            if os.access(siteid,os.F_OK):
                self.basedir=siteid
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
        self.basedir=tmp['Path']
        self.dataprefix=dataprefix
        self.zoo=zoo
        self.datematch=re.compile('\.[0-9]{8}\.[0-9]{6}\.')

    def search(self,start,end,*args,**kwargs):
        """ Librarian Generator function
        extra parameters given here will be passed to the returned narrator's init
        """
        return MarinemetNarrator(self,self.basedir,self.dataprefix,start,end,zoo=self.zoo,*args,**kwargs)

def main():
    from lg_dpl_toolbox.dpl.NetCDFZookeeper import GenericTemplateRemapNetCDFZookeeper 
    zoo=None#GenericTemplateRemapNetCDFZookeeper('kazr')
    lib=KAZRLibrarian('/data/mf2hsrldata','magkazrgeM1.a1.',zoo=zoo)
    zoo=None
    m=lib(start=datetime(2013,7,24,20,0,0),end=datetime(2013,7,25,0,0,0))

    for f in m:
        #print 'from librarian:',f
        if zoo:
            res=zoo(uri=f)
            print 'uri from zoo:',res
        else:
            print f
        #print 'content=',zoo.open(res)

if __name__=='__main__':
    if False:
        import cProfile
        cProfile.run('main()',sort='tottime')
    else:
        main()
