
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
import lg_dpl_toolbox.dpl.ARMNarrator as ARMNarrator

@dplkit.role.decorator.exposes_attrs_in_chain(['vdisType','vdisNativeTimeResolution'])
class VDisNarrator(ARMNarrator.ARMNarrator):
    """ Narrator Initialzation for PARS.  Shouldn't be created, should only be used by the PARSLibrarian

            :param host: source librarian object
            :param basedir: base directory to find PARS data
            :param dataprefix: filename prefixes for PARS data
            :param start: start time
            :param end: end time
            :param zoo: zookeeper, if narrator should return data, not file info
            :param kwargs: addtional parameters are passed to thezookeeper on read, if provided

         exposed attributes:
            - parsType : (string) identifying name
            - parsNativeTimeResolution : (timedelta) native time resolution (single typical resolution)
    """

    @property
    def vdisType(self):
        return self.instrumentname.upper()

    @property
    def vdisNativeTimeResolution(self):
        if self.timeres==None:
            if False:#self.zoo!=None:
                for fr in self:
                    if not hasattr(fr,'times') or fr.times.size<2:
                        continue
                    self.timeres=timedelta(seconds=(fr.times[-1]-fr.times[0]).total_seconds()/(fr.times.size-1))
                    break
            if self.timeres==None:
                self.timeres=timedelta(seconds=60)
        return self.timeres

    def __init__(self,host,basedir,instname,dataprefix,start,end):
        super(VDisNarrator,self).__init__(basedir,instname,dataprefix,'vdis',start,end)
        self.basedir=basedir
        self.instrumentname=instname
        self.timeres=None


class VDisLibrarian(dplkit.role.librarian.aLibrarian):
    """ Librarian Initialzation for PARS
  
            :param siteid: source site id for the data source. typically an hsrl instrument or base directory
            :param dataprefix: filename prefixes for PARS data
            :param zoo: zookeeper, if narrator should return data, not file info
    """
    def __init__(self, siteid,instrumentname,dataprefix,zoo=None):
        super(self.__class__,self).__init__()
        self.basedir=ARMNarrator.getBasedir(siteid)
        self.dataprefix=dataprefix
        self.instrumentname=instrumentname
        self.zoo=zoo

    def search(self,start,end,*args,**kwargs):
        """ Librarian Generator function
        extra parameters given here will be passed to the returned narrator's init
        """
        ret=VDisNarrator(self,self.basedir,self.instrumentname[self.instrumentname.find('vdis'):],self.dataprefix,start,end)
        if self.zoo is not None:
            ret=ARMNarrator.ARMFileNarrator(ret,self.zoo,*args,**kwargs)
        return ret

def main():
    from lg_dpl_toolbox.dpl.NetCDFZookeeper import GenericTemplateRemapNetCDFZookeeper 
    zoo=GenericTemplateRemapNetCDFZookeeper('vdis')
    lib=VDisLibrarian('/data/mf2hsrldata','mf2vdis','mf2vdisM1',zoo=zoo)
    zoo=None
    m=lib(start=datetime(2014,6,24,20,0,0),end=datetime(2014,7,24,22,0,0))

    for f in m:
        #print 'from librarian:',f
        if zoo:
            res=zoo(uri=f)
            print 'uri from zoo:',res
        else:
            print f
            print vars(f)
        #print 'content=',zoo.open(res)

if __name__=='__main__':
    if False:
        import cProfile
        cProfile.run('main()',sort='tottime')
    else:
        main()
