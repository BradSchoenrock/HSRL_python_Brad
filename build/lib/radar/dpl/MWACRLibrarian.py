
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
import RadarFilters

@dplkit.role.decorator.exposes_attrs_in_chain(['radarType','radarLambda','radarNativeTimeResolution'])
class MWACRNarrator(ARMNarrator.ARMNarrator):
    """ Narrator Initialzation for MWACR.  Shouldn't be created, should only be used by the MWACRLibrarian

            :param host: source librarian object
            :param basedir: base directory to find MWACR data
            :param dataprefix: filename prefixes for MWACR data
            :param start: start time
            :param end: end time
            :param zoo: zookeeper, if narrator should return data, not file info
            :param kwargs: addtional parameters are passed to thezookeeper on read, if provided

         exposed attributes:
            - radarType : (string) identifying name
            - radarLambda : (float) lambda of laser
            - radarNativeTimeResolution : (timedelta) native time resolution (single typical resolution)
    """

    @property
    def radarType(self):
        return self.instrumentname.upper()
    @property
    def radarLambda(self):
        return 3.154e-3

    @property
    def radarNativeTimeResolution(self):
        if self.timeres==None:
            if False:#self.zoo!=None:
                for fr in self:
                    if not hasattr(fr,'times') or fr.times.size<2:
                        continue
                    self.timeres=timedelta(seconds=(fr.times[-1]-fr.times[0]).total_seconds()/(fr.times.size-1))
                    break
            if self.timeres==None:
                self.timeres=timedelta(seconds=.5)
        return self.timeres

    def preYield(self,x,attrs,found):
        RadarFilters.addCommonAttributes(x,attrs)
        return True

    def __init__(self,host,basedir,instname,dataprefix,start,end):
        super(MWACRNarrator,self).__init__(basedir,instname,dataprefix,'radar/'+instname,start,end)
        self.host=host
        self.instrumentname=instname
        self.timeres=None

class MWACRLibrarian(dplkit.role.librarian.aLibrarian):
    """ Librarian Initialzation for MWACR
  
            :param siteid: source site id for the data source. typically an hsrl instrument or base directory
            :param dataprefix: filename prefixes for MWACR data
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
        ret=MWACRNarrator(self,self.basedir,self.instrumentname[self.instrumentname.find('mwacr'):],self.dataprefix,start,end)
        if self.zoo is not None:
            ret=ARMNarrator.ARMFileNarrator(ret,self.zoo,ret.preYield,*args,**kwargs)
        return ret

def main():
    from lg_dpl_toolbox.dpl.NetCDFZookeeper import GenericTemplateRemapNetCDFZookeeper 
    zoo=None#GenericTemplateRemapNetCDFZookeeper('mwacr')
    lib=MWACRLibrarian('/data/mf2hsrldata','magmwacr','magmwacrM1.a1.',zoo=zoo)
    zoo=GenericTemplateRemapNetCDFZookeeper('mwacr')
    m=lib(start=datetime(2013,6,21,20,0,0),end=datetime(2013,6,25,0,0,0))

    for f in m:
        #print 'from librarian:',f
        if zoo:
            res=zoo.open(zoo(uri=f))
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
