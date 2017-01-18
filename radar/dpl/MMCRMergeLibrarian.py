
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
class MMCRMergeNarrator(ARMNarrator.ARMNarrator):
    """ Narrator Initialzation for MMCR. Shouldn't be created, should only be used by the MMCRMergeLibrarian
  
            :param host: source librarian object
            :param basedir: base directory to find MMCR data
            :param dataprefix: filename prefixes for MMCR data
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
        return 'MMCR'

    @property
    def radarLambda(self):
        return 8.6e-3

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
                self.timeres=timedelta(seconds=10)
        return self.timeres

    def preYield(self,x,attrs,found):
        RadarFilters.addCommonAttributes(x,attrs)
        return True

    def __init__(self,host,basedir,dataprefix,start,end):
        super(MMCRMergeNarrator,self).__init__(basedir,'MMCR',dataprefix,'radar/merge',start,end)
        self.basedir=basedir
        self.timeres=None


class MMCRMergeLibrarian(dplkit.role.librarian.aLibrarian):
    """ Librarian Initialzation for MMCR
  
            :param siteid: source site id for the data source. typically an hsrl instrument or base directory
            :param dataprefix: filename prefixes for MMCR data
            :param zoo: zookeeper, if narrator should return data, not file info
    """
    def __init__(self, siteid,dataprefix,zoo=None):
        super(self.__class__,self).__init__()
        self.basedir=ARMNarrator.getBasedir(siteid)
        self.dataprefix=dataprefix
        #self.instrumentname=instrumentname
        self.zoo=zoo

    def search(self,start,end,*args,**kwargs):
        """ Librarian Generator function
        extra parameters given here will be passed to the returned narrator's init
        """
        ret=MMCRMergeNarrator(self,self.basedir,self.dataprefix,start,end)
        if self.zoo is not None:
            ret=ARMNarrator.ARMFileNarrator(ret,self.zoo,ret.preYield,*args,**kwargs)
        return ret

if __name__=='__main__':
    from lg_dpl_toolbox.dpl.NetCDFZookeeper import GenericTemplateRemapNetCDFZookeeper 
    zoo=GenericTemplateRemapNetCDFZookeeper('eurmmcrmerge')
    lib=MMCRMergeLibrarian('/data/ahsrldata','eurmmcrmerge.C1.c1.',zoo=zoo)
    zoo=None
    m=lib(start=datetime(2006,12,24,0,0,0),end=datetime(2006,12,25,0,0,0))

    for f in m:
        #print 'from librarian:',f
        if zoo:
            res=zoo(uri=f)
            print 'uri from zoo:',res
        else:
            print f
        #print 'content=',zoo.open(res)
