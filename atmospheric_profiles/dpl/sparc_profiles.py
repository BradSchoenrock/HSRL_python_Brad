import dplkit.role.librarian
import dplkit.role.narrator
import dplkit.role.filter
from datetime import datetime,timedelta
import numpy
import dplkit.role.decorator
import lg_base.core.array_utils as hau
import os,calendar
import lg_dpl_toolbox.dpl.ARMNarrator as ARMNarrator
import re

class SPARCTimeParser(object):
    def __init__(self):
        self.datematch=re.compile('_[0-9]{8}_[0-9]{6}\.')

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
        return datetime.strptime(tmp,'_%Y.%m%d_%H%M%S.')

class SPARCSondeNarrator(ARMNarrator.ARMNarrator):
    """ Narrator Initialzation for MWACR.  Shouldn't be created, should only be used by the MWACRLibrarian

            :param host: source librarian object
            :param basedir: base directory to find MWACR data
            :param dataprefix: filename prefixes for MWACR data
            :param start: start time
            :param end: end time
            :param zoo: zookeeper, if narrator should return data, not file info
            :param kwargs: addtional parameters are passed to thezookeeper on read, if provided

    """

    def __init__(self,host,basedir,start,end):
        super(SPARCSondeNarrator,self).__init__(basedir,"SPARCSonde",'sparc_radiosondes_','sondes',start,end,timeparse=SPARCTimeParser())
        self.host=host
        #self.instrumentname="SparcSonde"
        #self.timeres=None
        import atmospheric_profiles.soundings.sounding_utilities as su
        self.su=su

    def preYield(self,x,attrs,found):
        if not hasattr(x,'altitudes') or not hasattr(x,'times') or x.times.size==0:
            return False
        altmask=x.altitudes<90000 
        for n,v in vars(x).items():
            #print n,v
            if n.startswith('_'):
                continue
            v=v[altmask]
            if n=='altitudes':
                setattr(x,'_altitudevarname',n)
                setattr(x,n,hau.Z_Array(v).copy())
                continue
            if n in ('latitude','longitude','times'):
                if v.size>0:
                    setattr(x,n,v[0])
                continue
            setattr(x,n,hau.Z_Array(v.astype('double')))#.reshape([1]+list(v.shape))))
        x.temps+=273.15
        if not hasattr(x,'dew_points') and hasattr(x,'relative_humidity'):
            setattr(x,'dew_points',hau.Z_Array(self.su.cal_dew_point(x.relative_humidity,x.temps)))
        else:
            x.dew_points+=273.15
        if not hasattr(x,'frost_points') and hasattr(x,'dew_points'):
            setattr(x,'frost_points',hau.Z_Array(self.su.cal_frost_point(x.dew_points)))
        else:
            x.frost_points+=273.15

        x.sounding_type='sparc'
        x.sounding_id=attrs['location_code']
        x.station_id=attrs['location_code']
        x.top=max(x.altitudes)
        x.bot=min(x.altitudes)

        #print vars(x)
        #if x.times.size<=0:
        #    return False
        x.sample_latitude=x.latitude
        x.sample_longitude=x.longitude
        x.sample_time=x.times
        return True


class SPARCSondeLibrarian(dplkit.role.librarian.aLibrarian):
    """ Librarian Initialzation for MWACR
  
            :param siteid: source site id for the data source. typically an hsrl instrument or base directory
            :param dataprefix: filename prefixes for MWACR data
            :param zoo: zookeeper, if narrator should return data, not file info
    """
    def __init__(self, siteid,requested_altitudes=None,zoo=None):
        super(self.__class__,self).__init__()
        self.basedir=ARMNarrator.getBasedir(siteid)
        self.requested_altitudes=requested_altitudes
        self.zoo=zoo

    def search(self,start,end,*args,**kwargs):
        """ Librarian Generator function
        extra parameters given here will be passed to the returned narrator's init
        """
        ret=SPARCSondeNarrator(self,self.basedir,start,end)
        zoo=self.zoo
        if zoo is None and not kwargs.pop('filenames',False):
            from lg_dpl_toolbox.dpl.NetCDFZookeeper import GenericTemplateRemapNetCDFZookeeper 
            zoo=GenericTemplateRemapNetCDFZookeeper('sparcradiosondes')
        if zoo is not None:
            ret=ARMNarrator.ARMFileNarrator(ret,zoo,preYield=ret.preYield,*args,**kwargs)
        if self.requested_altitudes is not None:
            import atmospheric_profiles.dpl.dpl_temperature_profiles as dtp
            ret=dtp.dpl_radiosonderesample(ret,'altitudes',self.requested_altitudes,{'temps':[50,500],'pressures':[1,1000],
                'dew_points':[150,500],'frost_points':[50,500],'altitudes':None})
        return ret;


def main():
    sp=SPARCSondeLibrarian('bagohsrl')
    for x in sp(datetime(2015,6,25,1,0,0),datetime(2015,6,25,23,0,0)):
        print x
        #print vars(x)

if __name__ == '__main__':
    main()
