import dplkit.role.narrator
from datetime import datetime
import numpy as np
from parsing import flagParser

import dplkit.role.librarian

class dpl_hsrl_qa(dplkit.role.librarian.aLibrarian):
    """ HSRL_QA wrapper framestream

    :param instrument: instrument id string (eg. 'bagohsrl'). also supports naming the co-located HSRL to get a default radar source

    """

    def __init__(self, instrument,*args, **kwargs):
        super(self.__class__,self).__init__(None)
        import dpl_filters as qafilt
        import dpl_zookeeper as qazoo
        import dpl_librarian as qalib
        self.qafilt=qafilt
        self.instrument=instrument
        self.zoo=qazoo.QualityAssuranceZookeeper()
        self.lib=qalib.QualityAssuranceLibrarian(instrument=instrument)

    def __repr__(self):
        return 'DPL HSRL_QA Librarian (instrument="%s")' % (self.instrument)

    def search(self, start_time_datetime, end_time_datetime,hsrl_host=None,timealtsource=None,hostsource=None,hostsource_newframe=None,altname=None,constantAltitude=None,*args, **kwargs):
        """
        :param start_time_datetime: start time 
        :type start_time_datetime: datetime.datetime
        :param end_time_datetime: end time
        :type end_time_datetime: datetime.datetime
        """
        if len(args):
            print 'Unused dpl_hsrl_qa.search args = ',args
        if len(kwargs):
            print "Unused dpl_hsrl_qa.search kwargs = ",kwargs

        assert(timealtsource!=None or hostsource!=None)
        hsrl_host = hsrl_host or timealtsource or hostsource

        filecontent=self.zoo(self.lib(start_time=start_time_datetime,end_time=end_time_datetime))
        qcnarr=QCFlagNarrator(filecontent,altitude_axis=hsrl_host.altitudeAxis,binwidth=hsrl_host.hsrl_constants['binwidth']*1.5e8)#narrator streams the read flags as a single timestep per entry

        if timealtsource!=None:
            qcnarr=self.qafilt.QAFlagClonedSyncFilter(qasource=qcnarr,timealtsource=timealtsource,altname=altname,constantAltitude=constantAltitude)#takes the output of the narrator, and echos entries for each timestep (and potentially project range entries to altitude)
        else:
            qcnarr=self.qafilt.QAFlagClonedAttachFilter(qasource=qcnarr,hostsource=hostsource,hostsource_newframe=hostsource_newframe,altname=altname,constantAltitude=constantAltitude)#takes the output of the narrator, and echos entries for each timestep (and potentially project range entries to altitude)            

        return qcnarr

 
#this will contain narrators that use the zookeeper to create streams of status strings. those status strings are transformed by other narrators to
# make bit-arrays, matrixes, and such.  some arrays will be in altitude, others will be in range. at the matrix generation step, a platform altitude
# would be neccessary if range is used, for there the range will be converted to altitudes on every interval, and merged into the existing altitude array
# in accordance to the time/altitude stream given to that object (RuntimeError if range parameters are given and altitude isn't)

class QCFlagNarrator(dplkit.role.narrator.aNarrator):
    def __init__(self,content,altitude_axis=None,binwidth=None,start_time=None,end_time=None,translator=None):
        super(QCFlagNarrator,self).__init__(None)
        self.listArray=content['flags']
        self.starttime=start_time
        self.endtime=end_time
        self.flagparser=translator or flagParser(altitude_axis,binwidth=binwidth)
        self.provides={}
        self.provides['time']=dict(shortname='time',type=datetime)
        self.provides['flags']=dict(shortname='flags',type=np.ndarray)
        self.provides['range_flags']=dict(shortname='range_flags',type=np.ndarray)

    @property
    def flagmanager(self):
        return self.flagparser

    def makeFrame(self,entry):
        ret=dict(time=entry['header']['date'],flags=None,range_flags=None)
        ret['flags'],ret['range_flags']=self.flagparser.parseFlags(entry['content'])
        return ret

    def read(self):
        prior=None
        priorprior=None
        for x in self.listArray:
            priorprior=prior
            prior=x
            if self.starttime!=None and x['header']['date']<self.starttime:
                continue
            if priorprior!=None:
                yield self.makeFrame(priorprior)
            if self.endtime!=None and x['header']['date']>=self.endtime:
                break
        if prior!=None and (self.endtime==None or prior['header']['date']<self.endtime):
            yield self.makeFrame(prior)


class QCLogNarrator(dplkit.role.narrator.aNarrator):
    def __init__(self,content,starttime=None,endtime=None):
        super(QCLogNarrator,self).__init__(None)
        self.listArray=content['logs']
        self.starttime=starttime
        self.endtime=endtime

    def read(self):
        prior=None
        priorprior=None
        for x in self.listArray:
            priorprior=prior
            prior=x
            if self.starttime!=None and x['header']['date']<self.starttime:
                continue
            if priorprior!=None:
                yield priorprior
            if self.endtime!=None and x['header']['date']>=self.endtime:
                break
        if prior!=None and (self.endtime==None or prior['header']['date']<self.endtime):
            yield prior


def main():
    import sys
    l=dpl_hsrl_qa(sys.argv[1])
    s=datetime.strptime(sys.argv[2],"%Y%m%d")
    e=datetime.strptime(sys.argv[3],"%Y%m%d")
    for x in l(s,e):
        print x

if __name__ == '__main__':
    main()
