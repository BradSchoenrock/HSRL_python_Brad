#!/usr/bin/env python
from datetime import datetime,timedelta
#import hsrl.utils.json_config as jc
import sys,os
#import argparse
# adapted by Joe VanAndel from hsrl/dpl/dpl_live.py


import lg_dpl_toolbox.filters.substruct as frame_substruct
from hsrl.dpl.dpl_hsrl import dpl_hsrl
#import hsrl.dpl.dpl_artists as artists
import lg_dpl_toolbox.dpl.dpl_artists as tools_artists

#import wingdbstub

class netcdf_namer:
    """
    Create filenames with a specified path, prefix, suffix, and start-time
    :param path: directory name
    :param prefix
    :param suffix
    """
    def __init__(self,path,prefix,suffix):
        self.path=path
        self.prefix=prefix
        self.suffix=suffix

    def __call__(self,starttime,endtime,*args,**kwargs):
        return os.path.join(self.path,self.prefix+starttime.strftime('_%Y%m%dT%H%M_')
                            + self.suffix)


class SourceStream:
    """
        :param instruumnt ('gvhsrl', 'bagohsrl', etc)
        :param starttime: starttime 
        :type starttime: datetime
        :param endtime: endtime 
        :type endtime: datetime
        :param process_control: process control file or object for processing parameters
        :param minalt_km: minimum altitude in km. default is 0
        :param maxalt_km: maximum altitude in km. default is 15

        :
        """    
    def __init__(self, instrument, startTime, endTime, 
                 window=5*60 ,
                 process_control = None, minalt_km=0, 
                 maxalt_km = 15):
 
        self.instrument = instrument
        self.startTime = startTime
        self.endTime = endTime
        self.currentStart = startTime
        self.window = timedelta(seconds=window)
        self.process_control = process_control
        self.minalt_km = minalt_km
        self.maxalt_km = maxalt_km
        frame_substruct.SubstructBrancher.multiprocessable=False 
        self.dplobj=dpl_hsrl(instrument=self.instrument,
                        process_control=self.process_control)
        
    def __iter__(self):
        return self
    
    def next(self):
        while 1:
            if self.currentStart >= self.endTime:
                print "Completed processing of data from: ", self.startTime, \
                " to ", self.endTime
                raise StopIteration()
            currentStart = self.currentStart
            self.currentStart += self.window

            try:
                dplgen=self.dplobj(start_time_datetime=currentStart,
                                         end_time_datetime = self.currentStart,
                          reverse_padding=timedelta(seconds=60),
                          min_alt_m=self.minalt_km*1000.0,
                          max_alt_m=self.maxalt_km*1000.0,
                          with_profiles=False)
                break
            except KeyError, e:
                print 'Error processing ', currentStart, " thru ",self.currentStart
                print 'Trying next interval'
                continue

        print "SourceStream.next() start = ", currentStart, " end =",self.currentStart
        return dplgen, currentStart, self.currentStart


def livestream(instrument,template, starttime, endtime,
    process_control=None,maxalt_km=15,minalt_km=0, 
               window=5*60, outputDir='/tmp', filename_suffix='data_cset'):
    """
    Quick and simple livestream of HSRL data using DPL constructs.

    :param instrument: instrument name ('gvhsrl','bagohsrl', etc)
    :param starttime: starttime 
    :type starttime: datetime
    :param endtime: endtime 
    :type endtime: datetime
    :param process_control: process control file or object for processing parameters
    :param minalt_km: minimum altitude in km. default is 0
    :param maxalt_km: maximum altitude in km. default is 15
    """

    # create a generator that produces tuples of
    # streams, startTime, endTime
    srcStream = SourceStream(instrument,
                             starttime, endtime, window,
                             process_control,
                             minalt_km , maxalt_km)
    #dplgen = srcStream.dplgen()

    
    artist=tools_artists.dpl_multi_netcdf_artist( \
        sourcestreamstream = srcStream,  template=template,
        storagepath=outputDir, 
        filename_maker=
            netcdf_namer(outputDir, instrument, filename_suffix+".nc"))
         
        
    return artist

def main():
    #parser = argparse.ArgumentParser(description='create CfRadial files')
    addpath=os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]),os.path.pardir,os.path.pardir))
    sys.path.insert(0,addpath)
    idx=2
    parms=dict(template =  '/usr/local/hsrl/config/5min_hsrl_cfradial.cdl',
                process_control = 'process_control.json')
    while idx<len(sys.argv) and sys.argv[idx][0]=='-':
        sw=sys.argv[idx]
        parm =sys.argv[idx+1]
        idx+=2
        if sw in ('-a','--altmax','--maxalt'):
            parms['maxalt_km']=float(parm)
        elif sw in ('-h','--duration','--hours'):
            parms['hours']=float(parm)
        elif sw in ('-i','--delay','--ignore'):
            parms['now_delay']=timedelta(seconds=float(parm))
        elif sw in ('-s','--start'):
            try:
                # did user specify a relative start time in hours 
                parms['starttime']=datetime.utcnow()-timedelta(hours=float(parm))
            except ValueError:
                parms['starttime']=datetime.strptime(parm,'%Y%m%dT%H%M')
        elif sw in ('-e','--end'):
                parms['endtime']=datetime.strptime(parm,'%Y%m%dT%H%M')

        elif sw in ('-p','--processcontrol','--parameters','--process'):
            parms['process_control']=parm
        elif sw in ('-o','--output','--outputdir'):
            parms['outputDir']=parm
        else:
            RuntimeError('Unknown parameter %s : %s' % (sw,parm))

    if len(sys.argv)!=idx:
        print 'usage: %s instrument [-a maxalt_km] [-d display_json] [-h hours_duration]'
        print '\t -s --start :\tstart time as YYYYMMDDThhmm'
        print '\t -s --end :\tend time as YYYYMMDDThhmm'
        print '\t -a --altmax :\tspecify maximum altitude in km (default 15)'
        print '\t -i --ignore :\tignore x seconds prior to now (default 30)'
        print '\t -p --process :\tprocess control json filename'
        return 0
    instrument=sys.argv[1]
    artist=livestream(instrument,**parms)
    f=1
    for x in artist:
        print 'process frame %d' % f
        f+=1


if __name__ == '__main__':
    main()
