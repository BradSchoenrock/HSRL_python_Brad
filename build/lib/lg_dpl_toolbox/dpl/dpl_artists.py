import dplkit.role.artist
import dplkit.role.filter
import dplkit.role.decorator
import copy
from datetime import datetime,timedelta
import calendar
import os,types
import traceback
import sys

class dpl_window_caching_filter(dplkit.role.filter.aFilter):
    """ A caching filter that from the given frame stream, uses Time_Z_Group append() and trimTimeInterval() to maintain a specific amount of time in the frame, and yields that always

    :param framestream: input iterable framestream
    :param thewidth: time duration of the window
    :type thewindow: timedelta
    """
    def __init__(self,framestream,thewidth,includeIncompleteFrames=False):
        super(dpl_window_caching_filter,self).__init__(framestream)
        self.framestream=framestream
        self.width=thewidth
        self.includeIncompleteFrames=includeIncompleteFrames
 
    def process(self):
        frame=None
        for subframe in self.framestream:
            assert(subframe is not None)
            if frame==None:
                frame=copy.deepcopy(subframe)
            elif hasattr(frame,'append'):
                frame.append(subframe)
            else:
                assert('Dont know how to append')
            if frame.trimTimeInterval(width=self.width) or self.includeIncompleteFrames:
                yield frame

class default_multi_netcdf_namer:
    def __init__(self,path,prefix,suffix):
        self.path=path
        self.prefix=prefix
        self.suffix=suffix

    def __call__(self,starttime,endtime,*args,**kwargs):
        return os.path.join(self.path,self.prefix+starttime.strftime('_%Y%m%dT%H%M')+endtime.strftime('_%Y%m%dT%H%M')+self.suffix)

def stepsize_iterator(starttime,floorstep,staticstep=None,mode=None):
    if staticstep!=None:
        yield (floorstep+staticstep)-starttime
        while True:
            yield staticstep
    t=floorstep
    if mode=='month':
        info=calendar.monthrange(starttime.year,starttime.month)
        dt=timedelta(days=info[1])
        yield (floorstep+dt)-starttime
        while True:
            t=t+dt
            info=calendar.monthrange(t.year,t.month)
            dt=timedelta(days=info[1])
            yield dt
    

def multi_netcdf_filewindow(namestart,nameend,starttime,endtime,mode):
    stime=starttime
    etime=starttime
    if mode=='single':
        yield {namestart:starttime,nameend:endtime,'mode':mode},starttime,endtime
        return
    elif mode=='30minute':
        stepper=stepsize_iterator(starttime,starttime.replace(minute=starttime.minute-(starttime.minute%30),second=0,microsecond=0),timedelta(seconds=30*60))
    elif mode=='hour':
        stepper=stepsize_iterator(starttime,starttime.replace(minute=0,second=0,microsecond=0),timedelta(seconds=60*60))
    elif mode=='day':
        stepper=stepsize_iterator(starttime,starttime.replace(hour=0,minute=0,second=0,microsecond=0),timedelta(days=1))
    elif mode=='week':
        #info=calendar.monthrange(starttime.year,starttime.month)
        stepper=stepsize_iterator(starttime,starttime.replace(hour=0,minute=0,second=0,microsecond=0),timedelta(days=7))
    elif mode=='month':
        stepper=stepsize_iterator(starttime,starttime.replace(day=1,hour=0,minute=0,second=0,microsecond=0),mode="month")
    else:
        raise RuntimeError
    for dt in stepper:
        etime = stime+dt
        if etime>endtime:
            etime=endtime
        yield {namestart:stime,nameend:etime,'mode':mode},stime,etime
        if etime==endtime:
            return
        stime=etime


class dpl_multi_netcdf_artist(object):
    def __init__(self,sourcestreamstream,storagepath=None,filename_maker=None,yieldArtists=False,callFirst=None,*args,**kwargs):
        super(dpl_multi_netcdf_artist,self).__init__()
        self.sourcestreamstream=sourcestreamstream
        self.callFirst=callFirst
        if filename_maker:
            self.namer=filename_maker
        else:
            self.namer=default_multi_netcdf_namer(storagepath,"multi",".nc")
        self.args=args
        self.kwargs=kwargs
        self.yieldArtists=yieldArtists

    def artist_iter(self):
        for stream,stime,etime in self.sourcestreamstream:
                print self.namer(stime,etime)
                args=copy.copy(self.args)
                kwargs=copy.copy(self.kwargs)
                kwargs.update(dict(outputfilename=self.namer(stime,etime)))
                if self.callFirst!=None:
                    self.callFirst(stream,args,kwargs)
                try:
                    fileartist=dpl_netcdf_artist(stream,*args,**kwargs)
                except Exception, e:
                    print "artist_iter: dpl_netcdf_artist()  threw an exception", e
                    print traceback.format_exc()
            
                yield fileartist
  
    def __iter__(self):
        for fileartist in self.artist_iter():
            if self.yieldArtists:
                yield fileartist
            else:
                try:
                    for f in fileartist:
                        yield f
                except:
                    print 'error occurred'
                    print traceback.format_exc()
  
def datasetParametersFor(template,cfradial=None):
    if template.endswith(".json"):
        import json
        d=json.load(file(template))
        fmt=d.pop('format','NETCDF4')
        if fmt.lower()=='cfradial':
            return "NETCDF3_64BIT",True
        return fmt,False
    basetemplate=os.path.basename(template)
    if cfradial==None:
        cfradial=('cfradial' in basetemplate)
    if '3' in basetemplate and not cfradial:
        return 'NETCDF3_64BIT',cfradial
    return 'NETCDF4',cfradial
 
class dpl_netcdf_artist(dplkit.role.artist.aArtist):
    """ Artist for outputing a NetCDF from a named CDL template

    :param framestream: DPL iterable framestream
    :param template: CDL template filename
    :param outputfilename: output NetCDF filename
    :param format: NetCDF format type. Default is "NETCDF4" unless the template contians a "3", then "NETCDF3_64BIT". this will override the discovered default
    :param usecfradial: if true, will use cfradial output module. if false, will use uw hsrl module. if None, will use cfradial if 'cfradial' appears in the template name. this will override the discovered default
    :param selected_bindings: list of source object paths to include in the output file. if not specified, will include all available that are in the template.
    :param forModule: if the constructing entity prefer this search for templates in a config directory specific to another module, pass the module here
    :param withUnlimited: default (None) is to allow for an unlimited axis. If a fixed length is preferred instead, that length should be passed here
    :param basetime: if an explicit basetime is preferred over one that is discovered, pass it here.

    additional named args are passed to dpl_create_templatenetcdf.dpl_create_templatenetcdf
    """
    def __init__(self,framestream,template,outputfilename=None,format=None,usecfradial=None,
        selected_bindings=None,output=None, 
        basetime=None,**kwargs):
        super(dpl_netcdf_artist,self).__init__(framestream)
        self.template=template
        self.outputfilename=outputfilename
        #self.group=group
        #if forModule==None:
        #    fm=None
        #else:
        #    fm=[sys.modules[__name__]]
        #    if isinstance(forModule,types.ModuleType):
        #        fm.insert(0,forModule)
        #    else:
        #        for m in forModule:
        #            fm.insert(len(fm)-1,m)
        #self.modules=fm
        self.cfradial=usecfradial
        if output==None and format==None:
            format='NETCDF4'
            if '3' in template:
                format='NETCDF3_64BIT'
        self.format=format
        self.outnetcdf=output
        self.willclose=(output==None)
        self.basetime=basetime
        self.usecfradial=False #deprecated
        #self.attributes=addAttributes
        self.nctemplate=None
        self.framestream=framestream
        self.selectedVariables=selected_bindings
        #self.withUnlimited=withUnlimited
        self.templateargs=kwargs.copy()
        if 'cfradial' not in self.templateargs:
            self.templateargs['cfradial']=usecfradial
        if 'withBasetime' not in self.templateargs:
            self.templateargs['withBasetime']=basetime
        if self.acceptableMetaframe(framestream.provides):
            self.__opentemplate(framestream.provides)

    def acceptableMetaframe(self,metaframe):
        if not isinstance(metaframe,dict):
            return False
        if 'shape' in metaframe:
            return True
        keys=metaframe.keys()
        if len(keys)==0:
            return False
        return self.acceptableMetaframe(metaframe[keys[0]])

    def __opentemplate(self,metaframe):
        if self.outnetcdf==None:
            from netCDF4 import Dataset
            self.outnetcdf=Dataset(self.outputfilename,'w',clobber=True,format=self.format)
        selvar=self.selectedVariables if (self.selectedVariables!=None and len(self.selectedVariables)>0) else None
        if self.usecfradial:#deprecated
            import dpl_create_cfradial as dpl_ctnc
            self.nctemplate=dpl_ctnc.DplCreateCfradial(self.template,self.outnetcdf,metaframe,bindings_to_keep=selvar)
        else:
            import dpl_create_templatenetcdf as dpl_ctnc
            self.nctemplate=dpl_ctnc.dpl_create_templatenetcdf(self.template,self.outnetcdf,metaframe#,cfradial=self.cfradial
                ,bindings_to_keep=selvar,**self.templateargs)#withUnlimited=self.withUnlimited,withBasetime=self.basetime,addAttributes=self.attributes,
#                group=self.group,forModule=self.modules)

    def render(self):
        for frame in self.framestream:
            if frame==None:
                yield frame
                continue
            if self.nctemplate==None:
                self.__opentemplate(frame)
            if self.usecfradial:#deprecated
                self.nctemplate.append_data(frame)
            else:
                self.nctemplate.appendtemplatedata(frame)
            self.outnetcdf.sync()
            yield frame

    def __del__(self):
        # oddly enough, __del__ can be called with an incomplete object.
        # don't assume any attributes actually exist
        if  hasattr(self, 'nctemplate') and self.nctemplate!=None:
            self.nctemplate.close()
            
        if hasattr(self, 'outnetcdf') and self.outnetcdf and self.willclose:
            self.outnetcdf.close()



if __name__ == '__main__':
    for w,s,e in multi_netcdf_filewindow('start','end',datetime(2013,1,9,0,0,0),datetime(2013,5,1,1,0,0),'day'):
        print w
