import dplkit.role.artist
import dplkit.role.filter
import dplkit.role.decorator
import copy
from datetime import datetime,timedelta
import calendar
import os
import traceback
import lg_base.core.json_config as jc
import lg_base.core.locate_file as lf


class dpl_raman_images_artist(dplkit.role.artist.aArtist):
    """ Artist for rendering Raman Lidar Images

    The busy work on this artist is all in hsrl.data_stream.radar_display.show_radar()

    :param framestream: iterable dpl frame stream
    :param display_defaults: visual configuration file
    """
    def __init__(self,framestream,display_defaults,figurecontainer=None,subframe=None,streamname=None):
        super(dpl_raman_images_artist,self).__init__(framestream)
        self.display_defaults=display_defaults
        if isinstance(display_defaults,basestring):
            self.display_defaults = jc.json_config(lf.locate_file(display_defaults),'display_defaults',allow_missing_values=True)
        self.figcontainer=figurecontainer 
        self.instrument=self.platform+'-'+(streamname or self.ramanType)# if hasattr(framestream,'radarType') else 'Radar')
        self.framestream=framestream
        self.subframe=subframe
        import raman.graphics.raman_display as rd
        self.rd=rd
        import lg_base.graphics.graphics_toolkit as gt
        self.gt=gt
        #if isinstance(self.display_defaults,basestring):
        #    [self.display_defaults, dummy]= du.get_display_defaults(self.display_defaults,"new")
        self.stepper=None
        if self.figcontainer is None:
            self.figcontainer=self.gt.figurelist()

    def replace(self,framestream=None):
        if framestream is not None:
            self.framestream=framestream

    def renderframe(self,frame):
            print 'radar render frame called'
            myframe=frame
            allc=False
            #self.du.show_images(instrument=self.instrument,rs=myframe,
            #                processing_defaults=self.processing_defaults,
            #                display_defaults=self.display_defaults,
            #                max_alt=self.max_alt,auto_loop=None,figlist=self.figcontainer,**sondeparms)
            rs=myframe
            usetimes=None
            usealts=None
            nativeframe=rs
            if self.subframe is None:
                if usetimes is None and hasattr(nativeframe,'times'):
                    usetimes=nativeframe.times
                if usealts is None and hasattr(nativeframe,'altitudes'):
                    usealts=nativeframe.altitudes
                if usetimes is None or usealts is None:
                    print "Can't find native radar axes in artist."
            if self.subframe is not None:
                if not hasattr(rs,self.subframe):
                    print 'Subframe ',self.subframe,'is missing!!!'
                    return
                nativeframe=getattr(rs,self.subframe)
                if usealts is None and hasattr(nativeframe,'altitudes'):
                    usealts=nativeframe.altitudes
                if usetimes is None and hasattr(nativeframe,'times'):
                    usetimes=nativeframe.times
                if usetimes is None or usealts is None:
                    print "Can't find native Radar range axes in artist."
            if hasattr(rs,'rs_inv'):
                if usetimes is None and hasattr(rs.rs_inv,'times'):
                    usetimes=rs.rs_inv.times
                if usealts is None and hasattr(rs.rs_inv,'msl_altitudes'):
                    usealts=rs.rs_inv.msl_altitudes
            if hasattr(rs,'rs_mean'):
                if usetimes is None and hasattr(rs.rs_mean,'times'):
                    usetimes=rs.rs_mean.times
                if usealts is None and hasattr(rs.rs_mean,'msl_altitudes'):
                    usealts=rs.rs_mean.msl_altitudes

            #print nativeframe
            #print vars(nativeframe).keys()
            self.rendernativeframe(nativeframe,usetimes,usealts)

            self.figcontainer.shownew()
 
            if allc:
                del myframe

    def rendernativeframe(self,rs,usetimes,usealts):
            #toplevel
            if usealts is None and hasattr(rs,'altitudes'):
                usealts=rs.altitudes
            if usetimes is None and hasattr(rs,'times'):
                usetimes=rs.times

            gt=self.gt
            figs=self.figcontainer
            display_defaults=self.display_defaults
            instrument=self.instrument
            self.rd.show_raman(instrument,display_defaults,rs,usetimes,usealts,figs,consts=self.raman_constants_first if hasattr(self,'raman_constants_first') else None)

    def render(self):
        #if self.figcontainer!=None:
        #    self.figcontainer.close()
        #    self.figcontainer=None

        for frame in self.framestream:
            if frame is not None:
                self.renderframe(frame)
            yield frame

    @property
    def figs(self):
        return self.figcontainer

    def __iter__(self):
        return self.render()

    def __call__(self):
        for f in self:
            print 'loop'
            pass
        return self.figs

    def step(self):
        if self.stepper is None:
            self.stepper=self.render()
        try:
            r=self.stepper.next()
        except StopIteration:
            self.stepper=None
            r=None
        return r



class dpl_ramanmerge_images_artist(dplkit.role.artist.aArtist):
    """ Artist for rendering Raman Lidar Images

    The busy work on this artist is all in hsrl.data_stream.radar_display.show_radar()

    :param framestream: iterable dpl frame stream
    :param display_defaults: visual configuration file
    """
    def __init__(self,framestream,display_defaults,figurecontainer=None,subframe=None,streamname=None):
        super(dpl_ramanmerge_images_artist,self).__init__(framestream)
        self.display_defaults=display_defaults
        if isinstance(display_defaults,basestring):
            self.display_defaults = jc.json_config(lf.locate_file(display_defaults),'display_defaults',allow_missing_values=True)
        self.figcontainer=figurecontainer 
        self.instrument=self.platform+'-'+(streamname or self.ramanType)#self.ramanType
        self.framestream=framestream
        self.subframe=subframe
        import raman.graphics.raman_display as rd
        self.rd=rd
        import lg_base.graphics.graphics_toolkit as gt
        self.gt=gt
        #if isinstance(self.display_defaults,basestring):
        #    [self.display_defaults, dummy]= du.get_display_defaults(self.display_defaults,"new")
        self.stepper=None
        if self.figcontainer is None:
            self.figcontainer=self.gt.figurelist()

    def replace(self,framestream=None):
        if framestream is not None:
            self.framestream=framestream

    def renderframe(self,frame):
            print 'radar render frame called'
            myframe=frame
            allc=False
            #self.du.show_images(instrument=self.instrument,rs=myframe,
            #                processing_defaults=self.processing_defaults,
            #                display_defaults=self.display_defaults,
            #                max_alt=self.max_alt,auto_loop=None,figlist=self.figcontainer,**sondeparms)
            rs=myframe
            usetimes=None
            usealts=None
            nativeframe=rs
            if self.subframe is None:
                if usetimes is None and hasattr(nativeframe,'times'):
                    usetimes=nativeframe.times
                if usealts is None and hasattr(nativeframe,'altitudes'):
                    usealts=nativeframe.altitudes
                if usetimes is None or usealts is None:
                    print "Can't find native radar axes in artist."
            if self.subframe is not None:
                if not hasattr(rs,self.subframe):
                    print 'Subframe ',self.subframe,'is missing!!!'
                    return
                nativeframe=getattr(rs,self.subframe)
                if usealts is None and hasattr(nativeframe,'altitudes'):
                    usealts=nativeframe.altitudes
                if usetimes is None and hasattr(nativeframe,'times'):
                    usetimes=nativeframe.times
                if usetimes is None or usealts is None:
                    print "Can't find native Radar range axes in artist."
            if hasattr(rs,'rs_inv'):
                if usetimes is None and hasattr(rs.rs_inv,'times'):
                    usetimes=rs.rs_inv.times
                if usealts is None and hasattr(rs.rs_inv,'msl_altitudes'):
                    usealts=rs.rs_inv.msl_altitudes
            if hasattr(rs,'rs_mean'):
                if usetimes is None and hasattr(rs.rs_mean,'times'):
                    usetimes=rs.rs_mean.times
                if usealts is None and hasattr(rs.rs_mean,'msl_altitudes'):
                    usealts=rs.rs_mean.msl_altitudes

            #print nativeframe
            #print vars(nativeframe).keys()
            self.rendernativeframe(nativeframe,usetimes,usealts)

            self.figcontainer.shownew()
 
            if allc:
                del myframe

    def rendernativeframe(self,rs,usetimes,usealts):
            #toplevel
            if usealts is None and hasattr(rs,'altitudes'):
                usealts=rs.altitudes
            if usetimes is None and hasattr(rs,'times'):
                usetimes=rs.times

            gt=self.gt
            figs=self.figcontainer
            display_defaults=self.display_defaults
            instrument=self.instrument
            self.rd.show_ramanmerge(instrument,display_defaults,rs,usetimes,usealts,figs,consts=self.raman_constants_first if hasattr(self,'raman_constants_first') else None)

    def render(self):
        #if self.figcontainer!=None:
        #    self.figcontainer.close()
        #    self.figcontainer=None

        for frame in self.framestream:
            if frame is not None:
                self.renderframe(frame)
            yield frame

    @property
    def figs(self):
        return self.figcontainer

    def __iter__(self):
        return self.render()

    def __call__(self):
        for f in self:
            print 'loop'
            pass
        return self.figs

    def step(self):
        if self.stepper is None:
            self.stepper=self.render()
        try:
            r=self.stepper.next()
        except StopIteration:
            self.stepper=None
            r=None
        return r


class dpl_raman_profile_images_artist(dplkit.role.artist.aArtist):
    """ Artist for rendering Raman Lidar Images

    The busy work on this artist is all in hsrl.data_stream.radar_display.show_radar()

    :param framestream: iterable dpl frame stream
    :param display_defaults: visual configuration file
    """
    def __init__(self,framestream,display_defaults,figurecontainer=None,subframe=None,streamname=None):
        super(dpl_raman_profile_images_artist,self).__init__(framestream)
        self.display_defaults=display_defaults
        if isinstance(display_defaults,basestring):
            self.display_defaults = jc.json_config(lf.locate_file(display_defaults),'display_defaults',allow_missing_values=True)
        self.figcontainer=figurecontainer 
        self.instrument=self.platform+'-'+(streamname or self.ramanType)#self.ramanType
        self.framestream=framestream
        self.subframe=subframe
        import raman.graphics.raman_display as rd
        self.rd=rd
        import lg_base.graphics.graphics_toolkit as gt
        self.gt=gt
        #if isinstance(self.display_defaults,basestring):
        #    [self.display_defaults, dummy]= du.get_display_defaults(self.display_defaults,"new")
        self.stepper=None
        if self.figcontainer is None:
            self.figcontainer=self.gt.figurelist()

    def replace(self,framestream=None):
        if framestream is not None:
            self.framestream=framestream

    def renderframe(self,frame):
            print 'radar render frame called'
            myframe=frame
            allc=False
            #self.du.show_images(instrument=self.instrument,rs=myframe,
            #                processing_defaults=self.processing_defaults,
            #                display_defaults=self.display_defaults,
            #                max_alt=self.max_alt,auto_loop=None,figlist=self.figcontainer,**sondeparms)
            rs=myframe
            usetimes=None
            usedeltas=None
            usealts=None
            nativeframe=rs
            if self.subframe is None:
                if usetimes is None and hasattr(nativeframe,'times'):
                    usetimes=nativeframe.times
                if usetimes is None and hasattr(nativeframe,'delta_t'):
                    usedeltas=nativeframe.delta_t
                if usealts is None and hasattr(nativeframe,'altitudes'):
                    usealts=nativeframe.altitudes
                if usetimes is None or usealts is None:
                    print "Can't find native radar axes in artist."
            if self.subframe is not None:
                if not hasattr(rs,self.subframe):
                    print 'Subframe ',self.subframe,'is missing!!!'
                    return
                nativeframe=getattr(rs,self.subframe)
                if usealts is None and hasattr(nativeframe,'altitudes'):
                    usealts=nativeframe.altitudes
                if usetimes is None and hasattr(nativeframe,'delta_t'):
                    usedeltas=nativeframe.delta_t
                if usetimes is None and hasattr(nativeframe,'times'):
                    usetimes=nativeframe.times
                if usetimes is None or usealts is None:
                    print "Can't find native Radar range axes in artist."

            #print nativeframe
            #print vars(nativeframe).keys()
            self.rendernativeframe(nativeframe,usetimes,usedeltas,usealts)

            self.figcontainer.shownew()
 
            if allc:
                del myframe

    def rendernativeframe(self,rs,usetimes,usedeltas,usealts):
            #toplevel
            if usealts is None and hasattr(rs,'altitudes'):
                usealts=rs.altitudes
            if usetimes is None and hasattr(rs,'times'):
                usetimes=rs.times
            if usedeltas is None and hasattr(rs,'delta_t'):
                usedeltas=rs.delta_t
            #if usetimes is not None and usedeltas is not None:
            #    usetimes=[usetimes,usetimes+usedeltas]

            gt=self.gt
            figs=self.figcontainer
            display_defaults=self.display_defaults
            instrument=self.instrument
            self.rd.show_raman_profile(instrument,display_defaults,rs,usetimes,usealts,figs,consts=self.raman_constants_first if hasattr(self,'raman_constants_first') else None)

    def render(self):
        #if self.figcontainer!=None:
        #    self.figcontainer.close()
        #    self.figcontainer=None

        for frame in self.framestream:
            if frame is not None:
                self.renderframe(frame)
            yield frame

    @property
    def figs(self):
        return self.figcontainer

    def __iter__(self):
        return self.render()

    def __call__(self):
        for f in self:
            print 'loop'
            pass
        return self.figs

    def step(self):
        if self.stepper is None:
            self.stepper=self.render()
        try:
            r=self.stepper.next()
        except StopIteration:
            self.stepper=None
            r=None
        return r


class dpl_raman_inverted_profile_images_artist(dplkit.role.artist.aArtist):
    """ Artist for rendering Raman Lidar Images

    The busy work on this artist is all in hsrl.data_stream.radar_display.show_radar()

    :param framestream: iterable dpl frame stream
    :param display_defaults: visual configuration file
    """
    def __init__(self,framestream,display_defaults,figurecontainer=None,subframe=None,streamname=None):
        super(dpl_raman_inverted_profile_images_artist,self).__init__(framestream)
        self.display_defaults=display_defaults
        if isinstance(display_defaults,basestring):
            self.display_defaults = jc.json_config(lf.locate_file(display_defaults),'display_defaults',allow_missing_values=True)
        self.figcontainer=figurecontainer 
        self.instrument=self.platform+'-'+(streamname or self.ramanType)#self.ramanType
        self.framestream=framestream
        self.subframe=subframe
        import raman.graphics.raman_display as rd
        self.rd=rd
        import lg_base.graphics.graphics_toolkit as gt
        self.gt=gt
        #if isinstance(self.display_defaults,basestring):
        #    [self.display_defaults, dummy]= du.get_display_defaults(self.display_defaults,"new")
        self.stepper=None
        if self.figcontainer is None:
            self.figcontainer=self.gt.figurelist()

    def replace(self,framestream=None):
        if framestream is not None:
            self.framestream=framestream

    def renderframe(self,frame):
            print 'radar render frame called'
            myframe=frame
            allc=False
            #self.du.show_images(instrument=self.instrument,rs=myframe,
            #                processing_defaults=self.processing_defaults,
            #                display_defaults=self.display_defaults,
            #                max_alt=self.max_alt,auto_loop=None,figlist=self.figcontainer,**sondeparms)
            rs=myframe
            usetimes=None
            usedeltas=None
            usealts=None
            nativeframe=rs
            if self.subframe is None:
                if usetimes is None and hasattr(nativeframe,'times'):
                    usetimes=nativeframe.times
                if usetimes is None and hasattr(nativeframe,'delta_t'):
                    usedeltas=nativeframe.delta_t
                if usealts is None and hasattr(nativeframe,'altitudes'):
                    usealts=nativeframe.altitudes
                if usetimes is None or usealts is None:
                    print "Can't find native radar axes in artist."
            if self.subframe is not None:
                if not hasattr(rs,self.subframe):
                    print 'Subframe ',self.subframe,'is missing!!!'
                    return
                nativeframe=getattr(rs,self.subframe)
                if usealts is None and hasattr(nativeframe,'altitudes'):
                    usealts=nativeframe.altitudes
                if usetimes is None and hasattr(nativeframe,'delta_t'):
                    usedeltas=nativeframe.delta_t
                if usetimes is None and hasattr(nativeframe,'times'):
                    usetimes=nativeframe.times
                if usetimes is None or usealts is None:
                    print "Can't find native Radar range axes in artist."

            #print nativeframe
            #print vars(nativeframe).keys()
            self.rendernativeframe(nativeframe,usetimes,usedeltas,usealts)

            self.figcontainer.shownew()
 
            if allc:
                del myframe

    def rendernativeframe(self,rs,usetimes,usedeltas,usealts):
            #toplevel
            if usealts is None and hasattr(rs,'altitudes'):
                usealts=rs.altitudes
            if usetimes is None and hasattr(rs,'times'):
                usetimes=rs.times
            if usedeltas is None and hasattr(rs,'delta_t'):
                usedeltas=rs.delta_t
            #if usetimes is not None and usedeltas is not None:
            #    usetimes=[usetimes,usetimes+usedeltas]

            gt=self.gt
            figs=self.figcontainer
            display_defaults=self.display_defaults
            instrument=self.instrument
            self.rd.show_raman_inverted_profile(instrument,display_defaults,rs,usetimes,usealts,figs,consts=self.raman_constants_first if hasattr(self,'raman_constants_first') else None)

    def render(self):
        #if self.figcontainer!=None:
        #    self.figcontainer.close()
        #    self.figcontainer=None

        for frame in self.framestream:
            if frame is not None:
                self.renderframe(frame)
            yield frame

    @property
    def figs(self):
        return self.figcontainer

    def __iter__(self):
        return self.render()

    def __call__(self):
        for f in self:
            print 'loop'
            pass
        return self.figs

    def step(self):
        if self.stepper is None:
            self.stepper=self.render()
        try:
            r=self.stepper.next()
        except StopIteration:
            self.stepper=None
            r=None
        return r



if __name__ == '__main__':
    for w,s,e in multi_netcdf_filewindow('start','end',datetime(2013,1,9,0,0,0),datetime(2013,5,1,1,0,0),'day'):
        print w
