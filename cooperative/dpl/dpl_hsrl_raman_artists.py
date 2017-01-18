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



class dpl_raman_hsrl_profile_images_artist(dplkit.role.artist.aArtist):
    """ Artist for rendering HSRL/Radar Cooperative Particle Images

    :param framestream: iterable dpl frame stream
    :param display_defaults: visual configuration file
    """
    def __init__(self,framestream,display_defaults,subframe=None,figurecontainer=None,includenestedframes={}):
        super(dpl_raman_hsrl_profile_images_artist,self).__init__(framestream)
        self.display_defaults=display_defaults
        if isinstance(display_defaults,basestring):
            self.display_defaults = jc.json_config(lf.locate_file(display_defaults),'display_defaults',allow_missing_values=True)
        self.figcontainer=figurecontainer
        self.subframename=subframe
        self.framestream=framestream
        self.includenestedframes=includenestedframes
        import cooperative.graphics.hsrlraman_display as rd
        self.rd=rd
        import lg_base.graphics.graphics_toolkit as gt
        self.gt=gt
        #if isinstance(self.display_defaults,basestring):
        #    [self.display_defaults, dummy]= du.get_display_defaults(self.display_defaults,"new")
        self.stepper=None
        if self.figcontainer==None:
            self.figcontainer=self.gt.figurelist()

    def replace(self,framestream=None):
        if framestream!=None:
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
            wholeframe=None
            nativeframe=rs
            if self.subframename==None:
                if usetimes==None and hasattr(nativeframe,'start'):
                    usetimes=[nativeframe.start,nativeframe.start+nativeframe.width]
                if usealts==None and hasattr(nativeframe,'altitudes'):
                    usealts=nativeframe.altitudes
                if usealts==None and hasattr(nativeframe,'heights'):
                    usealts=nativeframe.heights
                if usetimes==None or usealts==None:
                    print "Can't find native raw ramanhsrl axes in artist."
            elif self.subframename!=None:
                wholeframe=rs
                if hasattr(rs,self.subframename):
                    nativeframe=getattr(rs,self.subframename)
                else:
                    nativeframe=None
                if usealts==None and hasattr(nativeframe,'altitudes'):
                    usealts=nativeframe.altitudes
                if usealts==None and hasattr(nativeframe,'heights'):
                    usealts=nativeframe.heights
                if usetimes==None and hasattr(nativeframe,'start'):
                    usetimes=[nativeframe.start,nativeframe.start+nativeframe.width]
                if usetimes==None or usealts==None:
                    print "Can't find native raw ramanhsrl",self.subframename,"range axes in artist."
                if hasattr(rs,'profiles'):
                    if usetimes==None and hasattr(rs.profiles,'start'):
                        usetimes=[rs.profiles.start,rs.profiles.start+rs.profiles.width]
                    if usealts==None and hasattr(rs.profiles,'msl_altitudes'):
                        usealts=rs.profiles.msl_altitudes

            #print nativeframe
            #print vars(nativeframe).keys()
            try:
                if usetimes!=None and usealts!=None:
                    self.rendernativeframe(nativeframe,usetimes,usealts,wholeframe=wholeframe)

                self.figcontainer.shownew()
            except AttributeError:
                print 'failed to render merge raman/hsrl'
                traceback.print_exc()
                pass
 
            if allc:
                del myframe

    def rendernativeframe(self,rs,usetimes,usealts,wholeframe=None):
            #toplevel

            gt=self.gt
            figs=self.figcontainer
            display_defaults=self.display_defaults
            instrument=self.instrument

            additionalparameters={}
            if wholeframe and self.includenestedframes:
                for paramname,framename in self.includenestedframes.items():
                    if hasattr(wholeframe,framename):
                        additionalparameters[paramname]=getattr(wholeframe,framename)

            self.rd.show_raman_hsrl_profile(display_defaults,rs,None,usetimes,usealts,figs,**additionalparameters)


    def render(self):
        for frame in self.framestream:
            if frame!=None:
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
        if self.stepper==None:
            self.stepper=self.render()
        try:
            r=self.stepper.next()
        except StopIteration:
            self.stepper=None
            r=None
        return r


class dpl_ramanmerge_hsrl_images_artist(dplkit.role.artist.aArtist):
    """ Artist for rendering HSRL/Radar Cooperative Particle Images

    :param framestream: iterable dpl frame stream
    :param display_defaults: visual configuration file
    """
    def __init__(self,framestream,display_defaults,subframe=None,figurecontainer=None,includenestedframes={}):
        super(dpl_ramanmerge_hsrl_images_artist,self).__init__(framestream)
        self.display_defaults=display_defaults
        if isinstance(display_defaults,basestring):
            self.display_defaults = jc.json_config(lf.locate_file(display_defaults),'display_defaults',allow_missing_values=True)
        self.figcontainer=figurecontainer
        self.subframename=subframe
        self.framestream=framestream
        self.includenestedframes=includenestedframes
        import cooperative.graphics.hsrlraman_display as rd
        self.rd=rd
        import lg_base.graphics.graphics_toolkit as gt
        self.gt=gt
        #if isinstance(self.display_defaults,basestring):
        #    [self.display_defaults, dummy]= du.get_display_defaults(self.display_defaults,"new")
        self.stepper=None
        if self.figcontainer==None:
            self.figcontainer=self.gt.figurelist()

    def replace(self,framestream=None):
        if framestream!=None:
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
            wholeframe=None
            nativeframe=rs
            if self.subframename==None:
                if usetimes==None and hasattr(nativeframe,'times'):
                    usetimes=nativeframe.times
                if usealts==None and hasattr(nativeframe,'altitudes'):
                    usealts=nativeframe.altitudes
                if usealts==None and hasattr(nativeframe,'heights'):
                    usealts=nativeframe.heights
                if usetimes==None or usealts==None:
                    print "Can't find native raw ramanhsrl axes in artist."
            elif self.subframename!=None:
                wholeframe=rs
                if hasattr(rs,self.subframename):
                    nativeframe=getattr(rs,self.subframename)
                else:
                    nativeframe=None
                if usealts==None and hasattr(nativeframe,'altitudes'):
                    usealts=nativeframe.altitudes
                if usealts==None and hasattr(nativeframe,'heights'):
                    usealts=nativeframe.heights
                if usetimes==None and hasattr(nativeframe,'times'):
                    usetimes=nativeframe.times
                if usetimes==None or usealts==None:
                    print "Can't find native raw ramanhsrl",self.subframename,"range axes in artist."
                if hasattr(rs,'rs_mmcr'):
                    if usealts==None and hasattr(rs.rs_mmcr,'heights'):
                        usealts=rs.rs_mmcr.heights
                    if usetimes==None and hasattr(rs.rs_mmcr,'times'):
                        usetimes=rs.rs_mmcr.times
                if hasattr(rs,'rs_inv'):
                    if usetimes==None and hasattr(rs.rs_inv,'times'):
                        usetimes=rs.rs_inv.times
                    if usealts==None and hasattr(rs.rs_inv,'msl_altitudes'):
                        usealts=rs.rs_inv.msl_altitudes
                if hasattr(rs,'rs_mean'):
                    if usetimes==None and hasattr(rs.rs_mean,'times'):
                        usetimes=rs.rs_mean.times
                    if usealts==None and hasattr(rs.rs_mean,'msl_altitudes'):
                        usealts=rs.rs_mean.msl_altitudes

            #print nativeframe
            #print vars(nativeframe).keys()
            try:
                if usetimes!=None and usealts!=None:
                    self.rendernativeframe(nativeframe,usetimes,usealts,wholeframe=wholeframe)

                self.figcontainer.shownew()
            except AttributeError:
                print 'failed to render merge raman/hsrl'
                traceback.print_exc()
                pass
 
            if allc:
                del myframe

    def rendernativeframe(self,rs,usetimes,usealts,wholeframe=None):
            #toplevel

            gt=self.gt
            figs=self.figcontainer
            display_defaults=self.display_defaults
            instrument=self.instrument

            additionalparameters={}
            if wholeframe and self.includenestedframes:
                for paramname,framename in self.includenestedframes.items():
                    if hasattr(wholeframe,framename):
                        additionalparameters[paramname]=getattr(wholeframe,framename)

            self.rd.show_ramanmerge_hsrl(display_defaults,rs,None,usetimes,usealts,figs,**additionalparameters)


    def render(self):
        for frame in self.framestream:
            if frame!=None:
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
        if self.stepper==None:
            self.stepper=self.render()
        try:
            r=self.stepper.next()
        except StopIteration:
            self.stepper=None
            r=None
        return r

class dpl_raman_hsrl_images_artist(dplkit.role.artist.aArtist):
    """ Artist for rendering HSRL/Radar Cooperative Particle Images

    :param framestream: iterable dpl frame stream
    :param display_defaults: visual configuration file
    """
    def __init__(self,framestream,display_defaults,subframe=None,figurecontainer=None,includenestedframes={}):
        super(dpl_raman_hsrl_images_artist,self).__init__(framestream)
        self.display_defaults=display_defaults
        if isinstance(display_defaults,basestring):
            self.display_defaults = jc.json_config(lf.locate_file(display_defaults),'display_defaults',allow_missing_values=True)
        self.figcontainer=figurecontainer
        self.subframename=subframe
        self.framestream=framestream
        self.includenestedframes=includenestedframes
        import cooperative.graphics.hsrlraman_display as rd
        self.rd=rd
        import lg_base.graphics.graphics_toolkit as gt
        self.gt=gt
        #if isinstance(self.display_defaults,basestring):
        #    [self.display_defaults, dummy]= du.get_display_defaults(self.display_defaults,"new")
        self.stepper=None
        if self.figcontainer==None:
            self.figcontainer=self.gt.figurelist()

    def replace(self,framestream=None):
        if framestream!=None:
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
            wholeframe=None
            nativeframe=rs
            if self.subframename is not None:
                if usetimes==None and hasattr(nativeframe,'times'):
                    usetimes=nativeframe.times
                if usealts==None and hasattr(nativeframe,'altitudes'):
                    usealts=nativeframe.altitudes
                if usealts==None and hasattr(nativeframe,'heights'):
                    usealts=nativeframe.heights
                if usetimes==None or usealts==None:
                    print "Can't find native ramanhsrl axes in artist."
            elif self.subframename is not None:
                wholeframe=rs
                if hasattr(rs,self.subframename):
                    nativeframe=getattr(rs,self.subframename)
                else:
                    nativeframe=None
                if usealts==None and hasattr(nativeframe,'altitudes'):
                    usealts=nativeframe.altitudes
                if usealts==None and hasattr(nativeframe,'heights'):
                    usealts=nativeframe.heights
                if usetimes==None and hasattr(nativeframe,'times'):
                    usetimes=nativeframe.times
                if usetimes==None or usealts==None:
                    print "Can't find native ramanhsrl",self.subframename,"range axes in artist."
                if hasattr(rs,'rs_mmcr'):
                    if usealts==None and hasattr(rs.rs_mmcr,'heights'):
                        usealts=rs.rs_mmcr.heights
                    if usetimes==None and hasattr(rs.rs_mmcr,'times'):
                        usetimes=rs.rs_mmcr.times
                if hasattr(rs,'rs_inv'):
                    if usetimes==None and hasattr(rs.rs_inv,'times'):
                        usetimes=rs.rs_inv.times
                    if usealts==None and hasattr(rs.rs_inv,'msl_altitudes'):
                        usealts=rs.rs_inv.msl_altitudes
                if hasattr(rs,'rs_mean'):
                    if usetimes==None and hasattr(rs.rs_mean,'times'):
                        usetimes=rs.rs_mean.times
                    if usealts==None and hasattr(rs.rs_mean,'msl_altitudes'):
                        usealts=rs.rs_mean.msl_altitudes

            #print nativeframe
            #print vars(nativeframe).keys()
            try:
                #if usetimes!=None and usealts!=None:
                self.rendernativeframe(nativeframe,usetimes,usealts,wholeframe=wholeframe)

                self.figcontainer.shownew()
            except AttributeError:
                print 'failed to render hsrl/raman inverted combined form'
                traceback.print_exc()
                pass
 
            if allc:
                del myframe

    def rendernativeframe(self,rs,usetimes,usealts,wholeframe=None):
            #toplevel

            gt=self.gt
            figs=self.figcontainer
            display_defaults=self.display_defaults
            instrument=self.instrument

            additionalparameters={}
            if wholeframe and self.includenestedframes:
                for paramname,framename in self.includenestedframes.items():
                    if hasattr(wholeframe,framename):
                        additionalparameters[paramname]=getattr(wholeframe,framename)

            self.rd.show_raman_hsrl(display_defaults,rs,None,usetimes,usealts,figs,**additionalparameters)


    def render(self):
        for frame in self.framestream:
            if frame!=None:
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
        if self.stepper==None:
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
