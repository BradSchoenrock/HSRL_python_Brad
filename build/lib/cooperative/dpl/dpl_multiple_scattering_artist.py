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

class dpl_multiple_scattering_artist(dplkit.role.artist.aArtist):
    """ Artist for rendering multiple scattering Images

    :param framestream: iterable dpl frame stream
    :param display_defaults: visual configuration file
    """
    def __init__(self,framestream,display_defaults,subframe='rs_multiple_scattering',figurecontainer=None,includenestedframes={}):
        super(dpl_multiple_scattering_artist,self).__init__(framestream)
        self.display_defaults=display_defaults
        if isinstance(display_defaults,basestring):
            self.display_defaults = jc.json_config(lf.locate_file(display_defaults),'display_defaults',allow_missing_values=True)
        self.figcontainer=figurecontainer
        try:
            self.instrument=getattr(self.framestream,'multiple_scattering_instrument')
        except:
            self.instrument='Multiple Scattering'
        self.subframename=subframe
        self.framestream=framestream
        self.includenestedframes=includenestedframes
        import cooperative.graphics.multiple_scattering_display as msd
        self.msd=msd
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
            nativeframe=rs
            if self.subframename==None:
                if usetimes==None and hasattr(nativeframe,'times'):
                    usetimes=nativeframe.times
                if usealts==None and hasattr(nativeframe,'msl_altitudes'):
                    usealts=nativeframe.msl_altitudes
                if usetimes==None or usealts==None:
                    print "Can't find native mscat axes in artist."
            if self.subframename!=None and hasattr(rs,self.subframename):
                nativeframe=getattr(rs,self.subframename)
                if usealts==None and hasattr(nativeframe,'msl_altitudes'):
                    usealts=nativeframe.msl_altitudes
                if usetimes==None and hasattr(nativeframe,'times'):
                    usetimes=nativeframe.times
                if usetimes==None or usealts==None:
                    print "Can't find native mscat",self.subframename,"range axes in artist."
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
            self.rendernativeframe(nativeframe,usetimes,usealts,wholeframe=rs)

            self.figcontainer.shownew()
 
            if allc:
                del myframe

    def rendernativeframe(self,rs,usetimes,usealts,wholeframe=None):
            #toplevel
            if usealts==None and hasattr(rs,'msl_altitudes'):
                usealts=rs.msl_altitudes
            if usetimes==None and hasattr(rs,'times'):
                usetimes=rs.times

            gt=self.gt
            figs=self.figcontainer
            display_defaults=self.display_defaults
            instrument=self.instrument

            additionalparameters={}
            if wholeframe and self.includenestedframes:
                for paramname,framename in self.includenestedframes.items():
                    additionalparameters[paramname]=getattr(wholeframe,framename)
            try:
                pp=self.spheroid_particle_parameters
            except:
                pp=None

            self.msd.show_multiple_scattering(self.instrument,display_defaults,rs,self.multiple_scattering_parameters,pp,usetimes,usealts,figs,**additionalparameters)


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
