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


class dpl_particle_images_artist(dplkit.role.artist.aArtist):
    """ Artist for rendering HSRL/Radar Cooperative Particle Images

    :param framestream: iterable dpl frame stream
    :param display_defaults: visual configuration file
    """
    def __init__(self,framestream,display_defaults,subframe=None,mode='mass_dimension',figurecontainer=None,includenestedframes={}):
        super(dpl_particle_images_artist,self).__init__(framestream)
        self.display_defaults=display_defaults
        if isinstance(display_defaults,basestring):
            self.display_defaults = jc.json_config(lf.locate_file(display_defaults),'display_defaults',allow_missing_values=True)
        self.figcontainer=figurecontainer
        self.subframename=subframe
        self.framestream=framestream
        self.includenestedframes=includenestedframes
        self.mode=mode
        if mode not in ('spheroid','mass_dimension'):
            raise NotImplementedError('Particle image generation mode for '+mode)
        instname=mode+'_instrument'
        if hasattr(self,instname):
            self.instrument=getattr(self,instname)
        else:
            self.instrument=mode
        import cooperative.graphics.coop_display as rd
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
            if self.subframename is None:
                if usetimes is None and hasattr(nativeframe,'times'):
                    usetimes=nativeframe.times
                if usetimes is not None and usetimes.size==0:
                    usetimes=None
                if usealts is None and hasattr(nativeframe,'msl_altitudes'):
                    usealts=nativeframe.msl_altitudes
                if usealts is None and hasattr(nativeframe,'heights'):
                    usealts=nativeframe.heights
                if usetimes is None or usealts is None:
                    print "Can't find native particle axes in artist."
            elif self.subframename is not None and hasattr(rs,self.subframename):
                wholeframe=rs
                nativeframe=getattr(rs,self.subframename)
                if usealts is None and hasattr(nativeframe,'msl_altitudes'):
                    usealts=nativeframe.msl_altitudes
                if usealts is None and hasattr(nativeframe,'heights'):
                    usealts=nativeframe.heights
                if usetimes is None and hasattr(nativeframe,'times'):
                    usetimes=nativeframe.times
                if usetimes is not None and usetimes.size==0:
                    usetimes=None
                if usetimes is None or usealts is None:
                    print "Can't find native Particle",self.subframename,"range axes in artist."
                if hasattr(rs,'rs_mmcr'):
                    if usealts is None and hasattr(rs.rs_mmcr,'heights'):
                        usealts=rs.rs_mmcr.heights
                    if usetimes is None and hasattr(rs.rs_mmcr,'times'):
                        usetimes=rs.rs_mmcr.times
                    if usetimes is not None and usetimes.size==0:
                        usetimes=None
                if hasattr(rs,'rs_inv'):
                    if usetimes is None and hasattr(rs.rs_inv,'times'):
                        usetimes=rs.rs_inv.times
                    if usetimes is not None and usetimes.size==0:
                        usetimes=None
                    if usealts is None and hasattr(rs.rs_inv,'msl_altitudes'):
                        usealts=rs.rs_inv.msl_altitudes
                if hasattr(rs,'rs_mean'):
                    if usetimes is None and hasattr(rs.rs_mean,'times'):
                        usetimes=rs.rs_mean.times
                    if usetimes is not None and usetimes.size==0:
                        usetimes=None
                    if usealts is None and hasattr(rs.rs_mean,'msl_altitudes'):
                        usealts=rs.rs_mean.msl_altitudes

            #print nativeframe
            #print vars(nativeframe).keys()
            try:
                if usetimes is not None and usealts is not None:
                    self.rendernativeframe(nativeframe,usetimes,usealts,wholeframe=wholeframe)

                self.figcontainer.shownew()
            except AttributeError:
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

            if self.mode=='mass_dimension':
                self.rd.show_mass_dimension_particle(self.instrument,display_defaults,rs,self.mass_dimension_particle_parameters,usetimes,usealts,figs)
            elif self.mode=='spheroid':
                self.rd.show_spheroid_particle(self.instrument,display_defaults,rs,self.spheroid_particle_parameters,usetimes,usealts,self.radarLambda,figs,**additionalparameters)
            else:
                raise NotImplementedError(self.mode)


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
