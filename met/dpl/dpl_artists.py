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

class dpl_met_images_artist(dplkit.role.artist.aArtist):
    """ Artist for rendering met Images

    :param framestream: iterable dpl frame stream
    :param display_defaults: visual configuration file
    """
    def __init__(self,framestream,display_defaults,framedomain=None,instrument=None,figurecontainer=None):
        super(dpl_met_images_artist,self).__init__(framestream)
        self.display_defaults=display_defaults
        if isinstance(display_defaults,basestring):
            self.display_defaults = jc.json_config(lf.locate_file(display_defaults),'display_defaults',allow_missing_values=True)
        self.figcontainer=figurecontainer
        self.instrument=instrument if instrument!=None else 'Met'
        self.framestream=framestream
        self.framedomain=framedomain
        import met.graphics.met_display as rd
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
            print 'met render frame called'
            myframe=frame
            allc=False
            #self.du.show_images(instrument=self.instrument,rs=myframe,
            #                processing_defaults=self.processing_defaults,
            #                display_defaults=self.display_defaults,
            #                max_alt=self.max_alt,auto_loop=None,figlist=self.figcontainer,**sondeparms)
            rs=myframe
            usetimes=None
            nativeframe=rs
            if self.framedomain is None:
                if usetimes==None and hasattr(nativeframe,'times'):
                    usetimes=nativeframe.times
                if usetimes==None:
                    print "Can't find native met axes in artist."
            if self.framedomain is not None:
                if not hasattr(rs,self.framedomain):
                    return
                nativeframe=getattr(rs,self.framedomain)
                if usetimes==None and hasattr(nativeframe,'times'):
                    usetimes=nativeframe.times
                if usetimes==None:
                    print "Can't find native met axes in artist."
            if hasattr(rs,'rs_particle'):
                if usetimes==None and hasattr(rs.rs_particle,'times'):
                    usetimes=rs.rs_particle.times
            if hasattr(rs,'rs_mmcr'):
                if usetimes==None and hasattr(rs.rs_mmcr,'times'):
                    usetimes=rs.rs_mmcr.times
            if hasattr(rs,'rs_inv'):
                if usetimes==None and hasattr(rs.rs_inv,'times'):
                    usetimes=rs.rs_inv.times
            if hasattr(rs,'rs_mean'):
                if usetimes==None and hasattr(rs.rs_mean,'times'):
                    usetimes=rs.rs_mean.times
 
            #print nativeframe
            #print vars(nativeframe).keys()
            self.rendernativeframe(nativeframe,usetimes)

            self.figcontainer.shownew()

            if allc:
                del myframe

    def rendernativeframe(self,rs,usetimes):
            #toplevel
            if usetimes==None and hasattr(rs,'times'):
                usetimes=rs.times

            gt=self.gt
            figs=self.figcontainer
            display_defaults=self.display_defaults
            instrument=self.instrument

            self.rd.show_met(display_defaults,rs,usetimes,figs)


    def render(self):
        #if self.figcontainer!=None:
        #    self.figcontainer.close()
        #    self.figcontainer=None
        #self.figcontainer=figcontainer #FIXME
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
