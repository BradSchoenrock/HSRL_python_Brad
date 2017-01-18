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

class dpl_profiles_images_artist(dplkit.role.artist.aArtist):
    """ Artist for rendering HSRL/Radar Cooperative Particle Images

    :param framestream: iterable dpl frame stream
    :param display_defaults: visual configuration file
    """
    def __init__(self,framestream,display_defaults,subframes,figurecontainer=None):
        super(dpl_profiles_images_artist,self).__init__(framestream)
        self.display_defaults=display_defaults
        if isinstance(display_defaults,basestring):
            self.display_defaults = jc.json_config(lf.locate_file(display_defaults),'display_defaults',allow_missing_values=True)
        self.figcontainer=figurecontainer
        self.subframenames=subframes
        self.framestream=framestream
        import cooperative.graphics.profiles_display as rd
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

            r=dict(wholeframe=frame)
            for x in self.subframenames:
                if hasattr(frame,x):
                    r[x]=getattr(frame,x)
            #print nativeframe
            #print vars(nativeframe).keys()
            try:
                #if usetimes!=None and usealts!=None:
                self.rendernativeframe(**r)

                self.figcontainer.shownew()
            except AttributeError:
                print 'failed to render all profiles form'
                traceback.print_exc()
                pass

    def rendernativeframe(self,*args,**kwargs):
            #toplevel

            gt=self.gt
            figs=self.figcontainer
            display_defaults=self.display_defaults
            instrument=self.instrument

            self.rd.show_all_profiles(display_defaults,figs,*args,**kwargs)


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
