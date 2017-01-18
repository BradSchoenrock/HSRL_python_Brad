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


class dpl_radar_images_artist(dplkit.role.artist.aArtist):
    """ Artist for rendering Radar Images

    The busy work on this artist is all in hsrl.data_stream.radar_display.show_radar()

    :param framestream: iterable dpl frame stream
    :param display_defaults: visual configuration file
    """
    def __init__(self,framestream,display_defaults,instrument=None,figurecontainer=None,subframe='rs_mmcr'):
        super(dpl_radar_images_artist,self).__init__(framestream)
        self.display_defaults=display_defaults
        if isinstance(display_defaults,basestring):
            self.display_defaults = jc.json_config(lf.locate_file(display_defaults),'display_defaults',allow_missing_values=True)
        self.figcontainer=figurecontainer 
        if instrument is None:
            print 'WARNING: RADAR Type not specified at init! Might not be right if multiple radar streams given'
        self.instrument=instrument or (framestream.radarType if hasattr(framestream,'radarType') else 'Radar')
        self.framestream=framestream
        self.subframe=subframe
        import radar.graphics.radar_display as rd
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
            nativeframe=rs
            if self.subframe is None:
                if usetimes is None and hasattr(nativeframe,'times'):
                    usetimes=nativeframe.times
                if usetimes is not None and usetimes.size==0:
                    usetimes=None
                if usealts is None and hasattr(nativeframe,'heights'):
                    usealts=nativeframe.heights
                if usetimes is None or usealts is None:
                    print "Can't find native radar axes in artist."
            if self.subframe is not None:
                if not hasattr(rs,self.subframe):
                    print 'Subframe ',self.subframe,'is missing!!!'
                    return
                nativeframe=getattr(rs,self.subframe)
                if usealts is None and hasattr(nativeframe,'heights'):
                    usealts=nativeframe.heights
                if usetimes is None and hasattr(nativeframe,'times'):
                    usetimes=nativeframe.times
                if usetimes is not None and usetimes.size==0:
                    usetimes=None
                if usetimes is None or usealts==None:
                    print "Can't find native Radar range axes in artist."
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
            self.rendernativeframe(nativeframe,usetimes,usealts)

            self.figcontainer.shownew()
 
            if allc:
                del myframe

    def rendernativeframe(self,rs,usetimes,usealts):
            #toplevel
            if usealts==None and hasattr(rs,'heights'):
                usealts=rs.heights
            if usetimes==None and hasattr(rs,'times'):
                usetimes=rs.times

            gt=self.gt
            figs=self.figcontainer
            display_defaults=self.display_defaults
            instrument=self.instrument
            self.rd.show_radar(instrument,display_defaults,rs,usetimes,usealts,figs)

    def render(self):
        #if self.figcontainer!=None:
        #    self.figcontainer.close()
        #    self.figcontainer=None

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
