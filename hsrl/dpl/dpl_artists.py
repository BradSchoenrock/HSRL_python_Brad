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

@dplkit.role.decorator.exposes_attrs_in_chain(['figs','display_defaults']) #this is really only useful if theres only one artist
class dpl_images_artist(dplkit.role.artist.aArtist):
    """ Artist for rendering HSRL Images

    :param framestream: iterable dpl frame stream
    :param instrument: instrument name
    :param max_alt: maximum altitude
    :param processing_defaults: processing parameters used with the frame stream
    :param display_defaults: visual configuration file
    """
    def __init__(self,framestream,max_alt,display_defaults,instrument=None
        ,processing_defaults=None,figurecontainer=None,limit_frame_to=None
        ,breakup_nesting=False,flat_frame=False,enable_masking=None):
        super(dpl_images_artist,self).__init__(framestream)
        self.framestream=framestream
        #self.provides=framestream.provides
        self.instrument=instrument or framestream.hsrl_instrument
        self.max_alt=max_alt
        self.enable_masking=enable_masking
        try:
            self.processing_defaults=processing_defaults or framestream.hsrl_process_control
        except:
            self.processing_defaults=None
        self.display_defaults=display_defaults
        if isinstance(display_defaults,basestring):
            self.display_defaults = jc.json_config(lf.locate_file(display_defaults),'display_defaults',allow_missing_values=True)
        self.limit_frame_to=limit_frame_to
        self.breakup_nesting=breakup_nesting
        self.flat_frame=flat_frame
        import lg_base.graphics.graphics_toolkit as gt
        self.gt=gt
        #import hsrl.data_stream.display_utilities as du
        import hsrl.graphics.hsrl_display as du
        self.du=du
        self.stepper=None
        self.figcontainer=figurecontainer or self.gt.figurelist()

    def replace(self,framestream=None):
        if framestream!=None:
            self.framestream=framestream

    def renderframe(self,frame,wholeframe=None):
            sondeparms={'sounding':None,'geo_corr':None,'rs_constants':None,'last_sounding_time':None,'rs_Cxx':None}
            myframe=frame

            #if hasattr(myframe,'profiles') and (not hasattr(myframe.profiles,'inv') or not hasattr(myframe.profiles.inv,'times') or myframe.profiles.inv.times.shape[0]==0):
            #    myframe=copy.copy(frame)
            #    delattr(myframe,'profiles')
            if hasattr(myframe,'profiles'):#only include sounding and geometry for profiles
                    if hasattr(self,'hsrl_Cxx'):
                        sondeparms['rs_Cxx']=self.hsrl_Cxx
                    elif hasattr(wholeframe,'rs_Cxx'):
                        sondeparms['rs_Cxx']=wholeframe.rs_Cxx
                    if hasattr(self,'hsrl_sounding'):
                        sondeparms['sounding']=self.hsrl_sounding
                    elif hasattr(wholeframe,'sounding'):
                        sondeparms['sounding']=wholeframe.sounding
                    if hasattr(self,'hsrl_cal'):
                        sondeparms['geo_corr']=self.hsrl_cal.geo.data
        
            if sondeparms['rs_constants'] is None:
                if hasattr(self,'hsrl_constants_first'):
                    sondeparms['rs_constants']=self.hsrl_constants_first
                elif hasattr(wholeframe,'rs_constants'):
                    sondeparms['rs_constants']=vars(wholeframe.rs_constants)
                else:
                    sondeparms['rs_constants']={}
            self.du.show_images(instrument=self.instrument,rs=myframe,
                            processing_defaults=self.processing_defaults,
                            display_defaults=self.display_defaults,
                            max_alt=self.max_alt,auto_loop=None,figlist=self.figcontainer,
                            enable_masking=self.enable_masking,
                            **sondeparms)
            self.figcontainer.reorder(self.display_defaults.get_attrs())


    def render(self):
        subframes=None
        for frame in self.framestream:
            if frame!=None:
                if self.limit_frame_to!=None or self.breakup_nesting:
                    if subframes is None:
                        if self.limit_frame_to is not None:
                            subframes=[self.limit_frame_to]
                        else:
                            subframes=[]
                            omit=None
                            include=None
                            if hasattr(frame,'rs_mean') or hasattr(frame,'rs_inv') or not hasattr(frame,'rs_raw'):
                                omit='raw'
                            else:
                                include='raw'
                            for x in vars(frame).keys():
                                if x.startswith('_'):
                                    continue
                                if omit is not None and omit in x:
                                    continue
                                if include is not None and include not in x:
                                    continue
                                subframes.append(x)
                    for k in subframes:
                        if hasattr(frame,k):
                            if self.flat_frame:
                                print 'Showing images for HSRL flat frame '+k
                                self.renderframe(getattr(frame,k),frame)
                            else:
                                import lg_base.core.array_utils as hau
                                sendframe=hau.Time_Z_Group(like=frame)
                                setattr(sendframe,k,getattr(frame,k))
                                print 'Showing images for HSRL subframe '+k
                                if k=='rs_raw':
                                    setattr(sendframe,'rs_mean',getattr(frame,k))
                                elif k=='rs_mean':
                                    setattr(sendframe,'rs_raw',getattr(frame,k))
                                self.renderframe(sendframe,frame)
                else:
                    self.renderframe(frame,frame)
                self.figcontainer.shownew()
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
