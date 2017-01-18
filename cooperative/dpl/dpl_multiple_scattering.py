
import dplkit.role.filter
import dplkit.role.decorator
import copy
from datetime import datetime,timedelta
import lg_base.core.array_utils as hau
import numpy as np
import os
import json
import traceback

from collections import namedtuple

@dplkit.role.decorator.exposes_attrs_in_chain(['multiple_scattering_parameters'])
@dplkit.role.decorator.autoprovidenested(nestedclasses=[hau.Time_Z_Group,hau.rs_xfer,dict])
class dpl_multiple_scattering(dplkit.role.filter.aFilter):
    """ HSRL/Particle multiple scattering frame adding filter. requires both scopes be on the same resolution/gridsize

    adds rs_multiple_scattering to frame

    :param framestream: input iterable framestream
    :param hsrlscope: nested subframe source for HSRL inverted data
    :param radarscope: nested subframe source for Radar
    :param outscope: change rs_particle to another subframe
    additional parameters are passed to the initializer of cooperative.core.multiple_scattering.multiple_scattering

    exposed attributes:

    - multiple_scattering_parameters
    """
    def __init__(self,framestream,ms_parameters,hsrlscope='rs_inv',particlescope='rs_particle',outputscope='rs_multiple_scattering',*args,**kwargs):
        super(dpl_multiple_scattering,self).__init__(framestream)
        self.framestream=framestream
        import cooperative.core.multiple_scattering as ms
        self.particlescope=particlescope
        self.hsrlscope=hsrlscope
        self.outputscope=outputscope
        try:
            pp=self.framestream.spheroid_particle_parameters
        except:
            pp=None
        #FIXME additional parameters. if calvals are needed, extract it from framestream
        self.callable=ms.multiple_scattering(ms_parameters,pp,*args,**kwargs)
        #print framestream.provides
 
    def triage(self,hsrl,particle,rs):
        timeaxis=None
        if (timeaxis==None or timeaxis.size==0) and hsrl!=None and hasattr(hsrl,hsrl._timevarname):
            timeaxis=getattr(hsrl,hsrl._timevarname)
        if (timeaxis==None or timeaxis.size==0) and particle!=None and hasattr(particle,particle._timevarname):
            timeaxis=getattr(particle,particle._timevarname)
        if timeaxis==None:
            return None
        return hau.Time_Z_Group(timeaxis.copy(),timevarname='times',altname='heights')

    
    @property
    def multiple_scattering_instrument(self):
        if self.framestream is not None:
            if hasattr(self.framestream,'spheroid_instrument'):
                ret=self.framestream.spheroid_instrument
            else:
                ret=self.framestream.hsrl_instrument
        else:
            if self.particlescope is not None and hasattr(self.particlescope,'spheroid_instrument'):
                ret=self.particlescope.spheroid_instrument
            else:
                ret=self.hsrlscope.hsrl_instrument
        return ret+' Multiple Scattering'

    @property
    def multiple_scattering_parameters(self):
        return self.callable.multiple_scatter_parameters

    def process(self):
        """ main dpl generator function
        """
        #lred=self.lred
        for rs in self.framestream:
            rs_inv=getattr(rs,self.hsrlscope) if hasattr(rs,self.hsrlscope) else None
            rs_particle=None if ( self.particlescope==None or not hasattr(rs,self.particlescope) ) else getattr(rs,self.particlescope)
            if rs_inv==None or not hasattr(rs_inv,'beta_a_backscat'):
                print 'MISSING INVERTED DATA'
                rs_multiple_scattering=self.triage(rs_inv,rs_particle,rs)
            else:
                #hsrl_sounding=self.framestream.hsrl_sounding

                # BEGIN WORK. extra parameters needed here should be added at __init__, like a json for particle parameters
                # rs should have rs.rs_inv, rs.rs_mmcr, and rs.rs_mean already on the same grid.
                # functions outside of this source can be called as needed, conventient, or for any reason.
                # parameters to those functions should be set up in __init__, discovered or extrapolated from the rs structure, or retrieved from the state of the self.framestream
                #   (instrument calibration tables, hsrl_constants,hsrl_sounding,hsrl_calvals,hsrl_Cxx, and hsrl_instrument are available here. e.g. "self.framestream.hsrl_calvals")
                # If neither are possible (runtime, depending on a global state not in the framestream), better design on that retrieval is NEEDED
                # additionally, if something OUTSIDE this object can get exposed attributes of the framestream's objects. e.g. "particle_parameters" is exposed on line 10 of this file
                #  this means that this DPL object can be used to get particle_parameters. exposure means that if a dpl_hsrl_radar object is anywhere in a DPL stream object,
                #  that object.particle_parameters will get the parameters from here.

                print 'CALLED Multiple Scattering PROCESS'

                
                print
                rs_multiple_scattering=self.callable(rs_inv,rs_particle,self.framestream.hsrl_constants_first)
                # END WORK
            if rs_multiple_scattering==None:
                pass
            elif self.outputscope!=None:
                setattr(rs,self.outputscope,rs_multiple_scattering)
                yield rs
            else:
                yield rs_multiple_scattering

