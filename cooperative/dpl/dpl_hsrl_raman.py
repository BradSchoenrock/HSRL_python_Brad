
import dplkit.role.filter
import dplkit.role.blender
import dplkit.role.decorator
import copy
from datetime import datetime,timedelta
import lg_base.core.array_utils as hau
import numpy as np
import os
import json
import traceback

from collections import namedtuple

#@dplkit.role.decorator.exposes_attrs_in_chain(['spheroid_particle_parameters'])
@dplkit.role.decorator.autoprovidenested(nestedclasses=[hau.Time_Z_Group,hau.rs_xfer,dict])
class dpl_hsrl_raman_profile(dplkit.role.blender.aBlender):
    """ HSRL/RamanMerge photon counts cooperative profile product test. requires both scopes be on the same resolution/gridsize

    adds ramanmerge_hsrl_test to frame

    :param framestream: input iterable framestream
    :param particle_parameters: None for default (as determined by calvals), for filename or dictionary of particle parameters/crystal distribution
    :param ramanscope: nested subframe source for Raman Merge (photon count data)
    :param hsrlscope: nested subframe source for HSRL inverted data
    :param outscope: change rs_particle to another subframe

    exposed attributes:

    - spheroid_particle_parameters
    """
    def __init__(self,framestream,hsrl_ramanmerge_parameters=None,ramanscope=None,hsrlscope=None\
                 ,outputscope=None,*args,**kwargs):
        super(dpl_hsrl_raman_profile,self).__init__([framestream] if framestream!=None else [hsrlscope,ramanscope])
        self.framestream=framestream
        self.ramanscope=ramanscope
        self.hsrlscope=hsrlscope
        self.load_parameters(hsrl_ramanmerge_parameters)
        #import cooperative.core.lidar_radar_eff_diameter as lred
        import cooperative.core.hsrl_raman_test as pp
        #self.lred=lred
        self.pp=pp
        self.outputscope=outputscope or 'ramanmerge_hsrl_test'
        #print framestream.provides
        
    def load_parameters(self,parms):
        if parms==None:
            self._parameters=None
            return 
        systemOnly=(parms==None)
        if parms==None:
            parms='spheroid_particle_parameters_default.json'
        if isinstance(parms,basestring):
            from lg_base.core.locate_file import locate_file
            parms=json.load(open(locate_file(parms,systemOnly=systemOnly),'r'))
        if not isinstance(parms,dict):
            raise RuntimeError('Particle Parameters need to be a json filename or a dictionary')
        self._parameters=parms

    @property
    def xxxxxspheroid_particle_parameters(self):
        return self._parameters

    def triage(self,hsrl,raman):
        timeaxis=None
        if (timeaxis==None or timeaxis.size==0) and hsrl!=None and hasattr(hsrl,hsrl._timevarname):
            timeaxis=getattr(hsrl,hsrl._timevarname)
        if (timeaxis==None or timeaxis.size==0) and raman!=None and hasattr(raman,raman._timevarname):
            timeaxis=getattr(raman,raman._timevarname)
        if timeaxis!=None:
            return hau.Time_Z_Group(timeaxis.copy(),timevarname='times',altname='altitudes',can_append=False)
        return None

    def combine(self):
        """ main dpl generator function
        """
        parameters=self._parameters
        #lred=self.lred
        pp=self.pp
        if self.framestream!=None:
          for rs in self.framestream:
            rs_mean=getattr(rs,self.hsrlscope) if hasattr(rs,self.hsrlscope) else None
            rs_merge=getattr(rs,self.ramanscope) if hasattr(rs,self.ramanscope) else None
            if rs_mean==None or rs_merge==None:
                print 'MISSING PARTICLE DATA'
                x=self.triage(rs_mean,rs_merge)
                if x!=None:
                    setattr(rs,self.outputscope,x)
                yield rs
                continue
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

            print 'CALLED HSRL RAMANMERGE PROCESS'

            if not hasattr(rs_mean,'molecular_counts') or not hasattr(rs_merge,'nitrogen_counts'):
                x=self.triage(rs_mean,rs_merge)
                if x!=None:
                    setattr(rs,self.outputscope,x)
                yield rs
                continue
            if not all([(rs_mean.molecular_counts.shape[x]==rs_merge.nitrogen_counts.shape[x]) for x in range(len(rs_merge.nitrogen_counts.shape))]):
                print 'Size Mismatch of variables!'
                print rs_mean.times
                print rs_merge.times
                raise RuntimeError('Array Size Mismatch in raman merge/hsrl processing. Check the log or contact your administrator with information on how to recreate this crash')

            #rs_particle=pp.process_spheroid_particle(rs_inv,rs_radar,particle_parameters\
            #        ,lambda_radar=self.framestream.radarLambda,entire_frame=rs\
            #        ,sounding=hsrl_sounding,bs_ratio_to_dmode=self.bs_ratio_to_dmode\
            #        ,size_dist=self.size_dist)
           
            mergemean=pp.process_hsrl_raman_profile(rs_mean,rs_merge,parameters,entire_frame=rs)
            # END WORK
            setattr(rs,self.outputscope,mergemean)
            yield rs
        else:
            meansrc=iter(self.hsrlscope)
            mergesrc=iter(self.ramanscope)
            while True:
                try:
                    if meansrc!=None:
                        rs_mean=meansrc.next()
                except StopIteration:
                    rs_mean=None
                    meansrc=None
                try:
                    if mergesrc!=None:
                        rs_merge=mergesrc.next()
                except StopIteration:
                    rs_merge=None
                    mergesrc=None
                if meansrc==None and mergesrc==None:
                    break #done!
                if rs_mean==None or rs_merge==None:
                    print 'MISSING HSRLRAMANAMERGE DATA'
                    x=self.triage(rs_mean,rs_merge)
                    if x!=None:
                        yield x
                    continue
                #hsrl_sounding=self.hsrlscope.hsrl_sounding

                # BEGIN WORK. extra parameters needed here should be added at __init__, like a json for particle parameters
                # rs should have rs.rs_inv, rs.rs_mmcr, and rs.rs_mean already on the same grid.
                # functions outside of this source can be called as needed, conventient, or for any reason.
                # parameters to those functions should be set up in __init__, discovered or extrapolated from the rs structure, or retrieved from the state of the self.framestream
                #   (instrument calibration tables, hsrl_constants,hsrl_sounding,hsrl_calvals,hsrl_Cxx, and hsrl_instrument are available here. e.g. "self.framestream.hsrl_calvals")
                # If neither are possible (runtime, depending on a global state not in the framestream), better design on that retrieval is NEEDED
                # additionally, if something OUTSIDE this object can get exposed attributes of the framestream's objects. e.g. "particle_parameters" is exposed on line 10 of this file
                #  this means that this DPL object can be used to get particle_parameters. exposure means that if a dpl_hsrl_radar object is anywhere in a DPL stream object,
                #  that object.particle_parameters will get the parameters from here.

                print 'CALLED SPHEROID PARTICLE PROCESS'

                if not hasattr(rs_merge,'nitrogen_counts') or not hasattr(rs_mean,'molecular_counts'):
                    x=self.triage(rs_mean,rs_merge)
                    if x!=None:
                        yield x
                    continue
                if not all([(rs_mean.molecular_counts.shape[x]==rs_merge.nitrogen_counts.shape[x]) for x in range(len(rs_merge.nitrogen_counts.shape))]):
                    print 'Size Mismatch of variables!'
                
                    raise RuntimeError('Array Size Mismatch in hsrl/raman merge processing. Check the log or contact your administrator with information on how to recreate this crash')

                mergemean=pp.process_hsrl_raman_profile(rs_mean,rs_merge,parameters,entire_frame=None)#,sounding=hsrl_sounding)
              
                # END WORK
                yield mergemean

#@dplkit.role.decorator.exposes_attrs_in_chain(['spheroid_particle_parameters'])
@dplkit.role.decorator.autoprovidenested(nestedclasses=[hau.Time_Z_Group,hau.rs_xfer,dict])
class dpl_hsrl_ramanmerge(dplkit.role.blender.aBlender):
    """ HSRL/RamanMerge photon counts cooperative product test. requires both scopes be on the same resolution/gridsize

    adds ramanmerge_hsrl_test to frame

    :param framestream: input iterable framestream
    :param particle_parameters: None for default (as determined by calvals), for filename or dictionary of particle parameters/crystal distribution
    :param ramanscope: nested subframe source for Raman Merge (photon count data)
    :param hsrlscope: nested subframe source for HSRL inverted data
    :param outscope: change rs_particle to another subframe

    exposed attributes:

    - spheroid_particle_parameters
    """
    def __init__(self,framestream,hsrl_ramanmerge_parameters=None,ramanscope='rlprofmerge',hsrlscope='rs_mean'\
                 ,outputscope=None,*args,**kwargs):
        super(dpl_hsrl_ramanmerge,self).__init__([framestream] if framestream!=None else [hsrlscope,ramanscope])
        self.framestream=framestream
        self.ramanscope=ramanscope
        self.hsrlscope=hsrlscope
        self.load_parameters(hsrl_ramanmerge_parameters)
        #import cooperative.core.lidar_radar_eff_diameter as lred
        import cooperative.core.hsrl_raman_test as pp
        #self.lred=lred
        self.pp=pp
        self.outputscope=outputscope or 'ramanmerge_hsrl_test'
        #print framestream.provides
        
    def load_parameters(self,parms):
        if parms==None:
            self._parameters=None
            return 
        systemOnly=(parms==None)
        if parms==None:
            parms='spheroid_particle_parameters_default.json'
        if isinstance(parms,basestring):
            from lg_base.core.locate_file import locate_file
            parms=json.load(open(locate_file(parms,systemOnly=systemOnly),'r'))
        if not isinstance(parms,dict):
            raise RuntimeError('Particle Parameters need to be a json filename or a dictionary')
        self._parameters=parms

    @property
    def xxxxxspheroid_particle_parameters(self):
        return self._parameters

    def triage(self,hsrl,raman):
        timeaxis=None
        if (timeaxis==None or timeaxis.size==0) and hsrl!=None and hasattr(hsrl,hsrl._timevarname):
            timeaxis=getattr(hsrl,hsrl._timevarname)
        if (timeaxis==None or timeaxis.size==0) and raman!=None and hasattr(raman,raman._timevarname):
            timeaxis=getattr(raman,raman._timevarname)
        if timeaxis!=None:
            return hau.Time_Z_Group(timeaxis.copy(),timevarname='times',altname='altitudes')
        return None

    def combine(self):
        """ main dpl generator function
        """
        parameters=self._parameters
        #lred=self.lred
        pp=self.pp
        if self.framestream!=None:
          for rs in self.framestream:
            rs_mean=getattr(rs,self.hsrlscope) if hasattr(rs,self.hsrlscope) else None
            rs_merge=getattr(rs,self.ramanscope) if hasattr(rs,self.ramanscope) else None
            if rs_mean==None or rs_merge==None:
                print 'MISSING PARTICLE DATA'
                x=self.triage(rs_mean,rs_merge)
                if x!=None:
                    setattr(rs,self.outputscope,x)
                yield rs
                continue
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

            print 'CALLED HSRL RAMANMERGE PROCESS'

            if not hasattr(rs_mean,'molecular_counts') or not hasattr(rs_merge,'nitrogen_counts'):
                x=self.triage(rs_mean,rs_merge)
                if x!=None:
                    setattr(rs,self.outputscope,x)
                yield rs
                continue
            if not all([(rs_mean.molecular_counts.shape[x]==rs_merge.nitrogen_counts.shape[x]) for x in range(len(rs_merge.nitrogen_counts.shape))]):
                print 'Size Mismatch of variables!'
                print rs_mean.times
                print rs_merge.times
                raise RuntimeError('Array Size Mismatch in raman merge/hsrl processing. Check the log or contact your administrator with information on how to recreate this crash')

            #rs_particle=pp.process_spheroid_particle(rs_inv,rs_radar,particle_parameters\
            #        ,lambda_radar=self.framestream.radarLambda,entire_frame=rs\
            #        ,sounding=hsrl_sounding,bs_ratio_to_dmode=self.bs_ratio_to_dmode\
            #        ,size_dist=self.size_dist)
            mergemean=pp.process_hsrl_ramanmerge(rs_mean,rs_merge,parameters,entire_frame=rs)
            # END WORK
            setattr(rs,self.outputscope,mergemean)
            yield rs
        else:
            meansrc=iter(self.hsrlscope)
            mergesrc=iter(self.ramanscope)
            while True:
                try:
                    if meansrc!=None:
                        rs_mean=meansrc.next()
                except StopIteration:
                    rs_mean=None
                    meansrc=None
                try:
                    if mergesrc!=None:
                        rs_merge=mergesrc.next()
                except StopIteration:
                    rs_merge=None
                    mergesrc=None
                if meansrc==None and mergesrc==None:
                    break #done!
                if rs_mean==None or rs_merge==None:
                    print 'MISSING HSRLRAMANAMERGE DATA'
                    x=self.triage(rs_mean,rs_merge)
                    if x!=None:
                        yield x
                    continue
                #hsrl_sounding=self.hsrlscope.hsrl_sounding

                # BEGIN WORK. extra parameters needed here should be added at __init__, like a json for particle parameters
                # rs should have rs.rs_inv, rs.rs_mmcr, and rs.rs_mean already on the same grid.
                # functions outside of this source can be called as needed, conventient, or for any reason.
                # parameters to those functions should be set up in __init__, discovered or extrapolated from the rs structure, or retrieved from the state of the self.framestream
                #   (instrument calibration tables, hsrl_constants,hsrl_sounding,hsrl_calvals,hsrl_Cxx, and hsrl_instrument are available here. e.g. "self.framestream.hsrl_calvals")
                # If neither are possible (runtime, depending on a global state not in the framestream), better design on that retrieval is NEEDED
                # additionally, if something OUTSIDE this object can get exposed attributes of the framestream's objects. e.g. "particle_parameters" is exposed on line 10 of this file
                #  this means that this DPL object can be used to get particle_parameters. exposure means that if a dpl_hsrl_radar object is anywhere in a DPL stream object,
                #  that object.particle_parameters will get the parameters from here.

                print 'CALLED SPHEROID PARTICLE PROCESS'

                if not hasattr(rs_merge,'nitrogen_counts') or not hasattr(rs_mean,'molecular_counts'):
                    x=self.triage(rs_mean,rs_merge)
                    if x!=None:
                        yield x
                    continue
                if not all([(rs_mean.molecular_counts.shape[x]==rs_merge.nitrogen_counts.shape[x]) for x in range(len(rs_merge.nitrogen_counts.shape))]):
                    print 'Size Mismatch of variables!'
                
                    raise RuntimeError('Array Size Mismatch in hsrl/raman merge processing. Check the log or contact your administrator with information on how to recreate this crash')

                mergemean=pp.process_hsrl_ramanmerge(rs_mean,rs_merge,parameters,entire_frame=None)#,sounding=hsrl_sounding)
              
                # END WORK
                yield mergemean

            
#@dplkit.role.decorator.exposes_attrs_in_chain(['mass_dimension_particle_parameters'])
@dplkit.role.decorator.autoprovidenested(nestedclasses=[hau.Time_Z_Group,hau.rs_xfer,dict])
class dpl_hsrl_raman(dplkit.role.blender.aBlender):
    """ HSRL/Raman Lidar Inverted Cooperative process. requires both scopes be on the same resolution/gridsize

    adds hsrl_raman_test to frame

    :param framestream: input iterable framestream
    :param particle_parameters: None for default (as determined by calvals), for filename or dictionary of particle parameters/crystal distribution
    :param ramanscope: nested subframe source for Raman Inverted data
    :param hsrlscope: nested subframe source for HSRL inverted data

    exposed attributes:

    - mass_dimension_particle_parameters
    """
    def __init__(self,framestream,parameters=None,ramanscope='rlprofdep',hsrlscope='rs_inv',outputscope=None,*args,**kwargs):
        super(dpl_hsrl_raman,self).__init__([framestream] if framestream!=None else [hsrlscope,ramanscope])
        self.framestream=framestream
        self.ramanscope=ramanscope
        self.hsrlscope=hsrlscope
        self.load_parameters(parameters)
        #import cooperative.core.lidar_radar_eff_diameter as lred
        import cooperative.core.hsrl_raman_test as pp
        #self.lred=lred
        self.pp=pp
        self.outputscope=outputscope or 'hsrl_raman_test'
        #print framestream.provides

    def triage(self,hsrl,raman):
        timeaxis=None
        if (timeaxis==None or timeaxis.size==0) and hsrl!=None and hasattr(hsrl,hsrl._timevarname):
            timeaxis=getattr(hsrl,hsrl._timevarname)
        if (timeaxis==None or timeaxis.size==0) and raman!=None and hasattr(raman,raman._timevarname):
            timeaxis=getattr(raman,raman._timevarname)
        if timeaxis!=None:
            return hau.Time_Z_Group(timeaxis.copy(),timevarname='times',altname='altitudes')
        return None
 
    def load_parameters(self,parms):
        if parms==None:
            self._parameters=None
            return
        if parms==None:
            try:
                cv=self.framestream.hsrl_constants_first if self.framestream!=None else self.hsrlscope.hsrl_constants_first
                if 'default_particle_parameters' in cv:
                    parms=cv['default_particle_parameters']
            except AttributeError:
                print 'WARNING: HSRL is missing. can\'t find hsrl_constants_first'
                #raise RuntimeError("Framestream doesn't have hsrl_constants. HSRL is missing?")
        systemOnly=(parms==None)
        if parms==None:
            parms='particle_parameters_default.json'
        if isinstance(parms,basestring):
            from lg_base.core.locate_file import locate_file
            parms=json.load(open(locate_file(parms,systemOnly=systemOnly),'r'))
        if not isinstance(parms,dict):
            raise RuntimeError('Particle Parameters need to be a json filename or a dictionary')
        
        self._particle_parameters=parms

    
    @property
    def xxxxmass_dimension_instrument(self):
        if self.framestream is not None:
            return self.framestream.hsrl_instrument+'-'+self.framestream.radarType
        return  self.hsrlscope.hsrl_instrument+'-'+self.ramanscope.radarType

    @property
    def xxxxmass_dimension_particle_parameters(self):
        return self._particle_parameters

    def combine(self):
        """ main dpl generator function
        """
        parameters=self._parameters
        #lred=self.lred
        pp=self.pp
        if self.framestream!=None:
          for rs in self.framestream:
            if not hasattr(rs,self.hsrlscope) or not hasattr(rs,self.ramanscope):
                print 'MISSING hsrlraman DATA'
                yield rs
                continue
            rs_inv=getattr(rs,self.hsrlscope)
            rs_raman=getattr(rs,self.ramanscope)

            # BEGIN WORK. extra parameters needed here should be added at __init__, like a json for particle parameters
            # rs should have rs.rs_inv, rs.rs_mmcr, and rs.rs_mean already on the same grid.
            # functions outside of this source can be called as needed, conventient, or for any reason.
            # parameters to those functions should be set up in __init__, discovered or extrapolated from the rs structure, or retrieved from the state of the self.framestream
            #   (instrument calibration tables, hsrl_constants,hsrl_sounding,hsrl_calvals,hsrl_Cxx, and hsrl_instrument are available here. e.g. "self.framestream.hsrl_calvals")
            # If neither are possible (runtime, depending on a global state not in the framestream), better design on that retrieval is NEEDED
            # additionally, if something OUTSIDE this object can get exposed attributes of the framestream's objects. e.g. "particle_parameters" is exposed on line 10 of this file
            #  this means that this DPL object can be used to get particle_parameters. exposure means that if a dpl_hsrl_radar object is anywhere in a DPL stream object,
            #  that object.particle_parameters will get the parameters from here.

            print 'CALLED hsrlraman PROCESS'

            if not hasattr(rs_raman,'backscatter') or not hasattr(rs_inv,'beta_a_backscat'):
                yield rs
                continue
            if not all([(rs_inv.beta_a_backscat.shape[x]==rs_raman.backscatter.shape[x]) for x in range(len(rs_raman.backscatter.shape))]):
                print 'Size Mismatch of variables!'
                print rs_inv.times
                print rs_raman.times
                raise RuntimeError('Array Size Mismatch in raman/hsrl processing. Check the log or contact your administrator with information on how to recreate this crash')


            hsrlraman=pp.process_hsrl_raman(rs_inv,rs_raman,parameters)#,lambda_radar=self.framestream.radarLambda,entire_frame=rs)

            # END WORK
            setattr(rs,self.outputscope,hsrlraman)
            yield rs
        else:
            invsrc=iter(self.hsrlscope)
            ramsrc=iter(self.ramanscope)
            while True:
                try:
                    if invsrc!=None:
                        rs_inv=invsrc.next()
                except StopIteration:
                    rs_inv=None
                    invsrc=None
                try:
                    if ramsrc!=None:
                        rs_raman=ramsrc.next()
                except StopIteration:
                    rs_raman=None
                    ramsrc=None
                if invsrc==None and ramsrc==None:
                    break #done!
                if rs_inv==None or rs_raman==None:
                    print 'MISSING HSRLRAMAN DATA'
                    x=self.triage(rs_inv,rs_raman)
                    if x!=None:
                        yield x
                    continue
                #hsrl_sounding=self.hsrlscope.hsrl_sounding

                # BEGIN WORK. extra parameters needed here should be added at __init__, like a json for particle parameters
                # rs should have rs.rs_inv, rs.rs_mmcr, and rs.rs_mean already on the same grid.
                # functions outside of this source can be called as needed, conventient, or for any reason.
                # parameters to those functions should be set up in __init__, discovered or extrapolated from the rs structure, or retrieved from the state of the self.framestream
                #   (instrument calibration tables, hsrl_constants,hsrl_sounding,hsrl_calvals,hsrl_Cxx, and hsrl_instrument are available here. e.g. "self.framestream.hsrl_calvals")
                # If neither are possible (runtime, depending on a global state not in the framestream), better design on that retrieval is NEEDED
                # additionally, if something OUTSIDE this object can get exposed attributes of the framestream's objects. e.g. "particle_parameters" is exposed on line 10 of this file
                #  this means that this DPL object can be used to get particle_parameters. exposure means that if a dpl_hsrl_radar object is anywhere in a DPL stream object,
                #  that object.particle_parameters will get the parameters from here.

                print 'CALLED HSRLRAMAN PROCESS'

                print 'Inverted HSRL:',rs_inv
                print 'Raman        :',rs_raman
                print '((((((((((((((((((((((((((((((((((((((((((((9'
                if not hasattr(rs_raman,'backscatter') or not hasattr(rs_inv,'beta_a_backscat'):
                    x=self.triage(rs_inv,rs_raman)
                    if x!=None:
                        yield x
                    continue
                if not all([(rs_inv.beta_a_backscat.shape[x]==rs_raman.backscatter.shape[x]) for x in range(len(rs_raman.backscatter.shape))]):
                    print 'Size Mismatch of variables!'
                    print rs_inv.times
                    print rs_raman.times
                    raise RuntimeError('Array Size Mismatch in hsrl/raman processing. Check the log or contact your administrator with information on how to recreate this crash')


                hsrlraman=pp.process_hsrl_raman(rs_inv,rs_raman,parameters)#,lambda_radar=self.radarscope.radarLambda,entire_frame=None)

                # END WORK
                yield hsrlraman
