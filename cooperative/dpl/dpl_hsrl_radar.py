
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

@dplkit.role.decorator.exposes_attrs_in_chain(['allradar_hsrl_parameters'])
@dplkit.role.decorator.autoprovidenested(nestedclasses=[hau.Time_Z_Group,hau.rs_xfer,dict])
class dpl_allradar_hsrl_cooperative(dplkit.role.blender.aBlender):
    """ HSRL/Radar Particle information filter using Spheroid Technique. requires both scopes be on the same resolution/gridsize

    adds /p outputscope to frame

    :param framestream: input iterable framestream
    :param particle_parameters: None for default (as determined by calvals), for filename or dictionary of particle parameters/crystal distribution
    :param radarscope: nested subframe source for Radar
    :param hsrlscope: nested subframe source for HSRL inverted data
    :param outscope: change rs_particle to another subframe

    exposed attributes:

    - spheroid_particle_parameters
    """
    def __init__(self,framestream,allradar_hsrl_parameters=None,scopes={},outputscope=None,*args,**kwargs):
        super(dpl_allradar_hsrl_cooperative,self).__init__([framestream] if framestream is not None else [hsrlscope,radarscope])
        self.framestream=framestream
        self.scopes=scopes
        self.load_parameters(allradar_hsrl_parameters)
        #import cooperative.core.lidar_radar_eff_diameter as lred
        import cooperative.hsrl_radar.radar_hsrl_processing as rhp
        #self.lred=lred
        self.rhp=rhp

        self.outputscope=outputscope or 'rs_radar_hsrl'
        #print framestream.provides
        
    def load_parameters(self,parms):
        if parms==None:
            try:
                cv=self.framestream.hsrl_constants_first if self.framestream!=None else self.hsrlscope.hsrl_constants_first
                if 'default_allradar_hsrl_parameters' in cv:
                    parms=cv['default_allradar_hsrl_parameters']
            except AttributeError:
                print 'WARNING: HSRL is missing. can\'t find hsrl_constants_first'
                #raise RuntimeError("Framestream doesn't have hsrl_constants. HSRL is missing?")
        systemOnly=(parms==None)
        if parms==None:
            parms='allradar_hsrl_parameters_default.json'
        if isinstance(parms,basestring):
            from lg_base.core.locate_file import locate_file
            parms=json.load(open(locate_file(parms,systemOnly=systemOnly),'r'))
        if not isinstance(parms,dict):
            raise RuntimeError('ALLRADAR_HSRL Parameters need to be a json filename or a dictionary')
        self._allradar_hsrl_parameters=parms
    
    @property
    def allradar_hsrl_instrument(self):
        return  '-'.join([v for k,v in self.scopes.keys()])

    @property
    def allradar_hsrl_parameters(self):
        return self._allradar_hsrl_parameters

    def triage(self,*fr):
        timeaxis=None
        for fram in fr:
            if (timeaxis==None or timeaxis.size==0) and fram!=None and hasattr(fram,fram._timevarname):
                timeaxis=getattr(fram,fram._timevarname)
        if timeaxis is None:
            return None
        return hau.Time_Z_Group(timeaxis.copy(),timevarname='times',altname='heights')

    def combine(self):
        """ main dpl generator function
        """
        allradar_hsrl_parameters=self._allradar_hsrl_parameters
        #lred=self.lred
        rhp=self.rhp
        if self.framestream!=None:
          for _rs in self.framestream:
            args=[allradar_hsrl_parameters]
            kwargs={}
            shapes={}
            willTriage=False
            rs=copy.copy(_rs)
            for k,v in self.scopes.items():
                kwargs[k]=getattr(rs,v) if hasattr(rs,v) else None
                if kwargs[k] is None:
                    print 'MISSING COOPERATIVE HSRL/RADAR DATA:',v
                    continue
                if 'rs_inv' in v or 'rs_mean' in v:
                    shapes[v]=kwargs[k].beta_a_backscat.shape if hasattr(kwargs[k],'beta_a_backscat') and hasattr(kwargs[k].beta_a_backscat,'shape') else None
                else:
                    shapes[v]=kwargs[k].Backscatter.shape if hasattr(kwargs[k],'Backscatter') and hasattr(kwargs[k].Backscatter,'shape') else None

            if not willTriage:
                startshape=None
                firstshape=None
                for k,v in shapes.items():
                    if v is None:
                        #willTriage=True
                        pass
                    elif startshape is None:
                        startshape=v
                        firstshape=k
                    elif len(v)!=len(startshape):
                        print self.scopes[k],firstshape,"has different size!",v,startshape
                        #willTriage=True
                    else:
                        for i in range(len(v)):
                            if startshape[i]!=v[i]:
                                print self.scopes[k],firstshape,'dimension',i,'differs',v,startshape
                                willTriage=True

            if willTriage:
                x=self.triage(*[v for k,v in kwargs.items()])
                if x!=None:
                    setattr(rs,self.outputscope,x)
                    for k,v in kwargs.items():
                        if v is None:
                            setattr(rs,self.scopes[k],copy.deepcopy(x))
            else:
                phrf=rhp.process_hsrl_radar(*args,**kwargs)
                # END WORK
                setattr(rs,self.outputscope,phrf)
            yield rs
        else:
            invsrc=iter(self.hsrlscope)
            radsrc=iter(self.radarscope)
            while True:
                try:
                    if invsrc!=None:
                        rs_inv=invsrc.next()
                except StopIteration:
                    rs_inv=None
                    invsrc=None
                try:
                    if radsrc!=None:
                        rs_radar=radsrc.next()
                except StopIteration:
                    rs_radar=None
                    radsrc=None
                if invsrc==None and radsrc==None:
                    break #done!
                if rs_inv==None or rs_radar==None:
                    print 'MISSING PARTICLE DATA'
                    x=self.triage(rs_inv,rs_radar)
                    if x!=None:
                        yield x
                    continue
                hsrl_sounding=self.hsrlscope.hsrl_sounding

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

                if not hasattr(rs_radar,'Backscatter') or not hasattr(rs_inv,'beta_a_backscat'):
                    x=self.triage(rs_inv,rs_radar)
                    if x!=None:
                        yield x
                    continue
                if not all([(rs_inv.beta_a_backscat.shape[x]==rs_radar.Backscatter.shape[x]) for x in range(len(rs_radar.Backscatter.shape))]):
                    print 'Size Mismatch of variables!'
                
                    raise RuntimeError('Array Size Mismatch in Cooperative processing. Check the log or contact your administrator with information on how to recreate this crash')

                rs_particle=pp.process_spheroid_particle(rs_inv,rs_radar,particle_parameters,lambda_radar=self.radarscope.radarLambda,entire_frame=None,sounding=hsrl_sounding)
              
                # END WORK
                yield rs_particle

  
@dplkit.role.decorator.exposes_attrs_in_chain(['spheroid_particle_parameters'])
@dplkit.role.decorator.autoprovidenested(nestedclasses=[hau.Time_Z_Group,hau.rs_xfer,dict])
class dpl_spheroid_particle(dplkit.role.blender.aBlender):
    """ HSRL/Radar Particle information filter using Spheroid Technique. requires both scopes be on the same resolution/gridsize

    adds rs_particle to frame

    :param framestream: input iterable framestream
    :param particle_parameters: None for default (as determined by calvals), for filename or dictionary of particle parameters/crystal distribution
    :param radarscope: nested subframe source for Radar
    :param hsrlscope: nested subframe source for HSRL inverted data
    :param outscope: change rs_particle to another subframe

    exposed attributes:

    - spheroid_particle_parameters
    """
    def __init__(self,framestream,particle_parameters=None,radarscope='rs_mmcr'\
                 ,hsrlscope='rs_inv',outputscope=None,*args,**kwargs):
        super(dpl_spheroid_particle,self).__init__([framestream] if framestream is not None else [hsrlscope,radarscope])
        self.framestream=framestream
        self.radarscope=radarscope
        self.hsrlscope=hsrlscope
        self.load_particle_parameters(particle_parameters)
        #import cooperative.core.lidar_radar_eff_diameter as lred
        import cooperative.core.spheroid_particle_processing as pp
        #self.lred=lred
        self.pp=pp


        """ 
        #compute tables of mode diameter vs radar_backscatter/lidar_backscatter
        import cooperative.core.radar_lidar_backscat_ratio as rlr
        self.backscat_ratio = rlr.radar_lidar_backscat_ratio(self.framestream.radarLambda,self._particle_parameters)
        """

        #define size distribution and methods for deff_prime, deff etc
        import cooperative.core.spheroid_utilities as su
        self.size_dist = su.dstar_table(self._particle_parameters,self.framestream.radarLambda)
        """
        #compute table converting d_eff_prime to d_mode using Mie theory, table index = particle diameter in microns
        import cooperative.core.bs_ratio_to_dmode as bsr
        #self.bs_ratio_to_dmode = bsr.bs_ratio_to_dmode(self.framestream.radarLambda,self._particle_parameters)
        self.bs_ratio_to_dmode = bsr.bs_ratio_to_dmode(self.framestream.radarLambda,self.size_dist)
        """
        
        #generates a correction to the assumed backscatter phase function based on first estimate of size
        #produces corrections needed only for drizzle and rain
        import cooperative.core.lidar_p180_water as lpw
        print 'particle_parameters'

        """
        print self._particle_parameters['p180_water']
        #need to get lidar wavelength from calvals to replace fixed 532 nm value
        self.lidar_p180_water=lpw.lidar_p180_water(532e-9,self._particle_parameters)
   
        #provides non-Rayleigh correction to the radar backscatter phase function
        import cooperative.core.radar_p180_water as rpw
        self.radar_p180_water=rpw.radar_p180_water(self.framestream.radarLambda,self._particle_parameters)
        """

        self.outputscope=outputscope or 'rs_particle'
        #print framestream.provides
        
    def load_particle_parameters(self,parms):
        if parms==None:
            try:
                cv=self.framestream.hsrl_constants_first if self.framestream!=None else self.hsrlscope.hsrl_constants_first
                if 'default_spheroid_particle_parameters' in cv:
                    parms=cv['default_spheroid_particle_parameters']
            except AttributeError:
                print 'WARNING: HSRL is missing. can\'t find hsrl_constants_first'
                #raise RuntimeError("Framestream doesn't have hsrl_constants. HSRL is missing?")
        systemOnly=(parms==None)
        if parms==None:
            parms='spheroid_particle_parameters_default.json'
        if isinstance(parms,basestring):
            from lg_base.core.locate_file import locate_file
            parms=json.load(open(locate_file(parms,systemOnly=systemOnly),'r'))
        if not isinstance(parms,dict):
            raise RuntimeError('Particle Parameters need to be a json filename or a dictionary')
        self._particle_parameters=parms
    
    @property
    def spheroid_instrument(self):
        if self.framestream is not None:
            return self.framestream.hsrl_instrument+'-'+self.framestream.radarType
        return  self.hsrlscope.hsrl_instrument+'-'+self.radarscope.radarType

    @property
    def spheroid_particle_parameters(self):
        return self._particle_parameters

    def triage(self,hsrl,radar):
        timeaxis=None
        if (timeaxis==None or timeaxis.size==0) and hsrl!=None and hasattr(hsrl,hsrl._timevarname):
            timeaxis=getattr(hsrl,hsrl._timevarname)
        if (timeaxis==None or timeaxis.size==0) and radar!=None and hasattr(radar,radar._timevarname):
            timeaxis=getattr(radar,radar._timevarname)
        if timeaxis!=None:
            return hau.Time_Z_Group(timeaxis.copy(),timevarname='times',altname='heights')
        return None

    def combine(self):
        """ main dpl generator function
        """
        particle_parameters=self._particle_parameters
        #lred=self.lred
        pp=self.pp
        if self.framestream!=None:
          for rs in self.framestream:
            rs_inv=getattr(rs,self.hsrlscope) if hasattr(rs,self.hsrlscope) else None
            rs_radar=getattr(rs,self.radarscope) if hasattr(rs,self.radarscope) else None
            if rs_inv is None or rs_radar is None:
                print 'MISSING PARTICLE DATA','inv_fail' if rs_inv is None else 'inv_ok','radar_fail' if rs_radar is None else 'radar_ok'
                x=self.triage(rs_inv,rs_radar)
                if x!=None:
                    setattr(rs,self.outputscope,x)
                    if rs_inv is None and rs_radar is not None:
                        setattr(rs,self.hsrlscope,copy.deepcopy(x))
                    if rs_inv is not None and rs_radar is None:
                        setattr(rs,self.radarscope,copy.deepcopy(x))
                yield rs
                continue
            hsrl_sounding=self.framestream.hsrl_sounding

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

            if not hasattr(rs_radar,'Backscatter') or not hasattr(rs_inv,'beta_a_backscat'):
                x=self.triage(rs_inv,rs_radar)
                if x!=None:
                    setattr(rs,self.outputscope,x)
                yield rs
                continue
            if not all([(rs_inv.beta_a_backscat.shape[x]==rs_radar.Backscatter.shape[x]) for x in range(len(rs_radar.Backscatter.shape))]):
                print 'Size Mismatch of variables!'
                print rs_inv.times
                print rs_radar.times
                raise RuntimeError('Array Size Mismatch in Cooperative processing. Check the log or contact your administrator with information on how to recreate this crash')

            #rs_particle=pp.process_spheroid_particle(rs_inv,rs_radar,particle_parameters\
            #        ,lambda_radar=self.framestream.radarLambda,entire_frame=rs\
            #        ,sounding=hsrl_sounding,bs_ratio_to_dmode=self.bs_ratio_to_dmode\
            #        ,size_dist=self.size_dist)
            rs_particle=pp.process_spheroid_particle(rs_inv,rs_radar,particle_parameters\
                    ,lambda_radar=self.framestream.radarLambda,entire_frame=rs\
                    ,sounding=hsrl_sounding,size_dist=self.size_dist)
            # END WORK
            setattr(rs,self.outputscope,rs_particle)
            yield rs
        else:
            invsrc=iter(self.hsrlscope)
            radsrc=iter(self.radarscope)
            while True:
                try:
                    if invsrc!=None:
                        rs_inv=invsrc.next()
                except StopIteration:
                    rs_inv=None
                    invsrc=None
                try:
                    if radsrc!=None:
                        rs_radar=radsrc.next()
                except StopIteration:
                    rs_radar=None
                    radsrc=None
                if invsrc==None and radsrc==None:
                    break #done!
                if rs_inv==None or rs_radar==None:
                    print 'MISSING PARTICLE DATA'
                    x=self.triage(rs_inv,rs_radar)
                    if x!=None:
                        yield x
                    continue
                hsrl_sounding=self.hsrlscope.hsrl_sounding

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

                if not hasattr(rs_radar,'Backscatter') or not hasattr(rs_inv,'beta_a_backscat'):
                    x=self.triage(rs_inv,rs_radar)
                    if x!=None:
                        yield x
                    continue
                if not all([(rs_inv.beta_a_backscat.shape[x]==rs_radar.Backscatter.shape[x]) for x in range(len(rs_radar.Backscatter.shape))]):
                    print 'Size Mismatch of variables!'
                
                    raise RuntimeError('Array Size Mismatch in Cooperative processing. Check the log or contact your administrator with information on how to recreate this crash')

                rs_particle=pp.process_spheroid_particle(rs_inv,rs_radar,particle_parameters,lambda_radar=self.radarscope.radarLambda,entire_frame=None,sounding=hsrl_sounding)
              
                # END WORK
                yield rs_particle

            
@dplkit.role.decorator.exposes_attrs_in_chain(['mass_dimension_particle_parameters'])
@dplkit.role.decorator.autoprovidenested(nestedclasses=[hau.Time_Z_Group,hau.rs_xfer,dict])
class dpl_mass_dimension_particle(dplkit.role.blender.aBlender):
    """ HSRL/Radar Particle information filter using Mass Dimension Technique. requires both scopes be on the same resolution/gridsize

    adds rs_particle to frame

    :param framestream: input iterable framestream
    :param particle_parameters: None for default (as determined by calvals), for filename or dictionary of particle parameters/crystal distribution
    :param radarscope: nested subframe source for Radar
    :param hsrlscope: nested subframe source for HSRL inverted data

    exposed attributes:

    - mass_dimension_particle_parameters
    """
    def __init__(self,framestream,particle_parameters=None,radarscope='rs_mmcr',hsrlscope='rs_inv',outputscope=None,*args,**kwargs):
        super(dpl_mass_dimension_particle,self).__init__([framestream] if framestream!=None else [hsrlscope,radarscope])
        self.framestream=framestream
        self.radarscope=radarscope
        self.hsrlscope=hsrlscope
        self.load_particle_parameters(particle_parameters)
        #import cooperative.core.lidar_radar_eff_diameter as lred
        import cooperative.core.mass_dimension_particle_processing as pp
        #self.lred=lred
        self.pp=pp
        self.outputscope=outputscope or 'rs_particle'
        #print framestream.provides

    def triage(self,hsrl,radar):
        timeaxis=None
        if (timeaxis==None or timeaxis.size==0) and hsrl!=None and hasattr(hsrl,hsrl._timevarname):
            timeaxis=getattr(hsrl,hsrl._timevarname)
        if (timeaxis==None or timeaxis.size==0) and radar!=None and hasattr(radar,radar._timevarname):
            timeaxis=getattr(radar,radar._timevarname)
        if timeaxis!=None:
            return hau.Time_Z_Group(timeaxis.copy(),timevarname='times',altname='heights')
        return None
 
    def load_particle_parameters(self,parms):
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
    def mass_dimension_instrument(self):
        if self.framestream is not None:
            return self.framestream.hsrl_instrument+'-'+self.framestream.radarType
        return  self.hsrlscope.hsrl_instrument+'-'+self.radarscope.radarType

    @property
    def mass_dimension_particle_parameters(self):
        return self._particle_parameters

    def combine(self):
        """ main dpl generator function
        """
        particle_parameters=self._particle_parameters
        #lred=self.lred
        pp=self.pp
        if self.framestream!=None:
          for rs in self.framestream:
            if not hasattr(rs,self.hsrlscope) or not hasattr(rs,self.radarscope):
                print 'MISSING PARTICLE DATA'
                yield rs
                continue
            rs_inv=getattr(rs,self.hsrlscope)
            rs_radar=getattr(rs,self.radarscope)

            # BEGIN WORK. extra parameters needed here should be added at __init__, like a json for particle parameters
            # rs should have rs.rs_inv, rs.rs_mmcr, and rs.rs_mean already on the same grid.
            # functions outside of this source can be called as needed, conventient, or for any reason.
            # parameters to those functions should be set up in __init__, discovered or extrapolated from the rs structure, or retrieved from the state of the self.framestream
            #   (instrument calibration tables, hsrl_constants,hsrl_sounding,hsrl_calvals,hsrl_Cxx, and hsrl_instrument are available here. e.g. "self.framestream.hsrl_calvals")
            # If neither are possible (runtime, depending on a global state not in the framestream), better design on that retrieval is NEEDED
            # additionally, if something OUTSIDE this object can get exposed attributes of the framestream's objects. e.g. "particle_parameters" is exposed on line 10 of this file
            #  this means that this DPL object can be used to get particle_parameters. exposure means that if a dpl_hsrl_radar object is anywhere in a DPL stream object,
            #  that object.particle_parameters will get the parameters from here.

            print 'CALLED MASS DIMENSION PARTICLE PROCESS'

            if not hasattr(rs_radar,'Backscatter') or not hasattr(rs_inv,'beta_a_backscat'):
                yield rs
                continue
            if not all([(rs_inv.beta_a_backscat.shape[x]==rs_radar.Backscatter.shape[x]) for x in range(len(rs_radar.Backscatter.shape))]):
                print 'Size Mismatch of variables!'
                print rs_inv.times
                print rs_radar.times
                raise RuntimeError('Array Size Mismatch in Cooperative processing. Check the log or contact your administrator with information on how to recreate this crash')


            rs_particle=pp.process_mass_dimension_particle(rs_inv,rs_radar,particle_parameters,lambda_radar=self.framestream.radarLambda,entire_frame=rs)

            # END WORK
            setattr(rs,self.outputscope,rs_particle)
            yield rs
        else:
            invsrc=iter(self.hsrlscope)
            radsrc=iter(self.radarscope)
            while True:
                try:
                    if invsrc!=None:
                        rs_inv=invsrc.next()
                except StopIteration:
                    rs_inv=None
                    invsrc=None
                try:
                    if radsrc!=None:
                        rs_radar=radsrc.next()
                except StopIteration:
                    rs_radar=None
                    radsrc=None
                if invsrc==None and radsrc==None:
                    break #done!
                if rs_inv==None or rs_radar==None:
                    print 'MISSING PARTICLE DATA'
                    x=self.triage(rs_inv,rs_radar)
                    if x!=None:
                        yield x
                    continue
                hsrl_sounding=self.hsrlscope.hsrl_sounding

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

                print 'Inverted HSRL:',rs_inv
                print 'Radar        :',rs_radar
                print '((((((((((((((((((((((((((((((((((((((((((((9'
                if not hasattr(rs_radar,'Backscatter') or not hasattr(rs_inv,'beta_a_backscat'):
                    x=self.triage(rs_inv,rs_radar)
                    if x!=None:
                        yield x
                    continue
                if not all([(rs_inv.beta_a_backscat.shape[x]==rs_radar.Backscatter.shape[x]) for x in range(len(rs_radar.Backscatter.shape))]):
                    print 'Size Mismatch of variables!'
                    print rs_inv.times
                    print rs_radar.times
                    raise RuntimeError('Array Size Mismatch in Cooperative processing. Check the log or contact your administrator with information on how to recreate this crash')


                rs_particle=pp.process_mass_dimension_particle(rs_inv,rs_radar,particle_parameters,lambda_radar=self.radarscope.radarLambda,entire_frame=None)

                # END WORK
                yield rs_particle
