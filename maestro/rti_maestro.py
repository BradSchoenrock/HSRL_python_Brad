from datetime import timedelta,datetime
import matplotlib.pyplot as plt
import numpy as np
import hsrl.calibration.cal_file_generation as cfg
import raman.core.calibration.cal_file_generation as rcfg
import hsrl.data_stream.get_scale_corrections as gsc
import hsrl.data_stream.performance_model as pm
import hsrl.data_stream.hsrl_read_utilities as hru
import subprocess as sp
from collections import OrderedDict
import traceback
import os,sys
import lg_base.core.json_config as jc
import json
import lg_base.core.locate_file as lf
import lg_base.core.open_config as oc
import lg_base.core.read_utilities as ru
import dplkit.role.filter
import dplkit.role.decorator
import logging
import lg_base.formats.json_dictionary_utils as jdu
import lg_dpl_toolbox.filters.altslice as aslice
LOG = logging.getLogger(__name__)

class PassThru(dplkit.role.filter.aFilter):
    def __init__(self,source,name):
        super(PassThru,self).__init__(source)
        self.source=source
        self.name=name

    def process(self):
        import lg_base.core.array_utils as hau
        for f in self.source:
            print 'PassThru',self.name,dir(f)
            for n,v in self.provides.items():
                if n not in dir(f):
                    raise RuntimeError('Missing '+n+' in '+self.name)
                if isinstance(getattr(f,n),hau.T_Array) and len(getattr(f,n).shape)==0:
                    raise RuntimeError('Dimensions of '+n+' are 0 in '+self.name+' type '+repr(type(getattr(f,n))))
            yield f

def getdictionary(obj):
    return obj.get_dict()

def formatted_json_print(D,levels=1,prefix=''):
    ignoreParts=('doc','docs','documentation','parameters','Parameters')
    for k in D.keys():
        if isinstance(D[k],dict):
            continue
        if k in ignoreParts:
            continue
        print prefix,k,'=',D[k]
    if levels>0:
        for k in D.keys():
            if not isinstance(D[k],dict):
                continue
            if k in ignoreParts:
                continue
            print prefix,k,'= {'
            formatted_json_print(D[k],levels=levels-1,prefix=prefix+'    ')
            print prefix,'}'

def formatted_print(D):       
    i=0
    line = '    '
    if isinstance(D,dict) or hasattr(D,'keys'):
        names = sorted(D.keys())
    else:
        names = sorted(vars(D).keys())
    for name in names:
        if not name[0] =='_' :
            if isinstance(D,dict):
                v=D[name]
            elif hasattr(D,name):
                v=getattr(D,name)
            else:
                raise RuntimeError('Structure '+str(D)+' advertises having field '+name+' but isn\'t there')
            if hasattr(v,'shape'):
                shape_text = str(v.shape)
            elif isinstance(v,dict):
                if 'shape' in v:
                    shape_text='array as dictionary shape '+str(v['shape'])
                elif 'fields' in v:
                    shape_text='cal table (as dict) '+(','.join(v['fields']))
                elif 'data' in v and 'shape' in v['data']:
                    shape_text='cal table (as dict) '+str(v['data']['shape'])
                else:
                    shape_text='dictionary value or object'
            elif isinstance(v,hru.calibration_vector_table):
                shape_text = 'cal table '
                if hasattr(v,'fields'):
                    shape_text=shape_text+(','.join(v.fields))
                else:
                    shape_text=shape_text+str(v.data.shape)
            else:
                shape_text = 'value or object'
            i=i+1
            line = line + name + ' ' + shape_text + '\t'
            if i>=2:
                print line.expandtabs(50)
                i = 0
                line = '    '
    if i > 0:
        #print line
        print line.expandtabs(50)       
    return

def deephasattr(f,fields):
    if isinstance(fields,basestring):
        return deephasattr(f,fields.split('.'))
    if f is None:
        return False
    if len(fields)==0:
        return True
    if isinstance(f,dict):
        if not fields[0] in f:
            return False
        return deephasattr(f[fields[0]],fields[1:])
    if not hasattr(f,fields[0]):
        return False
    return deephasattr(getattr(f,fields[0]),fields[1:])

def deepgetattr(f,fields):
    if isinstance(fields,basestring):
        return deepgetattr(f,fields.split('.'))
    if len(fields)==0:
        return f
    if isinstance(f,dict):
        return deepgetattr(f[fields[0]],fields[1:])
    return deepgetattr(getattr(f,fields[0]),fields[1:])

def list_fields(frame,fieldlist):
    for field in fieldlist:

        if deephasattr(frame,field):
            print field
            f=deepgetattr(frame,field)
            if not isinstance(f,dict):
                f=vars(f)
            formatted_print(f)

class Rti(object):

    """ HSRL Plotting Routines

  Examples: ::

    r=Rti('gvhsrl','18-Aug-11 00:00', 2, 0, 15)
    r=Rti('gvhsrl', '18-aug-11 00:00','18-aug-11 2:00' 0 15)
    r=Rti('ahsrl', '18-aug-11 1:00','23:21:05', 0, 15)
    r=Rti('ahsrl','1-jul-12 1:00','2:21',0,15,mol_norm_alt=2.1,display='all_plots.json'
              ,netcdf_defaults=None, process_control=None)

  optional cmd line parameters--defaults are provided in "process_control.json" 
    mol_norm_alt------altitude at which to normalize optical depth(km) 
    process_control------select process options other than "process_control.json"
    netcdf_defaults---selection and translation of netcdf variable names to read
                         will be controlled by contents of xxxx_netcdf_defaults.json
                         where xxxx is the instrument('eg ahsrl, mf2hsrl, bagohsrl ...)
    z_res,t_res --- altitude and time resolution structures
    mass_dimension_particle_parameters --- mass dimension particle parameter json filename
    spheroid_particle_parameters --- spheroid model particle parameter json filename
  r.next()

    """
    def __init__( self, instruments, start_time, plot_length, min_alt 
                ,max_alt, display, mol_norm_alt=None, process_control=None
                ,z_res=None,t_res=None,filterclasses=[],filterparams=[]
                ,fullfilterclasses=[],fullfilterparams=[],radar_process_control=None
                ,output_defaults=None,time_sourcename=None,alt_sourcename=None,process=True
                ,*args,**kwargs):

        """ constructs a new Rti instance

        instruments     = hsrl id string (eg. 'ahsrl','gvhsrl','nshsrl','mf2hsrl') 
                          and/or other instruments in a list or tuple ('mmcr','magkazr','nsakazr','magmwacr','spheroid_particle','mass_dimension_particle','magpars2S1','magpars2S2').
        start_time     = time string, first data to plot (eg.'5-may-11 12:30').
        plot_length    = time period to plot in hours (e.g. 2) or
                              end time of first interval(e.g. '6-may-10 1:30').
        min_alt        = minimum altitude (msl) to display (km).
        max_alt        = maximum altitude (msl) to display (km).
        mol_norm_alt    = normalization altitude for optical depth (km).
                              if not supplied mol_norm_alt is taken from
                              process_control.json                                
        display        = name of  .json file containing display directives
        process_control   = alternate to xxxx_process_control.json file containing
                              processing defaults. A dictionary containing the entries
                              in process_control.json may also be provided directly.
        z_res          = requested altitude resolution
                       = not supplied--->select number of alt bins = number of image y-pixels
                       = {'n_dr':n}--->select number of alt bins = n * instrument bin width
                       = {'manual':dz}--->set altitude resolution to dz meters
        t_res          = requested time resolution
                       = not supplied--->select number of time bins = number of image x-pixels
                       = {'n_dt':n}--->select number of time bins = n * instrument bin width
                       = {'manual':dt}--->set time resolution to dt seconds               
        mass_dimension_particle_parameters --- mass dimension particle parameter json filename, dictionary, or None
        spheroid_particle_parameters --- spheroid model particle parameter json filename, dictionary, or None

        expert parameters:
            filterclasses and filterparams --- in-order lists of classes and initialziation parameters respectively to append to the dpl frame stream before accumulation
            fullfilterclasses and fillfilterparams --- in-order lists of classes and initialziation parameters respectively to append to the dpl frame stream before accumulation
        """
        
        print "*****************hi******************"

        self.calibration_overrides=None
        self.libs=OrderedDict()
        if isinstance(instruments,basestring):
            instruments=[instruments]
        self.instrumentrequests=set(instruments)
        firstbase=None
        claimedkwargs=[]
        for inst in instruments:
            libbase=None
            if inst.endswith('hsrl_raw'):
                from hsrl.dpl.dpl_hsrl import dpl_hsrl
                self.libs['hsrl_raw']=dpl_hsrl(instrument=inst.replace('_raw',''))
                hsrlSource=self.libs['hsrl_raw']
                libbase=self.libs['hsrl_raw'].instrument
            elif inst.endswith('hsrl'):
                from hsrl.dpl.dpl_hsrl import dpl_hsrl
                self.libs['hsrl']=dpl_hsrl(process_control=process_control,instrument=inst)
                hsrlSource=self.libs['hsrl']
                libbase=self.libs['hsrl'].instrument
            elif inst.endswith('hsrl_profile') and len(inst.split('_'))==2:
                from hsrl.dpl.dpl_hsrl import dpl_hsrl
                self.libs['hsrl_profile']=dpl_hsrl(process_control=process_control,instrument=inst.replace('_profile',''))
                hsrlSource=self.libs['hsrl_profile']
                libbase=self.libs['hsrl_profile'].instrument
            elif 'mmcr' in inst or 'kazr' in inst or 'mwacr' in inst:
                from radar.dpl.dpl_radar import dpl_radar
                tmp=dpl_radar(instrument=inst,process_control=radar_process_control)
                self.libs[tmp.instrument]=tmp
                libbase=self.libs[tmp.instrument].instrumentbase
            elif 'met' in inst:
                from met.dpl.dpl_marinemet import dpl_marinemet
                self.libs['met']=dpl_marinemet(instrument=inst)
                libbase=self.libs['met'].instrumentbase
            elif 'rain' in inst:
                from precip.dpl.dpl_rain import dpl_rain
                self.libs['rain']=dpl_rain(instrument=inst)
                libbase=self.libs['rain'].instrumentbase
            elif 'vdis' in inst:
                from precip.dpl.dpl_vdis import dpl_vdis
                self.libs['vdis']=dpl_vdis(instrument=inst)
                libbase=self.libs['vdis'].instrumentbase
            elif 'pars' in inst:
                from pars.dpl.dpl_pars import dpl_pars
                if not ('pars' in self.libs):
                    self.libs['pars']=OrderedDict()
                self.libs['pars']['rs_'+inst[3:]]=dpl_pars(instrument=inst)
                libbase=self.libs['pars']['rs_'+inst[3:]].instrumentbase
            elif inst.startswith('rlprof'):
                from raman.dpl.raman_dpl import dpl_raman
                if not ('raman' in self.libs):
                    self.libs['raman']=OrderedDict()
                self.libs['raman'][inst]=dpl_raman('bagohsrl',inst)
                libbase=self.libs['raman'][inst].instrumentbase
            elif inst in ('spheroid_particle','mass_dimension_particle','multiple_scattering','raman_hsrl_test','allradar_hsrl_coop',\
                    'ramanmerge_hsrl_test','raman_inv','raman_hsrl_profile'):
                claimedkwargs.append(inst+'_parameters')
                pass # runtime only FIXME this should check for hsrl and radar
            else:
                raise NotImplementedError('Unknown instrument '+inst)
            if libbase is not None:
                if firstbase is None:
                    firstbase=libbase
                if firstbase!=libbase:
                    raise RuntimeError('Incongruence in libraries:',inst,'is congruent to the',libbase,'while all others are',firstbase)
        self.instrument=firstbase

        #if 'hsrl' not in self.libs and 'hsrl_raw' not in self.libs and 'hsrl_profile' not in self.libs:
        #    raise RuntimeError('Need an HSRL')
        import lg_base.graphics.graphics_toolkit as gt
        self.gt=gt
        self.figs=self.gt.figurelist()
        self.artistparams=dict()
        if kwargs.pop('singleSandbox',True):
            self.artistparams['figurecontainer']=self.figs
        self.dropcontent=kwargs.pop('DropFrames',None)
        self.config=None
        self.setdisplay(display)
        self.setconfig(output_defaults or 'rti_default_config.json')#used for optional default printer, format, and path info. older files also contain display info 
        self.z_res=z_res
        self.t_res=t_res
        self.alt_host=alt_sourcename
        self.time_host=time_sourcename
        self.initargs=args#particle_parameters=particle_parameters
        if len(args)>0:
            raise TypeError("Initialization doesn't support additional unnamed parameters")
        self.initkwargs=kwargs
        self.filterclasses=filterclasses
        self.filterparams=filterparams
        self.fullfilterclasses=fullfilterclasses
        self.fullfilterparams=fullfilterparams
        self.searchparms=dict(min_alt_m=min_alt*1000.0,max_alt_m=max_alt*1000.0,mol_norm_alt_m=mol_norm_alt*1000.0 if mol_norm_alt is not None else None)        
        if 'with_profiles' in kwargs:
            self.searchparms['with_profiles']=kwargs.pop('with_profiles')
        #import hsrl.dpl_experimental.dpl_artists as artists
        #self.artists=artists
        for k in kwargs.keys():
            if not k in claimedkwargs:
                print "WARNING: '"+k+"' is an invalid keyword argument for this function for the requested instruments"

        import lg_base.core.read_utilities as hru
        rs_date = hru.convert_date_str(start_time)
        self.start_time = rs_date['datetime']

        # if plot_length entered as end time string

        if isinstance(plot_length, basestring):

            # check for time without date--assume it is same day as start

            if plot_length.find(':') < 3:
                index = start_time.find(' ')
                plot_length = start_time[0:index + 1] + plot_length
            rs_date = hru.convert_date_str(plot_length)
            end_time = rs_date['datetime']
            self.delta_t = end_time - self.start_time
            if self.delta_t <= timedelta(seconds=0):
                print ' '
                print 'ERROR----end time must be later than start time'
                print ' '
                raise RuntimeError, 'ERROR----end time must be later than start time'
        else:
            self.delta_t = timedelta(days=(plot_length) / 24.0)  # plot length in hours converted to timedelta

        if process:
            self.update()
            self.display()
            print 'use process=False to have Rti initialize without processing'
        else:
            print 'RTI Initialized. To process, use r.reprocess(), or r.update() to just run processing without display'

    def reload_cals(self):
        for l,libr in self.libs.items():
            if hasattr(libr,'reload'):
                libr.reload()
                print 'Reloaded librarian object',l

    def next(self):
        """ Plot next time interval on top of last plot"""
        now = datetime.utcnow()
        self.start_time+=self.delta_t
        if now<(self.start_time+self.delta_t):
            self.start_time=now-self.delta_t
        self.update()
        self.display()

    def cnext(self):
        """ Clear data from plots, then process and plot next time interval"""
        now = datetime.utcnow()
        self.start_time+=self.delta_t
        if now<(self.start_time+self.delta_t):
            self.start_time=now-self.delta_t
        self.clear()
        self.update()
        self.display()

    def change_sounding_calvals(self,**kwargs):
        for k,v in self.libs.items():
            if not k.startswith('hsrl'):
                continue
            print 'updating calvals for',k,'to',kwargs
            v.cal.change_sounding_calvals(**kwargs)

    def setdisplay(self,display):
        #[self.disp,dummy]=du.get_display_defaults(display,"new")
        print 'fetching display info from ',display
        self.dispfilename=display
        self.display_defaults_multisource_generator=jdu.jsonListAccumulator(self.dispfilename,majorkey='display_defaults',allow_missing_values=True)
        #self.setconfig(display)

    def setconfig(self,conf=None):
        if conf is not None:
            self.confname=conf
        try:
            tmp=jc.json_config(lf.locate_file(self.confname),'config',allow_missing_values=False)
            self.config=tmp
            print 'using RTI config '+lf.locate_file(self.confname)
        except (IOError,KeyError):#doesn't exist on this scope, or doesn't have the key.
            if self.config is None:
                print 'WARNING : no defaults for printer or netcdf output'

    def replot(self,display=None):
        """ replot(display=display_defaults.json')
            clear figures, reload display defaults file, and replot last data
            display_defaults_file = optional--new display defaults file name
                                  = .json file containing display directives"""

        self.clear()
        self.setdisplay(display or self.dispfilename)#this triggers a reload
        self.display()

    def addplot(self,display=None):
        """ addplot(display=display_defaults.json')
            reload display defaults file, and replot last data using
            display_defaults_file = optional--new display defaults file name
                                  = .json file containing display directives"""
        if display is not None:
            self.setdisplay(display)
        self.display()


    def clear(self):
        """ clear data from all figure windows"""

    # loop through all existing figures
        if self.figs is not None:
            self.figs.clear()
        self.repaint()

    def close(self):
        """ close all figure windows """
        if self.figs is not None:
            if True:
                self.figs.close()
            else:
                for x in self.figs:
                    f=self.figs.figure(x)
                    f.clf()
                    del f


    def repaint(self):
        if self.figs is not None:
            self.figs.shownew()

    def print_fig(self,fig_num,printer=None):
        """print requested figure number"""
        print 'printing figure # = ',fig_num
        f=plt.figure(fig_num)
        f.savefig('myfilename.pdf', format='pdf')
        if printer is None:
            if self.config is None:
                raise RuntimeError('Need printer parameter')
            printer = self.config.get_value('printer_id','name')
        printer = '-d' + printer[:]
        sp.call(['/usr/bin/lp', printer, 'myfilename.pdf'])

    def save_fig(self,fig_num,filename=None,dots_per_inch=None):
        """save requested figure number"""

        print 'saving figure # =',fig_num
        if filename is None:
            #default filename was provided in display defaults
            try:
                fig_dir = self.config.get_value('figure_download_directory','fig_dir')
                print
                print 'figure directory selected from display_defaults.json file'
                print 'figure will be written to: "'+fig_dir + '/filename.xxx"'
                print '   xxx = ',
                for item in plt.gcf().canvas.get_supported_filetypes():
                    print item + ' | ',
                print    
                filename = raw_input('filename.xxx = ? ')
                if not filename.startswith(('.','/')):
                    print 'Prefixing',fig_dir
                    filename = os.path.join(fig_dir ,filename)
            except:
                print 'display.json file did not provide a default directory'
                print 'figure will be written to: "path/filename.xxx"'
                print '   xxx = ',
                for item in plt.gcf().canvas.get_supported_filetypes():
                    print item + ' | ',
                print    
                filename = raw_input('path/filename.xxx = ? ')

        if dots_per_inch is None:
            dots_per_inch = float(raw_input('resolution (dpi) =? '))
        if isinstance(fig_num,basestring):
            f=self.figs.figure(fig_num)
        else:
            f=plt.figure(fig_num)
        f.savefig(filename,dpi = dots_per_inch,bbox_inches='tight')

    def print_cmd_list(self):
        print '\n'
        print """cmds-->
        r.next()    r.cnext()       r.print_fig(#)   r.cal_gen('cal_type') 
        r.clear()   r.list_vars()   r.inspect_data() r.print_mscat_defaults()
        r.replot()  r.save_fig(#)   r.cmp_calfile()  r.print_particle_defaults()   
        r.addplot() r.plot_model()  r.reprocess()    r.print_display_defaults()             
        r.close()   r.cmp_aeronet() r.cal_scale()    r.print_hsrl_defaults()
        r.loop()                    r.write_netcdf()

        """
    def reprocess(self):
        """ Reloads calvals_xxhsrl.txt, process_control.json, and xxx_plot.json
            then recomputes last request with new profiles plotted.'
            on figures over old plots."""
        self.reload_cals()
        self.update()
        self.display()

    def update(self,initOnly=False):
        dplc=None
        try:
            dplc=self.makeNewFramestream()
            self.rs_dpl=dplc
            self.rs=None
            if not initOnly:
                for f in dplc():
                    if self.rs is not None:
                        print 'Extra call?'
                    self.rs=f
        except KeyboardInterrupt:
            raise
        except:
            self.rs_dpl=None
            del dplc
            traceback.print_exc() #raising here would be nice, but the exception retains the mess, so python can't exit

    def display(self):
        if self.rs_dpl is None:
            return
        art=None
        try:
            art=self.makeNewImageArtistFor(self.rs_dpl)
            self.display_defaults=art.display_defaults
            for x in art():
                pass
            self.repaint()
            self.print_cmd_list()
        except KeyboardInterrupt:
            raise
        except:
            self.rs_dpl=None
            del art
            traceback.print_exc() #raising here would be nice, but the exception retains the mess, so python can't exit
        if not plt.isinteractive():
            print 'Matplotlib not interactive mode. calling ion()'
            plt.ion()
            
    def plot_model(self):
        print
        print
        print 'compute theoretical combined hi return for the current data'
        print 'segment using system optical parameters in "xxhsrl_performance_specs.json".'
        print 'The Rayleigh scattering cross section is taken from rs_Cxx.'
        print 'The geometric correction is taken from rs_cal.'
        print 'Combined_counts are raw dark corrected counts, with no pileup'
        print 'or baseline corrections.  Comparisons are only valid in clear air with minimal'
        print ' low altitude extinction.'
        print

        
        #Cmm_scale =  np.nanmax(self.rs_dpl.hsrl_cal.i2scan.data[:,2]/self.rs_dpl.hsrl_cal.i2scan.data[:,1])
        #print 'cmm_scale = ',Cmm_scale
        import copy
        Cxx = copy.deepcopy(self.rs_dpl.hsrl_Cxx)
        #Cxx.Cmm = Cxx.Cmm/Cmm_scale
        molecular_counts =    self.rs.profiles.dc_molecular_counts.copy()
        combined_hi_counts  = self.rs.profiles.dc_combined_hi_counts.copy()
        
        pm.performance_model(self.rs_dpl.hsrl_instrument,Cxx,self.rs_dpl.hsrl_cal
                             ,combined_hi_counts,molecular_counts,self.rs.rs_mean.transmitted_energy
                             ,self.rs.rs_inv.seeded_shots,self.rs.rs_mean.times,self.rs_dpl.hsrl_constants_first)



    def cal_scale(self):
        """ rescale or turn off selected data corrections  """

        last_i2_scale = self.rs_dpl.hsrl_corr_adjusts['i2_corr']
        last_i2a_scale = self.rs_dpl.hsrl_corr_adjusts['i2a_corr']
        last_Cam_scale = self.rs_dpl.hsrl_corr_adjusts['Cam_corr']

        corr_adjusts = gsc.get_scale_corrections(self.rs_dpl.hsrl_corr_adjusts,self.rs_dpl.hsrl_instrument)

    # if new i2 correction was supplied, get new calibration constants
        #FIXME if rs_dpl was a live dpl object that was being stepped by hand, corr_adjusts might be able to be overridden directly
        # but it isn't. update() currently always makes a new one, and runs it.
        self.searchparms['corr_adjusts']=corr_adjusts #updates the search parameter corr_adjusts override

        self.update()       #(inCalibration=True)
        self.display()

    ## location of different structures:
    # multiple scattering : self.rs_dpl.multiple_scattering_parameters
    # spheroid particle   : self.rs_dpl.spheroid_particle_parameters
    # mass_dim particle   : self.rs_dpl.mass_dimension_particle_parameters
    # radar               : not exposed (may have to use Rti's and assume its good... which should be safe)

    def print_hsrl_defaults(self):
        aeD ={}
        eD ={}
        dD ={}
        p_defaults = self.rs_dpl.hsrl_process_control.get_dict()
        for key in p_defaults:
            if p_defaults[key].has_key('enable')and p_defaults[key]['enable'] > 0:
                eD[key] = p_defaults[key]
            elif not p_defaults[key].has_key('enable'):
                aeD[key] = p_defaults[key]
            else:
                dD[key] = p_defaults[key]
        print '----------------disabled processes------------------------'        
        formatted_json_print(dD)
        print
        print
        print '-----------------always enabled processes or settings-----'
        formatted_json_print(aeD)
        print
        print
        print '-----------------enabled processes-------------------------'
        formatted_json_print(eD)

    def print_particle_defaults(self):
        print 
        print '----------------particle parameters----------------------'
        p_defaults = self.rs_dpl.spheroid_particle_parameters
        for key in p_defaults:
            if not key == 'parameters':
                if key == 'size_distribution':
                   print '     Size distribution parameters'
                   for key2 in p_defaults['size_distribution']: 
                      print '         ', key2, ' : ',p_defaults['size_distribution'][key2]
                else:    
                   print '    ',key,' : ', p_defaults[key]


    def print_mscat_defaults(self):
        print 
        print '----------------multiple scattering parameters-------------'
        p_defaults = self.rs_dpl.multiple_scattering_parameters
        for key in p_defaults:
            if not key == 'parameters':     
                print '    ',key,' : ', p_defaults[key]

    def print_display_defaults(self,d=None,k="",i=""):
        if d is None:
            for k,i,d in self.display_defaults_multisource_generator:
                self.print_display_defaults(d,k,i)
            return
        eD ={}
        dD ={}
        aeD={}
        p_defaults,unused_defaults = d.get_filtereddict()
        auD={}
        uD={}
        for key in p_defaults:
            if p_defaults[key].has_key('enable')and p_defaults[key]['enable'] > 0:
                eD[key] = p_defaults[key]
            elif not p_defaults[key].has_key('enable'):
                aeD[key]=p_defaults[key]
            else:
                dD[key] = p_defaults[key]
        for key in unused_defaults:
            if unused_defaults[key].has_key('enable'):
                uD[key] = unused_defaults[key]
            else:
                auD[key]=unused_defaults[key]
        if len(i)>0:
            prefix='%s:%s ' % (k,i)
        else:
            prefix=self.disp+' '
        print prefix+'-----------------unused/deprecated figures---------------'
        formatted_json_print(uD)
        print
        print
        print prefix+'-----------------unused/deprecated figure settings--------'
        formatted_json_print(auD)
        print
        print
        print prefix+'----------------disabled figures------------------------'        
        formatted_json_print(dD)
        print
        print
        print prefix+'-----------------global figure settings------------------'
        formatted_json_print(aeD)
        print
        print
        print prefix+'-----------------enabled figures-------------------------'
        formatted_json_print(eD)

    def print_radar_defaults(self):
        if self.libs.has_key('radar'):
            print self.libs['radar']
        else:
            print
            print "Can't print radar defaults--no radar instrument in request"
            print

    def list_vars(self):
        """print list of data structures and with variables and array sizes"""
        fields=['rs_mean','rs_inv','rs_inv.qa_flags','rs_raw','rs_Cxx','profiles','profiles.inv','raw_profiles','rs_init','rs_mmcr',\
                             'rs_spheroid_particle','rs_particle','rs_multiple_scattering','rs_marinemet','rs_pars2S1','rs_pars2S2','rs_vdis','rs_rain',
                             'rlprofaerosol','rlprofext','raman_rawprofile','raman_profile','raman_profile.inv','rlprofext_low','rlprofmerge',
                             'rlprofmerge_low','rlproftemp','rlprofdep','raman_inv','raman_hsrl_test','ramanmerge_hsrl_test','raman_hsrl_profile','rs_allradar_hsrl_coop']
        for f in vars(self.rs).keys():
            if not f.startswith('rs_'):
                continue
            if f in fields:
                continue
            fields.append(f)
        list_fields(self.rs,fields)
        return


    def inspect_data(self,request_time=None,request_alt=None,var_str=None):
        """ get information about a vairable at a particular time and altitude,
            invoke special debug processing routines if requested"""

        import lg_base.core.read_utilities as hru
        import lg_base.core.array_utils as hau
        from time import mktime

        #D = None
        D = self.rs
        
        for f in ('rs_mean','rs_inv','rs_mmcr','rs_marinemet','raman_inv','rlprofaerosol','rlproffex1thor'):
            if hasattr(D,f):
                D = getattr(D,f)
                break
        if D is None:
            print
            print 'Structure containing times and msl_altitudes not found'

            return
        print
        print f
        print D
        print dir(D)
        print
  
        times = getattr(D,'times')
        if request_time is None:
            print 'available times: %s  ---> %s ' %(times[0],times[-1])
            request_time = raw_input('requested time, e.g.(2013-09-18 10:00:00) = ?  ').strip()
            start_time_day_str = times[0].strftime('%Y-%m-%d')
            #if fractional second round to next lowest second
            if request_time.find('.') >= 0:
                index = request_time.find('.')
                request_time=request_time[:index]

            if len(request_time)==0:
                request_time=times[0].strftime("%Y-%m-%d %H:%M:%S")
                print 'Empty input. Using '+request_time
            if request_time[3:].find(':') < 0 :
                #request_time does not include  seconds--add seconds
                request_time = request_time + ':00'    
            if not request_time.find('-') > 0:
                #time does not include year-month_day--add these
                request_time = start_time_day_str + ' ' +request_time

            request_time = datetime.strptime(request_time,'%Y-%m-%d %H:%M:%S')    

        #dt = hau.T_Array([ abs((n - request_time).total_seconds()) for n in D['times'] ])
        #for i in range(len(times)):
        dt = hau.T_Array([ abs((times[i] - request_time).total_seconds()) for i in range(len(times))])
        t_index = dt.argmin()

        altitudes=getattr(D,'msl_altitudes')
        if request_alt is None:
            print 'available alitudes (km):  %3.2f--->%3.2f'\
                  %(altitudes[0]/1000.0,altitudes[-1]/1000.0)
            request_alt  = raw_input('requested altitude (km) = ?   ').strip()
            if len(request_alt)==0:
                request_alt='5.0'
        dz = np.abs(altitudes-np.float(request_alt)*1000.0)
        z_index = dz.argmin()
        print 
        print 'altitude = %4.1f, time = ' % altitudes[z_index],
        print times[t_index],
        print ',  index[%i,%i]' %(t_index,z_index)


        #get requested variables

        while 1:
            if var_str is None:
                print 
                var_str = raw_input('requested variable (e.g. struct.var or sp for special process) = ? ').strip()
            #leave routine if cariage return
            if len(var_str) == 0:
                break
            #do special processing if 'sp'
            if var_str == 'sp':
                import special_process_utilities.special_process as spu 
                spu.special_process(self, t_index, z_index)
            vi=var_str.split('.')
            D = self.rs
            try:
                for x in vi[:-1]:
                    D = getattr(D,x)
                variable = vi[-1]
                var=getattr(D,variable)
                print 'Array type = ', type(var)

                if isinstance(var,hau.TZ_Array):
                    print '%s[%i,%i] = %g' %(variable,t_index,z_index
                                             , var[t_index,z_index])
                elif isinstance(var,hau.T_Array):
                    print '%s[%i] = %g' %(variable,t_index
                                          , var[t_index])
                elif isinstance(var,hau.Z_Array):
                    print '%s[%i] = %g' %(variable,z_index,var[z_index])
                else:
                    print 'Unknown array type'

            except KeyError:
                self.list_vars()
                print
                print 'variable ',var_str,' not found--see above variable list'
            var_str=None

    def cmp_calfile(self,calfile=None):
        """ plot comparison of current default 'baseline','geofile','diff_geofile'
            ,or 'i2_scan' new calfile"""
        print ' '
        print "Supply 'baseline','geofile','diff_geo','i2a_temp_table',or 'i2_scan' filename with complete path "

        if calfile is None:
            calfile = raw_input('file name to compare #? ')
        calfile.strip()
        if not os.path.isfile(calfile):
            print
            print 'calibration filename = "'+calfile+'" was not found'
            print
            return
        import matplotlib.pylab as plt
        if 'geo' in calfile:
            if '/i2a_diff_geofile' in calfile:
                pass
            elif 'diff_default_geofile' in calfile:
                geo = self.rs_dpl.hsrl_cal.diff_geo
            elif 'geofile' in calfile:
                geo = self.rs_dpl.hsrl_cal.geo
                new_geo = hru.read_geo_corr(None,None,max_range_bin=geo.data.shape[0],filename=calfile)
                range_sq_corr = geo.data[:,0]**2
                range_sq_corr_new = new_geo.data[:,0]**2
            elif 'nadir_default_geofile' in calfile:
                pass
            elif 'i2a_d_geo' in calfile:
                pass 
            else:
                print
                print 'calibration file = "',calfile,'" was not found'
                print
                return



            print
            print 'Currently active calfile header'
            print self.rs_dpl.hsrl_cal.geo.header
            print
            print 'Replacement calfile header'
            print new_geo.header
            print
            plt.figure(2000)
            plt.plot(geo.data[:,0]/1000.0,geo.data[:,1] * 1e6 / range_sq_corr,'c'
                     ,geo.data[:,0]/1000.0,geo.data[:,1],'c'
                     ,new_geo.data[:,0]/1000.0,new_geo.data[:,1] * 1e6 / range_sq_corr,'r'
                     ,new_geo.data[:,0]/1000.0,new_geo.data[:,1],'r')
            plt.title(calfile)
            plt.xlabel('range (km)')
            plt.grid(True)
            plt.legend(['current','current*r^2','new','new*r^2'])
            ax=plt.gca()
            ax.set_yscale('log')
            plt.ylabel('Correction')



            plt.figure(2001)
            plt.plot(new_geo.data[:,0]/1000.0,new_geo.data[:,1]/geo.data[:len(new_geo.data[:,1]),1],'r')
            plt.title(calfile)
            plt.xlabel('range (km)')
            plt.grid(True)
            ax=plt.gca()
            plt.ylabel('new / curent_active')





            #replace currently active geofile with contents of calfile
            self.update_calibration_overrides(geo=new_geo)

        elif 'baseline' in calfile:    
            new_blc = hru.read_baseline(None,self.start_time,filename=calfile)
            new_data=new_blc.data
            blc = self.rs_dpl.hsrl_cal.baseline.data
            print
            print 'Currently active calfile header'
            print self.rs_dpl.hsrl_cal.baseline.header
            print
            print 'Replacement calfile header'
            print  new_blc.header
            print 


            print new_data
            plt.figure(2000)
            plt.plot(blc[:,0],blc[:,1],'c'
                     ,new_data[:,0],new_data[:,1],'r')
            plt.title(calfile)
            plt.xlabel('bin #')
            plt.grid(True)
            plt.legend(['current','new'])
            ax=plt.gca()
            ax.set_yscale('log')
            plt.ylabel('chi Baseline')

            plt.figure(2001)
            plt.plot(blc[:,0],blc[:,2],'c'
                     ,new_data[:,0],new_data[:,2],'k')
            plt.title(calfile)
            plt.xlabel('bin #')
            plt.grid(True)
            plt.legend(['current','new'])
            ax=plt.gca()
            ax.set_yscale('log')
            plt.ylabel('clo Baseline')

            plt.figure(2002)
            plt.plot(blc[:,0],blc[:,3],'c'
                     ,new_data[:,0],new_data[:,3],'b')
            plt.title(calfile)
            plt.xlabel('bin #')
            plt.grid(True)
            plt.legend(['current','new'])
            ax=plt.gca()
            ax.set_yscale('log')
            plt.ylabel('mol Baseline')

            plt.figure(2003)
            plt.plot(blc[:,0],blc[:,4],'c'
                     ,new_data[:,0],new_data[:,4],'g')
            plt.title(calfile)
            plt.xlabel('bin #')
            plt.grid(True)
            plt.legend(['current','new'])
            ax=plt.gca()
            ax.set_yscale('log')
            plt.ylabel('cpol Baseline')

            plt.figure(2004)
            plt.plot(blc[:,0],blc[:,5],'c'
                     ,new_data[:,0],new_data[:,5],'m')
            plt.title(calfile)
            plt.xlabel('bin #')
            plt.grid(True)
            plt.legend(['current','new'])
            ax=plt.gca()
            ax.set_yscale('log')
            plt.ylabel('I2a Baseline')

            """
            plt.figure(2001)
            plt.plot(new_geo.data[:,0]/1000.0,new_geo.data[:,1]/geo.data[:len(new_geo.data[:,1]),1],'r')
            plt.title(calfile)
            plt.xlabel('range (km)')
            plt.grid(True)
            ax=plt.gca()
            plt.ylabel('new / curent_active')
            """
            self.update_calibration_overrides(baseline=new_blc)

        elif 'i2-scan-' in calfile:
            new_i2 = hru.read_i2_scan(None,None,filename=calfile)
            new_data=new_i2.data
            i2_scan = self.rs_dpl.hsrl_cal.i2scan.data
            print
            print 'Currently active calfile header'
            print self.rs_dpl.hsrl_cal.i2scan.header
            print
            print 'Replacement calfile header'
            print new_data.header
            print 
            plt.figure(2000)
            plt.plot(i2_scan[:,0],i2_scan[:,1] ,'m' \
                     ,i2_scan[:,0],i2_scan[:,2],'c'
                     ,new_data[:,0],new_data[:,1],'r'\
                     ,new_data[:,0],new_data[:,2],'b')      
            plt.title(calfile)
            plt.xlabel('range (km)')
            plt.grid(True)
            plt.legend(['old_chi','old_mol','new_chi','new_mol'])
            ax=plt.gca()
            ax.set_yscale('log')
            plt.ylabel('Correction')



            plt.figure(2001)
            plt.plot(new_data[:,0],new_data[:,1]/i2_scan[:len(new_data[:,1]),1],'r'
                     ,new_data[:,0],new_data[:,2]/i2_scan[:len(new_data[:,2]),2],'b')
            plt.title(calfile)
            plt.xlabel('range (km)')
            plt.grid(True)
            ax=plt.gca()
            plt.ylabel('new / old')

            self.update_calibration_overrides(i2scan=new_i2)
        elif 'i2a_temp_table' in calfile:
            pass
        else:
            print
            print 'calibration file "',calfile,'" was not found'
            print
            return
        self.update()
        self.display()
        return


    def cmp_aeronet(self,aeronet_file_name=None,cmp_alt=None,plot_limit=None):
        """cmp_aeronet(self,aeronet_file_name,cmp_alt=None,plot_limit=None):
           read aeronet file and make plot comparing aeronet optical depth with
           hsrl optical depth  on existing od_vs_time plot
           
           """
        from hsrl.utilities.compare_hsrl_aeronet import compare_hsrl_aeronet
        if aeronet_file_name == None:
            print
            print "usage:   cmp_aeronet('path/aeronet_filename'"
            print         ",cmp_alt=compare_altitude(km),plot_limt=OD_plot_limit)"
            print "where: 'aeronet_filename' is a aeronet data file downloaded"
            print "from 'http://aeronet.gsfc.nasa.gov/'"
            print "will plot aeronet aerosol od's on existing od_vs_time plot"
            print "along with hsrl_ods from alt=cmp_alt vs aeronet values."
            print "Optical depth comparisions will be limited to OD's <= plot_limit"
            print "If aeronet file is .zip file will *.lev10 file will be automatically extracted"
            return

        if os.path.isfile(aeronet_file_name):
            compare_hsrl_aeronet(self,aeronet_file_name,cmp_alt,plot_limit)
        else:
            print
            print "usage:   cmp_aeronet('path/aeronet_filename',compare_altitude(km),OD_plot_limit)"
            print 'Aeronet file "' +aeronet_file_name+'" was not found.'
            print
            print 'Download aeronet file for the time period of interest'
            print 'from "http://aeronet.gsfc.nasa.gov/"'
            print "if aeronet_file is *.zip code will automatically extract *.lev10 files"
            print
            return


        return

    def cal_gen(self,cal_type):
        """
        Generate baseline, geofile, or diff_geofile cal files.'
        Geofiles can be computed from surface or constant altitude'
        flight legs with zenith pointing telescope.
        cal_type = 'geo'         --generate standard geofile from mol channel
                 = 'wfov_geo'    --geofile from wide-fov-channel
                 = 'baseline'    --baseline file from data with telescope covered
                 = 'i2scan'      --i2_scanfile from calibration scan data
                 = 'i2scan_manual' --i2_scanfile from manual i2 scan
                 = 'i2scan_from_spectra'--i2_scanfile from external i2 spectrum
                                          file
                 = 'i2a_temp_table'--temp table for use with argon buffered
                                         i2 cell
                 = 'diff_geo'      --diff_geofile, differential geometry
                                         between combined and molecular channels
                                         from data with i2 cell removed
                 = 'i2a_diff_geo'  --i2a_diff_geofile from data with i2a cell
                                         removed
                 = '1064_chi_gain'--ratio combined_1064 to combined_hi channel gain
                 = 'raman_wfov_geo'---wfov geofile for raman.
                 """
        if cal_type not in ['geo','wfov_geo','composite_geo','i2a_diff_geo','baseline','i2scan'
                            ,'cw_i2scan','i2scan_manual','i2scan_from_spectra'
                            ,'i2a_temp_table','diff_geo','1064_chi_gain'
                            ,'diff_geo_test','raman_wfov_geo' 
                            ,'1064_532_diff_geo','cpol_diff_geo']:
            print
            print 'UNKNOWN cal_type'
            print
            print "Valid cal_types:"
            print  "  'geo'           --generate standard geofile from mol channel"
            print  "  'wfov_geo'      --geofile from wide-fov-channel"
            print  "  'composite_geo' --geofile from wfov_geo patched together with standard geo"
            print  "  'baseline'      --baseline file from data with telescope covered"
            print  "  'i2scan'        --i2_scanfile from calibration scan data"
            print  "  'cw_i2scan'     --i2_scanfile from seperate cw seed laser"
            print  "  'i2scan_manual' --i2_scanfile from manual i2 scan"
            print  "  'i2scan_from_spectra'--i2_scanfile from external i2 spectrum"
            print  "                            file"
            print  "  'i2a_temp_table' --temp table for use with argon buffered"
            print  "                       i2 cell"
            print  "  'diff_geo'       --diff_geofile, differential geometry"
            print  "                      between combined and molecular channels"
            print  "                      from data with i2 cell removed"
            print  "  'i2a_diff_geo'   --i2a_diff_geofile from data with i2a cell"
            print  "                      removed"
            print  "  '1064_chi_gain'  --ratio combined_1064 to combined_hi channel gain"
            print  "  'raman_wfov_geo' --raman geofile from wfov /channel ratios"
            print  "  '1064_532_diff_geo'--1064/532 differential geometry correction"
            print  "  'cpol_diff_geo'  --cross pol differential geometry correction"
            return

        if cal_type == 'geo':
            if self.rs_dpl.hsrl_instrument == 'gvhsrl' \
               and self.rs_dpl.hsrl_constants_first['installation'] == 'airborne'\
               and any(self.rs.rs_raw.telescope_pointing == 0):
                cfg.make_geofile_new(
                    self.rs_dpl.hsrl_instrument
                    ,self.rs.profiles
                    ,self.rs_dpl.hsrl_cal
                    ,self.rs_dpl.hsrl_Cxx
                    ,self.rs_dpl.hsrl_constants_first
                    ,self.rs_dpl.corr_adjusts           
                    ,self.rs_dpl.hsrl_process_control)
                return    
            else: #lidar point upward
                #reinitialize with altitudes covering full range of lidar
                # and 7.5 m resolution
                if self.rs_dpl.hsrl_constants_first['installation'] == 'ground':   
                    self.searchparms['min_alt_m'] = self.rs_dpl.hsrl_constants_first['lidar_altitude']
                else:
                    print
                    print 'altitude variation during flight segment = '\
                          ,(np.max(self.rs.rs_raw.GPS_MSL_Alt) -np.min(self.rs.rs_raw.GPS_MSL_Alt)),' m'
                    print
                    self.searchparms['min_alt_m'] = np.mean(self.rs.rs_raw.GPS_MSL_Alt)

                self.searchparms['max_alt_m'] = 30000.0+self.searchparms['min_alt_m']
                self.searchparms['mol_norm_alt_m'] = self.searchparms['min_alt_m'] + 200
                self.setdisplay('all_plots.json')
                self.searchparms['altres_m']=None#none means best for image, or native
                self.searchparms['forimage']=False#setting this to false means is native

                #reprocess data with new resolution and no geo_corr
                if not self.searchparms.has_key('corr_adjusts'):
                    self.searchparms['corr_adjusts']={}
                self.searchparms['corr_adjusts']['geo_corr'] = 0.0
                self.update()
                #reprocessed data begins at lidar and has lidar instrument resolution

                #use repocessed data to compute geo_corr
                if cal_type == 'geo':
                    cfg.make_geofile(
                        self.rs_dpl.hsrl_instrument
                        ,self.rs.profiles
                        ,self.rs_dpl.hsrl_cal
                        ,self.rs_dpl.hsrl_Cxx
                        ,self.rs_dpl.hsrl_constants_first
                        ,process_defaults=self.rs_dpl.hsrl_process_control)
                         
                return
            
       
        
        elif cal_type == 'composite_geo' and hasattr(self.rs.rs_raw,'molecular_wfov_counts'):

            write_flg = False
            print
            print 'begining wfov geo correction'
            print

            # copied_profiles = copy.deepcopy(profiles)
            # copied_raw      = copy.deepcopy(raw)


 
            [wfov_ranges,wfov_geo,wfov_assumptions] =cfg.make_wfov_geofile(
                self.rs_dpl.hsrl_instrument
                ,self.rs.rs_raw
                ,self.rs.profiles
                ,self.rs_dpl.hsrl_Cxx
                ,self.rs_dpl.hsrl_cal
                ,self.rs_dpl.hsrl_constants_first
                ,self.rs_dpl.hsrl_process_control
                ,corr_adjusts=self.rs_dpl.hsrl_corr_adjusts
                ,write_flg = False)
            print
            print 'begining standard geo correction'
            print

            #lidar point upward
            #reinitialize with altitudes covering full range of lidar
            # and 7.5 m resolution
            if self.rs_dpl.hsrl_constants_first['installation'] == 'ground':   
                self.searchparms['min_alt_m'] = self.rs_dpl.hsrl_constants_first['lidar_altitude']
            else:
                print
                print 'altitude variation during flight segment = '\
                      ,(np.max(self.rs.rs_raw.GPS_MSL_Alt) -np.min(self.rs.rs_raw.GPS_MSL_Alt)),' m'
                print
                self.searchparms['min_alt_m'] = np.mean(self.rs.rs_raw.GPS_MSL_Alt)

            self.searchparms['max_alt_m'] = 30000.0+self.searchparms['min_alt_m']
            self.searchparms['mol_norm_alt_m'] = self.searchparms['min_alt_m'] + 200
            self.setdisplay('all_plots.json')
            self.searchparms['altres_m']=None#none means best for image, or native
            self.searchparms['forimage']=False#setting this to false means is native

            #reprocess data with new resolution and no geo_corr
            if not self.searchparms.has_key('corr_adjusts'):
                self.searchparms['corr_adjusts']={}
            self.searchparms['corr_adjusts']['geo_corr'] = 0.0
            self.update()
            #reprocessed data begins at lidar and has lidar instrument resolution

            #compute geo_corr from molecular return and expected clear air return
            [std_ranges,std_geo,std_assumptions]=cfg.make_geofile(
                self.rs_dpl.hsrl_instrument
                ,self.rs.profiles
                ,self.rs_dpl.hsrl_cal
                ,self.rs_dpl.hsrl_Cxx
                ,self.rs_dpl.hsrl_constants_first
                ,write_flg = False, fit_params = wfov_assumptions
                ,process_defaults=self.rs_dpl.hsrl_process_control)
            #combine wfov and conventional geo correction into composite
            cfg.make_composite_geofile(
                self.rs_dpl.hsrl_instrument
                ,self.rs_dpl.hsrl_cal
                ,self.rs.rs_raw.times
                ,wfov_ranges
                ,wfov_geo
                ,wfov_assumptions
                ,std_ranges
                ,std_geo
                ,std_assumptions
                ,process_defaults=self.rs_dpl.hsrl_process_control)                
            return
        elif cal_type == 'baseline':
            cfg.make_baseline_file(
                self.rs_dpl.hsrl_instrument
                ,self.rs.rs_raw
                ,self.rs_dpl.hsrl_constants_first
                ,self.rs_dpl.hsrl_process_control
                ,self.rs_dpl.hsrl_corr_adjusts)
            return
        elif cal_type == 'cw_i2scan':
             cfg.make_cw_i2_scan_file(
                self.rs_dpl.hsrl_instrument
                ,self.rs
                ,self.rs_dpl.hsrl_constants_first
                ,process_defaults=self.rs_dpl.hsrl_process_control)
             return
        elif cal_type == 'i2scan':
            cfg.make_i2_scan_file(
                self.rs_dpl.hsrl_instrument
                ,self.rs
                ,self.rs_dpl.hsrl_constants_first
                ,process_defaults=self.rs_dpl.hsrl_process_control)
            return
        elif cal_type == 'i2scan_manual':
            #use when narrow scan is not present
            cfg.make_i2_scan_file_manual(
                self.rs_dpl.hsrl_instrument
                ,self.rs
                ,self.rs_dpl.hsrl_constants_first
                ,process_defaults=self.rs_dpl.hsrl_process_control)
            return
        elif cal_type == 'i2scan_from_spectrum':
            #uses i2 scan from configuration files rather than current scan'
            cfg.make_i2_scan_from_i2_spectrum_file(
                self.rs_dpl.hsrl_instrument
                ,self.rs
                ,self.rs_dpl.hsrl_constants_first
                ,process_defaults=self.rs_dpl.hsrl_process_control)
            return

        elif cal_type == 'i2a_temp_table':
            if self.rs_dpl.hsrl_instrument == 'bagohsrl':
                cfg.make_temperature_table(
                    self.rs_dpl.hsrl_instrument
                    ,self.rs_dpl.hsrl_cal
                    ,self.rs.profiles.times
                    ,self.rs_dpl.hsrl_constants_first
                    ,self.rs_dpl.hsrl_corr_adjusts
                    ,self.rs_dpl.hsrl_process_control)
                return
            else:
                print 'no temp measurement channel on ',self.rs_dpl.hsrl_instrument
                return
        elif cal_type == '1064_532_diff_geo':
                #reinitialize with altitudes covering full range of lidar
                # and 7.5 m resolution   
                self.searchparms['min_alt_m'] = self.rs_dpl.hsrl_constants_first['lidar_altitude']
                self.searchparms['max_alt_m'] = 30000.0+self.searchparms['min_alt_m']
                self.searchparms['mol_norm_alt_m'] = self.searchparms['min_alt_m'] + 200
                #self.setdisplay('all_plots.json')
                self.searchparms['altres_m'] = None#none means best for image, or native
                self.searchparms['forimage']=False #setting this to false means is native

                #reprocess data with new resolution and no corr_adjusts['diff_1064_532_geo_corr']
                if not self.searchparms.has_key('corr_adjusts'):
                    self.searchparms['corr_adjusts']={}

                #disable diff_1064_532_geo_corr    
                self.searchparms['corr_adjusts']['diff_1064_532_geo_corr'] = 0.0

                #reprocess data
                self.update()
                
                #reprocessed data begins at lidar and has lidar instrument resolution
                cfg.make_diff_geo_1064_532(
                    self.rs_dpl.hsrl_instrument
                    ,self.rs.profiles.inv
                    ,self.rs_dpl.hsrl_constants_first
                    ,self.rs_dpl.hsrl_process_control)
                #reactivate diff_1064_532_geo_corr
                self.searchparms['corr_adjusts']['diff_1064_532_geo_corr'] = 1.0
            
        elif cal_type == 'i2a_diff_geo':
            cfg.make_i2a_diff_geofile(
                self.rs_dpl.hsrl_instrument
                ,self.rs.rs_raw
                ,self.rs_dpl.hsrl_Cxx
                ,self.rs_dpl.hsrl_cal
                ,self.rs_dpl.hsrl_process_control    
                ,self.rs_dpl.hsrl_constants_first
                ,self.rs_dpl.hsrl_corr_adjusts)
            return

        elif 0: # and cal_type == 'diff_geo' and self.rs_dpl.hsrl_instrument == 'gvhsrl':

            try:
                if self.rs.profiles.times.shape[0] < 1:
                    raise RuntimeError, "cal_gen - not enough times in profile"
            except IndexError:
                raise RuntimeError, "cal_gen - no profile exists - select data with 'r.params()'"

            # generate results at native altitude resolution without geo correction
            # only processses the last data chunk--does not work across cal or raob changes

            self.rs_dpl.hsrl_constants_first['first_bin_to_process'] = 0 #FIXME
            self.searchparms['min_alt_m'] = 0
            self.searchparms['max_alt_m'] = 60000
            self.searchparms['mol_norm_alt_m'] = 10000
            self.searchparms['altres_m']=None#none means best for image, or native
            self.searchparms['forimage']=False#setting this to false means is native
            if not self.searchparms.has_key('corr_adjusts'):
                self.searchparms['corr_adjusts']={}
            self.searchparms['corr_adjusts']['geo_corr'] = 0
            self.searchparms['corr_adjusts']['i2a_dif_geo_corr']=0

            self.update()
            #form average profiles selecting desired telescope pointing direction

        elif cal_type == 'diff_geo':
            print
            print 'making diff geofile'
            print
            cfg.make_diff_geofile(self.rs_dpl.hsrl_instrument
                                  ,self.rs.raw_profiles 
                                  ,self.rs_dpl.hsrl_cal
                                  ,self.rs_dpl.hsrl_constants_first
                                  ,self.rs_dpl.hsrl_process_control 
                                  ,self.rs_dpl.hsrl_corr_adjusts)
            self.display()
        elif cal_type == 'cpol_diff_geo':
            print
            print 'making diff cpol geofile'
            print
            cfg.make_diff_cross_pol(self.rs_dpl.hsrl_instrument
                                  ,self.rs.raw_profiles 
                                  ,self.rs_dpl.hsrl_cal
                                  ,self.rs_dpl.hsrl_constants_first
                                  ,self.rs_dpl.hsrl_process_control 
                                  ,self.rs_dpl.hsrl_corr_adjusts)
            self.display()
        elif cal_type == '1064_chi_gain' :
            #compute the gain ratio between combined_1064 and combined_hi
            cfg.compute_combined_1064_to_combined_hi_gain(
                self.rs.profiles   
                ,self.rs_dpl.hsrl_process_control)                       
        elif cal_type == 'raman_wfov_geo':
            #print 'self.rs.raman_profile'
            #print dir(self.rs.raman_profile)
            #print 'self.rs.raman_rawprofile'
            #print dir(self.rs.raman_rawprofile)
            chan_sel_dict = dict(sum_elastic_counts_high='sum_elastic_counts_low'
                                   ,sum_nitrogen_counts_high='sum_nitrogen_counts_low'
                                   ,sum_water_counts_high= 'sum_water_counts_low')
            print 'calling raman_wfov_geofile'
            rcfg.make_wfov_geofile(chan_sel_dict,self.rs,self.rs_dpl.raman_constants_first,self.rs_dpl.raman_process_control)
        else:  
            print
            print "Usage: r.cal_gen(type)"
            print "where type = ('geo' |'wfov_geo' |'baseline' | 'diff_geo'"
            print "            | i2a_diff_geo' | 'i2scan' | 'composite_geo' | 'i2scan_manual'"
            print "            | 'i2scan_from_spectrum' | 'i2a_temp_table' | '1064_chi_gain' "
            print "            | '1064_532_diff_geo')"
            print "cover telescope for 'baseline' correction"
            print "Use clear air data for 'geo' correction"
            print "Remove I2 cell for 'diff_geo' correction"
            print "Remove Argon buffered I2 cell for 'i2a_diff_geo' correction"

        return

    def write_netcdf(self,tag=None,netcdf_format=None,output_dir=None):
        print '\n'
        print '\n'
        print 'entering write_netcdf------------------------------------'
        print 'usage: write_netcdf(tag=None,netcdf_format=None,output_dir=None)'
        print "tag is "
        print tag
        print "netcdf_format is "
        print netcdf_format
        print "output_dir is "
        print output_dir
        print '\n'
        print '\n'




        timerange=[self.start_time,self.start_time+self.delta_t]
        start_time_str = timerange[0].strftime("%Y%m%dT%H%M")
        end_time_str = timerange[-1].strftime("%Y%m%dT%H%M")
        if tag is None:
            tag = raw_input('file name tag?  ')
        if tag == '':
            filename=self.instrument +'_'+start_time_str \
                +'_' + end_time_str +'.nc'
        else:
            filename=self.instrument +'_'+start_time_str \
                +'_' + end_time_str+'_' +tag+'.nc'
            
                  
        if netcdf_format is None:
            if self.config is None:
                raise RuntimeError('Need netcdf_format parameter')
            netcdf_format = self.config.get_value('netcdf','format')
        if output_dir is None:
            if self.config is None:
                raise RuntimeError('Need output_dir parameter')
            output_dir = self.config.get_value('netcdf','output_dir')
        filename = os.path.join(output_dir,filename)
        print ' ',filename
        if netcdf_format in ['uwlidar','uwlidar3','uwlidar3raw','hsrl_nomenclature.cdl','hsrl3_processed.cdl','hsrl3_raw.cdl']:#FIXME 20121207 JPG why are there 2 classes? this is supposed to be 1 exporter, configured by the CDL/initialization!
            import lg_dpl_toolbox.dpl.dpl_create_templatenetcdf as dpl_ctnc
            print 'writing ' + netcdf_format + ' format netcdf '
            template='hsrl_nomenclature.cdl'
            if netcdf_format=='uwlidar3':
                template='hsrl3_processed.cdl'
            elif netcdf_format=='uwlidar3raw':
                template='hsrl3_raw.cdl'
            elif netcdf_format.endswith('cdl'):
                template=netcdf_format
        elif 'cfradial' in netcdf_format:
            if netcdf_format.endswith('cdl'):
                template=netcdf_format
            else:
                template='hsrl_cfradial.cdl'
            print 'writing cfradial format netcdf to location ',output_dir
        else:
            template=netcdf_format

        self.makeNewNCOutputFor(self.rs_dpl,template,filename)
        return filename

    def loop(self,maxcount=-1):
        """ Automatically advance through data"""
        plt.ion()
        try:
            while maxcount!=0:

                self.figs.clear()
                try:
                    self.next()
                except RuntimeError, exc:
                    print exc
                    break
                except OSError, exc:
                    print exc
                    break
                # some log plots get np.ma.MaskErrors - known bug, no simple fix
                except np.ma.MaskError, exc:
                    print exc
                    break
                if maxcount>0:
                    maxcount-=1
        finally:
            plt.ioff()

    def getDisplayGeneratorFor(self,streamnames):
        print self.display_defaults_multisource_generator
        if isinstance(streamnames,basestring):
            streamnames=[streamnames]
        elif isinstance(streamnames,tuple):
            streamnames=list(streamnames)
        assert(isinstance(streamnames,list))
        #newDisplay=json.load(oc.open_config(self.disp), object_pairs_hook=OrderedDict)
        #for _,k,v in jdu.jsonParameterListGenerator(newDisplay,majorkey='display_defaults',keylist=streamnames+['default']):
        #    yield k,v
        for _,k,v in self.display_defaults_multisource_generator(keylist=streamnames+['default']):
            yield k,v

    def makeNewImageArtistFor(self,framestream):
        self.artistlist=OrderedDict()
        #artists = self.artists
        artist=framestream
        def allprofiles(frametypes):
            ret=[]
            for f in frametypes:
                if 'profile' in f:
                    ret.append(f)
            return ret

        descretealts=self.config.get_value('select_plot_altitude','plot_altitude','display')
        try:
            len(descretealts) #deliberately try to get check scalar or list/tuple/whatever
            descretealts=np.array(descretealts) #FIXME
        except:
            pass
        descretemaskname="plot_alt_index"
        rangemin=self.config.get_value('select_plot_layer','plot_alt_min','display')
        rangemax=self.config.get_value('select_plot_layer','plot_alt_max','display')
        rangemaskname='layer_indices'
        print 'slice from display config',descretealts,rangemin,rangemax
        if descretealts is not None:
            artist=aslice.DescreteAltitudeMask(artist,alts=descretealts*1000.0,maskname=descretemaskname)
        if rangemin is not None or rangemax is not None:
            artist=aslice.RangeAltitudeMask(artist,minalt=rangemin*1000.0,maxalt=rangemax*1000.0,maskname=rangemaskname)

        import hsrl.dpl.dpl_artists as hsrl_artists
        if 'hsrl' in self.libs or 'hsrl_raw' in self.libs or 'hsrl_profile' in self.libs:
          for d,p in self.getDisplayGeneratorFor('hsrl'):
            artist=hsrl_artists.dpl_images_artist(framestream=artist,
                                                  max_alt=self.searchparms['max_alt_m'],breakup_nesting=True,
                                                  display_defaults=d,**jdu.mergedDictionaryCopy(self.artistparams,p))
            self.artistlist['hsrl']=artist#FIXME we dont use the individual artists yet... but may in the future as they have their own sandbox
            self.display_defaults_multisource_generator.storeStructure(artist.display_defaults)
        if 'met' in self.libs:
            import met.dpl.dpl_artists as met_artists
            for d,p in self.getDisplayGeneratorFor('met'):
                artist=met_artists.dpl_met_images_artist(framestream=artist,display_defaults=d,framedomain='rs_marinemet'
                    ,**jdu.mergedDictionaryCopy(self.artistparams,p))
                self.artistlist['met']=artist
                self.display_defaults_multisource_generator.storeStructure(artist.display_defaults)
        if 'rain' in self.libs:
            import precip.dpl.rain_artist as rain_artist
            for d,p in self.getDisplayGeneratorFor('rain'):
                artist=rain_artist.dpl_rain_images_artist(framestream=artist,display_defaults=d,framedomain='rs_rain'
                    ,**jdu.mergedDictionaryCopy(self.artistparams,p))
                self.artistlist['rain']=artist
                self.display_defaults_multisource_generator.storeStructure(artist.display_defaults)
        if 'vdis' in self.libs:
            import precip.dpl.vdis_artist as vdis_artist
            for d,p in self.getDisplayGeneratorFor('vdis'):
                artist=vdis_artist.dpl_vdis_images_artist(framestream=artist,display_defaults=d,framedomain='rs_vdis'
                    ,**jdu.mergedDictionaryCopy(self.artistparams,p))
                self.artistlist['vdis']=artist
                self.display_defaults_multisource_generator.storeStructure(artist.display_defaults)
        if 'pars' in self.libs:
            import pars.dpl.dpl_artists as pars_artists
            for k,lib in self.libs['pars'].items():
                for d,p in self.getDisplayGeneratorFor(['pars',k]):
                    artist=pars_artists.dpl_pars_images_artist(framestream=artist,display_defaults=d,instrument=k
                        ,framedomain=k,**jdu.mergedDictionaryCopy(self.artistparams,p))
                    self.artistlist[k]=artist
                    self.display_defaults_multisource_generator.storeStructure(artist.display_defaults)
                  
        if 'raman' in self.libs:
            import raman.dpl.dpl_artists as raman_artists
            import cooperative.dpl.dpl_hsrl_raman_artists as hsrlraman_artists
            mergescope=None
            invscope=None
            for k,lib in self.libs['raman'].items():
                for d,p in self.getDisplayGeneratorFor([k,'raman_merge' if 'merge' in k else 'raman_inv' ,'raman']):
                    if 'merge' in k:
                        mergescope=k
                        artist=raman_artists.dpl_ramanmerge_images_artist(framestream=artist,display_defaults=d,streamname=k
                            ,subframe=k,**jdu.mergedDictionaryCopy(self.artistparams,p))
                    else:
                        invscope=k
                        artist=raman_artists.dpl_raman_images_artist(framestream=artist,display_defaults=d,streamname=k
                            ,subframe=k,**jdu.mergedDictionaryCopy(self.artistparams,p))
                    self.artistlist[k]=artist
                    self.display_defaults_multisource_generator.storeStructure(artist.display_defaults)
                if 'merge' in k:
                    for d,p in self.getDisplayGeneratorFor(['raman_rawprofiles','raman',k]):#FIXME
                        artist=raman_artists.dpl_raman_profile_images_artist(framestream=artist,streamname='raman_rawprofile',
                            display_defaults=d,subframe='raman_rawprofile',**jdu.mergedDictionaryCopy(self.artistparams,p))                        
                        self.display_defaults_multisource_generator.storeStructure(artist.display_defaults)
                    for d,p in self.getDisplayGeneratorFor(['raman_profiles','raman',k]):#FIXME
                        artist=raman_artists.dpl_raman_profile_images_artist(framestream=artist,streamname='raman_profile',
                            display_defaults=d,subframe='raman_profile',**jdu.mergedDictionaryCopy(self.artistparams,p))
                        self.display_defaults_multisource_generator.storeStructure(artist.display_defaults)
                else:
                    for d,p in self.getDisplayGeneratorFor(['raman_profiles','raman',k+'_profile']):#FIXME
                        artist=raman_artists.dpl_raman_inverted_profile_images_artist(framestream=artist,streamname=k+'_profile',
                            display_defaults=d,subframe=k+'_profile',**jdu.mergedDictionaryCopy(self.artistparams,p))
                        self.display_defaults_multisource_generator.storeStructure(artist.display_defaults)
            if 'raman_inv' in self.instrumentrequests:
                invscope='raman_inv'
                for d,p in self.getDisplayGeneratorFor(['raman_inv','raman']):#FIXME
                    artist=raman_artists.dpl_raman_images_artist(framestream=artist,display_defaults=d,subframe='raman_inv',streamname='raman_inv_UW'
                        ,**jdu.mergedDictionaryCopy(self.artistparams,p))
                    self.artistlist['raman_inv']=artist
                    self.display_defaults_multisource_generator.storeStructure(artist.display_defaults)
            if 'raman_hsrl_profile' in self.instrumentrequests:
                for d,p in self.getDisplayGeneratorFor(['raman_hsrl_profile','raman']):#FIXME
                    artist=hsrlraman_artists.dpl_raman_hsrl_profile_images_artist(framestream=artist,display_defaults=d,subframe='raman_hsrl_profile',
                        includenestedframes=dict(hsrl='profiles',raman='raman_profile') 
                        ,**jdu.mergedDictionaryCopy(self.artistparams,p))
                    self.artistlist['raman_hsrl_profile']=artist
                    self.display_defaults_multisource_generator.storeStructure(artist.display_defaults)
            if 'ramanmerge_hsrl_test' in self.instrumentrequests:
                for d,p in self.getDisplayGeneratorFor(['ramanmerge_hsrl_test','raman']):#FIXME
                    artist=hsrlraman_artists.dpl_ramanmerge_hsrl_images_artist(framestream=artist,display_defaults=d,subframe='ramanmerge_hsrl_test',
                        includenestedframes=dict(hsrl='rs_mean',raman=mergescope) 
                        ,**jdu.mergedDictionaryCopy(self.artistparams,p))
                    self.artistlist['ramanmerge_hsrl_test']=artist
                    self.display_defaults_multisource_generator.storeStructure(artist.display_defaults)
            if 'raman_hsrl_test' in self.instrumentrequests:
                for d,p in self.getDisplayGeneratorFor(['raman_hsrl_test','raman']):#FIXME
                    artist=hsrlraman_artists.dpl_raman_hsrl_images_artist(framestream=artist,display_defaults=d,subframe='raman_hsrl_test',
                        includenestedframes=dict(hsrl='rs_inv',raman=invscope)
                        ,**jdu.mergedDictionaryCopy(self.artistparams,p))
                    self.artistlist['raman_hsrl_test']=artist
                    self.display_defaults_multisource_generator.storeStructure(artist.display_defaults)
        import cooperative.dpl.dpl_all_profiles_artist as apa
        for d,p in self.getDisplayGeneratorFor(['profiles_test','all_profiles']):#FIXME
            artist=apa.dpl_profiles_images_artist(framestream=artist,display_defaults=d,
                subframes=allprofiles(framestream.provides.keys()),**jdu.mergedDictionaryCopy(self.artistparams,p))
            self.artistlist['all_profiles']=artist
            self.display_defaults_multisource_generator.storeStructure(artist.display_defaults)
        radars={}
        for inst in self.libs.keys():
            if inst=='mmcr' or 'kazr' in inst or 'mwacr' in inst:
                basetyp=None
                if inst=='mmcr':
                    basetyp=['mmcr']
                elif 'kazr' in inst:
                    basetyp=[inst[3:],'kazr']
                elif 'mwacr' in inst:
                    basetyp=['mwacr']
                radars['rs_'+basetyp[0]]=inst
                import radar.dpl.dpl_artists as radar_artists
                import cooperative.dpl.dpl_artists as coop_artists
                for d,p in self.getDisplayGeneratorFor([inst]+(basetyp if basetyp else [])+['radar']):
                    artist=radar_artists.dpl_radar_images_artist(framestream=artist,instrument=inst,display_defaults=d,subframe='rs_'+inst
                        ,**jdu.mergedDictionaryCopy(self.artistparams,p))
                    self.artistlist[inst]=artist
                    self.display_defaults_multisource_generator.storeStructure(artist.display_defaults)
                if 'mass_dimension_particle' in self.instrumentrequests:
                    for d,p in self.getDisplayGeneratorFor('mass_dimension_particle'):#FIXME
                        artist=coop_artists.dpl_particle_images_artist(framestream=artist,display_defaults=d
                            ,mode='mass_dimension',subframe='rs_particle',**jdu.mergedDictionaryCopy(self.artistparams,p))
                        self.artistlist['mass_dimension_particle']=artist
                        self.display_defaults_multisource_generator.storeStructure(artist.display_defaults)
                if 'spheroid_particle' in self.instrumentrequests:
                    for d,p in self.getDisplayGeneratorFor('spheroid_particle'):#FIXME
                        artist=coop_artists.dpl_particle_images_artist(framestream=artist,display_defaults=d,mode='spheroid',subframe='rs_spheroid_particle',
                                                                       includenestedframes={'rs_radar':'rs_'+inst,'rs_inv':'rs_inv','rs_pars2S1':'rs_pars2S1','rs_pars2S2':'rs_pars2S2','rs_pars2S1':'rs_pars2S1','rs_marinemet':'rs_marinemet'}
                                                                       ,**jdu.mergedDictionaryCopy(self.artistparams,p))
                        self.artistlist['spheroid_particle']=artist
                        self.display_defaults_multisource_generator.storeStructure(artist.display_defaults)
        if len(radars)>0:
            import cooperative.dpl.dpl_hsrl_radar_artists as hr_artist
            if 'allradar_hsrl_coop' in self.instrumentrequests:
                pp=None
                if 'allradar_hsrl_coop_parameters' in self.initkwargs:
                    pp=self.initkwargs['allradar_hsrl_coop_parameters']
                scopes=OrderedDict()
                scopes['rs_inv']='rs_inv'
                scopes['rs_mean']='rs_mean'
                for k,r in radars.items():
                    scopes['rs_'+k]='rs_'+r #framestream,display_defaults,subframe=None,figurecontainer=None,includenestedframes={}):
                for d,p in self.getDisplayGeneratorFor('allradar_hsrl_coop'):#FIXME
                    artist=hr_artist.dpl_allradar_hsrl_images_artist(artist,display_defaults=d,subframe='rs_allradar_hsrl_coop',includenestedframes=scopes,**jdu.mergedDictionaryCopy(self.artistparams,p))
                    self.artistlist['allradar_hsrl_coop']=artist
                    self.display_defaults_multisource_generator.storeStructure(artist.display_defaults)

        if 'multiple_scattering' in self.instrumentrequests:
            import cooperative.dpl.dpl_multiple_scattering_artist as ms_artist
            for d,p in self.getDisplayGeneratorFor('multiple_scattering'):#FIXME
                artist=ms_artist.dpl_multiple_scattering_artist(framestream=artist,display_defaults=d,subframe='rs_multiple_scattering'
                    ,includenestedframes={'profiles':'profiles'},**jdu.mergedDictionaryCopy(self.artistparams,p))
                self.artistlist['multiple_scattering']=artist
                self.display_defaults_multisource_generator.storeStructure(artist.display_defaults)
        return artist

    def assembleAttributeDictionary(self,prefix,content,omitkeys=[]):
        if isinstance(content,dict):
            ret=OrderedDict()
            keys=[k for k in content.keys()]
            keys.sort()
            for k in keys:
                if k in omitkeys or k.startswith('#'):
                    continue
                v=content[k]
                ret.update(self.assembleAttributeDictionary(prefix+'__'+k,v,omitkeys))
            return ret
        else:
            return {prefix:content}

    def makeNewNCOutputFor(self,framestream,template,filename):

        print '\n'
        print '\n'
        print 'entering makeNewNCOutputFor '
        print "framestream "
        print framestream
        print "template "
        print template
        print "filename "
        print filename
        print '\n'
        print '\n'

        
        import lg_dpl_toolbox.dpl.dpl_artists as artists
        output=None
        attrs=OrderedDict()
        try:
            consts=framestream.hsrl_constants_first
            attrs['hsrl_wavelength_nm']=consts['wavelength']
            if not 'installation' in consts or consts['installation'] == 'ground':
                attrs['hsrl_altitude_m']=consts['lidar_altitude']
                attrs['hsrl_latitude_degN']=consts['latitude']
                attrs['hsrl_longitude_degE']=consts['longitude']
        except:
            print 'Exception occurred loading attributes'
            traceback.print_exc()
        additional_attributes=OrderedDict()
        additional_attributes['hsrl_process_control']=('hsrl_processing_parameter',getdictionary)
        additional_attributes['radarLambda']=('radar_wavelength_m',None)
        additional_attributes['spheroid_particle_parameters']=None
        additional_attributes['raman_hsrl_test_parameters']=None
        additional_attributes['ramanmerge_hsrl_test_parameters']=None
        additional_attributes['mass_dimension_particle_parameters']=None
        additional_attributes['multiple_scattering_parameters']=None
        additional_attributes['allradar_hsrl_coop_parameters']=None
        omitkeys=('doc','docs','documentation','parameters','Parameters')
        for a,k in additional_attributes.items():
            if hasattr(framestream,a):
                attr=getattr(framestream,a)
                attrs.update(self.assembleAttributeDictionary(k[0] if k and k[0] else a,k[1](attr) if k and k[1] else attr,omitkeys))
        try:
            fmt,cfradial=artists.datasetParametersFor(template)
            if True:#not cfradial:
                cals=[]
                for cal in framestream.hsrl_cal_stream:
                    cals.append(cal)
                if len(cals)>0:
                    from netCDF4 import Dataset
                    output=Dataset(filename,'w',clobber=True,format=fmt)
                    import lg_dpl_toolbox.filters.substruct as frame_substruct
                    for f in artists.dpl_netcdf_artist(frame_substruct.TupleNarrator(cals,framestream.hsrl_cal_stream.provides),template,basetime=cals[0]['chunk_start_time'],\
                                                       output=output,withUnlimited=len(cals) if '3' in fmt else None,\
                                                       forModule=[sys.modules[type(framestream.hsrl_cal_stream).__module__],sys.modules[__name__],artists],
                                                       addAttributes=attrs):
                        pass
        except:
            traceback.print_exc()
            print 'exception occurred when trying to store calibration info'
            pass

        for f in artists.dpl_netcdf_artist(framestream,template,outputfilename=filename,output=output,forModule=[sys.modules[__name__],artists],addAttributes=attrs):
            pass
        if output:
            output.close()

    def update_calibration_overrides(self,**kwargs):
        for k,v in kwargs.items():
            if v is None:
                if self.calibration_overrides is not None and k in self.calibration_overrides:
                    del self.calibration_overrides[k]
                print 'Cleared calibration override for rs_cal.'+k
                continue
            if self.calibration_overrides is None:
                self.calibration_overrides=dict()
            self.calibration_overrides[k]=v
            print 'Set calibration override for rs_cal.'+k+' =',v

    def replace_calibration_overrides(self,**kwargs):
        self.calibration_overrides=None
        self.update_calibration_overrides(**kwargs)

    def makeNewFramestream(self):
        import lg_dpl_toolbox.filters.substruct as frame_substruct
        frame_substruct.SubstructBrancher.multiprocessable=False #Rti only works properly in single-process mode, and multiprocess queue will sometimes hang or crash.

        hsrlSource=None
        usedhsrl=None
        for hsrln in ('hsrl','hsrl_raw','hsrl_profile'):
            if hsrln in self.libs:
                hsrlSource=self.libs[hsrln]
                usedhsrl=hsrln
                break
        raw_only= (usedhsrl=='hsrl_raw')
        doinversion=(usedhsrl!='hsrl_profile')
        import hsrl.data_stream.set_resolution as sr
        import lg_dpl_toolbox.filters.time_frame as time_slicing
        import lg_dpl_toolbox.filters.resample_altitude as altitude_resampling

        def updated(dictionary,**kwargs):
            ret=dictionary.copy()
            ret.update(kwargs)
            return ret

        searchparms=self.searchparms.copy()
        searchparms['start_time_datetime']=self.start_time
        searchparms['end_time_datetime']=self.start_time+self.delta_t
        
        altitudes=None
        alt_host=self.alt_host
        time_host=self.time_host
        timesource=None
        framesource=None
        tslib=None
        tsgparms=dict()
        import lg_dpl_toolbox.dpl.TimeSource as TimeSource
        import lg_base.core.canvas_info as ci
        canvas_info=ci.load_canvas_info()

        if hsrlSource is not None:
            hsrlconst=hsrlSource(cal_only=True,raw_only=True,**searchparms)
            hsrlcal=None
            #p=hsrlnar.provides
            

            rs_constants=hsrlconst.hsrl_constants_first
            #number_x_pixels = process_control.get_value('image_pixel_size','x_pixels')
            #number_y_pixels = process_control.get_value('image_pixel_size','y_pixels')

            if 'altres_m' not in searchparms:
                binwidth = rs_constants['binwidth']*1.5e8
                [alt_res,n_range_ave] = sr.set_z_resolution( 
                    min_alt=searchparms['min_alt_m'], max_alt=searchparms['max_alt_m'], binwidth=binwidth
                    ,number_y_pixels=canvas_info['canvas_pixels']['y'], **({} if not self.z_res else self.z_res))        
                searchparms['altres_m']=alt_res

            if 'timeres_timedelta' not in searchparms:
                bin_dt = float(rs_constants['integration_time'])
                if not 'quarter_wave_plate_rotation' in rs_constants or  rs_constants['quarter_wave_plate_rotation'] == 'fixed' or rs_constants['quarter_wave_plate_rotation']=='none':
                    frame_dt = 0.0
                else:
                    frame_dt = float(rs_constants['polarization_integration_time'])

                [time_res,ntime_ave] = sr.set_time_resolution(self.delta_t, canvas_info['canvas_pixels']['x']
                                                              ,frame_dt,bin_dt, **({} if not self.t_res else self.t_res))
                searchparms['timeres_timedelta']=time_res

            if time_host is None:
                time_host='hsrl'
            if alt_host is None:
                alt_host='hsrl'

            if time_host == 'hsrl' or alt_host == 'hsrl':
                hsrlnar=hsrlSource(constsrc=hsrlconst,calsrc=hsrlcal,raw_only=raw_only,do_inversion=doinversion,calibration_overrides=self.calibration_overrides \
                    ,requested_altitudes=altitudes,**searchparms)
                dplc=hsrlnar
            if alt_host=="hsrl":
                altitudes=hsrlnar.altitudeAxis


        def callnarrate(obj,*args,**kwargs):
            return obj.narrateSubstruct(*args,**kwargs)

        def ginsu_rs_inv_time(hsrlnarsplitter,*args,**kwargs):
            return time_slicing.JustTime(time_slicing.TimeGinsu(hsrlnarsplitter.narrateSubstruct('rs_inv'),timefield='times',dtfield='delta_t',onlyTime=True),*args,**kwargs)
        def compound_rs_inv(hsrlnarsplitter):
            ret=hsrlnarsplitter.narrateSubstruct('rs_inv')
            return ret

        def noop(obj):
            return obj

        class JITTY(object):
            def __init__(self,sourceclass,obj,callfunc,*args,**kwargs):
                self.sourceclass=sourceclass
                self.object=obj
                self.callfunc=callfunc
                self.jit=None
                self.keepargs=args
                self.keepkwargs=kwargs

            def __call__(self,*args,**kwargs):
                if self.jit is None:
                    self.jit=self.sourceclass(self.object)
                if len(args)==0 and len(kwargs)==0:
                    return self.callfunc(self.jit,*self.keepargs,**self.keepkwargs)
                else:
                    if len(self.keepargs)>0 or len(self.keepkwargs):
                        print 'ignoring init args in JITTY',self.keepargs,self.keepkwargs
                    return self.callfunc(self.jit,*args,**kwargs)

        hsrlnarsplitter=None

        narration=OrderedDict()
        countMatch=[]
        radars=[]
        radartypes=[]

        for inst in self.libs.keys():
            if inst=='mmcr' or 'kazr' in inst or 'mwacr' in inst:
                radars.append(inst)
                basetyp=None
                if inst=='mmcr':
                    basetyp='mmcr'
                elif 'kazr' in inst or 'mwacr' in inst:
                    basetyp=inst[3:]
                radartypes.append(basetyp)
                if time_host is None:
                    time_host = 'rs_'+inst
                curs=self.libs[inst](**updated(searchparms,forimage=True,inclusive=time_host!='rs_'+inst))
                ginsuparms=dict(timefield='times',dtfield='delta_t' if 'delta_t' in curs.provides else None)
                narration['rs_'+inst]=time_slicing.TimeGinsu(altitude_resampling.ResampleXd(curs,'heights',altitudes) if ('rs_'+inst)!=alt_host else curs,
                                                        **ginsuparms )
                countMatch.append('rs_'+inst)
                if ('rs_'+inst)==time_host:
                    tslib=self.libs[inst]
                    tsgparms=ginsuparms
        if 'met' in self.libs:#add merge to rs_mmcr, refit
            if time_host is None:
                time_host = 'rs_met'
            ginsuparms=dict(timefield='times',dtfield=None)
            narration['rs_marinemet']=time_slicing.TimeGinsu(self.libs['met'](**updated(searchparms,forimage=True,inclusive=time_host!='rs_marinemet')),**ginsuparms)
            countMatch.append('rs_marinemet')
            if 'rs_marinemet'==time_host:
                tslib=self.libs['met']
                tsgparms=ginsuparms
        if 'rain' in self.libs:#add merge to rs_mmcr, refit
            if time_host is None:
                time_host = 'rs_rain'
            ginsuparms=dict(timefield='times',dtfield=None)
            narration['rs_rain']=time_slicing.TimeGinsu(self.libs['rain'](**updated(searchparms,forimage=True,inclusive=time_host!='rs_rain')),**ginsuparms)
            countMatch.append('rs_rain')
            if 'rs_rain'==time_host:
                tslib=self.libs['rain']
                tsgparms=ginsuparms
        if 'vdis' in self.libs:#add merge to rs_mmcr, refit
            if time_host is None:
                time_host = 'rs_vdis'
            ginsuparms=dict(timefield='times',dtfield=None)
            narration['rs_vdis']=time_slicing.TimeGinsu(self.libs['vdis'](**updated(searchparms,forimage=True,inclusive=time_host!='rs_vdis')),**ginsuparms)
            countMatch.append('rs_vdis')
            if 'rs_vdis'==time_host:
                tslib=self.libs['vdis']
                tsgparms=ginsuparms
        if 'pars' in self.libs:#add merge to rs_mmcr, refit
            ginsuparms=dict(timefield='times',dtfield=None)
            for k,lib in self.libs['pars'].items():
                if time_host is None:
                    time_host = k
                narration[k]=time_slicing.TimeGinsu(lib(**updated(searchparms,forimage=True,inclusive=time_host!=k)),**ginsuparms)
                countMatch.append(k)
                if k==time_host:
                    tslib=lib
                    tsgparms=ginsuparms
        if 'raman' in self.libs:
            for k,lib in self.libs['raman'].items():
                if time_host is None:
                    time_host = k
                nar=lib(**updated(searchparms,forimage=True,inclusive=time_host!=k))
                pp=None
                if 'merge' in k:
                    import raman.dpl.raman_inv_dpl as rmprof
                    narration['raman_rawprofile']=rmprof.dpl_raman_profile_filter(nar)
                    import raman.dpl.raman_inv_dpl as raman_inv
                    pp=None
                    if 'raman_inv_parameters' in self.initkwargs:
                        pp=self.initkwargs['raman_inv_parameters']
                    nar=raman_inv.RamanRangedFilter(nar,process_control=pp)
                ginsuparms=dict(timefield='times',dtfield=None,getTimeFunction=time_slicing.getTimeMiddleTime)#FIXME use start/width
                narration[k]=time_slicing.TimeGinsu(altitude_resampling.ResampleXd(nar,'altitudes',altitudes) if k!=alt_host else nar,**ginsuparms)
                if "merge" in k:
                    print 'APPLYING ROLLING TO MERGE RANGED DATA'
                    import lg_dpl_toolbox.filters.dpl_rolling_window_filter as drwf
                    fparms=raman_inv.getRollingDescription(nar.raman_process_control)
                    narration[k]=drwf.dpl_rolling_window_filter(narration[k],fparms)
                countMatch.append(k)
                if k==time_host:
                    tslib=lib
                    tsgparms=ginsuparms
        if time_host!='hsrl':# and (len(narration)>1 or hsrlSource is not None):
            print 'restarting hsrl narrator for new time source'
            #timesource=JITTY(frame_substruct.SubstructBrancher,narration[time_host],callnarrate,None)
            if timesource is not None:
                pass
            elif tslib is not None:
                timesource=JITTY(noop,time_slicing.TimeGinsu(tslib(**updated(searchparms,inclusive=False)),**tsgparms),time_slicing.JustTime,searchparms['start_time_datetime'],searchparms['end_time_datetime'])
                framesource=JITTY(noop,tslib(**updated(searchparms,inclusive=False)),noop)
            else:
                raise RuntimeError("FAIL")
                timesource=JITTY(noop,narration[time_host],time_slicing.JustTime,searchparms['start_time_datetime'],searchparms['end_time_datetime'])
                framesource=JITTY(noop,narration[time_host],noop)
            if hsrlSource is not None:
                hsrlnar=hsrlSource(raw_only=raw_only,do_inversion=doinversion,calibration_overrides=self.calibration_overrides\
                   ,timesource=timesource(),
                   requested_altitudes=altitudes,**updated(searchparms,inclusive=True))
                dplc=hsrlnar
            if tslib is None and framesource is not None and time_host in narration:
                narration[time_host]=framesource()
            print 'done'


        if len(narration.keys())>0: #must merge them
            import lg_base.core.array_utils as hau

            if hsrlnarsplitter is None and hsrlSource is not None:
                print 'making splitter for grouping size'
                hsrlnarsplitter=frame_substruct.SubstructBrancher(hsrlnar)
                print 'done making splitter'
            if timesource is None and time_host=='hsrl':
                timesource=JITTY(noop,hsrlnarsplitter,ginsu_rs_inv_time,searchparms['start_time_datetime'],searchparms['end_time_datetime'])
                print 'splitter given to time source'
            if framesource is None and hsrlSource is not None:
                framesource=JITTY(noop,hsrlnarsplitter,compound_rs_inv)
            import functools
            from dplkit.simple.blender import TimeInterpolatedMerge
            import lg_dpl_toolbox.dpl.dpl_configured_callable as dcc
            def funky(x,y):
                try:
                    print y,x
                except:
                    traceback.print_exc()
                    print 'cowardly ignoring error'
                return x

            if time_host not in ('hsrl','common') and hsrlnarsplitter is not None: #need to slice up hsrl, which is already resampled using timesource above
                ginsuparms=dict(timefield='times',dtfield='delta_t',getTimeFunction=time_slicing.getTimeMiddleTime)#FIXME use start/width
                for f in ('rs_inv','rs_mean','rs_raw'):
                    if f in hsrlnarsplitter.provides:
                        print 'splitting up hsrl source ',f
                        narration[f]=time_slicing.TimeGinsu(hsrlnarsplitter.narrateSubstruct(f),**ginsuparms)
                        if 'raw' not in f:
                            countMatch.append(f)

            merging=OrderedDict()
            print 'merging begins'

            for k,narr in narration.items():
                if k.startswith(('rs_mmcr','rlprof')) or 'kazr' in k or 'mmcr' in k or 'mwacr' in k:
                    if timesource is not None:#k!=time_host:
                        timenar=timesource()
                        sp=TimeInterpolatedMerge(timenar,[narr], allow_nans=True)
                        sp=TimeSource.AddAppendableTime(sp,'times','delta_t')
                        sp=frame_substruct.Retyper(sp,functools.partial(hau.Time_Z_Group,timevarname='times',altname='heights' if 'rlprof' not in k else 'altitudes'))
                    else:
                        sp=TimeSource.AddAppendableTime(narr,'times','delta_t')
                    if framesource is not None:
                        if k in countMatch:
                            sp=frame_substruct.CountDeGinsu(frame_substruct.FrameLength(framesource(),'times'),sp)
                        else:
                            sp=frame_substruct.TimeDeGinsu(framesource(),sp)
                elif k.startswith(('rs_marinemet','rs_rain','rs_vdis')):
                    if timesource is not None:#k!=time_host:
                        timenar=timesource()
                        sp=TimeInterpolatedMerge(timenar,[narr], allow_nans=True)
                        sp=TimeSource.AddAppendableTime(sp,'times','delta_t')
                        sp=frame_substruct.Retyper(sp,functools.partial(hau.Time_Z_Group,timevarname='times'))
                    else:
                        sp=TimeSource.AddAppendableTime(narr,'times','delta_t')
                    if framesource is not None:
                        if k in countMatch:
                            sp=frame_substruct.CountDeGinsu(frame_substruct.FrameLength(framesource(),'times'),sp)
                        else:
                            sp=frame_substruct.TimeDeGinsu(framesource(),sp)
                elif k.startswith(('rs_pars',)):
                    if timesource is not None:
                        timenar=timesource()#,searchparms['start_time_datetime'],searchparms['end_time_datetime']) #time_slicing.TimeGinsu(hsrlnarsplitter.narrateSubstruct('rs_inv'),timefield='times',dtfield='delta_t',onlyTime=True)
                        if True and  'pars' in k:
                            narr=dcc.dpl_configured_callable(narr,funky,y=narr.parsType+' sliced')
                        sp=time_slicing.Nearest(narr,timenar)
                        sp=TimeSource.AddAppendableTime(sp,'times','delta_t')
                        x=sp.parsType
                        if True and 'pars' in k:
                            sp=dcc.dpl_configured_callable(sp,funky,y=sp.parsType+' appended')
                    else:
                        assert(0)
                        sp=TimeSource.AddAppendableTime(narr,'times','delta_t')
                    if framesource is not None:
                        if True and 'pars' in k:
                            sp=dcc.dpl_configured_callable(sp,funky,y=sp.parsType+' pre-regrouped')
                        if k in countMatch:
                            sp=frame_substruct.CountDeGinsu(frame_substruct.FrameLength(framesource(),'times'),sp)
                        else:
                            sp=frame_substruct.TimeDeGinsu(framesource(),sp)
                        if True and 'pars' in k:
                            sp=dcc.dpl_configured_callable(sp,funky,y=sp.parsType+' regrouped')
                elif 'qa_flags' in k or k in ('rs_raw','rs_mean','rs_inv'):
                    sp=narr
                    sp=TimeSource.AddAppendableTime(sp,'times','delta_t')
                    if framesource is not None:
                        if k in countMatch:
                            sp=frame_substruct.CountDeGinsu(frame_substruct.FrameLength(framesource(),'times'),sp)
                        else:
                            sp=frame_substruct.TimeDeGinsu(framesource(),sp)
                    if False:
                        tmp=frame_substruct.SubstructBrancher(sp)
                        sp=tmp.narrateSubstruct('qaflags')
                elif 'profile' in k:
                    sp=narr
                else:
                    raise NotImplementedError(k)
                merging[k]=sp

            if hsrlSource is not None and len(merging)>0:
                print 'combining streams in compositor'
                dplc=frame_substruct.NestingCompositer(hsrlnarsplitter.narrateSubstruct(None),merging)
            elif len(merging)>0:
                if len(merging)>1:
                    if framesource is not None:
                        dplc=frame_substruct.NestingCompositer(
                            frame_substruct.DropFrameContent(frame_substruct.Nester(framesource(),"garb",hau.Time_Z_Group()),'garb'),merging)
                    else:
                        raise NotImplementedError("OOPS more than one with no pace")
                else:
                    for k,v in merging.items():
                        dplc=frame_substruct.Nester(v,k,hau.Time_Z_Group())
            else:
                pass
            #print dplc.provides
            #dplc=frame_substruct.SubstructMerger('rs_inv',merging,hau.Time_Z_Group,countMatch=countMatch)
        if 'raman' in self.libs and 'raman_inv' in self.instrumentrequests:#add merge to rs_mmcr, refit
            import raman.dpl.raman_inv_dpl as raman_inv
            pp=None
            if 'raman_inv_parameters' in self.initkwargs:
                pp=self.initkwargs['raman_inv_parameters']
            found=False
            for f in self.libs['raman'].keys():
                if 'merge' in f:
                    inv_raman=raman_inv.RamanInvertFilter(dplc,process_control=pp,datasourcescope=f,destinationscope='raman_inv',requested_altitudes=altitudes)
                    dplc=inv_raman
                    found=True
                    dplc=raman_inv.dpl_raman_profile_filter(dplc,srcscopename=f,subscopename='raman_profile'
                        ,calsource=inv_raman.cals)
                else:
                    dplc=raman_inv.dpl_raman_inverted_profile_filter(dplc,srcscopename=f,subscopename=f+'_profile')
            assert(found)
        ms_particle=None
        if len(radars)>0 and 'hsrl' in self.libs:#add merge to rs_mmcr, refit
            import cooperative.dpl.dpl_hsrl_radar as dpl_hsrl_radar
            if 'allradar_hsrl_coop' in self.instrumentrequests:
                pp=None
                if 'allradar_hsrl_coop_parameters' in self.initkwargs:
                    pp=self.initkwargs['allradar_hsrl_coop_parameters']
                scopes=OrderedDict()
                scopes['rs_inv']='rs_inv'
                scopes['rs_mean']='rs_mean'
                for i,r in enumerate(radars):
                    scopes['rs_'+radartypes[i]]='rs_'+r
                dplc=dpl_hsrl_radar.dpl_allradar_hsrl_cooperative(dplc,pp,scopes=scopes,outputscope='rs_allradar_hsrl_coop')
            if 'mass_dimension_particle' in self.instrumentrequests:
                pp=None
                if 'mass_dimension_particle_parameters' in self.initkwargs:
                    pp=self.initkwargs['mass_dimension_particle_parameters']
                hsrl_merge=dpl_hsrl_radar.dpl_mass_dimension_particle(dplc,pp,radarscope='rs_'+radars[0],hsrlscope='rs_inv',outputscope='rs_particle')
                dplc=hsrl_merge
                ms_particle='rs_particle'
            if 'spheroid_particle' in self.instrumentrequests:
                pp=None
                if 'spheroid_particle_parameters' in self.initkwargs:
                    pp=self.initkwargs['spheroid_particle_parameters'] 
                hsrl_merge=dpl_hsrl_radar.dpl_spheroid_particle(dplc,pp,radarscope='rs_'+radars[0],hsrlscope='rs_inv',outputscope='rs_spheroid_particle')
                dplc=hsrl_merge
                ms_particle='rs_spheroid_particle'
        if 'raman' in self.libs and 'hsrl' in self.libs:#add merge to rs_mmcr, refit
            import cooperative.dpl.dpl_hsrl_raman as dpl_hsrl_raman
            if 'ramanmerge_hsrl_test' in self.instrumentrequests:
                pp=None
                if 'ramanmerge_hsrl_test_parameters' in self.initkwargs:
                    pp=self.initkwargs['ramanmerge_hsrl_test_parameters']
                found=False
                for f in self.libs['raman'].keys():
                    if 'merge' in f:
                        hsrl_raman=dpl_hsrl_raman.dpl_hsrl_ramanmerge(dplc,pp,ramanscope=f,hsrlscope='rs_mean',outputscope='ramanmerge_hsrl_test')
                        dplc=hsrl_raman
                        found=True
                if not found:
                    raise RuntimeError("ramanmerge_hsrl_test needs a merge product")
            if 'raman_hsrl_profile' in self.instrumentrequests:
                pp=None
                if 'ramanmerge_hsrl_test_parameters' in self.initkwargs:
                    pp=self.initkwargs['ramanmerge_hsrl_test_parameters']
                hsrl_raman=dpl_hsrl_raman.dpl_hsrl_raman_profile(dplc,pp,ramanscope='raman_profile',hsrlscope='profiles',outputscope='raman_hsrl_profile')
                dplc=hsrl_raman
            if 'raman_hsrl_test' in self.instrumentrequests:
                pp=None
                if 'raman_hsrl_test_parameters' in self.initkwargs:
                    pp=self.initkwargs['raman_hsrl_test_parameters']
                if 'raman_inv' in self.instrumentrequests:
                    hsrl_raman=dpl_hsrl_raman.dpl_hsrl_raman(dplc,pp,ramanscope='raman_inv',hsrlscope='rs_inv',outputscope='raman_hsrl_test')
                    dplc=hsrl_raman
                else:
                    found=False
                    for f in self.libs['raman'].keys():
                        if 'merge' not in f:
                            hsrl_raman=dpl_hsrl_raman.dpl_hsrl_raman(dplc,pp,ramanscope=f,hsrlscope='rs_inv',outputscope='raman_hsrl_test')
                            dplc=hsrl_raman
                            found=True
                    if not found:
                        raise RuntimeError("raman_hsrl_test needs an inverted (ext or dep) product")
        if 'multiple_scattering' in self.instrumentrequests:
            import cooperative.dpl.dpl_multiple_scattering as ms
            pp=None
            if 'multiple_scattering_parameters' in self.initkwargs:
                pp=self.initkwargs['multiple_scattering_parameters']
            hsrl_merge=ms.dpl_multiple_scattering(dplc,pp,particlescope=ms_particle,hsrlscope='rs_inv',outputscope='rs_multiple_scattering')
            dplc=hsrl_merge
        for i,x in enumerate(self.filterclasses):
            dplc=x(dplc,**self.filterparams[i])

        if self.dropcontent is not None:
            import lg_dpl_toolbox.filters.substruct as frame_substruct
            dplc=frame_substruct.DropFrameContent(dplc,self.dropcontent)

        dplc=time_slicing.FrameCachedConcatenate(dplc)

        if len(self.fullfilterclasses)>0:
            for i,x in enumerate(self.fullfilterclasses):
                dplc=x(dplc,**self.fullfilterparams[i])

            dplc=time_slicing.FrameCachedConcatenate(dplc)#just to make sure modifications of the above are saved too.

        if 'timeres_timedelta' in searchparms:
            ts=TimeSource.TimeGenerator(start_time=searchparms['start_time_datetime']
                                       ,end_time=searchparms['end_time_datetime']
                                       ,time_resolution=searchparms['timeres_timedelta'])
            import lg_dpl_toolbox.filters.fill as fill
            def ignoreSomeGroups(name):
                return 'raw' in name or 'profile' in name
            dplc=fill.FillIn(dplc,[np.array([x['start'] for x in ts]),altitudes],ignoreGroups=ignoreSomeGroups)

        return dplc

    @staticmethod
    def available(baseinstrument=None,site=None,start=None,end=None,duration=None,include_derivative=False):
        """ Discover what datasets are available
        :param baseinstrument: instrument name of the HSRL identifying the site (HSRL-centricity). specify this or site, not both
        :param site: NOTIMPLEMENTED site enumeration (agnostic). specify this or baseinstrument, not both.
        :param start: start time to explore (default now). can be a datetime or a string (optional)
        :param end: end time to explore. can be a datetime or a string (optional). don't specify this and duration
        :param duration: duration, if specified in preference over end. is a timedelta or hours as a float. default is 5 minutes if end and duration not specified
        :param include_derivative: (default false) include all names of compatible derivative streams (like particle measurements)
        :return: tuple of available instrument streams, including optionally derivative

        example:

          Rti.available('mf2hsrl',start='20-Sep-13 5:00',include_derivative=True)
        will return:
          ('mf2hsrl', 'magkazrge', 'magkazrmd', 'magpars2S1', 'magpars2S2', 'mass_dimension_particle', 'spheroid_particle', 'multiple_scattering')
        """
        ret=[]
        libs=OrderedDict()
        if site is not None:
            raise NotImplementedError
        elif baseinstrument is not None:
            possible_instruments={'ahsrl':['ahsrl','mmcr'],
                                  'gvhsrl':['gvhsrl'],
                                  'mf2hsrl':['mf2hsrl','mf2kazrge','mf2kazrmd','mf2mwacr','mf2pars2S1','mf2pars2S2','mf2vdis','mf2met','mf2marinemet','mf2rain'
                                             ,'magkazrge','magkazrmd','magmwacr','magpars2S1','magpars2S2','magmarinemet'
                                             ,'tmpkazrge','tmpkazrmd','tmpmwacr','tmpvdis','tmpmet','tmprain'],
                                  'nshsrl':['nshsrl','nskazrge','nskazrmd','nsmwacr'
                                            ,'nsakazrge','nsakazrmd','nsamwacr'],
                                  'bagohsrl':['bagohsrl','rlprofaerosol','rlprofext','rlprofext_low','rlprofmerge','rlprofmerge_low','rlprofmerge_high','rlproftemp','rlprofdep'],
                                  'rbhsrl':['rbhsrl']}
            if baseinstrument not in possible_instruments:
                raise RuntimeError(baseinstrument,'not supported. must be one of',possible_instruments.keys())
            DOEsuff= dict(mf2hsrl='M1',nshsrl='C1')
            DOEinf=DOEsuff[baseinstrument] if baseinstrument in DOEsuff else ''
            from radar.dpl.MWACRLibrarian import MWACRLibrarian
            from pars.dpl.PARSLibrarian import PARSLibrarian
            from radar.dpl.KAZRLibrarian import KAZRLibrarian
            from precip.dpl.RainLibrarian import RainLibrarian
            from precip.dpl.VDisLibrarian import VDisLibrarian
            armlist=dict(mwacr=MWACRLibrarian,pars=PARSLibrarian,kazr=KAZRLibrarian,rain=RainLibrarian,vdis=VDisLibrarian)
            for inst in possible_instruments[baseinstrument]:
                if inst.endswith('hsrl'):
                    from hsrl.dpl.HSRLLibrarian import HSRLLibrarian
                    libs[inst]=HSRLLibrarian(instrument=inst)
                elif inst=='mmcr':
                    from radar.dpl.MMCRMergeLibrarian import MMCRMergeLibrarian
                    libs[inst]=MMCRMergeLibrarian('ahsrl',['eurmmcrmerge.C1.c1.','nsaarscl1clothC1.c1.'])
                elif 'met' in inst:
                    from met.dpl.MarinemetLibrarian import MarinemetLibrarian
                    libs[inst]=MarinemetLibrarian(baseinstrument,[inst])
                else:
                    for pref,lib in armlist.items():
                        if pref in inst:
                            libs[inst]=lib(baseinstrument,inst,[inst])
                            break
                if inst not in libs:
                    raise NotImplementedError('instrument',inst)
        else:
            raise RuntimeError('Must specify baseinstrument or site')

        import lg_dpl_toolbox.filters.time_frame as time_frame
        import lg_base.core.read_utilities as hru
        if isinstance(start,basestring):
            rs_date = hru.convert_date_str(start)
            oldstart=start
            start = rs_date['datetime']
        if isinstance(end, basestring):

            # check for time without date--assume it is same day as start

            if end.find(':') < 3:
                index = oldstart.find(' ')
                end = oldstart[0:index + 1] + end
            rs_date = hru.convert_date_str(end)
            end = rs_date['datetime']
        if end is None and duration is None:
            duration=5.0/60
        if duration is not None and not isinstance(duration,timedelta):
            duration = timedelta(days=(duration) / 24.0)  # plot length in hours converted to timedelta

        windowinfo=time_frame.parse_timewindow(start,end,duration)
        if windowinfo['endtime'] is None:
            windowinfo['endtime']=windowinfo['starttime']+windowinfo['windowwidth']
        for instrument,lib in libs.items():
            if libraryHas(lib,windowinfo['starttime'],windowinfo['endtime']):
                ret.append(instrument)

        if include_derivative:
            havehsrl=False
            haveradar=False
            haveraman=False
            for inst in ret:
                if inst.endswith('hsrl'):
                    havehsrl=True
                elif inst=='mmcr' or 'kazr' in inst or 'mwacr' in inst:
                    haveradar=True
                elif inst.startswith('rlprof'):
                    haveraman=True
            if havehsrl and haveradar:
                ret.extend(['mass_dimension_particle','spheroid_particle','allradar_hsrl_coop'])
            if havehsrl and haveraman:
                ret.extend(['raman_hsrl_test','ramanmerge_hsrl_test','raman_hsrl_profile'])
            if havehsrl:
                ret.extend(['multiple_scattering'])
        return tuple(ret)
#this is ripped from Picnic.views.py
def libraryHas(lib,starttime,endtime):
    times=[]
    fn=None
    t=None
    srchres=lib(start=starttime,end=endtime)
    for x in srchres:
        times.append(x['start'])
        if len(times)==2:
            break
        fn=x
        t=times[0]

    success=False
    if len(times)==0:
        success=False
        #print 'no data'
    elif len(times)>=2:
        success=True
        #print 'more than 1'
    elif t>=starttime and t<=endtime:
        success=True
        #print 'time in range'
    elif starttime>datetime.utcnow():
        success=False
    elif (fn['start']+fn['width'])>starttime and fn['start']<endtime:
        success=True
        #print 'end times work'
    #elif 'data' in x and (starttime-t).total_seconds()<(60*60):
    #    success=True
    #    print 'data time may intersect'
    #elif 'data' not in x and (starttime-t).total_seconds()<(3*60*60):
    #    success=True
    #    print 'cal time may intersect'
    return success
