#!/usr/bin/python
# -*- coding: utf-8 -*-

print """
*******************************************************
WARNING - this module iS ***NOT*** Supported -- WARNING
replace your import with 

from maestro.rti_maestro import Rti
*******************************************************

"""



import sys
import os.path
from time import sleep
from subprocess import call
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime,timedelta
import lg_base.core.read_utilities as hru
import hsrl.data_stream.hsrl_read_utilities as hhru
import hsrl.calibration.calibration_utilities as cu
import hsrl.calibration.cal_read_utilities as cru
import hsrl.soundings.sounding_utilities as su
import hsrl.data_stream.processing_utilities as pu
import hsrl.data_stream.display_utilities as du
import lg_base.graphics.graphics_toolkit as gt
import hsrl.calibration.cal_file_generation as cfg
from hsrl.data_stream.get_scale_corrections import get_scale_corrections
from lg_base.core.open_config import open_config
import lg_base.core.array_utils as hau
import hsrl.data_stream.main_routines as mr
import hsrl.data_stream.performance_model as pm
import hsrl.data_stream.input_translators as it
import hsrl.data_stream.set_resolution as sr
import hsrl.data_stream.iodine_argon_temperatures as i2at
import json
#from hsrl.utils.json_config import json_config
import lg_base.core.json_config as jc
from netCDF4 import Dataset
import copy

from lg_base.core.locate_file import locate_file

class Rti:

    """ HSRL Plotting Routines

  Examples:
  r=Rti('gvhsrl','18-Aug-11 00:00', 2, 0, 15)
  r=Rti('gvhsrl', '18-aug-11 00:00','18-aug-11 2:00' 0 15)
  r=Rti('ahsrl', '18-aug-11 1:00','23:21:05', 0, 15)
  r=Rti('ahsrl','1-jul-12 1:00','2:21',0,15,mol_norm_alt=2.1,display='all_plots.json'
              ,netcdf_defaults=None, cmd_defaults=None)
  optional cmd line parameters--defaults are provided in "process_control.json" 
    mol_norm_alt------altitude at which to normalize optical depth 
    display-----------select display options other than "display_defaults.json" 
    cmd_defaults------select process options other than "process_control.json"
    netcdf_defaults---selection and translation of netcdf variable names to read
                         will be controlled by contents of xxxx_netcdf_defaults.json
                         where xxxx is the instrument('eg ahsrl, mf2hsrl, bagohsrl ...)
  r.next()

  """
    def __init__( self, instrument, start_time, plot_length, min_alt 
            ,max_alt, mol_norm_alt=None, display=None,cmd_defaults=None
            ,netcdf_defaults=None,z_res=None,t_res=None):
             
        """ constructs a new Rti instance

        instrument     = hsrl id string (eg. 'ahsrl','gvhsrl','nshsrl','mf2hsrl').
        start_time     = time string, first data to plot (eg.'5-may-11 12:30').
        plot_length    = time period to plot in hours (e.g. 2) or
                              end time of first interval(e.g. '6-may-10 1:30').
        min_alt        = minimum altitude (msl) to display (km).
        max_alt        = maximum altitude (msl) to display (km).
        mol_norm_alt    = normalization altitude for optical depth (km).
                              if not supplied mol_norm_alt is taken from
                              process_control.json                                
        display        = name of  .json file containing display directives
                              (default = 'flight_plots.json')
        cmd_defaults   = alternate to xxxx_process_control.json file containing
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
        """
   
        plt.ion()

        if not (instrument == 'ahsrl' or instrument == 'nshsrl' or instrument == 'bagohsrl' \
           or instrument == 'mf2hsrl' or instrument =='gvhsrl'):
             print 
             print 'ERROR: ****instrument',instrument, 'does not exist***'
             return
        
        self._time_travel_delta = None
        self._time_travel_str = ""
        self._verbose = 0
        self.profiles=None
        self.figs=None
        # for debugging , simulate that we're running at a previous time
        if os.environ.has_key('TIME_TRAVEL'):
            self.time_travel(os.environ['TIME_TRAVEL'])

#        gt.init_colorbar_status()
        # older versions of matplotlib can crash when re-using windows
       
        try:
           
            self.rti_init( instrument, start_time
                 ,plot_length,min_alt,max_alt,mol_norm_alt
                 ,display,cmd_defaults=cmd_defaults
                 ,netcdf_defaults=netcdf_defaults,z_res=z_res,t_res=t_res)
           
            self.update()
        except RuntimeError, emsg:
            print 'ERROR: Rti.__init__(): %s' % emsg
            return
        except OSError, emsg:
            print 'ERROR: Rti.__init__(): %s' % emsg
            return
        self.display()
    
    def get_netcdf_defaults(self,netcdf_default_file):
        if netcdf_default_file == None:
            netcdf_default_file = self.rs_static.instrument +'_netcdf_defaults.json'
        print 'opening netcdf defaults file ',netcdf_default_file
        fd = open_config(netcdf_default_file)
        dd = json.load(fd)            
        netcdf_defaults=dd
        fd.close()
        self._netcdf_defaults=netcdf_defaults
       
    def get_process_control_defaults(self,cmd_defaults=None,mol_norm_alt=None
                                 ,display=None):
    #def get_process_control_defaults(self,cmd_defaults,mol_norm_alt,ext_pts_to_ave
    #       ,m_smooth,ext_bin_delta,display):

        print 'cmd_def',cmd_defaults
        #get default options for unsupplied cmd line params from process_control.json
        if cmd_defaults == None:
            cmd_defaults = 'process_control.json'
            print 'cd',cmd_defaults
        self.cmd_defaults = cmd_defaults
        print 'cmd_defaults  =', cmd_defaults 
        #if it is a string assume that it is the name of the process cntrl structure.
        if isinstance(cmd_defaults,basestring):
            processing_defaults = jc.json_config(locate_file(cmd_defaults),'process_defaults')
           
       
   
       
        
            
        #if values are supplied on cmd line place them in processing_defaults        
        if mol_norm_alt !=None:
           # processing_defaults['mol_norm_alt']=1000.0*float(mol_norm_alt)
           value = 1000.0*float(mol_norm_alt)
           processing_defaults.set_value('mol_norm_alt','meters',value)
  
  
        self._processing_defaults=processing_defaults
        
        return 
        

    def processing_parameters(self):
        
        return self._processing_defaults
    def corr_adjusts(self):
        return self._corr_adjusts
    
    def verbose(self, i=0):
        self._verbose = i

    def time_travel(self, str=None):
        """ set current time for testing """
        if str == None:
            print "time_travel = %s" % self._time_travel_str
        elif str == 'off':
            self._time_travel_delta = None
            self._time_travel_str = ""
        else:
            self._time_travel_str = str
            print 'WARNING - TIME TRAVEL active - shifting time to %s' % str
            tt_days = hru.convert_date_str(str)['datetime']
            self._time_travel_delta =  datetime.utcnow() - tt_days

    def update(self):
        [self.rs, self.rs_init] = \
            mr.update_cal_and_process(self.rs_static, self.rs_init)
       
        sel_telescope_dir = self.rs_static.processing_defaults.get_value('averaged_profiles','telescope_pointing') 
        #this accumulates. set self.profiles to None to clear it in other functions
        self.profiles=mr.generate_ave_profiles(self.rs.rs_raw,self.rs.rs_mean,self.rs.rs_inv.qc_mask,
            self.rs.rs_Cxx,self.rs_init.rs_constants,self.rs_static.processing_defaults,sel_telescope_dir,
                                               self.rs_static.corr_adjusts,old_profiles=self.profiles)
       
        #compute temperatures if data is available
        if hasattr(self.profiles,'i2a_mol_ratio'):
            self.profiles.i2a_temperatures=i2at.compute_temperature_from_i2a(
                   self.rs_static.instrument,self.profiles.i2a_mol_ratio
                ,self.rs_init.rs_cal.i2a_temp_table,self.rs_init.sounding.pressures)
            self.profiles.i2a_temperatures[self.rs.rs_mean.msl_altitudes
                     <= self.rs_init.rs_constants['lidar_altitude']+200] = np.nan       
        self.rs.profiles=copy.deepcopy(self.profiles)
       
    def cmp_calfile(self,calfile=None):
        """ plot comparison of current default 'baseline','geofile','diff_geofile'
            ,or 'i2_scan' new calfile"""
        print ' '
        print "Supply 'baseline','geofile','diff_geo','i2a_temp_table',or 'i2_scan' filename with complete path "
        self.profiles=None
        
        if calfile==None:
            calfile = raw_input('file name to compare #? ')
        [header,data]=cu.compare_calfiles(self.rs_static.instrument,calfile,self.rs_init.rs_cal)
        print ' '
        print 'new calfile header'
        print header
        print ' '

  
        
        print 'calfile', calfile
        print calfile.find('i2a_temp_table_')

        #replace default calibrations with values from new file.
        if calfile.find('.blc')>0:
            #baseline correction file request
            print 'comparing baseline'
            self.rs_init.rs_cal.baseline.data=data
        elif calfile.find('geofile')>=0 and calfile.find('diff')<0:
            #standard geofile request
            print ' '
            print 'updating both zenith and nadir geo_corr with new data'
            print ' '
            self.rs_init.rs_cal.geo.header=header
            self.rs_init.rs_cal.geo.data=data
            self.rs_init.rs_cal.n_geo.header=header
            self.rs_init.rs_cal.n_geo.data=data
        elif calfile.find('i2a_mol_diff_geofile') >= 0:
            print 'comparing i2a_mol_diff_geofiles'
            self.rs_init.rs_cal.i2a_diff_geo.header = header
            self.rs_init.rs_cal.i2a_diff_geo.data = data
        elif calfile.find('geofile') > 0 and calfile.find('diff') >= 0:
            print 'replacing diff_geo data'
            self.rs_init.rs_cal.diff_geo.data=data
        elif calfile.find('i2-scan') > 0:
            print 'comparing i2_scans'
            self.rs_init.rs_cal.i2scan.data=data
       
        elif calfile.find('i2a_temp_table_')>=0:
            print 'comparing i2a_temp_tables'
            self.rs_init.rs_cal.i2a_temp_table.data = data
            self.rs_init.rs_cal.i2a_temp_table.header = header
            print '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
            print header
        else:
            print 'invalid file calfile name'

        #reprocess last data request with new calibration    
        self.update()
        self.display()

    def reprocess(self):  
     # read new process control defaults
        #use defaults or existing values for all inputs
        print "reading process control defaults"
        self.profiles=None
        self.get_process_control_defaults(self.cmd_defaults,self.rs_static.mol_norm_alt,self.display)
        #self.get_process_control_defaults(self.cmd_defaults,self.rs_static.mol_norm_alt
        #         ,ext_pts_to_ave,m_smooth
        #         ,ext_bin_delta)  #,self.display)
        self.rs_static.processing_defaults = self._processing_defaults
        self.rs_static.corr_adjusts = self._corr_adjusts
       
             
        #try:
        #    self.update()
        #except RuntimeError, emsg:
        #    print 'Rti.new_calvals(): %s', emsg
        #except OSError, emsg:
        #    print 'Rti.new_calvals(): %s', emsg
        #self.display()
        
        #def new_calvals(self):
        """ read new calibration values """

        print 'Reading new calibration values'

        # read new calvals file
        last_i2_scale = self.rs_static.corr_adjusts['i2_corr']
        if self.rs_static.instrument == 'bagohsrl':
           last_i2a_scale = self.rs_static.corr_adjusts['i2a_corr']
        else:
           last_i2a_scale = 1.0 
        last_Cam_scale = self.rs_static.corr_adjusts['Cam_corr']

        # reread calvals file
        calvals = cru.cal_file_reader(self.rs_static.instrument)

        print 'interval_start_time',self.rs_init.interval_start_time
        
        # select calvals for current time
        self.rs_init.rs_constants = cru.select_sys_constants(calvals,
                self.rs_init.interval_start_time)
        #first entry into list of calvals constants employed.
        self.rs_init.rs_constants_used =self.rs_init.rs_constants
               
                       

    

        # if new i2 or Cam correction was supplied, compute new calibration constants

        #note: i2 scan correction is product of values from calvals
        #      corr_adjusts['i2_corr']
        if self.rs_static.corr_adjusts['i2_corr'] != last_i2_scale:
            print 'i2 scaling changed from ', last_i2_scale, ' to ', \
                self.rs_static.corr_adjusts['i2_corr']
        if self.rs_static.corr_adjusts['Cam_corr'] != last_Cam_scale:
            print 'Cam scaling changed from ', last_Cam_scale, ' to ', \
                self.rs_static.corr_adjusts['Cam_corr']        
        if self.rs_init.rs_constants.has_key('i2a_scan_adjustment') \
               and self.rs_static.corr_adjusts['i2a_corr'] != last_i2a_scale:
            print 'i2a scaling changed from ', last_i2a_scale, 'to ',\
                self.rs_static.corr_adjusts['i2a_corr']
            i2a_scan_corr = self.rs_static.corr_adjusts['i2a_corr']\
                    *self.rs_init.rs_constants['i2a_scan_adjustment']
        else:
            i2a_scan_corr = 1.0
            
        self.rs_init.rs_Cxx = cu.quick_cal(
                self.rs_init.rs_cal.i2scan.data,
                self.rs_init.rs_cal.i2scan.Cam
                      * self.rs_static.corr_adjusts['Cam_corr']
                      *self.rs_init.rs_constants['Cam_adjustment'],
                self.rs_init.rs_cal.i2scan.Cam_i2a
                      * self.rs_static.corr_adjusts['Cam_corr']
                      *self.rs_init.rs_constants['Cam_adjustment'],
                self.rs_init.sounding,
                self.rs_init.rs_constants['wavelength'],
                self._processing_defaults.get_value('molecular_spectrum','model'),
                self.rs_static.corr_adjusts['i2_corr']
                *self.rs_init.rs_constants['i2_scan_adjustment']
                ,i2a_scan_corr)
                
        try:
            self.update()
        except RuntimeError, emsg:
            print 'Rti.new_calvals(): %s', emsg
        except OSError, emsg:
            print 'Rti.new_calvals(): %s', emsg
        self.display()

    def print_fig(self,fig_num,printer=None):
        """print requested figure number"""
        print 'printing figure # = ',fig_num
        f=plt.figure(fig_num)
        f.savefig('myfilename.pdf', format='pdf')
       
        if printer==None:
            printer = self.rs_static.display_defaults.get_value(
                'printer_id','name',key='config')
        printer = '-d' + printer[:]
        call(['/usr/bin/lp', printer, 'myfilename.pdf'])
        
    def save_fig(self,fig_num,filename=None,dots_per_inch=None):
        """save requested figure number"""

        print 'saving figure # =',fig_num
        if filename==None:
            filename = raw_input('file name (e.g. path/figame.png) = ? ')
        if dots_per_inch==None:
            dots_per_inch = np.float(raw_input('resolution (dpi) =? '))
        f=plt.figure(fig_num)
        f.savefig(filename,dpi = dots_per_inch)
    
    def write_netcdf(self,tag=None,netcdf_format=None,output_dir=None):
        
        start_time_str = \
              self.rs.rs_raw.times[0].strftime("%Y%m%dT%H%M")
        end_time_str = \
            self.rs.rs_raw.times[-1].strftime("%Y%m%dT%H%M")
        if tag==None:
            tag = raw_input('file name tag?  ')
        if tag == '':
            filename=self.rs_static.instrument +'_'+start_time_str \
               +'_' + end_time_str +'.nc'
        else:
            filename=self.rs_static.instrument +'_'+start_time_str \
               +'_' + end_time_str+'_' +tag+'.nc'
    
        if netcdf_format==None:
            netcdf_format = self.rs_static.display_defaults.get_value('netcdf','format',key='config')
        if output_dir==None:
            output_dir =self.rs_static.display_defaults.get_value('netcdf','output_dir',key='config')
        filename = os.path.join(output_dir,filename)
        print '      ',filename
        if netcdf_format in ['uwlidar','uwlidar3','uwlidar3raw','hsrl_nomenclature.cdl','hsrl3_processed.cdl','hsrl3_raw.cdl'] :#FIXME 20121207 JPG why are there 2 classes? this is supposed to be 1 exporter, configured by the CDL/initialization!
            import lg_dpl_toolbox.dpl.dpl_create_templatenetcdf as dpl_ctnc
            print 'writing ' + netcdf_format + ' format netcdf '
            template='hsrl_nomenclature.cdl'
            if netcdf_format=='uwlidar3':
                template='hsrl3_processed.cdl'
            elif netcdf_format=='uwlidar3raw':
                template='hsrl3_raw.cdl'
            elif netcdf_format.endswith('cdl'):
                template=netcdf_format
            format="NETCDF4"
            if '3' in template:
                format='NETCDF3_CLASSIC'
            n=Dataset(filename,'w',clobber=True,format=format)
            print 'writing template netcdf with template ',template
            v=dpl_ctnc.dpl_create_templatenetcdf(locate_file(template),n,self.rs)
            v.appendtemplatedata(self.rs)
        elif 'cfradial' in netcdf_format:
            if netcdf_format.endswith('cdl'):
                template=netcdf_format
            else:
                template='hsrl_cfradial.cdl'
            print 'writing cfradial format netcdf to location ',output_dir
            import lg_dpl_toolbox.dpl.dpl_create_cfradial as dpl_ctnc
            n=Dataset(filename,'w',clobber=True)
            v=dpl_ctnc.DplCreateCfradial(locate_file(template),n,self.rs)
            v.append_data(self.rs)
            v.close()
        n.sync()
        n.close()

    def read_processed_hsrl_data(self,netcdf_filename=None):
        if netcdf_filename==None:
            netcdf_filename = raw_input('path/netcdf_filename?    ')
        hhru.read_processed_hsrl_data(netcdf_filename)

    def plot_model(self):
        """plot profile computed from model specs"""

        pm.performance_model(
            self.rs_init.rs_Cxx,
            self.rs_init.rs_cal,
            self.rs.profiles.dc_combined_hi_counts,
            self.rs.rs_raw.transmitted_energy,
            self.rs.rs_raw.seeded_shots,
            self.rs.rs_raw.times,
            self.rs_init.rs_constants,
            )

    def cal_scale(self):
        """ rescale or turn off calibrations """
    
        self.profiles=None
        last_i2_scale = self.rs_static.corr_adjusts['i2_corr']
        last_i2a_scale = self.rs_static.corr_adjusts['i2a_corr']
        last_Cam_scale = self.rs_static.corr_adjusts['Cam_corr']
       
        self.rs_static.corr_adjusts = \
            get_scale_corrections(self.rs_static.corr_adjusts,self.rs_static.instrument)

    # if new i2 correction was supplied, get new calibration constants

        
        if self.rs_static.corr_adjusts['i2_corr'] != last_i2_scale \
               or self.rs_static.corr_adjusts['Cam_corr'] != last_Cam_scale\
               or self.rs_static.corr_adjusts['i2a_corr'] != last_i2a_scale :
            if self.rs_init.rs_constants.has_key('i2a_scan_adjustment'):
                print 'i2a scaling changed from ', last_i2a_scale, 'to ',\
                      self.rs_static.corr_adjusts['i2a_corr']
                i2a_scan_corr = self.rs_static.corr_adjusts['i2a_corr']\
                    *self.rs_init.rs_constants['i2a_scan_adjustment']
            else:
                i2a_scan_corr = 1.0
           
            self.rs_init.rs_Cxx = cu.quick_cal(self.rs_init.rs_cal.i2scan.data
                ,self.rs_init.rs_cal.i2scan.Cam
                     * self.rs_static.corr_adjusts['Cam_corr']
                     * self.rs_init.rs_constants['Cam_adjustment']
                ,self.rs_init.rs_cal.i2scan.Cam_i2a
                     * self.rs_static.corr_adjusts['Cam_corr']
                     * self.rs_init.rs_constants['Cam_adjustment']
                ,self.rs_init.sounding
                , self.rs_init.rs_constants['wavelength']
                ,self._processing_defaults.get_value('molecular_spectrum','model')
                ,self.rs_static.corr_adjusts['i2_corr']
                    *self.rs_init.rs_constants['i2_scan_adjustment']
                ,i2a_scan_corr)           
        try:
            self.update()       #(inCalibration=True)
        except RuntimeError, emsg:
            print 'Rti.cal_scale(): %s', emsg
        except OSError, emsg:
            print 'Rti.cal_scale(): %s', emsg
        self.display()




    
    def cal_gen(self,cal_type):
        self.profiles=None
        if cal_type == 'geo':
            if self.rs_static.instrument == 'gvhsrl':
                cfg.make_geofile_new(
                       self.rs_static.instrument
                      ,self.rs.profiles
                      ,self.rs_init.rs_cal
                      ,self.rs.rs_Cxx
                      ,self.rs_init.rs_constants
                      ,self.rs_static.corr_adjusts           
                      ,self.rs_static.processing_defaults)
                return    
            else:
                #reinitialize with altitudes covering full range of lidar
                # and 7.5 m resolution
                min_alt = self.rs_init.rs_constants['lidar_altitude']/1000.0
                max_alt = 30.0+min_alt
                mol_norm_alt = min_alt + 0.2
                display='all_plots.json'
                self.rti_init( self.rs_static.instrument, self.start_time
                   ,self.plot_length,min_alt,max_alt,mol_norm_alt
                   ,display,self.cmd_defaults
                   ,self._netcdf_defaults,z_res={'manual':7.5},t_res=None)

                #reprocess data with new resolution and no geo_corr
                self.rs_static.corr_adjusts['geo_corr'] = 0.0
                if self.rs_static.processing_defaults.get_value('averaged_profiles','apply_mask'):
                    print 
                    print '         you must turn off averaged_profiles in process_control.json'
                    print "         before running r.cal_gen('geo')"
                    print 
                    return
                self.update()

                #use reprocessed data to compute geo_corr
                cfg.make_geofile(
                       self.rs_static.instrument
                      ,self.rs.profiles
                      ,self.rs_init.rs_cal
                      ,self.rs.rs_Cxx
                      ,self.rs_init.rs_constants)
                return
        if cal_type == 'wfov_geo':
            cfg.make_wfov_geofile(
                       self.rs_static.instrument
                      ,self.rs.rs_raw
                      ,self.rs.profiles
                      ,self.rs_init.rs_cal
                      ,self.rs.rs_Cxx
                      ,self.rs_init.rs_constants
                      ,self.rs_static.processing_defaults)
            return
        elif cal_type == 'baseline':
            cfg.make_baseline_file(
                       self.rs_static.instrument
                      ,self.rs.rs_raw
                      ,self.rs_init.rs_constants
                      ,self.rs_static.corr_adjusts)
            return
        elif cal_type == 'i2scan':
            cfg.make_i2_scan_file(
                       self.rs_static.instrument
                      ,self.rs
                      ,self.rs_init.rs_constants)
            return
        elif cal_type == 'i2a_temp_table':
            if self.rs_static.instrument == 'bagohsrl':
                cfg.make_temperature_table(
                      self.rs_static.instrument
                      ,self.rs_init.rs_cal
                      ,self.rs.profiles.times
                      ,self.rs_init.rs_constants
                      ,self.rs_static.corr_adjusts
                      ,self.rs_static.processing_defaults)
                return
            else:
                print 'no temp measurement channel on ',self.rs_static.instrument
                return
        elif cal_type == 'i2a_diff_geo':
            cfg.make_i2a_diff_geofile(self.rs_static.instrument
                                      ,self.rs.rs_raw
                                      ,self.rs_init.rs_cal
                                      ,self._processing_defaults    
                                      ,self.rs_init.rs_constants
                                      ,self.rs_static.corr_adjusts)
        
        elif self.rs_static.instrument == 'gvhsrl' and cal_type == 'diff_geo':
            
            try:
                if self.rs.profiles.times.shape[0] < 1:
                    raise RuntimeError, "cal_gen - not enough times in profile"
            except IndexError:
                raise RuntimeError, "cal_gen - no profile exists - select data with 'r.params()'"

            # generate results at native altitude resolution without geo correction
            # only processses the last data chunk--does not work across cal or raob changes
            
            self.rs_init.rs_constants['first_bin_to_process'] = 0
            self.rs_static.mol_norm_alt = 10000
            self.rs_static.min_alt = 0
            self.rs_static.max_alt = 60000
            self.rs_static.alt_res = 7.5
            self.rs_static.n_range_ave = 1
            self.rs_static.corr_adjusts['geo_corr'] = 0
            self.rs_static.corr_adjusts['i2a_dif_geo_corr']=0
      
            #******************************************************
            #get soundings at native resolution
            requested_altitudes=np.arange(0,self.rs_static.max_alt,self.rs_static.alt_res)
            self.rs_init.rs_soundings = su.sounding_archive(
                 self.rs_static.instrument,
                 self.rs_init.rs_constants['sounding_type'],
                 self.rs_init.rs_constants['sounding_id'],
                 self.rs_init.interval_start_time,requested_altitudes)
        
            #get current sounding
            self.rs_init.sounding = self.rs_init.rs_soundings.profile(self.rs_init.interval_start_time,[],[])
            if self.rs_init.rs_constants.has_key('i2a_scan_adjustment'):
               i2a_scan_corr = self.rs_static.corr_adjusts['i2a_corr']\
                    *self.rs_init.rs_constants['i2a_scan_adjustment']
            else:
               i2a_scan_corr = 1.0
  
            #compute new calibrations at native resolution 
            self.rs_init.rs_Cxx = cu.quick_cal(self.rs_init.rs_cal.i2scan.data
                ,self.rs_init.rs_cal.i2scan.Cam
                    *self.rs_static.corr_adjusts['Cam_corr']
                    *self.rs_init.rs_constants['Cam_adjustment']
                ,self.rs_init.rs_cal.i2scan.Cam_i2a
                    *self.rs_static.corr_adjusts['Cam_corr']
                    *self.rs_init.rs_constants['Cam_adjustment']
                ,self.rs_init.sounding
                ,self.rs_init.rs_constants['wavelength']
                ,self._processing_defaults.get_value('molecular_spectrum','model')
                ,self.rs_static.corr_adjusts['i2_corr']
                     * self.rs_init.rs_constants['i2_scan_adjustment']
                ,i2a_scan_corr) 

            if self.rs_static.instrument == 'bagohsrl':
                self.rs_static.corr_adjusts['i2a_dif_geo_corr'] = 0.0
        
        
            #try:
            #    self.update()
            #except RuntimeError, emsg:
            #    print 'Rti.cal_gen(): %s', emsg
            #except OSError, emsg:
            #    print 'Rti.cal_gen(): %s', emsg
       
            if self.rs_init.rs_constants.has_key('installation') \
                   and self.rs_init.rs_constants['installation'] == 'airborne':
                #select telescope pointing directions to use in generating new cal files    
                sel_telescope_dir=raw_input('select telescope pointing as zenith or nadir ?  ')
            else:
                sel_telescope_dir = 'zenith'
            
            #form average profiles selecting desired telescope pointing direction
            self.rs.profiles=mr.generate_ave_profiles(self.rs.rs_raw,self.rs.rs_mean,self.rs.rs_inv.qc_mask
                                ,self.rs.rs_Cxx,self.rs_init.rs_constants
                                 ,self._processing_defaults
                                 ,sel_telescope_dir,self.rs_static.corr_adjusts)
       
            
        elif cal_type == 'diff_geo':
                cfg.make_diff_geofile(self.rs_static.instrument
                                     ,self.rs.rs_raw
                                     ,self.rs_init.rs_cal
                                     ,self.rs_static.corr_adjusts 
                                     ,self.rs_init.rs_constants)
              
                self.display()
        else:  
            print
            print "Usage: r.cal_gen(type)"
            print "where type = ('geo' |'wfov_geo' |'baseline' | 'diff_geo'| i2a_diff_geo' | 'i2scan' | 'i2a_temp_table')"
            print "cover telescope for 'baseline' correction"
            print "Use clear air data for 'geo' correction"
            print "Remove I2 cell for 'diff_geo' correction"
            print "Remove Argon buffered I2 cell for 'i2a_diff_geo' correction"
  
        return

    def replot(self,display=None):
        """ replot(display=display_defaults.json')
            clear figures, reload display defaults file, and replot last data
            display_defaults_file = optional--new display defaults file name
                                  = .json file containing display directives"""
 
        self.clear()
        [self.rs_static.display_defaults, config]= \
           du.get_display_defaults(self.rs_static.display_defaults_file,"new")
        self.display()

    def clear(self):
        """clear all figures"""

    # loop through all existing figures

        for x in self.figs:#plt._pylab_helpers.Gcf.get_all_fig_managers():
            f=self.figs.figure(x)#.num
            f.clear()
            #gt.del_colorbar(f.number)
            if hasattr(f,'fig_colorbars'):
                delattr(f,'fig_colorbars')
        self.repaint()

    def close(self):
        """ close all figures """
        if self.figs!=None:
          if True:
              self.figs.close()
          else:
            for x in self.figs:
              f=self.figs.figure(x)
              f.clf()
              del f
        

    def repaint(self):
        firstfig=None
        plt.show(block=False)

        for x in self.figs:
            fig = self.figs.figure(x)
            #help(fig.canvas)
            #fig.set_size_inches(fig.get_size_inches(),forward=True)
            i=fig.get_size_inches()*fig.get_dpi();
            fig.canvas.resize(i[0],i[1])
            if hasattr(fig.canvas,'repaint'):
                fig.canvas.repaint()
            else:
                fig.canvas.draw()
            if firstfig==None or firstfig.number>fig.number:
                firstfig=fig
        if firstfig!=None:
            #firstfig.set_size_inches(firstfig.get_size_inches(),forward=True)
            i=firstfig.get_size_inches()*firstfig.get_dpi()
            firstfig.canvas.resize(i[0],i[1])
            if hasattr(firstfig.canvas,'repaint'):
                firstfig.canvas.repaint()
            else:
                firstfig.canvas.draw()

    def display(self):
        isint=plt.isinteractive()
        prefint=True
        if prefint != isint:
            if prefint:
                plt.ion()
            else:
                plt.ioff()
        self.figs=du.show_images(
            self.rs_static.instrument,
            self.rs,
            self.rs_init.sounding,
            self.rs_init.rs_constants,
            self.rs_init.rs_cal.geo.data,
            self.rs_static.processing_defaults,
            self.rs_static.display_defaults,
            self.rs_init.last_sounding_time,
            self.rs_static.max_alt,
            self.rs_static.auto_loop,
            self.figs
            )

        self.repaint()
            
        if isint != plt.isinteractive():
            if isint:
                plt.ion()
            else:
                plt.ioff()
 
        print '\n'
        print """cmds-->
        r.next()       r.cnext()  r.cal_scale()    r.cal_gen('cal_type') 
        r.clear()      r.loop()   r.print_fig(#)   r.save_fig(#) 
        r.replot()     r.close()  r.cmp_calfile()  r.show_params()
        r.loop()       r.help()   r.reprocess()    r.write_netcdf()     
        r.plot_model()"""
    def help(self):
        print '\n'
        print 'r.next()-------- Plot next time interval on top of last plot.'
        print 'r.cnext()------- Clear plot and plot next time interval.'
        print 'r.clear()------- Clear old plots.'
        print 'r.close()------- Close all plots.'
        print 'r.replot()------ Replot last time interval.'
        print 'r.replot(display="display_defaults.json") for different plots of last data'
        
        print 'r.reprocess(0--- Reloads calvals_xxhsrl.txt and process_control.json',
        print '                 then recomputes last request with new profiles plotted.'
        print '                 along with old.'
        print 'r.cal_gen()----- Generate baseline, geofile, or diff_geofile cal files.'
        print '                 Geofiles can be computed from surface or constant altitude'
        print '                 flight legs with zenith pointing telescope.'
        print 'r.cal_scale()--- Scale or disable various data corrections.'
        print 'r.cmp_calfile()- Plot a comparison of a new baseline, geo, diff_geo, or i2_scan' 
        print '                 file with current default calibration values and recompute'
        print '                 the last request with the new calibration'
        print 'r.loop()-------- Automatically advance through data.'
       
        print 'r.print_fig(#)-- Print fig on printer designated in "display_defaults.json."'
        print 'r.save_fig(#)--- Save figure as graphics file'
        print 'r.plot_model()-- Plot signals projected from model specs provided in'
        print "                 'hsrl_config/xxhsrl_specs.json'."
        print 'r.show_params()--Show r.params call used to start current process.'
        print """r.params(start_time, plot_length, min_alt,
                    max_alt, mol_norm_alt[, display_defaults_file])"""
        print 'r.write_netcdf()-create netcdfile_tag.nc as per display_display.json directives.'   
        print ' '
        print 'r.help()------- Help.'

    def loop(self,maxcount=-1):
        done = False
        while not done:
            fig_numbers = [x for x in self.figs]
            #                           plt._pylab_helpers.Gcf.get_all_fig_managers()]

            for i in fig_numbers:#range(len(fig_numbers)):
                f=self.figs.figure(i)
                f.clf()
                #gt.del_colorbar(f.number)
                if hasattr(f,'fig_colorbars'):
                    delattr(f,'fig_colorbars')
            try:
                self._next()
            except RuntimeError, exc:
                print exc
                done = True
            except OSError, exc:
                print exc
                done = True
            # some log plots get np.ma.MaskErrors - known bug, no simple fix
            except np.ma.MaskError, exc:
                print exc
                done = True
            if maxcount>0:
                maxcount-=1
            if maxcount==0:
                done=True
    
    def next(self, repeat=0):
        """ plot the next time interval, reporting any exceptions"""
        try:
            self._next(repeat)
        except RuntimeError, exc:
            print exc
        except OSError, exc:
            print exc
        except np.ma.MaskError, exc:
            print exc


    def _next(self, repeat=0):
        """ plot next data """

        now = datetime.utcnow()
        warning = ''
        if self._time_travel_delta != None:
          now=now-self._time_travel_delta
          warning = '(TIME_TRAVEL!)'
        if self._verbose: 
            print 'Current time is= ', now, warning, \
            " plot length = %.2f hours "% (self.rs_static.delta_t * 24.0)

        # truncate the end of the interval to 'now' 
        # (there is not any data after 'now')
        start_adj = ''
        end_adj = ''
        if self.rs_init.interval_end_time > now:
            end_adj = '**'
            start_adj = '**'
            self.rs_init.interval_end_time = now
            self.rs_init.interval_start_time = now \
                - self.rs_static.delta_t
            if self._verbose: print "_next :  [now-delta, now]"
        elif repeat == 0:
            # move start and end times by delta
            if self._verbose: print '_next: [start + delta, end + delta] '
            self.rs_init.interval_start_time = \
                self.rs_init.interval_start_time \
                + self.rs_static.delta_t
            self.rs_init.interval_end_time = \
                self.rs_init.interval_start_time \
                + self.rs_static.delta_t
        if self.rs_init.interval_start_time > now:
            # interval stops at now
            if self._verbose: print "_next : [now-delta, now] " 
            end_adj = '**'
            start_adj = '**'
            self.rs_init.interval_start_time = now \
                - self.rs_static.delta_t
            self.rs_init.interval_end_time = now
        elif self.rs_init.interval_end_time > now:
            end_adj = '**'
            start_adj = '**'
            # interval ends at now
            if self._verbose: print "_next : [now-delta, now]"
            self.rs_init.interval_end_time = now
            self.rs_init.interval_start_time = \
                self.rs_init.interval_end_time - self.rs_static.delta_t
        if self._verbose:
            print '_next: plot [%s' %start_adj,  
            print self.rs_init.interval_start_time, ',', 
            print self.rs_init.interval_end_time, '%s]'% end_adj

    # don't catch the RuntimeError exception, so Rti.loop() will see it
    
        #bin_dt is the instrument bin time durration
        bin_dt = float(self.rs_init.rs_constants['integration_time'])
        if self.rs_init.rs_constants['quarter_wave_plate_rotation'] == 'fixed' \
                 or self.rs_init.rs_constants[
                 'quarter_wave_plate_rotation']=='none':
            frame_dt = 0.0
        else:
            frame_dt = float(self.rs_init.rs_constants['polarization_integration_time'])
        import lg_base.core.canvas_info as ci
        canvas_info=ci.load_canvas_info()
        number_x_pixels = canvas_info['canvas_pixels']['x']
        [time_res,ntime_ave] = sr.set_time_resolution(self.rs_static.delta_t, number_x_pixels
                        ,frame_dt,bin_dt, **({} if not self.rs_static.t_res else self.rs_static.t_res))
      

        requested_times = np.array( [(self.rs_init.interval_start_time + timedelta(seconds=x))\
                for x in np.arange(time_res.total_seconds()/2.0,self.rs_static.delta_t.total_seconds()\
                +time_res.total_seconds()/2.0,time_res.total_seconds()) ])
 
        self.rs_static.requested_times = requested_times
        self.rs_static.time_res = time_res 
        self.update()     #(inCalibration=True)
        self.display()

    def cnext(self, repeat=0):
        """ plot next data wihout overplotting graphs """

        self.clear()
        self.next(repeat)

    def rti_init( self, instrument, start_time, plot_length, min_alt
         ,max_alt,mol_norm_alt=None,display=None,cmd_defaults=None
         ,netcdf_defaults=None,z_res=None,t_res=None):
    #def rti_init( self, instrument, start_time, plot_length, min_alt
    #     ,max_alt,mol_norm_alt=None,ext_pts_to_ave=None,m_smooth=None
    #     ,ext_bin_delta=None,display=None,cmd_defaults=None
    #     ,netcdf_defaults=None,z_res=None,t_res=None):
        """ rti_init(instrument,start_time,plot_length,min_alt,max_alt,mol_norm_alt\
            display_defaults_file)
            Generates rti images of hsrl data.
            start_time          = time string, first data to plot (eg.'5-may-11 12:30').
            instrument          = hsrl id string (eg. 'ahsrl','gvhsrl','nshsrl','mf2hsrl').
            plot_length         = time period to plot in hours (e.g. 2) this may also be
                                  entered as end time of first interval(e.g. '6-may-10 1:30)
                                  or as 1:30 if end time is on same day as start time.
            min_alt             = minimum altitude (msl) to display.
            max_alt             = maximum altitude (msl) to display.
            mol_norm_alt         = normalization altitude for print mol_norm_alt
                                  if not supplied mol_norm_alt is taken from
                                  process_control.json
           display_defaults_file= name of .json file containing display directives
                                   """
        #get initial value of corr_adjust from process_control json
        self.get_process_control_defaults(cmd_defaults)
        self._corr_adjusts =self._processing_defaults.get_dict('corr_adjusts')
       
       
        # input parameters that may be changed by "update_cal_and_process"
        self.rs_init = hau.rs_xfer()
        self.rs_init.rs_constants_used=[]

        # input parameters that are not changed by "update_cal_and_process"
        self.rs_static = hau.rs_xfer()

        # save instrument name 
        self.rs_static.instrument = instrument
       

        self.__params__(start_time, plot_length, min_alt
            ,max_alt,mol_norm_alt,display,cmd_defaults
            ,z_res=z_res,t_res=t_res)
        #self.__params__(start_time, plot_length, min_alt
        #    ,max_alt,mol_norm_alt,ext_pts_to_ave,m_smooth,ext_bin_delta,display
        #,cmd_defaults,z_res=z_res,t_res=t_res)       
    
    def show_params(self):
        """ show last self.params() call"""
        
        print ' '
        print "r.params('%s', %s, %4.2f, %4.2f" % (self.start_time, repr(self.plot_length)
            ,self.min_alt, self.max_alt),
        for item in self.rs_static.processing_defaults:
            if item == 'display':
               print ", %s='%s'" % (item,self.rs_static.processing_defaults[item]),
            elif item == 'mol_norm_alt':
               print ", %s=%s" % (item,self.rs_static.processing_defaults[item]/1000.0),
            #don't print these items with r.show_params()   
            elif item == 'image_total_xy_pixels' \
                     or item == 'do_ave_profiles' \
                     or item == 'spectral_model':
                a=1
            elif item[0] != '#':
               print ", %s=%s" % (item,self.rs_static.processing_defaults[item]),
        print ')'
        print ' '
        
        
      

    def params(self, start_time, plot_length, min_alt, max_alt
           , mol_norm_alt=None,display = None
           ,z_res = None,t_res=None):
    #def params(self, start_time, plot_length, min_alt
    #       , max_alt, mol_norm_alt=None,ext_pts_to_ave=None,m_smooth=None
    #       ,ext_bin_delta=None,display = None,z_res = None,t_res=None):
        
        """ params(instrument,start_time,plot_length,min_alt,max_alt,mol_norm_alt\
            display_defaults_file)
            Generates rti images of hsrl data.
            start_time          = time string, first data to plot (eg.'5-may-11 12:30').
            plot_length         = time period to plot in hours (e.g. 2) this may also be
                                  entered as end time of first interval(e.g. '6-may-10 1:30)
                                  or as 1:30 if end time is on same day as start time.
            min_alt             = minimum altitude (msl) to display.
            max_alt             = maximum altitude (msl) to display.
            mol_norm_alt        = normalization altitude for optical depth (km).
                                  if not supplied mol_norm_alt is taken from
                                  calvals_xhsrl.txt
           display_defaults_file= name of .json file containing display directives
           """
       
        if display == None:
            display = self.rs_static.display_defaults_file
        if mol_norm_alt == None:
            mol_norm_alt = self.rs_static.processing_defaults['mol_norm_alt']/1000.0
        #if m_smooth == None:
        #    m_smooth = self.rs_static.processing_defaults['m_smooth']
        #if ext_pts_to_ave == None:
        #   ext_pts_to_ave = self.rs_static.processing_defaults['ext_pts_to_ave']
        #if ext_bin_delta == None:
        #   ext_bin_delta = self.rs_static.processing_defaults['ext_bin_delta']
      
    
        self.__params__(start_time, plot_length, min_alt, max_alt,mol_norm_alt
                 ,display,self.cmd_defaults,z_res=z_res,t_res=t_res)
        #self.__params__(start_time, plot_length, min_alt
        #     ,max_alt,mol_norm_alt, ext_pts_to_ave,m_smooth,ext_bin_delta
        #     ,display,self.cmd_defaults,z_res=z_res,t_res=t_res)
        try:
            self.update()
            self.display()
            #self.replot()
        except RuntimeError, exc:
            print exc
            done = True
        except OSError, exc:
            print exc

    def __params__(self, start_time, plot_length, min_alt,max_alt
           ,mol_norm_alt=None,display=None
           ,cmd_defaults=None,netcdf_defaults=None,z_res=None,t_res=None):

    #def __params__(self, start_time, plot_length, min_alt
    #         ,max_alt,mol_norm_alt=None,ext_pts_to_ave=None,m_smooth=None,ext_bin_delta=None
    #         ,display=None,cmd_defaults=None,netcdf_defaults=None,z_res=None,t_res=None):    

        
        # save all user params so we can show their last values
        self.start_time = start_time
        self.plot_length = plot_length
        self.min_alt = min_alt
        self.max_alt = max_alt
    
        #get processing defaults
        self.get_process_control_defaults(cmd_defaults,mol_norm_alt,display)
        #self.get_process_control_defaults(cmd_defaults,mol_norm_alt,ext_pts_to_ave
        #                    ,m_smooth,ext_bin_delta,display) 
        self.rs_static.processing_defaults = self._processing_defaults
        self.rs_static.corr_adjusts = self._corr_adjusts

        #get netcdf_defaults list of requested variables to read from netcdf
        self.get_netcdf_defaults(netcdf_defaults)
        self.rs_static.netcdf_defaults = self._netcdf_defaults
        self.rs_static.t_res = t_res
        self.rs_static.z_res = z_res
        
        #get display_defaults        import display_utilities as du
        #[display_defaults, config]= \
        #   du.get_display_defaults(self.rs_static.processing_defaults['display'],"new") 
        
        [display_defaults, config]= du.get_display_defaults(display,"new") 
       
       
        print '   printer selection = ', config.get_value('printer_id','name',key='config')

        # convert alt from km to meters

        min_alt = 1000 * float(min_alt)
        max_alt = 1000 * float(max_alt)
       
        # check if housekeeping data plots are requested
        # if so it will be necessary to read the housekeeping data

        data_request = 'images'
        for (name, value) in display_defaults.items():

           # these displays require housekeeping data

            if name.find('show') >= 0 and value[0] > 0:
                data_request = 'images housekeeping'
                print 'requesting housekeeping data'
                break

       
       
        max_range_bin = 4000
        #else:
        #    max_range_bin = np.int(max_alt / 7.5)
        #    if max_range_bin > 4000:
        #        max_range_bin = 4000
        #    elif max_range_bin < 2000:
        #        max_range_bin = 2000

        # convert start_time to python datetime

        interval_start_time_str = start_time
        rs_date = hru.convert_date_str(start_time)
        interval_start_time = rs_date['datetime']

        # if plot_length entered as end time string

        if isinstance(plot_length, basestring):

            # check for time without date--assume it is same day as start

            if plot_length.find(':') < 3:
                index = start_time.find(' ')
                plot_length = start_time[0:index + 1] + plot_length
            rs_date = hru.convert_date_str(plot_length)
            end_time = rs_date['datetime']
            delta_t = end_time - interval_start_time
            if delta_t <= timedelta(seconds=0):
                print ' '
                print 'ERROR----end time must be later than start time'
                print ' '
                raise RuntimeError, 'ERROR----end time must be later than start time'
        else:
            delta_t = timedelta(days=(plot_length) / 24.0)  # plot length in hours converted to timedelta



        rs_mem = hau.Time_Z_Group()
        rs_Cxx = hau.Time_Z_Group()# cu.Cxx()

        # loop setup

        interval_start_time_str = start_time
        rs_date = hru.convert_date_str(start_time)
        interval_start_time = rs_date['datetime']
        current_file_end_time = interval_start_time

        # read calvals file of system calibration constants
        # and select current value of constants

        self.rs_init.calvals = cu.calval_info(self.rs_static.instrument)
        self.rs_init.rs_constants = self.rs_init.calvals.select_time(interval_start_time)
        self.rs_init.rs_constants_used.append(self.rs_init.rs_constants)



        #set altitude resolution
        import lg_base.core.canvas_info as ci
        canvas_info=ci.load_canvas_info()
        number_y_pixels = canvas_info['canvas_pixels']['y']
        #binwidth in meters
        binwidth = self.rs_init.rs_constants['binwidth']*1.5e8
        print 'y pixels',number_y_pixels
        [alt_res,n_range_ave] = sr.set_z_resolution( 
               min_alt=min_alt, max_alt=max_alt, binwidth=binwidth
            ,number_y_pixels=number_y_pixels, **({} if not z_res else z_res))        
        
        # set time resolution
        
        #bin_dt is the instrument bin time durration
        bin_dt = float(self.rs_init.rs_constants['integration_time'])
        if self.rs_init.rs_constants['quarter_wave_plate_rotation'] == 'fixed' \
                 or self.rs_init.rs_constants[
                 'quarter_wave_plate_rotation']=='none':
            frame_dt = 0.0
        else:
            frame_dt = float(self.rs_init.rs_constants['polarization_integration_time'])
        number_x_pixels = canvas_info['canvas_pixels']['x']
        [time_res,ntime_ave] = sr.set_time_resolution(delta_t, number_x_pixels
                        ,frame_dt,bin_dt, **({} if not t_res else t_res))
        # Matt Edit:  Avoid preavaraging when QWP is rotating
        if self.rs_init.rs_constants['quarter_wave_plate_rotation'] == \
                'rotating':
            ntime_ave = 1

        self.rs_static.ntime_ave = ntime_ave
    
        print 
        print 'time_resolution = ',  time_res.total_seconds(), ' seconds'
        print 'number of time bins in preaverage = ', ntime_ave
        print

       
        requested_times = np.array( [(interval_start_time + timedelta(seconds=x))\
                for x in np.arange(time_res.total_seconds()/2.0,delta_t.total_seconds()\
                +time_res.total_seconds()/2.0,time_res.total_seconds()) ])

        #note 0.5 factor to prevent missing pixels
        #self.rs_static.ntime_ave = max(int(0.5*time_res.total_seconds()/integration_time),1)

        #self.rs_static.processing_defaults['do_ave_profiles'] = \
        #       display_defaults.enabled('backscat_profile') \
        #       or display_defaults.enabled('depol_profile') \
        #       or display_defaults.enabled('od_profile')
  
        
 
      #if display_defaults.enabled('raw_profiles') \
        #    or display_defaults.enabled('dif_geo_profiles'):
        #    self.rs_static.processing_defaults = 2
            

        # read initial range dependent cal_vectors, baseline, geo_corr, diff_geo_corr
        # and the i2_scan_file.

        # read initial values of the calibration vectors, baseline,geo,diff_geo, i2scan
        #if self.rs_static.processing_defaults.has_key('alternate_cal_dir'):
        #    alternate_cal_dir = self.rs_static.processing_defaults['alternate_cal_dir']
        #else:
        #    alternate_cal_dir = None

        alternate_cal_dir = self.rs_static.processing_defaults.get_value('alternate_cal_dir','full_dir_path') 
        if alternate_cal_dir == 'None':
             alternate_cal_dir = None          

        rs_cal = mr.cal_vectors(self.rs_static.instrument,interval_start_time,
                             max_range_bin,alternate_cal_dir)

        
       

        # get soundings file extended to 50 km if needed by climatology
        # soundings include one prior to start time to end of file
        # values returned at alt resolution given by 7.5*n_range_ave
     
        self.rs_static.requested_altitudes=np.arange(0,50000,7.5*n_range_ave)    
        rs_soundings = su.sounding_archive(
             self.rs_static.instrument,
             self.rs_init.rs_constants['sounding_type'],
             self.rs_init.rs_constants['sounding_id'],
             interval_start_time,
             self.rs_static.requested_altitudes)

        

        sounding = rs_soundings.profile(interval_start_time,[],[])    
       
        
                
        sounding_date = sounding.times
    
        print 'Initial sounding= ' + sounding.station_id[:] + '  ' \
            + sounding_date.strftime('%d-%b-%y %H:%M')
    
        last_sounding_time = 0
        
       
    
        # inversion coef are generated at the requested altitude resolution up to 50km
        if rs_cal.i2scan.Cam_i2a == None:
            rs_cal.i2scan.Cam_i2a = 0.0

        if self.rs_init.rs_constants.has_key('i2a_scan_adjustment'):
            i2a_scan_corr = self.rs_static.corr_adjusts['i2a_corr']\
                    *self.rs_init.rs_constants['i2a_scan_adjustment']
        else:
            i2a_scan_corr = 1.0
            
        rs_Cxx = cu.quick_cal(rs_cal.i2scan.data
            ,rs_cal.i2scan.Cam
                  * self.rs_static.corr_adjusts['Cam_corr']
                  * self.rs_init.rs_constants['Cam_adjustment']
            ,rs_cal.i2scan.Cam_i2a
                  * self.rs_static.corr_adjusts['Cam_corr']
                  * self.rs_init.rs_constants['Cam_adjustment']
            ,sounding
            ,self.rs_init.rs_constants['wavelength']
            ,self._processing_defaults.get_value('molecular_spectrum','model')
            ,self.rs_static.corr_adjusts['i2_corr']
            *self.rs_init.rs_constants['i2_scan_adjustment']
            ,i2a_scan_corr)
  
       

        # no data currently in memory

        rs_data_end_time = None
        interval_end_time = interval_start_time + delta_t

        auto_loop = 0  # wait for operator input before next set of plots

        #are photon counting stats requested
        #compute_stats = self._processing_defaults['compute_stats'][0] 
 
        compute_stats =self._processing_defaults.enabled('compute_stats')


        # package outputs into structures


        # input parameters that are not changed by "update_cal_and_process"

        self.rs_static.start_time = start_time
        self.rs_static.plot_length = plot_length
        self.rs_static.min_alt = min_alt
        self.rs_static.max_alt = max_alt
        self.rs_static.mol_norm_alt = mol_norm_alt
        self.rs_static.processing_defaults = self._processing_defaults
        self.rs_static.corr_adjusts = self._corr_adjusts      
        self.rs_static.display_defaults_file = display
        #self.rs_static.display_defaults_file = self._processing_defaults['display']
        #self.rs_static.spectral_model = self._processing_defaults['spectral_model']
        self.rs_static.delta_t = delta_t
        #self.rs_static.corr_adjusts =  corr_adjusts
        self.rs_static.auto_loop = auto_loop
        self.rs_static.display_defaults = display_defaults
        self.rs_static.config = config
        self.rs_static.config
        self.rs_static.time_res = time_res
        self.rs_static.requested_times=requested_times
        self.rs_static.alt_res = alt_res
        self.rs_static.max_range_bin = max_range_bin
        self.rs_static.n_range_ave = n_range_ave
        self.rs_static.compute_stats\
                = self.rs_static.processing_defaults.enabled('compute_stats')
        self.rs_static.data_request = data_request
        self.rs_static.display_defaults_file

       
        # These parameters may be changed by "update_cal_and_process" as
        # it steps through the data.

        self.rs_init.interval_start_time = interval_start_time
        self.rs_init.interval_end_time = interval_end_time
        self.rs_init.rs_data_end_time = rs_data_end_time
        self.rs_init.current_file_end_time = current_file_end_time

        # data and calibration objects

        self.rs_init.rs_Cxx = rs_Cxx
        self.rs_init.rs_cal = rs_cal
        self.rs_init.rs_soundings = rs_soundings
        self.rs_init.sounding = sounding
        self.rs_init.rs_mem = rs_mem
        self.rs_init.last_sounding_time = last_sounding_time

        return 


# for testing

if __name__ == '__main__':
    try:
        print "r = Rti('gvhsrl','now', 0.5, 0.0, 18, 3.0,'flight_plots.json')"
    except RuntimeError, msg:
        print 'WARNING: ', msg
    except OSError, msg:
        print 'WARNING: ', msg

