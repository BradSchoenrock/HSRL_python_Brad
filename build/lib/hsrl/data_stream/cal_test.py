#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os.path
from time import sleep
from subprocess import call
#from PyQt4.QtGui import QApplication

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from matplotlib.dates import date2num
from matplotlib.dates import num2date
import matplotlib.pyplot as plt
import hsrl.data_stream.hsrl_read_utilities as hru
import hsrl.calibration.calibration_utilities as cu
import lg_base.formats.calvals as cru
#import sounding_utilities as su
import processing_utilities as pu
import tz_utilities as tzu
import display_utilities as du
#import cal_file_generation as cfg
from get_scale_corrections import get_scale_corrections
#import rti_fig as rf
from main_routines import update_cal_and_process
import performance_model as pm
import json


class Rti:
  """ HSRL Plotting Routines

  Example:
  r = Rti('gvhsrl','18-Aug-11 00:00', 2, 1.5, 15, 2.1, 'alt_disp_def.json')

  r.next()

  """

  def __init__(
    self,
    instrument,
    start_time,
    plot_length,
    min_alt,
    max_alt,
    od_norm_alt,
    display_defaults_file = 'display_defaults.json'
    ):
    """ constructs a new Rti instance

        instrument          = hsrl id string (eg. 'ahsrl','gvhsrl','nshsrl','mf2hsrl').
        start_time          = time string, first data to plot (eg.'5-may-11 12:30').
        plot_length         = time period to plot in hours (e.g. 2) or
                              end time of first interval(e.g. '6-may-10 1:30').
        min_alt             = minimum altitude (msl) to display.
        max_alt             = maximum altitude (msl) to display.
        od_norm_alt         = normalization altitude for optical depth (m).
                              if not supplied od_norm_alt is taken from
                              calvals_xhsrl.txt
       display_defaults_file= name of alternate .json file containing display directives
                              this is optional, if not suppied values are taken
                              from 'display_defaults.json' 
    """
    plt.ion()
    (self.rs_static, self.rs_init) = rti_init(
      instrument,
      start_time,
      plot_length,
      min_alt,
      max_alt,
      od_norm_alt,
      display_defaults_file,
    )
    try:
      self.update()
    except RuntimeError, emsg:
      print 'Rti.__init__(): %s' % emsg
      return
    except OSError, emsg:
      print 'Rti.__init__(): %s' % emsg
      return
    self.display()

  def update(self, inCalibration=False):
    [self.rs, self.rs_init] = update_cal_and_process(self.rs_static,
                                                     self.rs_init, inCalibration)

  def new_calvals(self):
    """ read new calibration values """

    print 'Reading new calibration values'
    #read new calvals file
    last_i2_scale=self.rs_static.psel[4]
    #reread calvals file
    calvals=cru.cal_file_reader(self.rs_static.instrument)
    #select calvals for current time
    self.rs_init.rs_constants=\
      cru.select_sys_constants(calvals,self.rs_init.interval_start_time)
    #get half width,in bin numbers, of interval used to compute extinction cross section
    self.rs_init.rs_constants['extinction_bin_delta']=self.rs_static.display_defaults['extinction_bin_delta'][0]

    #update i2_scan adjustment from new calvals
    self.rs_static.psel[4]=self.rs_init.rs_constants['i2_scan_adjustment']
    #if new i2 correction was supplied, compute new calibration constants
    if self.rs_static.psel[4] != last_i2_scale:
      print 'i2 scaling changed from ',last_i2_scale,' to ',self.rs_static.psel[4]
      self.rs_init.rs_Cxx=cu.quick_cal(self.rs_init.rs_cal.i2scan.data\
                                       ,self.rs_init.rs_cal.i2scan.constant,self.rs_init.sounding\
                                       ,self.rs_init.rs_constants['wavelength']\
                                       ,self.rs_static.spectral_model,self.rs_static.psel)  
    try:
      self.update()
    except RuntimeError, emsg:
      print 'Rti.new_calvals(): %s', emsg
    except OSError, emsg:
      print 'Rti.new_calvals(): %s', emsg
    self.display()

  def print_fig(self):
    """print requested figure number"""
    fig_num=raw_input('figure #? ')
    plt.figure(fig_num)
    plt.savefig("myfilename.pdf", format="pdf")
    printer=self.rs_static.config.get_value('printer_id', 'name')
    printer='-d'+ printer[:]
    call(['/usr/bin/lp', printer,'myfilename.pdf'])

  def plot_model(self):
    """plot profile computed from model specs"""
    pm.performance_model(self.rs_init.rs_Cxx\
                         ,self.rs_init.rs_cal\
                         ,self.rs.profiles.dc_combined_hi_counts\
                         ,self.rs.rs_raw.transmitted_energy\
                         ,self.rs.rs_raw.seeded_shots\
                         ,self.rs.rs_raw.times\
                         ,self.rs_init.rs_constants)

  def cal_scale(self):
    """ rescale or turn off calibrations """

    last_i2_scale=self.rs_static.psel[4]
    self.rs_static.psel=get_scale_corrections(self.rs_static.psel,self.rs_static.config)
    #if new i2 correction was supplied, get new calibration constants
    if self.rs_static.psel[4] != last_i2_scale:
      self.rs_init.rs_Cxx=cu.quick_cal(self.rs_init.rs_cal.i2scan.data\
                                       ,self.rs_init.rs_cal.i2scan.constant,self.rs_init.sounding\
                                       ,self.rs_init.rs_constants['wavelength']\
                                       ,self.rs_static.spectral_model,self.rs_static.psel)  
    try:
      self.update()
    except RuntimeError, emsg:
      print 'Rti.cal_scale(): %s', emsg
    except OSError, emsg:
      print 'Rti.cal_scale(): %s', emsg
    self.display()


  def cal_gen(self):
    #generate new cal files: baseline, geo, and diff_geo files.
    print('\nProcess request at native altitude resolution with no geo corr for cal_gen')

    print ' '
    print 'you must make initial request for min_alt=0 and max_alt =35 km'
    print 'and also use "calibration.json" for display_defaults'
    print ' '

    #generate results at native altitude resolution without geo correction
    #only processses the last data chunk--does not work across cal or raob changes
    #calibration.json must be set to produce native altitude resolution


    self.rs_init.rs_constants['first_bin_to_process']=0
    self.rs_static.od_norm_alt=10000
    self.rs_static.min_alt=0
    self.rs_static.max_alt=35000
    self.rs_static.alt_res=7.5
    self.rs_static.psel[3]=0
    try:
      self.update()
    except RuntimeError, emsg:
      print 'Rti.cal_gen(): %s', emsg
    except OSError, emsg:
      print 'Rti.cal_gen(): %s', emsg
    cfg.generate_new_cal_file(self.rs_static.instrument,self.rs,self.rs_init.rs_cal\
                              ,self.rs_init.rs_constants)
    self.display()

  def replot(self):
    """ clear figures and replot last"""
    self.clear()
    self.display()

  def clear(self):
    """clear all figures"""
    # loop through all existing figures
    for x in plt._pylab_helpers.Gcf.get_all_fig_managers() :
      plt.figure(x.num)
      plt.clf()


  def display(self):
    du.show_images(
      self.rs_static.instrument,
      self.rs,
      self.rs_init.sounding,
      self.rs_init.rs_constants,
      self.rs_static.display_defaults,
      self.rs_init.last_sounding_time,
      self.rs_static.max_alt,
      self.rs_static.auto_loop,
    )

    # force a redraw
    for x in plt._pylab_helpers.Gcf.get_all_fig_managers(): 
#      print 'updating  %d' % x.num
      fig = plt.figure(x.num)
      #QApplication.processEvents()
      plt.draw()
    # JVA - this shouldn't be necessary, but without it, the last figure
    # stays behind the other figures
    plt.figure(plt._pylab_helpers.Gcf.get_all_fig_managers()[-1].num)
    #QApplication.processEvents()
    plt.draw()


    print '\n'
    print 'cmds-->r.next(), r.cnext(), r.loop(), r.new_calvals(), r.cal_scale()'
    print '       r.clear(), r.cal_gen(), r.print_fig(), r.replot(), r.plot_model()'
    print '       r.help()'

  def help(self):
    print '\n'
    print 'r.next()------- plot next time interval on top of last plot'
    print 'r.cnext()------ clear plot and plot next time interval'
    print 'r.clear()------ clear old plots'
    print 'r.replot()----- replot last time interval'
    print 'r.new_calvals-- reloads calvals_xxhsrl.txt before next set of plots'
    print 'r.cal_gen()---- generate baseline, geofile, or diff_geofile cal files'
    print 'r.cal_scale()-- scale or disable various data corrections'
    print 'r.loop()------- automatically advance through data'
    print 'r.print_fig()-- print figure on printer designated in "display_defaults.json"'
    print 'r.plot_model()- plot signals expected from model specifications' 
    print 'r.help()------ help'


  def loop(self):
    done = False
    while not done:
      fig_numbers = [x.num
                     for x in plt._pylab_helpers.Gcf.get_all_fig_managers()] 

      for i in range(len(fig_numbers)):
        plt.figure(fig_numbers[i])
        plt.clf()
      try:
        self.next()
      except RuntimeError as exc:
        print exc
        done = True
      except OSError as exc:
        print exc
        done = True

    # this seems to be required
    self.replot() 


  def next(self, repeat=0 ):
    """ plot next data """
    now = date2num(datetime.utcnow())
    print 'Current time is= ', num2date(now)

    if self.rs_init.interval_end_time > now:
      self.rs_init.interval_end_time = now
      self.rs_init.interval_start_time = now \
        - self.rs_static.delta_t
    elif repeat == 0:
      self.rs_init.interval_start_time = \
        self.rs_init.interval_start_time + self.rs_static.delta_t
      self.rs_init.interval_end_time = self.rs_init.interval_start_time \
        + self.rs_static.delta_t
      now = date2num(datetime.utcnow())
    if self.rs_init.interval_start_time > now:
      self.rs_init.interval_start_time = now - self.rs_static.delta_t
      self.rs_init.interval_end_time = now
    elif self.rs_init.interval_end_time > now:
      self.rs_init.interval_end_time = now
      self.rs_init.interval_start_time = self.rs_init.interval_end_time \
        - self.rs_static.delta_t
      print 'end time = now '

    # don't catch the RuntimeError exception, so Rti.loop() will see it
    self.update()
    self.display()

  def cnext(self, repeat=0):
    """ plot next data wihout overplotting graphs """
    self.clear()
    self.next(repeat)

# for testing 
if __name__ == '__main__':
  try:
    r = Rti('gvhsrl','now', 0.5, 0.0, 18, 3.0,'flight_plots.json')
  except RuntimeError, msg:
    print 'WARNING: ', msg
  except OSError, msg:
    print 'WARNING: ', msg
