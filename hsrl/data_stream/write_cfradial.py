#!/usr/bin/python
# -*- coding: utf-8 -*-

from hsrl.dpl.dpl_hsrl import dpl_hsrl
from netCDF4 import Dataset
import datetime
import lg_dpl_toolbox.dpl.dpl_create_cfradial as cfr
import lg_dpl_toolbox.dpl.dpl_artists as artists
import lg_dpl_toolbox as lgtb
from lg_base.core.locate_file import locate_file
import traceback

def write_cfradial(output_filename, start_dt, minutes, timeres_s = 5, 
                   altres_m =60, maxtimeslice_td=datetime.timedelta(seconds=30*60), 
                   instrument='gvhsrl', min_alt_m=0, max_alt_m=5000,store_calibrations=False):
    """ writes HSRL data in CfRadial netcdf format

          output_filename        = where to write the data
          start_td             = = datetime.datetime object first time to retrieve.
          minutes                = how many minutes to process
          timeres_s              = time resolution in seconds (native would be 2.5)
          altres_m               = altitude resolution in meters
          maxtimeslice_timedelta = datetime.timedelta object for amount of data processed (safe is 1 or 2 hours)
          instrument             = hsrl id string (eg. 'ahsrl','gvhsrl','nshsrl','mf2hsrl').
          min_alt_m              = minimum altitude in meters to display
          max_alt_m              = maximum altitude in meters to display.
    """

    cdl = locate_file('hsrl_cfradial.cdl', forModule=lgtb)
    print 'CDL = ', cdl
    timeres_td = datetime.timedelta(seconds=timeres_s)

    netcdf = Dataset(output_filename, 'w', clobber=True)
    delta = datetime.timedelta(minutes=minutes)
    timeres_delta = datetime.timedelta(seconds=timeres_s)
    end_dt = start_dt + delta

    gen = dpl_hsrl(instrument)

    if store_calibrations: # to store calibrations, newer actors are needed, as well as the precall methods (FIXME better design)
        import maestro.netcdf_precall as npc
        args=[]
        kwargs=dict(output=netcdf,template=cdl,usecfradial=True,basetime=start_dt)
        x=npc.addConstantsToParms(npc.addCalibrationsToNetCDF())
        hsrlnar=gen(start_dt, end_dt, timeres_timedelta=timeres_delta, min_alt_m=min_alt_m, max_alt_m=max_alt_m,  altres_m=altres_m)
        x(hsrlnar,args,kwargs)
        nar=artists.dpl_netcdf_artist(hsrlnar,*args,**kwargs)
        #framestream,template,outputfilename=None,format=None,usecfradial=None,selected_bindings=None,output=None,forModule=None,withUnlimited=None,basetime=None,addAttributes={}):
        for x in nar:
          pass
    else:
        v = None
        try:
            # store each lidar record
            for tzg in gen(start_dt, end_dt, timeres_timedelta=timeres_delta, min_alt_m=min_alt_m, max_alt_m=max_alt_m,  altres_m=altres_m):
                if v == None:
                    v = cfr.DplCreateCfradial(cdl, netcdf, tzg)
                v.append_data(tzg)

            v.close()

        except RuntimeError, msg:
            print msg
            traceback.print_exc()
            print 'write_cfradial: could not process data for %s starting at %s' % \
                  (instrument, start_dt.strftime('%Y-%m-%d %H:%M:%S'))

if __name__ == '__main__':
#  start_dt =datetime.datetime(2012, 6, 20, 1, 0, 0)
    # 2012/2/22 is a research flight
    start_dt = datetime.datetime(2012, 2, 22, 16, 0, 0)
    write_cfradial('cf_radialtest.nc', start_dt, 10)  
    write_cfradial('cf_radialtest_120_10.nc', start_dt, 10, altres_m=120, timeres_s=5)