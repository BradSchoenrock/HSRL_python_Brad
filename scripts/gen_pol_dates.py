#!/usr/bin/env python
from __future__ import print_function
import os
import time
import sys
import netCDF4
import numpy as np

DRIFT_THRESHOLD=5.0   # rotation rates smaller than this are considered fixed, rather than rotating

def list_files(path):
    directories=[]
    files = []

    data = os.path.abspath(os.path.expanduser(path))
    for root, dirs_o, files_o in os.walk(data):
        for name in dirs_o:
            directories.append(os.path.join(root, name))
        for name in files_o:
            file_path = os.path.join(root, name)
            if os.path.isfile(file_path) :
		froot, ext = os.path.splitext(name)
		# only process netcdf files with the correct prefix that are NOT in the badsum directory
		if ext == '.nc' and froot.find('gvhsrl',0) == 0 and root.find('badsum') == -1:
		    files.append(file_path)

    return directories, files

def scan_files(file_list):
    prev_state = 'unknown'
    output = []
    log = open('/tmp/gen_pol_dates.log','w+')
    drifting_log = open('/tmp/gen_pol_drifting.log','w+')

    for f in file_list:
    	print ('processing %s' % f, file = log)
	key = 'DATA_base_time'
	try:
	    ds = netCDF4.Dataset(f)
	    year, month, day, hour, minute, second =  ds.variables[key][0:6]
	    month_str = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'][month-1]
	except KeyError, e:
	    print ('could not get %s from %s' % (key, f), file = log)
	    continue
	except RuntimeError, e:
	    print ('could not open %s' % f, file = log)
	    continue

	cur_state = 'fixed'  # until proven otherwise
	if ds.variables.has_key('polarization'):
	  pol_var = ds.variables['polarization']
          if len(pol_var) > 1:
	    timeD = ds.variables['DATA_time']
	    timeDsec = timeD[:,3]*3600+timeD[:,4]*60+timeD[:,5]+timeD[:,6]*1e-3;
	    median_angle_rate = np.median(np.diff(pol_var)/np.diff(timeDsec)*180/np.pi)

	    # we have a polorization variable AND it changes value, then the polarizer is rotating
	    if median_angle_rate >DRIFT_THRESHOLD:
		cur_state = 'rotating'
	    elif median_angle_rate > 0.1:
		# don't flag measurement noise for the polarization angle as angle drift
		print("%s  %f.1 %d-%s-%d %02d:%02d:%02d '%s'" % ('drifting', median_angle_rate,
			    day, month_str, year-2000, hour, minute, second, f), file=drifting_log)
	
	# if our state changed, print out the time for the start of the file
	if prev_state != cur_state:
		output.append("%d-%s-%d %02d:%02d:%02d, '%s'" % (day,month_str,year-2000, hour, minute, second, cur_state),
		)
		print("%s %d-%s-%d %02d:%02d:%02d '%s'" % (cur_state, day, month_str, year-2000, hour, minute, second, f), file=log)
	    	prev_state = cur_state

	ds.close() # done with this dataset

    output.reverse()
    return output

 
if __name__ == '__main__':
    if len(sys.argv) != 2:
	print ("usage: %s directory" % (sys.argv[0]))
	sys.exit(1)

    dirs, files = list_files(sys.argv[1])
    if len(files) == 0:
	print("%s : no files found for %s" % (sys.argv[0], sys.argv[1]))
	sys.exit(2)
    files.sort()
    cal_vals = open('/tmp/calvals.txt', 'w+')
    output = scan_files(files)
    for l in output:
	print(l, file=cal_vals)
