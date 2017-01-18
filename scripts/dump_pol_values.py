#!/usr/bin/env python
from __future__ import print_function
import netCDF4
import sys
import numpy as np
DRIFT_THRESHOLD=5.0

def dump_pol_data(filename):
    np.set_printoptions(precision=3,suppress=True, threshold=7200)
    try:
	ds = netCDF4.Dataset(filename)
	year, month, day, hour, minute, second =  ds.variables['DATA_base_time'][0:6]
	month_str = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'][month-1]
    except KeyError, e:
            print ('could not get DATA_base_time from %s' % f)
	    return
    except RuntimeError, e:
            print ('could not open %s' % f)
	    return

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
	    cur_state = 'drifting'
	print("%s  %f.1 %d-%s-%d %02d:%02d:%02d '%s'" % (cur_state, median_angle_rate,
			day, month_str, year-2000, hour, minute, second, filename) )

	ds.close() # done with this dataset

if __name__ == '__main__':
    if len(sys.argv) != 2:
	print ('usage %s netcdf_file')
	sys.exit(1)

    dump_pol_data(sys.argv[1])

