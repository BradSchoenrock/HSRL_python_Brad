#!/opt/anaconda/envs/production/bin/python
from __future__ import print_function
from datetime import datetime, timedelta
from maestro.rti_maestro import Rti
from hsrl.data_stream.hsrl_read_utilities import get_path_to_data
from lg_base.core.locate_file import locate_file
import os,os.path
import sys

def create_cfradial(startStr, endStr, output_dir='/tmp'):
    # set location of hsrl_python.json
    os.environ['HSRL_CONFIG'] = '/usr/local/hsrl/config'
    os.environ['GRIB_CACHE'] = '/data/GFS_cache'
    # allow configuration files outside the hsrl_python tree
    os.environ['OVERRIDE_SYSTEM_DEFAULTS']='1'
    os.environ['LANG']='C'

    try:
	os.mkdir(output_dir)
    except OSError:
	pass

    cdl=os.path.join(os.environ['HSRL_CONFIG'],'5min_hsrl_cfradial.cdl')

    r = Rti(instruments='gvhsrl', start_time=startStr, plot_length = endStr, 
	    min_alt =0, max_alt=15, display='no_plots.json', 
	    cmd_defaults='process_control.json',
	    t_res={'manual':0.5},
	    z_res={'manual':7.5})
    print ('converting to CfRadial from  {0} to {1} '.format(startStr , endStr))
    r.write_netcdf(netcdf_format=cdl, tag='', output_dir=output_dir)

if __name__ == '__main__':
    usage = '''
    {0} start_time end_time [output_dir] 
e.g.:
"02-Jun-2015 16:45:00" "02-Jun-2015 16:50:00" /tmp
(Note quotes around each date/time)
'''
    startStr, endStr = '02-Jun-2015 16:45:00', '02-Jun-2015 16:50:00'
    outputDir = '/tmp'
    if len(sys.argv) == 2 or len(sys.argv) == 1:
        print(usage.format(sys.argv[0]))
	sys.exit(1)
    if len(sys.argv) == 4:
	outputDir=sys.argv[3]
    if len(sys.argv) == 3 or len(sys.argv)==4 :
	startStr = sys.argv[1]
	endStr = sys.argv[2]
	
    create_cfradial(startStr, endStr, outputDir)
