#!/opt/anaconda/envs/production/bin/python
from __future__ import print_function
from datetime import datetime, timedelta
from maestro.rti_maestro import Rti
from hsrl.data_stream.hsrl_read_utilities import get_path_to_data
from lg_base.core.locate_file import locate_file
import os,os.path
import subprocess

# set location of hsrl_python.json
os.environ['HSRL_CONFIG'] = '/usr/local/hsrl/config'
os.environ['GRIB_CACHE'] = '/data/GFS_cache'
# allow configuration files outside the hsrl_python tree
os.environ['OVERRIDE_SYSTEM_DEFAULTS']='1'
os.environ['LANG']='C'
now = datetime.now()
delta = timedelta(minutes=5)
start = now - delta

startStr = start.strftime("%d-%b-%y %X")
endStr = now.strftime("%X")
#data_root = get_path_to_data('gvhsrl', None)
data_root = '/data/hsrl/data.archiver'
#output_dir=os.path.join(data_root, start.strftime('%Y/%m/%d'), 'cfr_height')
date_dir = start.strftime('%Y%m%d')
output_dir=os.path.join(data_root, 'cfradial', 'uw', date_dir)
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
output = r.write_netcdf(netcdf_format=cdl, tag='', output_dir=output_dir)
print ('converting to CfRadial from  {0} to {1} '.format(startStr , endStr))
end = datetime.now()
duration = end - now
print ('conversion complete after', duration)

# create a latest data time string from a filename
# of the form gvhsrl_20150602T1606_20150602T1609.nc
bname = os.path.basename(output)
ldt = bname[bname.index('_')+1:]
ldt = ldt[:ldt.index('_')]  # yields 20150602T1606
ldt = ldt[0:8] + ldt[9:] + '00'
cmd = ('LdataWriter', '-dir', os.path.join(os.environ['DATA_DIR'], 'cfradial','uw'), '-dtype', 'nc', '-ext', 'nc', '-rpath', 
	os.path.join(date_dir, bname), '-ltime' ,ldt )
print ('invoking command = ', cmd)
subprocess.call(cmd)
print ('command completed')
