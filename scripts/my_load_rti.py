import sys
import os,os.path
hsrl=os.path.join(os.path.expanduser('~'),'code','hsrl_python')
sys.path.append(hsrl)
from maestro.rti_maestro import Rti
#print """
#sample run:


r = Rti('gvhsrl', '29-jul-15 19:30', '20:00:00', 0, 15, display='bm_plots.json',t_res={'manual':15.0},z_res={'manual':30})

r.write_netcdf(netcdf_format='bm_hsrl_cfradial.cdl', tag='', output_dir='/tmp')

#"""

# default 04-Feb-12 14:00
# matt 13-dec-16 17:00
# bruce 29-jul-15 21:45
