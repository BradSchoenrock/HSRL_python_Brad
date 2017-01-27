#!/usr/bin/env python
#from hsrl.data_stream.rti import Rti

import sys

sys.path.insert(0, '.')

from maestro.rti_maestro import Rti

r = Rti('gvhsrl', '29-jul-15 19:30', '20:00:00', 0, 15, display='bm_plots.json',t_res={'manual':15.0},z_res={'manual':30})
r.write_netcdf(netcdf_format='bm_hsrl_cfradial.cdl', tag='', output_dir='/tmp')

