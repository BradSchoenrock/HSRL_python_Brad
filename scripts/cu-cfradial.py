#!/usr/bin/env python
#from hsrl.data_stream.rti import Rti

from maestro.rti_maestro import Rti

r = Rti('gvhsrl', '14-feb-12 17:10:00', '19:40:00', 0, 15, display='no_plots.json',t_res={'manual':15.0},z_res={'manual':30})
r.write_netcdf(netcdf_format='bm_hsrl_cfradial.cdl', tag='', output_dir='/tmp')

