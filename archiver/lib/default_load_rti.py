from __future__ import print_function

import logging
#logging.basicConfig(level=logging.DEBUG)
from datetime import datetime, timedelta
from maestro.rti_maestro import Rti
now = datetime.now()
endStr = now.strftime("%X")
delta = timedelta(minutes=5)
start = now - delta
startStr = start.strftime("%d-%b-%y %X")

print ('plotting from {0} to {1}'.format(startStr, endStr))

r = Rti(instruments='gvhsrl',start_time=startStr, plot_length = endStr, 
        min_alt = 1.5, max_alt = 15, display= 'flight_plots.json', mol_norm_alt = 2.1, 
        cmd_defaults='process_control.json')
