import sys
import os,os.path
hsrl=os.path.join(os.path.expanduser('~'),'code','hsrl_python')
sys.path.append(hsrl)
from maestro.rti_maestro import Rti
#print """
#sample run:

r=Rti(instruments= 'gvhsrl',start_time = '29-jul-15 19:30', 
plot_length=0.5, min_alt=1.5, max_alt=15, mon_norm_alt=2.1, 
display='bm_plots.json')
#"""

# default 04-Feb-12 14:00
# matt 13-dec-16 17:00
# bruce 29-jul-15 21:45
