import sys
import os,os.path
# these 2 lines add the current directory to PYTHONPATH.  
# This should probably be done in the project attributes
#hsrl=os.path.join(os.path.expanduser('~'),'code','eol_hsrl_python')
#sys.path.append(hsrl)
from maestro.rti_maestro import Rti
#print """
#sample run:

r=Rti(instruments= 'gvhsrl',start_time = '17-Jul-12 2:00', 
plot_length=2.0, min_alt=1.5, max_alt=15, mol_norm_alt=2.1, 
display='mh_plots.json') #t_res={'n_dt':20}
#"""
