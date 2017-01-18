
import numpy as np
import lg_base.core.array_utils as hau
import scipy.special as ss
from datetime import datetime
import copy
try:
    from bottleneck import nanmean,nansum,anynan
except ImportError:
    print
    print 'No bottleneck.nanmean available! Falling back to SLOW scipy.stats.nanmean'
    print
    from scipy.stats import nanmean
    from numpy import nansum
    def anynan(x):
        return np.any(np.isnan(x))  

def process_hsrl_radar(process_parameters=None,rs_inv=None,rs_mean=None,rs_mwacr=None,rs_kazrge=None,**kwargs):
    #create timez group and add heights            
    rs_cooperative=hau.Time_Z_Group(rs_inv.times.copy(),timevarname='times',altname='heights')
    setattr(rs_cooperative,'heights',rs_inv.msl_altitudes.copy())



    print 'called process_hsrl_radar'
    print kwargs
    print rs_mwacr
    print rs_kazrge
    #raise RuntimeError('RUNNING HSRLRADAR')



    return rs_cooperative