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


#operates on count (and possibly inv) profiles
def process_hsrl_raman_profile(hsrl_profile,raman_profile,parameters,entire_frame=None):
    print
    print 'entering process_hsrl_raman_profile in cooperative-------------------'
    rs =hau.Time_Z_Group(like=hsrl_profile)
    rs.hsrl_profile =copy.deepcopy(hsrl_profile)
    rs.raman_profile = copy.deepcopy(raman_profile)
    print 'leaving process_hsrl_raman_profile'
    return rs


#operates on counts from merge files
def process_hsrl_ramanmerge(rs_mean,rs_merge,parameters,entire_frame=None):
    print 'entering process_hsrl_ramanmerge in cooperative-------------------'
    rs = copy.deepcopy(rs_mean)
    rs.raman_backscatt = rs_merge.beta_a_backscat
    return rs


#operates on inverted data
def process_hsrl_raman(rs_inv,rs_raman,parameters,entire_frame=None):
    print 'entering process_hsrl_raman in cooperative-------------------'
    return None
