import numpy as np
import lg_base.graphics.graphics_toolkit as gt
# Try to use the much faster nanmean from bottleneck, otherwise fall back
# to the scipy.stats version
try:
    from bottleneck import nanmean, nanmax
except ImportError:
    print
    print "No bottleneck.nanmean available! Falling back to SLOW scipy.stats.nanmean"
    print
    from scipy.stats import nanmean

#for plotting cooperate product hsrl and radars
def show_hsrl_radar(display_defaults,rs,parameters,usetimes,usealts,figs,rs_inv=None,rs_mean=None,rs_mwacr=None,rs_kazrge=None,**kwargs):
    pass