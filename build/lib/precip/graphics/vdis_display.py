import lg_base.graphics.graphics_toolkit as gt
import numpy as np
import lg_base.core.array_utils as hau
try:
    from bottleneck import nanmean,nansum
except:
    print
    print 'No bottleneck.nanmean available! Falling back to SLOW scipy.stats.nanmean'
    print
    from scipy.stats import nanmean
    from numpy import nansum

def show_vdis(instrument,display_defaults,rs,usetimes,usealts,figs):
                if hasattr(rs,'rain_rate'):
                    gt.plot_vs_time('vdis_precip_rate'
                       ,instrument
                       ,usetimes
                       ,[rs.rain_rate]  
                       ,['b']
                       ,[2]
                       ,[]
                       ,''    
                       ,'vdis precip rate (mm/hr)'
                       ,[]
                       ,'vdis precip rate vs time'   
                       ,False #FIXME
                       ,display_defaults
                       ,figs)
                    print 'rendering vdis precip rate'               
