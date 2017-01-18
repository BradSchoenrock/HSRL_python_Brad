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

def show_rain(instrument,display_defaults,rs,usetimes,usealts,figs):
                if hasattr(rs,'precip_rate'):
                    gt.plot_vs_time('rainguage_precip_rate'
                       ,instrument
                       ,usetimes
                       ,[rs.precip_rate]  
                       ,['b']
                       ,[2]
                       ,[]
                       ,''    
                       ,'raingauge precip rate (mm/hr)'
                       ,[]
                       ,'raingauge precip rate vs time'   
                       ,False #FIXME
                       ,display_defaults
                       ,figs)
                    print 'rendering rainguage precip rate'               
