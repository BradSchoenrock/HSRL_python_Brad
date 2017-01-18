
import lg_base.graphics.graphics_toolkit as gt


def show_met(display_defaults,rs,usetimes,figs):



    if display_defaults.enabled('sfcmet_rain_rate'):
        if hasattr(rs,'rain_intensity_mean_wx1')\
              and hasattr(rs,'rain_intensity_mean_wx2'):
            print 'rendering rain rate from wx files'        
            gt.plot_vs_time('sfcmet_rain_rate'
                 ,'rain gauge'
                 ,rs.times
                 ,[rs.rain_intensity_mean_wx1,rs.rain_intensity_mean_wx2]
                 ,['g','b']
                 ,[2,2]
                 ,['wx1','wx2']
                 ,'upper left'    
                 ,'Rain intensity'
                 ,'mm/hr'
                 ,'Rain rate vs time'   
                 ,False #FIXME
                 ,display_defaults
                 ,figs)
        else:
            print 'no rain intensity plot----rs.rain_intensity_wx1 or wx2 not found'
    if display_defaults.enabled('sfcmet_rain_accumulated'):
           if hasattr(rs,'rain_accumulated_wx1')\
               and hasattr(rs,'rain_accumulated_wx2'):
               print 'rendering accumulated rain from wx files'  
               gt.plot_vs_time('sfcmet_rain_accumulated'
                    ,'rain gauge'
                    ,rs.times
                    ,[rs.rain_accumulated_wx1,rs.rain_accumulated_wx2]
                    ,['g','b']
                    ,[2,2]
                    ,['wx1','wx2']
                    ,'upper left'    
                    ,'Rain Accumulated'
                    ,'mm'
                    ,'Rain Accumulated'   
                    ,False #FIXME
                    ,display_defaults
                    ,figs)
           else:
               print 'no accumulated rain plot--rs.rain_accumulated_wx1 or wx2 missing'
    return
