import lg_base.graphics.graphics_toolkit as gt
import numpy as np

# Try to use the much faster nanmean from bottleneck, otherwise fall back
# to the scipy.stats version
try:
    from bottleneck import nanmean,nansum
except ImportError:
    print
    print "No bottleneck.nanmean available! Falling back to SLOW scipy.stats.nanmean"
    print
    from scipy.stats import nanmean,nansum
    
def show_radar(instrument,display_defaults,rs,usetimes,usealts,figs):
            #toplevel
    if usealts==None and hasattr(rs,'heights'):
                usealts=rs.heights
    if usetimes==None and hasattr(rs,'times'):
                usetimes=rs.times

    plot_alt_index = rs.plot_alt_index if hasattr(rs,'plot_alt_index') else 0
    plot_altitude = usealts[plot_alt_index]
    layer_indices = rs.layer_indices if hasattr(rs,'layer_indices') else np.arange(len(usealts))
        
    print
    print
    print 'plot_alt_index = ',plot_alt_index
    print 'plot altitude =  ',plot_altitude
    print
            
    print display_defaults.enabled('radar_backscatter_image')
    if display_defaults.enabled('radar_backscatter_image') and hasattr(rs,'Backscatter'):
                title_str = instrument+' backscatter cross section '
      
                print 'rendering ',title_str
                try:
                    gt.rti_fig('radar_backscatter_image'
                    , instrument
                    , rs.Backscatter
                    , usetimes
                    , usealts
                    , title_str
                    , '1/(m sr)'
                    , None                 
                    , display_defaults
                    ,figs)
                except AttributeError:
                    raise RuntimeError, "show_images - no data for radar backscatter plot"

    if display_defaults.enabled('radar_backscatter_profile') and hasattr(rs,'Backscatter'):
                temp = rs.Backscatter.copy()
               
                temp[temp <= 0.0] = np.NaN
                radar_bs_profile = nanmean(temp,0)

               
                try:
                   gt.plot_vs_altitude('radar_backscatter_profile'                     
                    ,instrument                   
                    ,usetimes                        
                    ,usealts                    
                    ,[radar_bs_profile]
                    ,'r'                   
                    ,[2]
                    ,[]                
                    ,None          
                    ,'Backscatter cross section'                  
                    ,'1/(m sr)'                         
                    ,'backscatter'            
                    ,False    #clear figure                
                    ,display_defaults                    
                    ,figs)
                
                except: #AttributeError:
                    print
                    print "show_images - error plotting radar_backscatter_profile"
                    print
            
    if display_defaults.enabled('radar_reflectivity_image') and hasattr(rs,'Reflectivity'):
                title_str = instrument+' reflectivity '
                print 'rendering ',title_str
                try:
                    gt.rti_fig('radar_reflectivity_image'
                    , instrument
                    , rs.Reflectivity
                    , usetimes
                    , usealts
                    , title_str
                    , 'dBz'
                    , None             
                    , display_defaults
                    ,figs)
                except AttributeError:
                    #raise
                    raise RuntimeError, "show_images - no data for radar reflectivity plot"

    if display_defaults.enabled('reflectivity_vs_time') and hasattr(rs,'Reflectivity'):
                title_str = ' reflectivity vs time  '+str(plot_altitude/1000.0) +' km' 
                print 'Rendering ',title_str
                try:
                    gt.plot_vs_time('reflectivity_vs_time'
                       ,instrument
                       ,usetimes
                       ,[rs.Reflectivity[:,plot_alt_index]]  
                       ,['r']
                       ,[2]
                       ,[]
                       ,'upper left'    
                       ,'Reflectivity'
                       ,'dBz'
                       ,title_str   
                       ,False #FIXME
                       ,display_defaults
                       ,figs)
                except AttributeError:
                    raise RuntimeError, "Radar display - no data for radar reflectivity vs time plot"
                
    if display_defaults.enabled('radar_velocity_image') and hasattr(rs,'MeanDopplerVelocity'):

                title_str = instrument+' Doppler velocity '
                print 'rendering ',title_str     
                try:
                    gt.rti_fig('radar_velocity_image'
                    , instrument
                    , rs.MeanDopplerVelocity
                    , usetimes
                    , usealts
                    , title_str
                    , 'm/s'
                    , None                             
                    , display_defaults
                    ,figs)
                except AttributeError:
                    raise RuntimeError, "show_images - no data for radar velocity plot"
   
    if display_defaults.enabled('radar_spectralwidth_image') and hasattr(rs,'SpectralWidth'):
                title_str = instrument+' spectral width '
                print 'rendering',title_str
                
                try:
                    gt.rti_fig('radar_spectralwidth_image'
                    , instrument
                    , rs.SpectralWidth
                    , usetimes
                    , usealts
                    , title_str
                    , 'm/s'
                    , None             
                    , display_defaults
                    ,figs)
                except AttributeError:
                    raise RuntimeError, "show_images - no data for radar spectral width plot"
                
    if display_defaults.enabled('platform_vertical_velocity') \
                 and hasattr(rs,'vertically_averaged_doppler'):
              
        title_str = 'plattform_vertical velocity'
        print 'rendering ',title_str
        try:
            gt.plot_vs_time('platform_vertical_velocity'
                  ,instrument
                  ,usetimes
                  ,[rs.vertically_averaged_doppler]  
                  ,['r']
                  ,[2]
                  ,[]
                  ,'upper left'    
                  ,'Platform vertical Velocity '
                  ,'m/s'
                  ,title_str   
                  ,False #FIXME
                  ,display_defaults
                  ,figs)
        except:
            print 'unable to plot platform vertical velocity'
   
    if display_defaults.enabled('radar_backscatter_vs_time') \
                 and hasattr(rs,'Backscatter'):          
        title_str = 'radar_backscatter vs time'
        print 'rendering ',title_str
        try:
            gt.plot_vs_time('radar_backscatter_vs_time'
                  ,instrument
                  ,usetimes
                  ,[rs.Backscatter[:,plot_alt_index]]  
                  ,['r']
                  ,[2]
                  ,[]
                  ,'upper left'    
                  ,'backscatter cross-section '
                  ,'1/(m sr)'
                  ,title_str   
                  ,False #FIXME
                  ,display_defaults
                  ,figs)   
	except AttributeError:
            raise RuntimeError, "plot_vs_time - no data for radar backscatter "

   
    if display_defaults.enabled('radar_doppler_velocity_vs_time') \
                 and hasattr(rs,'MeanDopplerVelocity'):          
        title_str = 'Doppler velocity vs time'
        print 'rendering ',title_str
        try:
            gt.plot_vs_time('radar_doppler_velocity_vs_time'
                  ,instrument
                  ,usetimes
                  ,[rs.MeanDopplerVelocity[:,plot_alt_index]]  
                  ,['r']
                  ,[2]
                  ,[]
                  ,'upper left'    
                  ,'Doppler velocity '
                  ,'m/s'
                  ,title_str   
                  ,False #FIXME
                  ,display_defaults
                  ,figs)   
	except AttributeError:
            raise RuntimeError, "plot_vs_time - no data for radar Doppler velocity "
        
            #BEGIN RENDERING PROCEDURES FOR PARTICLE
            #NOTES:
            # rs is the rs_particle structure. NOT the base 'rs' composite frame of all hsrl/radar/everything  Anything combining those fields should only be done in processing, not visualization
            # a separate artist for specific comparisons may be made in the future for this if its really desired... or rs may be pushed here anyway.

