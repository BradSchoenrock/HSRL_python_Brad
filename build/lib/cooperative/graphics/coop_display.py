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

def define_multiple_alt_plots(var,altitudes,indices,lines
                        ,widths,colors,legends,this_width=None
                        ,legend_addon=None,alt_color_list=None):
    """
       [lines,colors,widths,legends]=define_multiple_alt_plots(var,altitudes \
                        ,indices,lines,widths,,colors,legends,this_width=None)
       provides list of variables at altitudes[indices] for plot_vs_time
       var = variable to plot at altitudes[indices[:]]
       altitudes      = altitude vector corresponding to variable altitudes scale
       indices        = list of altitude indices to plot
       lines          = list of variables at requested levels
       widths         = list of line widths, set = 2 unless this_width is specified
       colors         = list of line colors
       legends        = list of legends providing the altitude for each plot
       this_width     = width of plot lines for this call if included
       alt_color_list = if not == None, uses different list of colors
       make initial call with lines, widths, colors, and legends as empty lists
       subsequent calls can add additional varaibles to the plot. New varia

    """
     
    width = 2
    #override default value of 2 if this width is provided
    if not this_width == None:
        width = this_width
    if alt_color_list == None:
        color_list = ['r','b','g','k','c','m','r','b','g','k','c','m']
    else:
        color_list =  ['k','c','m','r','b','g','k','c','m','r','b','g']
        
    if not isinstance(indices,(list,tuple,np.ndarray)):
      indices=[indices]
    list_size = len(lines)
    for i in range(len(indices)):
        lines.append(var[:,indices[i]])
        widths.append(width)
        #i_color = np.mod(list_size+i,len(color_list))
        i_color = i
        colors.append(color_list[i_color])
        if legend_addon == None:
            legends.append('%4.2f km'%(altitudes[indices[i]]/1000.0))
        else:
            legends.append(('%4.2f km'%(altitudes[indices[i]]/1000.0))
                     +legend_addon)           
    return lines,colors,widths,legends





def show_spheroid_particle(instrument,display_defaults,rs,particle_parameters,usetimes,usealts,radarLambda,figs,rs_radar=None,rs_inv=None,rs_pars2S2=None,rs_pars2S1=None,rs_marinemet=None,size_dist=None):

    #fix me--single bit radar mask because intepolation not aware of summode
    mask=None
    if mask is None and hasattr(rs,'qc_mask'):
        mask=np.ones(rs.qc_mask.shape,dtype='uint16')
    if mask is None and hasattr(rs,'effective_diameter_prime'):
        mask = np.ones(rs.effective_diameter_prime.shape,dtype='uint16')
    if mask is None and hasattr(rs_radar,'qc_radar_mask'):
        mask = np.ones(rs_radar.qc_radar_mask.shape,dtype='uint16')

    if rs_radar!=None and display_defaults.get_value('mask_image','radar_sn_mask'):
        mask[rs_radar.qc_radar_mask >= 0.5] = 1
        mask[rs_radar.qc_radar_mask < 0.5] = 0

    plot_alt_index = rs.plot_alt_index if hasattr(rs,'plot_alt_index') else 0
    layer_indices = rs.layer_indices if hasattr(rs,'layer_indices') else np.arange(len(usealts))
        
    if display_defaults.get_value('mask_image','mol_sn_ratio') and hasattr(rs,'qc_mask'):    
        cmask=np.zeros_like(rs.qc_mask,dtype='uint32')
        #2**4 is the mol_sn_ratio bit in qc_mask
        cmask = rs.qc_mask & 16 > 0
        mask &= cmask
    
    if usealts is not None:
           
          print 'plots depicting layer properties include '+str(len(layer_indices)) \
                 +' bins '+ str(usealts[layer_indices[0]])+'-->'+str(usealts[layer_indices[-1]])+ ' km'
          print
          print

    #define particle parameters string for use in title of particle figures
    if particle_parameters['type'] == 'rain':
        p_params_str = ' %s , $\\alpha$=%3.1f,$\\gamma$=%3.1f' \
            %(particle_parameters['size_distribution']['form']
            ,particle_parameters['size_distribution']['alpha_water']\
            ,particle_parameters['size_distribution']['g_water'])
    else:
        p_params_str = ' %s, $\\alpha_w$=%3.1f,$\\gamma_w$=%3.1f,\\alpha_i$=%3.1f,$\\gamma_i$=%3.1f' \
            %(particle_parameters['size_distribution']['form']
            ,particle_parameters['size_distribution']['alpha_water']\
            ,particle_parameters['size_distribution']['g_water']\
            ,particle_parameters['size_distribtution']['alpha_ice']\
            ,particle_parameters['size_distribution']['g_ice'])

        
    if display_defaults.enabled('rain_extinction_image') and hasattr(rs_inv,'extinction_aerosol'):
        title_str = ' rain extinction '
        background_extinction = particle_parameters['background_aerosol_bs']
        if particle_parameters['ext_source'] == 'bs/p180':
            rain_extinction = rs_inv.beta_a_backscat/particle_parameters['p180_water'] \
                       - particle_parameters['background_aerosol_bs']
        else :
            rain_extinction = rs_inv.extinction_aerosol - particle_parameters['background_aerosol_bs']
            
        print 'rendering ',title_str
        units = '1/m'
        try:
	    gt.rti_fig('rain_extinction_image'
                    , instrument
                    , rain_extinction
                    , usetimes
                    , usealts
                    , title_str
                    , units
                    , mask             
                    , display_defaults
                    ,figs)
	except AttributeError:
	        raise RuntimeError, "show_images - rain extinction"



            
    if display_defaults.enabled('effective_diameter_prime_image') and hasattr(rs,'effective_diameter_prime'):
        title_str = ' effective diameter prime '
        if display_defaults.get_value('effective_diameter_prime_image','units') == 'mm':
            units = 'mm'
            scale = 1e3
            print
            print '****************effective diameter image units in mm'
            print
        else:
            units = 'microns'
            scale = 1e6
            print
            print '*****************effective diameter image units in microns'
            print
        print 'rendering ',title_str
        try:
	    gt.rti_fig('effective_diameter_prime_image'
                    , instrument
                    , rs.effective_diameter_prime * scale
                    , usetimes
                    , usealts
                    , title_str
                    , units
                    , mask             
                    , display_defaults
                    ,figs)
	except AttributeError:
	        raise RuntimeError, "show_images - no data for effective diameter prime plot"

    if display_defaults.enabled('effective_diameter_image') and hasattr(rs,'effective_diameter'):
       
        title_str = 'Eff diameter, '+ p_params_str
	print 'rendering ',title_str
        
        if display_defaults.get_value('effective_diameter_image','units') == 'mm':
            units = 'mm'
            scale = 1e3
        else:
            units = 'microns'
            scale = 1e6

	try:
	    gt.rti_fig('effective_diameter_image'
                    , instrument
                    , rs.effective_diameter * scale
                    , usetimes
                    , usealts
                    , title_str
                    , units
                    , mask            
                    , display_defaults
                    ,figs)
	except AttributeError:
                    raise RuntimeError, "show_images - no data for effective diameter plot"
                
    if display_defaults.enabled('mode_diameter_image') and hasattr(rs,'mode_diameter'):             
        
         title_str = 'Mode dia, '+ p_params_str       
         if display_defaults.get_value('mode_diameter_image','units') == 'mm':
             units = 'mm'
             scale = 1e3
         else:
             units = 'microns'
             scale = 1e6

         print 'rendering ',title_str
         try:
             gt.rti_fig('mode_diameter_image'
             , instrument
             , rs.mode_diameter * scale 
             , usetimes
             , usealts
             , title_str
             , units
             , mask            
             , display_defaults
             ,figs)
         except AttributeError:
             raise RuntimeError, "show_images - no data for mode diameter plot"
         
    if display_defaults.enabled('mean_diameter_image') and hasattr(rs,'mean_diameter'):             
        
         title_str = 'Mean dia, '+ p_params_str        
         if display_defaults.get_value('mean_diameter_image','units') == 'mm':
             units = 'mm'
             scale = 1e3
         else:
             units = 'microns'
             scale = 1e6

         print 'rendering ',title_str
         try:
             gt.rti_fig('mean_diameter_image'
             , instrument
             , rs.mode_diameter * scale
             , usetimes
             , usealts
             , title_str
             , units
             , mask            
             , display_defaults
             ,figs)
         except AttributeError:
             raise RuntimeError, "show_images - no data for mean diameter plot"

    if display_defaults.enabled('mean_mass_diameter_image') and hasattr(rs,'mean_mass_diameter'):             
         
         title_str = 'Mean mass dia, '+ p_params_str       
         if display_defaults.get_value('mean_mass_diameter_image','units') == 'mm':
             units = 'mm'
             scale = 1e3
         else:
             units = 'microns'
             scale = 1e6
         
         print 'rendering ',title_str
         try:
             gt.rti_fig('mean_mass_diameter_image'
             , instrument
             , rs.mode_diameter * scale
             , usetimes
             , usealts
             , title_str
             , units
             , mask            
             , display_defaults
             ,figs)
         except AttributeError:
             raise RuntimeError, "show_images - no data for mean mass diameter plot"
    
    if display_defaults.enabled('non_rayleigh_adjustment_image') and hasattr(rs,'non_rayleigh_adjustment'):
        
        title_str = \
	       'Non-Rayleigh radar adjustment, dist($\\alpha_w$=%3.1f, $\\gamma_w$=%3.1f)' \
               %(particle_parameters['size_distribution']['alpha_water']\
                 ,particle_parameters['size_distribution']['g_water'])            

        print 'rendering ',title_str
        if 1: #try:
	    gt.rti_fig('non_rayleigh_adjustment_image'
                    , instrument
                    , rs.non_rayleigh_adjustment
                    , usetimes
                    , usealts
                    , title_str
                    , ''
                    , None           
                    , display_defaults
                    ,figs)
	if 0: #except AttributeError:
	    raise RuntimeError, "show_images - no radar_p180_water  image"
    
    if display_defaults.enabled('particle_phase_image') and hasattr(rs,'phase'):
        title_str = 'particle phase'
        temp = rs.phase.copy()
        temp[np.isnan(temp)] = 3
        cb_labels=[gt.ClassColormapEntry(0,'water',(0,0,0.9)),   #blue
                      gt.ClassColormapEntry(1,'ice',(.9,.9,.9)),  #white
                      gt.ClassColormapEntry(3,'',(0,0,0))]         #black for undefined
	print 'rendering ',title_str
	try:
	    gt.rti_fig('particle_phase_image'
                    , instrument
                    , temp
                    , usetimes
                    , usealts
                    , title_str
                    , cb_labels
                    , None            
                    , display_defaults
                    ,figs)
	except AttributeError:
                    raise RuntimeError, "show_images - no data for effective diameter plot"           
	    
    if display_defaults.enabled('liquid_water_content_image') and hasattr(rs,'LWC'):
      
        title_str = 'Water content, '+ p_params_str        
        try:
	    gt.rti_fig('liquid_water_content_image'
                    , instrument
                    , rs.LWC * 1e3  #convert from kg/m^3 to gr/m^3
                    , usetimes
                    , usealts
                    , title_str
                    , '$g/m^3$'
                    , mask           
                    , display_defaults
                    ,figs)
	except AttributeError:
	    raise RuntimeError, "show_images - no data for liquid water content plot"

    
    if display_defaults.enabled('adjusted_p180_water_image') and hasattr(rs,'lidar_mie_p180_water'):
       
        title_str = \
	       'adjusted p180/4pi water,      dist($\\alpha_w$=%3.1f, $\\gamma_w$=%3.1f)' \
                %(particle_parameters['size_distribution']['alpha_water']
                  ,particle_parameters['size_distribution']['g_water']) 
        if 1: #try:
	    gt.rti_fig('adjusted_p180_water_image'
                    , instrument
                    , rs.lidar_mie_p180_water
                    , usetimes
                    , usealts
                    , title_str
                    , '1/sr'
                    , mask            
                    , display_defaults
                    ,figs)
	if 0: #except AttributeError:
	    raise RuntimeError, "show_images - no adjusted p180 water image"
        
    if display_defaults.enabled('precip_rate_image') and hasattr(rs,'hsrl_radar_precip_rate'):    
      
        title_str = \
	       'precip rate,      dist($\\alpha_w$=%3.1f, $\\gamma_w$=%3.1f)' \
                %(particle_parameters['size_distribution']['alpha_water']
                  ,particle_parameters['size_distribution']['g_water'])        
        title_str = 'precip rate, , '+ p_params_str
        print 'rendering ',title_str
        try:
	    gt.rti_fig('precip_rate_image'
                    , instrument
                    , rs.hsrl_radar_precip_rate * 3600 * 1e3  #convert from m/s to mm/hr
                    , usetimes
                    , usealts
                    , title_str
                    , '$mm/hr$'
                    , mask            
                    , display_defaults
                    ,figs)
	except AttributeError:
	    raise RuntimeError, "show_images - no data for precip rate image"

    if display_defaults.enabled('radar_weighted_fall_velocity_image') and hasattr(rs,'rw_fall_velocity'):           
        title_str = 'rw fall velocity, '+ p_params_str          
        print 'rendering ',title_str
        try:
	    gt.rti_fig('radar_weighted_fall_velocity_image'
                    , instrument
                    , rs.rw_fall_velocity 
                    , usetimes
                    , usealts
                    , title_str
                    , '$m/s$'
                    , mask            
                    , display_defaults
                    ,figs)
	except AttributeError:
	    raise RuntimeError, "show_images - no data for computed radar-weighted fall velocity image"
    
    if display_defaults.enabled('model_spectral_width_image') \
               and hasattr(rs,'model_spectral_width'):
                  
            title_str = 'Model spectral width, '+ p_params_str        
            print 'rendering ',title_str
            try:
                gt.rti_fig('model_spectral_width_image'
                    , instrument
                    , rs.model_spectral_width
                    , usetimes
                    , usealts
                    , title_str
                    , '$m/s$'
                    , mask
                    , display_defaults
                    ,figs)
            except AttributeError:
                raise RuntimeError, "show_images - no data for model spectral width image"
  
    if display_defaults.enabled('particle_diameter_profiles') and hasattr(rs,'mode_diameter'):
        gt.plot_vs_altitude('particle_diameter_profiles'                     
                 ,'lidar_radar'                   
                 ,usetimes                
                 ,usealts              
                 ,[nanmean(rs.effective_diameter,0)*1e6,nanmean(rs.mean_diameter,0)*1e6,nanmean(rs.mode_diameter,0)*1e6]
                 ,['b','g','k']               
                 ,[]
                 ,['eff','mean','mode']                 
                 ,'upper right'          
                 ,'particle diameter'                  
                 ,'microns'                         
                 ,'particle diameters'            
                 ,None              
                 ,display_defaults                    
                 ,figs)    
   





    if display_defaults.enabled('fall_velocity_vs_time') \
       and hasattr(rs,'rw_fall_velocity')\
       and rs_radar!=None and hasattr(rs_radar,'MeanDopplerVelocity'):

        Doppler_good = rs_radar.MeanDopplerVelocity.copy()
        Doppler_good[mask==0]=np.NaN 
        fall_velocity_good = rs.rw_fall_velocity.copy()
        fall_velocity_good[mask==0] = np.NaN

        lines = []
        widths = []
        colors = []
        legends = []
        # line width to 1 for this plot--plot all doppler velocities
        [lines,colors,widths,legends]=define_multiple_alt_plots(
                rs_radar.MeanDopplerVelocity
                ,usealts
                ,plot_alt_index
                ,lines
                ,widths
                ,colors
                ,legends
                ,1
                ,' all Rdar'
                ,alt_color_list=True)
        #append doppler_good to lines
        #default linewidth of 2 for good doppler velocities
        [lines,colors,widths,legends]=define_multiple_alt_plots(
                Doppler_good
                ,usealts
                ,plot_alt_index
                ,lines
                ,widths
                ,colors
                ,legends
                ,2
                ,' good-radar'
                ,alt_color_list=True)
        #append all modelded fall velocities with linewidth = 1
        [lines,colors,widths,legends]=define_multiple_alt_plots(
                rs.rw_fall_velocity
                ,usealts
                ,plot_alt_index
                ,lines
                ,widths
                ,colors
                ,legends
                ,1,' model')
        #append good fall velocities with default linewidth = 2
        [lines,colors,widths,legends]=define_multiple_alt_plots(
               fall_velocity_good
                ,usealts
                ,plot_alt_index
                ,lines
                ,widths
                ,colors
                ,legends)

      
        title_str = 'radar_weighted velocity'
        print 'rendering ',title_str
        
        try:
            gt.plot_vs_time('fall_velocity_vs_time'
                  ,''
                  ,usetimes
                  ,lines  
                  ,colors
                  ,widths
                  ,legends
                  ,'upper left'    
                  ,'Radar weighted vertical Velocity '
                  ,'m/s'
                  ,title_str   
                  ,False #FIXME
                  ,display_defaults
                  ,figs)
            
	except AttributeError:
            raise RuntimeError, "plot_vs_time - no data for doppler , fall velocity "
  
    if display_defaults.enabled('doppler_and_model_velocities_vs_time') \
       and hasattr(rs,'rw_fall_velocity')\
       and rs_radar!=None and hasattr(rs_radar,'MeanDopplerVelocity'):

        Doppler_good = rs_radar.MeanDopplerVelocity.copy()
        Doppler_good[mask==0]=np.NaN 
        fall_velocity_good = rs.rw_fall_velocity.copy()
        fall_velocity_good[mask==0] = np.NaN

        lines = []
        widths = []
        colors = []
        legends = []
        # line width to 1 for this plot--plot all doppler velocities
        [lines,colors,widths,legends]=define_multiple_alt_plots(
                rs.mw_fall_velocity
                ,usealts
                ,plot_alt_index
                ,lines
                ,widths
                ,colors
                ,legends
                ,1
                ,' mw wieghted')
       
        #append all modelded fall velocities with linewidth = 1
        [lines,colors,widths,legends]=define_multiple_alt_plots(
                rs.rw_fall_velocity
                ,usealts
                ,plot_alt_index
                ,lines
                ,widths
                ,colors
                ,legends
                ,3
                 ,' rw model')
        #append good fall velocities with default linewidth = 2
        [lines,colors,widths,legends]=define_multiple_alt_plots(
               rs_radar.MeanDopplerVelocity
                ,usealts
                ,plot_alt_index
                ,lines
                ,widths
                ,colors
                ,legends
                ,alt_color_list =True
                ,legend_addon =  ' Doppler')

      
        title_str = 'fall velocity comparison'
        print 'rendering ',title_str
        
        try:
            gt.plot_vs_time('doppler_and_model_velocity_vs_time'
                  ,''
                  ,usetimes
                  ,lines  
                  ,colors
                  ,widths
                  ,legends
                  ,'upper left'    
                  ,'vertical Velocity '
                  ,'m/s'
                  ,title_str   
                  ,False #FIXME
                  ,display_defaults
                  ,figs)
            
	except AttributeError:
            raise RuntimeError, "plot_vs_time - no data for doppler and model velocities "

        
    if display_defaults.enabled('depol_and_doppler_velocity_vs_time') \
       and rs_radar!=None and hasattr(rs_radar,'MeanDopplerVelocity'):

        
        lines = []
        widths = []
        colors = []
        legends = []
        # line width to 2 for this plot--plot doppler velocities
        [lines,colors,widths,legends]=define_multiple_alt_plots(
                rs_radar.MeanDopplerVelocity
                ,usealts
                ,plot_alt_index
                ,lines
                ,widths
                ,colors
                ,legends
                ,2
                ,' DopV'
                ,alt_color_list=True)
       
        # line width to 2 for this plot--plot depol *10

        [lines,colors,widths,legends]=define_multiple_alt_plots(
                rs_inv.linear_depol*10
                ,usealts
                ,plot_alt_index
                ,lines
                ,widths
                ,colors
                ,legends
                ,2
                ,' cpol*10')
        
        title_str = \
	    'depol and Doppler velocity'
        print 'rendering ',title_str
        try:
            gt.plot_vs_time('depol_and_doppler_velocity_vs_time'
                  ,''
                  ,usetimes
                  ,lines  
                  ,colors
                  ,widths
                  ,legends
                  ,'upper left'    
                  ,'vel, depolarization * 10'
                  ,'m/s'
                  ,title_str   
                  ,False #FIXME
                  ,display_defaults
                  ,figs)
            
	except AttributeError:
            raise RuntimeError, "plot_vs_time - no data for depol and doppler  "
    
    
    if display_defaults.enabled('deff_prime_and_approximation_vs_time'):
        units = display_defaults.get_value('deff_prime_and_approximation','units')
        """plot effective diameter prime along with approximate value using constant p180 
           and Rayleigh approximation for radar.
           """
        
        if units == 'microns':
            scale = 1e6
        else: #mm
            scale =1e3
            
        if radarLambda < 5.0e-3: #3 mm radar
             k_sq_water = 0.7
             k_sq_ice = 0.186
        elif radarLambda < 1e-2: # 8.6 mm radar
             k_sq_water = 0.88
             k_sq_ice = 0.176
        else:    #x-band
             k_sq_water = 0.93
             k_sq_ice = 0.176

        p180_ice = particle_parameters['p180_ice']
        p180_water = particle_parameters['p180_water']
        radar_const_water = radarLambda * (2.0 * p180_water /(np.pi**3 * k_sq_water))**0.25
        radar_const_ice = radarLambda *(2.0 * p180_ice / (np.pi**3 * k_sq_ice))**0.25

        aerosol_background = particle_parameters['background_aerosol_bs']
        defp_water = radar_const_water *(rs_radar.Backscatter\
                    /(rs_inv.beta_a_backscat - aerosol_background))**0.25
        defp_ice = radar_const_ice *(rs_radar.Backscatter\
                      /(rs_inv.beta_a_backscat - aerosol_background))**0.25

        defp_ice[rs.phase == 0] = np.NaN
        defp_water[rs.phase == 1] = np.NaN
        lines = []
        widths = []
        colors = []
        legends = []
        # line width to 3 for effective_diameter_prime plot
        [lines,colors,widths,legends]=define_multiple_alt_plots(
                rs.effective_diameter_prime * scale
                ,usealts
                ,plot_alt_index
                ,lines
                ,widths
                ,colors
                ,legends
                ,3
                ," d_eff'")
            
            
        # line width to 1 for constant p180 Rayleigh water approx plot
        [lines,colors,widths,legends]=define_multiple_alt_plots(
                defp_water*scale
                ,usealts
                ,plot_alt_index
                ,lines
                ,widths
                ,colors
                ,legends
                ,1
                ," approx-water")
            
        # line width to 1 for constant p180 Rayleigh ice approx plot
            
        [lines,colors,widths,legends]=define_multiple_alt_plots(
                defp_ice*scale
                ,usealts
                ,plot_alt_index
                ,lines
                ,widths
                ,colors
                ,legends
                ,1
                ," approx-ice"
                ,alt_color_list=True)    
            
        title_str = "deff',constant p180-Ray approx"                 
        print 'rendering ',title_str
        try:
            gt.plot_vs_time('deff_pime_and_approximation_vs_time'
                  ,''
                  ,usetimes
                  ,lines  
                  ,colors
                  ,widths
                  ,legends
                  ,'upper left'    
                  ,'deff prime'
                  ,units
                  ,title_str   
                  ,False #FIXME
                  ,display_defaults
                  ,figs)
            
        except AttributeError:
            raise RuntimeError, "plot_vs_time - no data for deff prime and approximation "

   
    if display_defaults.enabled('model_fall_vel_vs_mode_dia') \
            and hasattr(rs,'mode_diameter'):       
        title = 'model fall vel vs dia'
        
        #sphere fall speeds (m/s) vs d(mm) from literature 
        d=np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8])
        v=np.array([0.27,0.72,1.17,1.62,2.06,2.47,2.87,3.27,4.03,3.67,4.03,5.17,5.65,6.09])
        
        gt.plot_xy('model_fall_vel_vs_mode_dia'
                ,''
                ,usetimes
                ,[d[:],rs.mode_diameter*1000.0, rs.mode_diameter*1000.0]
                ,[v[:],rs.mw_fall_velocity, rs.rw_fall_velocity]
                ,['c','k','r']
                ,['+','*','.']
                ,[2,5,5]
                ,['None']
                ,[0,0,0]
                ,['sphere','mw_vel','rw_vel']
                ,'lower right'
                ,'Mode diameter'
                ,'mm'
                ,'model fall velocity'
                ,'m/s'
                ,title
                ,[]
                ,[]
                ,[]
                ,[]
                ,display_defaults
                ,figs)


        
    if display_defaults.enabled('depol_vs_doppler_velocity') \
            and hasattr(rs_radar,'MeanDopplerVelocity'):
        title = 'Depol vs doppler'
        temp_vel = rs_radar.MeanDopplerVelocity[:,layer_indices].copy()
    
        temp_depol = rs_inv.linear_depol[:,layer_indices].copy()
        
        gt.plot_xy('depol_vs_doppler_velocity'
                ,''
                ,usetimes
                ,[temp_vel,temp_vel[rs_radar.Backscatter[:,plot_alt_index] > 1e-7]]
                ,[temp_depol*100,temp_depol[rs_radar.Backscatter[:,plot_alt_index] >1e-7]*100]
                ,['k','r']
                ,['*','.']
                ,[5,5]
                ,['None']
                ,[0,0]
                ,[]
                ,'lower right'
                ,'Doppler velocity'
                ,'m/s'
                ,'Depolarization'
                ,'%  '
                ,title
                ,[]
                ,[]
                ,[]
                ,[]
                ,display_defaults
                ,figs)
        
    if display_defaults.enabled('depol_vs_deff_prime'):
        title = 'Depol vs deff_prime'
        temp_deff_prime = rs.effective_diameter_prime[:,layer_indices].copy()
    
        temp_depol = rs_inv.linear_depol[:,layer_indices].copy()
        
        gt.plot_xy('depol_vs_doppler_velocity'
                ,''
                ,usetimes
                ,[temp_deff_prime[rs_radar.Backscatter[:,plot_alt_index] < 1e7]*1e6\
                  ,temp_deff_prime[rs_radar.Backscatter[:,plot_alt_index] >= 1e-7]*1e6\
                  ,temp_deff_prime[rs_radar.Backscatter[:,plot_alt_index] >1e-6]*1e6]\
                ,[temp_depol[rs_radar.Backscatter[:,plot_alt_index]<1e7]*100\
                  ,temp_depol[rs_radar.Backscatter[:,plot_alt_index] >=1e-7]*100\
                  ,temp_depol[rs_radar.Backscatter[:,plot_alt_index] >1e-6]*100]  
                ,['c','k','r']
                ,['*','*','*']
                ,[5,5,5]
                ,['None']
                ,[0,0,0]
                ,[]
                ,'lower right'
                ,'deff prime'
                ,'mm'
                ,'Depolarization'
                ,'%  '
                ,title
                ,[]
                ,[]
                ,[]
                ,[]
                ,display_defaults
                ,figs)
        
    if display_defaults.enabled('doppler_vs_deff_prime') \
            and hasattr(rs_radar,'MeanDopplerVelocity'):
        title = 'Doppler vs deff prime'
        temp_vel = rs_radar.MeanDopplerVelocity[:,layer_indices].copy()
        temp_deff_prime = rs.effective_diameter_prime[:,layer_indices].copy()

        gt.plot_xy('doppler_velocity_vs_deff_prime %4.2f-->%4.2f km,'\
                   %(usealts[layer_indices[0]],usealts[layer_indices[-1]])
                ,''
                ,usetimes
                ,[temp_deff_prime[mask[:,layer_indices] > 0.5]*1000.0,temp_deff_prime[mask[:,layer_indices] > 0.5]*1000.0]
                ,[temp_vel[mask[:,layer_indices] > 0.5],rs.rw_fall_velocity[mask[:,layer_indices] > 0.5]]
                ,['k','r']
                ,['*','.']
                ,[5,5]
                ,['None']
                ,[0,0]
                ,['radar','model']
                ,'lower right'
                ,'Eff diameter prime'
                ,'mm'
                ,'Doppler velocity'
                ,'m/s'
                ,title
                ,[]
                ,[]
                ,[]
                ,[]
                ,display_defaults
                ,figs)
   
    if display_defaults.enabled('p180_x_bs_vs_extinction')\
            and hasattr(rs,'p180_extinction') \
            and hasattr(rs,'extinction_aerosol'):
        title ='p180/4pi*bs, %4.2f-->%4.2f km' \
                %(usealts[layer_indices[0]]/1000.0,usealts[layer_indices[-1]]/1000.0)
        print 'rendering '+ title
        temp_ext = rs.extinction_aerosol[:,layer_indices]
        temp_p180_ext = rs.p180_extinction[:,layer_indices]
        x =  np.array([1e-6,1e-2])
        print 'p180'
        print rs.p180_extinction[:,plot_alt_index]
        print 'aer ext'
        print rs.extinction_aerosol
        gt.plot_2d_histogram('p180_x_bs_vs_extinction'
                          ,''
                          ,usetimes
                          ,[temp_ext,x,x]
                          ,[temp_p180_ext,x,x/2.0]
                          ,['-r','-g']          # colors and linetypes for optional overplot
                          ,None               #legend for optional overplot
                          ,'upper right'   # position for optional legend, eg. 'upper left'   
                          ,'extinction'
                          ,'1/m'
                          ,'backscat * 4pi/p180'
                          ,'1/m'
                          ,title
                          ,display_defaults
                          ,figs)
       
    
    if display_defaults.enabled('model_vs_doppler_velocity') \
            and hasattr(rs_radar,'MeanDopplerVelocity')\
            and hasattr(rs,'rw_fall_velocity'):

    
        t_vel = rs_radar.MeanDopplerVelocity.copy()
        t_model = rs.rw_fall_velocity.copy()
        t_vel[t_vel > 4.0] = np.NaN
        t_vel[t_vel < 0.05] = np.NaN
        t_model[t_model < 0.05] = np.NaN
        t_vel[mask < 0.5] = np. NaN
       
        t_vel[:,:layer_indices[0]] = np.NaN
        t_model[:,:layer_indices[0]] = np.NaN
        if layer_indices[-1]== len(t_model[0,:])-1:
            t_model[:,layer_indices[-1]:] = np.NaN
            t_vel[:,layer_indices[-1]:] = np.NaN 
        else:    
            t_model[:,(layer_indices[-1]+1):] = np.NaN
            t_vel[:,(layer_indices[-1]+1):] = np.NaN
        
        t_model[np.isnan(t_vel)] = np.NaN
        t_vel[np.isnan(t_model)] = np.NaN

        title = 'Model vs Doppler Velocity, %4.2f-->%4.2f km' \
                %(usealts[layer_indices[0]]/1000.0,usealts[layer_indices[-1]]/1000.0)
             
        if len(t_vel.shape) >1:
            new_size = t_vel.shape[0] * t_vel.shape[1]
            t_vel = np.reshape(t_vel,new_size)
            t_model = np.reshape(t_model,new_size)

        indices = np.arange(t_vel.shape[0])
        indices = indices[~np.isnan(t_vel)]
        if len(indices)==0:
            x_vec=[]
            pfit=[]
        else:
            x_vec = np.array([0, 4.0])
            pfit  = np.array([0, 4.0])
            #px = np.polyfit(t_vel[indices],t_model[indices],1)
            #x_vec = np.arange(50)*0.1
            #pfit=np.polyval(px,x_vec)
        #title = 'Doppler vs Model Velocity'
        print 'rendering '+ title
        gt.plot_2d_histogram('model_vs_doppler_velocity'
                          ,''
                          ,usetimes
                          ,[t_vel,x_vec]
                          ,[t_model,pfit]
                          ,['-r',]          # colors and linetypes for optional overplot
                          ,None               #legend for optional overplot
                          ,'upper right'   # position for optional legend, eg. 'upper left'   
                          ,'Doppler velocity'
                          ,'m/s'
                          ,'modeled velocity'
                          ,'m/s'
                          ,title
                          ,display_defaults
                          ,figs)
      
    if display_defaults.enabled('doppler_velocity_vs_deff_prime_hist') \
            and hasattr(rs_radar,'MeanDopplerVelocity'):

        temp_vel = rs_radar.MeanDopplerVelocity.copy()
        temp_deff_prime = rs.effective_diameter_prime.copy()

        #limit points plotted to altitudes between lo_lmt and hi_lmt
        temp_vel = rs_radar.MeanDopplerVelocity[:,layer_indices].copy()
        temp_deff_prime = rs.effective_diameter_prime[:,layer_indices].copy()
        temp_mask = mask[:,layer_indices]
       
        temp_mask[temp_deff_prime < 1.0e-5] = 0
        temp_mask[np.isnan(temp_deff_prime)] = 0
        temp_rw_fall = rs.rw_fall_velocity[:,layer_indices].copy()
        title= 'Doppler vel vs deff_prime, z=%4.2f-->%4.2f km, '\
                    %(usealts[layer_indices[0]]/1000.0,usealts[layer_indices[-1]]/1000.0)
        xlabel='Effective diameter prime'
        x_units = 'mm'
        ylabel='Doppler velocity'
        y_units = 'm/s'
        
        coef = particle_parameters['size_distribution']['alpha_water'] \
                   /particle_parameters['size_distribution']['g_water']
    
        legend_list = 'D^{%4.1f}exp(-%4.2f\cdot(D/D_m)^{%4.1f})'\
                      %(particle_parameters['size_distribution']['alpha_water'] \
                      ,coef,particle_parameters['size_distribution']['g_water'])
       
       
        legend_list = "$"+legend_list+"$"
        gt.plot_2d_histogram('doppler_velocity_vs_deff_prime_hist'
                          ,''
                          ,usetimes
                          ,[temp_deff_prime[temp_mask >0]*1000.0,temp_deff_prime[temp_mask >0]*1000.0]
                          ,[temp_vel[temp_mask > 0.5],temp_rw_fall[temp_mask > 0.5]]
                          ,['.r']          # colors and linetypes for optional overplot
                          ,[legend_list]     #legend for optional overplot
                          ,'upper right'   # position for optional legend, eg. 'upper left'   
                          ,xlabel
                          ,x_units
                          ,ylabel
                          ,y_units
                          ,title
                          ,display_defaults
                          ,figs)
    else:
        print 'no data for doppler vs deff_prime histogram'

                       
    if display_defaults.enabled('doppler_velocity_vs_dstar_hist') \
            and hasattr(rs_radar,'MeanDopplerVelocity') and hasattr(rs,'dstar'):
       
        temp_vel = rs_radar.MeanDopplerVelocity.copy()
        temp_deff_prime = rs.dstar.copy()

        #limit points plotted to altitudes between lo_lmt and hi_lmt
        temp_vel = rs_radar.MeanDopplerVelocity[:,layer_indices].copy()
        temp_dstar = rs.dstar[:,layer_indices].copy()
        temp_mask = mask[:,layer_indices]
       
        temp_mask[temp_dstar < 1.0e-5] = 0
        temp_mask[np.isnan(temp_dstar)] = 0
        temp_rw_fall = rs.rw_fall_velocity[:,layer_indices].copy()
        title= 'Doppler vel vs dstar, z=%4.2f-->%4.2f km, '\
                    %(usealts[layer_indices[0]]/1000.0,usealts[layer_indices[-1]]/1000.0)
        xlabel='dstar'
        x_units = 'mm'
        ylabel='Doppler velocity'
        y_units = 'm/s'
        
        coef = particle_parameters['size_distribution']['alpha_water'] \
                   /particle_parameters['size_distribution']['g_water']
    
        legend_list = 'D^{%4.1f}exp(-%4.2f\cdot(D/D_m)^{%4.1f})'\
                      %(particle_parameters['size_distribution']['alpha_water'] \
                      ,coef,particle_parameters['size_distribution']['g_water'])
       
       
        legend_list = "$"+legend_list+"$"
        gt.plot_2d_histogram('doppler_velocity_vs_dstar_hist'
                          ,''
                          ,usetimes
                          ,[temp_dstar[temp_mask >0]/1000.0,temp_dstar[temp_mask >0]/1000.0]
                          ,[temp_vel[temp_mask > 0.5],temp_rw_fall[temp_mask > 0.5]]
                          ,['.r']          # colors and linetypes for optional overplot
                          ,[legend_list]     #legend for optional overplot
                          ,'upper right'   # position for optional legend, eg. 'upper left'   
                          ,xlabel
                          ,x_units
                          ,ylabel
                          ,y_units
                          ,title
                          ,display_defaults
                          ,figs)
    else:
        print 'no data for doppler vs dstar histogram'
        
    if display_defaults.enabled('spectral_width_vs_deff_prime_hist') \
            and hasattr(rs_radar,'SpectralWidth'):
        temp_spectral_width = rs_radar.SpectralWidth.copy()
        temp_deff_prime = rs.effective_diameter_prime.copy()
        
        #limit points plotted to altitudes between lo_lmt and hi_lmt
        temp_spectral_width = rs_radar.SpectralWidth[:,layer_indices].copy()
        temp_deff_prime = rs.effective_diameter_prime[:,layer_indices].copy()    
        temp_mask = mask[:,layer_indices]
        temp_model_width = rs.model_spectral_width[:,layer_indices].copy()
        title= 'Spectral width vs deff_prime, z=%4.2f-->%4.2f km, '\
                    %(usealts[layer_indices[0]]/1000.0,usealts[layer_indices[-1]]/1000.0)
        
        xlabel='Effective diameter prime'

        
        coef = particle_parameters['size_distribution']['alpha_water'] \
                   /particle_parameters['size_distribution']['g_water']
    
        legend_list = 'D^{%4.1f}exp(-%4.1f\cdot(D/D_m)^{%4.1f})'\
                      %(particle_parameters['size_distribution']['alpha_water'] \
                      ,coef,particle_parameters['size_distribution']['g_water'])
           
        legend_list = "$"+legend_list+"$"
        legend_position = 'upper right'
        x_units = 'mm'
        ylabel='Spectral width'
        y_units = 'm/s'
        gt.plot_2d_histogram('spectral_width_vs_deff_prime_hist'
                          ,''
                          ,usetimes
                          ,[temp_deff_prime[temp_mask >0]*1000.0,temp_deff_prime[temp_mask>0.5]*1000.0]
                          ,[temp_spectral_width[temp_mask > 0.5],temp_model_width[temp_mask > 0.5]]
                          ,['.r']
                          ,[legend_list]
                          ,'upper right'
                          ,xlabel
                          ,x_units
                          ,ylabel
                          ,y_units
                          ,title
                          ,display_defaults
                          ,figs)
    else:
        print 'no data for radar spectral width vs deff_prime histogram'
   
    if display_defaults.enabled('p180_vs_deff_prime_hist') \
           and hasattr(rs,'effective_diameter_prime')\
           and hasattr(rs_inv,'beta_a_backscat'):
      
       p180=rs_inv.p180[:,layer_indices].copy()
       temp_deff_prime = rs.effective_diameter_prime[:,layer_indices].copy()
               
       print 'qc_mask applied to p180 deff prime histogram plot' 
       p180[rs_inv.qc_mask[:,layer_indices] == 0]=np.NaN
       temp_deff_prime[rs_inv.qc_mask[:,layer_indices] == 0] = np.NaN
       
       #title= 'P180 vs deff prime, z=%4.2f-->%4.2f km, '\
       #             %(usealts[layer_indices[0]]/1000.0,usealts[layer_indices[-1]]/1000.0)
       title = 'P180 vs deff prime--fix me'

       xlabel='Effective diameter prime'
       y_units = '1/sr'
       x_units = 'm'
       ylabel='P(180)/4$\\pi$'
       gt.plot_2d_histogram('p180_vs_deff_prime_hist'
                          ,' '
                          ,usetimes
                          ,[temp_deff_prime]
                          ,[p180]
                          ,None                            #no optional second plot
                          ,None                            #no legend for second plot
                          ,None                            #no position for legend      
                          ,xlabel
                          ,x_units
                          ,ylabel
                          ,y_units
                          ,title
                          ,display_defaults
                          ,figs)                     
    else:
           print 'no data for p180 vs effective diameter histogram'
           
    if display_defaults.enabled('spectral_width_vs_time') \
       and hasattr(rs,'model_spectral_width')\
       and rs_radar!=None and hasattr(rs_radar,'SpectralWidth'):
        lines = []
        widths = []
        colors = []
        legends = []
        # all radar measued spectal widths plot
        [lines,colors,widths,legends]=define_multiple_alt_plots(
                    rs_radar.SpectralWidth
                    ,usealts
                    ,plot_alt_index
                    ,lines
                    ,widths
                    ,colors
                    ,legends
                    ,1
                    ,'SW_all'
                    ,alt_color_list=True)
        #unmasked portion of the radar measured spectral widths
        [lines,colors,widths,legends]=define_multiple_alt_plots(
                    rs_radar.SpectralWidth * mask
                    ,usealts
                    ,plot_alt_index
                    ,lines
                    ,widths
                    ,colors
                    ,legends
                    ,2
                    ,'SW_good'
                    ,alt_color_list=True)
        #model calculated spectral widths
        [lines,colors,widths,legends]=define_multiple_alt_plots(
                    rs.model_spectral_width
                    ,usealts
                    ,plot_alt_index
                    ,lines
                    ,widths
                    ,colors
                    ,legends
                    ,2
                    ,'model')   

        title_str = 'spectral_width'
        print 'rendering ',title_str
        

        try:
            gt.plot_vs_time('spectral_width_vs_time'
                  ,'Radar'
                  ,usetimes
                  ,lines  
                  ,colors
                  ,widths
                  ,legends
                  ,'upper left'    
                  ,'Spectral Width'
                  ,'m/s'
                  ,title_str   
                  ,False #FIXME
                  ,display_defaults
                  ,figs)
            
	except AttributeError:
	    raise RuntimeError, "plot_vs_time - no data for spectral width "
    if display_defaults.enabled('hsrl_radar_liquid_water_vs_time') and hasattr(rs,'LWC'):

                lines = []
                widths = []
                colors = []
                legends = []
                # line width to 3 for effective_diameter_prime plot
                [lines,colors,widths,legends]=define_multiple_alt_plots(
                    rs.LWC * 1e3  #convert from kg/m^3 to gr/m^3
                    ,usealts
                    ,plot_alt_index
                    ,lines
                    ,widths
                    ,colors
                    ,legends)
                
                title_str = 'liquid water'
                print 'rendering ',title_str
                  
                gt.plot_vs_time('hsrl_radar_liquid_water_vs_time'
                  ,' '
                  ,usetimes
                  ,lines
                  ,colors
                  ,widths
                  ,legends
                  ,None              
                  ,instrument+' liquid water '
                  ,'gr/m^3'
                  ,title_str   
                  ,False #FIXME
                  ,display_defaults
                  ,figs)
    if display_defaults.enabled('eff_diameter_vs_time') and hasattr(rs,'effective_diameter'):

                lines = []
                widths = []
                colors = []
                legends = []
                # line width to 3 for effective_diameter_prime plot
                [lines,colors,widths,legends]=define_multiple_alt_plots(
                    rs.effective_diameter * 1000  #convert from meters to mm
                    ,usealts
                    ,plot_alt_index
                    ,lines
                    ,widths
                    ,colors
                    ,legends)
                
                title_str = 'effective diameter'
                print 'rendering ',title_str
                  
                gt.plot_vs_time('eff_diameter_vs_time'
                  ,' '
                  ,usetimes
                  ,lines
                  ,colors
                  ,widths
                  ,legends
                  ,None              
                  ,instrument+' eff diameter '
                  ,'mm'
                  ,title_str   
                  ,False #FIXME
                  ,display_defaults
                  ,figs)
    if display_defaults.enabled('mode_diameter_vs_time') and hasattr(rs,'mode_diameter'):

                lines = []
                widths = []
                colors = []
                legends = []
                # line width to 3 for effective_diameter_prime plot
                [lines,colors,widths,legends]=define_multiple_alt_plots(
                    rs.mode_diameter * 1000  #convert from meters to mm
                    ,usealts
                    ,plot_alt_index
                    ,lines
                    ,widths
                    ,colors
                    ,legends)
                
                title_str = 'Mode diameter'
                print 'rendering ',title_str
                  
                gt.plot_vs_time('mode_diameter_vs_time'
                  ,' '
                  ,usetimes
                  ,lines
                  ,colors
                  ,widths
                  ,legends
                  ,None              
                  ,instrument+' mode diameter '
                  ,'mm'
                  ,title_str   
                  ,False #FIXME
                  ,display_defaults
                  ,figs)
        
    if display_defaults.enabled('Z-R_rainfall_image') and hasattr(rs_radar,'Reflectivity'):
        title_str = 'Z-R rainfall'
        ZR_rainfall = ((10**(rs_radar.Reflectivity/10.0))/200)**(5.0/8.0)
        if 1: #try:
	    gt.rti_fig('Z-R_rainfall_image'
                    , instrument
                    , ZR_rainfall
                    , usetimes
                    , usealts
                    , title_str
                    , 'mm/hour'
                    , mask           
                    , display_defaults
                    ,figs)
	if 0: #except AttributeError:
	    raise RuntimeError, "show_images - no data for Z-R rainfall plot"
    if display_defaults.enabled('hsrl_radar_precip_rate_vs_time') and hasattr(rs,'hsrl_radar_dv_precip_rate'):
                print 'rendering hsrl_radar_precip_rate_vs_time'
                lines = []
                widths = []
                colors = []
                legends = []
                # precip rate--note conversion from m/s to mm/hour
                [lines,colors,widths,legends]=define_multiple_alt_plots(
                    rs.hsrl_radar_dv_precip_rate * 3600 * 1e3
                    ,usealts
                    ,plot_alt_index
                    ,lines
                    ,widths
                    ,colors
                    ,legends
                    ,1
                    ,' Doppler * LWC')
                [lines,colors,widths,legends]=define_multiple_alt_plots(
                    rs.hsrl_radar_precip_rate * 3600 * 1e3
                    ,usealts
                    ,plot_alt_index
                    ,lines
                    ,widths
                    ,colors
                    ,legends
                    ,2
                    ,' Model_V * LWC')
                

                if (display_defaults.get_value('hsrl_radar_precip_rate_vs_time','insitu_source')=='both_pars'\
                      or display_defaults.get_value('hsrl_radar_precip_rate_vs_time','insitu_source')=='pars2S1') \
                      and not rs_pars2S1 == None:
                    print
                    print 'rs_pars'
                    for i in range(len(rs_pars2S1.precip_rate)):
                        print usetimes[i],rs_pars2S1.precip_rate[i]
                    print    
                    lines.append(rs_pars2S1.precip_rate) 
                    legends.append('pars2S1')
                    colors.append('m')
                    widths.append(4)
                if  (display_defaults.get_value('hsrl_radar_precip_rate_vs_time','insitu_source')=='both_pars'\
                      or display_defaults.get_value('hsrl_radar_precip_rate_vs_time','insitu_source')=='pars2S2') \
                      and not rs_pars2S2 == None:   
                    lines.append(rs_pars2S2.precip_rate)
                    legends.append('pars2S2')
                    colors.append('c')
                    widths.append(4)
                title_str = 'precip rate'
                print 'rendering ',title_str
                gt.plot_vs_time('hsrl_radar_precip_rate_vs_time'
                  ,instrument
                  ,usetimes
                  ,lines
                  ,colors
                  ,widths
                  ,legends
                  ,'upper left'    
                  ,instrument+' precip rate '
                  ,'mm/hr'
                  ,title_str   
                  ,False #FIXME
                  ,display_defaults
                  ,figs)
    if display_defaults.enabled('ZR_accumulated_precip') and hasattr(rs,'hsrl_radar_precip_rate'):
                title_str ='ZR accumulated precipitation'
                print 'rendering '+title_str
                #ZR
                ZR_rainfall = ((10**(rs_radar.Reflectivity/10.0))/200)**(5.0/8.0) #mm/hour
                ZR_rainfall[np.isnan(ZR_rainfall)]=0.0
                ZR_accumulated = ZR_rainfall * (np.ones_like(rs.hsrl_radar_precip_rate)*rs.delta_t[:,np.newaxis])/3600.0
                ZR_accumulated_1 = np.cumsum(ZR_accumulated,0)
    
                ZR_rainfall = ((10**(rs_radar.Reflectivity/10.0))/88.0)**(1.0/1.5) #mm/hour
                ZR_rainfall[np.isnan(ZR_rainfall)]=0.0
                ZR_accumulated = ZR_rainfall * (np.ones_like(rs.hsrl_radar_precip_rate)*rs.delta_t[:,np.newaxis])/3600.0
                ZR_accumulated_2 = np.cumsum(ZR_accumulated,0)
                
                #hsrl_radar
                temp = (rs.hsrl_radar_precip_rate * (np.ones_like(rs.hsrl_radar_precip_rate)*rs.delta_t[:,np.newaxis])).copy()
                temp[np.isnan(temp)]=0.0
                precip = np.cumsum(temp,0)
                lines = []
                widths = []
                colors = []
                legends = []
                [lines,colors,widths,legends]=define_multiple_alt_plots(
                     precip * 1000.0             #m to mm
                     ,usealts
                     ,plot_alt_index
                     ,lines
                     ,widths
                     ,colors
                     ,legends
                     ,2
                     ,' hsrl_radar')
                
                [lines,colors,widths,legends]=define_multiple_alt_plots(
                     ZR_accumulated_1           # mm
                     ,usealts
                     ,plot_alt_index
                     ,lines
                     ,widths
                     ,colors
                     ,legends
                     ,1
                     ,' Z=200R^0.625'
                     ,alt_color_list=True)
                [lines,colors,widths,legends]=define_multiple_alt_plots(
                     ZR_accumulated_2           # mm
                     ,usealts
                     ,plot_alt_index
                     ,lines
                     ,widths
                     ,colors
                     ,legends
                     ,1
                     ,' Z=88R^0.66'
                     )
                gt.plot_vs_time('ZR_accumulated_precip'
                  ,instrument
                  ,usetimes
                  ,lines
                  ,colors
                  ,widths
                  ,legends
                  ,'upper left'    
                  ,'accumulated precip'
                  ,'mm'
                  ,title_str   
                 ,False #FIXME
                  ,display_defaults
                  ,figs)


    if display_defaults.enabled('hsrl_radar_accumulated_precip_vs_time') and hasattr(rs,'hsrl_radar_dv_precip_rate'):

                #multiply by the time interval between profiles in seconds yields precipitation in meters per profile
                temp = (rs.hsrl_radar_dv_precip_rate * ( np.ones_like(rs.hsrl_radar_dv_precip_rate)*rs.delta_t[:,np.newaxis])).copy()
                temp[np.isnan(temp)]=0.0
                dv_precip = np.cumsum(temp,0)
                temp = (rs.hsrl_radar_precip_rate * (np.ones_like(rs.hsrl_radar_precip_rate)*rs.delta_t[:,np.newaxis])).copy()
                temp[np.isnan(temp)]=0.0
                precip = np.cumsum(temp,0)

                lines = []
                widths = []
                colors = []
                legends = []
                # line width to 1 for this plot--plot all doppler velocities
                [lines,colors,widths,legends]=define_multiple_alt_plots(
                     precip * 1000.0             #m to mm
                     ,usealts
                     ,plot_alt_index
                     ,lines
                     ,widths
                     ,colors
                     ,legends
                     ,2
                     ,' model MWV * LWC')
              
                if not rs_pars2S1 == None \
                        and (display_defaults.get_value('hsrl_radar_accumulated_precip_vs_time'\
                              ,'insitu_source') == 'pars2S1'\
                        or display_defaults.get_value('hsrl_radar_accumulated_precip_vs_time'\
                               ,'insitu_source')=='both_pars'):
                   
                  
                    #rs_pars2S1.precip_rate supplied as mm/hour
                    temp = (rs_pars2S1.precip_rate * rs.delta_t/3600.0).copy()
                    temp[np.isnan(temp)]=0.0
                    accumulated_precip = np.cumsum(temp,0)
                    lines.append(accumulated_precip)
                    legends.append('pars2S1')
                    colors.append('k')
                    widths.append(4)

                    if not rs_pars2S2 == None \
                          and (display_defaults.get_value('hsrl_radar_accumulated_precip_vs_time'\
                              ,'insitu_source') == 'both_pars'):       
                        temp = (rs_pars2S2.precip_rate * rs.delta_t/3600.0).copy()
                        temp[np.isnan(temp)]=0.0
                        accumulated_precip = np.cumsum(temp,0)  
                        lines.append(accumulated_precip)  
                        legends.append('pars2S2')
                        colors.append('m')
                        widths.append(4)
                     
                if not rs_pars2S2 == None \
                        and display_defaults.get_value('hsrl_radar_accumulated_precip_vs_time'\
                        ,'insitu_source') == 'pars2S2':
                    temp = (rs_pars2S2.precip_rate * rs.delta_t/3600.0).copy()
                    temp[np.isnan(temp)]=0.0
                    accumulated_precip = np.cumsum(temp,0)  
                    lines.append(accumulated_precip)  
                    legends.append('pars2S2')
                    colors.append('m')
                    widths.append(4)
                          
                if not rs_marinemet == None \
                        and display_defaults.get_value('hsrl_radar_accumulated_precip_vs_time'\
                        ,'insitu_source') == 'wx1':
                    lines.append(rs_marinemet.rain_accumulated_wx1)
                    colors.append('g')
                    widths.append(2)
                    legends.append('wx2')
                if not rs_marinemet == None \
                        and display_defaults.get_value('hsrl_radar_accumulated_precip_vs_time'\
                        ,'insitu_source') == 'wx2':
                    lines.append(rs_marinemet.rain_accumulated_wx2)
                    widths.append(4)
                    colors.append('g')
                    legends.append('wx2')
                title_str =instrument + ' accumulated precip'
                gt.plot_vs_time('hsrl_radar_accumulated_precip_vs_time'
                  ,instrument
                  ,usetimes
                  ,lines
                  ,colors
                  ,widths
                  ,legends
                  ,'upper left'    
                  ,'accumulated precip'
                  ,'mm'
                  ,title_str   
                 ,False #FIXME
                  ,display_defaults
                  ,figs)
    if display_defaults.enabled('mode_diameter_histogram') and hasattr(rs,'mode_diameter'):
        max_diameter = display_defaults.get_value('mode_diameter_histogram','x max')
        min_diameter = display_defaults.get_value('mode_diameter_histogram','x min')
        nbins = display_defaults.get_value('mode_diameter_histogram','nbins')
        diameters = np.arange(nbins) * (max_diameter-min_diameter)/nbins \
                              + min_diameter/nbins
        units_str = display_defaults.get_value('mode_diameter_histogram','x_units')
        #extract the plot layer designated in display_defaults.json
        layer_dia = rs.mode_diameter[:,layer_indices].copy()
        layer_mask = mask[:,layer_indices].copy()
        good_layer_dia = layer_dia[layer_mask >0]
        layer_phase = rs.phase[:,layer_indices]
        good_layer_phase = layer_phase[layer_mask>0]
        if units_str == 'mm':
            layer_dia = layer_dia * 1000.0
            good_layer_dia = good_layer_dia * 1000.0
        elif units_str == 'microns':
            layer_dia = layer_dia * 1e6
            good_layer_dia = good_layer_dia * 1e6
        else:
            print 'WARNING----ERROR in mode_diameter_histogram x_units specification'
            
        if display_defaults.get_value('mode_diameter_histogram','phase')=='water':
            mode_dia_hist,jnk = np.array(np.histogram(layer_dia[layer_phase ==0],nbins
                                ,range=(np.min(diameters),np.max(diameters))))
            good_mode_dia_hist,jnk = np.array(np.histogram(good_layer_dia[good_layer_phase ==0],nbins
                                ,range=(np.min(diameters),np.max(diameters))))
        elif display_defaults.get_value('mode_diameter_histogram','phase')=='ice':
            mode_dia_hist,jnk = np.array(np.histogram(layer_dia[layer_phase ==1],nbins
                                ,range=(np.min(diameters),np.max(diameters))))
            good_mode_dia_hist,jnk = np.array(np.histogram(good_layer_dia[good_layer_phase ==1],nbins
                                ,range=(np.min(diameters),np.max(diameters))))
        elif display_defaults.get_value('mode_diameter_histogram','phase')=='all':
            mode_dia_hist,jnk = np.array(np.histogram(layer_dia,nbins
                                ,range=(np.min(diameters),np.max(diameters))))
            good_mode_dia_hist,jnk = np.array(np.histogram(good_layer_dia,nbins
                                ,range=(np.min(diameters),np.max(diameters))))
        else:
            raise RuntimeError('unrecognized phase--coop_display.py--mode_diameter_histogram')
       
            
        gt.plot_xy('mode_diameter_histogram'     #plot name
                        ,' '
                        ,rs.times
                        ,[diameters,diameters]
                        ,[mode_dia_hist,good_mode_dia_hist]
                        ,['b','r']
                        ,['*','+']
                        ,[5,5]
                        ,['-','-']
                        ,[2,2]
                        ,['N','N & mask']                      #legend
                        ,'upper right'
                        ,'Mode diameter '
                        ,units_str
                        ,'number'
                        ,' '
                        ,'Dm histogram ' + str(int(usealts[layer_indices[0]])) \
                                           + '-->' + str(int(usealts[layer_indices[-1]])) +'(m)' 
                        ,'' #text_str
                        ,None #'text_position_x
                        ,None #text_position_y
                        ,None #text_angle
                        ,display_defaults
                        ,figs)
    if display_defaults.enabled('model_fall_velocity_vs_mode_diameter') and hasattr(rs,'mode_diameter'):
        title_str = 'model fall velocity'
        print 'rendering '+title_str
        diameters = rs.mode_diameter[:,layer_indices]
        mw_velocities= rs.mw_fall_velocity[:,layer_indices]
        rw_velocities= rs.rw_fall_velocity[:,layer_indices]
        #gunn values sea level--unweighted
        d = np.array([.01,.02,.03,.04,.05,.06,.07,.08,.09,.1,.12,.14,.16,.18,.20])/100.0
        v = np.array([27 ,72 ,117,162,206,247,287,327,367,403,464,517,565,609,649])/100.0
        gt.plot_xy('model_fall_velocity_vs_mode_diameter'     #plot name
                        ,' '
                        ,rs.times
                        ,[diameters*1e3,diameters*1e3,d*1e3]   #convert m to mm
                        ,[mw_velocities,rw_velocities,v]
                        ,['r','r','m']
                        ,['*','*','*']
                        ,[5,5]                           #marker size
                        ,['None','None','-']                 #line style
                        ,[1,1,3]                           #line width
                        ,['mass weighted','radar weighed','Gunn sea level']     #legend
                        ,'lower right'
                        ,'Mode diameter '
                        ,units_str
                        ,'fall velocity'
                        ,'m/s'
                        ,'fall velocity ' + str(int(usealts[layer_indices[0]])) \
                                           + '-->' + str(int(usealts[layer_indices[-1]])) +'(m)' 
                        ,'' #text_str
                        ,None #'text_position_x
                        ,None #text_position_y
                        ,None #text_angle
                        ,display_defaults
                        ,figs)
        
    if display_defaults.enabled('model_particle_size_distributions'): # and size_dist != None:
        D = np.logspace(-2,1,200)
        Dm = [.005,.01,.025,.05,0.1,0.25,0.5,1.0]
        colors =['b','r','g','k','c','m','b','r']
        
        
        alpha = particle_parameters['size_distribution']['alpha_water']
        gam   = particle_parameters['size_distribution']['g_water']
        #exponent = alpha.copy()
            
      
        
        wp = display_defaults.get_value('model_particle_size_distributions','weight_by_power')
        x = []
        y = []
        c = []
        legend = []
        for i in range(len(Dm)):
            x.append(D)
            y_value = (D**(alpha+wp)*np.exp(-(alpha/gam)*(D/Dm[i])**gam))
            y_value = y_value/np.max(y_value)
            #y_value = size_dist.number(D,Dm,'water')
            y_value[y_value <= 0] = np.NaN
            y.append(y_value)
            legend.append('Dm='+str(Dm[i]))
            c.append(colors[i])
        gt.plot_xy('model_particle_size_distributions'  #plot name
                        ,' '
                        ,rs.times
                        ,x
                        ,y
                        ,c
                        ,[]
                        ,[]
                        ,['-']
                        ,[2]
                        ,legend
                        ,'upper right'
                        ,'particle size '
                        ,'mm'
                        ,'number * D**'+str(wp)
                        ,[]
                        ,'model size distributions'
                        ,'' #text_str
                        ,None #'text_position_x
                        ,None #text_position_y
                        ,None #text_angle
                        ,display_defaults
                        ,figs)    
    return

def show_mass_dimension_particle(instrument,display_defaults,rs,particle_parameters,usetimes,usealts,figs):
            #toplevel
            if usealts==None and hasattr(rs,'heights'):
                usealts=rs.heights
            if usetimes==None and hasattr(rs,'times'):
                usetimes=rs.times


            #EXAMPLE CODE COPIED FROM RADAR RENDERING
            print 'CALLED MASS DIMENSION PARTICLE RENDER'

      
            if display_defaults.enabled('effective_diameter_prime_image') and hasattr(rs,'effective_diameter_prime'):
                title_str = ' effective diameter prime '
      
                print 'rendering ',title_str
                try:
                    gt.rti_fig('effective_diameter_prime_image'
                    , instrument
                    , rs.effective_diameter_prime * 1e6
                    , usetimes
                    , usealts
                    , title_str
                    , 'microns'
                    , None            
                    , display_defaults
                    ,figs)
                except AttributeError:
                    raise RuntimeError, "show_images - no data for effective diameter prime plot"
           
            if display_defaults.enabled('effective_diameter_image') and hasattr(rs,'effective_diameter'):
                
                if particle_parameters['type'] == 'rain':
                     title_str = ' Effective Diameter, gamma dist($\\alpha$=%3.1f,$\\gamma$=%3.1f %s)' \
                                 %(particle_parameters['alpha_water'],particle_parameters['g_water'] \
                                 ,particle_parameters['type'])
                else:     
                    title_str = ' Effective Diameter, gamma dist($\\alpha$=%3.1f,$\\gamma$=%3.1f %s)' \
                                 %(particle_parameters['alpha_ice'],particle_parameters['g_ice'] \
                                 ,particle_parameters['type'])
                    
                print 'rendering ',title_str
                try:
                    gt.rti_fig('effective_diameter_image'
                    , instrument
                    , rs.effective_diameter * 1e6
                    , usetimes
                    , usealts
                    , title_str
                    , 'microns'
                    , None            
                    , display_defaults
                    ,figs)
                except AttributeError:
                    raise RuntimeError, "show_images - no data for effective diameter plot"
            if display_defaults.enabled('particle_number_density_image') and hasattr(rs,'num_particles'):
                if particle_parameters['type'] == 'rain':
                     title_str = ' Number Density, gamma dist($\\alpha$=%3.1f,$\\gamma$=%3.1f %s)' \
                                 %(particle_parameters['alpha_water'],particle_parameters['g_water'] \
                                 ,particle_parameters['type'])
                else:     
                    title_str = ' Number Density, gamma dist($\\alpha$=%3.1f,$\\gamma$=%3.1f %s)' \
                                 %(particle_parameters['alpha_ice'],particle_parameters['g_ice'] \
                                 ,particle_parameters['type'])
                print 'rendering ',title_str
                try:
                    gt.rti_fig('particle_number_density_image'
                    , instrument
                    , rs.num_particles  * 1e6
                    , usetimes
                    , usealts
                    , title_str
                    , '1/liter'
                    , None            
                    , display_defaults
                    ,figs)
                except AttributeError:
                    raise RuntimeError, "show_images - no data for number density plot"
            if display_defaults.enabled('liquid_water_content_image') and hasattr(rs,'LWC'):
            
                if particle_parameters['type'] == 'rain':
                     title_str = ' Liquid Water, gamma dist($\\alpha$=%3.1f,$\\gamma$=%3.1f %s)' \
                                 %(particle_parameters['alpha_water'],particle_parameters['g_water'] \
                                 ,particle_parameters['type'])
                else:     
                    title_str = ' Liquid Water, gamma dist($\\alpha$=%3.1f,$\\gamma$=%3.1f %s)' \
                                 %(particle_parameters['alpha_ice'],particle_parameters['g_ice'] \
                                 ,particle_parameters['type'])
                  
                print 'rendering ',title_str
                try:
                    gt.rti_fig('liquid_water_content_image'
                    , instrument
                    , rs.LWC   * 1e3
                    , usetimes
                    , usealts
                    , title_str
                    , '$g/m^3$'
                    , None            
                    , display_defaults
                    ,figs)
                except AttributeError:
                    raise RuntimeError, "show_images - no data for liquid water content plot"
               
            if display_defaults.enabled('mode_diameter_image') and hasattr(rs,'mode_diameter'):
                
                if particle_parameters['type'] == 'rain':
                     title_str = ' Mode Diameter, gamma dist($\\alpha$=%3.1f,$\\gamma$=%3.1f %s)' \
                                 %(particle_parameters['alpha_water'],particle_parameters['g_water'] \
                                 ,particle_parameters['type'])
                else:     
                    title_str = ' Mode Diameter, gamma dist($\\alpha$=%3.1f,$\\gamma$=%3.1f %s)' \
                                 %(particle_parameters['alpha_ice'],particle_parameters['g_ice'] \
                                 ,particle_parameters['type'])
                    
                print 'rendering ',title_str
                try:
                    gt.rti_fig('mode_diameter_image'
                    , instrument
                    , rs.mode_diameter 
                    , usetimes
                    , usealts
                    , title_str
                    , 'microns'
                    , None            
                    , display_defaults
                    ,figs)
                except AttributeError:
                    raise RuntimeError, "show_images - no data for mode diameter plot"

            if display_defaults.enabled('precip_rate_image')\
                      and hasattr(rs,'LWC')\
                      and hasattr(rs,'hsrl_radar_precip_rate'):
                
                if particle_parameters['type'] == 'rain':
                     title_str = 'Precip Rate, gamma dist($\\alpha$=%3.1f,$\\gamma$=%3.1f %s)' \
                                 %(particle_parameters['alpha_water'],particle_parameters['g_water'] \
                                 ,particle_parameters['type'])
                else:     
                    title_str = ' Precip Rate, gamma dist($\\alpha$=%3.1f,$\\gamma$=%3.1f %s)' \
                                 %(particle_parameters['alpha_ice'],particle_parameters['g_ice'] \
                                 ,particle_parameters['type'])

                print 'rendering ',title_str

                try:
                    gt.rti_fig('precip_rate_image'
                    , instrument
                    , rs.hsrl_radar_precip_rate 
                    , usetimes
                    , usealts
                    , title_str
                    , 'mm/hr'
                    , None
                    , display_defaults
                    ,figs)
                except AttributeError:
                    raise RuntimeError, "show_images - no data for precip rate"

            if display_defaults.enabled('hsrl_radar_precip_vs_time'):
               
                legend = ['z= %4.2f km' %(alt1)]
                gt.plot_vs_time('hsrl_radar_precip_vs_time'
                  ,instrument
                  ,usetimes
                  ,[rs.hsrl_radar_precip_rate[:,plot_alt_index]]  
                  ,['r']
                  ,[2]
                  ,legend
                  ,'upper left'    
                  ,instrument+' precip rate (mm/hr)'
                  ,[]
                  ,instrument+' precip rate vs time'   
                  ,False #FIXME
                  ,display_defaults
                  ,figs)
                                
          #FIXME add line plots
            """
  
  %make average particle diameter profile figure
  figure(fig_offset+509)
    clf
    %find valid data points 
    mask=zeros(size(data.effective_diameter_prime));
    mask(data.effective_diameter_prime>0 & data.effective_diameter_prime<1000)=1;
    eff_dia_prime=data.effective_diameter_prime;
    eff_dia_prime(mask==0)=0;
    mean_eff_dia_prime=sum(eff_dia_prime,1)./sum(mask,1);
    mean_diameter=data.mean_diameter;
    mask=zeros(size(mean_diameter));
    mask(~isnan(mean_diameter))=1;
    mean_diameter(mask==0)=0;
    mean_mean_dia=sum(mean_diameter,1)./sum(mask,1);
    h509=plot(mean_eff_dia_prime,dcal.range/1000,'r',mean_eff_dia,dcal.range/1000,'k',mean_mean_dia,dcal.range/1000,'b');
    set(h509([1 2 3]),'LineWidth',3);
    legend('d_e_f_fprime',['d_e_f_f(\alpha=2,\gamma=1,',particle_type]...
          ,['mean d, \alpha=2,\gamma=1,',particle_type]);
    set(gcf,'name','particle diameter profile')
       
    labl=['Time averaged particle diameter  ',datestr(start_time_num,1),' '];
  
    if end_time_num-start_time_num< 1/25 %less than 1 hr
      labl=[labl, datestr(start_time_num,13),'--> ',datestr(end_time_num,13)];
    elseif floor(start_time_num)==floor(end_time_num)  %same day
      labl=[labl, datestr(start_time_num,15),'--> ',datestr(end_time_num,15)];
    else
      labl=[labl, datestr(start_time_num,15),'-->',datestr(end_time_num,1)...
    ,' ',datestr(end_time_num,15)];
    end  
    title(labl);
    
    xlabel('Diameter (microns)');
    ylabel('Altitude (km)');
    ax=axis;
    xmax1=max(mean_eff_dia_prime);
    xmax2=max(mean_mean_dia);
    xmax=min(1000,max(xmax1,xmax2));
    axis([0 xmax ax(3) ax(4)]);
    grid
    end    
"""
            # end graphing code. don't modify the below

