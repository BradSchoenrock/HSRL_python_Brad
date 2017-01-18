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



def show_raman_inverted_profile(instrument,display_defaults,rs,timerange,usealts,figs,consts=None):
    pass

    
def show_raman_profile(instrument,display_defaults,rs,timerange,usealts,figs,consts=None):
    """
       show_raman_profile(instrument,display_defaults,rs,timerange,usealts,figs)
       Generates plots from raman_profile, raman_profile.inv  and
       raman_rawprofile streams.
       These are averaged over the entire requested time interval
    """
       
    if usealts is None and hasattr(rs,'altitudes'):
                usealts=rs.altitudes
    #if usetimes is None and hasattr(rs,'times'):
    usetimes=np.array([rs.start_time,rs.end_time])
    print
    print
    print 'rl_raw_count_profiles'
    print 'dir(rs)'
    print dir(rs)
    print
    if display_defaults.enabled('rl_polarization_count_profiles') \
                   and hasattr(rs,'depolarization_counts_high'):

    
        gt.plot_vs_altitude('rl_polarization_count_profiles'
                    ,instrument               
                    ,usetimes                        
                    ,usealts                    
                    ,[rs.sum_elastic_counts_high[0,:]/rs.sum_mean_shots
                      ,rs.sum_depolarization_counts_high[0,:]/rs.sum_mean_shots
                      ,rs.sum_elastic_counts[0,:]/rs.sum_mean_shots]
                    ,['b','c','r']                   
                    ,[1,2,3,1,2,3,1,2,3]
                    ,['elastic_hi_par','elastic_hi_perp','elastic_total']          #legend list        
                    ,'upper right'                 #legend position
                    ,'Count rate'       #x lablel           
                    ,'Hz'                         
                    ,'rl_polarization_counts'        #title str    
                    ,False              #clear figure                
                    ,display_defaults                    
                    ,figs)
        
    if display_defaults.enabled('rl_raw_count_profiles') and hasattr(rs,'nitrogen_counts'):
   
                print 'Rendering rl_raw_count_profiles'
                if 1: #try:
                   k_n2  = consts['nitrogen_hi_to_lo_gain_ratio'][1]
                   k_e   = consts['elastic_hi_to_lo_gain_ratio'][1]
                   k_h2o = consts['water_hi_to_lo_gain_ratio'][1]
                 
                   gt.plot_vs_altitude('rl_raw_count_profiles'
                    ,instrument               
                    ,usetimes                        
                    ,usealts                    
                    ,[k_n2 * rs.sum_nitrogen_counts_low[0,:]/rs.sum_mean_shots
                      ,rs.sum_nitrogen_counts_high[0,:]/rs.sum_mean_shots
                      ,rs.sum_nitrogen_counts[0,:]/rs.sum_mean_shots
                      ,k_e * rs.sum_elastic_counts_low[0,:]/rs.sum_mean_shots
                      ,rs.sum_elastic_counts_high[0,:]/rs.sum_mean_shots
                      ,rs.sum_elastic_counts[0,:]/rs.sum_mean_shots
                      ,k_h2o * rs.sum_water_counts_low[0,:]/rs.sum_mean_shots
                      ,rs.sum_water_counts_high[0,:]/rs.sum_mean_shots
                      ,rs.sum_water_counts[0,:]/rs.sum_mean_shots]
                    ,['b','c','b','r','c','r','g','c','g']                   
                    ,[1,2,3,1,2,3,1,2,3]
                    ,['n2_lo*gain','n2_hi','n2_comp','el_lo*gain','el_hi','el_comp'
                            ,'h2o_lo*gain','h2o_hi','h2o_comp']          #legend list        
                    ,'lower left'                 #legend position
                    ,'Count rate'       #x lablel           
                    ,'Hz'                         
                    ,'rl_composite_counts'        #title str    
                    ,False              #clear figure                
                    ,display_defaults                    
                    ,figs)
                
                if 0: #except: #AttributeError:
                    #raise
                    raise RuntimeError, 'raman_display.py - error plotting merge_count_profile'


      

    if display_defaults.enabled('rl_merge_count_profiles') and\
               hasattr(rs,'nitrogen_counts'):
        print 'Rendering rl_merge_count_profiles'
        
        if 1: #try:
            gt.plot_vs_altitude('rl_merge_count_profiles'
                            ,instrument
                            ,usetimes
                            ,usealts
                            ,[rs.sum_nitrogen_counts[0,:]/rs.sum_mean_shots
                            ,rs.sum_elastic_counts[0,:]/rs.sum_mean_shots
                            ,rs.sum_water_counts[0,:]/rs.sum_mean_shots]
                            ,['b','r','g']
                            ,[2,2,2]
                            ,['N2','elastic','h20']
                            ,'lower left' ,
                            'rl_count rate'
                            ,'MHz'
                            ,'rl_counts'
                            ,False #clear figure
                            ,display_defaults
                            ,figs)
                
        if 0: #except: #AttributeError:
                    #raise
                    raise RuntimeError, "raman_display.py - error plotting merge count profile"
  
    if display_defaults.enabled('rl_backscatter_profile') and hasattr(rs,'inv'):
                print 'Rendering rl_backscatter_profile'
                
                if 1: #try:
                   gt.plot_vs_altitude('rl_backscatter_profile'                     
                    ,instrument                   
                    ,usetimes                        
                    ,usealts                    
                    ,[rs.inv.beta_r_355 * 3.0 / (8*np.pi),rs.inv.beta_a_backscat[0,:]]
                    ,['b','r']                   
                    ,[2,2]               
                    ,['mol355','aerosol']
                    ,'lower left'                   
                    ,'Aerosol Backscatter cross section '                  
                    ,'1/(m sr)'                         
                    ,'rl_Backscatter'            
                    ,False    #clear figure                
                    ,display_defaults                    
                    ,figs)
                
                if 0: #except: #AttributeError:
                    #raise
                    raise RuntimeError, "raman_display.py - error plotting backscatter profile"

    if display_defaults.enabled('rl_scattering_ratio_profile') and hasattr(rs,'inv'):
                print 'Rendering rl_scattering_ratio_profile'
                
                if 1: #try:
                   gt.plot_vs_altitude('rl_scattering_ratio_profile'                     
                    ,instrument                   
                    ,usetimes                        
                    ,usealts                    
                    ,[rs.inv.aerosol_backscatter_ratio[0,:],rs.inv.aerosol_backscatter_ratio_low[0,:],rs.inv.aerosol_backscatter_ratio_high[0,:]]
                    ,['r','c','m']                   
                    ,[3,1,1]               
                    ,['SR','SR_lo','SR_hi']
                    ,'upper right'                   
                    ,'rl_Aerosol Backscatter Ratio'                  
                    ,None                       
                    ,'SR'            
                    ,False    #clear figure                
                    ,display_defaults                    
                    ,figs)
                
                if 0: #except: #AttributeError:
                    #raise
                    raise RuntimeError, "raman_display.py - error plotting merge_count_profile"

    if display_defaults.enabled('rl_extinction_profile') and hasattr(rs,'inv'):
       
        if 1:  #try:
                   gt.plot_vs_altitude('rl_extinction_profile'                     
                    ,'hsrl_rl'                   
                    ,rs.times         #usetimes                        
                    ,rs.altitudes  #usealts                    
                    ,[rs.inv.extinction_aerosol[0,:]
                       ,rs.inv.extinction[0,:]]
                    ,['m','c']                   
                    ,[2,2]                              #linewidth list
                    ,['aerosol_ext','total_ext']                 #legend list, None = no legend
                    ,'lower left'                              #legend position eg. 'upper left'
                    ,'Aersol extinction cross section'        #xlabel          
                    ,'1/(m)'                         #x units
                    ,'rl_extinction'                      #string in title 
                    ,False                              #clear figure before plotting               
                    ,display_defaults                    
                    ,figs)
                
        if 0: #except: #AttributeError:
                    #raise
                    raise RuntimeError, "show_images - error plotting hsrl and Raman extinction profiles"
        
    if display_defaults.enabled('rl_od_profile') and hasattr(rs,'inv') \
             and (not ('installation' in consts) \
                  or consts['installation'] == 'ground' \
                  or consts['installation'] == 'shipborne'\
                  or (hasattr(rs.profiles,'telescope_pointing') and
                    (rs.telescope_pointing[0] > 0.9 \
                  or rs.telescope_pointing[0] < 0.1))): 
        
        #Rayliegh backscatter phase function time altitude step size
        #pdr=(8*np.pi/3)*(rs.rs_mean.msl_altitudes[3]-rs.rs_mean.msl_altitudes[2])
        #pdr = pdr/np.cos(np.pi*rs_consts['telescope_zenith_angle']/180.0)
        #if no installation in calvals assume ground based instrument
        #or if airborne based the telescope always pointed up
        ref_alt = rs.altitudes[rs.inv.mol_norm_index]/1000.0
        xlabel = 'Optical depth (ref alt = '+ str(ref_alt)+' km)'
       
        if not('installation' in consts) \
                or (consts['installation'] == 'ground')\
                or (consts['installation'] == 'shipborne')\
                or rs.telescope_pointing[0] > 0.9:
            print 'ground based with zenith pointing telescope'
            
            dr = (rs.altitudes[1]-rs.altitudes[0])
            dr = dr/np.cos(np.pi*consts['telescope_roll_angle_offset']/180.0)
           
            od_at_norm_alt= \
               dr*np.sum(rs.inv.beta_r_355[0:rs.inv.mol_norm_index])
            
            if hasattr(rs.inv,'wfov_corr_optical_depth'):

                od_total=rs.inv.wfov_corr_optical_depth[0,:] 
                od_aerosol = rs.inv.wfov_corr_optical_depth_aerosol[0,:]
                xlabel = 'WFOV corrected OD (ref alt = '+ str(ref_alt)+' km)'
            else:
                od_total=rs.inv.optical_depth[0,:] 
                od_aerosol = rs.inv.optical_depth_aerosol[0,:]
                title = 'Optical depth'
                              
        #if airborne with telescope always pointed at nadir
        elif (consts['installation']=='airborne') and hasattr(rs,'inv') \
            and hasattr(rs,'telescope_pointing') and rs.telescope_pointing[0]<0.1:
                print 'airborne with telescope always pointed at nadir'
                mean_alt=rs.mean_GPS_MSL_Alt
                aircraft_alt_index = len(rs.altitudes[rs.altitudes<=mean_alt])
                #Rayliegh backscatter phase function time altitude step size
                pdr=-(8*np.pi/3)*shared_dz
                pdr = pdr/np.cos(np.pi*consts['telescope_roll_angle_offset']/180.0)
                mol_optical_depth = np.zeros_like(rs.inv.beta_r_backscat)
        
                #make reverse list of indices for integration looking downward
                indices = range(len(mol_optical_depth[:aircraft_alt_index])-1,-1,-1)
                mol_optical_depth[:aircraft_alt_index]= pdr*np.cumsum(rs.inv.beta_r_backscat[indices])
                mol_optical_depth=mol_optical_depth \
                      -mol_optical_depth[rs.inv.mol_norm_index]
                 
                od_aerosol=rs.inv.optical_depth[0,:] - mol_optical_depth
                od_total = rs.inv.optical_depth[0,:]

        #plot if ground based or if nearly all telescope up or down pointing
        if not consts['installation']=='airborne'\
            or  rs.telescope_pointing[0] <0.1 \
            or  rs.telescope_pointing[0] >0.9:
            
            print 'rendering raman optical depth profile'
            gt.plot_vs_altitude('rl_od_profile'        #display defaults plot name
                 ,instrument                       #instrument name
                 ,timerange                  #python datetimes vector
                 ,rs.altitudes         #altitude vector (meters)
                 ,[od_aerosol,od_total]     #variables
                 ,['r','b']                        #colors, [] default colors
                 ,[]                               #widths, [] sets widths = 2
                 ,['aerosol','total']                   #legend list, [] = no legend
                 ,'lower right'                    #legend position, [] ok if list []
                 ,xlabel                           #xlabel
                 ,[]                               #x units
                 ,'rl_optical depth'               #plot title
                 ,False                       # =1, clear figure before new plot
                 ,display_defaults                    
                 ,figs)
        else:
            print ' '
            print 'mixed up and down pointing no OD plot'
            print ' '

def show_ramanmerge(instrument,display_defaults,rs,usetimes,usealts,figs,consts=None):
    """
      show_ramanmerge(instrument,display_defaults,rs,usetimes,usealts,figs)
      show plots from rlprofmerge data stream
    """
   

    plot_alt_index = rs.plot_alt_index if hasattr(rs,'plot_alt_index') else 0
    plot_altitude = usealts[plot_alt_index]
    layer_indices = rs.layer_indices if hasattr(rs,'layer_indices') else np.arange(len(usealts))

 

    if display_defaults.enabled('rl_h2o_to_n2_count_ratio_image') and hasattr(rs,'water_counts'):

                title_str = instrument+' h20_to_n2_counts_ratio '
                print 'rendering ',title_str
                try:
                    gt.rti_fig('rl_h2o_to_n2_count_ratio_image'
                    , instrument
                    , rs.water_counts / rs.nitrogen_counts
                    , usetimes
                    , usealts
                    , title_str
                    , ' '
                    , None
                    , display_defaults
                    ,figs)
                except AttributeError:
                    raise RuntimeError, "show_images - no data for raman h20_to_n2_count_ratio  plot"

    if display_defaults.enabled('rl_hi_to_lo_channel_ratios') and \
        hasattr(rs,'nitrogen_counts_high'):

        print 'Rendering rl_hi_to_lo_channel ratios'
        title ='rl nfov to wfov vs time Z = %4.2f km' %(rs.altitudes[plot_alt_index]/1000.0)
        legend = ['n2 %3.1f' %( np.nanmean(rs.nitrogen_counts_high[:,plot_alt_index] \
                          / rs.nitrogen_counts_low[:,plot_alt_index]))]
        legend.append('elastic %3.1f' %(np.nanmean(rs.elastic_counts_high[:,plot_alt_index] \
                          / rs.elastic_counts_low[:,plot_alt_index])))
                                
        legend.append('h2o %3.1f' %(np.nanmean(rs.water_counts_high[:,plot_alt_index] \
                          / rs.water_counts_low[:,plot_alt_index])))                                
        gt.plot_vs_time('rl_hi_to_lo_channel_ratios'
                 ,instrument
                 ,usetimes
                 ,[rs.nitrogen_counts_high[:,plot_alt_index] \
                          / rs.nitrogen_counts_low[:,plot_alt_index]
                     ,rs.elastic_counts_high[:,plot_alt_index] \
                          / rs.elastic_counts_low[:,plot_alt_index]
                     ,rs.water_counts_high[:,plot_alt_index] \
                          / rs.water_counts_low[:,plot_alt_index]]
                 ,['b','r','g']
                 ,[2,2,2]
                 ,legend
                 ,'upper left'    
                 ,'narrow to wfov count ratios'
                 ,''
                 ,title  
                 ,False
                 ,display_defaults
                 ,figs)
        
    if display_defaults.enabled('rl_filter_mode') and \
        hasattr(rs,'filter_mode'):

        print 'Rendering Raman filter_mode'
        try:  
            title ='Raman filter mode'  
            gt.plot_vs_time('rl_filter_mode'
                 ,instrument
                 ,usetimes
                 ,[rs.filter_mode]
                 ,['b']
                 ,[2]
                 ,None
                 ,'upper left'    
                 ,'filter mode'
                 ,''
                 ,title  
                 ,False
                 ,display_defaults
                 ,figs)             
        except AttributeError:
             raise RuntimeError, "show_images - Raman filter mode plot"

    if display_defaults.enabled('rl_count_rates_vs_time') and \
        hasattr(rs,'nitrogen_counts'):

        print 'Rendering Raman count rates vs time'

        try:
            title ='Raman counting rates,Z=%4.2f km' %(rs.altitudes[plot_alt_index]/1000.0) 
            gt.plot_vs_time('rl_counts_vs_time'
                 ,instrument
                 ,usetimes
                 ,[rs.nitrogen_counts[:,plot_alt_index]
                       ,rs.elastic_counts[:,plot_alt_index]
                       ,rs.water_counts[:,plot_alt_index]]
                 ,['b','r','g']
                 ,[2,2,2]
                 ,['n2','elastic','h2o']
                 ,'upper left'    
                 ,'Count rate'
                 ,'Hz'
                 ,title  
                 ,False
                 ,display_defaults
                 ,figs)         
        except AttributeError:
             raise RuntimeError, "show_images - Raman counting rates plot"
         
    if display_defaults.enabled('rl_dark_counts_vs_time') and \
        hasattr(rs,'nitrogen_counts'):

        print 'Rendering Raman dark counts vs time'

        if 1: #try:
            #low channel is wfov, high is narrow fov
            title ='Raman dark count rates'
            temp_n2_low = rs.nitrogen_low_dark_counts.copy()
            temp_n2_low[temp_n2_low <= 0] = np.NaN
            temp_elastic_low = rs.elastic_low_dark_counts.copy()
            temp_elastic_low[temp_elastic_low <= 0] = np.NaN
            temp_water_low = rs.water_low_dark_counts.copy()
            temp_water_low[temp_water_low<=0] = np.NaN
            temp_n2_high = rs.nitrogen_high_dark_counts.copy()
            temp_n2_high[temp_n2_high <= 0] = np.NaN
            temp_elastic_high = rs.elastic_high_dark_counts.copy()
            temp_elastic_high[temp_elastic_high <= 0] = np.NaN
            temp_water_high = rs.water_high_dark_counts.copy()
            temp_water_high[temp_water_high<=0] = np.NaN
            gt.plot_vs_time('rl_dark_counts_vs_time'
                 ,instrument
                 ,usetimes
                 ,[temp_elastic_low, temp_n2_low, temp_water_low
                   ,temp_elastic_high,temp_n2_high,temp_water_high]
                 ,['r','b','g','r','b','g']
                 ,[1,1,1,3,3,3]
                 ,['elastic_low','N2_low','H2O_low'
                     ,'elastic_high','N2_high','H2O_high']
                 ,'upper left'    
                 ,'Count rate'
                 ,'Hz'
                 ,title  
                 ,False
                 ,display_defaults
                 ,figs)         
        if 0: #except AttributeError:
             raise RuntimeError, "show_images - Raman counting rates plot"
         
def show_raman(instrument,display_defaults,rs,usetimes,usealts,figs,consts=None,qc_mask=None):
    """
        show_raman(instrument,display_defaults,rs,usetimes,usealts,figs)
        Display variables produced by raman_inversion
    """
    
            #toplevel
    if usealts is None and hasattr(rs,'altitudes'):
                usealts=rs.altitudes
    if usetimes is None and hasattr(rs,'times'):
                usetimes=rs.times
     
    plot_alt_index = rs.plot_alt_index if hasattr(rs,'plot_alt_index') else 0
    plot_altitude = usealts[plot_alt_index]
    layer_indices = rs.layer_indices if hasattr(rs,'layer_indices') else np.arange(len(usealts))
    
         
    if display_defaults.enabled('rl_backscatter_image') and (hasattr(rs,'beta_a_backscat') or hasattr(rs,'beta') or hasattr(rs,'particulate_backscatter_be')):

                title_str = instrument+' backscatter cross section '
                print 'rendering ',title_str
                try:
                    gt.rti_fig('rl_backscatter_image'
                    , instrument
                    , rs.beta_a_backscat if hasattr(rs,'beta_a_backscat') else (rs.beta if hasattr(rs,'beta') else rs.particulate_backscatter_be)
                    , usetimes
                    , usealts
                    , title_str
                    , '1/(m sr)'
                    , None                 
                    , display_defaults
                    ,figs)
                except AttributeError:
                    raise RuntimeError, "show_images - no data for raman backscatter plot"
            
    if display_defaults.enabled('rl_extinction_image') and (hasattr(rs,'extinction_aerosol') or hasattr(rs,'extinction_n2') or hasattr(rs,'extinction_be')):

                title_str = instrument+' Extinction '
                print 'rendering ',title_str     
                try:
                    gt.rti_fig('rl_extinction_image'
                    , instrument
                    , rs.extinction_aerosol if hasattr(rs,'extinction_aerosol') else (rs.extinction_n2 if hasattr(rs,'extinction_n2') else rs.extinction_be)
                    , usetimes
                    , usealts
                    , title_str
                    , '1/m'
                    , None                             
                    , display_defaults
                    ,figs)
                except AttributeError:
                    raise RuntimeError, "show_images - no data for raman extinction plot"

    
    if display_defaults.enabled('rl_h2o_to_n2_count_ratio_image') and hasattr(rs,'water_counts'):

                title_str = instrument+' h20_to_n2_counts_ratio '
                print 'rendering ',title_str
                try:
                    gt.rti_fig('rl_h2o_to_n2_count_ratio_image'
                    , instrument
                    , rs.water_counts / rs.nitrogen_counts
                    , usetimes
                    , usealts
                    , title_str
                    , ' '
                    , None
                    , display_defaults
                    ,figs)
                except AttributeError:
                    raise RuntimeError, "show_images - no data for raman h20_to_n2_count_ratio  plot"  
   

    if display_defaults.enabled('rl_depol_image') and (hasattr(rs,'linear_depol') or hasattr(rs,'depolarization_ratio')):

                title_str = instrument+' Linear Depolarization '
                print 'rendering ',title_str     
                try:
                    gt.rti_fig('rl_depol_image'
                    , instrument
                    , rs.linear_depol if hasattr(rs,'linear_depol') else rs.depolarization_ratio
                    , usetimes
                    , usealts
                    , title_str
                    , '%'
                    , None                             
                    , display_defaults
                    ,figs)
                except AttributeError:
                    raise RuntimeError, "show_images - no data for linear depolarization plot"
            
    if display_defaults.enabled('rl_aerosol_scattering_ratio_image') and (hasattr(rs,'aerosol_backscatter_ratio') or hasattr(rs,'scattering_ratio_e_n2')):

                if hasattr(rs,'aerosol_backscatter_ratio'):
                    title_str = instrument+' Aerosol Scattering Ratio '
                else:
                    title_str = instrument+' Scattering Ratio '
                print 'rendering ',title_str     
                try:
                    gt.rti_fig('rl_aerosol_scattering_ratio_image'
                    , instrument
                    , rs.aerosol_backscatter_ratio if hasattr(rs,'aerosol_backscatter_ratio') else rs.scattering_ratio_e_n2
                    , usetimes
                    , usealts
                    , title_str
                    , '%'
                    , None                             
                    , display_defaults
                    ,figs)
                except AttributeError:
                    raise RuntimeError, "show_images - no data aerosol scattering ratio plot"
    
    if display_defaults.enabled('rl_raw_scattering_ratio_image') and hasattr(rs,'raw_scattering_ratio'):

                title_str = instrument+' Raw Scattering Ratio '
                print 'rendering ',title_str     
                try:
                    gt.rti_fig('rl_raw_scattering_ratio_image'
                    , instrument
                    , rs.raw_scattering_ratio
                    , usetimes
                    , usealts
                    , title_str
                    , '%'
                    , None                             
                    , display_defaults
                    ,figs)
                except AttributeError:
                    raise RuntimeError, "show_images - no data aerosol scattering ratio plot"

    if display_defaults.enabled('rl_od_vs_time') and hasattr(rs,'optical_depth_aerosol'):
        """
        alt_res = shared_dz
        alt1=float(display_defaults.get_value('od_vs_time','altitude1'))
        alt2=float(display_defaults.get_value('od_vs_time','altitude2'))
       
        alt1_index = int(alt1* 1000.0 / alt_res)
        alt2_index = int(alt2* 1000.0 / alt_res)
        nalts = rs.rs_inv.optical_depth_aerosol.shape[1]
        if alt1_index >= nalts:
            alt1_index = nalts - 1
        if alt2_index >= nalts:
            alt2_index = nalts - 1
        """    
        lines = []
        legend = []
        if display_defaults.get_value('rl_od_vs_time','use_plot_layer_alts'):
            alt1_index = layer_indices[0]
            alt2_index = layer_indices[-1]
        else:
             alt_res = rs.altitudes[1]-rs.altitudes[0]
             alt1=float(display_defaults.get_value('rl_od_vs_time','altitude1'))
             alt2=float(display_defaults.get_value('rl_od_vs_time','altitude2'))
       
             alt1_index = int(alt1* 1000.0 / alt_res)
             alt2_index = int(alt2* 1000.0 / alt_res)
             nalts = rs.optical_depth_aerosol.shape[1]
             if alt1_index >= nalts:
                alt1_index = nalts - 1
             if alt2_index >= nalts:
                alt2_index = nalts - 1
            
        if display_defaults.get_value('rl_od_vs_time','show_od_1_and_2'):
            if display_defaults.get_value('rl_od_vs_time','show_od_1_and_2')==1:
                lines.append(rs.optical_depth_aerosol[:, alt1_index])
                legend.append('Z= %4.2f' %(rs.altitudes[alt1_index]))
            else: 
                lines.append(rs.rs_inv.optical_depth_aerosol[:, alt2_index])
                legend.append('Z= %4.2f' %(rs.altitudes[alt2_index])) 
        if display_defaults.get_value('rl_od_vs_time','show_od_difference'):
            od = rs.optical_depth_aerosol[:, alt2_index]\
                          - rs.optical_depth_aerosol[:, alt1_index]
            if qc_mask is not None:
              od[qc_mask[:,alt1_index]==0] =np.NaN
            lines.append(od)
            legend.append('d_OD')
        print figs
        gt.plot_vs_time('rl_od_vs_time'
                 ,instrument
                 ,rs.times
                 ,lines
                 ,['g','b','r']
                 ,[2,2,2]
                 ,legend
                 ,'upper left'    
                 ,'Aerosol optical depth'
                 ,None
                 ,'Aerosol optical depth vs time'   
                 ,False
                 ,display_defaults
                 ,figs)        
               
    return figs
