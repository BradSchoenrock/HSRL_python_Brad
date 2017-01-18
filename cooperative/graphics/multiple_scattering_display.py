import numpy as np
import lg_base.graphics.graphics_toolkit as gt
# Try to use the much faster nanmean from bottleneck, otherwise fall back
# to the scipy.stats version

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


def show_multiple_scattering(instrument,display_defaults,rs_multiple_scattering,multiple_scatter_parameters,particle_parameters,usetimes,usealts,figs,profiles=None):
    print 'entering show_multiple_scattering()'
    
    if display_defaults.enabled('multiple_scattering_ratios') \
	   and hasattr(rs_multiple_scattering,'ms_ratios_profile'):
        print 'Rendering mulitiple_scattering_ratios'
	[nalts,norders] =rs_multiple_scattering.ms_ratios_profile.shape
        legend_list =['2','3','4','5','6']
        linewidth_list = [1,1,1,1,1,1,1]
        color_list = ['g','r','c','m','k','g','y'] 
        legend_list = legend_list[:norders-2]
        linewidth_list = linewidth_list[:norders-2]
        color_list = color_list[:norders-2]
	plot_list =[]
	for i in range(2,norders):
            plot_list.append(rs_multiple_scattering.ms_ratios_profile[:,i])
        if norders>2:
            plot_list.append(nansum(rs_multiple_scattering.ms_ratios_profile[:,2:],1))
	    legend_list.append('total')
            linewidth_list.append(3)
            color_list.append('b')
	gt.plot_vs_altitude('multiple_scattering_ratios'       
                 ,instrument                   
                 ,usetimes                        
                 ,rs_multiple_scattering.msl_altitudes                   
                 ,plot_list
                 ,color_list               
                 ,linewidth_list
                 ,legend_list                
                 ,'upper right'          
                 ,'multiple/single scatter '                  
                 ,[]                         
                 ,'Multiple scatter'            
                 ,None              
                 ,display_defaults                    
                 ,figs)
   
    if display_defaults.enabled('multiple_scatter_weighted_diameter') \
	   and hasattr(rs_multiple_scattering,'weighted_diameter'):
	instrument = ''
	print 'Rendering---multiple scatter weighted diameter'    
	gt.plot_vs_altitude('multiple_scatter_weighted_diameter'       
                 ,instrument                   
                 ,usetimes
                 ,rs_multiple_scattering.msl_altitudes
                 ,[rs_multiple_scattering.weighted_diameter[0,:]*1e3]
                 ,['b']               
                 ,[2]
                 ,[]                
                 ,'upper right'          
                 ,'particle diameter '                  
                 ,'mm'                         
                 ,'water, ice weighted diameter'            
                 ,None            
                 ,display_defaults                    
                 ,figs)
        
    if display_defaults.enabled('multiple_scatter_extinction_profile') \
	   and hasattr(rs_multiple_scattering,'extinction_profile'):
	instrument = ''
	print 'Rendering---multiple scatter extinction profile'    
	gt.plot_vs_altitude('multiple_scatter_extinction_profile'       
                 ,'ms'                  
                 ,usetimes
                 ,rs_multiple_scattering.msl_altitudes
                 ,[rs_multiple_scattering.extinction_profile[0,:]]
                 ,['b']               
                 ,[2]
                 ,[]                
                 ,'upper right'          
                 ,'extinction profile for ms '                  
                 ,'1/m'                         
                 ,'extinction'            
                 ,None            
                 ,display_defaults                    
                 ,figs)
    if display_defaults.enabled('ms_color_ratio_profile') \
	   and hasattr(rs_multiple_scattering,'ms_ratios_profile_2'):
	instrument = ''
	print 'Rendering---multiple scatter color_ratio profile'
        ms_ratio_total_profile = nansum(rs_multiple_scattering.ms_ratios_profile[:,2:],1)
        ms_ratio_total_profile_2 = nansum(rs_multiple_scattering.ms_ratios_profile_2[:,2:],1) 
	gt.plot_vs_altitude('ms_color_ratio_profile'       
                 ,instrument                   
                 ,usetimes
                 ,rs_multiple_scattering.msl_altitudes
                 ,[ms_ratio_total_profile_2 / ms_ratio_total_profile]
                 ,['b']               
                 ,[2]
                 ,[]                
                 ,'upper right'          
                 ,str(rs_multiple_scattering.wavelength) + '/' + str(rs_multiple_scattering.second_wavelength) + 'ratio'                  
                 ,''                         
                 ,'multiple scatter color ratio'            
                 ,None            
                 ,display_defaults                    
                 ,figs)
        
    if display_defaults.enabled('multiple_scatter_od_correction') \
	   and hasattr(rs_multiple_scattering,'extinction_profile')\
           and hasattr(rs_multiple_scattering,'ms_ratios_profile')\
           and hasattr(profiles,'inv'):
        
        start_index = len(rs_multiple_scattering.ranges[rs_multiple_scattering.ranges < 0.1])
       
        beta = rs_multiple_scattering.extinction_profile[0,start_index:].copy()
        beta[np.isnan(beta)]=0.0
        dr =rs_multiple_scattering.ranges[start_index+1] \
                -rs_multiple_scattering.ranges[start_index]
        beta[np.isnan(beta)]=0.0
        od = np.cumsum(beta*dr)

        ms_ratio_total_profile = nansum(rs_multiple_scattering.ms_ratios_profile[:,2:],1)
        apparent_od = od - np.log(1.0 + ms_ratio_total_profile[start_index:])/2.0
	instrument = ''

        if particle_parameters['type']=='rain' \
                  or particle_parameters['type'] == 'water_cloud':
              p180 = particle_parameters['p180_water']
              legend =['beta/'+str(p180),'with ms','measured od']
        elif particle_parameters['type'] == 'ice':
              p180 = particle_parameters['p180_ice']
              legend =['beta/'+str(p180),'with ms','measured od']
        else:
              legend = ['beta/p180','with ms','measured od']
	print 'Rendering---multiple scatter od correction'
        
	gt.plot_vs_altitude('multiple_scatter_od_correction'       
                 ,instrument                   
                 ,usetimes
                 ,rs_multiple_scattering.msl_altitudes[start_index:]
                 ,[od,apparent_od,profiles.inv.optical_depth_aerosol[0,start_index:]]
                 ,['b','g','k']               
                 ,[2,2,2]
                 ,legend                
                 ,'lower right'          
                 ,'optical depth '                  
                 ,''                         
                 ,'optical depth corrections'            
                 ,None              
                 ,display_defaults                    
                 ,figs)
    if display_defaults.enabled('multiple_scatter_od_correction_2') \
	   and hasattr(rs_multiple_scattering,'extinction_profile')\
           and hasattr(rs_multiple_scattering,'ms_ratios_profile_2'):

        start_index = len(rs_multiple_scattering.ranges[rs_multiple_scattering.ranges < 0.1])
       
        beta = rs_multiple_scattering.extinction_profile[0,start_index:].copy()
        beta[np.isnan(beta)]=0.0
      
        dr =rs_multiple_scattering.ranges[start_index+1] \
                -rs_multiple_scattering.ranges[start_index]
        beta[np.isnan(beta)]=0.0
       
        od = np.cumsum(beta*dr)
        ms_ratio_total_profile_2 = nansum(rs_multiple_scattering.ms_ratios_profile_2[:,2:],1)
        apparent_od = od - np.log(1.0 + ms_ratio_total_profile_2[start_index:])/2.0
	instrument = ''
	print 'Rendering---multiple scatter od correction for second wavelength'
       
	gt.plot_vs_altitude('multiple_scatter_od_correction_2'       
                 ,instrument                   
                 ,usetimes
                 ,rs_multiple_scattering.msl_altitudes[start_index:]
                 ,[od,apparent_od]
                 ,['b','g']               
                 ,[2,2]
                 ,['beta/p180 od','+ ms']                
                 ,'lower right'          
                 ,'optical depth--'+str(rs_multiple_scattering.second_wavelength *1e9)                  
                 ,'nm'                         
                 ,'optical depth corrections'            
                 ,None              
                 ,display_defaults                    
                 ,figs)   
    
    if display_defaults.enabled('ms_ratio_image') \
           and hasattr(rs_multiple_scattering,'ms_ratio_total')\
           and rs_multiple_scattering.ms_ratio_total.shape[0]>10 :
        print 'Rendering ms_ratio_image'
        array_size = rs_multiple_scattering.ms_ratio_total.shape
        N2 = multiple_scatter_parameters['highest_order']
        title_str = 'multiple/single scatter up to order '+ str(N2) 
        
        try:
            gt.rti_fig('ms_ratio_image'
                    ,'ms'
                    , rs_multiple_scattering.ms_ratio_total
                    , rs_multiple_scattering.times 
                    , rs_multiple_scattering.msl_altitudes
                    , title_str
                    , ''
                    , None             
                    , display_defaults
                    ,figs)
        except AttributeError:
            raise RuntimeError, "show_images - no data for multiple scatter ratios plot"
    
