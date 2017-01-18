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

#any of the profile sections can be None. never assume availability without checking first. any not named here would be found in kwargs
def show_all_profiles(display_defaults,figs, wholeframe,rlprofaerosol_profile=None,  profiles=None, raman_rawprofile=None, raman_profile=None, raw_profiles=None,*args,**kwargs):
    #all frame parts are in kwargs dictionary
    import time

    print
    print
    print
    print '****************************************************************************'
    print 'entering show_all_profiles'
    print
    print
    print 'alprofiles artist has ',kwargs.keys()
    print

    if 0: #not rlprofaerosol_profile is None:     
        print 'rlprofaerosol'
        print dir(rlprofaerosol_profile)
        print   
    if 0: #not profiles is None :
        print 'profiles'
        print dir(profiles)
        print
        print dir(profiles.inv)
        print
    if 0: #not raman_rawprofile is None:
        print 'raman_rawprofile'
        print dir(raman_rawprofile)
        print
    if 0: #not raw_profiles is None:
        print 'raw_profiles'
        print dir(raw_profiles)
        print
       
    if 0: #display_defaults.enabled('hsrl_thor_extinction_profiles') and  kwargs.has_key('rlproffex1thor_profile'):
        for name in kwargs:
            print 'name',kwargs
        print '************************************************************************ '   
        for name in kwargs:
            print 'name',name
            print dir(kwargs[name])
            if name == 'rlproffex1thor_profile':
               print 
               print 'rlproffex1thor_profile'
               print 
            elif name == 'raman_hsrl_profile':
               print
               print 'raman_hsrl_profile'
               print dir(kwargs[name])
               print
               print    'kwargs[name].raman_profile'
               print dir(kwargs[name].raman_profile)
               print
               print    'kwargs[name].hsrl_profile'
               print dir(kwargs[name].hsrl_profile)
               print
               print       'kwargs[name].hsrl_profile.inv'
               print dir(kwargs[name].hsrl_profile.inv)
               print
               print 'dir(kwargs[name].raman_profile.inv)'
               print dir(kwargs[name].raman_profile.inv)
      
    if display_defaults.enabled('hsrl_thor_extinction_profiles'):
                   name =  'rlproffex1thor_profile'
                   gt.plot_vs_altitude('hsrl_thor_extinction_profiles'                     
                    ,'hsrl_thor'                   
                    ,profiles.times         #usetimes                        
                    ,kwargs[name].altitudes  #usealts                    
                    ,[kwargs[name].extinction_be[0,:]
                      ,kwargs['raman_hsrl_profile'].raman_profile.inv.extinction_aerosol[0,:]
                      ,profiles.inv.extinction_aerosol[0,:]]
                    ,['m','r','g']                   
                    ,[2,2,2]                              #linewidth list
                    ,['thor_be_355','uw_ext_355','hsrl_532']                 #legend list, None = no legend
                    ,'lower left'                              #legend position eg. 'upper left'
                    ,'Aersol extinction cross section'        #xlabel          
                    ,'1/(m)'                         #x units
                    ,'rl_extinction'                      #string in title 
                    ,False                              #clear figure before plotting               
                    ,display_defaults                    
                    ,figs)
                
    if 0: #except: #AttributeError:
                    #raise
                    raise RuntimeError, "show_images - error plotting  and Raman extinction_profiles"
                
    if display_defaults.enabled('aero_thor_extinction_profiles'):
                   thor =  'rlproffex1thor_profile'
                   gt.plot_vs_altitude('aero_thor_extinction_profiles'                     
                    ,'aero_thor'                   
                    ,profiles.times         #usetimes                        
                    ,kwargs[thor].altitudes  #usealts                    
                    ,[kwargs[thor].extinction_be[0,:]
                      ,rlprofaerosol_profile.extinction_n2[0,:]]
                    ,['m','r']                   
                    ,[2,2]                              #linewidth list
                    ,['thor_be','rlprofaerosol']                 #legend list, None = no legend
                    ,'lower left'                              #legend position eg. 'upper left'
                    ,'Aersol extinction cross section'        #xlabel          
                    ,'1/(m)'                         #x units
                    ,'rl_extinction'                      #string in title 
                    ,False                              #clear figure before plotting               
                    ,display_defaults                    
                    ,figs)
                
    if 0: #except: #AttributeError:
                    #raise
                    raise RuntimeError, "show_images - error plotting  and Raman extinction_profiles"
                
    if display_defaults.enabled('hsrl_thor_backscatter_profiles'):
                   name =  'rlproffex1thor_profile'
                   gt.plot_vs_altitude('hsrl_thor_backscatter_profiles'                     
                    ,'hsrl_thor'                   
                    ,profiles.times         #usetimes                        
                    ,kwargs[name].altitudes  #usealts                    
                    ,[kwargs[name].particulate_backscatter_be_no_smooth[0,:]
                      ,kwargs['raman_hsrl_profile'].raman_profile.inv.beta_a_backscat[0,:]
                      ,profiles.inv.beta_a_backscat[0,:]
                      ,profiles.inv.beta_a_1064_backscat[0,:]]
                    ,['m','r','g','r']                   
                    ,[2,2,2,2]                              #linewidth list
                    ,['thor_355','uw__355','hsrl_355','hsrl_a1064']                 #legend list, None = no legend
                    ,'lower left'                              #legend position eg. 'upper left'
                    ,'Aersol backscatter cross section'        #xlabel          
                    ,'1/(m sr)'                         #x units
                    ,'backscatter'                      #string in title 
                    ,False                              #clear figure before plotting               
                    ,display_defaults                    
                    ,figs)
    if display_defaults.enabled('thor_backscatter_profiles'):
                   thor =  'rlproffex1thor_profile'
                   gt.plot_vs_altitude('thor_backscatter_profiles'                     
                    ,'thor'                   
                    ,profiles.times         #usetimes                        
                    ,kwargs[thor].altitudes  #usealts                    
                    ,[kwargs[thor].particulate_backscatter_e_n2[0,:]
                     ,kwargs[thor].particulate_backscatter_e[0,:]
                     ,kwargs[thor].particulate_backscatter_be_no_smooth[0,:]
                     ,kwargs[thor].particulate_backscatter_be[0,:] 
                     ,kwargs[thor].particulate_backscatter_e_beS[0,:]
                     ,rlprofaerosol_profile.beta[0,:]]
                    ,['m','r','g','c','k','b']                   
                    ,[2,2,2,2,2,2]                              #linewidth list
                    ,['e_n2','e','be_noS','be','e_be_S','aer']                 #legend list, None = no legend
                    ,'lower left'                              #legend position eg. 'upper left'
                    ,'Aersol backscatter cross section'        #xlabel          
                    ,'1/(m sr)'                         #x units
                    ,'backscatter'                      #string in title 
                    ,False                              #clear figure before plotting               
                    ,display_defaults                    
                    ,figs)
                   
    if 0: #except: #AttributeError:
                    #raise
                    raise RuntimeError, "show_images - error plotting hsrl and Raman backscat_profiles"
                
    if display_defaults.enabled('rl_hsrl_backscat_ratio_profile'):
                print 'Rendering rl_hsrl_backscat_ratio_profile'
                name = 'raman_hsrl_profile'
                print kwargs[name].hsrl_profile.inv.beta_r_backscat.shape
                if 1: #try:
                   gt.plot_vs_altitude('rl_hsrl_backscat_ratio_profile'                     
                    ,'rl_hsrl'                   
                    ,profiles.times                      
                    ,kwargs[name].raman_profile.altitudes                   
                    ,[kwargs[name].raman_profile.inv.aerosol_backscatter_ratio_low[0,:]   
                         ,kwargs[name].raman_profile.inv.aerosol_backscatter_ratio_high[0,:]
                    ,kwargs[name].raman_profile.inv.aerosol_backscatter_ratio[0,:]
                    ,kwargs[name].hsrl_profile.inv.beta_a_backscat[0,:]
                              /kwargs[name].hsrl_profile.inv.beta_r_backscat]
                      
                    ,['r','m','k','g']                   
                    ,[2,2,2,2]               
                    ,['low_355','high_355','merge_355','hsrl_532']
                    ,'lower left'                   
                    ,'Aerosol/molecular backscatter ratio '                  
                    ,None                         
                    ,'backscatter_ratios'            
                    ,False    #clear figure                
                    ,display_defaults                    
                    ,figs)
                
                if 0: #except: #AttributeError:
                    #raise
                    raise RuntimeError, "raman_display.py - error plotting raman-hsrl scattering ratio profiles"
              
                
    time.sleep(10)
    return


