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

#for plotting cooperate product of rs_mean and raman photon counts (merge) in profile
def show_raman_hsrl_profile(display_defaults,rs,parameters,usetimes,usealts,figs,hsrl=None,raman=None):
    """
       plot results created by the  hsrl_raman_profile stream
    """
    if 1:
        print
        print
        print 'entering show_raman_hsrl_profile-------------------------------------------'
        print 
        print dir(rs)

                         
    if display_defaults.enabled('rl_hsrl_backscatter_profile'):
                print 'Rendering rl_hsrl_backscatter_profile'
                
                if 1: #try:
                   gt.plot_vs_altitude('rl_hsrl_backscatter_profile'                     
                    ,'rl_hsrl'                   
                    ,rs.hsrl_profile.times                      
                    ,rs.hsrl_profile.msl_altitudes                    
                    ,[rs.hsrl_profile.inv.beta_a_backscat[0,:]   
                         ,rs.raman_profile.inv.beta_a_backscat[0,:]
                      ,rs.hsrl_profile.inv.beta_a_1064_backscat[0,:]
                         ,rs.raman_profile.inv.beta_r_355*3.0/(8.0*np.pi)
                         ,rs.hsrl_profile.inv.beta_r_backscat
                         ,rs.hsrl_profile.inv.beta_r_1064*3.0/(8.0*np.pi)]
                    ,['g','m','r','b','b','b']                   
                    ,[2,2,2,1,2,3]               
                    ,['aerosol532','aerosol355','aerosol_1064','mol355','mol532','mol1064']
                    ,'lower left'                   
                    ,'Aerosol Backscatter cross section '                  
                    ,'1/(m sr)'                         
                    ,'rl_hsrl_Backscatter'            
                    ,False    #clear figure                
                    ,display_defaults                    
                    ,figs)
                
                if 0: #except: #AttributeError:
                    #raise
                    raise RuntimeError, "raman_display.py - error plotting raman backscatter profile"
                                    
    if display_defaults.enabled('rl_hsrl_extinction_profile'):
        
        if 1:  #try:
                   gt.plot_vs_altitude('rl_hsrl_extinction_profile'                     
                    ,''                   
                    ,rs.hsrl_profile.times         #usetimes                        
                    ,rs.hsrl_profile.msl_altitudes                    
                    ,[rs.hsrl_profile.inv.extinction_aerosol[0,:]
                       ,rs.raman_profile.inv.extinction_aerosol[0,:]]
                    ,['g','m']                   
                    ,[2,2,2,1,2,3]                              #linewidth list
                    ,['532nm','355nm']                 #legend list, None = no legend
                    ,'lower left'                              #legend position eg. 'upper left'
                    ,'Aersol extinction cross section'        #xlabel          
                    ,'1/(m)'                         #x units
                    ,'rl_hsrl_extinction'                      #string in title 
                    ,False                              #clear figure before plotting               
                    ,display_defaults                    
                    ,figs)
                
        if 0: #except: #AttributeError:
                    #raise
                    raise RuntimeError, "show_images - error plotting hsrl and Raman backscat_profiles"


#for plotting cooperate product of rs_mean and raman photon counts (merge)
def show_ramanmerge_hsrl(display_defaults,rs,parameters,usetimes,usealts,figs,hsrl=None,raman=None):
    pass


#for plotting cooperateive product of rs_inv and raman inverted data (dep or ext)
def show_raman_hsrl(display_defaults,rs,parameters,usetimes,usealts,figs,hsrl=None,raman=None):
    """
       show_raman_hsrl(display_defaults,rs,parameters,usetimes,usealts,figs,hsrl=None,raman=None)
       plot variables from hsrl rs_inv data stream with those from raman_inv data stream
    """
    if 0:
        print 'entering show_raman_hsrl-----------------------------------------------'
        print 'raman structure contains:'
        print dir(rs.raman_inv)
        print 'hsrl structure contains:'
        print dir(rs.rs_inv)
        print
  
        

   
    return
