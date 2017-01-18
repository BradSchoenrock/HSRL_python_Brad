#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import hsrl.filters.savitzky_golay as sg

import lg_base.core.decoratortools as nt

try:
    from bottleneck import allnan, nanmean
except ImportError:
    def allnan(x):
        return np.all(np.isnan(x))       


def filtered_extinction(inv,msl_altitudes,min_alt,window_od, t_window_length
                        ,max_z_window_length,order,adaptive,telescope_pointing=None):
    """filtered_extinction(inv,msl_altitudes,telescope_pointing, min_alt
                        window_od,max_window_length,order)
       Stravitzky-Golay fitting which can use a fixed window or
       an adaptive widow which restricts the length of the
       fit to a user specified optical depth interval estimated from the
       integrated backscatter cross section divided by a p180/4pi = 0.025
       inv                = dictionary containing beta_a_backscat_par, beta_a_backscat_perp
                            and beta_r_backscat, the Rayliegh backscatter
       msl_altitudes      = vector of bin altitudes (m)
       min_alt            = smooth alititudes > (min_alt + max_window_length/2.0)
                             (this is ignored  in current code)
       window_od          = max estimated od within fit window
       t_window_length    = length of time window (seconds)
       max_z_window_length= max length of altitude fit window (m)
       order              = order of polynomial to use in fit
       adaptive           = 0, fixed length window
                          = 1, window length can not exceed od limit"""

    if not hasattr(inv,'Nm'):
      print 'ERROR: XXXXXXXXXXXXXXXXXX Filtered extinction missing Nm in inv. Nothing to do.'
      return inv

    data_len = len(inv.Nm[0,:])    
    ntimes = len(inv.Nm[:,0])
    
    #distance between altitude bins
    dz=0
    for dzi in range(len(msl_altitudes)-1):
      ndz = msl_altitudes[dzi+1]-msl_altitudes[dzi]
      if ndz>0 and abs(dz-ndz)<.01:
        break
      dz=ndz
    z_window_pts = int(max_z_window_length/dz)
    #must be at least order +1
    if z_window_pts < order +1:
        z_window_pts = order +1
    #must be odd     
    if z_window_pts%2==0:
        z_window_pts = z_window_pts + 1

    if ntimes > order: 
        dt = inv.delta_t.copy()#(inv.times[2] - inv.times[1]).seconds
        dt[dt<1.0] = 1.0
        t_window_pts = np.array(t_window_length/dt,dtype='int')
        #must be at least order +1
        t_window_pts[t_window_pts < order +1]=order +1
        #must be odd     
        t_window_pts[t_window_pts%2==0]+=1

        filtered_Nm = inv.Nm.copy()
        integrated_backscat = inv.integrated_backscat.copy()
       
        #filter in time
        for k in range(data_len-1):
            if not allnan(filtered_Nm[:,k]) and  ntimes > t_window_pts:
               filtered_Nm[:,k] = sg.savitzky_golay(filtered_Nm[:,k],t_window_pts[k],order)
    else:
        filtered_Nm = inv.Nm.copy()
        integrated_backscat = inv.integrated_backscat.copy()

    if telescope_pointing is None :
        #if not provided assume zenith pointing
        t_pointing = np.ones_like(inv.Nm[:,0])
    else:    
        t_pointing = telescope_pointing.copy()
        t_pointing[t_pointing < 0.1] = -1.0
    extinction = np.zeros(filtered_Nm.shape)
    
    #if (data_len-min_alt/dz)  <  z_window_pts*5.0:
    if data_len  <  z_window_pts*5.0:
        print
        print 'WARNING---filtered_extinction--filter window length too long'
        print 'reseting z_window_pts from ',z_window_pts,
        #z_window_pts = int((data_len-min_alt/dz)/5.0)
        z_window_pts = int(data_len/5.0)
        z_window_pts = 2*(z_window_pts/2)+1
        print 'to ', z_window_pts
        if z_window_pts < order +2:
            print
            #raise ValueError, 'number of altitude resolution elements two few for filter length'
            print 'WARNING-----Savitzky_golay---number of altitude resolution elements two few for filter length'
            print '            inv.extintiction returned as NaN'
            print
            inv.extinction = np.NaN * inv.Nm
            return inv
        print
    #start_pt = int(np.ceil(z_window_pts/2 + min_alt/dz))
    start_pt = int(np.ceil(z_window_pts/2))
    end_pt   = data_len -z_window_pts/2 -1
    
    if adaptive == 0:
        slope_Nm = np.zeros_like(filtered_Nm[0,:])
        inv.p180 = np.zeros_like(integrated_backscat)
                             
        dbeta_dr = np.zeros_like(filtered_Nm[0,:])
        dbeta_dr[1:] =  + 0.5*(1/inv.beta_r_backscat[1:])\
                    *(inv.beta_r_backscat[1:]
                    -inv.beta_r_backscat[:-1])/dz 
        for i in range(ntimes):
            #compute extinction from the first derivative of a filtered Nm
            slope_Nm[start_pt:end_pt] = -sg.savitzky_golay(
                  filtered_Nm[i,start_pt:end_pt],z_window_pts,order,deriv = 1) / dz
            extinction[i,start_pt:end_pt]= \
               (-0.5*(1/filtered_Nm[i,start_pt:end_pt])
               *slope_Nm[start_pt:end_pt] 
               +dbeta_dr[start_pt:end_pt])*t_pointing[i]
        if 0:
          import matplotlib.pyplot as plt
          bin_vec = np.arange(len(filtered_Nm[0,:]))
          bin_vec2 = np.arange(len(slope_Nm))
          plt.figure(777777)
          plt.plot(slope_Nm,bin_vec2,'b',extinction[0,:],bin_vec2,'r',filtered_Nm[0,:],bin_vec2,'g'
                 ,-dbeta_dr,bin_vec2,'k',inv.p180[0,:],bin_vec2,'c')
          ax=plt.gca()
          ax.set_xscale('log')
          ax.grid(True)
         
    #when adaptive ==1,  window length derived from integrated backscat
    else:    
        #pick a intermediate value of p180/4pi for od estimate   
        p180 = 0.025

        #half of the maximum fit window in bins
        max_half_win =int((max_z_window_length/dz)/2)    
        wind_od = window_od/2.0
        
            
        #reflect data at end for padding
        yy=np.zeros((ntimes,data_len+max_half_win+1))
        end_range = range(data_len-1,(data_len-max_half_win-1),-1)
        yy[0:ntimes,0:data_len]=inv.Nm[0:ntimes,0:data_len]           
        yy[0:ntimes,data_len:data_len+max_half_win] = inv.Nm[0:ntimes,end_range]
        filtered_Nm = inv.Nm.copy()
        #extinction = np.zeros_like(inv.Nm)
        extinction = np.zeros(inv.Nm.shape)
        low_limit = np.zeros(data_len)
        high_limit = np.zeros_like(low_limit)
        beta_a=np.zeros(data_len+max_half_win+1)
    
        for k in range(ntimes):
            beta_a[0:data_len]=(inv.beta_a_backscat_par[k,0:data_len]\
                      +inv.beta_a_backscat_perp[k,0:data_len])
            beta_a[range((data_len-max_half_win),data_len)] \
                  =(inv.beta_a_backscat_par[k,end_range] \
                        +inv.beta_a_backscat_perp[k,end_range])
    
            #compute integrated backscatter
            beta_a[np.isnan(beta_a)]=0
            int_bs_od=np.cumsum(beta_a)*dz/(2*p180)
            int_bs_od[data_len:data_len+max_half_win+1]=int_bs_od[data_len]
            #find end points of fit for each data point
            #use the optical depth estimated from the integrated backscatter to compute
            #limits for polynomial fit at each data point
            derivative_coefs=np.arange(order,0,-1)
            order_lmt = order/2 +1
                        
            for i in range(start_pt,data_len):
            
                lo_lmt=np.max([i-max_half_win,start_pt])
                while (int_bs_od[i] -int_bs_od[lo_lmt] > wind_od) and i-lo_lmt > order_lmt:
                   lo_lmt = lo_lmt + 1

                hi_lmt = i+max_half_win
                while (int_bs_od[hi_lmt] - int_bs_od[i] > wind_od) and hi_lmt-i > order_lmt:
                    hi_lmt = hi_lmt -1
                
                #make fitting interval symetric around i
                if i - lo_lmt < hi_lmt-i :
                    hi_lmt = 2*i -lo_lmt
            
                elif i - lo_lmt > hi_lmt -i :
                   lo_lmt = 2*i - hi_lmt
            
    

                ylocal = yy[k,lo_lmt:hi_lmt+1]
                x=np.arange(len(ylocal))
                #print x
                pc = np.polyfit(x,ylocal,order)
                filtered_Nm[k,i]= np.polyval(pc,i-lo_lmt)
    
                slope_Nm = np.polyval(derivative_coefs*pc[0:order],i-lo_lmt)/dz

                #print 'ext',i,lo_lmt,hi_lmt,ylocal,slope_Nm,pc,pc[0:order]\
                #         ,dz,np.polyval(pc,[0,1,2,3]),t_pointing[k]
       
                extinction[k,i]=(-0.5*(1/filtered_Nm[k,i])*slope_Nm \
                    + 0.5*(1/inv.beta_r_backscat[i])\
                    *(inv.beta_r_backscat[i]-inv.beta_r_backscat[i-1])/dz)\
                    *t_pointing[k]

                #print 'ext',i,lo_lmt,hi_lmt,extinction[k,i],slope_Nm\
                #      ,(inv.beta_r_backscat[i]-inv.beta_r_backscat[i-1])/dz\
                #      ,inv.beta_r_backscat[i],inv.beta_r_backscat[i-1],t_pointing[k]
                      
                xx=np.arange(len(x)*10)/10.0
           
       
            
    inv.extinction = type(filtered_Nm)(extinction)
    time_vec = np.ones_like(inv.extinction[:,0])
    beta_r_array = 8 * np.pi * (time_vec[:,np.newaxis] * inv.beta_r_backscat[np.newaxis,:])/3.0 
    inv.extinction_aerosol = inv.extinction - beta_r_array

    #compute p180 from integrated backscatter and optical depth over window segments
    inv.p180 = np.zeros_like(inv.extinction_aerosol)
    start_pt = int(np.ceil(z_window_pts/2))
    delta_tau = np.zeros_like(inv.extinction)
   
    delta_tau[:,start_pt:] = -0.5 * np.log(
        inv.Nm[:,start_pt:]/inv.Nm[:,:-start_pt]\
        *(beta_r_array[:,:-start_pt]/beta_r_array[:,start_pt:]))
    
    delta_tau[:,start_pt:] = delta_tau[:,start_pt:] \
          - (beta_r_array[:,start_pt:]+beta_r_array[:,:-start_pt]) * msl_altitudes[start_pt]/2.0
    
    inv.p180[:,start_pt:] = (\
        integrated_backscat[:,start_pt:] - integrated_backscat[:,:-start_pt])\
        /delta_tau[:,start_pt:]
   
    return inv       
