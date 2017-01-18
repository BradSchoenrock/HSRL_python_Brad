import numpy as np
from scipy.optimize import leastsq, curve_fit
from scipy.interpolate import splrep, splev
import copy
#from scipy.signal import savgol_coeffs

from lg_dpl_toolbox.filters.dpl_rolling_window_filter import WindowedFilterDescription
import lg_base.core.array_utils as hau

import hsrl.filters.savitzky_golay as sg

import lg_base.core.decoratortools as nt

from datetime import datetime,timedelta
from hsrl.data_stream.processing_utilities import deltaz

try:
    from bottleneck import allnan,nanmin
except ImportError:
    def allnan(x):
        return np.all(np.isnan(x))       

def listMerge(timezlist):
  ret=None
  for x in timezlist:
    if ret is None:
      ret=copy.deepcopy(x)
    else:
      ret.append(x)
  return ret

class filtered_extinction(object):
  def __init__(self):
    self.ref={}
    self.priorTime=None

  def svm(self,window_size, order, deriv=-1):
    if not window_size in self.ref:
      self.ref[window_size]={}
    if not order in self.ref[window_size]:
      self.ref[window_size][order]={}
    if not deriv in self.ref[window_size][order]:
      self.ref[window_size][order][deriv]=None
    return self.ref[window_size][order][deriv]

  def setsvm(self,val,window_size, order, deriv=-1):
    self.ref[window_size][order][deriv]=val

  def __call__(self,_inv,min_alt,window_od
                        ,max_z_window_length,order,adaptive,wfov_corr_enabled=None):
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
       window_od          = max estimated od within fit window
       t_window_length    = length of time window (seconds)
       max_z_window_length= max length of altitude fit window (m)
       order              = order of polynomial to use in fit
       adaptive           = 0, fixed length window
                          = 1, window length can not exceed od limit"""

    #inv is a rolling time_window into the inv structure
    inv=listMerge(_inv)
    
    ntimes = inv.times.size
    print 'ext window ',inv.times[0],inv.times[-1],ntimes
    useindex=0
    if self.priorTime is None:
      self.priorTime=inv.times[0]-timedelta(seconds=10)
    while useindex<(ntimes-1) and self.priorTime>=inv.times[useindex]:
      useindex=useindex+1
    self.priorTime=inv.times[useindex]

    #r_inv is a structure containing a single profile selected at useindex of inv
    r_inv=copy.deepcopy(_inv[useindex])


    if not hasattr(r_inv,'Nm') or not hasattr(r_inv,'msl_altitudes'):
      print 'WARNING: Missing value for extinction calculation. skipped'
      print  hasattr(r_inv,'Nm'),hasattr(r_inv,'msl_altitudes')
      return r_inv
    data_len = len(r_inv.Nm[0,:])    
    msl_altitudes=r_inv.msl_altitudes
    #distance between altitude bins
    adz=deltaz(msl_altitudes)
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

    t_window_pts=ntimes
   
    t_pts=np.arange(len(inv.Nm[:,0]))

    #provide time domain filtering with Savitzky-Golay type filter
    #use polyfit/polyval since inv contains rolling window used for single time
    if order <> 1: 
        pc = np.polyfit(t_pts,inv.Nm,order)           
        filtered_molecular_counts = np.polyval(pc,useindex)
        pcw = np.polyfit(t_pts,inv.beta_a_backscat,order)
        filtered_beta_a_backscat = np.polyval(pcw,useindex)
    else: # simplified computation when fit is linear
        filtered_molecular_counts = np.nanmean(inv.Nm,0) 
        filtered_beta_a_backscat = np.nanmean(inv.beta_a_backscat,0)
          
    t_pointing = np.ones([ntimes])
    
    if hasattr(r_inv,'telescope_pointing'):
      t_pointing[:]=np.NAN
      for i in range(ntimes):
          if inv.telescope_pointing[i] < 0.1:
            t_pointing[i] = -1.0
          elif inv.telescope_pointing[i] > 0.9:
            t_pointing[i] = 1.0

    r_inv.extinction = type(r_inv.Nm)(np.zeros(r_inv.Nm.shape))
    r_inv.extinction_aerosol = type(r_inv.Nm)(np.zeros(r_inv.Nm.shape))
    r_inv.p180 = type(r_inv.Nm)(np.zeros(r_inv.Nm.shape))

    #if (data_len -min_alt/dz) <  z_window_pts*5.0:
    if data_len <  z_window_pts*5.0:
        print 'WARNING---filtered_extinction--filter window length too long'
        print '          reseting z_window_pts from ',z_window_pts,
        #z_window_pts = int((data_len - min_alt/dz)/5.0)
        z_window_pts = int(data_len/5.0)
        z_window_pts = 2*(z_window_pts/2)+1
        print 'to ', z_window_pts
        if z_window_pts < order +2:
            print
            #raise ValueError, 'number of altitude resolution elements two few for filter length'
            print 'WARNING-----Savitzky_golay---number of altitude resolution elements two few for filter length'
            print '            inv.extiniction returned as NaN'
            print
            r_inv.extinction[:,:] = np.NaN * r_inv.Nm
            return r_inv
        print
    #start_pt = int(np.ceil(z_window_pts/2 + min_alt/dz))
    start_pt = int(np.ceil(z_window_pts/2))
    end_pt   = data_len -z_window_pts/2 -1
  
    if adaptive == 0:
        extinction = np.zeros_like(filtered_molecular_counts)
        slope_Nm = np.zeros_like(filtered_molecular_counts)
        dbeta_dr = np.zeros_like(filtered_molecular_counts)
        dbeta_dr[1:] =  + 0.5*(1/r_inv.beta_r_backscat[1:])\
                    *(r_inv.beta_r_backscat[1:]
                    -r_inv.beta_r_backscat[:-1])/adz[1:]
        _m=self.svm(z_window_pts,order,1)
        
        for i in [useindex]:
            #compute the first derivative of a filtered Nm
            tmp,m = sg.savitzky_golay(
                  filtered_molecular_counts[start_pt:end_pt],z_window_pts,order,deriv = 1,m=_m,retval=True)

            #temp_bs,m = sg.savitzky_golay(
            #      filtered_beta_a_backscat[i,start_pt:end_pt],z_window_pts,order,m=_m,retval=True)
            slope_Nm[start_pt:end_pt]=-tmp/adz[start_pt:end_pt]
            #if hasattr(r_inv'telescope_pointing') and r_inv.telescope_pointing <= 0.5:
            #    slope_Nm = -slope_Nm
                    
            if _m is None:
              self.setsvm(m,z_window_pts,order,1)
              _m=m
            extinction[start_pt:end_pt]= \
               (-0.5*(1/filtered_molecular_counts[start_pt:end_pt])
               *slope_Nm[start_pt:end_pt] 
               +dbeta_dr[start_pt:end_pt])*t_pointing[i]     
            if 0:
                    temp = np.zeros_like(r_inv.beta_r_backscat)
                    temp[1:]= r_inv.beta_r_backscat[1:]-r_inv.beta_r_backscat[:-1]
                    import matplotlib.pylab as plt
                    plt.figure(2222)
                    plt.plot(temp,r_inv.msl_altitudes,'r')
                    plt.grid(True)
                    plt.show(False)
                    
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
        yy[0:ntimes,0:data_len]=Nm[0:ntimes,0:data_len]           
        yy[0:ntimes,data_len:data_len+max_half_win] = Nm[0:ntimes,end_range]
        #filtered_Nm = Nm.copy()
        #extinction = np.zeros_like(inv.Nm)
        low_limit = np.zeros(data_len)
        high_limit = np.zeros_like(low_limit)
        beta_a=np.zeros(data_len+max_half_win+1)
    
        for k in [useindex]:#range(ntimes):
            beta_a[0:data_len]=(inv.beta_a_backscat_par[k,0:data_len]\
                      +inv.beta_a_backscat_perp[k,0:data_len])
            beta_a[range((data_len-max_half_win),data_len)] \
                  =(inv.beta_a_backscat_par[k,end_range] \
                        +inv.beta_a_backscat_perp[k,end_range])
    
            #compute integrated backscatter
            beta_a[np.isnan(beta_a)]=0
            int_bs_od=np.cumsum(beta_a)*adz/(2*p180)
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
                
                pc = np.polyfit(x,ylocal,order)
                filtered_Nm[k,i]= np.polyval(pc,i-lo_lmt)
    
                slope_Nm = np.polyval(derivative_coefs*pc[0:order],i-lo_lmt)/adz[i]

                #print 'ext',i,lo_lmt,hi_lmt,ylocal,slope_Nm,pc,pc[0:order]\
                #         ,dz,np.polyval(pc,[0,1,2,3]),t_pointing[k]
       
                extinction = (-0.5*(1/filtered_Nm[k,i])*slope_Nm \
                    + 0.5*(1/r_inv.beta_r_backscat[i])\
                    *(r_inv.beta_r_backscat[i]-r_inv.beta_r_backscat[i-1])/adz[i])\
                    *t_pointing[k]
               
                #print 'ext',i,lo_lmt,hi_lmt,extinction[k,i],slope_Nm\
                #      ,(inv.beta_r_backscat[i]-inv.beta_r_backscat[i-1])/dz\
                #      ,inv.beta_r_backscat[i],inv.beta_r_backscat[i-1],t_pointing[k]
                      
                xx=np.arange(len(x)*10)/10.0
 

    #subtract molecular extinction from total to get aerosol extinction  
    extinction_aerosol = extinction - 8 * np.pi * r_inv.beta_r_backscat/3.0
    r_inv.extinction[0,:] = extinction
    r_inv.extinction_aerosol[0,:] = extinction_aerosol
   
    #compute normalize backscatter phase function for aerosol
    #from integrated_backscatter/optical depth using blocks sized the same as used for
    #extinction filtering
    int_backscat = filtered_beta_a_backscat.copy()
    int_backscat[np.isnan(filtered_beta_a_backscat)] = 0.0
    int_backscat = np.cumsum(int_backscat* adz)
    tau = np.zeros_like(r_inv.extinction)
    tau[0,start_pt:] = -0.5 * np.log(
           filtered_molecular_counts[start_pt:]/filtered_molecular_counts[:-start_pt]\
          *(r_inv.beta_r_backscat[:-start_pt]/r_inv.beta_r_backscat[start_pt:]))

    #subtract Rayliegh optical depth
    tau[0,start_pt:] -= 8.0 * np.pi * (r_inv.beta_r_backscat[start_pt:] \
                    + r_inv.beta_r_backscat[:-start_pt]) * msl_altitudes[start_pt]/6.0
    r_inv.p180 = np.zeros_like(r_inv.beta_a_backscat)
    r_inv.p180[0,start_pt:] = (int_backscat[start_pt:] - int_backscat[:-start_pt])\
                /(tau[0,start_pt:])
    
    """ 
    temp_ext = extinction_aerosol.copy()
    temp_ext[temp_ext == 0] = np.NaN
    r_inv.p180 = filtered_beta_a_backscat/temp_ext
    """
    return r_inv       


def extinction_filter_setup(processing_defaults,rs_constants,provides,mean_provides): 
        ground_based=not rs_constants.has_key('installation') \
                        or  (rs_constants['installation'] == 'ground'
                             or rs_constants['installation'] == 'shipborne')
        if (ground_based and not processing_defaults.enabled('extinction_processing')) or \
           (not ground_based and not processing_defaults.enabled('airborne_extinction_processing')):
          return None
        od_threshhold = processing_defaults.get_value('extinction_processing','od_threshhold')
        z_window_width = processing_defaults.get_value('extinction_processing','alt_window_length')
        t_window_width = processing_defaults.get_value('extinction_processing','time_window_length')
        min_filter_alt = processing_defaults.get_value('extinction_processing','min_alt')
        order = processing_defaults.get_value('extinction_processing','polynomial_order')
        adaptive =  processing_defaults.get_value('extinction_processing','adaptive')
        wfov_corr_enabled = processing_defaults.enabled('wfov_corr') and 'raw_molecular_wfov_counts' in mean_provides
        if ground_based:
            min_filter_alt = np.max([min_filter_alt
                ,rs_constants['lidar_altitude']+100])
            print 'min_filter_alt adjusted for ground based lider elevation'   
        else:
          od_threshhold = processing_defaults.get_value('airborne_extinction_processing','od_threshhold',missing=od_threshhold)
          z_window_width = processing_defaults.get_value('airborne_extinction_processing','alt_window_length',missing=z_window_width)
          t_window_width = processing_defaults.get_value('airborne_extinction_processing','time_window_length',missing=t_window_width)
          min_filter_alt = processing_defaults.get_value('airborne_extinction_processing','min_alt',missing=min_filter_alt)
          order = processing_defaults.get_value('airborne_extinction_processing','polynomial_order',missing=order)
          adaptive =  processing_defaults.get_value('airborne_extinction_processing','adaptive',missing=adaptive)

        print ' '
        if adaptive:
             print '    extinction computed with adaptive Savitzky-Golay filter'
             print '    od threshhold       = ',od_threshhold
        else:
             print '    extinction computed with fixed window Savitzky_Golay filter'
       
        print '    z_window width      = ',z_window_width, ' (m)'
        print '    t_window width      = ',t_window_width, ' (sec)'
        print '    polynomial order    = ',order
        print '    min altitude applied = ',min_filter_alt, ' (m)'

        kwcargs=dict(min_alt=min_filter_alt,window_od=od_threshhold
                        ,max_z_window_length=z_window_width,order=order,adaptive=adaptive
                        ,wfov_corr_enabled=wfov_corr_enabled)

        minwidth=order+1
        if (minwidth%2)!=1:
          minwidth+=1

        return WindowedFilterDescription(
            filter=filtered_extinction(),width=minwidth,modrequirement=(2,1),#ensures odd count
            time_width=timedelta(seconds=t_window_width),edgemode='fullduplicate'
           ,kwcargs=kwcargs)


class wfov_geo_extinction_filter_preprocess(object):
  def __init__(self):
    self.priorTime=None



  def __call__(self,_mean,t_filter_order,correct_below_range,min_fit_alt,z_norm_interval,enable_z_fit):
    """
    def func(x,a):
        #fourth order polynomial with last value = 1 and df/dx of last value = 0
        #return a[0] + a[1] * x + a[2] * x**2 + a[3] * x**3 + a[4] * x**4
        #val = a[0] + a[1] * x**2 + a[2] * x**3 + a[3] * x**4
        #val = val -(2*a[1]*x[-1] + 3*a[2]*x[-1]**2 + 4*a[3]*x[-1]**3) * x
        val = a[0] * x**2 + a[1] * x**3 + a[2] * x**4
        # linear term from df/dx = 0 at x[-1]
        b = -(2*a[0]*x[-1] + 3*a[1]*x[-1]**2 + 4*a[2]*x[-1]**3)
        # constant term from x[-1] = 1
        a = 1.0 - (b*x[-1] + a[0]*x[-1]**2 + a[1]*x[-1]**3 + a[2]*x[-1]**4)
        #full quadratic
        val = a + b*x + val
        return val


    def func2(x,a):
        #fourth order polynomial with last value = 1 and df/dx of first and last value = 0
        #return a[0] + a[1] * x + a[2] * x**2 + a[3] * x**3 + a[4] * x**4
        val = a[0] * x**3 + a[1] * x**4
        #second order term from df/dx = 0 at first and last points
        c = -(3 * a[0] * (x[-1]**2 - x[0]**2) + 4 * a[1] * (x[-1]**3 - x[0]**3))\
                        /(2 * (x[-1] - x[0]))
        # linear term from df/dx = 0 at x[-1]
        b = -(2* c *x[-1] + 3*a[0]*x[-1]**2 + 4*a[1]*x[-1]**3)
        # constant term from x[-1] = 1
        a = 1.0 - (b*x[-1] + c * x[-1]**2 + a[0]*x[-1]**3 + a[1]*x[-1]**4)
        val = a + b*x + c*x**2 + val
        return val
    """
    def func3(x,a):
        #fifth order polynomial with last value = 1 and df/dx of first and last value = 0
        #return a[0] + a[1] * x + a[2] * x**2 + a[3] * x**3 + a[4] * x**4 + a[5] * x**5
        val = a[0] * x**3 + a[1] * x**4 +a[2] * x**5
        #second order term from df/dx = 0 at first and last points
        c = -(3*a[0]*(x[-1]**2 - x[0]**2) + 4*a[1]*(x[-1]**3 - x[0]**3) +5*a[2]*(x[-1]**4-x[0]**4))\
                        /(2*(x[-1] - x[0]))
        # linear term from df/dx = 0 at x[-1]
        b = -(2* c *x[-1] + 3*a[0]*x[-1]**2 + 4*a[1]*x[-1]**3 + 5*a[2]*x[-1]**4)
        # constant term from x[-1] = 1
        a = 1.0 - (b*x[-1] + c * x[-1]**2 + a[0]*x[-1]**3 + a[1]*x[-1]**4 + a[2]*x[-1]**5)
        val = a + b*x + c*x**2 + val
        return val
    
    def residuals(p,x,y):
        res = y - func3(x,p)
        return res

    s_time = datetime.utcnow()
    mean=listMerge(_mean)
    
   
      
    ntimes = mean.times.size
    useindex=0
    if self.priorTime is None:
      self.priorTime=mean.times[0]-timedelta(seconds=10)
    while useindex<(ntimes-1) and self.priorTime>=mean.times[useindex]:
      useindex=useindex+1
    self.priorTime=mean.times[useindex]
   
    r_mean=copy.deepcopy(_mean[useindex])

    print 'wfov window ',r_mean.times[0],r_mean.times[-1],ntimes 
    if not hasattr(mean,'raw_molecular_counts') or not hasattr(mean,'raw_molecular_wfov_counts') or not hasattr(mean,'msl_altitudes'):
      print  
      print 'WARNING: Missing value for wfov calculation. skipped'
      print "hasattr(r_inv,'raw_molecular_wfov_counts') = ", hasattr(mean,'raw_molecular_wfov_counts')
      print "hasattr(r_inv,'msl_altitudes') = ",hasattr(mean,'msl_altitudes')
      return r_mean

    if not hasattr(r_mean,'raw_molecular_wfov_counts') or not hasattr(r_mean,'msl_altitudes'):
      return r_mean

    if not hasattr(r_mean,'calibration_wfov_mol_ratio'):
       print 'CALIBRATION WFOV mol ratio missing!'
       return r_mean

    if r_mean.msl_altitudes[-1] < (correct_below_range + z_norm_interval + r_mean.lidar_altitude) :
          print 'WARNING--max altitude two low for wfov correction--no corr applied'
          print '         wfov cor requires max alt >', (correct_below_range + z_norm_interval + r_mean.lidar_altitude)
          return r_mean

    
    data_len = len(mean.raw_molecular_counts[0,:])    
    msl_altitudes = mean.msl_altitudes

    #distance between altitude bins
    adz=hau.Z_Array(np.zeros(msl_altitudes.shape))
    adz[1:]=msl_altitudes[1:]-msl_altitudes[:-1]
    adz[0]=adz[1]
    
    t_window_pts=ntimes
    ones_array = np.ones_like(mean.raw_molecular_counts)  
    molecular_wfov_counts = (mean.raw_molecular_wfov_counts - mean.m_wfov_dark_counts *ones_array)
    molecular_counts = (mean.raw_molecular_counts - mean.mol_dark_counts * ones_array)
    
    if 0:
  
        import matplotlib.pylab as plt
        plt.figure(1111)
        plt.plot(msl_altitudes,np.nanmean(molecular_counts,0),'c'
           ,msl_altitudes,46*(np.nanmean(molecular_wfov_counts,0)),'r')
        ax=plt.gca()
        ax.set_yscale('log')
    #filter in time
    t_pts=np.arange(len(molecular_counts[:,0]))
    fit_time =  datetime.utcnow()

    top_norm_bin = len(msl_altitudes[msl_altitudes <= correct_below_range + z_norm_interval])
    top_bin = len(msl_altitudes[msl_altitudes <= correct_below_range])
    bot_bin = len(msl_altitudes[msl_altitudes <= min_fit_alt])
    
    if t_filter_order <> 1: 
        pc = np.polyfit(t_pts,molecular_counts,t_filter_order)           
        filtered_molecular_counts = np.polyval(pc,useindex)
        pcw = np.polyfit(t_pts,molecular_wfov_counts,t_filter_order)
        filtered_molecular_wfov_counts = np.polyval(pcw,useindex)
        pcw = np.polyfit(t_pts,mean.raw_molecular_wfov_counts,t_filter_order)
        filtered_raw_molecular_wfov_counts = np.polyval(pcw,useindex)
    else: # simplified computation when fit is linear
        filtered_molecular_counts = np.nanmean(molecular_counts[:,:top_norm_bin],0) 
        filtered_molecular_wfov_counts = np.nanmean(molecular_wfov_counts[:,:top_norm_bin],0)
        filtered_raw_molecular_wfov_counts = np.nanmean(mean.raw_molecular_wfov_counts[:,:top_norm_bin],0)    
   
    
    wfov_mol_ratio = filtered_molecular_wfov_counts /filtered_molecular_counts

    #wfov_mol_ratio is now ratioed to the value measured with the geometric correction 
    wfov_mol_ratio /= mean.calibration_wfov_mol_ratio[:top_norm_bin]
    if 0:
        import matplotlib.pylab as plt
        plt.figure(6666)
        plt.plot(np.arange(top_norm_bin),wfov_mol_ratio,'b'\
                 ,np.arange(top_norm_bin),wfov_mol_ratio[:top_norm_bin],'r')
        
                 
    # add functional code here.
    #   'inv' is a normal 2-D frame of data (read only)
    #   'r_inv' is the 1-D frame to be returned.
    #   'useindex' is the time index of the 'r_inv' frame within 'inv'

   

    
    #normalize wfov_correction by mean value between top_bin and 1.5 km above top_bin


    wfov_mol_ratio /= np.nanmean(wfov_mol_ratio[top_bin:top_norm_bin])
    #wfov_mol_ratio /= wfov_mol_ratio[top_bin]

    if 0:
        import matplotlib.pylab as plt
        plt.figure(6667)
        plt.plot(msl_altitudes[:top_norm_bin],wfov_mol_ratio,'g'\
                 ,msl_altitudes[:top_norm_bin],wfov_mol_ratio[:top_norm_bin],'r')
        plt.grid(True)
        
        
    #bin_vec = np.arange(len(wfov_mol_ratio))
    #smoothed_wfov_mol_ratio = np.ones_like(wfov_mol_ratio)
    bin_vec = np.arange(top_norm_bin)
    smoothed_wfov_mol_ratio = np.ones(len(msl_altitudes))
    if enable_z_fit == "spline":
        #wfov fit with spline
        #forces correction factor to 1.0 and slope to zero at top of layer
        #weights = np.sqrt(filtered_molecular_wfov_counts[bot_bin:top_bin])
        weights = filtered_molecular_wfov_counts[bot_bin:top_bin] \
                  / np.sqrt(filtered_raw_molecular_wfov_counts[bot_bin:top_bin])
        temp =wfov_mol_ratio[bot_bin:top_bin].copy()
        weights = weights[np.isnan(temp)==0]
        bins = bin_vec[bot_bin:top_bin].copy()
        bins = bins[np.isnan(temp)==0]
        temp = temp[np.isnan(temp)==0]
        
        #force correction factor to 1.0 and slope = 0.0 at top of corrected layer
        temp[-2:]=1.0
        weights[-2:]=1000.0
        
        npts = bins.shape[0] - 1-2 
        s = 0.05*(npts-np.sqrt(2*npts))
        tck = splrep(bins,temp,w=weights, s=s) 
        smoothed_wfov_mol_ratio[bot_bin:top_bin] =splev(bin_vec[bot_bin:top_bin],tck)
    elif enable_z_fit == "polynomial":      #smooth over ranges bot_bin to top_bin if enabled
        # 'wfov correction using polynomial fit'
        temp = wfov_mol_ratio[bot_bin:top_bin].copy()
        temp[np.isnan(temp)] = 0
        p = leastsq(residuals
                ,x0=np.array([-1.0e-5,0.0000,0.0])
                ,args=(bin_vec[bot_bin:top_bin] \
                ,temp))
           
        #use fitted values between min_fit_alt and fit_below_alt
        smoothed_wfov_mol_ratio[bot_bin:top_bin] = func3(bin_vec[bot_bin:top_bin]
                           ,p[0])
       
    else: #no smoothing
        # 'wfov correction with no range smoothing'
        smoothed_wfov_mol_ratio[:top_bin]=wfov_mol_ratio[:top_bin]
        
    #make value constant below bot_bin
    smoothed_wfov_mol_ratio[:bot_bin] = smoothed_wfov_mol_ratio[bot_bin]

    #smallestval=nanmin(smoothed_wfov_mol_ratio[smoothed_wfov_mol_ratio>0.0])
    #smoothed_wfov_mol_ratio[np.logical_not(smoothed_wfov_mol_ratio>0.0)]=smallestval
    
    if 0:  #plot fit as function of altitude
       import matplotlib.pylab as plt
       plt.figure(6668)
       plt.plot(msl_altitudes[:top_bin],wfov_mol_ratio[:top_bin],'c'
            ,msl_altitudes[bot_bin:top_bin],smoothed_wfov_mol_ratio[bot_bin:top_bin],'r'
            ,msl_altitudes[top_bin:],smoothed_wfov_mol_ratio[top_bin:],'g')
       plt.grid(True)
       plt.xlabel('Range (m)')
       plt.ylabel('Correction')
       ax=plt.gca()
       ax.set_ylim(0.7,1.1)
       ax.set_xlim(0,8000)
       
       #plt.show()
    r_mean.filtered_wfov_mol_ratio = hau.TZ_Array(np.ones(r_mean.raw_molecular_counts.shape))

    #no correction for NaN points
    smoothed_wfov_mol_ratio[np.isnan(smoothed_wfov_mol_ratio)] = 1.0
    #limit size of correction
    smoothed_wfov_mol_ratio[smoothed_wfov_mol_ratio < 0.5] = 0.5
    
    r_mean.filtered_wfov_mol_ratio[0,:] = smoothed_wfov_mol_ratio

    
    r_mean.wfov_extinction_corr = hau.TZ_Array(np.zeros(r_mean.raw_molecular_counts.shape))
    r_mean.wfov_extinction_corr[0,1:top_bin-1] = smoothed_wfov_mol_ratio[2:top_bin]-  smoothed_wfov_mol_ratio[:top_bin-2]
    r_mean.wfov_extinction_corr[0,1:top_bin-1] = -0.25 * r_mean.wfov_extinction_corr[0,1:top_bin-1] \
                                           /(smoothed_wfov_mol_ratio[1:top_bin-1] * adz[1:top_bin-1])

    # make wfov geo correction to mean counts.
    for chan in vars(r_mean).keys():
        if not chan.endswith('_counts'):
          continue
        applyFilter=True
        for badPart in ['dark','raw_','var_','wfov']:
          if badPart in chan:
            applyFilter=False
            break
        if not applyFilter:
          continue
        if hasattr(r_mean,chan):
            setattr(r_mean,chan,(getattr(r_mean,chan) *   r_mean.filtered_wfov_mol_ratio))

    return r_mean       


def wfov_geo_extinction_filter_setup(processing_defaults,rs_constants,provides,cal_provides):

        if not processing_defaults.enabled('wfov_corr') or 'raw_molecular_wfov_counts' not in provides:
            return None
      
        if 'rs_cal' not in cal_provides or 'geo' not in cal_provides['rs_cal'] or 'wfov_mol_ratio' not in cal_provides['rs_cal']['geo']:
            print 'NO WFOV MOL RATIO in calibration geo table'
            return None
        # figure out any parameters to the __call__ function above. add them to kwcargs dictionary
        #t_window_width = 0 #target number of seconds for the window FIXME
        #minwidth = 0 #minimum number of time steps in a window

        minwidth = 0
        
        t_window_width = processing_defaults.get_value('wfov_corr','window_durration')
        t_filter_order = processing_defaults.get_value('wfov_corr','time_filter_order')
        min_fit_alt = processing_defaults.get_value('wfov_corr','min_fit_range')\
                      * 1000.0 + rs_constants['lidar_altitude']
        correct_below_range = processing_defaults.get_value('wfov_corr','correct_below_range')\
                      * 1000.0 +rs_constants['lidar_altitude']
        try:
            z_norm_interval = processing_defaults.get_value('wfov_corr','z_norm_interval') \
                          *1000.0
        except:
            z_norm_interval = 1500.0
            print
            print "z_norm_interval not provided in 'wfov_corr' processing_defaults, default=1500 m"
            print
        try:
            enable_z_fit = processing_defaults.get_value('wfov_corr','enable_z_fit')
            print "wfov correction using " + enable_z_fit
        except:
            enable_z_fit = 1
            print
            print "z_fit_enable not provided in 'wfov_corr' processing_defaults, default=True"
            print         
        #kwcargs=dict(additionalparameterFIXME=None)
        kwcargs=dict(t_filter_order=t_filter_order
                     ,min_fit_alt=min_fit_alt,correct_below_range=correct_below_range,z_norm_interval=z_norm_interval,enable_z_fit=enable_z_fit)

        minwidth=t_filter_order+1
        if (minwidth%2)!=1:
          minwidth+=1

        return WindowedFilterDescription(
            filter=wfov_geo_extinction_filter_preprocess(),width=minwidth,modrequirement=(2,1),#ensures odd count
            time_width=timedelta(seconds=t_window_width),edgemode='fullduplicate'
           ,kwcargs=kwcargs)


def mean_filter_setup(processing_defaults,rs_constants,provides,cal_provides):
      
      ret=[]
      print
      print
      print 'entering mean_filter_setup'
      print
      wfovconfig=wfov_geo_extinction_filter_setup(processing_defaults,rs_constants,provides,cal_provides)
      if wfovconfig is not None:
        ret.append(wfovconfig)

      #if the wfov geo extinction filter should trigger a different set of code from the 'exctinction filter setup',
      # the below can be skipped in favor of another extinction filter that is deliberately intended for wfov data cases.

      if len(ret)==0:
        return None
      return ret


def inv_filter_setup(processing_defaults,rs_constants,provides,mean_provides):
      
      ret=[]

      #if the wfov geo extinction filter should trigger a different set of code from the 'exctinction filter setup',
      # the below can be skipped in favor of another extinction filter that is deliberately intended for wfov data cases.

      tmp=extinction_filter_setup(processing_defaults,rs_constants,provides,mean_provides)
      if tmp is not None:
        ret.append(tmp)

      if len(ret)==0:
        return None
      return ret

