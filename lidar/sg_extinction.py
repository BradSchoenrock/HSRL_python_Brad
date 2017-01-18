import numpy as np
import hsrl.filters.savitzky_golay as sg


def polyFilter(Nm,order):
    """ Input of TxA matrix
    returns a 1xA vector filtered using polynomial order given
    """
    t_pts = np.arange(Nm.shape[0])
    t_center_pt=Nm.shape[0]/2
    filtered_Nm = Nm[np.arange(1),:].copy()
    for x in range(filtered_Nm.shape[1]):
        p = np.polyfit(t_pts,Nm[:,x],order)
        #get z-array of fitted values at center pt in time dimension
        filtered_Nm[0,x] = np.polyval(p,t_center_pt)
    return filtered_Nm

class sg_extinction(object):
    def __init__(self,Filter):
        self.order = Filter['order']
        self.z_window_pts = Filter['z_window_pts']
        self.start_bin = Filter['start_bin']
        self.dz = Filter['dz']
        self.adaptive = False

    def __call__(self,times,delta_ts,Nm,beta_a_backscat,integ_backscat,telescope_pointing=None,beta_r=None,state=[]):
            """filtered_extinction(Nm,beta_a_backscat,beta_r_backscat,integ_backscat
               ,msl_altitudes,process_control,telescope_pointing=None)
               
          
               Nm(0:nt,0:nz)= corrected molecular counts--slice of nt time points
                            for first and last t_window_pts/2 time steps Nm  must
                            be padded to provide full width slices.
               nz       = number of points in full altitude profile
               beta_a_backscat = altitude vector of beta_a_backscat at center of slice
               int_backscat    = altitude vector of integrated backscat at center of slice
            """
            if 0: #Nm.shape[0]==1:
                bin_vec = np.arange(Nm.shape[1])
                import matplotlib.pylab as plt
                plt.figure(2000)
                plt.plot(bin_vec,Nm[0,:])
                
            if len(state)==0:
                dbeta_dr = np.zeros_like(beta_r)
                dbeta_dr[1:] =  +  0.5*(1/beta_r[1:])\
                    * (beta_r[1:] - beta_r[:-1])/self.dz
                state.append(dbeta_dr)
            dbeta_dr=state[0]
    
            if Nm.shape[0] >1:  #if Nm has time dimension greater that 1
                #fit in time dimension
                filtered_Nm = runTimeFilter(Nm,self.order)
            else:
                filtered_Nm = Nm.copy()
                
            if self.adaptive == False:
                slope_Nm = np.NaN * filtered_Nm
                extinction = np.NaN * filtered_Nm
            
                #compute the first derivative of a filtered Nm 
                slope_Nm[0,self.start_bin:] = -sg.savitzky_golay(
                      filtered_Nm[0,self.start_bin:]
                      ,self.z_window_pts,self.order,deriv = 1)/self.dz
                
                #compute total extinction (aerosol+molecular)
                if telescope_pointing == None or telescope_pointing > .9 : #zenith pointing tel
                    extinction[0,self.start_bin:]= \
                       (-0.5*(1/filtered_Nm[0,self.start_bin:])
                        *slope_Nm[0,self.start_bin:] 
                       +dbeta_dr[self.start_bin:])
                elif telescope_pointing<.1:  #downward pointing telescope
                    extinction[0,start_bin:end_bin]= \
                       -(-0.5*(1/filtered_Nm[0,self.start_bin:])*slope_Nm[0,self.start_bin:] 
                       +dbeta_dr[self.start_bin:])
    
                #compute aerosol extinction alone                
                extinction_aerosol = extinction - beta_r
               
                temp_ext = extinction_aerosol.copy()
                temp_ext[temp_ext ==0] = np.NaN
                p180 = beta_a_backscat/temp_ext
               
                return extinction, extinction_aerosol, p180
            else:
               raise NotImplemented('adaptive extinction filter not present') 
        

def filter_setup(altitudes,processing_defaults,rs_constants,delta_t=None):
        """filter_setup(altitudes,beta_r_backscat,delta_t,processing_defaults,rs_constants
           sets up filter parameters
           altitudes = altitude vector
           delta_t   = time spacing between profiles
        """
        
        dz = altitudes[1] - altitudes[0]
        order = processing_defaults.get_value('extinction_processing'
                                                    ,'polynomial_order')
        od_threshhold = processing_defaults.get_value('extinction_processing'
                                                      ,'od_threshhold')

        #calculate number of points in altitude window
        z_window_pts = int(processing_defaults.get_value('extinction_processing'
                ,'alt_window_length')/dz)
        z_window_pts = max(z_window_pts , order +1)
        if (z_window_pts%2)!=1:
            z_window_pts = z_window_pts +1
        z_window_width = z_window_pts * dz
        t_window_time=processing_defaults.get_value('extinction_processing'
                 ,'time_window_length')
        if delta_t is not None and delta_t>0:
            #calculate number of points in time window
            t_window_pts = int(t_window_time/delta_t)
            t_window_pts = max(t_window_pts , order +1)
            if (t_window_pts%2)!=1:
                t_window_pts = t_window_pts +1
            t_window_width = t_window_pts * delta_t
        else:
            t_window_pts=None
            t_window_width=None

        #find lowest bin to filter
        min_filter_alt = processing_defaults.get_value('extinction_processing','min_alt')
        first_bin_to_process = processing_defaults.get_value('first_bin_to_process','bin_number')
        binwidth = rs_constants['binwidth'] * 1.5e8  #width(m) = width(sec) * speed_of_light/2
        if not rs_constants.has_key('installation') \
                        or  (rs_constants['installation'] == 'ground'
                             or rs_constants['installation'] == 'shipborne'):
            min_filter_alt = np.max([min_filter_alt
                ,rs_constants['lidar_altitude']+binwidth * first_bin_to_process]) 
            print 'min_filter_alt adjusted for ground based lidar elevation'
        start_bin = int(min_filter_alt / dz) + 1    
        min_filter_alt = start_bin * dz 
        adaptive =  processing_defaults.get_value('extinction_processing','adaptive')
       
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

        WindowedFilterDescription = dict(order=order,t_window_pts=t_window_pts,dz=dz
             ,z_window_pts=z_window_pts,start_bin=start_bin,min_filter_alt=min_filter_alt
             ,z_window_width=z_window_width,t_window_width=t_window_width
             ,adaptive=adaptive,max_z_window_length=z_window_width
             ,window_od=od_threshhold,t_window_time=t_window_time)

       

        return WindowedFilterDescription

    


