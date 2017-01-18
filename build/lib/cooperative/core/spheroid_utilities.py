import numpy as np
import math as math
import matplotlib.pylab as plt
import scipy.special as gm
import os
import lg_base.core.read_utilities as ru
from lg_base.core.locate_file import locate_file
from datetime import datetime,timedelta


def K_squared(wavelength,temperature = None):
    if wavelength < 5.0e-3: #3 mm radar
             k_sq_water = 0.7
    elif wavelength < 1e-3: # 8.6 mm radar
             k_sq_water = 0.88
    else:    #x-band
             k_sq_water = 0.93
    k_sq_ice = 0.176
	
    return k_sq_water,k_sq_ice



def read_mie_files(radar_wavelength):
        """[D,qsca_lidar,qback_lidar,qsca_radar,qback_radar]
                 =read_mie_files(radar_wavelength)
           read scattering and backscatter efficiencies from mie files
	   and provide values scattering efficeinies as funtion of diameter.
           radar_wavelength in meters (either 3.154e-3 or 8.6e-3)
           D = vector of diameters (microns) at which the efficiencies are provided.
        """
	scattering_model = 'Mie'   # 'Rayleigh' OR  'Mie'
	
        if scattering_model == 'Mie' :  #for mie cross sections
	    #reading mie theory for lidar	
            mie_filename = locate_file('mie_coefficients_0_18000.txt')
            if os.path.isfile(mie_filename):
                print 'reading  ' + mie_filename
                [hh_lidar,dd_lidar]=ru.readascii(mie_filename)
                lambda_lidar = 0.532  #in microns
            else:
                print "unable to compute radar/lidar backscatter ratio"
                raise RuntimeError, ("Can't find mie results in '",mie_filename,"'")
            print
            print 'lambda', radar_wavelength
            print
            if radar_wavelength == 3.154e-3:
                radar_mie_filename = locate_file('mie_data_3.2mm_285K.txt')
                lambda_radar = radar_wavelength * 1e6 #in microns
            elif radar_wavelength == 8.6e-3:
                lambda_radar = radar_wavelength * 1e6 #wavelength in m converted to microns
                radar_mie_filename = locate_file('mie_data_8.6mm_285K.txt')
                
            else:
                print "unable to compute radar/lidar backscatter ratio"
                raise RuntimeError, ("Can't find mie results in '",radar_mie_filename,"'")
	if scattering_model == 'Mie':
            print 'reading  '+ radar_mie_filename
            
            [hh_radar,dd_radar]=ru.readascii(radar_mie_filename)

            #convert wavenumber to diameter in microns
            dmax_lidar = dd_lidar[-1,0]*lambda_lidar/np.pi
            #find largest diameter in microns for this radar wavelength
            dmax_radar = dd_radar[-1,0] * lambda_radar /np.pi

            #limit max size to smallest max of lidar or radar
            d_max = np.min([dmax_lidar,dmax_radar])
            D = np.arange(d_max,dtype='float')
            #place lidar and radar backscat eff on common scale--1 micron/pt
	    #and divide backscatter by 4pi to conform to Deirmenjian normalized phase function
            qsca_lidar = np.interp(D,dd_lidar[:,0]*lambda_lidar/np.pi,dd_lidar[:,2])
            qback_lidar = np.interp(D,dd_lidar[:,0]*lambda_lidar/np.pi,dd_lidar[:,3])
	    qback_lidar = qback_lidar/(4.0 *np.pi)
            qsca_radar = np.interp(D,dd_radar[:,0]*lambda_radar/np.pi,dd_radar[:,2])
            qback_radar = np.interp(D,dd_radar[:,0]*lambda_radar/np.pi,dd_radar[:,3])
            qback_radar = qback_radar/(4.0*np.pi)
	    
            if 0:
                 import matplotlib.pylab as plt
                
                 
                 plt.plot(D,size_dist.area_weighted_number(D,D[200],'water'))
                 plt.ylabel('number density')
                 plt.xlabel('Diameter (microns)')
                 ax=plt.gca()
                 ax.set_xscale('log')
                 ax.grid(True)
                 
	#for rayleigh radar cross section and geometeric lidar with p180/4pi=0.65	 
	if scattering_model == 'Rayleigh':  
	    D = np.arange(10000,dtype='float')
	    k_sq_water,k_sq_ice = K_squared(radar_wavelength)
	    qsca_lidar = 2.0 * np.ones_like(D)
	    qback_lidar = 0.065 * qsca_lidar
	    #qsca_radar from van de Hulst page 70
	    qsca_radar = (8.0/3.0) * np.pi**4 * k_sq_water * (D*1e-6 / radar_wavelength)**4
	    qback_radar = 3 * qsca_radar/(8.0 * np.pi)

	return D,qsca_lidar,qback_lidar,qsca_radar,qback_radar


 
    
class size_dist(object):
    def __init__(self,particle_parameters,radar_wavelength):
	  """  
          __init__(self,particle_parameters,radar_wavelength)
          Calls to size_dist methods express particle sizes in meters
	  Internally size_dist converts particle sizes to microns.
	  """
          #read the mie files for lidar and radar scattering efficiencies
          [self.D,self.qsca_lidar,self.qback_lidar,self.qsca_radar,self.qback_radar] \
                    = read_mie_files(radar_wavelength)

          #internally--size_dist operates with sizes and radar wavelength in microns
	  #external calls to methods supply sizes in meters

	  #convert radar wavelength from meters to microns
	  self.radar_wavelength = radar_wavelength * 1e6
          
          #get size distribution discription from particle parameters
          self.distribution = particle_parameters['size_distribution']['form']
          self.alpha_water = particle_parameters['size_distribution']['alpha_water']
          self.alpha_ice = particle_parameters['size_distribution']['alpha_water']
          self.g_water = particle_parameters['size_distribution']['g_water']
          self.g_ice = particle_parameters['size_distribution']['g_water']

         
    def area_weighted_number(self,D,Dm,phase_name):
         if self.distribution == 'modified_gamma':
             """modified_gamma(x,Dm,alpha,gam):
                D = vector of sizes 
                Dm = mode diameter
                N ~ D^alpha *exp(-(alpha/gam)*(D/Dm)^gam)
                N = Relative number of particles as function of D
                    for a given Dm,alpha,gam
		output is independent of units as long as D and Dm are consistant.    
                """
             
             if phase_name == 'water':  #for water
                 alpha = self.alpha_water
                 gam = self.g_water
             elif phase_name == 'ice':
                 alpha  = self.alpha_ice
                 gam   = self.g_ice
             else:
                 raise RuntimeError('unrecognized phase requested') 
             #area weighted size distribution
             N = D**(alpha + 2) * np.exp(-(alpha/gam)*(D/Dm)**gam)
             N = N/np.nansum(N)
             return N
         elif self.distribution == 'inverse_modified_gamma':
             """inverse_modified_gamma(x,Dm,alpha,gam):
                D = vector of sizes 
                Dm = mode diameter
                N ~ D^alpha *exp(-(alpha/gam)*(D/Dm)^-gam)
                N = Relative number of particles as function of D
                    for a given Dm,alpha,gam
                """
         
             if phase_name == 'water':  
                 alpha = self.alpha_water
                 gam = self.g_water
             else:
                 alpha  = self.alpha_ice
                 gam   = self.g_ice
             #area weighted size distribution
             N = D**(2 - alpha) * np.exp(-(alpha/gam)*(D/Dm)**-gam)
             N = N/np.nansum(N)
             return N
             
         elif self.distribution == 'oconnor_gamma':
                  
             """oconnor_gamma(D,D0,mu)
                O'Connors normalized gamma distribution
                D = Volume weighted Diameter (microns)
                D0 = median of volume weighted diameter
                alpha = shape factor = mu in O'Connor

                Notice (alpha + 2) to provide particle area weighted distribution
                """
             if phase_name == 0:  #for water
                mu = self.alpha_water
             else:
                mu  = self.alpha_ice

             #convert mode diameter to median volume diameter
             D0 = Dm *(3.67+mu)/mu
             N = (D/D0)**(mu + 2) * np.exp(-(3.67+mu)*D/D0)
             N = N/np.sum(N)
             return N
         
    def number(self,D,Dm,phase_name):
        """number(self,D,Dm,phase_name)
           convert area weighted number density to number density
        """
	N = self.area_weighted_number(D,Dm,phase_name)*D**-2.0
        N = N/np.nansum(N)
        return N
    
    def radar_weighted_sizes(self,D,dD,Dm,phase_name,zeta=None):
        """radar_weighted_sizes(self,D,Dm,phase,zeta=None)
           compute radar weighted number density distribution
           D = vector of diameters (m)
           dD = seperation between D values for integration
           Dm = mode diameter (m)
           phase_name =  'ice' or  'water' 
           zeta = ellipticity of ice spheroid (sphere=1,zero thickness sheet=0)
        """
        dist = np.zeros_like(D)
        
        if phase_name == 'water':
                #radar weighted size distribution water
                #dist = self.number(D,Dm,phase_name) * D**6.0
                index = D.copy() * 1e6
                index[index >= len(self.qback_radar)] = len(self.qback_radar)-1
                index = index.astype('int')
                #dist = self.number(D,Dm,phase_name) * D**6.0
                dist = self.number(D,Dm,phase_name) * D**2 * self.qback_radar[index]       
                dist /= np.nansum(dist * dD)

        elif phase_name == 'ice':
                 #radar weighted size distribution for ice
                 dist = self.number(D,Dm,phase_name) * D**(4.0+2.0*zeta)
                 dist /= np.nansum(dist * dD)
        else:
                raise RuntimeError('unrecognized phase requested')
       
        return dist
    
    def mass_weighted_sizes(self,D,dD,Dm,phase_name,zeta=None):
        """mass_weighted_sizes(self,D,Dm,phase,zeta=None)
           compute mass weighted number density distribution
           D = vector of diameters
           dD = seperation between D values for integration
           Dm = mode diameter
           phase =  'ice' or  'water' 
           zeta = ellipticity of ice spheroid (sphere=1,zero thickness sheet=0)
        """
        dist = np.zeros_like(D)
        if self.distribution == 'modified_gamma':
            if phase_name == 'water':
                #mass weighted size distribution water
                dist = self.number(D,dD,Dm,phase_name) * D**3.0
                dist /= np.nansum(dist * dD)
            elif phase_name == 'ice':
                 #mass weighted size distribution for ice
                 dist = self.number(D,dD,Dm,phase_name,zeta) * D**(2.0+zeta)
                 dist /= np.nansum(dist * dD)
                          
            else:
                raise RuntimeError('unrecognized phase requested')
        
        elif self.distribution == 'inverse_modified_gamma':
             pass
        elif self.distribution == 'oconnor_gamma':
             pass
        else:
             raise RuntimeError('unrecognized size distribution')
        return dist

    def eff_diameter(self,Dm_array,phase_array,zeta_array=None):
        """eff_diameter(self,Dm_array,phase_array,zeta_array=None)
           compute eff_diameter from array of mode diameters
           Dm_array = array of mode diameters
           phase_array = array of phase, (ice = 1, water =0)
           zeta_array = ellipticity of ice spheroids (sphere=1,zero thickness sheet=0)
        """
        eff_diameter = np.NaN * np.zeros_like(Dm_array)
        if self.distribution == 'modified_gamma':
             #do all water points (zeta==0)
             factor = (self.g_water/self.alpha_water)**(1.0/self.g_water) \
                       *gm.gamma((self.alpha_water + 4)/self.g_water)\
                       /gm.gamma((self.alpha_water +3)/self.g_water)     
             eff_diameter[phase_array==0] = factor * Dm_array[phase_array==0]
             if not zeta_array == None:                         
                 eff_diameter[phase_array==1] = Dm_array[phase_array==1]\
                       *(self.g_ice/self.alpha_ice)**(zeta_array[phase_array==1]/self.g_ice) \
                       *gm.gamma((self.alpha_ice + zeta_array[phase_array==1] +3)/self.g_ice)\
                       /gm.gamma((self.alpha_ice +3)/self.g_ice)                          
        
        elif self.distribution == 'inverse_modified_gamma':
            #do all water points (zeta==0)
             factor = (self.alpha_water/self.g_water)**(1.0/self.g_water) \
                       *gm.gamma((self.alpha_water - 4)/self.g_water)\
                       /gm.gamma((self.alpha_water - 3)/self.g_water)     
             eff_diameter[phase_array==0] = factor * Dm_array[phase_array==0]
             if not zeta_array == None:                         
                 eff_diameter[phase_array==1] = Dm_array[phase_array==1]\
                       *(self.alpha_ice/self.g_ice)**(zeta_array[phase_array==1]/self.g_ice) \
                       *gm.gamma((self.alpha_ice - zeta_array[phase_array==1] - 3)/-self.g_ice)\
                       /gm.gamma((self.alpha_ice - 3)/self.g_ice)                          
          
        elif self.distribution == 'oconnor_gamma':
             pass
        else:
             raise RuntimeError('unrecognized size distribution')
        return eff_diameter
   
    def deff_prime(self,Dm_array,phase_array,zeta_array):
        """deff_prime(self,Dm_array,phase_array,zeta_array = None)
        create deff_prime array the same size as Dm
        compute deff_prime from Dm for those points with phase == 0 (water)
        If zeta and Dm are given compute deff_prime for points with phase==1(ice)
        If no zeta array phase==1 points are left as NaN's
        """

        deffprime  = np.NaN * np.zeros_like(Dm_array)
        if self.distribution == 'modified_gamma':
            #do all water points (zeta==0)                       
            factor = (self.g_water/self.alpha_water)**(1/self.g_water)\
                     * (gm.gamma((self.alpha_water + 7)/self.g_water) \
                     /(gm.gamma((self.alpha_water + 3)/self.g_water)))**0.25
            deffprime[phase_array==0] = factor * Dm_array[phase_array==0]
            if not zeta_array == None:  #do all ice points
                deffprime[phase_array==1] = Dm_array[phase_array==1]\
                  *(self.g_ice/self.alpha_ice)**(zeta_array[phase_array==1]/self.g_ice)\
                  *(gm.gamma((self.alpha_ice + zeta_array[phase_array==1] + 6)/self.g_ice) \
                  /(gm.gamma((self.alpha_ice + 3)/self.g_water)))**0.25
            
        elif self.distribution == 'inverse_modified_gamma':
            #do all water points (zeta==0)                       
            factor = (self.alpha_water/self.g_water)**(1/self.g_water)\
                     * (gm.gamma((self.alpha_water - 7)/self.g_water) \
                     /(gm.gamma((self.alpha_water - 3)/self.g_water)))**0.25
            deffprime[phase_array==0] = factor * Dm_array[phase_array==0]
            if not zeta_array == None:  #do all ice points
                deffprime[phase_array==1] = Dm_array[phase_array==1]\
                  *(self.alpha_ice/self.g_ice)**(zeta_array[phase_array==1]/self.g_ice)\
                  *(gm.gamma((self.alpha_ice - 2 * zeta_array[phase_array==1] - 5)/self.g_ice) \
                  /(gm.gamma((self.alpha_ice - 3)/self.g_water)))**0.25
             
        elif self.distribution == 'oconnor_gamma':
             pass
        else:
             raise RuntimeError('unrecognized size distribution')


        return deffprime

    def mean_diameter(self,Dm_array,phase_array,zeta_array=None):
        """mean_diameter(self,Dm_array,phase_array,zeta_array=None)
           compute eff_diameter from array of mode diameters
           Dm_array = array of mode diameters
           phase_array = array of phase, (ice = 1, water =0)
           zeta_array = ellipticity of ice spheroids (sphere=1,zero thickness sheet=0)
        """
        #mean diameter from mode diameter
        mean_diameter = np.NaN * np.zeros_like(Dm_array)
        if self.distribution == 'modified_gamma':
             #do all water points (zeta==0)
             factor = (self.g_water/self.alpha_water)**(1.0/self.g_water) \
                       *gm.gamma((self.alpha_water + 2)/self.g_water)\
                       /gm.gamma((self.alpha_water +1)/self.g_water)     
             mean_diameter[phase_array==0] = factor * Dm_array[phase_array==0]
             if not zeta_array == None:                         
                 mean_diameter[phase_array==1] = Dm_array[phase_array==1]\
                    *(self.g_ice/self.alpha_ice)**(zeta_array[phase_array==1]/self.g_ice)\
                    *gm.gamma((self.alpha_ice + zeta_array[phase_array==1] +1)/self.g_ice)\
                    /gm.gamma((self.alpha_ice +1)/self.g_ice)                          
        
        elif self.distribution == 'inverse_modified_gamma':
            factor = (self.alpha_water/self.g_water)**(1/self.g_water)\
                     * (gm.gamma((self.alpha_water - 2)/self.g_water)\
                     /gm.gamma((self.alpha_water - 1.0)/self.g_water)) 
            mean_diameter[phase_array==0] = factor * Dm_array[phase_array==0]
            if not zeta_array == None:  #do all ice points
                deffprime[phase_array==1] = Dm_array[phase_array==1]\
                  (self.alpha_ice/self.g_ice)**(zeta_array[phase_array==1]/self.g_ice)\
                  *(gm.gamma((self.alpha_ice - zeta_array[phase_array==1] - 1)/self.g_ice))/gm.gamma((self.alpha_ice - 1)/self.g_ice) 
        elif self.distribution == 'oconnor_gamma':
             pass
        else:
             raise RuntimeError('unrecognized size distribution')
        return mean_diameter

    def dmode_from_lidar_radar_rayleigh(self,dmode,beta_ext_lidar,radar_backscatter\
                           ,zeta,phase):
         """dmode_from_lidar_radar(beta_ext_lidar,radar_backscatter \
                           ,zeta,phase)
            compute the mode diameter directly from the lidar extinction
            and radar backscatter cross sections this function assumes geometric
            optics cross section for lidar and Rayleigh scattering for radar"""

    
         k_sq_water,k_sq_ice = K_squared(self.radar_wavelength) 
    
         if self.distribution == 'modified_gamma':
             ag_water = (self.alpha_water/self.g_water)**(1.0/self.g_water)
             ag_ice = (self.alpha_ice/self.g_ice)**(1.0/self.g_ice)
             dmode[phase==0] = ag_water * self.radar_wavelength \
                 * (2.0 * radar_backscatter[phase==0] \
                  / (np.pi**3 * beta_ext_lidar[phase==0])\
                  * gm.gamma((self.alpha_water +3.0)/self.g_water)
                  /(gm.gamma((self.alpha_water +7.0)/self.g_water)*k_sq_water))**0.25
    
             dmode[phase>0] = ag_ice \
                 *(2.0 * radar_backscatter[phase ==1] \
                 /(np.pi**3 * self.radar_wavelength**4 *beta_ext_lidar[phase==1] \
                 * gm.gamma((self.alpha_ice +3.0)/self.g_ice)\
                 / (gm.gamma((2*zeta[phase>0] + self.alpha_ice +5.0)/self.g_ice) * k_sq_ice)) \
                 **(1.0 / (2.0 * zeta[phase>0] +2.0)))
         elif self.distribution == 'inverse_modified_gamma':
              ag_water = (self.alpha_water/self.g_water)**(1.0/self.g_water)
              ag_ice = (self.alpha_ice/self.g_ice)**(1.0/self.g_ice)
              dmode[phase==0] = ag_water * self.radar_wavelength \
                 * (2.0 * radar_backscatter[phase==0] \
                  / (np.pi**3 * beta_ext_lidar[phase==0])\
                  * gm.gamma((self.alpha_water - 3.0)/self.g_water)
                  /(gm.gamma((self.alpha_water - 7.0)/self.g_water)*k_sq_water))**0.25
              dmode[phase>0] = ag_ice \
                 *(2.0 * radar_backscatter[phase ==1] \
                 /(np.pi**3 * self.radar_wavelength**4 *beta_ext_lidar[phase==1] \
                 * gm.gamma((self.alpha_ice - 3.0)/self.g_ice)\
                 / (gm.gamma((-2*zeta[phase>0] + self.alpha_ice -5.0)/self.g_ice) * k_sq_ice)) \
                 **(1.0 / (2.0 * zeta[phase>0] +2.0))) 
         elif self.distribution == 'oconnor_gamma':
              pass
         else:
             raise RuntimeError('unkown size distribution requested')

         #convert particle size in microns to meters
         dmode = dmode *1e-6    
         return dmode
 
class dstar_table(size_dist):
   def __init__(self,particle_parameters,radar_wavelength):
      super(dstar_table,self).__init__(particle_parameters,radar_wavelength)	   
      """class dstar_table(object):    
           def __init__(self,particle_parameters,radar_wavelength)
	   Prepare table converting:
		   dstar = radar_wavelength *(qback_radar/qback_lidar)**0.25
	   to mode diameter for particle_parameters['size_distribution'],
	   using Mie theory computed radar and lidar q-factors defined in
	   size_dist object.

           Table indexed in 1-micron steps from 0 to length of table.
           Table max index at max value of dstar--note Non-Rayleigh scattering
           may produce lower dstar values for larger mode diameters--solutions
           not unique near top end of table.
           """
     

      if 1:
            #dstar is defined as:  radar_wavelength * ( bs_radar / bs_lidar )**0.25
            # self.D is vector of diameters in microns, in 1-micron linear steps
            dstar = np.NaN * np.zeros_like(self.D)
            self.mean_qback = np.NaN * np.zeros_like(self.D)
            bs_ratio = np.NaN * np.zeros_like(self.D)
            self.p180_vs_mode_dia = np.NaN * np.zeros_like(self.D)
            self.dstar_vs_mode_dia = np.NaN * np.zeros_like(self.D)
            
      
            if 1:
                 print 'number'
                 print self.number(self.D,self.D[1000],'water')
                 import matplotlib.pylab as plt
                 plt.figure(1994)
                 plt.plot(self.D,self.number(self.D,self.D[1000],'water')
                          /np.nanmax(self.number(self.D,self.D[1000],'water')),'r'
                     ,self.D,self.number(self.D,self.D[500],'water')
                          /np.nanmax(self.number(self.D,self.D[500],'water')),'b'
                     ,self.D,self.number(self.D,self.D[200],'water')
                          /np.nanmax(self.number(self.D,self.D[200],'water')),'g')
                 plt.ylabel('number density')
                 plt.xlabel('Diameter (microns)')
                 plt.legend(('Dm=1mm','Dm=0.5mm','Dm=0.2mm'), loc='lower left')
                 ax=plt.gca()
                 ax.set_yscale('log')
                 ax.set_xscale('linear')
                 ax.set_ylim(1e-6,1.1)
                 ax.grid(True)





                 
                 plt.figure(1995)
                 
                 plt.plot(self.D,self.area_weighted_number(self.D,self.D[1000],'water')
                          /np.max(self.area_weighted_number(self.D,self.D[1000],'water')),'r'
                     ,self.D,self.area_weighted_number(self.D,self.D[500],'water')
                          /np.max(self.area_weighted_number(self.D,self.D[500],'water')),'b'
                     ,self.D,self.area_weighted_number(self.D,self.D[200],'water')
                          /np.max(self.area_weighted_number(self.D,self.D[200],'water')),'g')
                 plt.ylabel('number density * D^2')
                 plt.xlabel('Diameter (microns)')
                 plt.legend(('Dm=1mm','Dm=0.5mm','Dm=0.2mm'), loc='lower left')
                 ax=plt.gca()
                 ax.set_yscale('log')
                 ax.set_xscale('linear')
                 ax.set_ylim(1e-6,1.1)
                 ax.grid(True)
                

            #loop over values of the mode diameter
            print
            print 'starting to compute dmode vs dstar table'
            stime = datetime.utcnow()
            valid_entry = 1
            for i in range(len(self.D)):
                #compute dstar as function of mode diameter from mie theory
                #for gamma size distribution
                number_density = self.area_weighted_number(self.D,self.D[i],'water')
            
                bs_ratio[i] = np.nansum(number_density * self.qback_radar) \
                           /np.nansum(number_density * self.qback_lidar)             
      
                
                #note: radar_wavelength in microns, yeilds dstar results in microns
                self.dstar_vs_mode_dia[i] = self.radar_wavelength * bs_ratio[i]**0.25

                #compute p180/4pi as a function of mode diameter
                self.p180_vs_mode_dia[i] = np.nansum(number_density * self.qback_lidar) \
                                     /(np.nansum( number_density * self.qsca_lidar))

                #compute mean backscatter efficeincy as fucntion of mode diameter 
                self.mean_qback[i] = np.nansum(number_density * self.qback_lidar) \
                                      /np.nansum(number_density)                    
                                     
                #if more than 2% of size distribution is larger than 2 mm dstar-->mode is not valid
                stop_table_at = 3000 #don't extend table beyound this size in microns if distribution
                stop_table_limit = 0.1
                if number_density[stop_table_at] > np.nanmax(number_density) * stop_table_limit :
                     if valid_entry == 1:
                         print
                         print 'dstar_vs_mode_dia table terminated at mode_dia = ',self.D[i]/1000.0,' mm'
                         print ' more than ',stop_table_limit *100,'% of size distribution >'\
                                      , stop_table_limit /1000.0,' mm in diameter'
                         print
                         valid_entry = 0
                     self.dstar_vs_mode_dia[i] = np.NaN
                     self.p180_vs_mode_dia[i]  = np.NaN
                     self.mean_qback[i] = np.NaN                
      
            print 'time to create dmode vs dstar table ',datetime.utcnow()-stime


            
            if 1:
                 """
                 import matplotlib.pylab as plt
                 plt.figure(1996)
                 plt.plot(D/1000.0,radar_weighted_diameter)
                 plt.ylabel('radar weighted diamter')
                 plt.xlabel('Diameter (mm)')
                 ax=plt.gca()
                 ax.set_xscale('log')
                 ax.set_yscale('log')
                 ax.grid(True)
                 """
                 alpha = str(particle_parameters['size_distribution']['alpha_water'])
                 gam = str(particle_parameters['size_distribution']['g_water'])
                 labl= 'a='+alpha+' g='+gam
                 plt.figure(1996)
                 plt.plot(self.D/1000.0,self.p180_vs_mode_dia,'b'
                               ,self.D/1000.0,self.mean_qback,'r'
                               ,self.D/1000.0,2.0*self.p180_vs_mode_dia,'c')
                 plt.ylabel('')
                 plt.xlabel('Mode diameter (mm)')
                 ax=plt.gca()
                 ax.set_xscale('log')
                 plt.legend('p180','<Qbs>')
                 ax.grid(True)
                 
                 plt.figure(1997)
                 plt.plot(self.D,self.qback_radar/self.qback_lidar,label=labl)
                 plt.ylabel('qback_radar/qback_lidar')
                 plt.xlabel('Diameter (microns)')
                
                 ax=plt.gca()
                 plt.xlim(10,1e4)
                 plt.legend(loc='lower right')
                 ax.set_xscale('log')
                 ax.set_yscale('log')
                 ax.grid(True)
                 
                 plt.figure(1998)
                 plt.plot(self.D/1000.0,bs_ratio,label=labl)
                 plt.ylabel('qback_radar/qback_lidar')
                 plt.xlabel('Mode Diameter (mm)')
                 plt.legend(loc='lower right')
                 ax=plt.gca()
                 ax.set_xscale('log')
                 ax.set_yscale('log')
                 plt.xlim(8e-3,10)
                 plt.ylim(1e-10,10)
                 ax.grid(True)
                 
                 plt.figure(1999)
                 plt.plot(self.D,self.dstar_vs_mode_dia,label=labl)
                 plt.xlabel('mode diameter (microns)')
                 plt.ylabel('Dstar')
                 plt.legend(loc='lower right')
                 plt.xlim(5,1e4)
                 plt.ylim(10,1e4)
                 ax=plt.gca()
                 ax.set_xscale('log')
                 ax.set_yscale('log')
                 ax.grid(True)
                
            max_mode_dia = np.nanargmax(self.dstar_vs_mode_dia)
          
            dstar = np.arange(np.int(self.dstar_vs_mode_dia[max_mode_dia]))
              
            self.dstar_to_dmode_table = np.interp(
	          dstar,self.dstar_vs_mode_dia[:max_mode_dia],self.D[:max_mode_dia])
            self.max_index = len(self.dstar_to_dmode_table)-1
	   
	    
            if 1:
                 import matplotlib.pylab as plt
                 f=plt.figure(2000)
                 ax1 = f.add_subplot(111)
                 plt.plot(dstar,self.dstar_to_dmode_table,label=labl)
                 plt.xlabel('Dstar (microns)')
                 plt.ylabel('Mode diameter')
                 plt.legend(loc='upper left')
                 plt.xlim(10,1e4)
                 ax = plt.gca()
                 ax.set_xscale('log')
                 ax.set_yscale('log')
                 ax.grid(True)
                 
                 ax2 = ax1.twinx()
                 print dir(ax2)             


                 #deff_prime = d_mode * (g/alpha) * (Gamma((alpha+7)/g)/Gamma((alpha+3)/g)**0.25
                 y_limits = ax1.get_ylim()
                 alpha = particle_parameters['size_distribution']['alpha_water']
                 gam   = particle_parameters['size_distribution']['g_water']
            
               
                 #convert mode diameters to effective_diameters prime
                 
                 dm_to_defp = (gam/alpha)**(1.0/gam)\
                           * (gm.gamma((alpha + 7.0)/gam)\
                           /gm.gamma((alpha + 3.0)/gam))**0.25
                 print 'dmode to deff prime = ',dm_to_defp          
                 y_limits = np.asarray(y_limits) * dm_to_defp
                 ax2.set_ylim(y_limits)

                 ax2.set_yscale('log')
                 ax2.set_ylabel('Effective diameter prime')
                 # ax2.set_yticklabels(ax.get_yticklabels())

     
                 plt.figure(2001)
                 plt.plot(self.D,self.p180_vs_mode_dia)
                 plt.ylabel('P(180)/4pi')
                 plt.xlabel('Mode diameter')
                 ax=plt.gca()
                 ax.set_xscale('log')
                 ax.grid(True)
                 
          
            

   def dmode_from_radar_lidar_mie(self,bs_radar,bs_lidar):
       """dmode_from_radar_lidar_mie(self,bs_radar,bs_lidar)
          use radar/lidar backscatter ratio and mie theory
	  to compute the mode diameter, lidar extinction and backscatter efficeincy for
	  the drop size distribution described in particle_parameters
          returns---mode diameter, backscatter cross section, backscatter efficeincy,
          and dstar.
	  """
     
       #convert dstar in microns to an integer index into dstar_to_dmode_table

       dstar = self.radar_wavelength * (bs_radar/bs_lidar)**0.25
       index = dstar.copy() # * 1e6 #index in microns
       index[np.isnan(index)]=0
       index[index <0] = 0
   
       #if dstar is too large for table set index to largest value in table
       index[index > len(self.dstar_to_dmode_table)] = len(self.dstar_to_dmode_table)-1
       index = index.astype(int)
    
        
       #read mode_diameters from dstar_to_dmode_table and covert micron values to meters
       mode_diameter = self.dstar_to_dmode_table[index]*1e-6

       #convert index zeros back to NaN's
       mode_diameter[index == 0] = np.NaN
       
       #redefine index as mode_diameter to address p180 table
       index = (mode_diameter*1e6).copy()
       index[np.isnan(index)]= 0
       index = index.astype(int)
       index[index > len(self.p180_vs_mode_dia)] = len(self.p180_vs_mode_dia)-1
       #compute lidar ext
       ext_lidar = bs_lidar/self.p180_vs_mode_dia[index]
       mode_diameter[index == 0] = np.NaN
       #return mode diameter, lidar extinction, and mean backscatter efficeincy
      
       return mode_diameter, ext_lidar,self.mean_qback[index],dstar
   


def liquid_water_content_mie(effective_diameter,beta_a_backscat,mean_qback):
    """
       liquid_water_content_mie(effective_diameter,beta_a_backscat,mean_qback)
       for liquid parts of scene.

       effective_diameter = array of values calculated from size distribution and radar/lidar ratio
       beta_a_backscat    = array of lidar backscatter cross sections, 1/(m sr)
       mean_qback         = array of backscatter efficeincies computed from mie and size distribution
                            as function of mode diameter, mode_diameter=1-->len(mean_qback) in microns 
       
       volume of droplet = (effective_diameter/24) * area of droplet
       volume of water per unit volume air = effective_diameter * total_area_of_droplets per unit volume
                                           = effective_diameter * beta_a_backscat/mean_qback
                                           
    """
    LWC = np.NaN * np.zeros_like(beta_a_backscat)
      
    density_of_water = 1000    # kg/m^3
    """
    p180 = 1/(4pi) * < Q_back >/< Q_sca >
    # this because of radar 4pi backscatter definition used in mie calculations
    """
    LWC = 2.0 * density_of_water * (effective_diameter / 3.0) \
                     * (beta_a_backscat / mean_qback)
    #i = int(beta_a_backscat.shape[0]/2.0)
    return LWC

def liquid_water_content_ext_approx(LWC,effective_diameter,extinction,phase):
       """liquid_water_content(deff,extinction,phase)
       effective_diameter = 3rd moment over 2nd moment of distribution (m)
                            can be an array 
       extinction         = particle extinction cross section (1/m)
                            shape must match effective_diameter
       phase              = particle phase (0 = water, 1 = ice)
                            shape must match effective_diameter
       returns liquid water content (gr/m^3)"""

       specific_gravity_of_ice = 0.91
       specific_gravity_of_water = 1.0
       #extinction = 2 * (number of particles per unit volume) * (average area of a particle)
       #extinction = beta = 2 * N * pi/4 *<D^2>
       #N * <D^2> = 2 * beta /pi  
       #particle mass per unit volume = particle_density *  * effective_diameter * extinction / 2.0
       #effective diameter = Deff = <D^3>/<D^2>
       # N * <D^3> = N * Deff *<D^2> = Deff * 2 * beta /pi
       #LWC = rho * N * pi/6.0 * <D^3> = rho * pi/6 * Deff * 2 * beta /pi = 1/3 * Deff * beta

       #for water density of 1000 kg/m^3, the liquid water content in kg/m^3 
  
       LWC[phase==0] = 333.33 * specific_gravity_of_water \
             * effective_diameter[phase==0] * extinction[phase==0]
       
       LWC[phase > 0] = 333.33 * specific_gravity_of_ice \
             * effective_diameter[phase > 0] * extinction[phase > 0]   #kg/m^3
       return
   
def ice_water_content(LWC,effective_diameter,extinction,phase):
    """liquid_water_content(deff,extinction,phase)
       effective_diameter = 3rd moment over 2nd moment of distribution (m)
                            can be an array 
       extinction         = particle extinction cross section (1/m)
                            shape must match effective_diameter
       phase              = particle phase (0 = water, 1 = ice)
                            shape must match effective_diameter
       returns liquid water content (gr/m^3)"""

    specific_gravity_of_ice = 0.91
       
    #extinction = 2 * (number of particles per unit volume) * (average area of a particle)
    #extinction = beta = 2 * N * pi/4 *<D^2>
    #N * <D^2> = 2 * beta /pi  
    #particle mass per unit volume = particle_density *  * effective_diameter * extinction / 2.0
    #effective diameter = Deff = <D^3>/<D^2>
    # N * <D^3> = N * Deff *<D^2> = Deff * 2 * beta /pi
    #LWC = rho * N * pi/6.0 * <D^3> = rho * pi/6 * Deff * 2 * beta /pi = 1/3 * Deff * beta

    #for water density of 1000 kg/m^3, the liquid water content in kg/m^3 
  

    LWC[phase > 0] = 333.33 * specific_gravity_of_ice \
                     * effective_diameter[phase > 0] * extinction[phase > 0]   #kg/m^3


    return LWC

       



def make_fall_velocity_vs_deff_prime_tables(deff_prime,zeta,alpha,gam,air_temperature,air_density,phase):
    """make_fall_velocity_vs_deff_prime_tables(deff_prime,alpha,gam,air_temperature,air_density)
       make table of fall_velocity vs deff_prime as function of particle aspect
       ratio parameter zeta
       deff_prime = input vector providing deff_prime axis for table
       zeta       = input vector providing zeta axis for table
       alpha, gam = size distribution parameters N=No*D^alpha*exp(-bD^gam)
       returns:
           mw_table = mass-weighted fall velocity vs deff_prime and zeta
           rw_table = radar-weighted fall velocity vs deff_prime and zeta
           one dimensional tables vs deff_prime returned for phase=='water'
       """

    #if water cloud
    if phase == 'ice':
        for i in range(len(deff_prime)):
            for k in range(len(zeta)):
                [mw_table[i,k],rw_table[i,k]] = weighted_fall_velocity(deff_prime,alpha,gam\
                           ,zeta,air_temperature,rho_air,'ice')
    else:   #water
         for i in range(len(deff_prime)):
              [mw_table[i],rw_table[i,k]]=weighted_fall_velocity(deff_prime,alpha,gam\
                           ,zeta,air_temperature,rho_air,'water')
    return mw_table,rw_table       
            
def fall_velocity_vs_deff_prime(deff_prime,alpha,gam,zeta,phase):
    """computed_fall_velocity_image(deff_prime,alpha,gam,zeta,temperature,pressure,phase)
       computes the radar weighted fall velocity from measured deff_prime and and the
       assumed size distribution prarameters
       
       deff_prime  = deff_prime[time,altitude], array of measured deff_prime (m)
       alpha,gam   = size dist parameters, N(D) = No *D^alpha*exp(-(D/Dm)^gam)
                     where Dm is the mode diameter
       zeta        = aspect ratio power law  h = D^zeta
       temperature = air temperature(kelvin)
       pressure    = air pressure (mb)
       phase       = 'water' or 'ice'
       
       """
    

    radar_weighted_fall_velocity = np.zeros_like(deff_prime)
    mass_weighted_fall_velocity = np.zeros_like(deff_prime)
    
    for i in range(len(deff_prime)):
       [radar_weighted_fall_velocity[i],mass_weighted_fall_velocity[i]]= \
                   weighted_fall_velocity(deff_prime[i],alpha,gam,zeta,temperature,rho_air,phase)
    if 1:
        plt.figure(1004)
        plt.plot(deff_prime*1e6,mass_weighted_fall_velocity \
                 ,deff_prime*1e6,radar_weighted_fall_velocity)
        plt.grid(True)
        plt.xlabel('deff_prime (microns)')
        plt.ylabel('fall velocity (m/s)')
        ax=plt.gca()
        ax.set_xscale('log')

    return
                                                   
       
def weighted_fall_velocity(Dm,particle_parameters,zeta,temperature,pressure,phase,size_dist):
    """radar_weighted_fall_velocity(Dm,particle_parameters,zeta,temperature,pressure,phase,size_dist)
       inputs:
         Dm                      = mode diameter from hsrl and radar (m)
         zeta                    = power in aspect ratio power law
         particle_parameters     = dictionary, contianing particle parameters
         temperature             = air temperature (deg-K)
         pressure                = air pressure
         phase                   = 'ice' or 'water'
         size_dist               = Class contianing size distributions and methods
         
       returns:
         nw_fall_velocity        = number weighted fall velocity(m/s)  
         rw_fall_velocity        = radar cross section weighted fall velocity (m/s)
         mw_fall_velocity        = mass weighted fall velocity (m/s)
    """


   
    #vector of diameters in m from 10 microns to 1 cm
    D = np.logspace(1.0,4.0,400)/1e6
   
    dD = np.zeros_like(D)
    dD[1:-1]= (D[2:]-D[:-2])/2.0
   
    Dr = 1.0
    #constants, program works in mks 
    rho_ice=0.91 * 1000.0 #kg/m^3
    rho_water = 1000.0 #kg/m^3
    gas_constant = 287 #J/(kg K)
          
    rw_fall_velocity = np.NaN * np.zeros_like(Dm)
    mw_fall_velocity = np.NaN * np.zeros_like(Dm)
    model_spectral_width = np.NaN * np.zeros_like(Dm)
    nw_fall_velocity = np.NaN * np.zeros_like(Dm)
    
    #loop over altitudes
    for k in range(Dm.shape[1]):        
        eta = dynamic_viscosity(temperature[k])
        #note pressures are supplied in mb
        rho_air = 100.0 *pressure[k]/(gas_constant * temperature[k]) #in kg/m^3
        X_const = best_constant(rho_air,eta)
        #loop over times
        for i in range(Dm.shape[0]):
                if not np.isnan(Dm[i,k]) and phase[i,k] == 1: #for ice phase
                    #gamma_ratio = math.gamma((alpha_ice + 3.0) / g_ice)\
                    #         / math.gamma((2.0 * zeta[i,k] + alpha_ice + 5.0) / g_ice)
                                  
                    #best number---equation 8 Mitchell and Heymsfiled 2005
                    #this is vector of len(D)
                    #note inputs supplied in mks units 
                    #X = best_number(D,Dr,zeta[i,k],rho_ice,rho_air,eta)
                    X = best_number(D,Dr,X_const,zeta[i,k],rho_ice)
   
                    #vector of fall velocities converting from cm/s to m/s
                    Vf = spheroid_fall_velocity(D,X,zeta[i,k],rho_air,eta,1.0)
                  
                    #radar weighted fall velocity
                    dist = size_dist.radar_weighted_sizes(D,dD,Dm[i,k],'ice',zeta[i,k])
                    rw_fall_velocity[i,k] = np.nansum(Vf * dist *dD)

                    #spectral width
                    model_spectral_width[i,k] = 2.0 * np.nansum((Vf - rw_fall_velocity[i,k])**2 \
                          * dist * dD)
                
                    #mass weighted fall velocity
                    dist = size_dist.mass_weighted_sizes(D,dD,Dm[i,k],'ice',zeta[i,k]) 
                    mw_fall_velocity[i,k] = np.nansum(Vf * dist *dD)

                    #number weighted fall velocity
                    #convert mass weighted distribution to number weighted
                    dist = dist * D**-(2.0 +zeta[i,k])
                    dist = dist /np.nansum(dist * dD)
                    nw_fall_velocity[i,k] = np.nansum(Vf * dist * dD)

                elif not np.isnan(Dm[i,k]): #for water phase
                                  
                    #best number---equation 8 Mitchell and Heymsfiled 2005
                    #this is vector of len(D)
                    #X = best_number(D,Dr,zeta[i,k],rho_water,rho_air,eta)
                    X = best_number(D,Dr,X_const,1.0,rho_water)
                    #vector of fall velocities for water drops
                    Vf = spheroid_fall_velocity(D,X,rho_water,rho_air,eta,0.0)
                   
                    #radar weighted fall velocity
                   
                    dist = size_dist.radar_weighted_sizes(D,dD,Dm[i,k],'water')
                    rw_fall_velocity[i,k] = np.nansum(Vf * dist * dD)
        
                    #spectral width
                    model_spectral_width[i,k] = 2*(np.sqrt(np.nanmean((Vf**2 * dist * dD)) \
                                    - (np.mean(Vf * dist * dD))**2))

                    #convert radar weighted sizes to mass weighted sizes
                    #from 6th to 3rd power weighting
                    dist = dist * D**-3.0
                    dist = dist / np.nansum(dist *dD)
                    
                    #mass weighted fall velocity
                    mw_fall_velocity[i,k] = np.nansum(Vf * dist * dD)
                    
                    #number weighted fall velocity
                    #convert mass weighted distribution to number weighted
                    dist = dist * D**-3.0
                    dist = dist /np.nansum(dist * dD)
                    nw_fall_velocity[i,k] = np.nansum(Vf * dist * dD)
                   
                else:
                   rw_fall_velocity[i,k]=np.NaN
                   mw_fall_velocity[i,k]=np.NaN
                   model_spectral_width[i,k]=np.NaN
                   nw_fall_velocity[i,k] = np.NaN
        
    return rw_fall_velocity,mw_fall_velocity,model_spectral_width,nw_fall_velocity



def spheroid_fall_velocity(D,X,zeta,rho_air,eta,phase):
    """spheriod_fall_velocity_vs_D(D,X,zeta,rho_air,phase)
       calculate the fall velocity of oblate spheroids whose
       height is described by power law h=sigma_h*D_r*(D/D_r)^eta
             where we assume D_r = 1 micron
             and sigma_h = 1 (i.e. 1 micron particles are spheres)
       D                       = vector of diameters, D (m)
       X                       = best number
       zeta                    = power in aspect ratio power law
       temperature             = air temperature (deg-K)
       rho_air                 = air density (kg/m^3)
       eta                     = dynamic viscosity (kg/(m sec))
       phase                   = 'ice' or 'water'                  
       fall_velocity           = particle fall velocity (m/s) """


    
    if phase=='ice':
        delta_0=5.83
        C_0=0.6 #page 4345 Khvorostyanov and Curry, Dec 2005 JAS
    else:     #phase == 'water':
        delta_0 = 9.06
        C_0 = 0.292
       
    C_1 = 4 / (np.sqrt(C_0) * delta_0**2)
    C_2 = 0.25 * delta_0**2
    a_0 = 1.7e-3                  #page 4345 K&C
    b_0 = 0.8

    # X=best number vector--from equation 8 Mitchell and Heymsfiled 2005
    #X_const=4 * aspect_ratio * rho_particle * rho_air * g / (3 * eta**2)
    #X=X_const * D**3

    #X_const = 4 * Dr**(1-zeta) * rho_particle * rho_air * g / (3 * eta**2)
    #X = X_const * D **(2 + zeta) 
    
    #fall velocities for single particles
    #a_0=0
    #fall_velocity = (eta / (rho_air * D)) \
    #      * (C_2 * (((1 + C_1 * X**0.5)**0.5) - 1)**2 - a_0 * X**b_0)
    
    fall_velocity = (eta / (rho_air * D)) \
          * (C_2 * (((1 + C_1 * X**0.5)**0.5) - 1)**2)
    return fall_velocity #in m/sec

def best_number(D,Dr,X_const,zeta,rho_particle):
    """ X=best number vector--from equation 8 Mitchell and Heymsfiled 2005
        X_const=4 * aspect_ratio * rho_particle * rho_air * g / (3 * eta**2)
        X=X_const * D**3
         
       
        D            = particle diameter (m)
        zeta         = aspect ratio power,  h = (D/Dr)^zeta
        Dr           = reference diameter where particles are spheres (m)
        rho_particle = particle density(kg/m^3)
        rho_air      = air density (kg/m^3)
        eta          = dynamic viscosity (gr/(cm sec))
        """
    g = 9.80  #acceleration gravity (m/s^2)
    # X=best number vector--from equation 8 Mitchell and Heymsfiled 2005
    #X_const=4 * aspect_ratio * rho_air * g / (3 * eta**2)

    
    #best number 
    X = X_const * rho_particle * D**(2 + zeta) * Dr**(1-zeta)
    return X

def best_constant(rho_air,eta):
    """ X=best number vector--from equation 8 Mitchell and Heymsfiled 2005
        X_const=4 * aspect_ratio * rho_particle * rho_air * g / (3 * eta**2)
        X=X_const * D**3
         
       
        D            = particle diameter (m)
        zeta         = aspect ratio power,  h = (D/Dr)^zeta
        Dr           = reference diameter where particles are spheres (m)
        rho_particle = particle density(kg/m^3)
        rho_air      = air density (kg/m^3)
        eta          = dynamic viscosity (gr/(cm sec))
        """

    g = 9.80  #acceleration gravity (m/s^2)
   
  
    # X=best number vector--from equation 8 Mitchell and Heymsfiled 2005
    #X_const=4 * aspect_ratio * rho_particle * rho_air * g / (3 * eta**2)
    #X=X_const * D**3

    #best number with diameters and densities converted from mks to cgs
    X_const = 4 * rho_air * g / (3 * eta**2)
    return X_const

def dynamic_viscosity(temperature):
    #dynamic viscosity in poise, gr/(cm sec), K&C page 4347
    #temperature in Kelvin
    if temperature-273 < 0:
        eta = 1.718e-4*(1+0.00285*(temperature-273)-6.9e-6*(temperature-273)**2)
    else:
        eta = 1.718e-4*(1+0.00285*(temperature-273))
    return eta * 0.1 #convert from gr/(cm sec) to kg/(m sec)

if __name__ == '__main__':

    
    #D = np.arange(1000)/1e6  #diameters in microns
    #zeta = .6
    #temperature = 280.
    #rho_air = 1.2
    phase = 'ice'

    #fall_velocity = spheroid_fall_velocity(D,zeta,temperature,rho_air,phase)
    

    
    alpha = 2.0
    gam =1.0
    #deff_prime = 1.2e-3
    
    #[rw_fall_velocity,mw_fall_velocity] = weighted_fall_velocity(deff_prime,alpha,gam,zeta,temperature,rho_air,phase)

    #print 'radar weighted fall velocity = ',rw_fall_velocity
    #print 'mass weighted fall velocity = ',mw_fall_velocity
    if 0:
        #particle diameters
        D = np.logspace(-6.0,-2.0,5)
        Dr = 1e-6
        zeta=0.6
        gamma_ratio = math.gamma((alpha + 3.0) / gam) / math.gamma((2.0 * zeta + alpha + 5.0) / gam)

        deff_prime = D * (Dr**((1.0-zeta)/(zeta +1.0))*(alpha / gam)**(-1.0/gam)* gamma_ratio**(-1.0/(2.0*zeta+2.0)))**((zeta+1.0)/2.0)

        Dm = Dr**((zeta - 1) / (zeta + 1)) * (alpha / gam)**(1/gam) \
           *gamma_ratio**(1 / (2 * zeta + 2))*deff_prime**(2 / (zeta + 1))
       
    
    
        fall_velocity_vs_deff_prime(deff_prime,alpha,gam,zeta,phase)
       
  

    if 1:   #weighted fall velocities
        zeta = 1.0
        alpha = 2.0
        gam  = 1.0
        phase = 'water'
        deff_prime =np.ones((2,4)) * 1e-3
        Dr = 1e-6
        rho_air =1.2 #kg/m^3
        temperature =280  #deg kelvin
        
        [ms_fv,rw_fv,sp_w]=weighted_fall_velocity(deff_prime,alpha,gam,zeta,temperature,rho_air,phase)

