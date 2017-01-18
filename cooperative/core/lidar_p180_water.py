import numpy as np
import os
import lg_base.core.read_utilities as ru
from lg_base.core.locate_file import locate_file

class lidar_p180_water(object):
    def __init__(self,wavelength,particle_parameters):
        """calculate size distribution weighted p180/4pi from mie q-factors contained in
        'cooparative/config/mie_coefficients.txt' file"""
        
        #if not particle_parameters['p180_water']['size_distribution'] == 'None':
        if not particle_parameters['size_distribution'] == 'None': 
            mie_filename = locate_file('mie_coefficients_0_18000.txt')
            if os.path.isfile(mie_filename):
                [hh,dd]=ru.readascii(mie_filename)
            else:
                print "unable to preform adaptive p180 operation"
                raise RuntimeError, ("Can't find mie results in '",mie_filename,"'")
            """
            alpha = particle_parameters['alpha_water']
            gam   = particle_parameters['g_water']
            """
        
            if particle_parameters['size_distribution'].has_key('modified_gamma'):
                 alpha = particle_parameters['size_distribution']['modified_gamma']['alpha_water']
                 gam   = particle_parameters['size_distribution']['modified_gamma']['g_water']
            elif particle_parameters['size_distribution'].has_key('oconnor_gamma'):
                 alpha = particle_parameters['size_distribution']['oconnor_gamma']['alpha_water']
            else:
                 raise RuntimeError("unknown particle size distribution")
             
            #convert wavenumber to diameter in microns
            self.D = dd[:,0]*wavelength*1e6/np.pi
            #q mie scattering
            qsca  = dd[:,2]
            #q mie backscatter
            qback = dd[:,3]

            #p180/4pi vs diameter table 
            self.p180_vs_d_table =np.zeros_like(self.D)

            #if particle_parameters['p180_water']['size_distribution'] == 'modified_gamma':
            if particle_parameters['size_distribution'].has_key('modified_gamma'):
                for i in range(1,len(self.D)-1):    
                    self.p180_vs_d_table[i] = np.sum(self.modified_gamma(self.D,self.D[i],alpha,gam)\
                        * qback)/(4.0*np.pi)\
                        /np.sum(self.modified_gamma(self.D,self.D[i],alpha,gam)*qsca)
            elif particle_parameters['size_distribution'].has_key('oconnor_gamma'):
                for i in range(len(self.D)):    
                   self.p180_vs_d_table[i] = np.sum(self.oconnor_gamma(self.D,self.D[i],alpha)\
                      * qback)/(4.0*np.pi)\
                      /np.sum(self.oconnor_gamma(self.D,self.D[i],alpha)*qsca)
            else:
                raise RuntimeError("unrecongized size distribution type"\
                      + particle_parameters.p180_water)
           
        else:  #constant value of p180/4pi requested
            self.D = None
            self.p180_vs_d_table = particle_parameters['p180_water']['value']
        return 

    
    def oconnor_gamma(D,Dm,alpha):
        """oconnor_gamma(D,D0,mu)
           O'Connors normalized gamma distribution
           D = Diameter (microns)
           D0 = median of volume weighted diameter
           Dm = mode diameter
           alpha = shape factor = mu in O'Connor

           Notice (alpha + 2) to provide particle area weighted distribution
           
    
           N = (D/D0)**(alpha + 2) * np.exp(-(3.67+alpha)*D/D0)
        """
        #convert mode diameter to median volume weighted diameter
        D0 = Dm * (3.67+alpha)/alpha
        
        N = (D/D0)**(alpha + 2) * np.exp(-(3.67+alpha)*D/D0)
        N = N/np.sum(N)
        return N
       

    def modified_gamma(self,D,Dm,alpha,gam):
       """modified_gamma(x,Dm,alpha,gam):
          D = vector of sizes 
          Dm = mode diameter
          alpha, gam = D^alpha *exp(-(alpha/gam)*(D/Dm)^gam)
          N = vector of sizes for a given Dm,alpha,gam

          Notice that distribution has been multiplied by d**2
          and qsca in order to provide a distribution
          wieghted by scattered intensity.
       """
       #area weighted size distribution
       N = D**(alpha+2) * np.exp(-(alpha/gam)*(D/Dm)**gam)
       N = N/np.sum(N)
       return N        
        
    def __call__(self,mode_diameter):
       """If constant p180/4pi was requested, mode_diameter can be any
          array with the dimensions of the data array--mode_diameter is
          not used in this case except to dimension the output"""
       
       if self.D == None:  #requested constant value
           lidar_p180_water = np.ones_like(mode_diameter) * self.p180_vs_d_table
       else:
           index = mode_diameter.copy()*1e6 #index in microns
           index[np.isnan(index)] = 0
           index[index < 0] = 0
           index[index >= self.D[-1]] = self.D[-1]
           index = index.astype(int)
           lidar_p180_water = np.zeros_like(mode_diameter)
           lidar_p180_water = self.p180_vs_d_table[index]
           lidar_p180_water[index == 0] = np.NaN
       return lidar_p180_water
