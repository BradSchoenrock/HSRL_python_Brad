import numpy as np
import os
import lg_base.core.read_utilities as ru
from lg_base.core.locate_file import locate_file


#Try to use the much faster nanmean from bottleneck, otherwise fall back
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
    
def modified_gamma(D,Dm,alpha,gam):
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

def oconnor_gamma(D,D0,alpha):
        """oconnor_gamma(D,D0,mu)
           O'Connors normalized gamma distribution
           D = Volume weighted Diameter (microns)
           D0 = median of volume weighted diameter
           alpha = shape factor = mu in O'Connor

           Notice (alpha + 2) to provide particle area weighted distribution
           
           """
        #convert mode diameter to median volume diameter
        D0 = Dm *(3.67+alpha)/alpha

        N = (D/D0)**(alpha + 2) * np.exp(-(3.67+alpha)*D/D0)
        N = N/np.sum(N)
        return N
    
class radar_p180_water(object):
    def __init__(self,wavelength,particle_parameters):
        """calculate size distribution weighted non-rayliegh radar backscatter phase
        function from mie q-factors contained in 'cooperative/config/mie_xx_temp.txt' file
        Returns a table of p180/4pi at 1 micron intervals"""
        
        if (particle_parameters['non-Rayleigh_radar'] == 'True' \
               or particle_parameters['non-Rayleigh_radar'] == 'true'):
            if wavelength == 8.6e-3:
                mie_filename = locate_file('mie_data_8.6mm_285K.txt')
                print 'reading  ' + mie_filename
                index_refraction = 4.72-2.72j
            elif wavelength == 3.154e-3:

                #hack--approximate wavelength used for mie calculation
                wavelength = 3.2e-3
            
                mie_filename = locate_file('mie_data_3.2mm_285K.txt')
                print 'reading  '+ mie_filename
                index_refraction = 3.08 - 1.78j
            else:
                print 'unknown radar wavelength = ',wavelength
                raise RuntimeError('Unknown radar wavelength')
            if os.path.isfile(mie_filename):
                [hh,dd]=ru.readascii(mie_filename)
            else:
                print "unable to preform adaptive p180 operation"
                raise RuntimeError, ("Can't find mie results in '",mie_filename,"'")


            #convert wavelength from meters to microns
            lambd = wavelength * 1e6
           
            #alpha = particle_parameters['alpha_water']
            #gam   = particle_parameters['g_water']

            if particle_parameters['size_distribution'].has_key('modified_gamma'):  
                alpha = particle_parameters['size_distribution']['modified_gamma']['alpha_water']
                gam   = particle_parameters['size_distribution']['modified_gamma']['g_water']
            elif particle_parameters['size_distribution'].has_key('oconnor_gamma'): 
                alpha = particle_parameters['size_distribution']['oconnor_gamma']['alpha_water']
                
            #find largest diameter in microns for this wavelength
            dmax = np.int(dd[-1,0] * lambd /np.pi)
            self.D = np.arange(dmax).astype('float')

            qsca = np.interp(self.D,dd[:,0]*lambd/np.pi,dd[:,2])
            qback = np.interp(self.D,dd[:,0]*lambd/np.pi,dd[:,3])

            k_squared = (index_refraction**2 -1)/(index_refraction**2 +2)
            k_squared = np.real(k_squared * np.conj(k_squared))

           
            print
            print '**********************************************************************************************'
            print
            print 'k_squared  ',   k_squared, wavelength
            print
            print
            
            #qback_theory = np.pi**5 * lambd**-4 * k_squared * self.D**4
            qback_theory = np.pi**3 * lambd**-4 * k_squared * self.D**4
            
            #hack
            qback_theory = qback_theory * 4.0*np.pi
            print
            print
            print '****** hack----qback_theory adjustment by factor of 4*pi'
            print
            
            if 1:
                 import matplotlib.pylab as plt
                 plt.figure(1999)
                 plt.plot(self.D,qback,'k',self.D,qback_theory,'r')
                 plt.xlabel('diameter (microns)')
                 plt.ylabel('q backscatter')
                 ax=plt.gca()
                 ax.set_xscale('log')
                 ax.set_yscale('log')
                 ax.grid(True)
            #if particle_parameters['p180_water']['size_distribution'] == 'modified_gamma':
            if particle_parameters['size_distribution'].has_key('modified_gamma'):
                print 'Non-Rayleigh radar corrections calcualted with modified gamma size distribution'
                D_eval = np.logspace(np.log10(self.D[15]),np.log10(self.D[-1]),500)        
                table =np.ones(len(D_eval))
                scale     = qback[100]/(self.D[100]*np.pi/lambd)**4   
                rayleigh = scale * (self.D*np.pi/lambd)**4
                qsca_rayleigh = rayleigh
                for i in range(len(D_eval)-1):
                    #table[i] = np.nansum(modified_gamma(self.D,D_eval[i],alpha,gam)* qback) \
                    #    /(4.0* np.pi * np.nansum(modified_gamma(self.D,D_eval[i],alpha,gam) * qsca))
                    
                    #this is the ratio of the actual backscatter to the Rayleigh computed backscatter
                    #table[i] = np.nansum(modified_gamma(self.D,D_eval[i],alpha,gam) * qsca *  qback) \
                    #    /(np.nansum(modified_gamma(self.D,D_eval[i],alpha,gam) * qsca *  qback_theory))

                    #is this the correct way to average?
                    table[i] =np.nansum(modified_gamma(self.D,D_eval[i],alpha,gam) * qsca)\
                       /np.nansum(modified_gamma(self.D,D_eval[i],alpha,gam) * qsca *  qback/qback_theory) 

                    if i%50 == 0:
                        import matplotlib.pylab as plt
                        plt.figure(3000)
                        plt.plot(self.D,modified_gamma(self.D,D_eval[i],alpha,gam))
                        plt.xlabel('diameter (microns)')
                        plt.ylabel('distribution')
                        ax=plt.gca()
                        ax.set_xscale('log')
                        ax.grid(True)
                      
            elif particle_parameters['size_distribution'].has_key('oconnor_gamma'):
                print 'Non-Rayleigh radar corrections calcualted with oconnor_gamma size distribution'
                for i in range(len(self.D)):    
                   table[i] = np.nansum(self.oconnor_gamma(self.D,D[i],alpha)*rayleigh)\
                      /np.nansum(self.oconnor_gamma(self.D,D[i],alpha))
            else:
                raise RuntimeError("unrecongized size distribution type")

            #fill table for all sizes using interpolation    
            self.non_rayleigh_vs_d_table = np.interp(self.D,D_eval,table) 

        else:  #no non-rayliegh radar correction
            self.D = None
            self.non_rayleigh_vs_d_table = 1.0
            
       
        if particle_parameters['non-Rayleigh_radar'] == 'True':
            import matplotlib.pylab as plt
            plt.figure(2000)
            plt.plot(dd[:,0]*lambd/np.pi,dd[:,3],self.D,qback,self.D,1.5*rayleigh,'r')
            plt.xlabel('diameter')
            plt.ylabel('qback')
            ax=plt.gca()
            ax.set_xscale('log')
            ax.grid(True)

            plt.figure(2001)
            plt.plot(dd[:,0]*lambd/np.pi,dd[:,2],'r',self.D,qsca,self.D,rayleigh,'g')
            plt.xlabel('diameter (microns)')
            plt.ylabel('qscat')
            ax=plt.gca()
            ax.set_xscale('log')
            ax.grid(True)

            plt.figure(2002)
            plt.plot(self.D,self.non_rayleigh_vs_d_table,dd[:,0]*lambd/np.pi,dd[:,5]*8*np.pi/3.0)
            plt.xlabel('diameter (microns)')
            plt.ylabel('non-Rayleigh correction factor')
            ax=plt.gca()
            ax.set_xscale('log')
            ax.grid(True)

            plt.figure(2003)
            plt.plot(dd[:,0]*lambd/np.pi,dd[:,5]*8*np.pi/3.0,'k'
                     ,self.D,modified_gamma(self.D,20.0,alpha,gam)/np.max(modified_gamma(self.D,20.0,alpha,gam)),'m'
                     ,self.D,modified_gamma(self.D,100.0,alpha,gam)/np.max(modified_gamma(self.D,100.0,alpha,gam)),'m'
                     ,self.D,modified_gamma(self.D,550.0,alpha,gam)/np.max(modified_gamma(self.D,550.0,alpha,gam)),'m'
                     ,self.D,2.0/3.0*qback/qsca,'c'
                     ,self.D,self.non_rayleigh_vs_d_table,'r')
            plt.xlabel('diameter (microns)')
            plt.ylabel('non-Rayleigh correction factor')
            ax=plt.gca()
            ax.set_xscale('log')
            ax.grid(True)
            
        return 


        
    def __call__(self,mode_diameter):
       """Returns ratio of radar backscatter to rayleigh approximation
          If non_rayleigh correction was not requested, mode_diameter can be any
          array with the dimensions of the data array--mode_diameter is
          not used in this case except to dimension the output"""
       
       if self.D == None:  #requested constant value
           #radar_p180_water = 3.0 * np.ones_like(mode_diameter)/(8.0*np.pi)
           non_rayleigh_adjustment = 1.0
       else:
           index = mode_diameter.copy()*1e6 #index in microns
           index[np.isnan(index)] = 0
           index[index < 0] = 0
           index[index >= self.D[-1]] = self.D[-1]
           index = index.astype(int)
           non_rayleigh_adjustment = self.non_rayleigh_vs_d_table[index]
           non_rayleigh_adjustment[index == 0] = np.NaN
           non_rayleigh_adjustment =1.0/non_rayleigh_adjustment
       return non_rayleigh_adjustment
