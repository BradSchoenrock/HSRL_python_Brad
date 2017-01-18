import numpy as np
import matplotlib.pylab as plt
import lg_base.core.array_utils as hau
import math
import cooperative.core.spheroid_utilities as spu
import lg_base.graphics.graphics_toolkit as gt
import cooperative.core.effective_diameter_prime as edp



    
def precip_computation(alpha, gam, D, extinction, beta_a_backscat, eff_dia_prime, phase
                       ,temperature, pressure):
    
                    alpha_over_g = alpha/gam
                    one_over_g = 1.0/gam
                    try:
                        mode_dia = Dr**((zeta-1)/(zeta+1.0)) * (alpha_over_g)**one_over_g \
                           *(math.gamma((alpha+3.0)/gam)
                           /math.gamma((2*zeta+alpha+5.0)/gam))**(1/(2*zeta+2.0))\
                           *eff_dia_prime**(2.0/(zeta+1.0))
                    except:
                        mode_dia = np.NaN
                    print
                    print 'mode diameter = %4.0f (microns)' %(mode_dia * 1e6)  
                    alpha_over_g = alpha/gam
                    one_over_g = 1.0/gam
                    # size distribution is number per dD
                    size_dist = D**alpha * np.exp(-alpha_over_g *(D/mode_dia)**gam)
                    dist_norm = np.sum(size_dist *dD)
                    size_dist /= dist_norm

                    #mass weighted
                    mw_size_dist = D**(2+zeta+alpha)*np.exp(-alpha_over_g * (D/mode_dia)**gam)
                    mw_norm = np.sum(mw_size_dist * dD)
                    ave_vol = (4.0 *np.pi)/3.0 * mw_size_dist/dist_norm
                    mw_size_dist /= mw_norm

                    #lidar extinction weighted--area weighted
                    aw_size_dist =D**(2+alpha)*np.exp(-alpha_over_g * (D/mode_dia)**gam)
                    aw_norm = np.sum(aw_size_dist* dD)
                    ave_area = (np.pi/4.0) * aw_norm/dist_norm
                    aw_size_dist /= aw_norm
            

                    #radar backscatter weigthed --volume-sqrd weighting
                    rw_size_dist = D**(4+2*zeta+alpha)*np.exp(-alpha_over_g * (D/mode_dia)**gam)
                    rw_norm = np.sum(rw_size_dist * dD)
                    ave_vol_sqrd = (np.pi**2)/36.0 * rw_norm/dist_norm
                    rw_size_dist /= rw_norm

                    #calculate fall velocity
                    eta = spu.dynamic_viscosity(temperature) # kg/(m sec)

                    #note pressures are supplied in mb
                    rho_air = 100.0 * pressure/(gas_constant * temperature) 


                    #best number---equation 8 Mitchell and Heymsfiled 2005
                    #this is vector of len(D) --best # computed in cgs units
                    #X = spu.best_number(D ,Dr ,zeta,rho_ice,rho_air,eta)

                    X_constant = best_constant(rho_air,eta)
                    if phase == 1:
                         X = best_number(D,Dr,X_const,zeta,rho_ice)
                    else:
                         X = best_number(D,Dr,X_const,zeta,rho_water)

                    #vector of fall velocities (m/s)
                    Vf = spu.spheroid_fall_velocity(D,X,zeta,rho_air,eta,1.0)
           
                    #radar weighted fall velocity
                    exponent = 4.0+ 2 * zeta +alpha
                    dist = np.exp(-alpha_over_g * (D / mode_dia)**gam) \
                        * D**exponent
                    dist /= np.sum(np.exp(-alpha_over_g * (D / mode_dia)**gam) \
                        * D**exponent * dD)
                    rw_fall_vel = np.sum(Vf * dist *dD)
            
                    model_spectral_width = 2.0 * np.sum((Vf - rw_fall_vel)**2 \
                          * dist * dD)
                    
                    #mass weighted fall velocity
                    exponent = 2.0+ 2 * zeta +alpha
                    dist = np.exp(-alpha_over_g * (D / mode_dia)**gam) \
                        * D**exponent
                    dist /= np.sum(np.exp(-alpha_over_g * (D / mode_dia)**gam) \
                        * D**exponent * dD)
                    mw_fall_vel = np.sum(Vf * dist *dD)

                    #effective diameter from effective diameter prime
                    try:
                       effective_dia = \
                              Dr**(1-zeta) * dmode**zeta \
                              * (alpha / g_ice)**(zeta/gam) \
                              *gm.gamma((alpha + zeta +3)/gam)\
                              /gm.gamma((alpha +3)/gam)
                    except:
                        effective_dia = np.NaN
                        
                    #liquid water content
                    #for water density of 1000 kg/m^3, the liquid water content in kg/m^3 
                    #LWC= (2.0/3.0) * 1e3 * effective_diameter * (extinction / 2.0)  #kg/m^3
                    #LWC = (1.0/3.0) * 1e3 * effective_diameter * extinction   #kg/m^3   
                    if phase:
                        specific_gravity = 0.91  #ice
                        p180 = .038
                    else:
                         specific_gravity = 1.0
                         p180 = 0.05
                    LWC_ext = 333.33 * effective_dia * extinction *specific_gravity  # kg/m^3
                    LWC_p180= 333.33 * effective_dia * beta_a_backscat * specific_gravity / p180 

                    return mode_dia,effective_dia,size_dist,mw_size_dist,aw_size_dist,rw_size_dist\
                           ,mw_fall_vel,rw_fall_vel,LWC_ext,LWC_p180

def special_process(struct,t_index,z_index):
    """special_process(struct,t_index,z_index)
       sandbox for various testing and debug functions
       struct = structure containing self from the maestro
                this should provide access to all varaibles
       t_index = requested time index
       z_index = requested altitude index"""

    def get_data_at_index(struct,t_index,z_index):
        """get_data_at_point(struct,t_index,z_index):
           extract values of variables at a given array index"""

        rs_init              = struct['rs_init'].__dict__
        rs_inv               = struct['rs_inv'].__dict__
        rs_mmcr              = struct['rs_mmcr'].__dict__
        rs_spheroid_particle = struct['rs_spheroid_particle'].__dict__    

       
        sounding             = rs_init['sounding']
        
    
        ptv = hau.Time_Z_Group()

        particle_parameters = rs_spheroid_particle['particle_parameters']
    
        ptv.beta_ext_radar  = np.NaN * np.zeros((1,1))
        ptv.beta_ext_lidar  = np.NaN * np.zeros((1,1))
        ptv.beta_a_backscat = np.NaN * np.zeros((1,1))
        ptv.radar_spectral_width   = np.NaN * np.zeros((1,1))
        ptv.MeanDopplerVelocity = np.NaN * np.zeros((1,1))
        ptv.zeta            = np.NaN * np.zeros((1,1))
        ptv.phase           = np.NaN * np.zeros((1,1))
        ptv.temps           = np.NaN * np.zeros(1)
        ptv.pressures       = np.NaN * np.zeros(1)
        ptv.eff_diameter_prime\
                               = np.NaN * np.zeros((1,1))
        ptv.dmode              = np.NaN * np.zeros((1,1))
        ptv.effective_diameter = np.NaN * np.zeros((1,1))
        ptv.mean_diameter      = np.NaN * np.zeros((1,1))
        ptv.mean_mass_diameter = np.NaN * np.zeros((1,1))
        ptv.LWC_ext            = np.NaN * np.zeros((1,1))
        ptv.LWC_p180           = np.NaN * np.zeros((1,1))
        ptv.hsrl_radar_precip_rate=np.NaN * np.zeros((1,1))
        ptv.rw_fall_velocity   = np.NaN * np.zeros((1,1))
        ptv.mw_fall_velocity   = np.NaN * np.zeros((1,1))
        ptv.model_spectral_width\
                               = np.NaN * np.zeros((1,1))                    

        ptv.beta_ext_lidar[0,0]  = rs_inv['beta_a'][t_index,z_index]
        ptv.beta_a_backscat[0,0] = rs_inv['beta_a_backscat'][t_index,z_index]
        ptv.beta_ext_radar[0,0] \
                        = rs_mmcr['Backscatter'][t_index,z_index] * 8 * np.pi/3.0
        ptv.radar_spectral_width[0,0]   = rs_mmcr['SpectralWidth'][t_index,z_index]
        ptv.MeanDopplerVelocity[0,0] \
                                 = rs_mmcr['MeanDopplerVelocity'][t_index,z_index]
        ptv.zeta[0,0]            = rs_spheroid_particle['zeta'][t_index,z_index]
        ptv.phase[0,0]           = rs_spheroid_particle['phase'][t_index,z_index]
        ptv.eff_diameter_prime[0,0] \
               = rs_spheroid_particle['effective_diameter_prime'][t_index,z_index]
        ptv.temps[0]             = sounding.temps[z_index]
        ptv.pressures[0]         = sounding.pressures[z_index]

        ptv.dmode[0,0]              = rs_spheroid_particle['mode_diameter'][t_index,z_index]
        ptv.effective_diameter[0,0] = rs_spheroid_particle['effective_diameter'][t_index,z_index]
        ptv.mean_diameter[0,0]      = rs_spheroid_particle['mean_diameter'][t_index,z_index]
        ptv.mean_mass_diameter[0,0] =rs_spheroid_particle['mean_mass_diameter'][t_index,z_index]
        ptv.LWC_ext[0,0]            = rs_spheroid_particle['LWC'][t_index,z_index]
        ptv.hsrl_radar_precip_rate[0,0]=rs_spheroid_particle['hsrl_radar_precip_rate'][t_index,z_index]
        ptv.rw_fall_velocity[0,0]   = rs_spheroid_particle['rw_fall_velocity'][t_index,z_index]
        ptv.mw_fall_velocity[0,0]   = rs_spheroid_particle['mw_fall_velocity'][t_index,z_index]
        ptv.model_spectral_width[0,0]= rs_spheroid_particle['model_spectral_width'][t_index,z_index]
                                          
        return ptv,particle_parameters
    
    def particle_size_distribution(struct,t_index,z_index):
        """plot particle size distribution"""

        ptv,particle_parameters = get_data_at_index(struct,t_index,z_index)


        alpha_water =  particle_parameters['alpha_water']
        alpha_ice   =  particle_parameters['alpha_ice']
        g_ice = particle_parameters['g_ice']
        g_water = particle_parameters['g_water']
        Dr      = particle_parameters['Dr']
        print 'Default processing'
        print 'effective_diameter_prime = %4i (microns)' %(ptv.eff_diameter_prime[0,0]*1e6)
        print 'effective_diameter       = %4i (microns)' %(ptv.effective_diameter[0,0]*1e6)
        print 'mode diameter            = %4i (microns)' %(ptv.dmode[0,0]*1e6)
        print 'mean diameter            = %4i (microns)' %(ptv.mean_diameter[0,0]*1e6)
        print 'Doppler velocity         = %4.2f (m/s)' %ptv.MeanDopplerVelocity[0,0]
        print 'rw_fall_vel              = %4.2f (m/s)' %ptv.rw_fall_velocity[0,0]
        print 'mw_fall_vel              = %4.2f (m/s)' %ptv.mw_fall_velocity[0,0]
        print 'model_spectral_width     = %4.2f (m/s)' %ptv.model_spectral_width[0,0]
        print 'radar_spectral_width     = %4.2f (m/s)' %ptv.radar_spectral_width[0,0]
        print 'alpha_water = ' ,alpha_water
        print 'alpha_ice   = ',alpha_ice
        print 'gamma_water = ',g_water
        print 'gamma_ice   = ',g_ice
        print 'zeta        = ',ptv.zeta[0,0]
    
        rho_ice = 0.91 *1000.0 # kg/m^3
        gas_constant = 287 #J/(kg K)
        #vector of diameters in m from 1 micron to 1 cm
        D = np.logspace(0.0,4.0,200)/1e6
        dD = np.zeros_like(D)
        dD[1:] = D[1:]-D[:-1]

    
        if ptv.phase[0,0]:
            alpha = alpha_ice
            gam = g_ice        
        else:    
             alpha = alpha_water
             gam = g_water
             
        colors = ['b','g','r','c']
        lwidth = 1.0
        while 1:
             #recalculate size distribution with new alpha?
            alpha_gamma_str = raw_input('new params:  alpha,gamma,zeta = ? or map ')
            #strip all blanks
            #alpha_gamma_str.lstrip()
            
            alpha_gamma_str.replace(" ","")
            n_alpha = 1
            n_gam = 1
            print alpha_gamma_str
            if not alpha_gamma_str == 'map':
                n_times=[2,2]
                if len(alpha_gamma_str) == 0:
                    break                
                elif alpha_gamma_str.find(',')>= 0:
                    index = alpha_gamma_str.find(',')
                    index2 = alpha_gamma_str[(index+1):].find(',')
                    index2 = index +index2+1
                    alpha = np.float(alpha_gamma_str[:index])
                    gam = np.float(alpha_gamma_str[index +1:index2])
                    ptv.zeta[0,0] = np.float(alpha_gamma_str[index2+1:len(alpha_gamma_str)])
            else:
                n_alpha = 20
                n_gam   = 40
                n_times = [n_alpha,n_gam]
            cost_func = np.zeros((n_alpha,n_gam))
            print n_times, n_times[0]
            print n_alpha,n_gam,ptv.zeta[0,0]
            if alpha_gamma_str == 'map':
               for i in range(1,n_times[0]):
                  for k in range(1,n_times[1]):
                    alpha = np.float(i)/2.0
                    gam =   .05 * np.float(k)
                    print
                    print 'alpha = ',alpha,' gamma= ',gam

                    particle_parameters['alpha_water'] = alpha
                    particle_parameters['alpha_ice'] = alpha
                    particle_parameters['g_water'] = gam
                    particle_parameters['g_ice'] = gam

                    beta_ext_radar = np.NaN * np.zeros((1,1))
                    beta_ext_lidar = np.NaN * np.zeros((1,1))
                    zeta =np.NaN * np.zeros((1,1))
                    beta_ext_lidar = rs_inv['beta_a_backscat'][t_index,z_index]/0.05

                    beta_ext_radar = rs_mmcr['Backscatter'][t_index,z_index] * 8 * np.pi/3.0
                    zeta            = rs_spheroid_particle['zeta'][t_index,z_index]

                    print beta_ext_lidar.shape
                    dmode = spu.dmode_from_lidar_radar(beta_ext_lidar,beta_ext_radar,particle_parameters
                           ,ptv.zeta,ptv.phase,8.6e-3)
                    print dmode
                    
                    
                    mode_dia,effective_dia,size_dist,mw_size_dist,aw_size_dist,rw_size_dist\
                      ,mw_fall_vel,rw_fall_vel,LWC_ext,LWC_p180 = precip_computation(\
                      alpha
                      ,gam
                      ,D
                      ,extinction
                      ,beta_a_backscat
                      ,eff_dia_prime
                      ,phase  
                      ,temperature
                      ,pressure)
            else: #single value not mapping cost function
                particle_parameters['alpha_water'] = alpha
                particle_parameters['alpha_ice'] = alpha
                particle_parameters['g_water'] = gam
                particle_parameters['g_ice'] = gam

               
                
                
                #eff dia prime from radar and lidar
                #ptv.eff_diameter_prime[0,0] = edp.lidar_radar_eff_diameter_prime(
                #      ptv.beta_ext_lidar
                #     ,ptv.beta_ext_radar
                #     ,ptv.phase
                #     ,8.6e-3)
               
                print 'source',particle_parameters['ext_source']
                print ptv.beta_a_backscat
                print particle_parameters['p180_ice']
                print ptv.beta_ext_radar
                print ptv.zeta
                print ptv.phase
                if ptv.phase[0,0] == 1:
                    p180 = particle_parameters['p180_ice']
                else:
                    p180 = 0.05
                print 'beta_ext' , ptv.beta_a_backscat/p180
                
                if particle_parameters['ext_source'] == 'bs/p180':
                    ptv.dmode = spu.dmode_from_lidar_radar(
                             ptv.beta_a_backscat/p180
                            ,ptv.beta_ext_radar
                            ,particle_parameters
                            ,ptv.zeta
                            ,ptv.phase
                            ,8.6e-3)
                else:    
                    ptv.dmode = spu.dmode_from_lidar_radar(
                            ptv.beta_ext_lidar
                           ,ptv.beta_ext_radar
                           ,particle_parameters
                           ,ptv.zeta
                           ,ptv.phase
                           ,8.6e-3)
               
                

                #compute effective diameter from mode_diameter
                ptv.effective_diameter[0,0],ptv.mean_diameter[0,0],ptv.mean_mass_diameter[0,0]\
                    = spu.effective_diameter(
                      ptv.dmode
                     ,particle_parameters
                     ,ptv.zeta
                     ,ptv.phase)

                #compute liquid water content (kg/m^3)
                LWC_ext = spu.liquid_water_content(
                      ptv.effective_diameter
                     ,ptv.beta_ext_lidar,ptv.phase)
                if ptv.phase[0,0] == 1:
                    ptv.LWC_p180 = spu.liquid_water_content(
                        ptv.effective_diameter
                        ,ptv.beta_a_backscat/particle_parameters['p180_ice']
                        ,ptv.phase)
                elif ptv.phase[0,0] == 0 :
                    ptv.LWC_p180 = spu.liquid_water_content(
                        ptv.effective_diameter
                        ,ptv.beta_a_backscat/0.05
                        ,ptv.phase)
                else:
                    ptv.LWC_p180 = np.NaN
                    
                
                #compute precip rate (m/s)
                #precip_rate = 1/density (m^3/kg) * LWC (kg/m^3) * fall_velocity (m/s)
                ptv.hsrl_radar_precip_rate[0,0] = 0.001 * ptv.LWC_p180[0,0] * ptv.MeanDopplerVelocity[0,0] 
                #LR_precip_rate[0,0]= 0.001 * LWC_p180[0,0] * MeanDopplerVelocity[0,0 

               
                ptv.rw_fall_velocity[0,0],ptv.mw_fall_velocity[0,0],ptv.model_spectral_width[0,0]\
                     = spu.weighted_fall_velocity( 
                      ptv.dmode
                      ,particle_parameters
                      ,ptv.zeta
                      ,ptv.temps
                      ,ptv.pressures
                      ,ptv.phase)
     
          


            size_dist = D**alpha * np.exp(-(alpha/gam) * (D/ptv.dmode[0,0])**gam)
            norm      = np.sum(size_dist * dD)
            size_dist /=norm
            print 'norm',norm
            #area weighted distribution  
            aw_size_dist = D**(alpha+2.0) * np.exp(-(alpha/gam) * (D/ptv.dmode[0,0])**gam)
            norm = np.sum(aw_size_dist * dD)
            aw_size_dist /= norm
            print norm
            #mass weighted distribution
            mw_size_dist = D**(alpha+ ptv.zeta + 2.0) * np.exp(-(alpha/gam) * (D/ptv.dmode[0,0])**gam)
            norm = np.sum(mw_size_dist * dD)
            mw_size_dist /= norm
            print norm
            #radar weighted distribution
            rw_size_dist = D**(alpha+ 2.0 * ptv.zeta + 4.0) * np.exp(-(alpha/gam) * (D/ptv.dmode[0,0])**gam)
            norm = np.sum(rw_size_dist * dD)
            rw_size_dist /= norm
           
            #compute model fall velocity at this point
             #constants, program works in mks 
            rho_ice=0.91 * 1000.0 #kg/m^3
            rho_water = 1000.0 #kg/m^3
            gas_constant = 287 #J/(kg K)
            rho_air = 100.0 *ptv.pressures/(gas_constant * ptv.temps) #in kg/m^3
            eta = spu.dynamic_viscosity(ptv.temps)
            
            #X = spu.best_number(D,Dr,ptv.zeta,rho_ice,rho_air,eta)

            X_const = spu.best_constant(rho_air,eta)
            if ptv.phase == 1:
                X = spu.best_number(D,Dr,X_const,ptv.zeta,rho_ice)
            else:
                X = spu.best_number(D,Dr,X_const,1.0,rho_water)
                
            Vf = spu.spheroid_fall_velocity(D,X,ptv.zeta,rho_air,eta,ptv.phase)
            
            #cost_func[i,k] = (doppler_vel-rw_fall_vel)**2 + (spectral_width-model_spectral_width)**2
            #print 'cost function = ',cost_func[i,k]
            if alpha_gamma_str == 'map':
                 plt.figure(999)
                 imgplot = plt.imshow(np.log(cost_func))

            
            plt.figure(1000)
            #plt.ioff()
            lines=plt.plot(D*1e6  ,np.transpose(size_dist), colors[0]
                     ,D*1e6,np.transpose(aw_size_dist),colors[1]
                     ,D*1e6,np.transpose(mw_size_dist),colors[2]
                     ,D*1e6,np.transpose(rw_size_dist),colors[3])
            plt.grid(True)             
            ax=plt.gca()
            plt.setp(lines,linewidth=lwidth)
            plt.legend(['number','area','mass','radar'])
            ax.set_xlim(0,4500)
            plt.ylabel('weighted size distributions')
            plt.xlabel('Diameter (microns)')
        

            #plt.ion()
            plt.figure(1001)
            #plt.ioff()
            lines1=plt.plot(D*1e6  ,np.transpose(size_dist), colors[0]
                     ,D*1e6,np.transpose(aw_size_dist),colors[1]
                     ,D*1e6,np.transpose(mw_size_dist),colors[2]
                     ,D*1e6,np.transpose(rw_size_dist),colors[3])
            plt.grid(True)             
            ax=plt.gca()
            plt.setp(lines1,linewidth=lwidth)
            plt.legend(['number','area','mass','radar'])
            plt.xlabel('Diameter (microns)')
            ax.set_xscale('log')
            #plt.show(False)
            
            plt.figure(1002)
            #plt.ioff()
            lines2=plt.plot(D*1e6,np.transpose(Vf),'r')
            plt.grid(True)
            ax=plt.gca()
            plt.setp(lines2,linewidth=lwidth)
            plt.xlabel('diameter (microns)')
            plt.ylabel('fall velocity (m/s)')
            ax.set_xscale('log')
            
            
           

            print
            print 'effective_diameter_prime = %4.0f (microns)' %(1e6 * ptv.eff_diameter_prime) 
            print 'effective diameter       = %4.0f (micons)' %(ptv.effective_diameter * 1e6)
            print 'mode diameter            = %4.0f (microns)' %(ptv.dmode * 1e6)
            print 'mean diameter            = %4.0f (microns)' %(ptv.mean_diameter *1e6)
            print 'LWC , extinction         = %6.4f (gr/m^3)' %(ptv.LWC_ext * 1e3)
            print 'LWC ,p(180)/4pi = %4.2f   = %6.4f (gr/m^3)' \
                             %(particle_parameters['p180_ice'],ptv.LWC_p180 * 1e3)
            if ptv.phase[0,0] == 1:
                print 'phase                    = ice'
            elif ptv.phase[0,0] == 0:
                print 'phase                    = water'
            else:
                print 'phase                    = NaN'
            print 'Doppler velocity         = %4.2f (m/s)' %ptv.MeanDopplerVelocity
            print 'rw_fall_vel              = %4.2f (m/s)' %ptv.rw_fall_velocity
            print 'mw_fall_vel              = %4.2f (m/s)' %ptv.mw_fall_velocity
            print 'model_spectral_width     = %4.2f (m/s)' %ptv.model_spectral_width
            print 'radar_spectral_width     = %4.2f (m/s)' %ptv.radar_spectral_width
            print 'alpha = ' ,alpha
            print 'gamma = ',gam
            print 'zeta  = ',ptv.zeta[0,0]
        
           
            
        
            lwidth = lwidth + 1.0
            if not plt.isinteractive():   
               plt.show(False)
            #plt.ion()
           
            
                
        return

    
    #uncover the data structures
    struct = struct.rs.__dict__

    print 
    print "available data structures"
    print struct.keys()
    print

    particle_size_distribution(struct,t_index,z_index)
    return
