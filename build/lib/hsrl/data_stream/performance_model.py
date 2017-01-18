import matplotlib.pyplot as plt
import numpy as np
from matplotlib.dates import num2date, date2num
import os
import json
from lg_base.core.open_config import open_config

def performance_model(instrument,rs_Cxx,rs_cal,combined_counts,mol_counts,transmitted_energy,num_seeded_shots,times,constants):
    """performance_model(instrument,rs_Cxx,rs_cal,dc_combined_counts,dc_mol_counts\
                  ,transmitted_energy,num_seeded_shots,times,constants):
    compute theoretical combined hi return for the current data
    segment using system optical parameters in 'xxhsrl_performance_specs.json'.
    The Rayleigh scattering cross section is taken from rs_Cxx.
    The geometric correction is taken from rs_cal.
    Combined_counts are raw dark corrected counts, with no pileup
    or baseline corrections.  Comparisons are only valid in clear air with minimal
    low altitude extinction.
    """

    #extract calibration info for current data from rs_cal
    instrument=rs_cal.instrument
    geo_cor=rs_cal.geo.data
    print ' '
    print rs_cal.geo.header
    print ' '

   
    
    #convert the Rayleigh scattering cross section, beta_r, to lidar's 
    #native range resolution which is supplied in geo_cor[0,:]
    #this assumes that lidar is ground based
    ranges=rs_Cxx.msl_altitudes-constants['lidar_altitude']
    beta_r=np.zeros_like(geo_cor[:,0])
    Cmm = np.zeros_like(geo_cor[:,0])
    beta_r[:]=np.interp(geo_cor[:,0],ranges[:],rs_Cxx.beta_r[:])
    Cmm[:] = np.interp(geo_cor[:,0],ranges[:],rs_Cxx.Cmm[:])
   
    #fetch component specifications from json file
    dict_spec_file={'gvhsrl' :      'gvhsrl_specs.json',\
                    'ahsrl'  :      'ahsrl_specs.json',\
                    'nshsrl' :      'nshsrl_specs.json',\
                    'mf2hsrl':      'mf2hsrl_specs.json',\
                    'bagohsrl' :      'bagohsrl_specs.json',\
                    'rbhsrl' :        'rbhsrl_specs.json'}

    fd=open_config(dict_spec_file[instrument],verbose=True)
    print fd
    if 0:
     if os.path.isfile(dict_spec_file[instrument]):
       fd = open(dict_spec_file[instrument], 'r')
       print 'Preparing theoretical profiles using '\
             +dict_spec_file[instrument]
       print 'See file ' +dict_spec_file[instrument] +' for component specifications'
     else:
       print ' '
       print "Performance spec file: "+dict_spec_file[instrument]\
             +" not found"
       return
   
    dd = json.load(fd)
    specs=dd['specs']

    #these dictionaries contain the optical transmission of components in different
    #sections of the optical path
    general_specs     = specs['general']
    telescope_trans   = specs['telescope_trans']
    fore_optics_trans = specs['fore_optics_trans']
    sky_filter_trans = specs['sky_filter_trans']
    combined_trans   = specs['combined_trans']
    molecular_trans  = specs['molecular_trans']
    detector_eff     = specs['detector_eff']

    
   
    #constants
    plancks_constant  = 6.6e-34;        #planck's constant J s
    c                 = 3.0e8;          #speed of light m/s



    #calibrated transmitted energy per accumulation interval
    cal_transmitted_energy=transmitted_energy

    
    #get power as measured after telescope turning mirror
    #reported by transmitted energy monitior in current data 
    # if transmitted_energy not in specs dictionary
    print
    print str(constants['wavelength']) + 'nm wavelength used to compute scattering cross section'
    print
    if not general_specs.has_key('transmitted_pulse_energy'):
        #compute total energy transmitted in requested time interval
        total_transmitted_energy=np.nansum(cal_transmitted_energy)
        print 'total_energy',total_transmitted_energy
        #note addition of one accumulation interval for first accumulation
        #delta_t=3600*24*(times[-1]-times[0]+times[1]-times[0])
        delta_t =  (times[-1]-times[0]+times[1]-times[0])
        delta_t = delta_t.total_seconds()
        print 'delta_t= ',delta_t
        #compute average transmitted power
        transmitted_pwr= total_transmitted_energy/delta_t

        #compute transmitted energy per laser pulse 
        pulse_energy=total_transmitted_energy/np.nansum(num_seeded_shots)
        print ' '
        print 'Power and shot energy reported by tranmitted energy monitor'
    else:
        #pulse_energy is in mJ
        pulse_energy = general_specs['transmitted_pulse_energy']*1000.0
        transmitted_pwr = pulse_energy * general_specs['repetition_rate'] 
        print ' '
        print 'Power and shot energy as reported by specs.json file'
    print 'transmitted power=           ',transmitted_pwr, ' mW'
    print 'transmitted energy per shot= ',pulse_energy, ' mJ'

   
    
    r_primary=float(general_specs['primary_mirror_dia'])/2.0
    r_obstruction=float(general_specs['secondary_obstruction_dia'])/2.0
    r_e_width=float(general_specs['e_output_beam_dia'])/2.0
    #part of the outgoing energy is blocked by secondary
    #area under Gaussian with distribution of output power on primary
    weighted_beam_area=1-np.exp(-(r_primary/r_e_width)**2)
    weighted_obs_area=1-np.exp(-(r_obstruction/r_e_width)**2)
    pulse_energy=pulse_energy*(weighted_beam_area-weighted_obs_area)\
         /weighted_beam_area

    #transmission of telescope and output window.
    telescope = compute_trans('telescope and output windows ',telescope_trans)
    
    
    max_power_density = telescope\
         *transmitted_pwr/(np.pi*r_e_width**2)\
         *np.exp(-(r_obstruction/r_e_width)**2)
    safe_limit=5e-7*1000   #single pulse mJ/cm^2
    pwr_limit=safe_limit*general_specs['repetition_rate']\
               /(0.5*general_specs['repetition_rate'])**0.25 #multi pulse
    print pwr_limit
    print 'Power density at obstruction edge                 =     ',max_power_density/1e4, 'mW/cm^2'
    print 'Power density divided by eye safe limit           =     ',max_power_density/(1e4*pwr_limit)
    print 'Energy transmitted after secondary blockage       =     ',pulse_energy, ' mJ'
    pulse_energy = pulse_energy*telescope;
    print 'Output optics transmission                        =     ',telescope
    print 'Energy transmitted after window and mirror losses =     ',pulse_energy, ' mJ'
    print ' '

    receiver_area=(np.pi/4.0)*(general_specs['primary_mirror_dia']**2\
              -general_specs['secondary_obstruction_dia']**2)


    #compute counts/laser pulse incident on reciever area with no window losses
    #note conversion of energy in mJ to Joules and correction of beta_r for conversion
    #from reference system wavelength to modeled system wavelength


    incident_power=np.zeros_like(beta_r)
    incident_counts=np.zeros_like(beta_r)

    #remove any NaN's in beta_r
    beta_r[np.isnan(beta_r)]=0
    
    incident_power[:]=pulse_energy/1000.0\
          *c/2\
          *receiver_area\
          *3/(8.0*np.pi)\
          *1/(geo_cor[:,0]**2)\
          *beta_r[:]\
          *np.exp(-np.cumsum(beta_r[:])*c*constants['binwidth'])
   
    
    #convert power to counts 
    incident_counts[:]=incident_power[:]\
          *constants['binwidth']\
          *general_specs['transmit_wavelength']/(plancks_constant*c)
    print
    print 'incident counts computed in binwidth = ', constants['binwidth']
    print 'beta_r[0]',beta_r[0]
    print
  
    #remove 1/r^2 component from geo correction
    geo_cor[:,1]=1e6*geo_cor[:,1]/geo_cor[:,0]**2
    #geo cor can't amplify signal
    geo_norm=np.min(geo_cor[100:1000,1])
    #rescale to geo_cor to minimum value
    geo_cor[:,1]=geo_cor[:,1]/geo_norm
    
    
    incident_counts[~np.isfinite(incident_counts)]=0
    incident_geo_cor_counts=np.zeros_like(beta_r)
    incident_geo_cor_counts[:]=incident_counts[:]/geo_cor[:,1]
    incident_geo_cor_counts[~np.isfinite(incident_geo_cor_counts)]=0

   


      

    fore_optics = compute_trans('optics between telescope and field stop',fore_optics_trans)
    
    sky_filter = compute_trans('skylight filter',sky_filter_trans)

    combined_channel = compute_trans('etalon to combined detector',combined_trans)               

    molecular_channel = compute_trans('combined/mol beamsplitter to mol detector',molecular_trans)      

   
    #combined channel all components
    detector_eff = specs['detector_eff']
    combined_eff=telescope\
          *fore_optics\
          *sky_filter\
          *combined_channel\
          *detector_eff['QE_combined']
    print 
    print 'combined channel efficiency                  = ', combined_eff 

    #molecular_eff = combined_eff\
          
    #      *molecular_channel\
    #      *specs['QE_combined']
    #print 'combined channel efficiency                  = ', combined_eff 
    


    start_time_str=times[0].strftime("%d-%b-%y %H:%M")
    end_time_str=times[-1].strftime("%d-%b-%y %H:%M")
    fig = 2998
    plt.figure(fig)
    plt.clf()
    plt.figure(fig).canvas.set_window_title('Fig '+str(fig)+'        Geometric correction')
    #plt.rcParams['figure.figsize'] = display_defaults['profile_graph_size']
    plt.plot(geo_cor[:,1],geo_cor[:,0],'b')
    # ,interp_geo_cor,msl_altitudes[:],'r')
    plt.grid(True)
    plt.xlabel('Geometric correction')
    plt.ylabel('Range bin number') 
    ax=plt.gca()
    ax.set_xscale('log')
    plt.title(instrument+'   geo correction '+start_time_str+'-->'+end_time_str)
    


    laser_pulse_bin=np.int(constants['apd_pulse_timing'][1]/constants['binwidth'])



    #do a pileup correction on mean profile
    [nshots,nalts]=combined_counts.shape;
    ones_array=np.ones((nshots,nalts))
    ones_array_t=np.transpose(ones_array)
    p_corr=constants['combined_hi_dead_time']\
         *combined_counts/constants['binwidth']
    p_corr[p_corr>.99]=.95
    pileup_cor_hi_counts=combined_counts/(ones_array-p_corr)

    #plot results
    nbins=incident_counts.shape[0]
    bin_vec=np.arange(nbins)
   
    start_time_str=times[0].strftime("%d-%b-%y %H:%M")
    end_time_str=times[-1].strftime("%d-%b-%y %H:%M")
    #incident counts per laser pulse
    fig = 2999
    plt.figure(fig)
    plt.clf()
    plt.figure(fig).canvas.set_window_title('Fig '+str(fig)+'        modeled vs measured')
    #plt.rcParams['figure.figsize'] = display_defaults['profile_graph_size']
    print
    print
    print "instrument =" ,instrument
    alt_scale = constants['binwidth'] * 1.5e5
    if not instrument == 'rbhsrl':
        print 'not rbhsrl, instrument =',instrument
        lines =plt.plot(incident_counts[:],bin_vec,'b'\
           ,incident_geo_cor_counts[:],bin_vec,'g'\
           ,incident_geo_cor_counts[:]*telescope,bin_vec,'c'\
           ,incident_geo_cor_counts[:]*telescope*fore_optics,bin_vec,'m'\
           ,incident_geo_cor_counts[:]*telescope*fore_optics\
             *sky_filter,bin_vec,'k'\
           ,incident_geo_cor_counts[:]*telescope*fore_optics\
             *sky_filter*combined_channel,bin_vec,'k'\
           ,incident_geo_cor_counts[:]*combined_eff,bin_vec,'r'\
           ,incident_geo_cor_counts[:]*Cmm[:]*combined_eff*detector_eff['QE_molecular'] \
                        /detector_eff['QE_combined'],bin_vec,'b')             
         
        plt.setp(lines[6],linewidth=3)
        plt.setp(lines[7],linewidth=3)
        plt.legend(('incident','in FOV','after tel','field stop','after etalon'\
                ,'at detector','comb_model','mol_model'),'upper right')
        plt.grid(True)
        plt.xlabel('Combined high (counts/shot/bin)')
        ax=plt.gca()
        ax.set_xscale('log')
    
        plt.title(instrument+' modeled at source binwidth '+start_time_str+'-->'+end_time_str)
    else:

        #plot at input systems binwidth 
        lines =plt.plot(incident_counts[:],bin_vec,'b'\
              ,incident_geo_cor_counts[:],bin_vec,'g'\
              ,incident_geo_cor_counts[:]*telescope,bin_vec,'c'\
              ,incident_geo_cor_counts[:]*telescope*fore_optics,bin_vec,'m'\
              ,incident_geo_cor_counts[:]*telescope*fore_optics\
               *sky_filter,bin_vec,'k'\
              ,incident_geo_cor_counts[:]*telescope*fore_optics\
              *sky_filter*combined_channel,bin_vec,'k'\
              ,incident_geo_cor_counts[:]*combined_eff,bin_vec,'r'\
              ,incident_geo_cor_counts[:]*Cmm[:]*combined_eff*detector_eff['QE_molecular'] \
                        /detector_eff['QE_combined'],bin_vec,'b')
        plt.setp(lines[6],linewidth=3)
        plt.setp(lines[7],linewidth=3)
        plt.legend(('incident','in FOV','after tel','after fore optics','field stop','after etalon'\
                ,'comb_model','mol_model'),'upper right')
        plt.axis([1e-5, 1e2,0,2000])
        plt.ylabel('Range bin')
        plt.grid(True)
        ax=plt.gca()
        ax.set_xscale('log')
        
        if 1: #plot at binwidth of modeled system
           plt.figure(3000) 
           binwidth_scale = general_specs['binwidth'] / constants['binwidth']  #scaling for larger bin size
          
           lines =plt.plot(binwidth_scale * incident_counts[:],bin_vec * alt_scale,'b'\
              ,binwidth_scale * incident_geo_cor_counts[:],bin_vec * alt_scale,'g'\
              ,binwidth_scale * incident_geo_cor_counts[:]*telescope,bin_vec * alt_scale,'c'\
              ,binwidth_scale * incident_geo_cor_counts[:]*telescope*fore_optics,bin_vec * alt_scale,'m'\
              ,binwidth_scale * incident_geo_cor_counts[:]*telescope*fore_optics\
              *sky_filter,bin_vec * alt_scale,'k'\
              ,binwidth_scale * incident_geo_cor_counts[:]*telescope*fore_optics\
              *sky_filter*combined_channel,bin_vec * alt_scale,'k'\
              ,binwidth_scale * incident_geo_cor_counts[:]*Cmm[:]*combined_eff*detector_eff['QE_molecular'] \
                        /detector_eff['QE_combined'],bin_vec * alt_scale,'r')
           plt.setp(lines[6],linewidth=3)
           plt.legend(('incident','in FOV','after tel','after fore optics','field stop','after etalon'\
                ,'mol_model'),'upper right')
           plt.axis([1e-3, 1e2,0,15]) 
           plt.ylabel('Altitude (km)')
           plt.grid(True)
           ax=plt.gca()
           ax.set_xscale('log')
           integration_time = general_specs
           plt.xlabel('Counts/pulse/'+str(general_specs['binwidth'])+'s binwidth')


    fig = 3001
    plt.figure(fig)
    plt.clf()
    plt.figure(fig).canvas.set_window_title('Fig '+str(fig)+'        modeled counts/sec')
    #plt.rcParams['figure.figsize'] = display_defaults['profile_graph_size']
   

    #convert to binwidth of modeled system
    binwidth_scale = general_specs['binwidth'] / constants['binwidth']
    
    if 0:   #normal plot with combined
        lines =plt.plot(
           incident_geo_cor_counts[:]*combined_eff * general_specs['repetition_rate']\
           *binwidth_scale,bin_vec/binwidth_scale,'r'\
           ,incident_geo_cor_counts[:]*Cmm[:]*combined_eff *binwidth_scale\
           *general_specs['repetition_rate']\
           *detector_eff['QE_molecular']/detector_eff['QE_combined'],bin_vec ,'b')
        plt.axis([1e-1, 1e4,0,bin_vec[-1]/binwidth_scale])
        plt.ylabel('Range bin number')
        plt.xlabel('Count rate (counts/bin/sec)  ,  bin = '+str(general_specs['binwidth']*1e6)+' us')
        plt.grid(True)
        ax=plt.gca()
        ax.set_xscale('log')
    else:
         detector_noise = constants['mol_detector_dark_count']
         shots_summed = general_specs['integration_time_seconds']\
                         * general_specs['fraction_of_pulses_in_mol']\
                         * general_specs['repetition_rate']
         detector_noise = np.sqrt(constants['mol_detector_dark_count'] *shots_summed * binwidth_scale)
         print 'number of mol shots summed in integration period = ',shots_summed
         if general_specs.has_key('bright_sky_noise'):
             sky_noise=np.sqrt(general_specs['bright_sky_noise'] * shots_summed * binwidth_scale)
             signal_to_noise_ratio =(incident_geo_cor_counts[:] *Cmm[:]*combined_eff*binwidth_scale*np.sqrt(shots_summed))\
                        /np.sqrt(general_specs['bright_sky_noise']
                        + constants['mol_detector_dark_count']
                        + incident_geo_cor_counts[:] *Cmm[:]*combined_eff*binwidth_scale)
             lines =plt.plot([detector_noise,detector_noise],[0,bin_vec[-1]*alt_scale],':k'
               ,shots_summed * incident_geo_cor_counts[:]*Cmm[:]*combined_eff *binwidth_scale\
               *detector_eff['QE_molecular']/detector_eff['QE_combined'],bin_vec * alt_scale,'r'
               ,[sky_noise,sky_noise],[0,bin_vec[-1]*alt_scale],':g'
               ,signal_to_noise_ratio,bin_vec *alt_scale,'c')

         else:    
             lines =plt.plot([detector_noise,detector_noise],[0,bin_vec[-1]*alt_scale],'-c'
               ,shots_summed * incident_geo_cor_counts[:]*Cmm[:]*combined_eff *binwidth_scale\
               *detector_eff['QE_molecular']/detector_eff['QE_combined'],bin_vec * alt_scale,':r')
         plt.axis([1e-1, 1e4,0,bin_vec[-1] * alt_scale])
         plt.ylabel('Altitude (km)')
         plt.setp(lines[:],linewidth=3)
         plt.grid(True)
         plt.xlabel('Counts/bin (bin = '+str(int(general_specs['binwidth']*1e9)) + 'ns, dt= ' +str(general_specs['integration_time_seconds'])+'s)')
         ax=plt.gca()
         ax.set_xscale('log')
         plt.title(instrument+' modeled count rate in model bins '+start_time_str+'-->'+end_time_str)
         plt.axis([1, 2e5,0,15])
   
    print ' '
    print 'incident signal at pt 2000 = ', incident_geo_cor_counts[2000] 
    print 'modeled signal at pt 2000  = ', incident_geo_cor_counts[2000]*combined_eff
    plt.ion()
    plt.show()








def compute_trans(optics_section_name,components):
    """multiply transmission of individual components in optics session
       to calculate total through put"""
    trans =1.0
    print 
    print optics_section_name
    print '       component          throughput'
    for item in components:
        component_trans = components[item]
        print '      ', item,' = ',component_trans
        trans = trans * component_trans
    print '       section thoughput = ' , trans    
    return trans



"""
 #compute window and telescope transmission for received photons
    if instrument.find('ahsrl')>=0:
       telescope_trans=specs['shelter_window_trans']\
            *specs['telescope_window_trans']\
            *specs['primary_mirror_refl']\
            *specs['secondary_mirror_refl']
       print 'telsecope efficiency                          =  ',telescope_trans

       #optics between telescope and field stop
       fore_optics_trans=specs['turning_mirror_refl']\
         *specs['Q_wave_plate_trans']\
         *specs['cpol_pickup_trans']\
         *specs['turning_mirror_refl']\
         *specs['transmit_receive_pol_refl']\
         *specs['turning_mirror_refl']\
         *specs['half_wave_plate_trans']\
         *specs['cpol_pol_combiner_trans']\
         *specs['turning_mirror_refl']\
         *specs['half_wave_plate_trans']\
         *specs['combined_hi_field_lens_trans']
       print 'optics between telescope and field stop trans = ',fore_optics_trans

       #optics between field stop and output of etalon
       sky_filter_trans=specs['field_stop_lens_trans']\
           *specs['interf_filter_trans']\
           *specs['turning_mirror_refl']\
           *specs['etalon_trans']
       print 'skylight filter trans                          =  ',sky_filter_trans   
      
       #optics between etalon and combined detector
       combined_trans=specs['cpol_beamsplitter_refl']\
         *specs['half_wave_plate_trans']\
         *specs['combined_hi_mol_beamsplitter_refl']\
         *specs['half_wave_plate_trans']\
         *specs['combined_lo_bs_trans']\
         *specs['combined_hi_field_lens_trans']
       print 'combined channel efficiency                  = ', combined_trans
    
    elif instrument.find('gvhsrl')>=0:

       telescope_trans=specs['shelter_window_trans']\
            *specs['telescope_window_trans']\
            *specs['primary_mirror_refl']\
            *specs['secondary_mirror_refl']
       print 'telsecope efficiency                          =  ',telescope_trans
       
       #optics between telescope and field stop
       fore_optics_trans=specs['turning_mirror_refl']\
         *specs['Q_wave_plate_trans']\
         *specs['turning_mirror_refl']\
         *specs['turning_mirror_refl']\
         *specs['transmit_receive_pol_refl']\
         *specs['turning_mirror_refl']\
         *specs['half_wave_plate_trans']\
         *specs['turning_mirror_refl']\
         *specs['cpol_pol_combiner_refl']\
         *specs['turning_mirror_refl']\
         *specs['field_stop_lens_trans']
       print 'optics between telescope and field stop trans = ',fore_optics_trans

   

       #optics between field stop and output of etalon
       sky_filter_trans=specs['field_stop_lens_trans']\
           *specs['interf_filter_trans']\
           *specs['etalon_trans']    
       print 'skylight filter trans                          =  ',sky_filter_trans             
    
       telescope_trans=specs['shelter_window_trans']\
            *specs['telescope_window_trans']\
            *specs['primary_mirror_refl']\
            *specs['secondary_mirror_refl']\
            *specs['turning_mirror_refl'] 
       print 'telsecope efficiency                          =  ',telescope_trans
    
       combined_trans=specs['cpol_beamsplitter_trans']\
         *specs['turning_mirror_refl']\
         *specs['combined_hi_mol_beamsplitter_trans']\
         *specs['combined_lo_bs_trans']\
         *specs['turning_mirror_refl']\
         *specs['combined_hi_field_lens_trans']
        
       print 'combined_channel optics trans                 = ', combined_trans

       #combined QE
       combined_eff=telescope_trans\
         *fore_optics_trans\
         *sky_filter_trans\
         *combined_trans\
         *specs['QE_combined']
       print 'combined channel efficiency                  = ', combined_eff 


    #elif instrument.find('mf2hsrl')>=0:"""
