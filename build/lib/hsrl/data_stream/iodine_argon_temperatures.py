import os
import numpy as np
from bottleneck import nanmean,nansum
from copy import copy
from lg_base.core.fmt_mpl_date import fmt_mpl_datetime, fmt_mpl_time
from datetime import datetime
import lg_base.core.array_utils as hau
import hsrl.data_stream.processing_utilities as pu
import hsrl.calibration.calibration_utilities as cu
from scipy import interpolate
import hsrl.filters.savitzky_golay as sg

def compute_temperature_from_i2a(instrument,i2a_mol_ratio,i2a_temp_table,sounding_pressures,corr_adjusts):
    if not instrument == 'bagohsrl':
        print 'no temperature measurements from ', instrument
        return
    if i2a_temp_table.header == None:
        print
        print 'no i2a_temp_table found in month dir'
        print
        return
    else:
       print 'i2a_temp_table_header = ',i2a_temp_table.header
    table_pressures = i2a_temp_table.data[0,1:]
    table_i2a_i2_ratios = i2a_temp_table.data[1:,0]
    table_temps = i2a_temp_table.data[1:,1:]
    ratios=i2a_mol_ratio[0,:] * corr_adjusts['i2a_corr']/corr_adjusts['i2_corr']
    #ratios[np.isnan(ratios)]=1.0
    pressures = sounding_pressures[:len(i2a_mol_ratio[0,:])]
    
   
    i2a_temperatures=np.nan*np.zeros_like(pressures)

    spl = interpolate.RectBivariateSpline(table_i2a_i2_ratios,table_pressures,table_temps)

    for i in range(len(pressures)):
        i2a_temperatures[i]= spl(ratios[i],pressures[i])

        #print ratios[i],pressures[i],i2a_temperatures[i]

    return i2a_temperatures    



def make_i2a_mol_ratio_table(wavelength,i2_scan,process_defaults):
    table_temperatures = np.arange(200.0,300.0,5.0)        
    table_pressures = np.arange(1020.0,100.0,-5.0)
    i2a_i2_table = np.zeros((len(table_temperatures),len(table_pressures)))

    if (process_defaults.get_value('molecular_spectrum','model') == 'tenti_s6'):
        from tenti_s6 import tenti_s6
        for i in range(len(table_temperatures)):
            for j in range(len(table_pressures)):
                spectrum = tenti_s6(wavelength * 1e-9,table_temperatures[i],
                    table_pressures[j],
                    i2_scan[:, 0] * 1e9)
                #ratio of i2a to i2 signals as function of temp and press
                i2a_i2_table[i,j] = (np.sum(spectrum * i2_scan[:,6])) \
                         /(np.sum(spectrum *i2_scan[:,2]))
    elif   process_defaults.get_value('molecular_spectrum','model') == 'witschas':
        for i in range(len(table_temperatures)):
            for j in range(len(table_pressures)):
                spectrum = cu.witschas_spectrum(table_temperatures[i]
                    ,table_pressures[j],wavelength*1e-9
                    ,i2_scan[:, 0] * 1e9)
                spectrum = spectrum/np.sum(spectrum)
                #ratio of i2a to i2 signals as function of temp and press
                i2a_i2_table[i,j] = np.sum(spectrum * i2_scan[:,2]/i2_scan[:,6])
    elif   process_defaults.get_value('molecular_spectrum','model') == 'maxwellian':
        # spectral width of molecular scattering
        m_bar = 28.97 * 1.65978e-27  # average mass of an air molecule
        sigma_0 = 1 / (wavelength * 1e-9)  # number in 1/meters
        kb = 1.38044e-23  # Boltzmans constant J/(K deg)
        c = 3e8  # speed of light in m/s
        sigma = i2_scan[:, 0] * 1e9 / c  # wavenumber vector
        
        for i in range(len(table_temperatures)):
            for j in range(len(table_pressures)):
                norm = m_bar * c ** 2 / (8 * sigma_0 ** 2 * kb
                    * sounding.temps[ i])
                spectrum = exp(-norm * sigma ** 2) 
                spectrum = spectrum/np.sum(spectrum)
                #ratio of i2a to i2 signals as function of temp and press
                i2a_i2_table[i,j] = np.sum(spectrum * i2_scan[:,2]/i2_scan[:,6])        
    import matplotlib.pyplot as plt
    #print table_pressures[10]
    #print i2a_i2_table.shape
    plt.figure(199)
    plt.plot(i2_scan[:,0],i2_scan[:,2],'b',i2_scan[:,0],i2_scan[:,6],'r')
    ax=plt.gca()
    ax.grid(True)

    plt.figure(200)

    for i in range(0,180,10):
        plt.plot(table_temperatures,i2a_i2_table[:,i],'b')
    #plt.plot(table_temperatures,i2a_i2_table[:,0],'b'
    #        ,table_temperatures,i2a_i2_table[:,10],'r'
    #        ,table_temperatures,i2a_i2_table[:,180],'g')
    ax=plt.gca()
    ax.grid(True) 
    plt.show()

    #convert table to temperatures vs i2a_ratio and pressure

    ratio_min=np.min(i2a_i2_table)
    ratio_max=np.max(i2a_i2_table)
    ratios = np.arange(ratio_min,ratio_max,(ratio_max-ratio_min)/100)
    temperatures = np.Nan*np.zeros((len(ratios),len(table_pressures)))
    for i in range(len(table_pressures)):
        for j in range(len(ratios)):
             temperature[i,j] = np.interp(ratios(i),i2a_i2_table[:,i],table_temps)      
    
    return i2a_i2_table,table_temperatures,table_pressures 






def compute_i2a_mol_ratio(profiles,Cxx,corr_adjusts,process_defaults):
    """compute_i2a_mol_ratio(i2a,mol,process_defaults):
       compute ratio of argon broadened i2 filtered profile to normal
       i2 filtered profile
       i2a = corrected argon broadened profile
       mol = normal i2 filtered profile
       process_defaults = processing instructions from process_control.json"""


    i2a = profiles.molecular_i2a_counts[0,:].copy()
    i2a = i2a * corr_adjusts['i2a_ratio']
    mol = profiles.molecular_counts[0,:].copy()
    
    
    if process_defaults.enabled('i2a_mol_ratio'):
        #compute filter length in bins
        window = np.int(np.float(process_defaults.get_value('i2a_mol_ratio','filter_window')\
                 /(profiles.msl_altitudes[2]-profiles.msl_altitudes[1])))
        print 'applying 3-order savitzky_golay filter before computing i2a/i2 ratio',
        print ' window = ',window * (
              profiles.msl_altitudes[2]-profiles.msl_altitudes[1]), ' meters'
        #window must be odd number
        if window == 2*(window/2):
            window = window +1
    
        if window >=5:
            i2a=sg.savitzky_golay(i2a,window,3)
            mol=sg.savitzky_golay(mol,window,3)
        else:
            print 'Window too short---no filtering applied to i2a_mol_ratio'
   
    i2a_mol_ratio = (i2a 
         - Cxx.Cam_i2a*corr_adjusts['Cam_corr']
         * profiles.combined_hi_counts[0,:])\
         /(mol - Cxx.Cam*corr_adjusts['Cam_corr'] 
         * profiles.combined_hi_counts[0,:])
    i2a_mol_ratio = hau.Z_Array(i2a_mol_ratio[np.newaxis,:])
   
    return i2a_mol_ratio

    
def compute_i2a_diff_geo_corr(instrument,rs,rs_cal,process_defaults,constants,corr_adjusts):
    """compute_i2a_diff_geo_corr(i2a_molecular_counts,molecular_counts,time,rs_cal
       ,processing defaults,constants)
       from data acquired with the i2a filter removed.
    """
   
    import matplotlib.pyplot as plt
    import hsrl.data_stream.hsrl_read_utilities as hru


    #find bin # at range = 0
    [dark_end_time,pulse_time,end_cal_time]=constants['apd_pulse_timing']
 
    native_res=constants['binwidth']*3e8/2

    lidar_bin = np.int((pulse_time *3e8)/native_res)+1
    
    #restore all channels to uncorrected signals
    rs.molecular_counts   = rs.raw_molecular_counts.copy()
    rs.combined_hi_counts = rs.raw_combined_hi_counts.copy()
    rs.cross_pol_counts   = rs.raw_cross_pol_counts.copy()
    if hasattr(rs,'molecular_i2a_counts'):
        rs.molecular_i2a_counts = rs.raw_molecular_i2a_counts.copy()
    if hasattr(rs,'combined_lo_counts'):
        rs.combined_lo_counts = rs.raw_combined_lo_counts.copy()

    nalts = rs.molecular_counts.shape[1]
    
    if hasattr(rs,'raw_molecular_counts'):
        import matplotlib.pylab as plt
        bins=np.arange(rs.raw_molecular_counts.shape[1])
        plt.figure(5000)
        plt.plot(nanmean(rs.raw_molecular_counts,0),bins,'b',nanmean(rs.raw_molecular_i2a_counts,0),bins,'c')
        ax=plt.gca()
        ax.set_xscale('log')
    
    #redo dark and baseline corrections
    rs = pu.dark_count_correction(instrument,rs,Cxx,corr_adjusts,constants)
    rs = pu.baseline_correction(rs,rs_cal,nalts,corr_adjusts,constants)

    plt.figure(5001)
    plt.plot(nanmean(rs.molecular_counts,0),bins,'b',nanmean(rs.molecular_i2a_counts,0),bins,'c')
    ax=plt.gca()
    ax.grid(True)
    ax.set_xscale('log')

    
    norm_range = np.float(raw_input('Normalization range (m) ?'))
    
    norm_bin   = norm_range/(3e8*constants['binwidth'])
    norm_bin = np.int(norm_bin +constants['apd_pulse_timing'][1]/constants['binwidth'])+1
   
    i2a=nanmean(rs.molecular_i2a_counts,0)
    mol=nanmean(rs.molecular_counts,0)
    i2a_raw = i2a.copy()
    mol_raw = mol.copy()

    plt.figure(5002)
    plt.plot((i2a_raw/mol_raw)*mol_raw[norm_bin]/i2a_raw[norm_bin],bins,'c')
    ax=plt.gca()
    ax.set_xlim((.90,1.05)) 
    ax.grid(True)

    
    #apply third order, 21 bin filter (150 m filter)
    window=21
    i2a=sg.savitzky_golay(i2a,window,3)
    mol=sg.savitzky_golay(mol,window,3)

    plt.figure(5003)
    plt.plot((i2a/mol)*mol[norm_bin]/i2a[norm_bin],bins,'c')
    ax=plt.gca()
    ax.set_xlim((.90,1.05)) 
    ax.grid(True)



    i2a[:70] = i2a_raw[:70]
    mol[:70] = mol_raw[:70]

    plt.figure(5004)
    plt.plot((i2a/mol)*mol[norm_bin]/i2a[norm_bin],bins,'c')
    ax=plt.gca()
    ax.set_xlim((.90,1.05)) 
    ax.grid(True)

    i2a_2=i2a.copy()
    mol_2=mol.copy()
    
    #apply third order, 81 bin filter (600 m filter)
    window=81
    i2a=sg.savitzky_golay(i2a,window,3)
    mol=sg.savitzky_golay(mol,window,3)

    #leave bottom points unfiltered
    i2a[:100]=i2a_2[:100]
    mol[:100]=mol_2[:100]
    
    bin_vec = np.arange(mol.shape[0])
    #print bin_vec.shape,i2a.shape,mol.shape
    plt.figure(900)
    plt.plot(i2a,bin_vec,'b',mol,bin_vec,'r')
    ax=plt.gca()
    plt.xlabel('Counts')
    plt.ylabel('bin number')
    ax.legend(['mol_i2a','mol'],'upper right')
    ax.set_xscale('log')
    ax.grid(True)
    
    i2a_mol_ratio=i2a/mol
    i2a_mol_ratio=i2a_mol_ratio/i2a_mol_ratio[norm_bin]


    plt.figure(1000)
    plt.plot(i2a_raw/mol_raw/i2a_mol_ratio[norm_bin],'c',i2a_mol_ratio,bin_vec,'r')
    ax=plt.gca()
    ax.set_xlim((.90,1.05)) 
    ax.grid(True)

    raw_i2a_mol_ratio = i2a_mol_ratio.copy()

    #801 bin 3rd order filter
    window =801
    i2a_mol_ratio = sg.savitzky_golay(i2a_mol_ratio,window,3)
    i2a_mol_ratio[:600]=raw_i2a_mol_ratio[:600]
    
    plt.figure(1001)
    plt.plot(i2a_mol_ratio,bin_vec,'r')
    ax=plt.gca()
    ax.set_xlim((.90,1.05)) 
    ax.grid(True)
    
    i2a_mol_ratio[:lidar_bin]=np.NaN    
    inx = np.int(raw_input('make constant above bin #?'))    
    i2a_diff_geo = i2a_mol_ratio.copy()
    i2a_diff_geo[inx:]=i2a_diff_geo[inx]

    plt.plot(i2a_mol_ratio,bin_vec,'c',i2a_diff_geo,bin_vec,'r')
    ax=plt.gca()
    ax.set_xlim((.90,1.05))
    ax.grid(True)
    plt.ylabel('bin number')
    plt.xlabel('i2a_dif_geo')

    #write i2a_mol_diff_geofile
    #get path to directory
    import hsrl.calibration.cal_file_generation as cfg
    dir_path=cfg.calibration_path_for(instrument,rs.times[0],process_defaults)
    #define start time for use in filename
    start_time_str=rs.times[0].strftime("%Y%m%dT%H%M")          
    filename=os.path.join(dir_path,'i2a_mol_diff_geofile_'+start_time_str+'.geo')
    fileid=open(filename,'w')
    os.chmod(filename,0664)
    print >>fileid, '#profile data from %s -->%s UTC' %(rs.times[0],rs.times[-1])
    print >>fileid, '# Gains normalized range bin = %i ' %(norm_bin)
    print >>fileid, '#bin i2a_mol/mol'
    for i in range(len(mol)):
        print >>fileid, '%i   %10.4e'\
        %(i ,i2a_diff_geo[i])
    fileid.close()
    
    print '\nnew i2a_mol_diff_geofile=\n '+filename

