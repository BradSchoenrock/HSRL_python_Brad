

import numpy as np
import os
import scipy.signal as scipysig
import matplotlib.pyplot as plt
import copy
from datetime import datetime, timedelta
import scipy.optimize as sco
from scipy.interpolate import UnivariateSpline
import lg_base.core.array_utils as hau
import lidar.lidar_utilities as lu
import hsrl.calibration.calibration_utilities as cu
from lg_base.core.locate_file import locate_file
from hsrl.filters.polynomial_smoothing import polynomial_smoothing
# Try to use the much faster nanmean from bottleneck, otherwise fall back
# to the scipy.stats version
try:
    from bottleneck import nanmean, nansum
except ImportError:
    print
    print "No bottleneck.nanmean available! Falling back to SLOW scipy.stats.nanmean"
    print
    from scipy.stats import nanmean, nansum
import hsrl.data_stream.processing_utilities as pu
import atmospheric_profiles.soundings.sounding_utilities as su
import hsrl.filters.savitzky_golay as sg
import lg_base.core.read_utilities as ru
 
def calibration_path_for(instrument,date,process_defaults,alternate_cal_dir=None):
    import lg_dpl_toolbox.core.archival as hru
    basedir=hru.get_path_to_data(instrument,date)
    if alternate_cal_dir==None:
        if process_defaults!=None:
            alternate_cal_dir = process_defaults.get_value('alternate_cal_dir','full_dir_path') 
            if alternate_cal_dir == 'None':
                alternate_cal_dir = None
        if alternate_cal_dir==None:
            alternate_cal_dir=os.getenv('HSRL_ALTERNATE_CAL_DIR',None)
    if alternate_cal_dir!=None:
        basedir=os.path.join(alternate_cal_dir,instrument)
    ret=os.path.join(basedir,'%04i' % date.year,'%02i' % date.month,'%02i' % date.day,"calibration")
    try:
      os.makedirs(ret)
      os.chmod(ret,0775)
    except:
      pass
    return ret

class Ray_fit:
    """Fit signal with Rayleigh return"""
    
    def __init__(self,bin_offset,beta_r,bin_vec):
           self.bin_offset=bin_offset
           self.beta_r=beta_r
           
           
    def sig(self,k,x):
           return k*self.beta_r[x+self.bin_offset]/(x**2)
       
def make_i2a_diff_geofile(instrument,raw,rs_Cxx,rs_cal,process_defaults,constants,corr_adjusts):
    """compute_i2a_diff_geo_corr(instrument,raw,rs_cal,processing defaults,constants)
       1) Takes pileup corrected count profiles from data acquired with the i2a filter removed.
       2) Repeats dark and baseline corrections.
       3) Computes the ratio of corrected i2a and i2 counts from time averaged profiles.
    """
  
    if not hasattr(raw,'molecular_i2a_counts'):
        print
        print 'WARNING:*****no molecular_i2a_counts in record****************'
        print '             can not compute i2a_diff_geo_corr'
        return
    
    import matplotlib.pyplot as plt
   
    #make vector containing altitudes of the range bins
    nbins = len(raw.molecular_counts[0,:])
    bin_vec = np.arange(nbins)
    binwidth = constants['binwidth'] * 1.5e8  #binwidth in meters
    sbin = np.int(constants['apd_pulse_timing'][1]/constants['binwidth'])
    lidar_altitude = constants['lidar_altitude']
    bin_ranges = binwidth * np.arange(nbins)

    rs = hau.Time_Z_Group()

    #introduce i2a channel clock skew
              
                                                
    # generate dark and baseline corrected data.
    rs = pu.dark_count_correction(instrument,raw,raw,rs_Cxx,corr_adjusts,process_defaults,constants)
    rs = pu.baseline_correction(rs,rs_cal,nbins,corr_adjusts,constants)
    raw_combined_lo = np.nan*np.ones(nbins)
    raw_combined_hi = np.nan*np.ones(nbins)
    raw_molecular   = np.nan*np.ones(nbins)
    raw_molecular_i2a =np.nan*np.ones(nbins)
    raw_combined_lo[:nbins-sbin] = np.mean(rs.combined_lo_counts[:,sbin:],0)
    raw_combined_hi[:nbins-sbin] = np.mean(rs.combined_hi_counts[:,sbin:],0)
    raw_molecular[:nbins-sbin]   = np.mean(rs.molecular_counts[:,sbin:],0) * corr_adjusts['i2_corr']
    raw_molecular_i2a[:nbins-sbin]   = np.mean(rs.molecular_i2a_counts[:,sbin:],0) * corr_adjusts['i2a_corr']

    nranges = raw_molecular.shape
    print 'raw_molecular.shape =',nranges
    
    while 1:
        string = raw_input('i2a channel clock skew  (bins)  ')
        if len(string) == 0:
           raw.molecular_i2a_counts = temp
           break
        else:
           clk_skew = np.float(string)
           bin_vec = np.arange(len(raw_molecular_i2a))
           temp = raw_molecular_i2a.copy()                    
           temp = np.interp(bin_vec,bin_vec + clk_skew,temp)
           plt.figure(4900)
           plt.plot(raw_molecular_i2a/raw_molecular,bin_ranges/1000,temp/raw_molecular,bin_ranges/1000)
           plt.xlabel('molecular_i2a_counts/molecular_counts')
           plt.ylabel('range (km)')
           plt.grid(True)
    
    plt.figure(5000)
    plt.plot(raw_molecular,bin_ranges/1000.0,'b'
                 ,raw_molecular_i2a,bin_ranges/1000.0,'c')
    ax=plt.gca()
    ax.set_xscale('log')
    ax.legend(['i2','i2a'])
    ax.grid(True)
    plt.title('raw profiles for i2a diff geo corr')
    plt.ylabel('range (km)')
    plt.xlabel('counts/bin/shot')
    #plt.show(block=False)

    
    while 1:
        string = raw_input('i2a normalization range (km),CR to continue, c to clear fig ? ')
        if string == 'c':
           f=plt.figure(5001)
           f.clear()
           string = raw_input('i2a normalization range (km) ? ')    
        if len(string) == 0:
            break 
        else:
            try: 
                       norm_range = np.float(string)
                       norm_bin = 1000.0 * norm_range/(1.5e8*constants['binwidth'])
                       norm_bin = np.int(norm_bin)
   
                       i2a_raw = raw_molecular_i2a
                       mol_raw = raw_molecular

                       i2a = i2a_raw.copy()
                       mol = mol_raw.copy()

                       plt.figure(5001)
                       plt.plot((i2a_raw/mol_raw)*mol_raw[norm_bin]/i2a_raw[norm_bin],bin_ranges/1000.0,'c')
                       ax=plt.gca()
                       ax.set_xlim((.95,1.05)) 
                       ax.grid(True)
                       plt.title('normalized raw i2a/i2 ratio')
                       plt.ylabel('range (km)')
                       plt.xlabel('i2a/i2 ratio')    
                      
    
                       #apply third order, 21 bin filter (150 m filter)
                       window=21
                       i2a=sg.savitzky_golay(i2a,window,3)
                       mol=sg.savitzky_golay(mol,window,3)

                       plt.figure(5002)
                       plt.plot((i2a/mol)*mol[norm_bin]/i2a[norm_bin],bin_ranges/1000.0,'c')
                       ax=plt.gca()
                       ax.set_xlim((.95,1.05)) 
                       ax.grid(True)
                       plt.title('150m filtered i2a/i2 ratio')
                       plt.ylabel('ranges (km)')
                       plt.xlabel('i2a/i2 ratio') 
                       #plt.show(block=False)

                       i2a[:70] = i2a_raw[:70]
                       mol[:70] = mol_raw[:70]

                       plt.figure(5003)
                       plt.plot((i2a/mol)*mol[norm_bin]/i2a[norm_bin],bin_ranges/1000.0,'c')
                       ax=plt.gca()
                       ax.set_xlim((.95,1.05)) 
                       ax.grid(True)
                       plt.title('150m filtered i2a/i2 ratio[70:]')
                       plt.xlabel('i2a/i2 ratio')
                       plt.ylabel('range (km)')
    
    
                       i2a_2=i2a.copy()
                       mol_2=mol.copy()
    
                       #apply third order, 81 bin filter (600 m filter)
                       window=81
                       i2a=sg.savitzky_golay(i2a,window,3)
                       mol=sg.savitzky_golay(mol,window,3)

                       #leave bottom points unfiltered
                       i2a[:100]=i2a_2[:100]
                       mol[:100]=mol_2[:100]
    
                       #bin_vec = np.arange(mol.shape[0])
                       #print bin_vec.shape,i2a.shape,mol.shape

                       plt.figure(5004)
                       plt.plot(i2a,bin_ranges/1000.0,'b',mol,bin_ranges/1000.0,'r')
                       ax=plt.gca()
                       plt.xlabel('Counts')
                       plt.ylabel('range (km)')
                       ax.legend(['mol_i2a','mol'],'upper right')
                       ax.set_xscale('log')
                       ax.grid(True)
                       plt.title('600m filtered profiles')
    
    
                       i2a_mol_ratio=i2a/mol
                       i2a_mol_ratio=i2a_mol_ratio/i2a_mol_ratio[norm_bin]


                       plt.figure(5005)
                       plt.plot((i2a_raw/mol_raw)/i2a_mol_ratio[norm_bin],bin_ranges/1000.0,'c',i2a_mol_ratio,bin_ranges/1000.0,'r')
                       ax=plt.gca()
                       ax.set_xlim((.90,1.05)) 
                       ax.grid(True)
                       plt.xlabel('i2a/i2')
                       plt.ylabel('Range (km)')
                       plt.title('600m filtered i2a/i2 ratio')
                      

                       raw_i2a_mol_ratio = i2a_mol_ratio.copy()

                       #161 bin 1.2km 3rd order filter
                       window =161
                       i2a_mol_ratio = sg.savitzky_golay(i2a_mol_ratio,window,3)
                       i2a_mol_ratio[:150]=raw_i2a_mol_ratio[:150]
    
                       plt.figure(5006)
                       plt.plot(i2a_mol_ratio/i2a_mol_ratio[norm_bin],bin_ranges/1000.0,'r')
                       ax=plt.gca()
                       ax.set_xlim((.95,1.05)) 
                       ax.grid(True)
                       plt.title('1.2km filtered i2a/i2 ratio')    
                       plt.ylabel('range (km)')
                       plt.xlabel('i2a/i2')
    

                       raw_i2a_mol_ratio = i2a_mol_ratio.copy()

                       #500 bin 3.75km 3rd order filter
                       window =501
                       i2a_mol_ratio = sg.savitzky_golay(i2a_mol_ratio,window,3)
                       i2a_mol_ratio[:666]=raw_i2a_mol_ratio[:666]
    
                       plt.figure(5007)
                       plt.plot(i2a_mol_ratio/i2a_mol_ratio[norm_bin],bin_ranges/1000.0,'r')
                       ax=plt.gca()
                       ax.set_xlim((.95,1.05)) 
                       ax.grid(True)
                       plt.title('3.6km filtered i2a/i2 ratio[666:]')    
                       plt.ylabel('range (km)')
                       plt.xlabel('i2a/i2')
                       #plt.show(block=False)
            except:
                print ' '

    while 1:
        string = raw_input('make constant above (km),CR to continue, c to clear fig ? ')
        if string == 'c':
           f=plt.figure(5001)
           f.clear()
           string = raw_input('make constant above (km) ? ')    
        if len(string) == 0:
            break 
        else:
            try:             
                inx = np.int(np.float(string)*1000.0/binwidth)
                i2a_diff_geo = i2a_mol_ratio.copy()
                i2a_diff_geo[inx:]=i2a_diff_geo[inx]
    

                plt.figure(5008)
                plt.plot(i2a_mol_ratio,bin_ranges/1000.0,'c',i2a_diff_geo,bin_ranges/1000.0,'r')
                ax=plt.gca()
                ax.set_xlim((.95,1.05))
                ax.grid(True)
                plt.title('3.75km filtered i2a/i2 ratio[666:]')  
                plt.ylabel('range (km)')
                plt.xlabel('i2a_dif_geo')
                #plt.show(block=False)
            except:
                print ' '
    #write i2a_mol_diff_geofile
    #get path to directory
    dir_path=calibration_path_for(instrument,rs.times[0],process_defaults)
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
    
def make_i2a_diff_geofile_old(instrument,rs,rs_cal,process_defaults,constants,corr_adjusts):
    """compute_i2a_diff_geo_corr(i2a_molecular_counts,molecular_counts,time,rs_cal
       ,processing defaults,constants)
       1) Takes pileup corrected count profiles from data acquired with the i2a filter removed.
       2) Repeats dark and baseline corrections.
       3) Computes the ratio of corrected i2a and i2 counts from time averaged profiles.
    """
    
    import matplotlib.pyplot as plt
   
    #find bin # at range = 0
    [dark_end_time,pulse_time,end_cal_time]=constants['apd_pulse_timing']
    native_res=constants['binwidth']*3e8/2
 
    lidar_bin = np.int((pulse_time *3e8)/native_res)+1
    #restore all channels to uncorrected signals
    rs.molecular_counts   = rs.raw_molecular_counts.copy()
    rs.combined_hi_counts = rs.raw_combined_hi_counts.copy()
    rs.cross_pol_counts   = rs.raw_cross_pol_counts.copy()
    if hasattr(rs,'raw_combined_wfov_counts'):
        rs.combined_wfov_counts = rs.raw_combined_wfov_counts.copy()
    if hasattr(rs,'molecular_i2a_counts'):
        rs.molecular_i2a_counts = rs.raw_molecular_i2a_counts.copy()
    if hasattr(rs,'combined_lo_counts'):
        rs.combined_lo_counts = rs.raw_combined_lo_counts.copy()

    nalts = rs.molecular_counts.shape[1]
    
    if hasattr(rs,'raw_molecular_counts'):
        import matplotlib.pylab as plt
        bins=np.arange(rs.raw_molecular_counts.shape[1])
        plt.figure(5000)
        plt.plot(nanmean(rs.raw_molecular_counts,0),bins,'b'
                 ,nanmean(rs.raw_molecular_i2a_counts,0),bins,'c')
        ax=plt.gca()
        ax.set_xscale('log')
        ax.legend(['i2','i2a'])
        ax.grid(True)
        plt.title('raw profiles for i2a diff geo corr')
        plt.ylabel('range bins')
        plt.xlabel('counts/bin/shot')
        #plt.show(block=False)
        
    #redo dark and baseline corrections
    rs = pu.dark_count_correction(rs,nalts,corr_adjusts,process_defaults,constants)
    print 'energy',rs.transmitted_energy
    rs = pu.baseline_correction(rs,rs_cal,nalts,corr_adjusts,constants)

    plt.figure(5001)
    plt.plot(nanmean(rs.molecular_counts,0),bins,'b'
             ,nanmean(rs.molecular_i2a_counts,0),bins,'c')
    ax=plt.gca()
    ax.grid(True)
    ax.set_xscale('log')
    plt.title('dark and baseline corrected')
    ax.legend(['i2','i2a'])
    plt.ylabel('range bins')
    plt.xlabel('counts/bin/shot')    
    #plt.show(block=False)
    
    norm_range = np.float(raw_input('i2a normalization range (km) ? '))
    
    norm_bin   = 1000.0 * norm_range/(3e8*constants['binwidth'])
    norm_bin = np.int(norm_bin +constants['apd_pulse_timing'][1]/constants['binwidth'])+1
   
    i2a=nanmean(rs.molecular_i2a_counts,0)
    mol=nanmean(rs.molecular_counts,0)
    i2a_raw = i2a.copy()
    mol_raw = mol.copy()

    plt.figure(5002)
    plt.plot((i2a_raw/mol_raw)*mol_raw[norm_bin]/i2a_raw[norm_bin],bins,'c')
    ax=plt.gca()
    ax.set_xlim((.95,1.05)) 
    ax.grid(True)
    plt.title('normalized raw i2a/i2 ratio')
    plt.ylabel('i2a/i2 ratio')
    plt.xlabel('counts/bin/shot')
    #plt.show(block=False)
    
    #apply third order, 21 bin filter (150 m filter)
    window=21
    i2a=sg.savitzky_golay(i2a,window,3)
    mol=sg.savitzky_golay(mol,window,3)

    plt.figure(5003)
    plt.plot((i2a/mol)*mol[norm_bin]/i2a[norm_bin],bins,'c')
    ax=plt.gca()
    ax.set_xlim((.95,1.05)) 
    ax.grid(True)
    plt.title('150m filtered i2a/i2 ratio')
    plt.ylabel('i2a/i2 ratio')
    plt.xlabel('counts/bin/shot') 
    #plt.show(block=False)

    i2a[:70] = i2a_raw[:70]
    mol[:70] = mol_raw[:70]

    plt.figure(5004)
    plt.plot((i2a/mol)*mol[norm_bin]/i2a[norm_bin],bins,'c')
    ax=plt.gca()
    ax.set_xlim((.95,1.05)) 
    ax.grid(True)
    plt.title('150m filtered i2a/i2 ratio[70:]')
    plt.ylabel('i2a/i2 ratio')
    plt.xlabel('counts/bin/shot')
    #plt.show(block=False)
    
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

    plt.figure(5005)
    plt.plot(i2a,bin_vec,'b',mol,bin_vec,'r')
    ax=plt.gca()
    plt.xlabel('Counts')
    plt.ylabel('bin number')
    ax.legend(['mol_i2a','mol'],'upper right')
    ax.set_xscale('log')
    ax.grid(True)
    plt.title('600m filtered profiles')
    #plt.show(block=False)
    
    i2a_mol_ratio=i2a/mol
    i2a_mol_ratio=i2a_mol_ratio/i2a_mol_ratio[norm_bin]


    plt.figure(5006)
    plt.plot((i2a_raw/mol_raw)/i2a_mol_ratio[norm_bin],'c',i2a_mol_ratio,bin_vec,'r')
    ax=plt.gca()
    ax.set_xlim((.90,1.05)) 
    ax.grid(True)
    plt.xlabel('i2a/i2')
    plt.ylabel('Range bin')
    plt.title('600m filtered i2a/i2 ratio')
    #plt.show(block=False)

    raw_i2a_mol_ratio = i2a_mol_ratio.copy()

    #161 bin 1.2km 3rd order filter
    window =161
    i2a_mol_ratio = sg.savitzky_golay(i2a_mol_ratio,window,3)
    i2a_mol_ratio[:150]=raw_i2a_mol_ratio[:150]
    
    plt.figure(5007)
    plt.plot(i2a_mol_ratio,bin_vec,'r')
    ax=plt.gca()
    ax.set_xlim((.95,1.05)) 
    ax.grid(True)
    plt.title('1.2km filtered i2a/i2 ratio')    
    plt.ylabel('range bin')
    plt.xlabel('i2a/i2')
    #plt.show(block=False)
    
    i2a_mol_ratio[:lidar_bin]=np.NaN    
    inx = np.int(raw_input('make constant above bin #? '))    
    i2a_diff_geo = i2a_mol_ratio.copy()
    i2a_diff_geo[inx:]=i2a_diff_geo[inx]
    

    plt.figure(5008)
    plt.plot(i2a_mol_ratio,bin_vec,'c',i2a_diff_geo,bin_vec,'r')
    ax=plt.gca()
    ax.set_xlim((.95,1.05))
    ax.grid(True)
    plt.ylabel('bin number')
    plt.xlabel('i2a_dif_geo')
    #plt.show(block=False)
    
    #write i2a_mol_diff_geofile
    #get path to directory
    dir_path=calibration_path_for(instrument,rs.times[0],process_defaults)
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


def make_temperature_table(instrument,rs_cal,times
           ,rs_constants,corr_adjusts,process_defaults):
    import matplotlib.pyplot as plt
    
    wavelength = np.float(rs_constants['wavelength'])
    table_temperatures = np.arange(160.0,330.0,5.0)
    #table_pressures = np.arange(20.0,1020,20.0)
    table_pressures = np.arange(0.0,1001.,200.0)
    i2a_i2_table = np.zeros((len(table_temperatures),len(table_pressures)))
    print 'computing i2a/i2 temperature table with ',process_defaults.get_value(
                           'molecular_spectrum','model',key='process_defaults')
    print 'i2-scan-file = ',rs_cal.i2scan.filename
    done = False
    freq_offset = 0.0
    while not done:
        if (process_defaults.get_value('molecular_spectrum','model'
                 ,key='process_defaults') == 'tenti_s6'):
            from tenti_s6 import tenti_s6
            plt.figure(197)
            for j in range(len(table_pressures)):
                for i in range(len(table_temperatures)):
                    spectrum = tenti_s6(wavelength * 1e-9,table_temperatures[i],
                           table_pressures[j],
                           (rs_cal.i2scan.data[:, 0]+freq_offset) * 1e9)
                   
                    #ratio of i2a to i2 signals as function of temp and press
                    i2a_i2_table[i,j] = (np.sum(spectrum * rs_cal.i2scan.data[:,6])) \
                         /(np.sum(spectrum *rs_cal.i2scan.data[:,2]))\
                         *(rs_constants['i2a_scan_adjustment']/rs_constants['i2_scan_adjustment'])\
                         *(corr_adjusts['i2a_corr']/corr_adjusts['i2_corr'])
                if j==0 :
                    plt.plot(rs_cal.i2scan.data[:,0],spectrum*
                         rs_cal.i2scan.data[:,6]*rs_constants['i2a_scan_adjustment']
                         *corr_adjusts['i2a_corr'],'r')
               
                    plt.plot(rs_cal.i2scan.data[:,0],spectrum
                         * rs_cal.i2scan.data[:,2]*rs_constants['i2_scan_adjustment']
                         * corr_adjusts['i2_corr'],'c')
                    t_string = ' T= '+repr(table_temperatures[i]) + ' P= ' +repr(table_pressures[j])
                    print 'red and cyan -- ' +t_string
                    
                    print 'ratio = ', i2a_i2_table[i,j]
                    ax=plt.gca()
                    ax.grid(True)
                    #plt.show(block=False)
                   
                if j == len(table_pressures)-1:
                    plt.plot(rs_cal.i2scan.data[:,0]
                            ,spectrum*rs_cal.i2scan.data[:,6] * rs_constants['i2a_scan_adjustment']
                            *corr_adjusts['i2a_corr'] ,'g')         
                    plt.plot(rs_cal.i2scan.data[:,0]
                            ,spectrum*rs_cal.i2scan.data[:,2] * rs_constants['i2_scan_adjustment']
                            * corr_adjusts['i2_corr'],'k')
                    #plt.show(block=False)
                    
                    t_string = 'T= '+repr(table_temperatures[i]) + 'P= ' +repr(table_pressures[j])
                    print 'green and black -- '+t_string
                    print 'ratio = ', i2a_i2_table[i,j]
        elif   process_defaults.get_value('molecular_spectrum','model'
                 ,key='process_defaults') == 'witschas':
            plt.figure(198)
            for j in range(len(table_pressures)):
                for i in range(len(table_temperatures)):
                    spectrum = cu.witschas_spectrum(table_temperatures[i]
                        ,table_pressures[j],wavelength*1e-9
                        ,(rs_cal.i2scan.data[:, 0]+freq_offset) * 1e9)
                    spectrum = spectrum/np.sum(spectrum)
                    #ratio of i2a to i2 signals as function of temp and press
                    i2a_i2_table[i,j] = np.sum(spectrum * rs_cal.i2scan.data[:,6])\
                        /np.sum(spectrum * rs_cal.i2scan.data[:,2])\
                        * rs_constants['i2a_scan_adjustment']/rs_constants['i2_scan_adjustment']\
                        * corr_adjusts['i2a_corr']/corr_adjusts['i2_corr']            
            if j==0 :
                    plt.plot(rs_cal.i2scan.data[:,-1],spectrum,'r')
                    
            elif j == len(table_pressures)-1:    
                    plt.plot(rs_cal.i2scan.data[:,-1],spectrum,'k')
                    ax=plt.gca()
                    ax.grid(True)
            #plt.show(block=False)
            
        elif   process_defaults.get_value('molecular_spectrum','model'
                 ,key='process_defaults') == 'maxwellian':
            # spectral width of molecular scattering
            m_bar = 28.97 * 1.65978e-27  # average mass of an air molecule
            sigma_0 = 1 / (wavelength * 1e-9)  # number in 1/meters
            kb = 1.38044e-23  # Boltzmans constant J/(K deg)
            c = 3e8  # speed of light in m/s
            sigma = (rs_cal.i2scan.data[:, 0]+freq_offset) * 1e9 / c  # wavenumber vector

            for i in range(len(table_temperatures)):
                for j in range(len(table_pressures)):
                    norm = m_bar * c ** 2 / (8 * sigma_0 ** 2 * kb
                        * table_temperatures[ i])
                    spectrum = np.exp(-norm * sigma ** 2)
                    spectrum = spectrum/np.sum(spectrum)
                    #ratio of i2a to i2 signals as function of temp and press
                    i2a_i2_table[i,j] = np.sum(spectrum * rs_cal.i2scan.data[:,6])\
                         /np.sum(spectrum *rs_cal.i2scan.data[:,2])\
                         * rs_constants['i2a_scan_adjustment']/rs_constants['i2_scan_adjustment']\
                         * corr_adjusts['i2a_corr']/corr_adjusts['i2_corr']
   
   
    
        plt.figure(199)
        plt.plot(rs_cal.i2scan.data[:,0],rs_cal.i2scan.data[:,2]
             * rs_constants['i2_scan_adjustment']* corr_adjusts['i2_corr']
             ,'b',rs_cal.i2scan.data[:,0],rs_cal.i2scan.data[:,6]
             * rs_constants['i2a_scan_adjustment'] * corr_adjusts['i2a_corr']
             ,'r')
        ax=plt.gca()
        ax.grid(True)
        plt.title('i2 scan data')
        
        plt.figure(200)
        for i in range(len(table_pressures)):
            plt.plot(table_temperatures,i2a_i2_table[:,i],'b')
        plt.plot(table_temperatures,i2a_i2_table[:,0],'r')
        plt.plot(table_temperatures,i2a_i2_table[:,-1],'k')
        ax=plt.gca()
        ax.grid(True)
        plt.ylabel('i2a/i2 ratio')
        plt.xlabel('temperature (C)')
        plt.title('i2a_i2_table')
        
        #frequency_offset = freq_offset
        string = raw_input('recompute with frequency offset (GHz)?  ')
        if len(string)== 0:
            done = True          
        else:
            freq_offset = np.float(string)
            
    #convert table to temperatures vs i2a_ratio and pressure

    ratio_min=np.max(i2a_i2_table[0,:])
    ratio_max=np.min(i2a_i2_table[-1,:])
    ratios = np.arange(ratio_min,ratio_max,(ratio_max-ratio_min)/200)
    temperatures = np.nan*np.zeros((len(table_pressures),len(ratios)))
    
    #convert table i2a_ratios[temps,press] to temperatures[press,i2a_ratios]
    for i in range(len(table_pressures)):
        for j in range(len(ratios)-1):
             temperatures[i,j] = np.interp(ratios[j],i2a_i2_table[:,i]
                      ,table_temperatures[:])

    text_str = 'pressures='
   
    plt.figure(201)
    colors = ['r','g','b','m','c','k']
    for i in range(len(table_pressures)):
        labl = "%5i"%(np.int(table_pressures[i]))+'mb'
        plt.plot(ratios[:],temperatures[i,:],colors[i],label=labl)
    ax=plt.gca()
    ax.grid(True)
    legend = ax.legend(loc= 'upper left')
    plt.title('temp vs i2a_ratio')
    plt.xlabel('i2a/i2 ratio')
    plt.ylabel('temperature (C)')
    #plt.show(block=False)
    
    print '# i2scan data from--->%s' % rs_cal.i2scan.filename
    print '# table computed using the %s model' \
          %process_defaults.get_value('molecular_spectrum','model',key='process_defaults')
    print '# Frequency offset = %5.2e (GHz)' %(freq_offset)
    print '# temperatures as function of i2a_ratio and pressure'
    print '#ratios    pressures (mb)--->'
    print '%5.0f' %(0.0),
    for i in range(len(table_pressures)):
         print '%5.0f' %(table_pressures[i]),
    print
    for j in range(len(ratios)-1):
        print '%6.4f' %(ratios[j]),        
        for i in range(len(table_pressures)):
           print '%5.1f' %(temperatures[i,j]),
        print
    #write i2a_temp_table
    #get path to directory
    dir_path=calibration_path_for(instrument,times[0],process_defaults)
    #define start time for use in filename from i2scan fileless
    index = rs_cal.i2scan.filename.find('T')
    str = rs_cal.i2scan.filename[index-8:index+5]
    filename=os.path.join(dir_path,'i2a_temp_table_' + str + '.cal')
    fileid=open(filename,'w')
    os.chmod(filename,0664)
    
    print >>fileid, '# file created  %s' %datetime.now()
    print >>fileid, '# table computed using the %s model' \
          %process_defaults.get_value('molecular_spectrum','model',key='process_defaults')
    print >>fileid, '# i2scan data from--->%s' % rs_cal.i2scan.filename
    print >>fileid, '# Frequency offset = %5.2e (GHz)' %(freq_offset)
    print >>fileid, '# temperatures as function of i2a_ratio and pressure'
    print >>fileid, '# i2a_i2_ratios    pressures (mb)--->'

    #print header line
    print >>fileid, '#i2a/i2_ratio',
    for i in range(len(table_pressures)):
        print >>fileid, '%5.0f' %table_pressures[i],
    print >>fileid
    
    print >>fileid, '%5.0f' %(0.0),
    for i in range(len(table_pressures)):
        print >>fileid, '%5.0f' %table_pressures[i],
    print >>fileid
    for j in range(len(ratios)-1):
        print >>fileid, '%6.4f' %(ratios[j]),
        for i in range(len(table_pressures)):
            print >>fileid, '%5.1f' %(temperatures[i,j]),
        print >>fileid    
    fileid.close()
    
    print '\nnew i2a_temp_table=\n '+filename

     
    
    return i2a_i2_table,table_temperatures,table_pressures


 

def make_geofile_new(instrument,profiles,rs_cal,rs_Cxx,rs_constants,corr_adjusts
                     ,process_defaults):
    """
    make_geofile_new(instrument,profiles,rs_cal,rs_Cxx,rs_constants,corr_adjusts,process_defaults)
    instrument         = 'ahsrl','gvhsrl',mf2hsrl,'nshsrl','bagohsrl'
    profiles           = output from generate_profiles selected for zenith of nadir pointing
    rs_cal             = contains cal info from cal files
    rs_Cxx             = HSRL inversion results
    rs_consants        = system constants
    corr_adjusts       = scaling constants for tweeking calibrations
    process_defaults   = structure containing contents of process_control.json file
    """
   
   
    
    #make vector containing altitudes of the range bins
    nbins = len(profiles.raw_molecular_counts[0,:])
    bin_vec = np.arange(nbins)
    binwidth = rs_constants['binwidth'] * 1.5e8  #binwidth in meters
    sbin = np.int(rs_constants['apd_pulse_timing'][1]/rs_constants['binwidth'])
    
    if rs_constants['installation'] == 'ground':
        lidar_altitude = rs_constants['lidar_altitude']
        telescope_dir = 1.0
        
    else:   #airborne
        lidar_altitude =profiles.mean_GPS_MSL_Alt
        print
        print 'min aircraft altitude = ', profiles.min_GPS_MSL_Alt/1000.0, ' km'
        print 'max aircraft altitude = ', profiles.max_GPS_MSL_Alt/1000.0, ' km'

        if profiles.telescope_pointing == 1.0:
            telescope_dir = 1.0
        elif profiles.telescope_pointing == 0:
            telescope_dir = -1.0
        else:
            print
            print '********can not compute geo file with mixed zenith and nadir pointing profiles****'
            print '        no geofile created'
            print
            return
        
    
    bin_altitudes = lidar_altitude + (bin_vec-sbin) * telescope_dir *binwidth\
                * np.cos(np.pi*rs_constants['telescope_roll_angle_offset']/180.0)
        

    #possible fix me---compute new calibrations rather than interpolate existing
    n_cal_alts =  len(rs_Cxx.msl_altitudes)
  
    
    start_time_str = profiles.start_time.strftime('%d-%b-%y %H:%M')
    end_time_str   = profiles.end_time.strftime('%H:%M')
    
    
    #interpolate calibration to bin altitudes

    if telescope_dir > 0:
        cal_indices = np.arange(n_cal_alts)
    else:
        cal_indices = np.arange(n_cal_alts,0,-1)-1
        #bin_altitudes = bin_altitudes[bin_altitudes>=0]
     
   
    Cxx = cu.Cxx()
    #interpolate calibrations to new altitudes

  
    
    Cxx.Cam = rs_Cxx.Cam
    Cxx.Cmc = np.interp(bin_altitudes
           ,rs_Cxx.msl_altitudes, rs_Cxx.Cmm[:n_cal_alts])
    Cxx.Cmm = np.interp(bin_altitudes
           ,rs_Cxx.msl_altitudes, rs_Cxx.Cmm[:n_cal_alts])
    Cxx.beta_r =np.interp(bin_altitudes
            ,rs_Cxx.msl_altitudes, rs_Cxx.beta_r[:n_cal_alts])
    Cxx.times = rs_Cxx.times
    Cxx.msl_altitudes = bin_altitudes
    nalts = len(bin_altitudes)
   
      
   
    
    print ' '
    print 'select normalization altitude where geo_corr == 1'
    print'set approximately 10 km above/below lidar'
    norm_alt=np.float(raw_input('Normalization altitude (km) ?  '))
    norm_alt = norm_alt *1000.0


    
    


    rs = hau.Time_Z_Group()

    
    #get raw profiles and reprocess with range_process
    rs.transmitted_energy =profiles.transmitted_energy/profiles.seeded_shots
    rs.seeded_shots = profiles.seeded_shots.copy()
    rs.molecular_counts = profiles.raw_molecular_counts.copy()
    rs.raw_molecular_counts = profiles.raw_molecular_counts.copy()
    rs.combined_hi_counts = profiles.raw_combined_hi_counts.copy()
    rs.raw_combined_hi_counts = profiles.raw_combined_hi_counts.copy()
    if hasattr(profiles,'combined_lo_counts'):
        rs.combined_lo_counts = profiles.raw_combined_lo_counts.copy()
        rs.raw_combined_lo_counts = profiles.raw_combined_lo_counts.copy()
    rs.cross_pol_counts = profiles.raw_cross_pol_counts.copy()
    rs.raw_cross_pol_counts = profiles.raw_cross_pol_counts.copy()
    if hasattr(profiles,'telescope_pointing'):
        rs.telescope_pointing = profiles.telescope_pointing.copy()
    if hasattr(profiles,'raw_combined_wfov_counts'):
        rs.combined_wfov_counts=profiles.raw_combined_wfov_counts.copy()
        rs.raw_combined_wfov_counts \
                 = profiles.raw_combined_wfov_counts.copy()
    if hasattr(profiles,'raw_molecular_wfov_counts'):
        rs.molecular_wfov_counts = profiles.raw_molecular_wfov_counts.copy()
        rs.raw_molecular_wfov_counts \
                 = profiles.raw_molecular_wfov_counts.copy()
    rs.times=profiles.times.copy()
    rs.delta_t=profiles.delta_t.copy()
    #setting corr_adjusts['geo_corr'] =0 provides r-sqrd geo correction
    corr_adjusts['geo_corr'] = 0.0
   
    if hasattr(rs,'raw_molecular_wfov_counts'):
        plt.figure(294)
        plt.plot(rs.raw_molecular_wfov_counts[0,sbin:]
                 ,bin_altitudes[:nbins-sbin]/1000.0,'b'
                 ,rs.raw_molecular_counts[0,sbin:]
                 ,bin_altitudes[:nbins-sbin]/1000.0,'k')
        plt.grid(True)
        plt.xlabel('mol and mwfov counts')
        plt.ylabel('Altitudes (km)') 
        ax=plt.gca()
        ax.set_xscale('log')
        plt.title(instrument+'  molecular counts  '
              +start_time_str+'-->'+end_time_str)
        #plt.show(block=False)
        
    rs = pu.range_process( instrument, rs, 30.0, rs_constants
        ,rs_cal, False , corr_adjusts ,process_defaults)

    bin_altitudes = bin_vec * rs_constants['binwidth']*1.5e8
    if hasattr(profiles,'molecular_wfov_counts'):
        plt.figure(295)
        plt.plot(rs.molecular_wfov_counts[0,sbin:]
                 * rs_constants['molecular_to_wfov_gain_ratio']
                 ,bin_altitudes[:nbins-sbin]/1000.0,'c'
                 ,rs.molecular_counts[0,sbin:]
                 ,bin_altitudes[:nbins-sbin]/1000.0,'b')
        plt.grid(True)
        plt.xlabel('mol and mwfov counts')
        plt.ylabel('Altitudes (km)') 
        ax=plt.gca()
        ax.set_xscale('log')
        plt.title(instrument+'  molecular counts  '
              +start_time_str+'-->'+end_time_str)
        #plt.show(block=False)
        
    #compute alternate geo correction from wide fov telescope
   
    # ratio= rs.combined_wfov_counts/rs.combined_hi_counts
   
    if hasattr(profiles, 'combined_wfov_counts') \
           or hasattr(profiles,'molecular_wfov_counts'):
        nbins = rs.combined_hi_counts.shape[1]
        wfov_geo_corr = np.zeros(nbins)
        raw_wfov_geo_corr = np.zeros(nbins)
        filtered_mol_wfov = np.zeros(nbins)
        if hasattr(profiles,'combined_wfov_counts'):
            print 'wfov geo_corr computed from combined channels'
            wfov_geo_corr[:nbins-sbin] \
                 = rs.combined_wfov_counts[0,sbin:nbins]\
                 /rs.combined_hi_counts[0,sbin:nbins]
        else:
            print 'wfov geo_corr computed from molecular channels'
            filtered_mol_wfov[:nbins-sbin]\
                  = rs.molecular_wfov_counts[0,sbin:nbins].copy()
            filtered_mol_wfov[200:550] \
                  =sg.savitzky_golay(rs.molecular_wfov_counts[0,200:550]
                  ,11,4, deriv = 0) 
            filtered_mol_wfov[350:] \
                  =sg.savitzky_golay(rs.molecular_wfov_counts[0,350:]
                  ,301,3, deriv = 0)



            #plot wfov counts
            plt.figure(295)
            plt.plot(rs.molecular_wfov_counts[0,sbin:]
                 * rs_constants['molecular_to_wfov_gain_ratio']
                 ,bin_altitudes[:nbins-sbin]/1000.0,'c'
                 ,filtered_mol_wfov[:nbins-sbin]
                 * rs_constants['molecular_to_wfov_gain_ratio']
                 ,bin_altitudes[:nbins-sbin]/1000.0,'r'
                 ,rs.molecular_counts[0,sbin:]
                 ,bin_altitudes[:nbins-sbin]/1000.0,'b')
            plt.grid(True)
            plt.xlabel('mol and mwfov counts')
            plt.ylabel('Altitudes (km)') 
            ax=plt.gca()
            ax.set_xscale('log')
            plt.title(instrument+'  molecular counts  '
                +start_time_str+'-->'+end_time_str)
            #plt.show(block=False)
            
            raw_wfov_geo_corr[:nbins-sbin] \
                  = rs.molecular_wfov_counts[0,sbin:nbins] \
                  * rs_constants['molecular_to_wfov_gain_ratio']\
                                /rs.molecular_counts[0,sbin:nbins]
            wfov_geo_corr[:nbins-sbin] = filtered_mol_wfov[:nbins-sbin]\
                  * rs_constants['molecular_to_wfov_gain_ratio']\
                  /rs.molecular_counts[0,sbin:nbins]

       
        

        if 1:
            temp = np.ones_like(wfov_geo_corr)
            #try smooth ratio rather than wfov counts
            filter_start_bin = 100
            filter_width =201
            temp[filter_start_bin:] \
               =rs_cal.geo.data[filter_start_bin:,1] \
               *sg.savitzky_golay(raw_wfov_geo_corr[filter_start_bin:]\
               /rs_cal.geo.data[filter_start_bin:,1]\
               ,filter_width,3, deriv = 0) 

        #avoid glitch at begining of filter output
        wfov_geo_corr[filter_start_bin+int(filter_width/2.0):] \
              = temp[filter_start_bin+int(filter_width/2.0):]
        #plot wfov geo correction          
        plt.figure(296)
        plt.plot(bin_altitudes/1000.0,raw_wfov_geo_corr,'c'
                 ,bin_altitudes/1000.0,wfov_geo_corr,'b')
        plt.grid(True)
        plt.ylabel('wfov geo correction')
        plt.xlabel('Altitudes (km)') 
        ax=plt.gca()
        ax.set_yscale('log')
        plt.title(instrument+'  wfov geo correction  '
              +start_time_str+'-->'+end_time_str)
  
    #plot raw and range_corrected molecular return           
    plt.figure(297)
    plt.plot(profiles.raw_molecular_counts[0,:],bin_vec
             ,'k',rs.molecular_counts[0,:],bin_vec,'b')
    plt.grid(True)
    plt.xlabel('raw mol, range corr')
    plt.ylabel('Range bin number') 
    ax=plt.gca()
    ax.set_xscale('log')
    plt.title(instrument+'  raw and range corr mol profile  '
              +start_time_str+'-->'+end_time_str)

    
    rs.molecular_counts = rs.molecular_counts[0,0:nalts]
    rs.molecular_counts = rs.molecular_counts[np.newaxis,:]
  


    #plot range_corrected molecular return           
    plt.figure(298)
    plt.plot(rs.molecular_counts[0,:],bin_vec,'b'
               ,rs.combined_counts[0,:],bin_vec,'r'
              ,rs.cross_pol_counts[0,:],bin_vec,'g')
    plt.grid(True)
    plt.xlabel('range corrected')
    plt.ylabel('Range bin number')
    ax=plt.gca()
    ax.set_xscale('log')
    plt.title(instrument+'  range corr mol profile  '+start_time_str+'-->'+end_time_str)

   


    inv = cu.hsrl_inversion(rs, Cxx, rs_constants,corr_adjusts,process_defaults)

   
    #plot inverted returns           
    plt.figure(299)
    plt.plot(inv.Nm[0,:],bin_vec,'b'
               ,Cxx.beta_r,bin_vec,'k')
    plt.grid(True)
    plt.xlabel('Nm, beta_r')
    plt.ylabel('Range bin number')
    ax=plt.gca()
    ax.set_xscale('log')
    plt.title(instrument+'  inverted profiles  '+start_time_str+'-->'+end_time_str)
    #plt.show(block=False)

    #at this point bin_vec has sbin bins before lidar pulse--remove these bins
   

     
    Nm = np.NaN * np.zeros_like(bin_vec)
    beta_r = np.NaN * np.zeros_like(bin_vec)
    beta_a_backscat_par = np.NaN * np.zeros_like(bin_vec)
    linear_depol = np.NaN * np.zeros_like(bin_vec)
    circular_depol = np.NaN * np.zeros_like(bin_vec)


   
    Nm[:(nbins-sbin)] = inv.Nm[0,sbin:]
    beta_r[:(nbins-sbin)] = Cxx.beta_r[sbin:]
    beta_a_backscat_par[:(nbins-sbin)] = inv.beta_a_backscat_par[0,sbin:]
    linear_depol[:(nbins-sbin)] = inv.linear_depol[0,sbin:]
    if hasattr(inv,'circular_depol'):
        circular_depol[:(nbins-sbin)] = inv.circular_depol[0,sbin:]
    
    if telescope_dir >0:
        norm_bin = np.size(bin_altitudes[bin_altitudes <= norm_alt])
    else:
        norm_bin = np.size(bin_altitudes[bin_altitudes >= norm_alt])
    norm_bin = norm_bin - sbin
    
 

    #get beta_r relative to height above/below lidar, beta_r is molecular scattering cross section
  
    beta_r=Cxx.beta_r.copy()

    #compute attenuated beta_r
    att_beta_r=np.zeros_like(beta_r)

    #path length for one bin, with lidar pointed off of zenith
    delta_r = binwidth
    
    #attenuated Rayleigh scattering
    att_beta_r = beta_r * np.exp(-2*np.cumsum(beta_r)*delta_r)
    att_beta_r[0]=Cxx.beta_r[0]
    
    #compute geo correction normalized at norm_alt
    norm=att_beta_r[norm_bin]/Nm[norm_bin]   
  
   
    #plot molecular return and Rayleigh scattering cross section          
    plt.figure(300)
    plt.plot(norm*Nm,bin_vec,'b',beta_r,bin_vec,'c',att_beta_r,bin_vec,'r')
    plt.grid(True)
    plt.xlabel('beta_r, att_beta_r, normalized Nm')
    plt.ylabel('Range bin number') 
    ax=plt.gca()
    ax.set_xscale('log')
    plt.title(instrument+'  normalized mol profile  '+start_time_str+'-->'+end_time_str) 
    #plt.show(block=False)

    
    #compute geo correction normalized at norm_alt
    norm=att_beta_r[norm_bin]/Nm[norm_bin]   
    
    #plot molecular return and Rayleigh scattering cross section          
    plt.figure(300)
    plt.plot(norm*Nm,bin_vec,'b',att_beta_r,bin_vec,'k')
    plt.grid(True)
    plt.xlabel('att_beta_r, normalized Nm')
    plt.ylabel('Range bin number') 
    ax=plt.gca()
    ax.set_xscale('log')
    plt.title(instrument+'  normalized mol profile  '+start_time_str+'-->'+end_time_str)
    #plt.show(block=False)
    
   


    #geo_corr is relative to range from lidar
    geo_corr=att_beta_r/(norm*Nm)
    native_res=rs_constants['binwidth']*3e8/2

    s_Nm=Nm[0:len(geo_corr)]
    s_geo_corr=np.zeros_like(geo_corr)
    s_Nm = polynomial_smoothing(Nm,native_res,[1,3,100,10000])
    s_geo_corr=att_beta_r/(norm*s_Nm)
    s_geo_corr[0:500]=geo_corr[0:500]

   
   
    
    plt.figure(301)
    plt.plot(bin_vec[:4000],geo_corr[:4000],'k')
    plt.grid(True)
    plt.ylabel('geo correction')
    plt.xlabel('Range bin number') 
    #ax=plt.gca()
    #ax.set_xscale('log')
    plt.axis([0,4000,0,6])
    plt.title(instrument+'  raw geo corr profile  '+start_time_str+'-->'+end_time_str)

    #plot final pieced together profile
    plt.figure(302)
    plt.plot(bin_vec,geo_corr,'k',bin_vec,s_geo_corr,'r')
    plt.grid(True)
    plt.ylabel('geo correction')
    plt.xlabel('Range bin number') 
    #ax=plt.gca()
    #ax.set_xscale('log')
    plt.axis([0,4000,0,6])
    plt.title(instrument+'  geometry correction profile  '+start_time_str+'-->'+end_time_str)
    #plt.show(block=False)
    
    print
    print 'Correct for aerosol attenution below the normalization altitude?'
    print 'type <CR> for exit'

    bad_input = True
    while bad_input == True: 
       ref_value=raw_input('Input---lidar_ratio=xx or od=xx?  ')
       photometer_od = 0.0
       lidar_ratio = 0.0
       index = ref_value.find('=')+1
       if ref_value.find('od')==0 and index > 0:
           photometer_od = float(ref_value[index:])
           bad_input = False
       elif ref_value.find('lidar_ratio')==0 and index > 0 :
           lidar_ratio = float(ref_value[index:])
           bad_input = False
       elif len(ref_value) ==  0:
           break
       else:
           print ' '
           print '***************error input must be "lidar_ratio= " or "od= "'
           print ' '
    
    if lidar_ratio > 0 or photometer_od > 0:
        
      
        #get beta_a total as function of hieght above lidar
        if rs_constants['polarization_is_linear']:
            beta_a_backscat = beta_a_backscat_par * (1+profiles.inv.linear_depol)
        else:
            beta_a_backscat = beta_a_backscat_par \
                        * (1 + circular_depol)
        #assume lowest 20 bins are constant    
        beta_a_backscat[:20]=beta_a_backscat[20]   
        

        beta_a_backscat[np.isnan(beta_a_backscat)] = 0.0
        int_backscat = delta_r \
                    * beta_a_backscat.cumsum()
        print 'integrated backscatter = ', int_backscat
        if photometer_od > 0:
            lidar_ratio =  photometer_od/int_backscat[norm_bin]
            print ' '
            print 'Lidar ratio (based on integrated backscatter and supplied OD) = ',lidar_ratio
            print' '
            
        nbins=len(att_beta_r)
        aerosol_extinction = int_backscat[:nbins] * lidar_ratio
       
        #aerosol corrected geo_correction
        plt.figure(303)
        plt.plot(range(0,300),beta_a_backscat[:300],'g',range(0,300),aerosol_extinction[:300],'r',range(0,300),att_beta_r[:300],'b')
        plt.grid(True)
        plt.ylabel('')
        plt.xlabel('Range bin number')
        plt.legend(('beta_a_b','aer ext','att_b'),'lower left')
        #ax=plt.gca()
        #ax.set_xscale('log')
        #plt.axis([0,100,0,6])
        plt.title(instrument+'  aerosol ext, beta_r  '+start_time_str+'-->'+
             end_time_str)

        #aerosol corrected geo_correction
        plt.figure(304)
        plt.plot(bin_vec,aerosol_extinction,'r')
        plt.grid(True)
        plt.ylabel('aerosol_optical depth')
        plt.xlabel('Range bin number')
       
        #ax=plt.gca()
        #ax.set_xscale('log')
        plt.axis([0,4000,0,6])
        plt.title(instrument+'  aerosol optical depth  '+start_time_str+'-->'+
             end_time_str)
        #plt.show(block=False)
       
        att_mol = att_beta_r * np.exp(-2 * aerosol_extinction)

        #corr is relative to height above lidar
        ac_geo_corr = att_mol/(norm*Nm)
        ac_geo_corr = ac_geo_corr/ac_geo_corr[norm_bin]

        #aerosol corrected geo_correction
        plt.figure(305)
        plt.plot(bin_vec,s_geo_corr,'k',bin_vec,ac_geo_corr,'r')
        plt.grid(True)
        plt.ylabel('geo correction')
        plt.xlabel('Range bin number')
        #ax=plt.gca()
        #ax.set_xscale('log')
        plt.axis([0,4000,0,6])
        plt.title(instrument+'  geometry correction profile  '+start_time_str+'-->'+
             end_time_str)
        #plt.show(block=False)
        
        s_geo_corr=ac_geo_corr
        
    #Do you want to make correction constant above a given bin?
    print ' '
    print 'select bin # beyond which geo_corr is constant'
    make_constant_bin=int(raw_input('geo_corr=constant above bin#=?  '))
    #s_geo_corr[make_constant_bin:4001]=np.mean(s_geo_corr[make_constant_bin-10:make_constant_bin])
    mean_value_ext=np.mean(s_geo_corr[make_constant_bin-10:make_constant_bin])
    s_geo_corr=np.hstack((s_geo_corr[0:make_constant_bin],mean_value_ext*np.ones(4001-make_constant_bin)))
                       

    plt.figure(306)
    plt.plot(range(0,make_constant_bin),s_geo_corr[0:make_constant_bin],'k',range(make_constant_bin,3999)\
             ,s_geo_corr[make_constant_bin:3999],'r')
    plt.grid(True)
    plt.ylabel('geo correction')
    plt.xlabel('Range bin number') 
    #ax=plt.gca() ax.set_xscale('log')
    plt.axis([0,4000,0,6])
    plt.title(instrument+'  geo_corr  '+start_time_str+'-->'+end_time_str)

    
    plt.figure(307)
    plt.plot(bin_vec*0.0075,wfov_geo_corr,'k'
             ,bin_vec*.0075,s_geo_corr[:nbins],'r')
    plt.grid(True)
    plt.ylabel('geo correction')
    plt.xlabel('Range (km)') 
    ax=plt.gca()
    ax.set_yscale('log')
    plt.axis([0,30,0,6])
    plt.title(instrument+'  compare normal, wfov geo_corr '+start_time_str+'-->'+end_time_str)
    #plt.show(block=False)
        

   
    #write geofile
    #get path to directory
    dir_path=calibration_path_for(instrument,profiles.times[0],process_defaults)
    #define start time for use in filename
    str=profiles.start_time.strftime("%Y%m%dT%H%M")
    if telescope_dir < 0:
        filename=os.path.join(dir_path,'nadir_geofile_'+str+'.geo')
    else:
        filename=os.path.join(dir_path,'geofile_'+str+'.geo')

    fileid=open(filename,'w')
    os.chmod(filename,0664)
    print >>fileid, '#profile data from %s -->%s UTC' %(start_time_str,end_time_str)
    print >>fileid, '# and cal data from--->%s' % rs_cal.i2scan.filename
    print >>fileid, '# Gains normalized z = %i m' %(norm_alt)
    print >>fileid, '# Geo_corr computed using attenuated beta_r values'
    print >>fileid, '#Range  geo_correction'
    for i in range(np.size(s_geo_corr)):
        print >>fileid, '%10.4e   %10.3e '\
	%(i * binwidth ,s_geo_corr[i])
    fileid.close()
    
    print '\nnew geofile=\n '+filename

     
    return





def make_wfov_geofile(instrument,raw,profiles,rs_Cxx,rs_cal,rs_constants,process_defaults
                      ,corr_adjusts,write_flg = None):
        
    """
    make_wfov_geofile(instrument,profiles,rs_Cxx,rs_cal,rs_Cxx,rs_constants,corr_adjusts,process_defaults)
    Compute geometric correction using the wide field of view channel
    Write correction in /data/instrument/yyyy/mm/dd/calibration/geofile....in yyyy/mm/dd of data used. 
    Reconstruct range bin mol and wfov profiles at raw resolution from altitude bined average raw profiles.
    Do dark count and baseline corrections on these profiles.
    Assume that mol and wfov i2 filter transmissions are identical except for a constant factor.
    Assume that aerosol extinction is given by a constant lidar rato times the aerosol backscatter
    cross section for correction of wfov overlap in lowest layer.

    
    
    instrument         = 'ahsrl','gvhsrl',mf2hsrl,'nshsrl','bagohsrl'
    raw                = raw data structure containing
    profiles           = output from generate_profiles selected for zenith or nadir pointing
    rs_Cxx             = contains molecular extinction profile computed from sounding.
    rs_cal             = contains baseline correction data
    rs_constants       = system constants
    corr_adjusts       = scaling constants for tweeking calibrations
    """
    from scipy.optimize import minimize

    def lin_fit(slope,x,data):
        model = data[0] + (x -x[0]) * slope
        err = np.sum((model-data)**2)
        return err
    
    if not hasattr(profiles,'molecular_wfov_counts'):
        print 'Warning:***********no molecular_wfov channel---can not compute wfov_geo_corr*******'
        return
   
    #make vector containing altitudes of the range bins
    nbins = 4000
    bin_vec = np.arange(nbins)
    binwidth = rs_constants['binwidth'] * 1.5e8  #binwidth in meters
    sbin = np.int(rs_constants['apd_pulse_timing'][1]/rs_constants['binwidth'])
    zenith_angle = np.abs(rs_constants['telescope_roll_angle_offset'])
    lidar_altitude = rs_constants['lidar_altitude']
    bin_ranges = binwidth * np.arange(nbins)
    agl_bin_altitudes = binwidth * np.arange(nbins) * np.cos(np.pi * zenith_angle/180.0)

    rs = hau.Time_Z_Group()


    beta_a_backscat = profiles.inv.beta_a_backscat.copy()
    beta_a_backscat = np.interp(agl_bin_altitudes
                        ,profiles.inv.msl_altitudes-lidar_altitude
                        ,profiles.inv.beta_a_backscat[0,:])
    beta_a_backscat[np.isnan(beta_a_backscat)] = 0.0    
    beta_r = np.interp(agl_bin_altitudes,rs_Cxx.msl_altitudes-lidar_altitude
                       ,rs_Cxx.beta_r)
    beta_r[np.isnan(beta_r)] = 0.0
    mol_atten = np.exp(-2.0 * np.cumsum(beta_r) *binwidth)
    
                          
   
    start_time_str = profiles.start_time.strftime('%d-%b-%y %H:%M')
    end_time_str   = profiles.end_time.strftime('%H:%M')
   

    # generate corrected returns without geo correction.
    # adjust dark correction if desired
    loop = True
    plt.ion()
    while loop == True:
         c_raw = copy.deepcopy(raw)
         rs = copy.deepcopy(raw)
         rs = pu.dark_count_correction(instrument,c_raw,rs,rs_Cxx,corr_adjusts,process_defaults,rs_constants)
         if hasattr(rs_cal,'baseline'):
             rs = pu.baseline_correction(rs,rs_cal,nbins,corr_adjusts,rs_constants)
         else:
             print
             print
             print 'baseline correction was not found in rs_cal---no baseline correction'
             print
         #generate mean profiles, change orgin to timing of laser pulse
         raw_molecular_wfov_counts = np.nan*np.ones(nbins)
         raw_molecular_counts = np.nan*np.ones(nbins)
         raw_molecular_wfov_counts[:nbins-sbin] = nansum(rs.molecular_wfov_counts[:,sbin:],0)     
         raw_molecular_counts[:nbins-sbin] = nansum(rs.molecular_counts[:,sbin:],0)
         nbins = len(raw_molecular_counts)
         corr_molecular_counts = raw_molecular_counts * (binwidth * np.arange(nbins))**2 *1e-6                   
         #from this point onward all variables expressed in terms of range (not altitude)

         #compute alternate geo correction from wide fov telescope   
         wfov_geo_corr = np.zeros(nbins)
         raw_wfov_geo_corr = np.zeros(nbins)
         filtered_mol_wfov = np.zeros(nbins)
         molecular_wfov_counts = raw_molecular_wfov_counts.copy()
         molecular_counts = raw_molecular_counts.copy()
         filtered_mol_wfov = molecular_wfov_counts.copy()
         filtered_mol = molecular_counts.copy()
         dz        = binwidth  #binwidth (m)

         print 'Input data---mol_counts.shape',molecular_counts.shape

        
         plt.figure(193)
         plt.plot(raw_molecular_wfov_counts
                 * rs_constants['molecular_to_wfov_gain_ratio']
                 ,bin_ranges/1000.0,'c'
                 ,raw_molecular_counts
                 ,bin_ranges/1000.0,'b')
         plt.grid(True)
         plt.xlabel('mol and (mwfov * '+str(rs_constants['molecular_to_wfov_gain_ratio']) +') counts')
         plt.ylabel('Range (km)')
         plt.legend(('wfov_mol','mol'),'lower left')
         ax=plt.gca()
         ax.set_xscale('log')
         plt.title(instrument+'  molecular counts  '
              +start_time_str+'-->'+end_time_str)
         plt.show(block=False)
         

       
         filter_width =101
         temp_wfov = sg.savitzky_golay(raw_molecular_wfov_counts,filter_width,2, deriv = 0)
         temp_mol = sg.savitzky_golay(raw_molecular_counts,filter_width,2,deriv = 0)

         start_raw_smooth_bin = 1000
         raw_molecular_wfov_counts[start_raw_smooth_bin:] = temp_wfov[start_raw_smooth_bin:]
         raw_molecular_counts[start_raw_smooth_bin:] = temp_mol[start_raw_smooth_bin:]

         plt.figure(194)
         plt.plot(raw_molecular_wfov_counts
                 * rs_constants['molecular_to_wfov_gain_ratio']
                 ,bin_ranges/1000.0,'c'
                 ,raw_molecular_counts
                 ,bin_ranges/1000.0,'b')
         plt.grid(True)
         plt.xlabel('mol and (mwfov * '+str(rs_constants['molecular_to_wfov_gain_ratio']) +') counts')
         plt.ylabel('Range (km)')
         plt.legend(('wfov_mol','mol'),'lower left')
         ax=plt.gca()
         ax.set_xscale('log')
         plt.title(instrument+'  molecular counts  '
              +start_time_str+'-->'+end_time_str)
         plt.show(block=False)

         loop = False
         """
         print 
         print 'Dark count correction = ', corr_adjusts['mol_wfov_dark_count']
         response =raw_input('new mwfov dark correction? CR for none  ')
         if len(response) == 0:
             loop = False
         else:
             corr_adjusts['mol_wfov_dark_count'] = np.float(response)
             print 'new mol wfov dark = ',corr_adjusts['mol_wfov_dark_count']
         """

    #add overlap correction for wfov channel 
    print
    print 'select range where wfov_mol_counts are matched to integrated backscat * lidar ratio'
    print 'below this altitude geo_corr determined from density profile corrected for extinction'
    print 'using integrated backscatter cross section'
    print
    colors = ['r','b','g','k']
    color_index = 0
    legend_list = []
    while 1:
        
       fit_string = raw_input('fit altitude (km) ? ')
       if len(fit_string)==0:
           break
       if fit_string == 'c':
           print 'clearing figure 197'
           f=plt.figure(197)
           f.clear()
           fit_string = raw_input('fit altitude (km) ? ')
       ratio_string = raw_input('lidar ratio ? ')
       lg_str = 'LR='+ratio_string
       legend_list.append(lg_str)
       lidar_ratio = np.float(ratio_string)
       native_res=rs_constants['binwidth']*3e8/2
    
       #bin altitudes for raw data
       bin_ranges = native_res *np.arange(len(molecular_counts))
       index = len(rs_Cxx.msl_altitudes[rs_Cxx.msl_altitudes <= rs_constants['lidar_altitude']]) 
       
       #hack****************************************************
       #rs_Cxx altitudes and Cmm should have the same length,but don't always
       top_index = rs_Cxx.msl_altitudes.shape[0]

       
       #intepolate Cmm to range bin scale
       Cmm_bins = np.interp(bin_ranges,rs_Cxx.msl_altitudes[index:]-rs_Cxx.msl_altitudes[index],rs_Cxx.Cmm[index:top_index])
       

       try:
                     
          #index of bin_altitudes at fit altitude
           index = len(bin_ranges[bin_ranges/1000.0 \
                                           <= (np.float(fit_string))])
           #low_index = index - 30
           #hi_index = index + 30
           low_index = index - 10
           hi_index = index + 10
           print 'fit between ranges of ' \
                         ,bin_ranges[low_index]/1000.0,'-->'\
                         ,bin_ranges[hi_index]/1000.0, ' km'
                                       

           raw_wfov_counts = molecular_wfov_counts.copy()
           
           expected_wfov_return = Cmm_bins * beta_r *(np.exp(-2.0 * lidar_ratio
                                  * np.cumsum(beta_a_backscat) *binwidth))* mol_atten
           
           corr_wfov_return = raw_wfov_counts * (binwidth * np.arange(nbins))**2 *1e-6

           
                                                     
           scale = np.sum(corr_wfov_return[low_index:hi_index])\
                                 /np.sum(expected_wfov_return[low_index:hi_index])
           expected_wfov_return = expected_wfov_return * scale

           fbtp = process_defaults.get_value('first_bin_to_process','bin_number')
       
           plt.figure(196)
           plt.plot(corr_wfov_return
                    ,bin_ranges/1000.0,'c'
                    ,expected_wfov_return[fbtp:]
                    ,bin_ranges[fbtp:]/1000.0,'r')
           title_str = 'wfov return'
           ax = plt.gca()
           ax.set_xlim(np.min(expected_wfov_return)\
                                     ,np.max(expected_wfov_return))
           plt.grid(True)
           plt.title(title_str)
           plt.ylabel('Range (km)')
           ax.set_xscale('log')
           plt.show(block =False)
           
           #repeat expanded plot near surface
           plt.figure(197)
           plt.plot(corr_wfov_return
                    ,bin_ranges/1000.0,'c'
                    ,expected_wfov_return[fbtp:]
                    ,bin_ranges[fbtp:]/1000.0,colors[color_index])
           title_str = 'Molecular return'
           ax = plt.gca()
           ax.set_ylim(0,bin_ranges[2*hi_index]/1000.0)
           xmax = np.max(expected_wfov_return)
           ax.set_xlim(expected_wfov_return[hi_index],xmax)
           plt.grid(True)
           plt.title(title_str)
           plt.legend(['wfov','geo corr wfov'])
           plt.ylabel('Range (km)')
           #ax.set_xscale('log')
           plt.show(block =False)

           color_index = color_index + 1
           color_index = min(color_index,3)
           
         
        
           print
           print 'enter CR if fit altitude and lidar ratio are ok'
           print 'enter c to clear figure before recalculation'
       except:
           print 'input error--try again'

    raw_wfov_geo_corr = corr_wfov_return/corr_molecular_counts \
                 * rs_constants['molecular_to_wfov_gain_ratio']

    #set to raw ratio
    wfov_geo_corr = raw_wfov_geo_corr.copy()
    #copy expected return computed with integrated backscatter over lowest bins to correct for wfov overlap
    wfov_geo_corr[:index] = expected_wfov_return[:index]/corr_molecular_counts[:index]\
                            * rs_constants['molecular_to_wfov_gain_ratio']
    
    #create a filtered version of the correction that will be used all except closest ranges
    #compress dynamic range with approximate power law for first stage filter--expect corr to be ~ bin_ranges**2
    compression_exponent =1.9
    filter_width =101
    filtered_wfov_corr = sg.savitzky_golay(wfov_geo_corr*bin_ranges**compression_exponent,filter_width,2, deriv = 0)
    filtered_wfov_corr /= bin_ranges**compression_exponent

    #second state of filtering for longer ranges
    start_bin = 800 #6km
    filter_width = 101 #half-width =375m
    filtered_wfov_corr[start_bin:] = sg.savitzky_golay(filtered_wfov_corr[start_bin:],filter_width,1, deriv = 0)
    
    
    #plot raw wfov and corrected geo     
    plt.figure(198)
    plt.plot(bin_ranges[:index]/1000.0,raw_wfov_geo_corr[:index],'c'
             ,bin_ranges[index:]/1000.0,wfov_geo_corr[index:],'b',bin_ranges[index:]/1000.0,filtered_wfov_corr[index:],'r')
    plt.grid(True)
    plt.ylabel('wfov geo corr with overlap corr')
    plt.xlabel('Ranges (km)')
    plt.legend(('raw','overlap corr'),'upper right')
    ax=plt.gca()
    ax.set_ylim(-.5,10)
    #ax.set_xlim(-.05,2.0*bin_ranges[index]/1000.0)
    plt.title(instrument+'  wfov geo correction  '
                 +start_time_str+'-->'+end_time_str)
    plt.show(block =False)
       
    
    
    #if write_flg == None or write_flg == True:
    while 1:
      print ' '
      print 'enter <CR> if done'
      print 'select range to begin Savitzky-Golay filtering'
      print
      start_range = raw_input('start filtering at (km)=? ')
      if len(start_range) == 0:
          #exit on CR
          break
      start_range = np.float(start_range)
      filtered_start_bin = np.int(start_range*1000.0/binwidth)
      print 'start bin = ',filtered_start_bin
      
      print 'select range to begin linear fit to geo_correctin'
      start_fit_range = np.float(raw_input('fit after (km)=?  '))
      start_fit_bin = np.int(start_fit_range * 1000.0/binwidth)
      filtered_ref_bin =666
      print 'start_fit_bin = ',start_fit_bin
      print
      last_good_point = 200
      while not np.isnan(filtered_wfov_corr[last_good_point]):
          last_good_point = last_good_point +1
      last_good_point =1500

      mask = np.isfinite(filtered_wfov_corr)
      x = bin_ranges[start_fit_bin:last_good_point]-bin_ranges[start_fit_bin]
      data = filtered_wfov_corr[start_fit_bin:last_good_point]
      mask = np.isfinite(data)
      print 'ends  ',start_fit_bin, last_good_point
    
      out=minimize(lin_fit,0.0,args=(x[mask],data[mask]),method='Powell')
      slope = out.x
      
      print 'slope = ', slope

      print 'change slope?   --CR for no change'
      string = raw_input('new slope = ? ')
      if not len(string) == 0: 
          slope = np.float(string)    
      extrap = np.ones_like(bin_ranges)
      extrap[start_fit_bin:] = filtered_wfov_corr[start_fit_bin]\
                       + slope * (bin_ranges[start_fit_bin:]-bin_ranges[start_fit_bin])
      #extrap = extrap / filtered_wfov_corr[filtered_start_bin]
      try:
           #plot raw wfov geo correction     
           plt.figure(199)
           plt.plot(bin_ranges/1000.0,wfov_geo_corr,'c'\
             ,bin_ranges[filtered_start_bin:]/1000.0\
             ,filtered_wfov_corr[filtered_start_bin:],'k'\
             ,bin_ranges[start_fit_bin:]/1000.0
             ,extrap[start_fit_bin:],'r' )
           plt.grid(True)
           plt.ylabel('wfov geo corr with overlap corr')
           plt.xlabel('Ranges (km)')
           ax=plt.gca()
           ax.set_ylim(-.5,10)
           plt.title(instrument+'  wfov geo correction  '
                 +start_time_str+'-->'+end_time_str)
           plt.show(block =False)

           #plot corrections with expanded scale     
           plt.figure(200)
           plt.plot(bin_ranges/1000.0,wfov_geo_corr/filtered_wfov_corr[filtered_start_bin],'c'\
             ,bin_ranges[filtered_start_bin:]/1000.0
             ,filtered_wfov_corr[filtered_start_bin:]/filtered_wfov_corr[filtered_start_bin],'g'\
             ,bin_ranges[start_fit_bin:]/1000.0,extrap[start_start_bin:]\
             * filtered_wfov_corr[filtered_start_bin],'r' )
           plt.grid(True)
           plt.ylabel('wfov geo corr with overlap corr')
           plt.xlabel('Ranges (km)')
           plt.legend(('raw','overlap corr'),'upper right')
           ax=plt.gca()
           ax.set_ylim(0,2)
           ax.set_xlim(0.9*bin_ranges[index]/1000.0,1.1*bin_ranges[filtered_start_bin]/1000.0)
           plt.title(instrument+'  wfov geo correction  '
                 +start_time_str+'-->'+end_time_str)
           plt.show(block =False)
    
      except:
           print 'Input error--try again' 
            
     
   

    #add extrapolation to geo correction
    wfov_geo_corr = wfov_geo_corr/filtered_wfov_corr[filtered_ref_bin]
    #filtered_wfov_corr = filtered_wfov_corr/filtered_wfov_corr[lf_start_bin]
    wfov_geo_corr[filtered_start_bin:] = filtered_wfov_corr[filtered_start_bin:] \
               /filtered_wfov_corr[filtered_ref_bin]                          
    wfov_geo_corr[start_fit_bin:] = extrap[start_fit_bin:]

    
    plt.figure(201)
    plt.plot(bin_ranges/1000.0,wfov_geo_corr/raw_wfov_geo_corr,'b')
    plt.grid(True)
    plt.ylabel('wfov / raw wfov correction')
    plt.xlabel('Range (km)')
    ax=plt.gca()
    ax.set_ylim((0.5,1.5))
    plt.title(instrument+' wfov/ raw_wfov '+start_time_str+'-->'+end_time_str)
    
    plt.figure(202)
    plt.plot(bin_ranges/1000.0,wfov_geo_corr,'b')
    plt.grid(True)
    plt.ylabel('wfov geo correction')
    plt.xlabel('Range (km)')
    ax=plt.gca()
    ax.set_yscale('log')
    #plt.axis([0,make_constant_alt+1.0,0,6])
    plt.title(instrument+' wfov geo_corr '+start_time_str+'-->'+end_time_str)
    
   
    plt.figure(203)
    plt.plot(bin_ranges[:4000]/1000.0
        ,wfov_geo_corr[:4000]/(1e6 * rs_cal.geo.data[:4000,1]/rs_cal.geo.data[:4000,0]**2),'r')
    plt.grid(True)
    plt.ylabel('wfov_geo/old')
    plt.xlabel('Range (km)')
    ax=plt.gca()
    ax.set_ylim((.5,1.5))
    #ax.set_yscale('log')
    #plt.axis([0,make_constant_alt+1.0,0,6])
    plt.title(instrument+'  ratio wfov / current geo '+start_time_str+'-->'+end_time_str)
    plt.show(block =False)

    assumptions ={}
    assumptions['lidar_ratio'] =lidar_ratio
    assumptions['lidar_ratio_alt_lmt'] = agl_bin_altitudes[index]+lidar_altitude
    assumptions['lidar_ratio_bin_lmt'] = index

    if write_flg == None or write_flg == True:    
        assumptions['const_above_bin'] = start_fit_bin
        write_geofile(instrument,[profiles.start_time,profiles.end_time],bin_ranges,wfov_geo_corr,'wfov',assumptions,process_defaults=process_defaults)
      
    return bin_ranges, wfov_geo_corr, assumptions


def write_geofile(instrument,times,bin_ranges,geo_corr,geo_type,assumptions = None,process_defaults=None,wfov_mol_ratio = None):
    """write_geofile(instrument,times,bin_ranges,geo_corr,type,assumptions = None, process_defaults= None
                            wfov_mol_ratio = None)
        write geofile in "..../year/mm/dd/calibration/" directory
        instrument = instrument name (eg. bagohsrl)
        times = times used in generating geo correction data
        bin_ranges = vector of bin ranges (m)
        geo_corr = vector of geo corrections at specified ranges
        geo_type = 'std'--computed from ratio of mol return and Rayleigh profile
                 = 'wfov'--computed from ratio of wfov and molecular return 
                 = 'composite'--uses combination of 'std' and 'wfov'
        assumptions = dict['lidar_ratio'] = lidar_ratio
                          ['lidar_ratio_alt_lmt'] = apply lidar_ratio below this alt
                         ['lidar_ratio_bin_lmt'] = bin number for above alt limit
                          ['linear_extrap_bin_lmt'] = do linear extrap above this bin
                          ['linear_extrap_rate']    = slope of linear extrapolation
                          ['wfov_corr_exp'] = wfov needs multiplication  by exp(-2 * range * wfov_corr_exp)
                                              in order to match slope of expected molecular return.
      """


    start_time_str = times[0].strftime('%d-%b-%y %H:%M')
    end_time_str   = times[-1].strftime('%H:%M')

    #get path to directory

    dir_path=calibration_path_for(instrument,times[0],process_defaults)
    #define start time for use in filename
    str=times[0].strftime("%Y%m%dT%H%M")

    filename=os.path.join(dir_path,'geofile_'+str+'.geo')
    fileid=open(filename,'w')
    os.chmod(filename,0664)
    for item in assumptions:
       print item
    print >>fileid, '#'+geo_type+' geo_corr, data from %s -->%s UTC' %(start_time_str,end_time_str)
    print >>fileid, '# file created  %s' %datetime.now()

    print >>fileid, '#lidar ratio assumed = %2.1f below %2.1f km (bin_number=%i)'\
            %(assumptions['lidar_ratio'],assumptions['lidar_ratio_alt_lmt']\
               ,assumptions['lidar_ratio_bin_lmt'])
    if assumptions.has_key('linear_extrap_above_bin'):
        print >>fileid, '#linear extrapolation at bins greater than %i '\
             %(assumptions['linear_extrap_above_bin'])
        print >>fileid, '#linear_extrapolation_rate = %6.4f'%(assumptions['linear_extrap_rate'])
    if assumptions.has_key('wfov_corr_exp'):
        print >>fileid, '# wfov needs * exp(-2*range * %10.2e) '\
             %(assumptions['wfov_corr_exp'])
        
    print '#'
    
    if not wfov_mol_ratio == None:
        print >>fileid, '#Range       geo_correction    wfov_mol_ratio'
        for i in range(np.size(geo_corr)):
            print >>fileid, '%11.5e   %12.6e   %12.6e ' %(bin_ranges[i] ,geo_corr[i], wfov_mol_ratio[i])       
    else:
        print >>fileid, '#range       geo_correction'
        for i in range(np.size(geo_corr)):
            print >>fileid, '%11.5e   %12.6e'\
                      %(bin_ranges[i] ,geo_corr[i])

    fileid.close()

    print "\nnew geofile created with '" +geo_type+"'"
    print filename
      
    return

def make_composite_geofile(instrument,rs_cal,times,wfov_ranges,wfov_geo,wfov_assumptions
                           ,std_ranges,std_geo,std_assumptions,process_defaults=None):                     
    """
    make_composite_geofile(instrument,rs_cal,times,wfov_ranges,wfov_geo,wfov_assumptions
                           ,std_ranges,std_geo,std_assumptions)
    generate a composite geofile using wfov corrected for extinction estimated from integrate
    backscatter for ranges shorter than a specified splice point
    and a comparison of the molecular return to the expected clear air return beyound the
    splice point.
    """

    start_time_str = times[0].strftime('%d-%b-%y %H:%M')
    end_time_str   = times[-1].strftime('%H:%M')
    
    plt.figure(401)
    plt.plot(wfov_ranges[:4000]/1000.0,wfov_geo[:4000],'c'
             ,std_ranges[:4000]/1000.0,std_geo[:4000],'b')
    plt.grid(True)
    plt.ylabel('geo')
    plt.xlabel('Range (km)')
    ax=plt.gca()
    ax.set_ylim((0,5))
    plt.legend(('wfov','std'),'upper right')
    #ax.set_yscale('log')
    #plt.axis([0,make_constant_alt+1.0,0,6])
    plt.title(instrument+'  wfov and geo_corr '+start_time_str+'-->'+end_time_str)

    #filter_width =101
    #temp0=sg.savitzky_golay(wfov_geo,filter_width,4, deriv = 0)
    #temp0[:filter_width] = wfov_geo[:filter_width]

    plt.figure(402)
    plt.plot(wfov_ranges[:4000]/1000.0,wfov_geo[:4000]/std_geo[:4000],'c')
     #        ,wfov_ranges[:4000]/1000.0,temp0[:4000]/std_geo[:4000],'b')
    plt.grid(True)
    plt.ylabel('wfov geo / std geo')
    plt.xlabel('Range (km)')
    ax=plt.gca()
    ax.set_ylim((0,2))
    #ax.set_yscale('log')
    #plt.axis([0,make_constant_alt+1.0,0,6])
    plt.title(instrument+'  wfov_geo/std_corr '+start_time_str+'-->'+end_time_str)

    
    #combine wfov and std geo files into composite
    print 'select splice point between wfov_geo and geo computed from comparison'
    print 'expected clear air return'
    while 1:
  
        string = raw_input('splice at (km) ? ')
        if string == 'c':
                   f=plt.figure(296)
                   f.clear()
                   string = raw_input('splice at (km) ? ')    
        if len(string) == 0:
                   break
        splice_range = np.float(string)
                
                
                    
        #index of bin_altitudes at fit altitude
        index = len(wfov_ranges[wfov_ranges/1000.0 <= splice_range])


        plt.figure(403)
        plt.plot(wfov_ranges[:index]/1000.0,wfov_geo[:index],'c'
             ,std_ranges[index:]/1000.0,std_geo[index:]*wfov_geo[index+1]/std_geo[index+1],'b')
        plt.grid(True)
        plt.ylabel('wfov geo / std geo')
        plt.xlabel('Range (km)')
        ax=plt.gca()
        ax.set_ylim((0,2))
        #ax.set_yscale('log')
        #plt.axis([0,make_constant_alt+1.0,0,6])
        plt.title(instrument+'  spliced geofiles'+start_time_str+'-->'+end_time_str)

    
    #geo_corr = np.ones(4000)
    #ranges = np.arange(4000)
    geo_corr = wfov_geo.copy()
    geo_corr[index:4000]=std_geo[index:4000]*wfov_geo[index+1]/std_geo[index+1]

   


    #plot ratio of new file to currently active geofile
    plt.figure(404)
    npts = 4000
    plt.plot(wfov_ranges[:npts]/1000.0,geo_corr[:npts] \
       /(1e6*rs_cal.geo.data[:npts,1]/rs_cal.geo.data[:npts,0]**2),'r')
    plt.grid(True)
    plt.ylabel('newgeo/oldgeo')
    plt.xlabel('Range (km)')
    ax=plt.gca()
    ax.set_ylim((0.9,1.1))
    #ax.set_yscale('log')
    #plt.axis([0,make_constant_alt+1.0,0,6])
    plt.title(instrument+'  ratio: new/default geo '+start_time_str+'-->'+end_time_str)


    #write new geofile
    geo_type = 'composite'    
    write_geofile(instrument,times,wfov_ranges,geo_corr,geo_type,assumptions = wfov_assumptions,process_defaults=process_defaults)
    return    
    """
    make_wfov_geofile(instrument,profiles,rs_Cxx,rs_cal,rs_Cxx,rs_constants,corr_adjusts,process_defaults)
    Compute geometric correction using the wide field of view channel
    Write correction in /data/instrument/yyyy/mm/dd/calibration/geofile....in yyyy/mm/dd of data used. 
    Reconstruct range bin mol and wfov profiles at raw resolution from altitude bined average raw profiles.
    Do dark count and baseline corrections on these profiles.
    Assume that mol and wfov i2 filter transmissions are identical except for a constant factor.
    Assume that aerosol extinction is given by a constant lidar rato times the aerosol backscatter
    cross section for correction of wfov overlap in lowest layer.
    
    instrument         = 'ahsrl','gvhsrl',mf2hsrl,'nshsrl','bagohsrl'
    raw                = raw data structure containing
    profiles           = output from generate_profiles selected for zenith or nadir pointing
    rs_Cxx             = contains molecular extinction profile computed from sounding.
    rs_cal             = contains baseline correction data
    rs_constants       = system constants
    corr_adjusts       = scaling constants for tweeking calibrations
    """
   
    
def make_geofile(instrument,profiles,rs_cal,rs_Cxx,rs_constants,write_flg = None,fit_params = None,process_defaults=None):                     
    """
    make_geofile(instrument,profiles,rs_cal,rs_Cxx,rs_constants)

    Restore input vectors with raw count profiles, set range resolution to raw value, redo raw
    range processing with geo correction providing only 1/rsqd correction. Compare returns to that
    expected from atmosphere containing aerosols with an assumed lidar_ratio to compute a new geo_correction.
    
    instrument         = 'ahsrl','gvhsrl',mf2hsrl,'nshsrl','bagohsrl'
    profiles           = output from generate_profiles selected for zenith of nadir pointing
    rs_cal             = contains cal info from cal files
    rs_Cxx             = HSRL inversion results
    rs_consants        = system constants
    """
    def input_wfov_ext_correction(molecular_wfov_counts,mol_ext,mol_optical_depth,Cmm):
    
         print
         print
         print 'The wfov channel often decreases more slowly than expected'
         print 'the wfov extinction correction adds a exponential attenuation'
         print 'to correct for this. Adjust wfov_ext correction to match'
         print 'clear air modeled return'
         print
    
         while True:
            #if write_flg == None or write_flg == True:
            try:    
                string = raw_input('wfov ext correction (e.g. 1e-6)? ')
                if string == 'c':
                   f=plt.figure(296)
                   f.clear()
                   string = raw_input('wfov ext correction(e.g. 1e-6)? ')    
                if len(string) == 0:
                   break
                filtered_wfov = filter_profile(molecular_wfov_counts)    
                low_index = 950
                hi_index = 1050
                #kluge factor to make wfov and mol slope match
            
                wfov_ext_coef = np.float(string)  #5.0e-6
                wfov_ext = np.exp(-2 * ranges * wfov_ext_coef)\
                                    /np.exp(-2 * ranges[low_index] * wfov_ext_coef)
                print
                print
                print
                print 'KLUDGE---'+str(wfov_ext_coef)+' 1/m extinction applied to wfov in fig 297, 298'
                print
 
                scaling_ext = np.nanmean(filtered_wfov[low_index:hi_index])\
                             /np.nanmean(mol_ext[low_index:hi_index] * Cmm[low_index:hi_index]
                                     * np.exp(-2.0 * mol_optical_depth[low_index:hi_index]))
                                  
                #molecular wfov counts and return computed from sounding with no extinction
                plt.figure(296)
                plt.plot(ranges/1000.0
                               ,filtered_wfov,'c'
                               ,ranges/1000.0
                               ,filtered_wfov * wfov_ext, 'r'
                               ,ranges/1000.0
                               ,mol_ext * Cmm * np.exp(-2.0 * mol_optical_depth)\
                                   *scaling_ext,'k')
                ax = plt.gca()
                #x.set_xlim((0,10))
                ymax = np.nanmax(filtered_wfov)
         
                ax.set_ylim((ymax/100.0, 1.2*ymax))
                plt.grid(True)
                plt.title('wfov and model return')
                plt.xlabel('range (km)')
                plt.legend(('wfov','wfov * corr','model'),loc='lower left')
                plt.ylabel('molecular counts/bin/laser pulse')
                ax.set_yscale('log')
                plt.show(block=False)
            except:
                print 'Input error--try again'
                print
                
         #wfov_ext_coef is the extinction coefficient used to correct the wfov
         #wfov is a range vector, exp(-range*wfov_ext_coef)
         return wfov_ext_coef,wfov_ext
     
    def filter_profile(counts):
                      filtered_mol = counts.copy()
                      #filter with succesively wider windows

                      if 0: #savitzky_golay filter
                          n_filter_width =11
                          m_filter_width =21
                          w_filter_width =61
                          temp0=sg.savitzky_golay(filtered_mol,n_filter_width,4, deriv = 0)
                          temp1=sg.savitzky_golay(temp0,m_filter_width,4, deriv = 0)
                          temp2=sg.savitzky_golay(temp1,w_filter_width,3, deriv = 0)

                          #assemble filtered output from the three filters
                          filtered_mol[index:]=temp0[index:]
                          filtered_mol[(index+8*m_filter_width):] = temp1[(index+8*m_filter_width):]
                          filtered_mol[(index+16*w_filter_width):] = temp2[(index+16*w_filter_width):]
                      if 0: #polynomial filter
                           dz = ranges[1]- ranges[0]
                           smoothing =np.zeros(5)
                           smoothing[0] = True
                           smoothing[1] = 2        #filter order
                           smoothing[2] = 75.0     #filter width at start bin
                           smoothing[3] = 3000.0   #filter width at end bin
                           smoothing[4] = 400      #start bin
                           filtered_mol = polynomial_smoothing(filtered_mol,dz,smoothing)
                      if 1: #spline filter
                          #use Univariated Spline fit
                          bins = np.arange(len(filtered_mol))
                         
                          weights= np.sqrt(profiles.var_raw_molecular_counts[0,:] * profiles.seeded_shots)
                          #error_bounds = nansum(profiles.var_raw_molecular_counts[0,start_bin:end_bin])
                          #weights = np.sqrt(filtered_mol * profiles.seeded_shots)/profiles.seeded_shots
                          #weights = 1e7/bins**0.8 with s=100*len(good_bins)
           
                          error_sum = 0.002*nanmean(np.abs(filtered_mol[1:]-filtered_mol[:-1]))
                          
                          indices = np.arange(len(filtered_mol))
                          indices = indices[np.isnan(filtered_mol) == 0]
                          indices = indices[np.isnan(weights[indices]) == 0]
                          filtered_mol = filtered_mol[indices]
                          weights = weights[indices]
                          #weights = 1.0/weights
                          #weights[:] = 1.0/indices**1.5 + 0.5/len(indices)
                          weights =(500.0/indices**1.5) +0.5/len(indices)
                          weights[:500]=1.0
                          weights = weights/np.mean(weights)
                          good_bins = bins[indices]
                          
                          
                          s = UnivariateSpline(good_bins
                                , filtered_mol,w=weights, k=3,s=error_sum)
                          #s = UnivariateSpline(good_bins[start_bin:end_bin]
                          #      , filtered_mol[start_bin:end_bin],w=weights[start_bin:end_bin], k=3)
                          filtered_mol = s(bins)
                          #print 'knots = ', s.get_knots()
                      
                          
                          #filtered_mol  =  UnivariateSpline(np.arange(len(filtered_mol)), filtered_mol
                          #                 , w=None, bbox=[None, None], k=3, s=None)

                      return filtered_mol
                      
    #make vector containing range bins
    nbins = len(profiles.msl_altitudes)
    ranges = profiles.msl_altitudes - profiles.msl_altitudes[0]
    lidar_altitude = rs_constants['lidar_altitude']
   


    #make vector of mol scattering cross section vs range
    mol_ext = rs_Cxx.beta_r  

    #make vector of mol ext cross section vs bin_altitudes
    beta_a_backscat = profiles.inv.beta_a_backscat[0,:].copy()
    beta_a_backscat[np.isnan(beta_a_backscat)] = 0.0
    Cmm = rs_Cxx.Cmm[:len(ranges)]
    molecular_counts = profiles.molecular_counts[0,:].copy()
    if hasattr(profiles,'molecular_wfov_counts'):
        molecular_wfov_counts = profiles.molecular_wfov_counts[0,:].copy()

    #this integral must be multiplied by the lidar ratio to get optical depth
    #between the lidar and the scattering volume
    integral_beta_a_backscat =np.cumsum(beta_a_backscat) \
                    *(ranges[2]-ranges[1])
    mol_optical_depth = np.cumsum(mol_ext) * (ranges[2]-ranges[1])
    start_time_str = profiles.start_time.strftime('%d-%b-%y %H:%M')
    end_time_str   = profiles.end_time.strftime('%H:%M')
   
    plt.figure(294)
    if fit_params == None:
       plt.plot(beta_a_backscat,ranges/1000.0,'c',integral_beta_a_backscat,ranges/1000.0,'r' \
             ,mol_optical_depth,ranges/1000.0,'b')
       plt.legend(('$\\beta_a$','$\\int \\beta_a dr$)','mol_od'),'upper left')
    else:
       plt.plot(fit_params['lidar_ratio'] * integral_beta_a_backscat,ranges/1000.0
                         ,'r',beta_a_backscat,ranges/1000.0,'c'
                         ,integral_beta_a_backscat,ranges/1000.0,'r' \
                         ,mol_optical_depth,ranges/1000.0,'b')
       legend_str = '${%i}\\int\\beta_a dr$'%(fit_params['lidar_ratio'])        
       plt.legend(('$\\beta_a$',legend_str,'mol_od'),'upper left')
    
    ax=plt.gca()
    ax.set_xscale('log')
    plt.grid(True)
    plt.ylabel('Range (m)')
    plt.show(block =False)
    
    #plot molecular counts
    plt.figure(295)
    plt.plot(molecular_counts
                 ,ranges/1000.0,'c')
    plt.grid(True)
    plt.xlabel('mol counts')
    plt.ylabel('Range (km)') 
    ax=plt.gca()
    ax.set_xscale('log')
    plt.title(instrument+'  molecular counts  '
              +start_time_str+'-->'+end_time_str)
        
    plt.show(block=False)

    if hasattr(profiles,'molecular_wfov_counts'):
        filtered_wfov = filter_profile(molecular_wfov_counts)    
        wfov_ext_coef,wfov_ext = input_wfov_ext_correction(filtered_wfov,mol_ext,mol_optical_depth,Cmm)

    
    #generate a smoothed version of molecular_counts 
    filtered_mol = filter_profile(molecular_counts)
    
    if 1: 
        #compute overlap correction  
        print
        print 'select range where mol_counts is matched to beta_r'
        #print 'aerosol contribution to OD above fit altitude is assumed == 0'

        loop = True 
        while loop:
            if write_flg == None or write_flg == True:
                string = raw_input('fit range (km) ? ')
                if string == 'c':
                   f=plt.figure(297)
                   f.clear()
                   string = raw_input('fit range (km) ? ')    
                if len(string) == 0:
                   break
                fit_range = np.float(string)
                lidar_ratio =np.float(raw_input('lidar_ratio ? '))
            else:
                fit_range = ranges[fit_params['lidar_ratio_bin_lmt']]/1000.0
                lidar_ratio = fit_params['lidar_ratio']
                loop = False
            if 1: #try:
                      #index of bin_altitudes at fit altitude
                      index = len(ranges[ranges/1000.0 \
                                <= fit_range])
                      low_index = index - 10
                      hi_index = index + 10
                      #if write_flg == None or write_flg == True:
                      #    lidar_ratio =np.float(raw_input('lidar_ratio ? '))
                      print 'fitting beta_r between ranges of ' \
                         ,ranges[low_index]/1000.0,'-->'\
                         ,ranges[hi_index]/1000.0
                      print ' lidar ratio = ',lidar_ratio
                      if lidar_ratio >0:
                         print ' P180/4pi = ',1.0/lidar_ratio

                      
                   
                      #scaling is ratio of mol_counts to attenuated Rayleigh scattering 
                      scaling_ext = np.sum(molecular_counts[low_index:hi_index])\
                          / np.sum(mol_ext[low_index:hi_index]
                          * Cmm[low_index:hi_index]        
                          * np.exp(-2.0 * (lidar_ratio
                          *integral_beta_a_backscat[low_index:hi_index]+mol_optical_depth[low_index:hi_index])))
                      scaling_no_ext = np.sum(molecular_counts[low_index:hi_index])\
                          / np.sum(mol_ext[low_index:hi_index]*Cmm[low_index:hi_index]
                          * np.exp(-2.0 * mol_optical_depth[low_index:hi_index]))
                      print 'scaling_no_ext = ',scaling_no_ext,' scaling_ext = ',scaling_ext

                     
                        
                      #molecular counts = cyan
                      #filtered molecular counts =  blue
                      #molecular computed from sounding and integrated backscatter = red
                      #molecularcomputed from sounding with only mol extinction = black
                     
                      plt.figure(297)
                      plt.plot(ranges/1000.0
                           ,molecular_counts,'c'    
                           ,ranges/1000.0
                           ,filtered_mol,'b'
                           ,ranges/1000.0
                           ,mol_ext *Cmm* np.exp(-2 * (lidar_ratio
                           *integral_beta_a_backscat +mol_optical_depth)) * scaling_ext,'r'
                           ,ranges/1000.0
                           ,mol_ext * Cmm * np.exp(-2.0 * mol_optical_depth)
                           *scaling_no_ext,'k')
                      ax = plt.gca()
                      #x.set_xlim((0,10))
                      ax.set_ylim((1e-4, 10))
                      plt.grid(True)
                      plt.title('beta_r fit to mol counts')
                      plt.xlabel('range (km)')
                      plt.ylabel('molecular counts/bin/laser pulse')
                      plt.legend(('mol','filtered mol','model','model_no_ext'),loc='lower left')
                      ax.set_yscale('log')
                      plt.show(block=False)

                      #compare extinction corrected molecular profile with wfov counts
                      #use this to determine lidar ratio
                      if hasattr(profiles,'molecular_wfov_counts'):
                          corr_mol = mol_ext *Cmm* np.exp(-2 * (lidar_ratio
                                *integral_beta_a_backscat +mol_optical_depth))
                          scal = np.nanmean(corr_mol[low_index:hi_index]) \
                                 /np.nanmean(molecular_wfov_counts[low_index:hi_index]
                                             *wfov_ext[low_index:hi_index])
                        
                          print 'figure exists',plt.fignum_exists(298)
                          if plt.fignum_exists(298):
                             print 'clear fig',lidar_ratio 
                             f=plt.figure(298)
                             f.clear()
                             
                          plt.figure(298)
                          plt.plot(ranges/1000.0,molecular_wfov_counts *scal*wfov_ext,'c'
                                ,ranges/1000.0,corr_mol,'r')
                          plt.grid(True)
                          ax = plt.gca()
                          ax.set_xlim(0,ranges[hi_index]/1000.0)
                          ax.set_ylim(corr_mol[hi_index],10*corr_mol[hi_index])
                          plt.legend(('mol_wfov','ext_corr_mol'))
                          ax.set_yscale('log')
                          plt.title('ext corr mol and wfov counts, LR = '+str(lidar_ratio))
                          plt.xlabel('range (km)')
                          plt.show(block=False)

                      if write_flg == None or write_flg == True:
                          print
                          print '    enter CR if fit altitude and lidar_ratio are ok'
                          print '    enter c to clear figure before recalculation'
            if 0: #except:
                      print
                      print '************** input error--try again'
            if 1:        
                   ind =3800
                   overlap_corr = np.ones_like(mol_ext)
                   overlap_corr[:ind] = mol_ext[:ind] \
                           * Cmm[:ind]*np.exp(-2.0\
                           * (lidar_ratio * integral_beta_a_backscat[:ind]\
                           +mol_optical_depth[:ind])) * scaling_ext \
                           /filtered_mol[:ind]
                   no_ext_overlap_corr = np.ones_like(mol_ext)
                   no_ext_overlap_corr[:ind] = mol_ext[:ind]\
                       * np.exp(-2.0*mol_optical_depth[:ind]) \
                       * Cmm[:ind] * scaling_no_ext / molecular_counts[:ind]
                                
                   plt.figure(299)
                   plt.plot(ranges/1000.0,overlap_corr,'b'
                           ,ranges/1000.0,no_ext_overlap_corr,'r')
                   plt.grid(True)
                   plt.title('Geometry correction')
                   plt.xlabel('Range (km)')
                   plt.legend(('geo','no_ext_geo'))
                   ax = plt.gca()
                   ax.set_yscale('log')
                   plt.show(block=False) 

    temp=overlap_corr
    while 1:
        #check if user wants geo_corr = constant above some alt
        string2 = raw_input('make constant at ranges greater than (km) =  ? ')
        print 'CR for return to quit with no change'
        if len(string2) == 0:
            break
        index2 = len(ranges[ranges/1000.0 <= np.float(string2)])
        ave_val = np.mean(overlap_corr[index2-20:index2+20])
        temp = overlap_corr.copy()
        temp[index2:] = ave_val
        
        plt.figure(300)
        plt.plot(ranges/1000.0,temp,'b')
        plt.grid(True)
        plt.title('Geometry correction with constant section')
        plt.xlabel('Range (km)')
        ax = plt.gca()
        ax.set_yscale('log')
        plt.show(block=False)

    overlap_corr = temp
    
    #compare with default overlap correction
    plt.figure(301)
    plt.plot(ranges/1000.0,overlap_corr,'b'
       ,rs_cal.geo.data[:,0]/1000.0
       ,1e6*rs_cal.geo.data[:,1]/rs_cal.geo.data[:,0]**2,'r')
    plt.grid(True)
    plt.ylabel('geo correction')
    plt.xlabel('Range (km)')
    plt.legend(('new','default'))
    ax=plt.gca()
    ax.set_yscale('log')
    #plt.axis([0,make_constant_alt+1.0,0,6])
    plt.title(instrument+' compare with default geo_corr '+start_time_str+'-->'+end_time_str)
    plt.show(block=False)
    
    #ratio overlap/default correction
    plt.figure(302)
    npts = min(len(ranges),len(overlap_corr),4000)
    
    plt.plot(ranges[:npts]/1000.0,overlap_corr[:npts]
       /(1e6*rs_cal.geo.data[:npts,1]/rs_cal.geo.data[:npts,0]**2),'r')
    plt.grid(True)
    plt.ylabel('new_geo_corr/default_geo_corr')
    plt.xlabel('Range (km)')
    plt.show(block=False)
    #ax=plt.gca()
    #ax.set_yscale('log')
    #plt.axis([0,make_constant_alt+1.0,0,6])
    plt.title(instrument+' compare with default geo_corr '+start_time_str+'-->'+end_time_str)    


    assumptions ={}
    assumptions['lidar_ratio'] =lidar_ratio
    assumptions['lidar_ratio_alt_lmt'] = ranges[index]+lidar_altitude
    assumptions['lidar_ratio_bin_lmt'] = index
    assumptions['linear_extrap_above_bin'] = index2
    assumptions['linear_extrap_rate'] = 0.0
    assumptions['wfov_corr_exp'] = wfov_ext_coef
    #compute wfov/mol ratio if wfov channel is available
    if hasattr(profiles,'molecular_wfov_counts'):

        wfov_mol_ratio = profiles.molecular_wfov_counts[0,:] / profiles.molecular_counts[0,:]
        if 1:
            wfov_gain = rs_constants['molecular_to_wfov_gain_ratio']
            plt.figure(303)
            #plt.plot(ranges/1000.0,wfov_mol_ratio,'r')
            plt.plot(ranges/1000.0,profiles.molecular_wfov_counts[0,:],'c',ranges/1000.0,profiles.molecular_counts[0,:],'b'
                     ,ranges/1000.0,wfov_gain * profiles.molecular_wfov_counts[0,:],'r')
            plt.grid(True)
            plt.title('')
            plt.xlabel('Range (km)')
            ax = plt.gca()
            ax.set_yscale('log')
            plt.show(block=False)
        if 0:
            n_filter_width = 101
            temp0=sg.savitzky_golay(profiles.molecular_wfov_counts[0,:],n_filter_width,4, deriv = 0)
            smoothed_wfov_mol_ratio = temp0 / profiles.molecular_counts[0,:]
            temp0[:400] = wfov_mol_ratio[:400]
            plt.figure(303)
            plt.plot(ranges/1000.0,wfov_mol_ratio,'r'
                 ,ranges/1000.0,smoothed_wfov_mol_ratio,'k')
            plt.grid(True)
            plt.title('')
            plt.xlabel('Range (km)')
            ax = plt.gca()
            ax.set_yscale('log')
            plt.show(block=False)

    
    if (write_flg == None or write_flg == True) and 'wfov_mol_ratio' in locals():    
        write_geofile(instrument,[profiles.start_time,profiles.end_time]
                      ,ranges,overlap_corr,'std_geo',assumptions,process_defaults=process_defaults
                      ,wfov_mol_ratio=wfov_mol_ratio)
    elif write_flg == None or write_flg == True:                    #no wfov channel
        write_geofile(instrument,[profiles.start_time,profiles.end_time]
                      ,ranges,overlap_corr,'std_geo',assumptions,process_defaults=process_defaults)
                      
    return ranges,overlap_corr,assumptions





    
   
 

def make_cxx_file(instrument,rs_Cxx,time,process_defaults=None):
    """Write an ascii file of Cmm,Cma,Cam as a function of raw resolution bin #
       in the calibration directory. This is not used in processing and
       is provided for off-line computations. Time must be provided in
       matplotlib days. """

   

    #get path to directory
    start_time_str=time.strftime("%d-%b-%y %H:%M")
    dir_path=calibration_path_for(instrument,time,process_defaults)
    #define start time for use in filename
    str=time.strftime("%Y%m%dT%H%M")          
    filename=os.path.join(dir_path,'Cxx_'+str+'.cal')
    print 'writing %s in calibration directory' %(filename)
    sounding_time_str=rs_Cxx.sounding_time[0].strftime("%Y%m%dT%H%M")
    fileid=open(filename,'w')
    os.chmod(filename,0664)
    print >>fileid, '#Cxx profiles for %s  UTC' %(start_time_str)
    print >>fileid, '#sounding from %s at %s UTC' %(rs_Cxx.sounding_id,sounding_time_str)
    print >>fileid, '#bin      Cmm      Cmc      Cam      beta_r(1/m)'
    for i in range(len(rs_Cxx.Cmm[0,:])):
        print >>fileid, '%6.1   %10.4e        %10.4e        %10.4e      %10.4e'\
	%(rs_Cxx.msl_altitudes[i] ,rs_Cxx.Cmm[0,i],rs_Cxx.Cmc[0,i],rs_Cxx.Cam,rs_Cxx.beta_r[0,i])
    fileid.close()

    plt.figure(290)
    indices=np.arange(len(rs_Cxx.Cmm[0,:]))
    lines=plt.plot(rs_Cxx.Cmm[0,:],indices,'b',rs_Cxx.Cmc[0,:],indices,'r')
    plt.setp(lines[1],linewidth=3)
    plt.grid(True)
    ax=plt.gca()
    plt.legend(('Cmm','Cmc'),'upper right')
    plt.xlabel('Cal coef)')         
    plt.ylabel('Bin number')

    plt.title('Cal coeficients   '+str)
    plt.show(block=False)
    


def make_diff_geofile(instrument,profiles,rs_cal,rs_constants,process_defaults=None,corr_adjusts=None):
    """
    make_diff_geofile(instrument,profiles,rs_cal,rs_constants)
    This uses raw_profile stream as input
    Both iodine cells must be removed
    """
   
    start_time_str=profiles.start_time.strftime("%d-%b-%y %H:%M")
    end_time_str=profiles.end_time.strftime("%H:%M")
   
    #note raw_profiles have pileup correction applied--no other corrections 
    #signals
    comb_hi = np.zeros(len(profiles.sum_molecular_counts[0,:]))
    mol =              np.zeros_like(comb_hi)
    mol_i2a =          np.zeros_like(comb_hi)
    comb_lo =          np.zeros_like(comb_hi)
    #default differential geometry corrections
    diff_geo_hi =      np.ones_like(comb_hi)
    diff_geo_lo =      np.ones_like(comb_hi)
    diff_geo_i2a =     np.ones_like(comb_hi)
    raw_diff_geo_hi =  np.ones_like(comb_hi)
    raw_diff_geo_lo =  np.ones_like(comb_hi)
    raw_diff_geo_i2a = np.ones_like(comb_hi)
    temp_diff_geo_hi = np.ones_like(comb_hi)
    temp_diff_geo_lo = np.ones_like(comb_hi)
    temp_diff_geo_i2a =np.ones_like(comb_hi)
   
   
          
    ref_bin = 1000  #normalize at bin 1000 as first guess
    dark_scale = 1.0
    make_constant_bin =3800
    norm_index = 500
    import matplotlib.pylab as plt
    plt.ion()
    bin_vec=np.arange(len(mol))

    first_pass = True
    while 1:

        #baseline correction
        comb_hi[:] = profiles.sum_combined_hi_counts[0,:]- rs_cal.baseline.data[:,1]\
                *profiles.sum_transmitted_energy
        mol[:] = profiles.sum_molecular_counts[0,:] -rs_cal.baseline.data[:,3]\
                *profiles.sum_transmitted_energy
        if hasattr(profiles,'sum_combined_lo_counts'):
            comb_lo[:] = profiles.sum_combined_lo_counts[0,:]- rs_cal.baseline.data[:,2]\
                 *profiles.sum_transmitted_energy
        if hasattr(profiles,'sum_molecular_i2a_counts'):
             mol_i2a[:] = profiles.sum_molecular_i2a_counts[0,:] -rs_cal.baseline.data[:,5]\
                   *profiles.sum_transmitted_energy
        
        #dark count correction--no signal in dark correction applied
        mol = mol - profiles.sum_mol_dark_counts[0,0]
        comb_hi = comb_hi - profiles.sum_c_hi_dark_counts[0,0] * dark_scale
        raw_diff_geo_hi = comb_hi/mol
        if hasattr(profiles,'sum_combined_lo_counts'):
             comb_lo = comb_lo - profiles.sum_c_lo_dark_counts[0,0]
             raw_diff_geo_lo = comb_lo/mol 
        if hasattr(profiles,'sum_molecular_i2a_counts'):
            mol_i2a =  mol_i2a - profiles.sum_mol_i2a_dark_counts[0,0]
            raw_diff_geo_i2a = mol_i2a/mol 
            

        start_filter_index = 50
        temp_diff_geo_hi[start_filter_index:] \
                  =(sg.savitzky_golay(raw_diff_geo_hi[start_filter_index:]
                  ,51,3, deriv = 0)).copy()
        start_filter_index = 1000
        temp_diff_geo_hi[start_filter_index:] \
                  =sg.savitzky_golay(temp_diff_geo_hi[start_filter_index:]
                  ,201,3, deriv = 0)
        
        if hasattr(profiles,'sum_combined_lo_counts'):
           start_filter_index = 50
           temp_diff_geo_lo[start_filter_index:] \
                  =(sg.savitzky_golay(raw_diff_geo_lo[start_filter_index:]
                  ,51,3, deriv = 0)).copy()
           start_filter_index = 1000
           temp_diff_geo_lo[start_filter_index:] \
                  =sg.savitzky_golay(temp_diff_geo_lo[start_filter_index:]
                       ,201,3, deriv = 0)              
           
        if hasattr(profiles,'sum_molecular_i2a_counts'):
           start_filter_index = 50
           temp_diff_geo_i2a[start_filter_index:] \
                  =sg.savitzky_golay(raw_diff_geo_i2a[start_filter_index:]
                  ,51,3, deriv = 0)
           start_filter_index = 1000
           temp_diff_geo_i2a[start_filter_index:] \
                  =(sg.savitzky_golay(temp_diff_geo_i2a[start_filter_index:]
                  ,201,3, deriv = 0)).copy()
                                     
        plt.figure(191)
        plt.clf()
        plt.plot(raw_diff_geo_hi/temp_diff_geo_hi[norm_index],bin_vec
                 ,'c',temp_diff_geo_hi/temp_diff_geo_hi[norm_index],bin_vec,'r')
        plt.grid(True)
        ax=plt.gca()
        #ax.set_xscale('log')
        #plt.xlabel('mol')
        #ax.set_xlim((1e-7, 1))
        ax.set_xlim((.5,1.5))
        plt.ylabel('Lidar bin number')
        plt.title('diff geo comb hi')

        if hasattr(profiles,'sum_combined_lo_counts'):
            plt.figure(192)
            plt.clf()
            plt.plot(raw_diff_geo_lo/temp_diff_geo_lo[norm_index],bin_vec,'c'
                     ,temp_diff_geo_lo/temp_diff_geo_lo[norm_index],bin_vec,'r')
            plt.grid(True)
            ax=plt.gca()
            ax.set_xlim((.5,1.5))
            plt.ylabel('Lidar bin number')
            plt.title('diff geo comb lo')
            
        if hasattr(profiles,'sum_molecular_i2a_counts'):
            plt.figure(193)
            plt.clf()
            plt.plot(raw_diff_geo_i2a/temp_diff_geo_i2a[norm_index]
                     ,bin_vec,'c',temp_diff_geo_i2a/temp_diff_geo_i2a[norm_index],bin_vec,'r')
            plt.grid(True)
            ax=plt.gca()
            ax.set_xlim((.5,1.5))
            plt.ylabel('Lidar bin number')
            plt.title('diff_geo_i2a')
            
        diff_geo_hi = (temp_diff_geo_hi/temp_diff_geo_hi[norm_index]).copy()
        diff_geo_hi[make_constant_bin:] = diff_geo_hi[make_constant_bin]
        #diff_geo_hi = diff_geo_hi/diff_geo_hi[make_constant_bin]
        temp_diff_geo_hi = temp_diff_geo_hi/temp_diff_geo_hi[norm_index]
 
        if hasattr(profiles,'sum_combined_lo_counts'):
            diff_geo_lo = (temp_diff_geo_lo/temp_diff_geo_lo[norm_index]).copy()
            diff_geo_lo[make_constant_bin:] = diff_geo_lo[make_constant_bin]
            #diff_geo_lo = diff_geo_lo/diff_geo_lo[make_constant_bin]
            temp_diff_geo_lo = temp_diff_geo_lo/temp_diff_geo_lo[norm_index]

        if hasattr(profiles,'sum_molecular_i2a_counts'):
            diff_geo_i2a = (temp_diff_geo_i2a/temp_diff_geo_i2a[norm_index]).copy()
            diff_geo_i2a[make_constant_bin:] = diff_geo_i2a[make_constant_bin]
            #diff_geo_i2a = diff_geo_lo/diff_geo_i2a[make_constant_bin]
            temp_diff_geo_i2a = temp_diff_geo_i2a/temp_diff_geo_i2a[norm_index]
        if not first_pass:    
            plt.figure(194)
            plt.clf()
            plt.plot(temp_diff_geo_hi,bin_vec,'c',diff_geo_hi,bin_vec,'r')
            plt.grid(True)
            ax=plt.gca()
            ax.set_xlim((.5,1.5))
            plt.ylabel('Lidar bin number')
            plt.title('diff geo comb hi')

        if not first_pass and hasattr(profiles,'sum_combined_lo_counts'):
            plt.figure(195)
            plt.clf()
            plt.plot(temp_diff_geo_lo,bin_vec,'c',diff_geo_lo,bin_vec,'r')
            plt.grid(True)
            ax=plt.gca()
            ax.set_xlim((.5,1.5))
            plt.ylabel('Lidar bin number')
            plt.title('diff geo comb lo')
        
        if not first_pass and hasattr(profiles,'sum_molecular_i2a_counts'):
            plt.figure(196)
            plt.clf()
            plt.plot(temp_diff_geo_i2a,bin_vec,'c',diff_geo_i2a,bin_vec,'r')
            plt.grid(True)
            ax=plt.gca()
            ax.set_xlim((.5,1.5))
            plt.ylabel('Lidar bin number')
            plt.title('diff geo molecular i2a')

        first_pass = False
        print
        print 
        print 'normalize differential geometry at specified bin (or CR to continue)'
        ref_bin=raw_input('Set diff_geo==1 at bin#= ?,   ')
        if len(ref_bin)==0:
            break
        else:
            #make_constant_bin = int(ref_bin)
            dark_scale = raw_input('select scaling for combined dark count (1.0=no change) =? ')
            dark_scale = float(dark_scale)
            #make_constant_bin = int(make_constant_bin)
            #norm_index = make_constant_bin

    while 1:        
       
                

        input_str = raw_input('Make constant beyound bin #? (or CR to exit)   ')
        if len(input_str)==0:
            break

        #make constant after requested bin #
        make_constant_bin = int(input_str)
        temp_diff_geo_hi = diff_geo_hi.copy()
        temp_diff_geo_lo = diff_geo_lo.copy()
        temp_diff_geo_i2a = diff_geo_i2a.copy()
        temp_diff_geo_hi[make_constant_bin:] = temp_diff_geo_hi[make_constant_bin]
        temp_diff_geo_lo[make_constant_bin:] = temp_diff_geo_lo[make_constant_bin]
        temp_diff_geo_i2a[make_constant_bin:] = temp_diff_geo_i2a[make_constant_bin]
        
        plt.figure(197)
        plt.plot(diff_geo_hi,'c',temp_diff_geo_hi,bin_vec,'r')
        plt.grid(True)
        ax=plt.gca()
        ax.set_xlim((.5,1.5))
        plt.ylabel('Lidar bin number')
        plt.title('diff geo comb hi')
    
        if hasattr(profiles,'sum_combined_lo_counts'):
            plt.figure(198)
            plt.plot(temp_diff_geo_lo,bin_vec,'r')
            plt.grid(True)
            ax=plt.gca()
            ax.set_xlim((.5,1.5))
            plt.ylabel('Lidar bin number')
            plt.title('diff geo comb lo')
        
        if hasattr(profiles,'sum_molecular_i2a_counts'):
            plt.figure(198)
            plt.plot(temp_diff_geo_i2a,bin_vec,'r')
            plt.grid(True)
            ax=plt.gca()
            ax.set_xlim((.5,1.5))
            plt.ylabel('Lidar bin number')
            plt.title('diff geo mol i2a')

    diff_geo_hi = temp_diff_geo_hi
    diff_geo_lo = temp_diff_geo_lo
    diff_geo_i2a = temp_diff_geo_i2a
    
    #write diff_geofile
    #get path to directory
    dir_path=calibration_path_for(instrument,profiles.times[0],process_defaults)
    #define start time for use in filename
    str=profiles.start_time.strftime("%Y%m%dT%H%M")          
    filename=os.path.join(dir_path,'diff_geofile_'+str+'.geo')
    fileid=open(filename,'w')
    #os.chmod(fileid,0664)
    print >>fileid, '#profile data from %s -->%s UTC' %(start_time_str,end_time_str)
    print >>fileid, '# Gains normalized range bin = %i ' %(make_constant_bin)
    print >>fileid, '#bin #     diff_hi/mol diff_lo/mol diff_i2a/mol'
    for i in range(len(mol)):
        print >>fileid, '%10i   %10.4e   %10.4e   %10.4e'\
	%(i ,diff_geo_hi[i],diff_geo_lo[i],diff_geo_i2a[i])
    fileid.close()
    
    print '\nnew diff_geofile=\n '+filename
    plt.show()
    return


def make_diff_geo_1064_532(instrument, rs, rs_constants,process_defaults):
    """
       make_diff_geo_1064_532(instrument, rs, rs_constants,process_defaults)
       look for Rayleigh layer in overlap region
    """
    #make vector containing altitudes of the range bins
    nbins = len(rs.msl_altitudes)
    diff_geo_1064_532 = np.ones((nbins,2))   
    ranges = rs.msl_altitudes -rs_constants['lidar_altitude']   
    diff_geo_1064_532[:,0] =ranges
    #ratio of 1064/532 nm backscatter corrected for molecular attenuation
    diff_1064_532 = rs.raw_color_ratio[0,:].copy()
    #correct for aerosol attenuation
    A = process_defaults.get_value('color_ratio','angstrom_coef')
    diff_1064_532 *= np.exp( 2 * rs.optical_depth_aerosol[0,:] * (1.0 -0.5**A))
    plt.ion()
    plt.figure(2000)
    plt.plot(diff_1064_532,ranges,'r' \
             ,rs.raw_color_ratio[0,:],ranges,'c')
    plt.grid(True)
    ax=plt.gca()
    ax.set_xlim([0,2.0])
    plt.ylabel('Altitude (m)')
    plt.xlabel('1064 / 532 ratio corrected for mol attenuation')
    plt.title('raw 1064/532 ratio')
        
    
    make_constant_bin = np.int(7000.0/(ranges[2]-ranges[1]))
    get_str=False
    while 1:
        if get_str == True:
            try:
                input_str = raw_input('Set corr = 1  beyound range (km)? (or CR to exit)   ')
            except:
                print 'input err--try again'
            if len(input_str)==0:
               break
            make_constant_bin = np.int(np.float(input_str)*1000.0/(rs.msl_altitudes[2]-rs.msl_altitudes[1]))
        get_str = True   
        f=plt.figure(2001)
        f.clear()
        plt.plot(diff_1064_532[:]/diff_1064_532[make_constant_bin],ranges,'c' \
                 ,diff_1064_532[:make_constant_bin]/diff_1064_532[make_constant_bin],ranges[:make_constant_bin],'r' \
                 ,np.ones_like(ranges[make_constant_bin:]),ranges[make_constant_bin:],'r')
        plt.grid(True)
        ax=plt.gca()
        ax.set_xlim([0,2.0])
        plt.xlabel('1064 / 532 differential geometry')
        plt.ylabel('ranges')
    

    get_str = False
    low_bin = 0
    print
    print
    while 1:
        if get_str == True:
            try:
                input_str = raw_input('Make constant below bin (km)? (or CR to exit)   ')
            except:
                print 'input err--try again'
            if len(input_str)==0:
               break
            low_bin = np.int(np.float(input_str)*1000.0/(ranges[2]-ranges[1]))
        get_str=True   
        f=plt.figure(2002)
        f.clear()
        norm = diff_1064_532[make_constant_bin]
        plt.plot(diff_1064_532[:]/diff_1064_532[make_constant_bin]
                    ,ranges,'c' \
                 ,diff_1064_532[low_bin:make_constant_bin]/norm
                     ,ranges[low_bin:make_constant_bin],'r' \
                 ,np.ones_like(ranges[make_constant_bin:])
                     ,ranges[make_constant_bin:],'r'
                 ,diff_1064_532[low_bin]/norm*np.ones_like(ranges[:low_bin]),ranges[:low_bin],'r')
        plt.grid(True)
        ax=plt.gca()
        ax.set_xlim([0,2.0])
        plt.xlabel('1064 / 532 differential geometry')
        plt.ylabel('ranges')
    
  

    
    diff_geo_1064_532[low_bin:make_constant_bin,1] = \
                 diff_1064_532[low_bin:make_constant_bin]/diff_1064_532[make_constant_bin]

    diff_geo_1064_532[:low_bin,1] = diff_1064_532[low_bin]/norm * np.ones_like(ranges[:low_bin])
                 
    #get path to directory
    dir_path=calibration_path_for(instrument,rs.times[0],process_defaults)
    #define start time for use in filename
    start_time_str=rs.times[0].strftime("%Y%m%dT%H%M")          
    filename=os.path.join(dir_path,'diff_1064_532_geofile_'+start_time_str+'.geo')
    print
    print 'diff_1064_532_geofile written to:  ', filename
    print
    fileid=open(filename,'w')
    os.chmod(filename,0664)

    print >>fileid, '# '+instrument +' profile data from %s -->%s UTC' %(rs.times[0],rs.times[-1])
    print >>fileid, '# diff_geo set to one beyound range = %f5.1 ' %(ranges[make_constant_bin])
    print >>fileid, '# diff_geo set to one below range = %f5.1 ' %(ranges[low_bin])
    #print >>fileid, '#bin i2a_mol/mol'
    print >>fileid, '#range(m) diff_1064_532_geo'

    for i in range(len(diff_geo_1064_532[:,0])):
        print >>fileid, '%7.1f   %10.4e'\
        %(ranges[i] ,diff_geo_1064_532[i,1])
    fileid.close()
    if 1:
        plt.figure(2222)
        plt.plot(diff_geo_1064_532[:,1],diff_geo_1064_532[:,0],'r')
        plt.grid(True)
        ax=plt.gca()
        ax.set_xlim([0,1.5])
        plt.xlabel('1064/532 differential geo correction')
        plt.ylabel('range (m)')
        plt.show()

    
def make_diff_geofile_old(instrument,profiles,rs_cal,rs_constants,process_defaults=None):
    """
    make_diff_geofile(instrument,profiles,rs_cal,rs_constants)"""

    #select fit type
    fitting = '2nd_poly'


    #fit seperate returns and then ratio (if==1 fit ratios)
    fit_seperate_returns = 0


    
    print 'Select bin # at which to normalize Diff_geo correction'
    ref_bin=raw_input('Diff_geo==1 at bin#= ?,   ')
    print 'select bin # above which Diff_geo is constant'
    make_constant_bin=raw_input('Diff_geo=constant above bin#=?')
    ref_bin=int(ref_bin)
    make_constant_bin=int(make_constant_bin)
    start_time_str=profiles.start_time.strftime("%d-%b-%y %H:%M")
    end_time_str=profiles.end_time.strftime("%H:%M")

    #get soundings file extended to 50 km if needed by climatology
    #soundings include one prior to start time to end of file
    #values returned at alt resolution given by  raw resolution of instrument
    interval_start_time=profiles.times[0]
    alt_res=float(rs_constants['binwidth'])*1.5e8
    #rs_soundings = su.sounding_archive(
    #         self.rs_static.instrument,
    #         self.rs_init.rs_constants['sounding_type'],
    #         self.rs_init.rs_constants['sounding_id'],
    #         interval_start_time,
    #         7.5 * n_range_ave,
    #         50000)
    native_res=rs_constants['binwidth']*3e8/2
    rs_soundings = su.sounding_archive(
             instrument,
             rs_constants['sounding_type'],
             rs_constants['sounding_id'],
             interval_start_time,
             np.arange(0,50000,native_res))
    sounding = rs_soundings.profile(interval_start_time,[],[])
    
   
    beta_r=3.78e-6*sounding.pressures[:6000]/sounding.temps[:6000]
    beta_r=beta_r*(532/rs_constants['wavelength'])**4
    #note: beta_r is defined in terms of bins from the surface


  
    
    #bin offset to correct for dark_count interval
    bin_offset=float(rs_constants['apd_pulse_timing'][1])\
        /float(rs_constants['binwidth'])
    bin_offset=int(bin_offset)
    

   
    #note raw profiles have pileup correction applied--no other corrections 
                                                   
    #apply baseline correction 
    comb_hi = profiles.raw_combined_hi_counts[0,:]- rs_cal.baseline.data[:,1]\
              *profiles.transmitted_energy/profiles.seeded_shots
    comb_lo = profiles.raw_combined_lo_counts[0,:] -rs_cal.baseline.data[:,2]\
              *profiles.transmitted_energy/profiles.seeded_shots
    mol     = profiles.raw_molecular_counts[0,:] -rs_cal.baseline.data[:,3]\
              *profiles.transmitted_energy/profiles.seeded_shots


    
    #Dark count correction subtracts too much due to signal
    #from the previous laser pulse in dark interval.
    #Restore the missing signal to mol and combined hi channels
    range_at_dark_interval=1.5e8/float(rs_constants['laser_rep_rate'])
    binwidth_meters=1.5e8*float(rs_constants['binwidth'])
    dark_interval_meters=float(rs_constants['apd_pulse_timing'][0])*1.5e8
    laser_pulse_bin=int(rs_constants['apd_pulse_timing'][1]*1.5e8/binwidth_meters)
    lidar_altitude=float(rs_constants['lidar_altitude'])
    dark_interval_end=int(dark_interval_meters/binwidth_meters)
    

    #dark count correction--no signal in dark correction
    mol_dark_count=np.mean(mol[:dark_interval_end])
    mol=mol - mol_dark_count
    comb_hi_dark_count=np.mean(comb_hi[:dark_interval_end])
    comb_hi=comb_hi - comb_hi_dark_count
    comb_lo_dark_count=np.mean(comb_lo[:dark_interval_end]) 
    comb_lo=comb_lo - comb_lo_dark_count 




    #compute correction for signal from previous laser pulse in dark count
    #fit dark corrected molecular signal with Rayleigh return
    bin_vec=np.arange(mol.shape[0])
   
    
    #plot sounding up to dark interval vs bin_number
    plt.figure(189)

    sounding_bin_offset=int(lidar_altitude/binwidth_meters)-laser_pulse_bin
    indices=np.arange(len(sounding.altitudes[:]))
    top_index= np.max(indices[sounding.altitudes[:]<=sounding.top[0]])
    if sounding_bin_offset >= 0:
        lines=plt.plot(sounding.temps[sounding_bin_offset:5000+sounding_bin_offset]\
              ,np.arange(len(sounding.temps[0:5000])),'k'\
              ,sounding.temps[sounding_bin_offset:top_index]\
              ,np.arange(len(sounding.temps[sounding_bin_offset:top_index])),'r')
    else:
        lines=plt.plot(sounding.temps[abs(sounding_bin_offset):5000]\
              ,np.arange(len(sounding.temps[abs(sounding_bin_offset):5000])),'k'\
              ,sounding.temps[abs(sounding_bin_offset):top_index]\
              ,np.arange(len(sounding.temps[abs(sounding_bin_offset):top_index])),'r')   
    plt.setp(lines[1],linewidth=3)
    plt.grid(True)
    ax=plt.gca()
    plt.legend(('climate','measured'),'upper right')
    plt.xlabel('Temperature (K)')         
    plt.ylabel('lidar bin number')
    #time_str=sounding.times[0].strftime("%d-%m-%y %H:%M")
    time_str=sounding.times.strftime("%d-%m-%y %H:%M")
    plt.title('Sounding   '+time_str)
    plt.show(block=False)
    
    if 1:
        X=Ray_fit(sounding_bin_offset,beta_r,bin_vec)
        #distance between model and measurement
        errfunc = lambda p,x, y: X.sig(p,x) - y
        k0_mol=10000
        start_fit=3000
        end_fit=3500
        k1_mol,success = sco.leastsq(errfunc,k0_mol\
                       ,(bin_vec[start_fit:end_fit] ,mol[start_fit:end_fit]))

        #fit comb_hi return with molecular 
        k0_comb_hi=10000
        k1_comb_hi,success = sco.leastsq(errfunc,k0_comb_hi,(bin_vec[start_fit:end_fit] \
                  ,comb_hi[start_fit:end_fit]))

        plt.figure(191)
        plt.plot(mol,bin_vec,'b',X.sig(k1_mol,bin_vec),bin_vec,'c'\
              ,comb_hi,bin_vec,'r',X.sig(k1_comb_hi,bin_vec),bin_vec,'k')
    
        



        plt.grid(True)
        ax=plt.gca()
        ax.set_xscale('log')
        plt.xlabel('mol')
        ax.set_xlim((1e-7, 1))         
        plt.ylabel('Lidar bin number')
        plt.title('High alt fit to Rayliegh')
        plt.show(block=False)
        
        #extended bin number when return appears in middle of dark count
        mid_dark_bin=1/float(rs_constants['laser_rep_rate']\
                   *rs_constants['binwidth'])
        mid_dark_bin=int(mid_dark_bin-bin_offset/2.0)
    
        #add correction due to signal in dark count 
        print 'mol sig in dark    = ',X.sig(k1_mol,mid_dark_bin), 'photons/bin/shot'
        print 'comb_hi sig in dark= ',X.sig(k1_comb_hi,mid_dark_bin),'photons/bin/shot'
        mol=mol+X.sig(k1_mol,mid_dark_bin)
        comb_hi=comb_hi+X.sig(k1_comb_hi,mid_dark_bin)

   
  
    w0=10
    wm=400
    delta_w=(wm-w0)/float(len(mol))
    native_res=rs_constants['binwidth']*3e8/2

    s_mol=polynomial_smoothing(mol,native_res,[1,3,75,3000])
    s_comb_hi=polynomial_smoothing(comb_hi,native_res,[1,3,75,3000])
    s_comb_lo=polynomial_smoothing(comb_lo,native_res,[1,3,75,3000])
    

   
              
              
   
        
    #compute raw diff geometry corrections
    raw_comb_hi_to_mol=comb_hi/mol
    raw_comb_lo_to_mol=comb_lo/mol
    
    #compute smoothed diff geo corrections
    comb_hi_to_mol=s_comb_hi/s_mol
    comb_lo_to_mol=s_comb_lo/s_mol

    binlength=1.5e8*rs_constants['binwidth']



    #normalize to ref height
    comb_hi_to_mol=comb_hi_to_mol/comb_hi_to_mol[ref_bin]
    comb_lo_to_mol=comb_lo_to_mol/comb_lo_to_mol[ref_bin]

    #normalize to a range of bins around ref height for raw corrections
    w=(w0+ref_bin*delta_w)
    raw_comb_hi_to_mol=raw_comb_hi_to_mol/np.mean(raw_comb_hi_to_mol[(ref_bin-w):(ref_bin+w)])
    raw_comb_lo_to_mol=raw_comb_lo_to_mol/np.mean(raw_comb_lo_to_mol[(ref_bin-w):(ref_bin+w)])

    bin_vec=np.arange(len(mol))
    if fitting == '2nd_poly':
      s_raw_comb_hi_to_mol=raw_comb_hi_to_mol.copy()
      s_raw_comb_lo_to_mol=raw_comb_lo_to_mol.copy()
      for i in range(w0+75,len(mol)):
          w=int(w0+i*delta_w)
          if (w+i)>len(mol):
              start=int(len(mol)-2*w)
              s_raw_comb_hi_to_mol[i]=np.mean(raw_comb_hi_to_mol[start:])
              s_raw_comb_lo_to_mol[i]=np.mean(raw_comb_lo_to_mol[start:])
          else:
              start=int(i-w)
              end=int(i+w)
              x=range(0,2*w)
              p=np.polyfit(x,raw_comb_hi_to_mol[start:end],2)
              s_raw_comb_hi_to_mol[i]=np.polyval(p,w)
              p=np.polyfit(x,raw_comb_lo_to_mol[start:end],2)
              s_raw_comb_lo_to_mol[i]=np.polyval(p,w)   
       
      plt.figure(400)
      plt.plot(raw_comb_hi_to_mol,bin_vec,'r',s_raw_comb_hi_to_mol,bin_vec,'k')
      plt.grid(True)
      plt.axis([0.5, 1.5,0,bin_vec[-1]])
      plt.xlabel('smoothed and raw comb hi to mol ratio')
      ax.set_xlim((1e-7, 1))  
      plt.show(block=False)
      
    #use case with fit after creating ratios
    if fit_seperate_returns:
        diff_geo_hi=comb_hi_to_mol.copy()
        diff_geo_lo=comb_lo_to_mol.copy()      
    else:
        diff_geo_hi=s_raw_comb_hi_to_mol.copy()
        diff_geo_lo=s_raw_comb_lo_to_mol.copy()

    # provide constant value diff_geo corrections above specified bin #
    diff_geo_hi[make_constant_bin:]=diff_geo_hi[make_constant_bin]
    diff_geo_lo[make_constant_bin:]=diff_geo_lo[make_constant_bin]
    
    plt.figure(200)
    plt.plot(raw_comb_hi_to_mol,bin_vec,'r',comb_hi_to_mol,bin_vec,'g'\
             ,diff_geo_hi,bin_vec,'k')
    plt.grid(True)
    plt.xlabel('Diff geo correction')
    plt.ylabel('Range bin number') 
    ax=plt.gca()
    plt.axis([0.5, 1.5,0,bin_vec[-1]])
    start_time_str=profiles.start_time.strftime("%d-%b-%y %H:%M")
    end_time_str=profiles.end_time.strftime("%H:%M")
    plt.title(instrument+'   comb_hi/mol  '+start_time_str+'-->'+end_time_str)
    
    plt.figure(201)
    bin_vec=np.arange(len(mol))
    plt.plot(raw_comb_lo_to_mol,bin_vec,'c',comb_lo_to_mol,bin_vec,'g'\
             ,diff_geo_lo,bin_vec,'k')
    plt.grid(True)
    plt.xlabel('Diff geo correction')
    plt.ylabel('Range bin number') 
    ax=plt.gca()
    #ax.set_xscale('log')
    plt.axis([0.5, 1.5,0,bin_vec[-1]])
    plt.title(instrument+'   comb_lo/mol '+start_time_str+'-->'+end_time_str)
    plt.show(block=False) 

    #write diff_geofile
    #get path to directory
    dir_path=calibration_path_for(instrument,profiles.times[0],process_defaults)
    #define start time for use in filename
    str=profiles.start_time.strftime("%Y%m%dT%H%M")          
    filename=os.path.join(dir_path,'diff_geofile_'+str+'.geo')
    fileid=open(filename,'w')
    os.chmod(filename,0664)
    print >>fileid, '#profile data from %s -->%s UTC' %(start_time_str,end_time_str)
    print >>fileid, '# Gains normalized range bin = %i ' %(ref_bin)
    print >>fileid, '#bin #   diff_hi/mol diff_lo/mol raw_hi/mol raw_lo/mol'
    for i in range(len(mol)):
        print >>fileid, '%10i   %10.4e   %10.4e   %10.4e   %10.4e'\
	%(i ,diff_geo_hi[i],diff_geo_lo[i],raw_comb_hi_to_mol[i],raw_comb_lo_to_mol[i])
    fileid.close()
    
    print '\nnew diff_geofile=\n '+filename
    

    
def make_baseline_file_old(instrument,profiles,process_defaults=None):
    """make_baseline_file(instrument,profiles)"""
   
    mol=profiles.dc_molecular_counts[0,:]
    comb_hi=profiles.dc_combined_hi_counts[0,:]
    comb_lo=profiles.dc_combined_lo_counts[0,:]
    cpol=profiles.dc_cross_pol_counts[0,:]

    #use non-linear regression to fit baseline curves with function fp   

    # Parametric function: 'v' is the parameter vector, 'x' the independent variable
    fp = lambda v, x: v[0]*np.exp(-v[1]*x**v[2])
    # fp2 will have it's amplitude set by fp 1
    fp2= lambda v,a,x0,x: a*np.exp(-v[0]*(x-x0)**v[1])
    
    x=np.arange(0,mol.shape[0])
    #range bin to start fitting
    start_fit=200
    break_pt=1999
    
    # Error function
    e = lambda v, x, y: ((fp(v,x)-y)**2).sum()
    # error function for second fit
    e2= lambda v,a,x0,x,y:((fp2(v,a,x0,x)-y)**2).sum()

    #the bins below start fit will be copied from orginal arrays
    s_mol=mol.copy()
    s_comb_hi=comb_hi.copy()
    #s_comb_lo=comb_lo.copy()
    s_cpol=cpol.copy()
 
    #fit molecular with two functions matched at bin_number = break_pt
    # Initial parameter value
    v0 = [mol[start_fit], .001, 0.5]
    v= sco.fmin(e, v0, args=(x[start_fit:],mol[start_fit:]),maxiter=10000, maxfun=10000)
    s_mol[start_fit:(break_pt+1)]=fp(v,x[start_fit:(break_pt+1)])
    v0=[.001,.5]
    v = sco.fmin(e2,v0, args=(s_mol[break_pt],break_pt,x[break_pt:],s_mol[break_pt:]),maxiter=10000,maxfun=10000)
    s_mol[break_pt:]=fp2(v,s_mol[break_pt],break_pt,x[break_pt:])
    
    
    #fit comb_hi with two functions matched at bin_number = break_pt
    #initail values      
    v0 =[ comb_hi[start_fit],.001,0.5]
    # Fitting
    v = sco.fmin(e, v0, args=(x[start_fit:],comb_hi[start_fit:]),maxiter=10000, maxfun=10000)
    s_comb_hi[start_fit:(break_pt+1)]=fp(v,x[start_fit:(break_pt+1)])
    
    v0=[.001,.5]
    v = sco.fmin(e2,v0, args=(s_comb_hi[break_pt],break_pt,x[break_pt:],comb_hi[break_pt:]),maxiter=10000,maxfun=10000)
    s_comb_hi[break_pt:]=fp2(v,s_comb_hi[break_pt],break_pt,x[break_pt:])
    



    #fit cpol with two functions matched at bin_number = break_pt
    #initail values      
    v0 =[ cpol[start_fit],.001,0.5]
    v = sco.fmin(e, v0, args=(x[start_fit:4000],cpol[start_fit:4000]),maxiter=10000, maxfun=10000)
    s_cpol[start_fit:(break_pt+1)]=fp(v,x[start_fit:(break_pt+1)])
    v0=[.001,0.5]
    v = sco.fmin(e2,v0, args=(s_cpol[break_pt],break_pt,x[break_pt:],cpol[break_pt:]),maxiter=10000,maxfun=10000)
    s_cpol[break_pt:]=fp2(v,s_cpol[break_pt],break_pt,x[break_pt:])


    

    #fit comb_lo with single function--the comb_lo channel is never used at long ranges
    #initail values      
    v0 =[ comb_lo[start_fit],.001,0.5]
    v = sco.fmin(e, v0, args=(x[start_fit:4000],comb_lo[start_fit:4000]),maxiter=10000, maxfun=10000)
    s_comb_lo=fp(v,x)

    s_mol_med=scipysig.medfilt(mol,201)
    s_comb_hi_med=scipysig.medfilt(comb_hi,201)
    s_comb_lo_med=scipysig.medfilt(comb_lo,201)
    s_cpol_med=scipysig.medfilt(cpol,201)
    
    #baseline correction plots

    plt.figure(100)
    plt.rcParams['figure.figsize'] = [9,9]
    bin_vec=np.arange(profiles.raw_molecular_counts.shape[1])
    plt.plot(mol,bin_vec,'b',comb_hi,bin_vec,'r'
          ,comb_lo,bin_vec,'c',cpol,bin_vec,'g'\
         ,s_mol[start_fit:],bin_vec[start_fit:],'r',s_comb_hi[start_fit:],bin_vec[start_fit:],'k'\
         ,s_comb_lo[start_fit:],bin_vec[start_fit:],'k',s_cpol[start_fit:],bin_vec[start_fit:],'k') 
    plt.grid(True)
    plt.xlabel('Counts/shot/bin')
    plt.ylabel('Range bin number') 
    ax=plt.gca()
    ax.set_xscale('log')
    plt.axis([1e-7, 1e-3,0,4000])
    plt.legend(('mol','comb_hi','comb_lo','cpol','s_mol','s_comb_hi','s_comb_lo','s_cpol'),'upper right')

    plt.figure(101)
    plt.plot(mol,bin_vec,'b',s_mol[start_fit:],bin_vec[start_fit:],'r'\
             ,s_mol_med[start_fit:],bin_vec[start_fit:],'b')
    plt.grid(True)
    plt.xlabel('Counts/shot/bin')
    plt.ylabel('Range bin number') 
    ax=plt.gca()
    ax.set_xscale('log')
    plt.axis([1e-7, 1e-3,0,bin_vec[-1]])
    start_time_str=profiles.times[0].strftime("%d-%b-%y %H:%M")
    end_time_str=profiles.times[-1].strftime("%H:%M")
    plt.title(instrument+'   mol baseline corr '+start_time_str+'-->'+end_time_str)
    plt.legend(('mol','mol fit','mol median'),'upper right')

    plt.figure(102)
    plt.plot(comb_hi,bin_vec,'r',s_comb_hi[start_fit:],bin_vec[start_fit:],'k'\
             ,s_comb_hi_med[start_fit:],bin_vec[start_fit:],'b')
    plt.grid(True)
    plt.xlabel('Counts/shot/bin')
    plt.ylabel('Range bin number') 
    ax=plt.gca()
    ax.set_xscale('log')
    plt.axis([1e-7, 1e-3,0,bin_vec[-1]])
    start_time_str=profiles.times[0].strftime("%d-%b-%y %H:%M")
    end_time_str=profiles.times[-1].strftime("%H:%M")
    plt.title(instrument+'  comb_hi baseline corr '+start_time_str+'-->'+end_time_str)
    plt.legend(('comb hi','c hi fit','c hi median'),'upper right')
    
    plt.figure(103)
    plt.plot(cpol,bin_vec,'g'
         ,s_cpol[start_fit:],bin_vec[start_fit:],'k',s_cpol_med[start_fit:],bin_vec[start_fit:],'b')    
    plt.grid(True)
    plt.xlabel('Counts/shot/bin')
    plt.ylabel('Range bin number') 
    ax=plt.gca()
    ax.set_xscale('log')
    plt.axis([1e-7, 1e-3,0,bin_vec[-1]])
    start_time_str=profiles.times[0].strftime("%d-%b-%y %H:%M")
    end_time_str=profiles.times[-1].strftime("%H:%M")
    plt.title(instrument+'  cpol baseline corr '+start_time_str+'-->'+end_time_str)
    plt.legend(('cpol','c pol fit','cpol median'),'upper right')

    plt.figure(104)
    plt.plot(comb_lo,bin_vec,'c'
         ,s_comb_lo[start_fit:],bin_vec[start_fit:],'k',s_comb_lo_med[start_fit:],bin_vec[start_fit:],'b')    
    plt.grid(True)
    plt.xlabel('Counts/shot/bin')
    plt.ylabel('Range bin number') 
    ax=plt.gca()
    ax.set_xscale('log')
    plt.axis([1e-7, 1e-3,0,bin_vec[-1]])
    start_time_str=profiles.start_time.strftime("%d-%b-%y %H:%M")
    end_time_str=profiles.end_time.strftime("%H:%M")
    plt.title(instrument+'  comb_lo baseline corr '+start_time_str+'-->'+end_time_str)
    plt.legend(('comb lo','c lo fit','c lo median'),'upper right')
    
    
    #used filtered profiles at long ranges
    b_index=start_fit
    mol[b_index:]=s_mol[b_index:]
    comb_hi[b_index:]=s_comb_hi[b_index:]
    comb_lo[b_index:]=s_comb_lo[b_index:]
    cpol[b_index:]=s_cpol[b_index:]






    plt.figure(105)
    plt.plot(mol,bin_vec,'b',comb_hi,bin_vec,'r'\
          ,comb_lo,bin_vec,'c',cpol,bin_vec,'g') 
    plt.grid(True)
    plt.xlabel('Counts/shot/bin')
    plt.ylabel('Range bin number') 
    ax=plt.gca()
    ax.set_xscale('log')
    plt.axis([1e-7, 2,0,bin_vec[-1]])
    #start_time_str=profiles.times[0].strftime("%d-%b-%y %H:%M")
    #end_time_str=profiles.times[-1].strftime("%H:%M")
    plt.title(instrument+'   baseline corr '+start_time_str+'-->'+end_time_str)
    plt.legend(('mol','c hi','c lo','c pol'),'upper right')
    plt.show(block=False)

    #write baseline file
    #get path to directory
    dir_path=calibration_path_for(instrument,profiles.times[0],process_defaults)
    ave_energy_per_shot=profiles.transmitted_energy/profiles.seeded_shots
    #define start time for use in filename
    str=profiles.start_time.strftime("%Y%m%dT%H%M")          
    filename=os.path.join(dir_path,'baseline_correction_'+str+'.blc')
    fileid=open(filename,'w')
    os.chmod(filename,0664)
    print >>fileid, '#profile data from %s -->%s UTC' %(start_time_str,end_time_str)
    print >>fileid, '#ave energy per shot= %6.5f    mJ' %(ave_energy_per_shot)
    print >>fileid, '#bin_num  combined_hi   combined_lo    molecular  crosspol'
    nbins=len(mol)
    for i in range(len(mol)):
      print >> fileid, '%6i   %10.3e    %10.3e   %10.3e   %10.3e'\
	 %(i ,comb_hi[i],comb_lo[i],mol[i],cpol[i])
    fileid.close()
    
    print '\nnew baseline file=\n '+filename
    

     

def make_baseline_file(instrument,raw,rs_constants,processing_defaults,corr_adjusts):
    """make_baseline_file(instrument,profiles)"""
   
    start_time_str=raw.times[0].strftime("%d-%b-%y %H:%M")
    end_time_str=raw.times[-1].strftime("%H:%M")
    print
    print 'Baseline computed starting at ',start_time_str
    print 'Baseline end time ',raw.times[-1].strftime("%d-%b-%y %H:%M")
    print
    
    import matplotlib.pylab as plt
    [nshots,nbins]=raw.molecular_counts.shape

    #generate dark corrected profiles--turn off signal_in_dark correction temporarily
    #if it is present

    temp = corr_adjusts['signal_in_dark']
    corr_adjusts['signal_in_dark'] = 0

    rs=copy.deepcopy(raw)
    
    rs = pu.dark_count_correction(instrument,raw,rs,nbins,corr_adjusts
               ,processing_defaults,rs_constants)
        
    corr_adjusts['signal_in_dark'] = temp
    
    #make mean dark corrected profiles
    # and first short window smoothed version
    window =21
    start = 80
    order = 4
    mol = np.mean(rs.molecular_counts,0)/np.mean(rs.seeded_shots)
    s_mol = np.copy(mol)
    s_mol[start:] = sg.savitzky_golay(s_mol[start:]
                  ,window,order, deriv = 0)  
    comb_hi = np.mean(rs.combined_hi_counts,0)/np.mean(rs.seeded_shots)
    s_comb_hi = np.copy(comb_hi)
    s_comb_hi[start:] = sg.savitzky_golay(s_comb_hi[start:]
                  ,window,order, deriv = 0)
    if hasattr(rs,'combined_lo_counts'):
        comb_lo = np.mean(rs.combined_lo_counts,0)/np.mean(rs.seeded_shots)
        s_comb_lo = np.copy(comb_lo)
        s_comb_lo[start:] = sg.savitzky_golay(s_comb_lo[start:]
                  ,window,order, deriv = 0)
    else:
        comb_lo = np.zeros_like(mol)
        s_comb_lo = comb_lo
    if hasattr(rs,'cross_pol_counts'):
        cpol = np.mean(rs.cross_pol_counts,0)/np.mean(rs.seeded_shots)
        s_cpol = np.copy(cpol)
        s_cpol[start:] = sg.savitzky_golay(s_cpol[start:]
                  ,window,order, deriv = 0)
    else:
        cpol = np.zeros_like(mol)
        s_cpol = cpol
    if hasattr(rs,'molecular_i2a_counts'):
        mol_i2a = np.mean(rs.molecular_i2a_counts,0)/np.mean(rs.seeded_shots)
        s_mol_i2a = np.copy(mol_i2a)
        s_mol_i2a[start:] = sg.savitzky_golay(s_mol_i2a[start:]
                  ,window,order, deriv = 0)
    else:
        mol_i2a = np.zeros_like(mol)
        s_mol_i2a = mol_i2a

    if hasattr(rs,'combined_1064_counts'):
        comb_1064 = np.mean(rs.combined_1064_counts,0)/np.mean(rs.seeded_shots)
        s_comb_1064 = np.copy(comb_1064)
        s_comb_1064[start:] = sg.savitzky_golay(s_comb_1064[start:]
                  ,window,order, deriv = 0)
    else:
        comb_1064 = np.zeros_like(mol)
        s_comb_1064 = mol_i2a
        
    bin_vec = np.arange(len(mol))
    
    #baseline correction plots

    plt.figure(100)
    plt.rcParams['figure.figsize'] = [9,9]   
    plt.plot(mol,bin_vec,'b',comb_hi,bin_vec,'r'
          ,comb_lo,bin_vec,'c',cpol,bin_vec,'g',comb_1064,bin_vec,'m'\
         ,s_mol,bin_vec,'r',s_comb_hi,bin_vec,'k'\
         ,s_comb_lo,bin_vec,'k',s_cpol,bin_vec,'k'
         ,s_mol_i2a,bin_vec,'k',s_comb_1064,bin_vec,'m') 
    plt.grid(True)
    plt.xlabel('Counts/shot/bin')
    plt.ylabel('Range bin number') 
    ax=plt.gca()
    ax.set_xscale('log')
    plt.axis([1e-7, 1e-3,0,4000])
    plt.legend(('mol','comb_hi','comb_lo','cpol','mol_i2a','c_1064','s_mol','s_comb_hi','s_comb_lo','s_cpol'
                ,'s_mol_i2a,','s_c_1064'),'upper right')
    plt.show(block=False)

    #second longer window smoothing
    start = 200
    window = 201
    s_mol[start:] = sg.savitzky_golay(s_mol[start:]
                  ,window,order, deriv = 0)
    s_comb_hi[start:] = sg.savitzky_golay(s_comb_hi[start:]
                  ,window,order, deriv = 0)
    s_comb_lo[start:] = sg.savitzky_golay(s_comb_lo[start:]
                  ,window,order, deriv = 0)
    s_cpol[start:] = sg.savitzky_golay(s_cpol[start:]
                  ,window,order, deriv = 0)
    s_mol_i2a[start:] = sg.savitzky_golay(s_mol_i2a[start:]
                  ,window,order, deriv = 0)
    s_comb_1064[start:]=sg.savitzky_golay(s_comb_1064[start:]
                  ,window,order, deriv = 0)                        
    plt.figure(101)
    plt.plot(s_mol,bin_vec,'b'\
             ,s_comb_hi,bin_vec,'r'
             ,s_comb_lo,bin_vec,'c'
             ,s_cpol,bin_vec,'g'
             ,s_mol_i2a,bin_vec,'m'
             ,s_comb_1064,bin_vec,'k')
    plt.grid(True)
    plt.xlabel('Counts/shot/bin')
    plt.ylabel('Range bin number') 
    ax=plt.gca()
    ax.set_xscale('log')
    plt.axis([1e-7, 1e-3,0,bin_vec[-1]])
    plt.title(instrument+'   baseline corr '+start_time_str+'-->'+end_time_str)
    plt.legend(('mol','comb_hi','comb_lo','cpol','mol_i2a','s_1064'),'upper right')
    plt.show(block=False)
    
    #third still longer window smoothing
    start = 600
    window = 301
    s_mol[start:] = sg.savitzky_golay(s_mol[start:]
                  ,window,order, deriv = 0)
    s_comb_hi[start:] = sg.savitzky_golay(s_comb_hi[start:]
                  ,window,order, deriv = 0)
    s_comb_lo[start:] = sg.savitzky_golay(s_comb_lo[start:]
                  ,window,order, deriv = 0)
    s_cpol[start:] = sg.savitzky_golay(s_cpol[start:]
                  ,window,order, deriv = 0)
    s_mol_i2a[start:] = sg.savitzky_golay(s_mol_i2a[start:]
                  ,window,order, deriv = 0)
    s_comb_1064[start:]=sg.savitzky_golay(s_comb_1064[start:]
                  ,window,order, deriv = 0) 
    plt.figure(102)
    plt.plot(s_mol,bin_vec,'b'\
             ,s_comb_hi,bin_vec,'r'
             ,s_comb_lo,bin_vec,'c'
             ,s_cpol,bin_vec,'g'
             ,s_mol_i2a,bin_vec,'m'
             ,s_comb_1064,bin_vec,'k')
    plt.grid(True)
    plt.xlabel('Counts/shot/bin')
    plt.ylabel('Range bin number')
    ax=plt.gca()
    ax.set_xscale('log')
    plt.axis([1e-7, 1e-3,0,bin_vec[-1]])
    plt.title(instrument+'   baseline corr '+start_time_str+'-->'+end_time_str)
    plt.legend(('mol','comb_hi','comb_lo','cpol','mol_i2a','1064'),'upper right')
    plt.show(block=False)


    first_pt = 100
    mol[first_pt:]     = s_mol[first_pt:]
    comb_hi[first_pt:] = s_comb_hi[first_pt:]
    comb_lo[first_pt:] = s_comb_lo[first_pt:]
    cpol[first_pt:]    = s_cpol[first_pt:] 
    mol_i2a[first_pt:] = s_mol_i2a[first_pt:]
    comb_1064[first_pt:]=s_comb_1064[first_pt:]

    if 0:
      plt.figure(105)
      plt.plot(mol,bin_vec,'b',comb_hi,bin_vec,'r'\
          ,comb_lo,bin_vec,'c',cpol,bin_vec,'g') 
      plt.grid(True)
      plt.xlabel('Counts/shot/bin')
      plt.ylabel('Range bin number') 
      ax=plt.gca()
      ax.set_xscale('log')
      plt.axis([1e-7, 2,0,bin_vec[-1]])
      plt.title(instrument+'   baseline corr '+start_time_str+'-->'+end_time_str)
      plt.legend(('mol','c hi','c lo','c pol'),'upper right')


    #write baseline file
    #get path to directory
    dir_path=calibration_path_for(instrument,raw.times[0],processing_defaults)
    ave_energy_per_shot=np.mean(raw.transmitted_energy)\
            /np.mean(raw.seeded_shots)
    if hasattr(raw,'transmitted_1064_energy'):
        ave_1064_energy_per_shot = np.mean(raw.transmitted_1064_energy)\
            /np.mean(raw.seeded_shots)

    #define start time for use in filename
    str=raw.start_time.strftime("%Y%m%dT%H%M")
    filename=os.path.join(dir_path,'baseline_correction_'+str+'.blc')
    try:
        fileid=open(filename,'w')
        os.chmod(filename,0664)
    except:
        print
        print
        print 'unable to open ', filename
        print 'Probably a file protection problem'
        dir_path = os.getcwd()
        print 'writing to filename = '+dir_path+'/baseline_correction_'+str+'.blc'
        print
        print
        filename = dir_path + '/baseline_correction_'+str+'.blc'
        fileid=open(filename,'w')
        os.chmod(filename,0664)
    print >>fileid, '#profile data from %s -->%s UTC' %(start_time_str,end_time_str)
    print >>fileid, '#ave energy per shot= %6.5f    mJ' %(ave_energy_per_shot)
    if hasattr(raw,'transmitted_1064_energy'):
        print >>fileid, '#ave 1064 energy per shot= %6.5f    mJ' %(ave_1064_energy_per_shot)
    print >>fileid, '#bin_num  combined_hi   combined_lo    molecular  crosspol    mol_I2A   comb_1064'
    nbins=len(mol)
    for i in range(len(mol)):
      print >> fileid, '%6i   %10.3e    %10.3e   %10.3e   %10.3e   %10.3e    %10.3e'\
	 %(i ,comb_hi[i],comb_lo[i],mol[i],cpol[i],mol_i2a[i],comb_1064[i])
    fileid.close()
    
    print '\nnew baseline file=\n '+filename
    
    return
     
def make_cw_i2_scan_file(instrument,rs,rs_constants,process_defaults=None):
    import matplotlib.pylab as plt
    print 'entering cw_i2scan'

    #in this mode, cal pulse begins at normal end of cal pulse and extends to end of buffer
    [dark_interval_end_time, laser_pulse_time, cal_pulse_end_time] = \
        rs_constants['apd_pulse_timing']
    bin_duration = rs_constants['binwidth']
    dark_interval_end = int(dark_interval_end_time /bin_duration)
      #read theoretical i2 transmission file
    if 'lock_point_freq_offset' in rs_constants:
        i2_offset_freq = rs_constants['lock_point_freq_offset']
        print 'theoretical I2 absorption frequency increased by ',i2_offset_freq, ' GHz'
    else:
        i2_offset_freq = 0.0
    try:
         [header,i2t]=ru.readascii(locate_file('dfb_i2_spectra.txt'))
         dfb_to_i2_scale = rs_constants['i2_absorption_scale_factor']
         i2t_trans=10**(dfb_to_i2_scale*i2t[:,1])
         i2t_freq = (i2t[:,0] + i2_offset_freq) * 1e9
    except IOError: #file either not found by locate, or error occurred
         print
         print
         print 'i2-default-scan  not found--using "I2cell_272_31_extended_freq.txt"'
         print
         print
         [header,i2t]=ru.readascii(locate_file('I2cell_272_31_extended_freq.txt'))
         i2t_freq  = (i2t[:,0] + i2_offset_freq) * 1e9
         i2t_trans = i2t[:,2]
         
    #place on freq scale with smaller spacing     
    new_freqs = np.arange(np.nanmin(i2t_freq),np.nanmax(i2t_freq),25e6)
    i2t_trans = np.interp(new_freqs,i2t_freq,i2t_trans)
    i2t_freq = new_freqs
    
    try:
        use_superseed_freq = True
        print
        print 'Using superseed temperature to compute frequency shift, %g Ghz/K'\
                                   %(rs_constants['seedlaser_temp_to_freq'])
        print
        if rs.rs_raw.superseedlasercontrollog.shape[1]==11:
          index=7
        elif rs.rs_raw.superseedlasercontrollog.shape[1]==10:
          index=2
        else:
          print
          print
          raise NotImplemented("Can't get frequency from superseed board. %i element not supported" % (rs.rs_raw.superseedlasercontrollog.shape[1]))
          
        interf_freq = (rs.rs_raw.superseedlasercontrollog[:,index] \
                 - rs.rs_raw.superseedlasercontrollog[0,index])\
                 * rs_constants['seedlaser_temp_to_freq'] * 1e9
        
    except:  #use interferometer measured frequency
        print
        print 'no superseeder--using interferometer freq scale'
        print
        use_superseed_freq = False
        interf_freq = rs.rs_raw.interf_freq.copy()
       
   
    cal_with_seed_light = True
    times = rs.rs_raw.times.copy()

    print 'Using direct injection of doulbled seeder light for calibration'


   
    
    #select time interval
    date_str = datetime.strftime(times[0], '%d-%b-%y ')
    start_time_str = raw_input("Start time=? H:M:S or H:M ")
    start_time_str = date_str + start_time_str
    if start_time_str.count(':') == 2:
      start_time = datetime.strptime(start_time_str, '%d-%b-%y %H:%M:%S')
    else:
      start_time = datetime.strptime(start_time_str, '%d-%b-%y %H:%M')
    end_time_str = raw_input("End time H:M:S or H:M ")
    end_time_str = date_str + end_time_str
    if end_time_str.count(':') == 2:
       end_time = datetime.strptime(end_time_str, '%d-%b-%y %H:%M:%S')
    else:
       end_time = datetime.strptime(end_time_str, '%d-%b-%y %H:%M')
   #restrict times via user input
    time_mask = np.array([ ((x >= start_time) and (x<end_time)) 
                             for x in times[:] ])

    
    interf_freq = interf_freq[time_mask]
    molecular_cal_pulse = rs.rs_raw.molecular_counts[time_mask,:]
    combined_hi_cal_pulse = rs.rs_raw.combined_hi_counts[time_mask,:]
    shot_count = rs.rs_raw.shot_count[time_mask]
    short_cell_ratio = rs.rs_raw.filtered_energy[time_mask,0]/rs.rs_raw.nonfiltered_energy[time_mask,0]
    max_short_cell_ratio = rs.rs_raw.filtered_energy[time_mask,2]/rs.rs_raw.nonfiltered_energy[time_mask,0]
    times = times[time_mask]

    if hasattr(rs.rs_raw,'molecular_i2a_counts'):
        molecular_i2a_cal_pulse = rs.rs_raw.molecular_i2a_counts[time_mask,:]


    delta = interf_freq[1] - interf_freq[0]
   
        
    plt.figure(9999)
    plt.plot(np.arange(len(molecular_cal_pulse[0,:])),np.nanmean(molecular_cal_pulse,0))

    
    #remove dark counts
    
    if 0:
       molecular_cal_pulse = np.nanmean(molecular_cal_pulse[:,60:],1) \
                          -2e-5 *shot_count 
       #                   - rs_constants['molecular_detector_dark_count'] * shot_count
       combined_hi_cal_pulse = np.nanmean(combined_hi_cal_pulse[:,60:],1)\
       #                    -1e-5*shot_count
       #                   - rs_constants['comb_hi_detector_dark_count'] * shot_count
    else:
       molecular_cal_pulse = np.nanmean(molecular_cal_pulse[:,60:],1) 
       combined_hi_cal_pulse = np.nanmean(combined_hi_cal_pulse[:,60:],1)
       
    if hasattr(rs.rs_raw, 'molecular_i2a_counts'):
        if 0:
           molecular_i2a_cal_pulse = np.nanmean(molecular_i2a_cal_pulse[:,60:],1) \
                          - rs_constants['molecular_i2a_detector_dark_count'] * shot_count
        else:
           molecular_i2a_cal_pulse = np.nanmean(molecular_i2a_cal_pulse[:,60:],1) 
    if 1:
        import matplotlib.pylab as plt
        plt.figure(9000)
        plt.plot(interf_freq/1e9,molecular_cal_pulse/shot_count,'b'
                 ,interf_freq/1e9,combined_hi_cal_pulse/shot_count,'r')
        if hasattr(rs.rs_raw,'molecular_i2a_counts'):
            plt.plot(interf_freq/1e9,molecular_i2a_cal_pulse/shot_count,'m')
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('Counts/bin/laser_pulse')
        plt.grid(True)
 
        plt.figure(9001)
        plt.plot(interf_freq/1e9,molecular_cal_pulse,'b',interf_freq/1e9,combined_hi_cal_pulse,'r')
        if hasattr(rs.rs_raw,'molecular_i2a_counts'):
            plt.plot(interf_freq/1e9,molecular_i2a_cal_pulse,'m')
        plt.xlabel('Frequency (GHz)')
        plt.grid(True)
        ax = plt.gca()
        ax.set_yscale('log')

        i2t_trans = i2t_trans[i2t_freq <= max(interf_freq)]
        i2t_freq  = i2t_freq[i2t_freq <= max(interf_freq)]
        i2t_trans = i2t_trans[i2t_freq >= min(interf_freq)]
        i2t_freq = i2t_freq[i2t_freq >= min(interf_freq)]
        plt.figure(9002)
        plt.plot(interf_freq/1e9,(molecular_cal_pulse/combined_hi_cal_pulse)\
                / max(molecular_cal_pulse/combined_hi_cal_pulse),'.r')
        plt.plot(i2t_freq/1e9,i2t_trans,'k')
        if hasattr(rs.rs_raw,'molecular_i2a_counts'):
            plt.plot(interf_freq/1e9,molecular_i2a_cal_pulse/combined_hi_cal_pulse \
                     /max(molecular_i2a_cal_pulse/combined_hi_cal_pulse),'m')
        plt.xlabel('Frequency (GHz)')
        plt.grid(True)
        ax = plt.gca()
        #ax.set_yscale('log')

    freq = i2t_freq

 #compute and normalize i2 transmission
    #plot and allow adjustment of frequency offset between measurement and model
    print
    print 'Do you want to offset the zero point of the measured frequency?'
    print '   --Use this to center measured line on zero freq shift.'  
    offset_measured = 0.0
    
    #will present all results on the theoretical scan frequency axis
    while 1:

        #remove nan's
        i_freq=np.zeros_like(interf_freq)
        combined_hi = np.zeros_like(interf_freq)
        mol_i2 = np.zeros_like(interf_freq)
        scell_ratio = np.zeros_like(interf_freq)
        max_scell_ratio = np.zeros_like(interf_freq)
        if hasattr(rs.rs_raw,'molecular_i2a_cal_pulse'):
            mol_i2a = np.zeros_like(interf_freq)
        j=0    
        for i in range(len(interf_freq)):
            if (not hasattr(rs.rs_raw,'molecular_i2a_cal_pulse') \
                   and not np.isnan(combined_hi_cal_pulse[i]) \
                   and not np.isnan(molecular_cal_pulse[i])\
                   and not np.isnan(interf_freq[i])) :
                 combined_hi[j] = combined_hi_cal_pulse[i].copy()
                 mol_i2[j] = molecular_cal_pulse[i].copy()
                 scell_ratio[j] = short_cell_ratio[i].copy()
                 max_scell_ratio[j] = max_short_cell_ratio[i].copy()
                 i_freq[j] = interf_freq[i].copy()
                 j = j + 1                     
            elif (not np.isnan(molecular_i2a_cal_pulse[i])\
                   and not np.isnan(combined_hi_cal_pulse[i]) \
                   and not np.isnan(molecular_cal_pulse[i])\
                   and not np.isnan(interf_freq[i])) :
                 combined_hi[j] = combined_hi_cal_pulse[i].copy()
                 mol_i2[j] = molecular_cal_pulse[i].copy()
                 mol_i2a[j]= molecular_i2a_cal_pulse[i].copy()
                 scell_ratio[j] = short_cell_ratio[i].copy()
                 max_scell_ratio[j] = max_short_cell_ratio[i].copy()
                 i_freq[j]= interf_freq[i].copy()
                 j = j + 1

        combined_hi = combined_hi[:j-1]
        mol_i2 = mol_i2[:j-1]
        i_freq = i_freq[:j-1]
        scell_ratio = scell_ratio[:j-1]
        max_scell_ratio = max_scell_ratio[:j-1]
        if hasattr(rs.rs_raw,'molecular_i2a_cal_pulse'):
            mol_i2a = mol_i2a[:j-1]

        #make sure frequencies are monotonic
        for i in range(len(i_freq)-1):
            if  (i_freq[i+1]-i_freq[i])/delta > 0 :
                pass
            else:
                print
                print
                print
                print 'ERROR: frequency scale must be monotonic'
                print '     ---you probably requested points outside limits of scan'
                print
                print
                print jjjjjj
      
        combined_hi = np.interp(freq,i_freq+offset_measured,combined_hi)
        mol_i2 = np.interp(freq,i_freq+offset_measured,mol_i2)
        scell_ratio = np.interp(freq,i_freq+offset_measured,scell_ratio)
        max_scell_ratio = np.interp(freq,i_freq+offset_measured,max_scell_ratio)
        if hasattr(rs.rs_raw,'molecular_i2a_cal_pulse'):
                mol_i2a = np.interp(freq,i_freq+offset_measured,mol_i2a)

        
        f=plt.figure(4005)
        f.clear()
        plt.plot(freq/1e9,scell_ratio,'r',label='mean')
        plt.plot(freq/1e9,max_scell_ratio,'b',label='max')
        plt.xlabel('freq (GHz)')
        plt.ylabel('short_cell_ratio')
        ax=plt.gca()
        legend=ax.legend(loc='lower right')
        plt.grid(True)
            
        lo_freq_lmt = 6.0e9
        hi_freq_lmt = 9.0e9        
        norm = 1.0

        i2_trans = mol_i2/combined_hi
        #normalize to region between +2 and +3 GHz
        mask = np.array([ ((f >= lo_freq_lmt)
                     and (f < hi_freq_lmt)) for f in freq])
        i2_gain_ratio = nanmean(i2_trans[mask])

        #i2_trans=i2_trans/(i2_gain_ratio * norm)
        i2_trans = i2_trans/i2_gain_ratio
        
        if hasattr(rs.rs_raw,'molecular_i2a_counts'):
            #compute and normalize i2a transmission 
            i2a_trans = mol_i2a/combined_hi
            i2a_gain_ratio = nanmean(i2a_trans[mask])
            i2a_trans=i2a_trans/i2a_gain_ratio
            
        plt.figure(4006)
        plt.clf()
    
        #plot i2 transmission and theoretical values
        if hasattr(rs.rs_raw,'molecular_i2a_counts'):
            plt.plot(freq/1e9, i2_trans,'r'
                ,freq/1e9,i2t_trans,'g'
                ,freq/1e9,i2a_trans,'k')
        else:
            plt.plot(freq/1e9, i2_trans,'r'
                 ,freq/1e9, i2t_trans,'k')  

        ax=plt.gca()
        ax.set_yscale('log')
        plt.grid(True)
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('Counts')
        if hasattr(rs.rs_raw,'molecular_i2a_counts'):
            legend_str = ['i2','i2_theory','i2a']
        else:
            legend_str = ['i2','i2_theory']
        plt.legend(legend_str,loc='lower right')    
        title_str='I2 Transmission, gain_ratio= %4.3f' %(i2_gain_ratio)
        plt.title(title_str)
        plt.show(block=False)
        
       
        
        try:
           string = raw_input('offset measured frequency in GHz,(CR = exit) ? ')
           if len(string) == 0:
               break
           else:
               print 'offset by ',np.float(string),' GHz'
               #convert offset from GHz to Hz
               offset_measured=np.float(string)*1e9
               print 'offset (Hz) = ',offset_measured
            
           print 'shrink filter width for laser dither compansation ?'
           print 'shrink = 0.97 ---> 20 MHz at half power point'
           string = raw_input('fractional shrinkage = ?  ')
           if len(string)==0:
               break
           else:
               shrink = np.float(string)
               interf_freq = interf_freq * shrink
               print 'shrinking freq scale by ',shrink     
        except:
           print 'input error--try again'

            
    
#find min of transmission curve for molecular channel
    mask = np.array([ ((f >= -0.5e9) and (f < 0.5e9)) for f in freq])
    min_i2_trans = np.nanmin(i2_trans[mask])
   
   
    print 'gain_ratio= ', i2_gain_ratio
    print 'Cam = %8.2e ' %(min_i2_trans/i2_gain_ratio)
    print 'min_i2_transmission= %8.2e' %(min_i2_trans)
    print '1/i2_trans = %i' %(int(1.0/min_i2_trans))

    if hasattr(rs.rs_raw,'molecular_i2a_counts'):        
            #find min of i2a transmission curve
            min_i2a_trans = np.nanmin(i2a_trans[mask])
            print 'i2a_gain_ratio= ', i2a_gain_ratio
            print 'Cam_i2a = %8.2e ' %(min_i2a_trans/i2a_gain_ratio)
            print 'min_i2a_transmission= %8.2e' %(min_i2a_trans)
            print '1/i2a_trans = %i' %(int(1.0/min_i2a_trans))


    #normalize to max of combined channel
    
    
    p=np.polyfit(freq[np.abs(freq)<1.5e9],combined_hi[np.abs(freq)<1.5e9],2)
    center_comb_hi=np.polyval(p,freq[np.abs(freq)<1.5e9])
    norm_factor = np.max(center_comb_hi)
    print 'norm factor(max) = ',norm_factor
    norm_factor =np.polyval(p,0)
    print 'norm factor(0) =', norm_factor

    #add combined high on plot
    plt.figure(4006)
    plt.plot(freq/1e9,combined_hi/norm_factor,'r'
             ,freq[np.abs(freq)<1.5e9]/1e9,center_comb_hi/norm_factor,'k')
    plt.grid(True)
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Counts')
    plt.title('fit to combined_hi peak')
    plt.show(block=False)
    legend_str.append('chi')
    legend_str.append('chi_fit')
    plt.legend(legend_str, loc='lower right')

    #dither = np.zeros_like(freq)
    #dither[np.abs(freq)<= 20e6] = 1.0
        

    mol_i2 = mol_i2/norm_factor
    combined_hi = combined_hi/norm_factor
    if hasattr(rs.rs_raw,'molecular_i2a_counts'):
        mol_i2a=mol_i2a/norm_factor


    [cal_start,pulse,cal_end] = rs_constants['apd_pulse_timing'] 

        #write i2 scan file
    #get path to directory
    dir_path=calibration_path_for(instrument,rs.rs_raw.times[0],process_defaults)
    #define start time for use in filename
    str = start_time.strftime("%Y%m%dT%H%M")
    filename=os.path.join(dir_path,'i2-scan-' +str+'.cal')
    fileid=open(filename,'w')
    os.chmod(filename,0664)
    print >>fileid, '# Calibration scan as function of frequency offset from I2 line'
    print >>fileid, '# calibration scan data aquired on ' + start_time.strftime("%d-%b-%y")\
                 +' at '+ start_time.strftime("%H:%M") +' UT'
    print >>fileid, '# file created on '+ datetime.now().strftime("%d-%b-%y") +' at ' \
                 +datetime.now().strftime("%H:%M") + ' UT'
    print >>fileid, '# t_begin_cal_pulse= %4.2e ;  start time of cal pulse (sec).' %(cal_start)
    print >>fileid, '# t_end_cal_pulse=   %4.2e ;  end time of cal pulse (sec).' %(cal_end)
    print >>fileid, '# pulse_durration=   %4.2e ;  laser pulse durration (sec).' \
          %(rs_constants['binwidth'])
    print >>fileid, '# ratio of mol to combined channel gain = %6.3f ' %(i2_gain_ratio)
    print >>fileid, '# Cam = %8.3e ' %(min_i2_trans/i2_gain_ratio) 
    print >>fileid, '# Min iodine transmission =   %8.1e,  1/(min_trans) = %8.1e'\
                        %(min_i2_trans, 1/min_i2_trans)
    if hasattr(rs.rs_raw,'molecular_i2a_counts'):
        print >>fileid, '# ratio of mol_i2a to combined channel gain = %5.1f ' %(i2a_gain_ratio)
        print >>fileid, '# Cam_i2a = %8.3e ' %(min_i2a_trans/i2a_gain_ratio) 
        print >>fileid, '# Min iodine_argon transmission =   %8.3e,  1/(min_trans) = %8.3e'\
                        %(min_i2a_trans, 1/min_i2a_trans) 
    print >>fileid, '# '
    if hasattr(rs.rs_raw,'molecular_i2a_counts'):
        print >>fileid,\
           '#freq(GHz)    combined    molecular   i2_theory   i2_trans    i2a_trans    molecular_i2a '
        for i in range(np.size(freq)):
            print >>fileid, '%8.5f    %9.5f   %9.5f   %9.5f   %9.5f   %9.5f   %9.5f'\
            %(freq[i]/1e9 ,combined_hi[i], mol_i2[i]
                 ,i2t_trans[i],i2_trans[i],i2a_trans[i], mol_i2a[i])
    else:    
        print >>fileid, '#freq(GHz)  combined  molecular i2_measured i2_theory'
        for i in range(np.size(freq)):
            print >>fileid, '%8.5f   %9.5f  %9.5f  %9.5f  %9.5f'\
              %(freq[i]/1e9 ,combined_hi[i], mol_i2[i]
                 ,i2t_trans[i],i2_trans[i])
    fileid.close()
    
    print '\nnew i2 scan file=\n '+filename

    
def make_i2_scan_file(instrument,rs,rs_constants,process_defaults=None):
    """make_i2_scan_file(instrument,rs,rs_constants)
       make an i2_scan file from a calibration scan and store it in the
       calibration directory--instrument temperatures must be stable
       during measurement to minimize interferometer drift between the
       wide scan and the narrow scan"""

    #read theoretical i2 transmission file
    try:
         [header,i2t]=ru.readascii(locate_file('dfb_i2_spectra.txt'))
         dfb_to_i2_scale = rs_constants['i2_absorption_scale_factor']
         i2t_trans=10**(dfb_to_i2_scale*i2t[:,1])
         i2t_freq = i2t[:,0]*1e9
    except IOError: #file either not found by locate, or error occurred 
         [header,i2t]=ru.readascii(locate_file('I2cell_272_31_extended_freq.txt'))
         i2t_freq  = i2t[:,0]*1e9
         i2t_trans = i2t[:,2]

    try:
        use_superseed_freq = True
        print
        print 'Using superseed temperature to compute frequency shift, %g Ghz/K'\
                                   %(rs_constants['seedlaser_temp_to_freq'])
        print
        if rs.rs_raw.superseedlasercontrollog.shape[1]==11:
          index=7
        elif rs.rs_raw.superseedlasercontrollog.shape[1]==10:
          index=2
        else:
          print
          print
          raise NotImplemented("Can't get frequency from superseed board. %i element not supported" % (rs.rs_raw.superseedlasercontrollog.shape[1]))
          
        interf_freq = (rs.rs_raw.superseedlasercontrollog[:,index] \
                 - rs.rs_raw.superseedlasercontrollog[0,index])\
                 * rs_constants['seedlaser_temp_to_freq'] * 1e9
    except:  #use interferometer measured frequency
        print
        print 'no superseeder--using interferometer freq scale'
        print
        use_superseed_freq = False
        interf_freq = rs.rs_raw.interf_freq.copy()
    times = rs.rs_raw.times.copy()
    molecular_cal_pulse = rs.rs_raw.molecular_cal_pulse.copy()
    combined_hi_cal_pulse = rs.rs_raw.combined_hi_cal_pulse.copy()
    #assume system is locked at start of record--set this frequency to zero
    #nterf_freq=interf_freq-np.nanmean(interf_freq[:5])
    
    if hasattr(rs.rs_raw, 'molecular_i2a_counts'):
        molecular_i2a_cal_pulse = rs.rs_raw.molecular_i2a_cal_pulse.copy()

    #select time interval
    date_str = datetime.strftime(times[0], '%d-%b-%y ')
    start_time_str = raw_input("Start time=? H:M:S or H:M ")
    start_time_str = date_str + start_time_str
    if start_time_str.count(':') == 2:
      start_time = datetime.strptime(start_time_str, '%d-%b-%y %H:%M:%S')
    else:
      start_time = datetime.strptime(start_time_str, '%d-%b-%y %H:%M')
    end_time_str = raw_input("End time H:M:S or H:M ")
    end_time_str = date_str + end_time_str
    if end_time_str.count(':') == 2:
       end_time = datetime.strptime(end_time_str, '%d-%b-%y %H:%M:%S')
    else:
       end_time = datetime.strptime(end_time_str, '%d-%b-%y %H:%M')


    #restrict times via user input
    time_mask = np.array([ ((x >= start_time) and (x<end_time)) 
                             for x in times[:] ])
    interf_freq = interf_freq[time_mask]
    molecular_cal_pulse = molecular_cal_pulse[time_mask]
    combined_hi_cal_pulse = combined_hi_cal_pulse[time_mask]
    times = times[time_mask]
    
    if hasattr(rs.rs_raw,'molecular_i2a_counts'):
        molecular_i2a_cal_pulse = molecular_i2a_cal_pulse[time_mask]

    np.set_printoptions(threshold=np.NaN)

    #select only those points where frequency is being scanned
    i2_scan = (rs.rs_raw.op_mode[time_mask] & 1) == 1

    #select bins when frequency is scanning and only one of the calibration filters is inserted
    filter_1 = (rs.rs_raw.op_mode[time_mask] & 128) > 0
    filter_2  = (rs.rs_raw.op_mode[time_mask] & 32) > 0
    filter_1k = np.logical_and(np.logical_xor(filter_2 , filter_1) , i2_scan)
    #select center 5 GHz of scan
    filter_1k = np.logical_and(filter_1k, np.abs(interf_freq) < 2.5e9)
      
    #filter factor when only one of the calibration filters is inserted
    filter_1k_factor = 10.0**rs_constants['calibration_nd_filters'][0]


    #select bins when frequency is scanning and both filters are inserted.
    filter_1000k = np.logical_and(filter_1,i2_scan)
    
    plt.figure(4000)
    if hasattr(rs.rs_raw,'molecular_i2a_cal_pulse'):
        plt.plot(interf_freq[filter_1k]/1e9
             ,molecular_cal_pulse[filter_1k]/filter_1k_factor,'m'
             ,interf_freq[filter_1000k]/1e9
             ,molecular_cal_pulse[filter_1000k],'b'   
             ,interf_freq[filter_1000k]/1e9
             ,combined_hi_cal_pulse[filter_1000k],'r'
             ,interf_freq[filter_1k]/1e9
             ,molecular_i2a_cal_pulse[filter_1k]/filter_1k_factor,'k'
             ,interf_freq[filter_1000k]/1e9
             ,molecular_i2a_cal_pulse[filter_1000k],'c')
        plt.legend(['mi2_1k','mi2_1000k','c_1000k','mi2a_1k','mi2a_1000k'])
    else:
        plt.plot(interf_freq[filter_1k]/1e9
             ,molecular_cal_pulse[filter_1k]/1000,'m'
             ,interf_freq[filter_1000k]/1e9
             ,molecular_cal_pulse[filter_1000k] ,'b'
             ,interf_freq[filter_1000k]/1e9
             ,combined_hi_cal_pulse[filter_1000k],'r')
        plt.legend(('mi2_1k','mi2_1000k','c_1000k','mi2a_1k','mi2a_1000k'),'lower left')
    plt.grid(True)
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Counts')
    plt.title('Raw scan--linear '+start_time_str)

    plt.show(block=False)

    #logrithmic plot
    plt.figure(4001)
    if hasattr(rs.rs_raw,'molecular_i2a_cal_pulse'):
        plt.plot(interf_freq[filter_1k]/1e9
             ,molecular_cal_pulse[filter_1k]/filter_1k_factor,'g'
             ,interf_freq[filter_1000k]/1e9
             ,molecular_cal_pulse[filter_1000k],'b'   
             ,interf_freq[filter_1000k]/1e9
             ,combined_hi_cal_pulse[filter_1000k],'r'
             ,interf_freq[filter_1k]/1e9
                 ,molecular_i2a_cal_pulse[filter_1k]/filter_1k_factor,'k'
             ,interf_freq[filter_1000k]/1e9
             ,molecular_i2a_cal_pulse[filter_1000k],'c')
        plt.legend(('mi2_1k','mi2_1000k','c_1000k','mi2a_1k','mi2a_1000k'),'lower left')
    else:
        plt.plot(interf_freq[filter_1k]/1e9
             ,molecular_cal_pulse[filter_1k]/filter_1k_factor,'g'
             ,interf_freq[filter_1000k]/1e9
             ,molecular_cal_pulse[filter_1000k],'b'
             ,interf_freq[filter_1000k]/1e9
             ,combined_hi_cal_pulse[filter_1000k],'r')
        plt.legend(('mi2_1k','mi2_1000k','c_1000k'),'lower left')
    ax=plt.gca()
    try:
        ax.set_yscale('log')
    except:
        print
        print 'unable to do log conversion on cal pulse plot'
        print 'possibly due to empty array--figure 4001'
        print 'filter_1000k shape',filter_1000k.shape
        print 'filter_1k.shape',filter_1k.shape
        ax.set_yscale('linear')
        print
    plt.grid(True)
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Counts')
    plt.title('Raw scan-log '+start_time_str)
    plt.show(block=False)
    

    mol_1k = molecular_cal_pulse[filter_1k]
    mol_1000k = molecular_cal_pulse[filter_1000k]
    freq_1k = interf_freq[filter_1k]
    freq_1000k =interf_freq[filter_1000k]

    #adjust combined channel offset if needed
    chi_offset = 0
    combined_hi_max = np.nanmax(combined_hi_cal_pulse[filter_1000k])
    ratio_max = np.nanmax(mol_1000k/combined_hi_cal_pulse[filter_1000k])
    print
    print 'max_min comb_hi ', combined_hi_max \
                               , np.nanmin(combined_hi_cal_pulse[filter_1000k])
    while 1:
       plt.ion()
       f=plt.figure(5000)
       measured_trans = mol_1000k/(combined_hi_cal_pulse[filter_1000k]-chi_offset)
       mask = np.array([ ((f >= 2e9) and (f < 3e9)) for f in freq_1000k])
       measured_trans = measured_trans/ np.nanmean(measured_trans[mask])
       plt.plot(freq_1000k,measured_trans,'r'
                ,i2t_freq,i2t_trans,'k')
       plt.ylabel('transmission_i2')
       plt.xlabel('frequency offset (GHz)')
       plt.grid(True)
       plt.ylim((0,1.5))
       plt.xlim((np.nanmin(freq_1000k),np.nanmax(freq_1000k)))
       plt.show(block=False)
       offset_string = raw_input('Input combined_hi offset if needed, <CR> = no change    ')
       if len(offset_string) == 0:
            combined_hi_cal_pulse = combined_hi_cal_pulse - chi_offset
            break
       else:
         chi_offset = float(offset_string)
         #plt.show(block=False)
    #threshhold used to splice wide and narrow scans
    #this was set to 12---need to understand why the change!

    splice_level = 0.02
    print
    print 'splice level used to merge wide and narrow scans = ',splice_level
    print
    print 'Do you want to change the splice level used to form composite scan?'
    print 'when the number of photons in the narrow scan exceeds this number'
    print 'the composite output scan will switch from the narrow to wide scan'
    splice_str = raw_input('<CR> = no change)   ')
    if len(splice_str) == 0:
        pass
    else:
        splice_level = float(splice_str)
        print 'new splice level  = ',splice_level
    
    #select only those points with less than splice_level photons in cal pulse
    #from the narrow scan       
    mask_1k= np.array([((s < splice_level * filter_1k_factor)) for s in mol_1k[:] ])
    if len(mask_1k) == 0:   #not mask.any():
        print
        print 'ERROR---Frequency range must contain wide and narrow scans'
        print

    freq_1k_i2 = freq_1k[mask_1k]
    mol_1k_i2  = mol_1k[mask_1k]
    
     
    #select all other frequencys from the wide scan 
    maxfreq= np.nanmax(freq_1k_i2)
    minfreq=  np.nanmin(freq_1k_i2)
    
    mask = np.array([((f  > maxfreq              
             or f <= minfreq))for f in freq_1000k[:]])

    mol_1000k_i2 = mol_1000k[mask]
    freq_1000k_i2 = freq_1000k[mask]
    
    #concatenate mol i2 filter signals--frequency will be out of order 
    mol_i2 = np.concatenate((mol_1000k_i2,mol_1k_i2/filter_1k_factor),0)
    freq_i2 = np.concatenate((freq_1000k_i2,freq_1k_i2),0)
    
    #mol_i2 = np.concatenate((mol_1000k,mol_1k/filter_1k_factor),0)
    #freq_i2 = np.concatenate((freq_1000k,freq_1k),0)
    

    #repeat for i2a signals if they are present
    if hasattr(rs.rs_raw,'molecular_i2a_cal_pulse'):
        mol_1k_i2a = molecular_i2a_cal_pulse[filter_1k]
        mol_1000k_i2a = molecular_i2a_cal_pulse[filter_1000k]
        mask= np.array([((s < splice_level * filter_1k_factor)) for s in mol_1k_i2a[:] ])
        freq_1k_i2a = freq_1k[mask_1k]
        mol_1k_i2a  = mol_1k_i2a[mask_1k]
        mask = np.array([((f > np.nanmax(freq_1k_i2a)                 
             or f <= np.nanmin(freq_1k_i2a)))for f in freq_1000k[:]])
        mol_1000k_i2a = mol_1000k_i2a[mask]
        freq_1000k_i2a = freq_1000k[mask]
        
        #concatenate mol i2a filter signals
        mol_i2a = np.concatenate((mol_1000k_i2a,mol_1k_i2a/filter_1k_factor),0)
        freq_i2a = np.concatenate((freq_1000k_i2a,freq_1k_i2a),0)


   


    #rearange mol in freq sorted order
    indices =np.argsort(freq_i2)
    freq_mol_i2 = freq_i2[indices]
    molecular_i2 = mol_i2[indices]
    if hasattr(rs.rs_raw,'molecular_i2a_cal_pulse'):
       indices = np.argsort(freq_i2a)
       molecular_i2a = mol_i2a[indices]  
       freq_mol_i2a = freq_i2a[indices]

    #rearange combined hi in freq sorted order 
    indices = np.argsort(freq_1000k)
    combined_hi_cal_pulse = combined_hi_cal_pulse[indices].copy()
    freq_1000k = freq_1000k[indices].copy()
    """
    #read theoretical i2 transmission file
    try:
         [header,i2t]=ru.readascii(locate_file('dfb_i2_spectra.txt'))
         dfb_to_i2_scale = rs_constants['i2_absorption_scale_factor']
         i2t_trans=10**(dfb_to_i2_scale*i2t[:,1])
         i2t_freq = i2t[:,0]*1e9         
    except IOError: #file either not found by locate, or error occurred 
         [header,i2t]=ru.readascii(locate_file('I2cell_272_31_extended_freq.txt'))
         i2t_freq  = i2t[:,0]*1e9
         i2t_trans = i2t[:,2]
    """

    #normalize theoretical spectrum to 1--ie eliminate continuum absorption
    lo_freq_lmt = 2e9
    hi_freq_lmt = 4e9
    mask = np.array([ ((f >= lo_freq_lmt) and (f < hi_freq_lmt)) for f in i2t_freq])
    print 'normalizing i2 transmission between ', lo_freq_lmt/1e9, ' and ',hi_freq_lmt/1e9, ' GHz'
    norm = nanmean(i2t_trans[mask])
    i2t_trans = i2t_trans/norm

    
    plt.figure(4002)
    if hasattr(rs.rs_raw,'molecular_i2a_cal_pulse'):
        plt.plot(freq_mol_i2/1e9
             ,molecular_i2,'b'
             ,freq_1k_i2/1e9
             ,mol_1k_i2/filter_1k_factor,'g'    
             ,freq_1000k/1e9
             ,combined_hi_cal_pulse,'r'
             ,freq_mol_i2a/1e9
             ,molecular_i2a,'k')
          
        plt.legend(('mol','mol_1k','comb_hi','mol_i2a'),'lower left')
    else:
        plt.plot(freq_mol_i2/1e9,molecular_i2,'b'
                 ,freq_1k_i2/1e9,mol_1k_i2/filter_1k_factor,'g'    
             ,freq_1000k/1e9,combined_hi_cal_pulse,'r')
        plt.legend(('mol','combined'),'lower left')
        
    ax=plt.gca()
    ax.set_yscale('log')
    plt.grid(True)
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('counts')
    plt.title('composite scan ' +start_time_str)



    plt.figure(4003)
    plt.plot(freq_mol_i2/1e9,molecular_i2,'b'
             ,freq_1000k/1e9,combined_hi_cal_pulse,'r')
    plt.legend(('mol','combined'),'lower left')

    ax=plt.gca()
    ax.set_yscale('log')
    plt.grid(True)
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('counts')
    plt.title('composite I2 scan' +start_time_str)


    if hasattr(rs.rs_raw,'molecular_i2a_cal_pulse'):
        plt.figure(4004)
        plt.plot(freq_mol_i2a/1e9,molecular_i2a,'b'
             ,freq_1000k/1e9,combined_hi_cal_pulse,'r')
        plt.legend(('mol_i2a','combined'),'lower left')

        ax=plt.gca()
        ax.set_yscale('log')
        plt.grid(True)
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('counts')
        plt.title('composite I2A scan' +start_time_str)
    plt.show(block=False)



    
    #move set point frequency to another line?
    print
    print 'Do you want to move to a different freq set point?'
    print '    --Use this only if you want to move set point away from line 1109'
    offset_setpoint = rs_constants['lock_point_freq_offset']
    #will present all results on the theoretical scan frequency axis
    freq = i2t_freq
    while 1:
    
       
        plt.ion()
        plt.figure(4005)
        plt.clf()
        plt.plot(freq/1e9, i2t_trans,'k')  

        ax=plt.gca()
        ax.set_yscale('log')
        plt.grid(True)
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('Counts')
        
        title_str='I2 Transmission,lockpoint offset= %4.3f GHz' %(offset_setpoint)
        plt.title(title_str)
        plt.show(block=False)
        
        try:
           string = raw_input('Lock point frequency offset (GHz), (CR = exit) ? ')
           if len(string) == 0:
               break
           else:
               #convert offset from GHz to Hz
               offset_setpoint = np.float(string)*1e9
               freq=freq+offset_setpoint
        except:
               print 'input error--try again'
        

    


           

    #compute and normalize i2 transmission
    #plot and allow adjustment of frequency offset between measurement and model
    print
    print 'Do you want to offset the zero point of the measured frequency?'
    print '   --Use this to center measured line on zero freq shift.'  
    offset_measured = 0.0
    
    #will present all results on the theoretical scan frequency axis
    while 1:
    
        #interpolate to frequencys in theoretical i2 spectrum file
        combined_hi = np.interp(freq,freq_1000k+offset_measured,combined_hi_cal_pulse)
        mol_i2 = np.interp(freq,freq_mol_i2+offset_measured,molecular_i2)
       
        if hasattr(rs.rs_raw,'molecular_i2a_cal_pulse'):
            mol_i2a = np.interp(freq,freq_mol_i2a+offset_measured,molecular_i2a) 
        

        i2_trans = mol_i2/combined_hi
        #normalize to region between +2 and +3 GHz
        mask = np.array([ ((f >= lo_freq_lmt)
                     and (f < hi_freq_lmt)) for f in freq])
        i2_gain_ratio = nanmean(i2_trans[mask])
        i2_trans=i2_trans/(i2_gain_ratio * norm)
        if hasattr(rs.rs_raw,'molecular_i2a_counts'):
            #compute and normalize i2a transmission 
            i2a_trans = mol_i2a/combined_hi
            i2a_gain_ratio = nanmean(i2a_trans[mask])
            i2a_trans=i2a_trans/i2a_gain_ratio
       
        plt.ion()
        plt.figure(4006)
        plt.clf()
    
        #plot i2 transmission and theoretical values
        if hasattr(rs.rs_raw,'molecular_i2a_counts'):
            plt.plot(freq/1e9, i2_trans,'r'
                ,freq/1e9,i2t_trans,'g'
                ,freq/1e9,i2a_trans,'k')
        else:
            plt.plot(freq/1e9, i2_trans,'r'
                 ,freq/1e9, i2t_trans,'k')  

        ax=plt.gca()
        ax.set_yscale('log')
        plt.grid(True)
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('Counts')
        if hasattr(rs.rs_raw,'molecular_i2a_counts'):
            legend_str = ['i2','i2_theory','i2a']
        else:
            legend_str = ['i2','i2_theory']
        plt.legend(legend_str,'lower right')    
        title_str='I2 Transmission, gain_ratio= %4.3f' %(i2_gain_ratio)
        plt.title(title_str)
        plt.show(block=False)
        
        try:
           string = raw_input('offset measured frequency in GHz,(CR = exit) ? ')
           if len(string) == 0:
               break
           else:
               #convert offset from GHz to Hz
               offset_measured=np.float(string)*1e9
        except:
           print 'input error--try again'
       

    #find min of transmission curve for molecular channel
    mask = np.array([ ((f >= -0.5e9) and (f < 0.5e9)) for f in freq])
    min_i2_trans = np.nanmin(i2_trans[mask])
   
   
    print 'gain_ratio= ', i2_gain_ratio
    print 'Cam = %8.2e ' %(min_i2_trans/i2_gain_ratio)
    print 'min_i2_transmission= %8.2e' %(min_i2_trans)
    print '1/i2_trans = %i' %(int(1.0/min_i2_trans))

    if hasattr(rs.rs_raw,'molecular_i2a_counts'):        
            #find min of i2a transmission curve
            min_i2a_trans = np.nanmin(i2a_trans[mask])
            print 'i2a_gain_ratio= ', i2a_gain_ratio
            print 'Cam_i2a = %8.2e ' %(min_i2a_trans/i2a_gain_ratio)
            print 'min_i2a_transmission= %8.2e' %(min_i2a_trans)
            print '1/i2a_trans = %i' %(int(1.0/min_i2a_trans))


    #normalize to max of combined channel
    
    
    p=np.polyfit(freq[np.abs(freq)<1.5e9],combined_hi[np.abs(freq)<1.5e9],2)
    center_comb_hi=np.polyval(p,freq[np.abs(freq)<1.5e9])
    norm_factor = np.max(center_comb_hi)
    print 'norm factor(max) = ',norm_factor
    norm_factor =np.polyval(p,0)
    print 'norm factor(0) =', norm_factor

    #add combined high on plot
    plt.figure(4006)
    plt.plot(freq/1e9,combined_hi,'r',freq[np.abs(freq)<1.5e9]/1e9,center_comb_hi,'k')
    plt.grid(True)
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Counts')
    plt.title('fit to combined_hi peak')
    plt.show(block=False)
    legend_str.append('chi')
    legend_str.append('chi_fit')
    plt.legend(legend_str,'lower right')
    mol_i2= mol_i2/norm_factor
    combined_hi = combined_hi/norm_factor
    if hasattr(rs.rs_raw,'molecular_i2a_counts'):
        mol_i2a=mol_i2a/norm_factor

    
           
    [cal_start,pulse,cal_end] = rs_constants['apd_pulse_timing'] 
    
    #write i2 scan file
    #get path to directory
    dir_path=calibration_path_for(instrument,rs.rs_raw.times[0],process_defaults)
    #define start time for use in filename
    str = start_time.strftime("%Y%m%dT%H%M")
    filename=os.path.join(dir_path,'i2-scan-' +str+'.cal')
    fileid=open(filename,'w')
    os.chmod(filename,0664)
    print >>fileid, '# Calibration scan as function of frequency offset from I2 line'
    print >>fileid, '# calibration scan data aquired on ' + start_time.strftime("%d-%b-%y")\
                 +' at '+ start_time.strftime("%H:%M") +' UT'
    print >>fileid, '# file created on '+ datetime.now().strftime("%d-%b-%y") +' at ' \
                 +datetime.now().strftime("%H:%M") + ' UT'
    print >>fileid, '# t_begin_cal_pulse= %4.2e ;  start time of cal pulse (sec).' %(cal_start)
    print >>fileid, '# t_end_cal_pulse=   %4.2e ;  end time of cal pulse (sec).' %(cal_end)
    print >>fileid, '# pulse_durration=   %4.2e ;  laser pulse durration (sec).' \
          %(rs_constants['binwidth'])
    print >>fileid, '# ratio of mol to combined channel gain = %6.3f ' %(i2_gain_ratio)
    print >>fileid, '# Cam = %8.3e ' %(min_i2_trans/i2_gain_ratio) 
    print >>fileid, '# Min iodine transmission =   %8.1e,  1/(min_trans) = %8.1e'\
                        %(min_i2_trans, 1/min_i2_trans)
    if hasattr(rs.rs_raw,'molecular_i2a_counts'):
        print >>fileid, '# ratio of mol_i2a to combined channel gain = %5.1f ' %(i2a_gain_ratio)
        print >>fileid, '# Cam_i2a = %8.3e ' %(min_i2a_trans/i2a_gain_ratio) 
        print >>fileid, '# Min iodine_argon transmission =   %8.3e,  1/(min_trans) = %8.3e'\
                        %(min_i2a_trans, 1/min_i2a_trans) 
    print >>fileid, '# '
    if hasattr(rs.rs_raw,'molecular_i2a_counts'):
        print >>fileid,\
           '#freq(GHz)    combined    molecular   i2_theory   i2_trans    i2a_trans    molecular_i2a '
        for i in range(np.size(freq)):
            print >>fileid, '%8.5f    %9.5f   %9.5f   %9.5f   %9.5f   %9.5f   %9.5f'\
	    %(freq[i]/1e9 ,combined_hi[i], mol_i2[i]
                 ,i2t_trans[i],i2_trans[i],i2a_trans[i], mol_i2a[i])
    else:    
        print >>fileid, '#freq(GHz)  combined  molecular i2_measured i2_theory'
        for i in range(np.size(freq)):
            print >>fileid, '%8.5f   %9.5f  %9.5f  %9.5f  %9.5f'\
	      %(freq[i]/1e9 ,combined_hi[i], mol_i2[i]
                 ,i2t_trans[i],i2_trans[i])
    fileid.close()
    
    print '\nnew i2 scan file=\n '+filename

def make_i2_scan_from_i2_spectrum_file(instrument,rs,rs_constants,process_defaults=None):
    """make_i2_scan_file_from_dfb_scan(instrument,rs,rs_constants)
       make an i2_scan file from the combination of a calibration scan
       and a distributed feedback diode scan of the I2 spectrum.
       Store the result in the calibration directory--instrument
       temperatures must be stable during measurement to minimize
       interferometer drift between the wide scan and the narrow scan"""

    def select_times(times):
        #select time interval
        date_str = datetime.strftime(times[0], '%d-%b-%y ')
        start_time_str = raw_input("Start time=? H:M:S or H:M ")
        start_time_str = date_str + start_time_str
        if start_time_str.count(':') == 2:
            start_time = datetime.strptime(start_time_str, '%d-%b-%y %H:%M:%S')
        else:
            start_time = datetime.strptime(start_time_str, '%d-%b-%y %H:%M')
        end_time_str = raw_input("End time H:M:S or H:M ")
        end_time_str = date_str + end_time_str
        if end_time_str.count(':') == 2:
            end_time = datetime.strptime(end_time_str, '%d-%b-%y %H:%M:%S')
        else:
            end_time = datetime.strptime(end_time_str, '%d-%b-%y %H:%M')
        return date_str,start_time,start_time_str,end_time,end_time_str

    def plot_raw_counts(rs,start_time_str):
        #plot raw counts per laser pulse
        plt.figure(3999)
        if hasattr(rs,'molecular_i2a_cal_pulse'):
            plt.plot(rs.interf_freq/1e9
                 ,rs.molecular_cal_pulse/rs.seeded_shots,'b'   
                 ,rs.interf_freq/1e9
                 ,rs.combined_hi_cal_pulse/rs.seeded_shots,'r'
                 ,rs.interf_freq/1e9
                 ,rs.molecular_i2a_cal_pulse/rs.seeded_shots,'k')
            plt.legend(('i2','comb','i2a'),'lower left')
        else:
            plt.plot(rs.interf_freq/1e9
               ,rs.molecular_cal_pulse/rs.seeded_shots,'b'
               ,rs.interf_freq/1e9
               ,rs.combined_hi_cal_pulse/rs.seeded_shots,'r')
            plt.legend(('i2''comb'),'lower left')
        plt.grid(True)
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('Counts/laser_pulse')
        plt.title('Raw scan '+start_time_str)

        
    def make_manual_scan(rs,start_time,end_time,start_time_str):
        interf_freq = rs.interf_freq.copy()
  
        #energy normalization
        molecular_cal_pulse = rs.molecular_cal_pulse/rs.transmitted_energy
        combined_hi_cal_pulse = rs.combined_hi_cal_pulse/rs.transmitted_energy
    
        #assume system is locked at start of record--set this frequency to zero
        interf_freq=interf_freq-np.mean(interf_freq[:5])
    
        if hasattr(rs, 'molecular_i2a_counts'):
            molecular_i2a_cal_pulse = rs.molecular_i2a_cal_pulse/rs.transmitted_energy

        time_mask = np.array([ ((x >= start_time) and (x<end_time)) 
                             for x in rs.times[:] ])
        interf_freq = interf_freq[time_mask]
        molecular_cal_pulse = molecular_cal_pulse[time_mask]
        combined_hi_cal_pulse = combined_hi_cal_pulse[time_mask]
        if hasattr(rs,'molecular_i2a_counts'):
            molecular_i2a_cal_pulse = molecular_i2a_cal_pulse[time_mask]

        np.set_printoptions(threshold=np.NaN)  
        #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



        
        #rearange mol in freq sorted order
        indices =np.argsort(interf_freq)
        freq_mol_i2 = interf_freq[indices].copy()
        molecular_i2 = molecular_cal_pulse[indices].copy()
        if hasattr(rs,'molecular_i2a_cal_pulse'):
           molecular_i2a = molecular_i2a_cal_pulse[indices]  
           freq_mol_i2a = freq_mol_i2.copy()
        else:
            molecular_i2a = []
            freq_mol_i2a  = []
        
        #rearange combined hi in freq sorted order
        combined_hi = combined_hi_cal_pulse[indices].copy()
        combined_freq = freq_mol_i2.copy()
        """
        combined_hi = combined_hi_cal_pulse.copy()
        combined_freq = freq_1000k.copy()
        """
        plt.figure(4002)
        if hasattr(rs,'molecular_i2a_cal_pulse'):
            plt.plot(freq_mol_i2/1e9
               ,molecular_i2,'b'
               ,combined_freq/1e9
               ,combined_hi,'r'
               ,freq_mol_i2a/1e9
               ,molecular_i2a,'k')

            plt.legend(('mol','comb_hi','mol_i2a'),'lower left')
        else:
            plt.plot(freq_mol_i2/1e9,molecular_i2,'b'
             ,combined_freq/1e9,combined_hi,'r')
            plt.legend(('mol','combined'),'lower left')

        ax=plt.gca()
        ax.set_yscale('log')
        plt.grid(True)
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('counts')
        plt.title('composite scan ' +start_time_str)

        plt.figure(4003)
        plt.plot(freq_mol_i2/1e9,molecular_i2,'b'
             ,combined_freq/1e9,combined_hi,'r')
        plt.legend(('mol','combined'),'lower left')
        ax=plt.gca()
        ax.set_yscale('log')
        plt.grid(True)
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('counts')
        plt.title('composite I2 scan' +start_time_str)


        if hasattr(rs,'molecular_i2a_cal_pulse'):
            plt.figure(4004)
            plt.plot(freq_mol_i2a/1e9,molecular_i2a,'b'
                ,combined_freq/1e9,combined_hi,'r')
            plt.legend(('mol_i2a','combined'),'lower left')
            ax=plt.gca()
            ax.set_yscale('log')
            plt.grid(True)
            plt.xlabel('Frequency (GHz)')
            plt.ylabel('counts')
            plt.title('composite I2A scan' +start_time_str)
            plt.show(block=False)
        return combined_freq,combined_hi,freq_mol_i2,molecular_i2,freq_mol_i2a,molecular_i2a
    
    def make_composite_scans(rs,start_time,end_time,start_time_str ):
        """ combines wide and narrow scans to increase the dynamic range
            produces seperate combined, i2 and i2a scans with a different
            frequency scale for the i2a scan. Frequencies are sorted in
            ascending order but are not equally spaced."""

        interf_freq = rs.interf_freq.copy()
  
        #energy normalization
        molecular_cal_pulse = rs.molecular_cal_pulse/rs.transmitted_energy
        combined_hi_cal_pulse = rs.combined_hi_cal_pulse/rs.transmitted_energy
    
        #assume system is locked at start of record--set this frequency to zero
        interf_freq=interf_freq-np.mean(interf_freq[:5])
    
        if hasattr(rs, 'molecular_i2a_counts'):
            molecular_i2a_cal_pulse = rs.molecular_i2a_cal_pulse/rs.transmitted_energy

        time_mask = np.array([ ((x >= start_time) and (x<end_time)) 
                             for x in rs.times[:] ])
        interf_freq = interf_freq[time_mask]
        molecular_cal_pulse = molecular_cal_pulse[time_mask]
        combined_hi_cal_pulse = combined_hi_cal_pulse[time_mask]
        if hasattr(rs,'molecular_i2a_counts'):
            molecular_i2a_cal_pulse = molecular_i2a_cal_pulse[time_mask]

        np.set_printoptions(threshold=np.NaN)

        #select only those points where frequency is being scanned
        i2_scan = (rs.op_mode[time_mask] & 1) == 1

        #select bins when frequency is scanning and only one of the 1/1000 filters is inserted
        filter_1 = (rs.op_mode[time_mask] & 128) > 0
        filter_2  = (rs.op_mode[time_mask] & 32) > 0
        filter_1k  = (filter_2 ^ filter_1) * i2_scan

        #select bins when frequency is scanning and both 1/1000 filters are inserted.
        filter_1000k = filter_1*i2_scan
  
        plt.figure(4000)
        if hasattr(rs,'molecular_i2a_cal_pulse'):
            plt.plot(interf_freq[filter_1k]/1e9
             ,molecular_cal_pulse[filter_1k]/1000,'m'
             ,interf_freq[filter_1000k]/1e9
             ,molecular_cal_pulse[filter_1000k],'b'   
             ,interf_freq[filter_1000k]/1e9
             ,combined_hi_cal_pulse[filter_1000k],'r'
             ,interf_freq[filter_1k]/1e9
             ,molecular_i2a_cal_pulse[filter_1k]/1000,'k'
             ,interf_freq[filter_1000k]/1e9
             ,molecular_i2a_cal_pulse[filter_1000k],'c')
            plt.legend(['mi2_1k','mi2_1000k','c_1000k','mi2a_1k','mi2a_1000k'])
        else:
            plt.plot(interf_freq[filter_1k]/1e9
             ,molecular_cal_pulse[filter_1k]/1000,'m'
             ,interf_freq[filter_1000k]/1e9
             ,molecular_cal_pulse[filter_1000k] ,'b'
             ,interf_freq[filter_1000k]/1e9
             ,combined_hi_cal_pulse[filter_1000k],'r')
            plt.legend(('mi2_1k','mi2_1000k','c_1000k','mi2a_1k','mi2a_1000k'),'lower left')
        plt.grid(True)
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('Counts')
        plt.title('Energy normalized scan--linear '+start_time_str)

        plt.show(block=False)

        #logrithmic plot
        plt.figure(4001)
        if hasattr(rs,'molecular_i2a_cal_pulse'):
            plt.plot(interf_freq[filter_1k]/1e9
              ,molecular_cal_pulse[filter_1k]/1000,'g'
              ,interf_freq[filter_1000k]/1e9
              ,molecular_cal_pulse[filter_1000k],'b'   
              ,interf_freq[filter_1000k]/1e9
              ,combined_hi_cal_pulse[filter_1000k],'r'
              ,interf_freq[filter_1k]/1e9
              ,molecular_i2a_cal_pulse[filter_1k]/1000,'k'
              ,interf_freq[filter_1000k]/1e9
              ,molecular_i2a_cal_pulse[filter_1000k],'c')
            plt.legend(('mi2_1k','mi2_1000k','c_1000k','mi2a_1k','mi2a_1000k'),'lower left')
        else:
            plt.plot(interf_freq[filter_1k]/1e9
              ,molecular_cal_pulse[filter_1k]/1000,'g'
              ,interf_freq[filter_1000k]/1e9
              ,molecular_cal_pulse[filter_1000k],'b'
              ,interf_freq[filter_1000k]/1e9
              ,combined_hi_cal_pulse[filter_1000k],'r')
            plt.legend(('mi2_1k','mi2_1000k','c_1000k'),'lower left')
        ax=plt.gca()
        ax.set_yscale('log')
        plt.grid(True)
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('Counts/energy')
        plt.title('Energy normalized scan-log '+start_time_str)
        plt.show(block=False)
    
        mol_1k = molecular_cal_pulse[filter_1k]
        mol_1000k = molecular_cal_pulse[filter_1000k]
        freq_1k = interf_freq[filter_1k]
        freq_1000k =interf_freq[filter_1000k]
        combined_hi_cal_pulse=combined_hi_cal_pulse[filter_1000k]
 
        #select only those points with less than 12 photons in cal pulse
        #from the narrow scan
        mask= np.array([ ((s < 12)) for s in mol_1k[:] ])

        
        freq_1k_i2 = freq_1k[mask]
        mol_1k_i2  = mol_1k[mask]
        #select all other frequencys from the wide scan
        mask = np.array([((f > np.nanmax(freq_1k_i2)              
                 or f <= np.nanmin(freq_1k_i2)))for f in freq_1000k[:]])
    
        mol_1000k_i2 = mol_1000k[mask]
        freq_1000k_i2 = freq_1000k[mask]

        #concatenate mol i2 filter signals--frequencys will be out of order 
        mol_i2 = np.concatenate((mol_1000k_i2,mol_1k_i2/1000.0),0)
        freq_i2 = np.concatenate((freq_1000k_i2,freq_1k_i2),0)

        #repeat for i2a signals if they are present
        if hasattr(rs,'molecular_i2a_cal_pulse'):
                mol_1k_i2a = molecular_i2a_cal_pulse[filter_1k]
                mol_1000k_i2a = molecular_i2a_cal_pulse[filter_1000k]
                mask= np.array([((s < 12)) for s in mol_1k_i2a[:] ])
                freq_1k_i2a = freq_1k[mask]
                mol_1k_i2a  = mol_1k_i2a[mask]
                mask = np.array([((f > np.nanmax(freq_1k_i2a)                 
                   or f <= np.nanmin(freq_1k_i2a)))for f in freq_1000k[:]])
                mol_1000k_i2a = mol_1000k_i2a[mask]
                freq_1000k_i2a = freq_1000k[mask]
        
                #concatenate mol i2a filter signals
                mol_i2a = np.concatenate((mol_1000k_i2a,mol_1k_i2a/1000.0),0)
                freq_i2a = np.concatenate((freq_1000k_i2a,freq_1k_i2a),0)
        
        #rearange mol in freq sorted order
        indices =np.argsort(freq_i2)
        freq_mol_i2 = freq_i2[indices]
        molecular_i2 = mol_i2[indices]
        if hasattr(rs,'molecular_i2a_cal_pulse'):
           indices_i2a = np.argsort(freq_i2a)
           molecular_i2a = mol_i2a[indices_i2a]  
           freq_mol_i2a = freq_i2a[indices_i2a]
        else:
            molecular_i2a = []
            freq_mol_i2a  = []
        
        #rearange combined hi in freq sorted order
        
        indices = np.argsort(freq_1000k)
        combined_hi = combined_hi_cal_pulse[indices].copy()
        combined_freq = freq_1000k[indices].copy()
        """
        combined_hi = combined_hi_cal_pulse.copy()
        combined_freq = freq_1000k.copy()
        """
        plt.figure(4002)
        if hasattr(rs,'molecular_i2a_cal_pulse'):
            plt.plot(freq_mol_i2/1e9
               ,molecular_i2,'b'
               ,combined_freq/1e9
               ,combined_hi,'r'
               ,freq_mol_i2a/1e9
               ,molecular_i2a,'k')

            plt.legend(('mol','comb_hi','mol_i2a'),'lower left')
        else:
            plt.plot(freq_mol_i2/1e9,molecular_i2,'b'
             ,combined_freq/1e9,combined_hi,'r')
            plt.legend(('mol','combined'),'lower left')

        ax=plt.gca()
        ax.set_yscale('log')
        plt.grid(True)
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('counts')
        plt.title('composite scan ' +start_time_str)

        plt.figure(4003)
        plt.plot(freq_mol_i2/1e9,molecular_i2,'b'
             ,combined_freq/1e9,combined_hi,'r')
        plt.legend(('mol','combined'),'lower left')
        ax=plt.gca()
        ax.set_yscale('log')
        plt.grid(True)
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('counts')
        plt.title('composite I2 scan' +start_time_str)


        if hasattr(rs,'molecular_i2a_cal_pulse'):
            plt.figure(4004)
            plt.plot(freq_mol_i2a/1e9,molecular_i2a,'b'
                ,combined_freq/1e9,combined_hi,'r')
            plt.legend(('mol_i2a','combined'),'lower left')
            ax=plt.gca()
            ax.set_yscale('log')
            plt.grid(True)
            plt.xlabel('Frequency (GHz)')
            plt.ylabel('counts')
            plt.title('composite I2A scan' +start_time_str)
            plt.show(block=False)


        return combined_freq,combined_hi,freq_mol_i2,molecular_i2,freq_mol_i2a,molecular_i2a
    

    def load_i2_spectrum():
        #read theoretical i2 transmission file
        try:
            [header,i2t]=ru.readascii(locate_file('dfb_i2_spectra.txt'))
            dfb_to_i2_scale = rs_constants['i2_absorption_scale_factor']
            i2t_trans=10**(dfb_to_i2_scale*i2t[:,1])
            i2t_freq = i2t[:,0]*1e9
            print
            print 'DFB I2 scan used for basis of calibration'
            print
        except IOError: # dfb not found, or error reading 
            [header,i2t]=ru.readascii(locate_file('I2cell_272_31_extended_freq.txt'))
            i2t_freq  = i2t[:,0]*1e9
            i2t_trans = i2t[:,2]
            print
            print 'Theoretical scan used as basis for calibration'
            print
         
        #mask = np.array([ ((f >= freq_1000k[0]) and (f < freq_1000k[-1])) for f in i2t_freq[:]])
        #i2t_trans = i2t_trans[mask]
        #i2t_freq =i2t_freq[mask]
        #normalize I2 spectrum to 1--ie eliminate continuum absorption
        #lo_freq_lmt = -3e9
        #hi_freq_lmt = -2e9
        #mask = np.array([ ((f >= 2e9) and (f < 3e9)) for f in i2t_freq])
        #mask = np.array([ ((f >= lo_freq_lmt) and (f < hi_freq_lmt)) for f in i2t_freq])
        
        
        #move set point frequency to another line?
        print 'Do you want to move to a different freq set point?'
        offset_setpoint = rs_constants['lock_point_freq_offset']
        #will present all results on the theoretical scan frequency axis
        freq = i2t_freq
        while 1:
            plt.figure(4005)
            plt.clf()
            plt.plot(freq/1e9, i2t_trans,'k')

            ax=plt.gca()
            ax.set_yscale('log')
            plt.grid(True)
            plt.xlabel('Frequency (GHz)')
            plt.ylabel('Counts')

            title_str='I2 Transmission,lockpoint offset= %4.3f GHz' %(offset_setpoint)
            plt.title(title_str)
            plt.show(block=False)

            try:
               string = raw_input('Lock point frequency offset (GHz), (CR = exit) ? ')
               if len(string) == 0:
                   break
               else:
                   #convert offset from GHz to Hz
                   offset_setpoint = np.float(string)*1e9
                   freq=freq+offset_setpoint
            except:
               print 'input error--try again'



        norm = np.max(i2t_trans)
        i2t_trans = i2t_trans/norm
        return i2t_trans,freq
  

    #main code starts here 
    #********************************************************************************

    times =rs.rs_raw.times
    #select time interval
    [date_str,start_time,start_time_str,end_time,end_time_str]=select_times(times)

    #plot raw counts
    plot_raw_counts(rs.rs_raw,start_time_str)
    print 
    manual_scan_flag = raw_input("enter 'm' if this a manual scan, else CR for normal cal scan ")
    manual_scan_flag = manual_scan_flag.find('m') >= 0
    if manual_scan_flag:
        [frequency_comb,combined_hi,frequency_i2,molecular_i2
          ,frequency_i2a,molecular_i2a] = make_manual_scan(rs.rs_raw
          ,start_time,end_time,start_time_str )  
    else:    
        #return composite scans in frequency sorted order
        [frequency_comb,combined_hi,frequency_i2,molecular_i2
          ,frequency_i2a,molecular_i2a] = make_composite_scans(rs.rs_raw
          ,start_time,end_time,start_time_str )

    #selected and rearanged buffers are now:
    #       combined_hi[freq_combined]
    #       molecular_i2a[frequency_i2a]
    #       molecular_i2[frequency_mol_i2]

    #load the model i2 transmission spectrum and it's freq scale
    [i2t_transmission,i2t_freq] = load_i2_spectrum()


  
    #compute and normalize i2 transmission
    #plot and allow adjustment of frequency offset between measurement and model
    print 'Do you want to offset the zero point of the measured frequencys?'
    offset_measured = 0.0
    
    #will present all results on the theoretical scan frequency axis
    while 1:
        freq = i2t_freq.copy()
        freq_comb = frequency_comb.copy()
        freq_i2 = frequency_i2.copy()
        if len(molecular_i2a):
            freq_i2a= frequency_i2a.copy()
           
        #interpolate to frequencys in theoretical i2 spectrum file
        comb_hi = np.interp(freq,freq_comb+offset_measured,combined_hi)
        mol_i2 = np.interp(freq,freq_i2+offset_measured,molecular_i2)
        freq = freq + offset_measured
        if len(molecular_i2a):
            mol_i2a = np.interp(freq,freq_i2a+offset_measured,molecular_i2a) 

        #buffers interpolated to freq scale provided with spectrum are now:
        #     comb_hi[freq]
        #     mol_i2[freq]
        #     mol_i2a[freq]


        f=plt.figure(4007)
        f.clear()
        if len(molecular_i2a):
            plt.plot(freq,comb_hi,'r',freq,mol_i2,'b',freq,mol_i2a,'g')
            plt.legend(('comb_hi','mol_i2','mol_i2a'))
        else:
            plt.plot(freq,comb_hi,'r',freq,mol_i2,'b')
            plt.legend(('comb_hi','mol_i2'))
        plt.xlabel('frequency (GHz)')
        plt.grid(True)
        plt.title('Interpolated spectra')


        #trim to region around absorption line
        f_lmt = 4e9
        comb_hi = comb_hi[np.abs(freq)<f_lmt]
        mol_i2  = mol_i2[np.abs(freq)<f_lmt]
        if len(molecular_i2a):
            mol_i2a = mol_i2a[np.abs(freq)<f_lmt]
        i2t_trans = i2t_transmission[np.abs(freq)<f_lmt]
        freq = freq[np.abs(freq)<f_lmt]        
        
        #generate polynomial fit to combined channel transmission
        #this gives the transmission function of the etalon
        p=np.polyfit(freq,comb_hi,3)
        fit_comb_hi=np.polyval(p,freq)

        
        zero_index = np.argmin(np.abs(freq))
      
        norm_factor = fit_comb_hi[zero_index]
       
        #normalize to lock_pt value of combined channel
        print 'norm factor =', norm_factor
        fit_comb_hi=fit_comb_hi/norm_factor

        mol_i2= mol_i2/norm_factor
        comb_hi = comb_hi/norm_factor
        if  len(molecular_i2a):
           mol_i2a=mol_i2a/norm_factor
           i2a_trans = mol_i2a/comb_hi
       
        #multiply fitted comb_hi by i2 transmission to get fitted molecular 
        fit_mol = fit_comb_hi * i2t_trans
        
        i2_gain_ratio = np.mean(mol_i2[fit_mol>0.1] / fit_mol[fit_mol>0.1])
        fit_mol = fit_mol * i2_gain_ratio

        #this is the measured i2 transmission before any fitting
        i2_trans = mol_i2 /comb_hi
        i2_trans = i2_trans / i2_gain_ratio
   
        #find min of transmission curve for molecular channel
        #mask = np.array([ ((f >= -0.5e9) and (f < 0.5e9)) for f in freq])
        #Cam = np.nanmin(mol_i2[mask])
        Cam = fit_mol[zero_index] / fit_comb_hi[zero_index] 


        f=plt.figure(4008)
        f.clear()
        plt.plot(freq/1e9,comb_hi,'r',freq/1e9,fit_comb_hi,'k'
                 ,freq/1e9,mol_i2,'b',freq/1e9,fit_mol,'k')
        plt.grid(True)
        ax=plt.gca()
        ax.set_xlim((-6,6))
        ax.set_yscale('log')
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('Normalized counts/J')
        plt.title('fit to combined_hi peak')
        plt.legend(('comb_hi','fit_comb_hi','mol_i2','fit_mol_i2'))
        plt.show(block=False)
       
        try:
           string = raw_input('offset measured frequency in GHz,(CR = done) ? ')
           if len(string) == 0:
               break  
           else:
               #convert offset from GHz to Hz
               offset_measured=np.float(string)*1e9
        except:
           print 'input error--try again'
       




  
   
    print 'gain_ratio= ', i2_gain_ratio
    print 'Cam = %8.2e ' %(Cam)
    print 'min_i2_transmission= %8.2e' %(Cam/i2_gain_ratio)
    print '1/i2_trans = %i' %(int(1.0/(Cam/i2_gain_ratio)))

  
    if len(molecular_i2a):        
        #find i2a transmission at lock_pt
        min_i2a_trans = i2a_trans[zero_index]
        inx = np.array([ ((f >= 2.5e9) and (f < 3.5e9)) for f in freq])
        i2a_gain_ratio = np.mean(mol_i2a[inx] / fit_comb_hi[inx] )
        print 'i2a_gain_ratio= ', i2a_gain_ratio
        print 'Cam_i2a = %8.2e ' %(min_i2a_trans/i2a_gain_ratio)
        print 'min_i2a_transmission= %8.2e' %(min_i2a_trans)
        print '1/i2a_trans = %i' %(int(1.0/min_i2a_trans))
    
        plt.figure(4009)
        plt.plot(freq/1e9,fit_comb_hi*i2_gain_ratio,'c',freq/1e9,fit_comb_hi,'r'
                 ,freq/1e9,fit_mol,'b',freq/1e9,mol_i2a,'g')
        plt.grid(True)
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('Normalized counts/J')
        plt.title('I2 scan')
        plt.show(block=False)
        plt.legend(('fit_comb*gain','fit_comb','fit_mol','mol_i2a'))
    else:    
        plt.figure(4009)
        plt.plot(freq/1e9,fit_comb_hi*i2_gain_ratio,'c',freq/1e9,fit_comb_hi,'r',freq/1e9,fit_mol,'b')
        plt.grid(True)
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('Normalized counts/J')
        plt.title('I2 scan')
        plt.show(block=False)
        plt.legend(('fit_comb*gain','fit_comb','fit_mol'))
   

    
           
    [cal_start,pulse,cal_end] = rs_constants['apd_pulse_timing'] 
    
    #write i2 scan file
    #get path to directory
    dir_path=calibration_path_for(instrument,rs.rs_raw.times[0],process_defaults)
    #define start time for use in filename
    str = start_time.strftime("%Y%m%dT%H%M")
    filename=os.path.join(dir_path,'i2-scan-' +str+'.cal')
    fileid=open(filename,'w')
    os.chmod(filename,0664)
    print >>fileid, '# Calibration scan as function of frequency offset from I2 line'
    print >>fileid, '# calibration scan data aquired on ' + start_time.strftime("%d-%b-%y")\
                 +' at '+ start_time.strftime("%H:%M") +' UT'
    print >>fileid, '# file created on '+ datetime.now().strftime("%d-%b-%y") +' at ' \
                 +datetime.now().strftime("%H:%M") + ' UT'
    print >>fileid, '# t_begin_cal_pulse= %4.2e ;  start time of cal pulse (sec).' %(cal_start)
    print >>fileid, '# t_end_cal_pulse=   %4.2e ;  end time of cal pulse (sec).' %(cal_end)
    print >>fileid, '# pulse_durration=   %4.2e ;  laser pulse durration (sec).' \
          %(rs_constants['binwidth'])
    print >>fileid, '# ratio of mol to combined channel gain = %5.1f ' %(i2_gain_ratio)
    print >>fileid, '# Cam = %8.3e ' %(Cam) 
    print >>fileid, '# Min iodine transmission =   %8.1e,  1/(min_trans) = %8.1e'\
                        %(Cam*i2_gain_ratio, 1/(Cam*i2_gain_ratio))
    if hasattr(rs.rs_raw,'molecular_i2a_counts'):
        print >>fileid, '# ratio of mol_i2a to combined channel gain = %5.1f ' %(i2a_gain_ratio)
        print >>fileid, '# Cam_i2a = %8.3e ' %(min_i2a_trans/i2a_gain_ratio) 
        print >>fileid, '# Min iodine_argon transmission =   %8.3e,  1/(min_trans) = %8.3e'\
                        %(min_i2a_trans, 1/min_i2a_trans) 
    print >>fileid, '# '
    if hasattr(rs.rs_raw,'molecular_i2a_counts'):
        print >>fileid,\
           '#freq(GHz)    combined    molecular   i2_theory   i2_trans    i2a_trans    molecular_i2a '
        for i in range(np.size(freq)):
            print >>fileid, '%8.5f    %9.5f   %9.5f   %9.5f   %9.5f   %9.5f   %9.5f'\
	    %(freq[i]/1e9 ,fit_comb_hi[i], fit_mol[i]
                 ,i2t_trans[i],i2_trans[i],i2a_trans[i], mol_i2a[i])
    else:    
        print >>fileid, '#freq(GHz)  combined  molecular i2_measured i2_theory'
        for i in range(np.size(freq)):
            print >>fileid, '%8.5f   %9.5f  %9.5f  %9.5f  %9.5f'\
	      %(freq[i]/1e9 ,fit_comb_hi[i], fit_mol[i]
                 ,i2t_trans[i],i2_trans[i])
    fileid.close()
    
    print '\nnew i2 scan file=\n '+filename
    
def make_i2_scan_file_manual(instrument,rs,rs_constants,process_defaults=None):
    """make_i2_scan_file_manual(instrument,rs,rs_constants)
       make an i2_scan file from a manual calibration scan and store it in the
       calibration directory--instrument temperatures must be stable
       during measurement to minimize interferometer drift"""
   
    interf_freq = rs.rs_raw.interf_freq.copy()
    times = rs.rs_raw.times.copy()
    molecular_cal_pulse = rs.rs_raw.molecular_cal_pulse.copy()
    combined_hi_cal_pulse = rs.rs_raw.combined_hi_cal_pulse.copy()
    #assume system is locked at start of record--set this frequency to zero
    interf_freq=interf_freq-np.mean(interf_freq[:5])
    
    if hasattr(rs.rs_raw, 'molecular_i2a_counts'):
        molecular_i2a_cal_pulse = rs.rs_raw.molecular_i2a_cal_pulse.copy()

    #select time interval
    date_str = datetime.strftime(times[0], '%d-%b-%y ')
    start_time_str = raw_input("Start time=? H:M:S or H:M ")
    start_time_str = date_str + start_time_str
    if start_time_str.count(':') == 2:
      start_time = datetime.strptime(start_time_str, '%d-%b-%y %H:%M:%S')
    else:
      start_time = datetime.strptime(start_time_str, '%d-%b-%y %H:%M')
    end_time_str = raw_input("End time H:M:S or H:M ")
    end_time_str = date_str + end_time_str
    if end_time_str.count(':') == 2:
       end_time = datetime.strptime(end_time_str, '%d-%b-%y %H:%M:%S')
    else:
       end_time = datetime.strptime(end_time_str, '%d-%b-%y %H:%M')

    
    time_mask = np.array([ ((x >= start_time) and (x<end_time)) 
                             for x in times[:] ])
    
    interf_f = interf_freq[time_mask].copy()

    molecular_cal_pulse = molecular_cal_pulse[time_mask].copy()
    combined_hi_cal_pulse = combined_hi_cal_pulse[time_mask].copy()
    if hasattr(rs.rs_raw,'molecular_i2a_counts'):
        molecular_i2a_cal_pulse = molecular_i2a_cal_pulse[time_mask].copy()

    np.set_printoptions(threshold=np.NaN)
   

    plt.figure(4000)
    if hasattr(rs.rs_raw,'molecular_i2a_cal_pulse'):
        plt.plot(interf_f/1e9
             ,molecular_cal_pulse,'b'   
             ,interf_f/1e9
             ,combined_hi_cal_pulse,'r'
             ,interf_f/1e9
             ,molecular_i2a_cal_pulse,'k')
        plt.legend(['mol','comb','mol_i2a'],'lower left')
    else:
        plt.plot(interf_f/1e9
             ,molecular_cal_pulse,'b'
             ,interf_f/1e9
             ,combined_hi_cal_pulse,'r')
        plt.legend(('mol','comb'),'lower left')
    plt.grid(True)
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Counts')
    plt.title('Raw scan--linear '+start_time_str)

    plt.show(block=False)
    
    #logrithmic plot
    plt.figure(4001)
    if hasattr(rs.rs_raw,'molecular_i2a_cal_pulse'):
        plt.plot(interf_f/1e9
             ,molecular_cal_pulse,'b'   
             ,interf_f/1e9
             ,combined_hi_cal_pulse,'r'
             ,interf_f/1e9
             ,molecular_i2a_cal_pulse,'k')
        plt.legend(('mol','comb','mol_i2a'),'lower left')
    else:
        plt.plot(interf_f/1e9
             ,molecular_cal_pulse,'b'
             ,interf_f/1e9
             ,combined_hi_cal_pulse,'r')
        plt.legend(('mol','comb'),'lower left')
    ax=plt.gca()
    ax.set_yscale('log')
    plt.grid(True)
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Counts')
    plt.title('Raw scan-log '+start_time_str)
    plt.show(block=False)
   


    lo_mol_count = raw_input("Lowest valid mol count = ?  ")
    lo_mol_count = np.float(lo_mol_count)
    molecular_cal_pulse[molecular_cal_pulse <= lo_mol_count] = np.max(molecular_cal_pulse)/1e4
    if hasattr(rs.rs_raw,'molecular_i2a_cal_pulse'):
        molecular_i2a_cal_pulse[molecular_i2a_cal_pulse <= lo_mol_count] = np.max(molecular_i2a_cal_pulse)/1e4

    #logrithmic plot with bottom of line replaced with max/1e4
    plt.figure(4002)
    if hasattr(rs.rs_raw,'molecular_i2a_cal_pulse'):
        plt.plot(interf_f/1e9
             ,molecular_cal_pulse,'b'   
             ,interf_f/1e9
             ,combined_hi_cal_pulse,'r'
             ,interf_f/1e9
             ,molecular_i2a_cal_pulse,'k')
        plt.legend(('mol','comb','mol_i2a'),'lower left')
    else:
        plt.plot(interf_f/1e9
             ,molecular_cal_pulse,'b'
             ,interf_f/1e9
             ,combined_hi_cal_pulse,'r')
        plt.legend(('mol','comb'),'lower left')
    ax=plt.gca()
    ax.set_yscale('log')
    plt.grid(True)
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Counts')
    plt.title('scan with low limit applied-log '+start_time_str)
    plt.show(block=False)


    #rearange mol in freq sorted order
    indices =np.argsort(interf_f)
    freq_mol_i2 = interf_f[indices].copy()
    molecular_i2 = molecular_cal_pulse[indices].copy()
    combined_hi_cal_pulse = combined_hi_cal_pulse[indices].copy()
    if hasattr(rs.rs_raw,'molecular_i2a_cal_pulse'):
       molecular_i2a = molecular_i2a_cal_pulse[indices].copy()  

    #read theoretical i2 transmission file
    try:
        [header,i2t]=ru.readascii(locate_file('dfb_i2_spectra.txt'))
        dfb_to_i2_scale = rs_constants['i2_absorption_scale_factor']
        i2t_trans=10**(dfb_to_i2_scale*i2t[:,1])
        i2t_freq = i2t[:,0]*1e9
    except IOError:    
        [header,i2t]=ru.readascii(locate_file('I2cell_272_31_extended_freq.txt'))
        i2t_freq  = i2t[:,0]*1e9
        i2t_trans = i2t[:,2]

    #mask = np.array([ ((f >= freq_mol_i2[0]) and (f < freq_mol_i2[-1])) for f in i2t_freq[:]])
    #i2t_trans = i2t_trans[mask]
    #i2t_freq =i2t_freq[mask]
    #normalize theoretical spectrum to 1--ie eliminate continuum absorption
    mask = np.array([ ((f >= 2e9) and (f < 3e9)) for f in i2t_freq])
    norm = nanmean(i2t_trans[mask])
    i2t_trans = i2t_trans/norm




    #move set point frequency to another line?
    print 'Do you want to move to a different freq set point in theoretical spectrum?'
    offset_setpoint = rs_constants['lock_point_freq_offset']
    #will present all results on the theoretical scan frequency axis
    freq = i2t_freq
    while 1: 
        plt.figure(4005)
        plt.clf()
        plt.plot(freq/1e9, i2t_trans,'k')  

        ax=plt.gca()
        ax.set_yscale('log')
        plt.grid(True)
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('Counts')
        
        title_str='I2 Transmission,lockpoint offset= %4.3f GHz' %(offset_setpoint/1e9)
        plt.title(title_str)
        plt.show(block=False)
        
        try:
           string = raw_input('Lock point frequency offset (GHz), (CR = exit) ? ')
           if len(string) == 0:
               break
           else:
               #convert offset from GHz to Hz
               offset_setpoint = np.float(string)*1e9
               freq=freq+offset_setpoint
        except:
               print 'input error--try again'
   
    
    #compute and normalize i2 transmission
    #plot and allow adjustment of frequency offset between measurement and model
    print 'Do you want to offset the zero point of the measured frequencys?'
    offset = 0.0
    #will present all results on the theoretical scan frequency axis
    while 1:
    
        #interpolate to frequencys in theoretical i2 spectrum file
        combined_hi = np.interp(freq,freq_mol_i2+offset,combined_hi_cal_pulse)
        mol_i2 = np.interp(freq,freq_mol_i2+offset,molecular_i2)
        if hasattr(rs.rs_raw,'molecular_i2a_cal_pulse'):
            mol_i2a = np.interp(freq,freq_mol_i2+offset,molecular_i2a) 


        i2_trans = mol_i2/combined_hi
        #normalize to region between +2 and +3 GHz
        mask = np.array([ ((f >= 2e9) and (f < 3e9)) for f in freq])
        i2_gain_ratio = nanmean(i2_trans[mask])
        i2_trans=i2_trans/i2_gain_ratio
        if hasattr(rs.rs_raw,'molecular_i2a_counts'):
            #compute and normalize i2a transmission 
            i2a_trans = mol_i2a/combined_hi
            i2a_gain_ratio = nanmean(i2a_trans[mask])
            i2a_trans=i2a_trans/i2a_gain_ratio
       
        
        plt.figure(4006)
        plt.clf()
        freq_offset_theory = rs_constants['lock_point_freq_offset']
        #plot i2 transmission and theoretical values
        if hasattr(rs.rs_raw,'molecular_i2a_counts'):
            plt.plot(freq/1e9-freq_offset_theory, i2_trans,'r'
                ,freq/1e9,i2t_trans,'g'
                ,freq/1e9-freq_offset_theory,i2a_trans,'k')
        else:
            plt.plot(freq/1e9-freq_offset_theory, i2_trans,'r'
                 ,freq/1e9, i2t_trans,'k')  

        ax=plt.gca()
        ax.set_yscale('log')
        plt.grid(True)
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('Counts')
        if hasattr(rs.rs_raw,'molecular_i2a_counts'):
            plt.legend(('i2','i2_th','i2a'),'lower right')
        else:
            plt.legend(('i2','i2_filt','i2_th'),'lower right')
        title_str='I2 Transmission, gain_ratio= %4.3f' %(i2_gain_ratio)
        plt.title(title_str)
        plt.show(block=False)
        
        try:
           string = raw_input('offset frequency in GHz,(CR = exit) ? ')
           if len(string) == 0:
               break
           else:
               #convert offset from GHz to Hz
               offset=np.float(string)*1e9
        except:
           print 'input error--try again'

    mask = np.array([ ((f >= -3e9) and (f < 3e9)) for f in freq])
    freq = freq[mask].copy()
    i2_trans = i2_trans[mask].copy()
    i2t_trans = i2t_trans[mask].copy()
    combined_hi = combined_hi[mask].copy()
    mol_i2 = mol_i2[mask].copy()
    if hasattr(rs.rs_raw,'molecular_i2a_counts'):
        i2a_trans = i2a_trans[mask].copy()
    
    plt.figure(4007)
    plt.clf()
    #plot i2 transmission and theoretical values
    if hasattr(rs.rs_raw,'molecular_i2a_counts'):
        plt.plot(freq/1e9, i2_trans,'r'
                ,freq/1e9,i2t_trans,'g'
                ,freq/1e9,i2a_trans,'k')
    else:
        plt.plot(freq/1e9, i2_trans,'r',freq/1e9, i2t_trans,'k')  

    ax=plt.gca()
    ax.set_yscale('log')
    plt.grid(True)
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Counts')
    if hasattr(rs.rs_raw,'molecular_i2a_counts'):
            plt.legend(('i2','i2_th','i2a'),'lower right')
    else:
            plt.legend(('i2','i2_filt','i2_th'),'lower right')
    title_str='I2 Transmission, gain_ratio= %4.3f' %(i2_gain_ratio)
    plt.title(title_str)
    plt.show(block=False)
        
    #find min of transmission curve for molecular channel
    mask = np.array([ ((f >= -0.5e9) and (f < 0.5e9)) for f in freq])
    min_i2_trans = np.nanmin(i2_trans[mask])
   
   
    print 'gain_ratio= ', i2_gain_ratio
    print 'Cam = %8.2e ' %(min_i2_trans/i2_gain_ratio)
    print 'min_i2_transmission= %8.2e' %(min_i2_trans)
    print '1/i2_trans = %i' %(int(1.0/min_i2_trans))

    if hasattr(rs.rs_raw,'molecular_i2a_counts'):        
            #find min of i2a transmission curve
            min_i2a_trans = np.nanmin(i2a_trans[mask])
            print 'i2a_gain_ratio= ', i2a_gain_ratio
            print 'Cam_i2a = %8.2e ' %(min_i2a_trans/i2a_gain_ratio)
            print 'min_i2a_transmission= %8.2e' %(min_i2a_trans)
            print '1/i2a_trans = %i' %(int(1.0/min_i2a_trans))


    #normalize to max of combined channel
    
    
    p=np.polyfit(freq[np.abs(freq)<1.5e9],combined_hi[np.abs(freq)<1.5e9],2)
    center_comb_hi=np.polyval(p,freq[np.abs(freq)<1.5e9])
    norm_factor = np.max(center_comb_hi)
    print 'norm factor(max) = ',norm_factor
    norm_factor =np.polyval(p,0)
    print 'norm factor(0) =', norm_factor
    plt.figure(4008)
    plt.plot(freq/1e9,combined_hi,'r',freq[np.abs(freq)<1.5e9]/1e9,center_comb_hi,'k')
    plt.grid(True)
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Counts')
    plt.title('fit to combined_hi peak')
    ax=plt.gca()
    plt.show(block=False)

    mol_i2= mol_i2/norm_factor
    combined_hi = combined_hi/norm_factor
    if hasattr(rs.rs_raw,'molecular_i2a_counts'):
        mol_i2a=mol_i2a/norm_factor

    plt.figure(4008)
    plt.clf()
    #plot i2 transmission and theoretical values
    if hasattr(rs.rs_raw,'molecular_i2a_counts'):
        plt.plot(freq/1e9,i2t_trans,'c'
                ,freq/1e9, combined_hi,'r'
                ,freq/1e9,mol_i2,'b'
                ,freq/1e9,mol_i2a,'k')
                 
    else:
        plt.plot(freq/1e9, i2t_trans,'c'
                ,freq/1e9, combined_hi,'r'
                ,freq/1e9, mol_i2,'b')

    ax=plt.gca()
    ax.set_yscale('log')
    plt.grid(True)
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Normalized counts')
    if hasattr(rs.rs_raw,'molecular_i2a_counts'):
            plt.legend(('i2_th','comb_hi','mol_i2','mol_i2a'),'lower right')
    else:
            plt.legend(('i2_th','comb_hi','mol)i2'),'lower right')
    title_str='Final I2 scan, gain_ratio= %4.3f' %(i2_gain_ratio)
    plt.title(title_str)
    plt.show(block=False)

           
    [cal_start,pulse,cal_end] = rs_constants['apd_pulse_timing'] 
    
    #write i2 scan file
    #get path to directory
    dir_path=calibration_path_for(instrument,rs.rs_raw.times[0],process_defaults)
    #define start time for use in filename
    str = start_time.strftime("%Y%m%dT%H%M")
    filename=os.path.join(dir_path,'i2-scan-' +str+'.cal')
    fileid=open(filename,'w')
    os.chmod(filename,0664)
    print >>fileid, '# Calibration scan as function of frequency offset from I2 line'
    print >>fileid, '# calibration scan data aquired on ' + start_time.strftime("%d-%b-%y")\
                 +' at '+ start_time.strftime("%H:%M") +' UT'
    print >>fileid, '# file created on '+ datetime.now().strftime("%d-%b-%y") +' at ' \
                 +datetime.now().strftime("%H:%M") + ' UT'
    print >>fileid, '# t_begin_cal_pulse= %4.2e ;  start time of cal pulse (sec).' %(cal_start)
    print >>fileid, '# t_end_cal_pulse=   %4.2e ;  end time of cal pulse (sec).' %(cal_end)
    print >>fileid, '# pulse_durration=   %4.2e ;  laser pulse durration (sec).' \
          %(rs_constants['binwidth'])
    print >>fileid, '# ratio of mol to combined channel gain = %5.1f ' %(i2_gain_ratio)
    print >>fileid, '# Cam = %8.3e ' %(min_i2_trans/i2_gain_ratio) 
    print >>fileid, '# Min iodine transmission =   %8.1e,  1/(min_trans) = %8.1e'\
                        %(min_i2_trans, 1/min_i2_trans)
    if hasattr(rs.rs_raw,'molecular_i2a_counts'):
        print >>fileid, '# ratio of mol_i2a to combined channel gain = %5.1f ' %(i2a_gain_ratio)
        print >>fileid, '# Cam_i2a = %8.3e ' %(min_i2a_trans/i2a_gain_ratio) 
        print >>fileid, '# Min iodine_argon transmission =   %8.3e,  1/(min_trans) = %8.3e'\
                        %(min_i2a_trans, 1/min_i2a_trans) 
    print >>fileid, '# '
    if hasattr(rs.rs_raw,'molecular_i2a_counts'):
        print >>fileid,\
           '#freq(GHz)    combined    molecular   i2_theory   i2_trans    i2a_trans    molecular_i2a '
        for i in range(np.size(freq)):
            print >>fileid, '%8.5f    %9.5f   %9.5f   %9.5f   %9.5f   %9.5f   %9.5f'\
	    %(freq[i]/1e9 ,combined_hi[i], mol_i2[i]
                 ,i2t_trans[i],i2_trans[i],i2a_trans[i], mol_i2a[i])
    else:    
        print >>fileid, '#freq(GHz)  combined  molecular i2_measured i2_theory'
        for i in range(np.size(freq)):
            print >>fileid, '%8.5f   %9.5f  %9.5f  %9.5f  %9.5f'\
	      %(freq[i]/1e9 ,combined_hi[i], mol_i2[i]
                 ,i2t_trans[i],i2_trans[i])
    fileid.close()

    print '\nnew i2 scan file=\n '+filename


def compute_combined_1064_to_combined_hi_gain(profiles,processing_defaults):
    """compute_combined_1064_to_combined_hi_gain(profiles,processing_defaults)
       Use a layer with large particles where the backscatter cross section at
       1064 is approximately equal to the backscatter cross section at 532 to
       compute the gain ratio of the combined_1064 channel to the combined_532
       channel, corrected for the portion of the layer backscatter expected from
       molecules and the reduced attenuation at 1064 in the layers below the
       calibration layer, for details see:
                  'hsrl/derivations/calibration_1064_channel.pdf')"""
    alts = profiles.msl_altitudes
    count_ratio = profiles.combined_1064_counts[0,:] / profiles.combined_counts[0,:] 
    inv_scattering_ratio = profiles.inv.beta_r_backscat / profiles.inv.beta_a_backscat[0,:]
    mol_backscat_corr = (1.0 + inv_scattering_ratio / 16.0)/(1.0+inv_scattering_ratio)
    mol_od_corr = (15.0/16.0) * (profiles.inv.optical_depth[0,:]\
                                 - profiles.inv.optical_depth_aerosol[0,:])
    mol_ext_corr = np.exp(2.0*mol_od_corr)
    angstrom_coef = processing_defaults.get_value('color_ratio','angstrom_coef')
    aerosol_od_corr = (1-0.5**angstrom_coef)*profiles.inv.optical_depth_aerosol[0,:]
    aerosol_ext_corr = np.exp(2.0*aerosol_od_corr)
    gain_ratio = count_ratio /( mol_backscat_corr * mol_ext_corr * aerosol_ext_corr)
    
    plt.figure(1000)
    lines=plt.plot(count_ratio,alts,'k',gain_ratio,alts,'r',mol_backscat_corr,alts,'b'\
             ,mol_ext_corr,alts,'g',aerosol_ext_corr,alts,'c')
    plt.setp(lines[:],linewidth=2.0)
    ax=plt.gca()
    ax.set_xscale('log')
    plt.legend(('c_ratio','g_ratio','m_bs','m_ext','a_ext'),'upper left')
    plt.ylabel('Altitude(km)')
    plt.xlabel('Gain ratio')
    plt.grid('on')
    ax.set_xlim((1e-2,2.0))
    plt.title('gain ratio, angstrom_coef = ' +str(angstrom_coef))
    print 'altitude  gain_ratio  count_ratio  mol_back_corr  mol_ext_corr  aerosol_ext_corr'
    for i in range(len(alts)):
        print '%6.3f %13.2e %13.2e %13.2e %13.2e %13.2e'\
            %(alts[i]/1000.0,gain_ratio[i],count_ratio[i]\
            ,mol_backscat_corr[i],mol_ext_corr[i],aerosol_ext_corr[i])
    
    plt.show(block=False)





def make_diff_cross_pol(instrument,profiles,rs_cal,rs_constants,process_defaults=None,corr_adjusts=None):
    """
    make_diff_cross_pol(instrument,profiles,rs_cal,rs_constants)
    This uses raw_profile stream as input
    Must be very clean conditions with esentially no cross polarization present
    """
   
    start_time_str=profiles.start_time.strftime("%d-%b-%y %H:%M")
    end_time_str=profiles.end_time.strftime("%H:%M")
   
    #note raw_profiles have pileup correction applied--no other corrections 
    #signals
    comb_hi = np.zeros(len(profiles.sum_molecular_counts[0,:]))
    c_pol =              np.zeros_like(comb_hi)

    #differential geometry correction for cross polarization channel
    diff_geo_cpol =      np.ones_like(comb_hi)
    temp_diff_geo_cpol = np.ones_like(diff_geo_cpol)
    ref_bin = 1000  #normalize at bin 1000 as first guess
    dark_scale = 1.0
    make_constant_bin = len(comb_hi)-10
    norm_index = 500
    import matplotlib.pylab as plt
    plt.ion()
    bin_vec=np.arange(len(comb_hi))

    first_pass = True
    while 1:

        #baseline correction
        comb_hi[:] = profiles.sum_combined_hi_counts[0,:]- rs_cal.baseline.data[:,1]\
                *profiles.sum_transmitted_energy
        if hasattr(profiles,'sum_cross_pol_counts'):
            c_pol[:] = profiles.sum_cross_pol_counts[0,:]- rs_cal.baseline.data[:,4]\
                 *profiles.sum_transmitted_energy
        else:
            print "ERROR----no cross pol in data file----can not make 'cp_d_geo' file."
            return
        
        #dark count correction--no signal in dark correction applied
        c_pol = c_pol - profiles.sum_c_pol_dark_counts[0,0] * dark_scale
        comb_hi = comb_hi - profiles.sum_c_hi_dark_counts[0,0] * dark_scale

        raw_diff_geo_cpol = c_pol/comb_hi
      
            

        start_filter_index = 50
        temp_diff_geo_cpol[start_filter_index:] \
                  =(sg.savitzky_golay(raw_diff_geo_cpol[start_filter_index:]
                  ,51,3, deriv = 0)).copy()
        start_filter_index = 1000
        temp_diff_geo_cpol[start_filter_index:] \
                  =sg.savitzky_golay(temp_diff_geo_cpol[start_filter_index:]
                  ,201,3, deriv = 0)
                                     
        plt.figure(191)
        plt.clf()
        plt.plot(raw_diff_geo_cpol/temp_diff_geo_cpol[norm_index],bin_vec,'c')
        plt.grid(True)
        ax=plt.gca()
        #ax.set_xscale('log')
        #plt.xlabel('mol')
        #ax.set_xlim((1e-7, 1))
        ax.set_xlim((.5,1.5))
        plt.ylabel('Lidar bin number')
        plt.title('diff cpol geo')

      
        
        diff_geo_cpol = (temp_diff_geo_cpol/temp_diff_geo_cpol[norm_index]).copy()
        diff_geo_cpol[make_constant_bin:] = diff_geo_cpol[make_constant_bin]
        
        temp_diff_geo_cpol = temp_diff_geo_cpol/temp_diff_geo_cpol[norm_index]
 
       
        if not first_pass:    
            plt.figure(194)
            plt.clf()
            plt.plot(temp_diff_geo_cpol,bin_vec,'c',diff_geo_cpol,bin_vec,'r')
            plt.grid(True)
            ax=plt.gca()
            ax.set_xlim((.5,1.5))
            plt.ylabel('Lidar bin number')
            plt.title('diff geo comb cpol')

        first_pass = False
        print
        print 
        print 'normalize differential geometry at specified bin (or CR to continue)'
        ref_bin=raw_input('Set diff_geo==1 at bin#= ?,   ')
        if len(ref_bin)==0:
            break
        else:
            #make_constant_bin = int(ref_bin)
            dark_scale = raw_input('select scaling for combined dark count (1.0=no change) =? ')
            if len(dark_scale)==0:
                dark_scale = 1.0
            dark_scale = float(dark_scale)
            #make_constant_bin = int(make_constant_bin)
            #norm_index = make_constant_bin

    while 1:        
       
                

        input_str = raw_input('Make constant beyound bin #? (or CR to exit)   ')
        if len(input_str)==0:
            break

        #make constant after requested bin #
        make_constant_bin = int(input_str)
        temp_diff_geo_cpol = diff_geo_cpol.copy()
      
        temp_diff_geo_cpol[make_constant_bin:] = temp_diff_geo_cpol[make_constant_bin]
    
        plt.figure(197)
        plt.plot(diff_geo_cpol,bin_vec,'c',temp_diff_geo_cpol,bin_vec,'r')
        plt.grid(True)
        ax=plt.gca()
        ax.set_xlim((.5,1.5))
        plt.ylabel('Lidar bin number')
        plt.title('diff geo comb cpol')
    
      

    diff_geo_cpol = temp_diff_geo_cpol

    #don't allow  values less that 1.0
    #diff_geo_cpol[diff_geo_cpol<1.0] = 1.0
    
    #write diff_geofile
    #get path to directory
    dir_path=calibration_path_for(instrument,profiles.times[0],process_defaults)
    #define start time for use in filename
    str=profiles.start_time.strftime("%Y%m%dT%H%M")          
    filename=os.path.join(dir_path,'cpol_diff_geo_'+str+'.geo')
    fileid=open(filename,'w')
    #os.chmod(fileid,0664)
    print >>fileid, '#profile data from %s -->%s UTC' %(start_time_str,end_time_str)
    print >>fileid, '# Gains normalized range bin = %i ' %(make_constant_bin)
    print >>fileid, '#bin #  cpol/combined_hi'
    for i in range(len(comb_hi)):
        print >>fileid, '%5i   %10.4e'  \
	%(i ,diff_geo_cpol[i])
    fileid.close()
    
    print '\nnew diff_cpol_geofile=\n '+filename
    plt.show()
    return
