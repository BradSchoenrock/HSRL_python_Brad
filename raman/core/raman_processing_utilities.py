import copy
import numpy as np
import lg_base.core.array_utils as hau
import lidar.lidar_utilities as lu
import lidar.sg_extinction as lsge
from lidar.filtered_extinction import filtered_extinction
from lg_base.core.locate_file import locate_file
from datetime import datetime,timedelta
#import lg_base.core.read_utilities as ru
from lg_base.formats.vector_table import calibration_vector_table
# Try to use the much faster nanmean from bottleneck, otherwise fall back
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

import lg_base.core.array_utils as hau
from datetime import datetime

callist=dict(
    geo=dict(prefix='geofile_',suffix='.geo'),
    )

def calDataInfo(name):
    if name not in callist:
        r='Unknown cal file type '+name+'. Add it to callist in TableLibrarian.py'
        print r
        raise RuntimeError(r)
    return callist[name].copy()

import lg_dpl_toolbox.formats.VectorTableLibrarian

class TableLibrarian(lg_dpl_toolbox.formats.VectorTableLibrarian.VectorTableLibrarian):
    """ Librarian Initialzation for Calibration Tables
  
            :param siteid: source site id for the data source. typically an hsrl instrument or base directory
            :param datatype: file type to list
    """
    def __init__(self, basedir, instrument,datatype,**kwargs):
        super(TableLibrarian,self).__init__(instrumentname=instrument,datatype=calDataInfo(datatype),basedir=basedir,subdir='raman_cal',**kwargs)

def findFile(basedir,instrument,caltype,moment,*args,**kwargs):
    lib=TableLibrarian(basedir,instrument,caltype,*args,**kwargs)#,completeList=True)
    nar=lib(moment,moment+timedelta(days=31))#TableNarrator(getBasedir(datatype),datatype,caltype,caltype,st,et,completeList=True)
    for f in nar:
        return f
    return None

def findCalFile(*args,**kwargs):
    v=findFile(*args,**kwargs)
    if v is None:
        altp=kwargs.pop('alternativePath',None)
        r="Couldn't find cal file "+args[2]+' for '+args[1]+' time '+args[3].strftime('%Y.%m.%dT%H:%M:%S'+('' if altp is None else (' using alt path '+altp)))
        return (None, datetime(2100,1,1,0,0,0))
        #raise RuntimeError(r)
    return (v['path'],v['start']+v['width'])

def find_cal_file(basedir, isntname, file_type, start_time,alternate_cal_dir=None,filename=None):
    """returns None if cal file was not found
       this is used to find active calibration directory for this instrument
       and time"""

 # raise RuntimeError =[] for case when no file is found

 # routines requires both name and extension for match
    if filename is not None:
        return (filename,datetime(2100,1,1,0,0,0))

    return findCalFile(basedir,isntname,file_type,start_time,alternativePath=alternate_cal_dir)

class cal_vector_object():
    pass
def load_raman_tables(basedir,instname,consts,time=None,alternate_cal_dir=None):
    rs_cal = cal_vector_object()

    usegeofile=None 
    #usegeofile=locate_file('geofile_20150906T0830.geo')
    #usegeofile=locate_file('geofile_20150902T0759.geo')
    if usegeofile is not None:
        print 'WARNING using located geo file instead of data-archive version'
        import time
        time.sleep(15)
    geofile,geoexp=find_cal_file(basedir,instname,'geo',time,alternate_cal_dir=alternate_cal_dir,filename=usegeofile)
    #rs_cal.geo = calibration_vector_table(locate_file('geofile_raman_test.geo'),expire_time=None)
    rs_cal.geo = calibration_vector_table(geofile,expire_time=geoexp)
    print 'using RAMAN GEOFILE ' + geofile
    rs_cal.expire_time=rs_cal.geo.expire_time

    return rs_cal

def raman_quick_cal(consts,rs_cal,sounding):
    Cxx = cal_vector_object()   
    wavelength_elastic = consts['wavelength_elastic']
    wavelength_n2_raman = consts['wavelength_nitrogen_raman']
    """
    #get nitrogen Raman scattering cross section profile
    beta_raman_nitrogen = Raman_nitrogen_cross_section(
            wavelength_nitrogen,sounding.pressures,sounding.temperature,sounding.altitudes)
    """
   
    Cxx.altitudes = sounding.altitudes.copy()
    #get elastic backscatter and extinction
    Cxx.beta_r_355=lu.Rayleigh_cross_section(
         wavelength_elastic,sounding.pressures,sounding.temps,sounding.altitudes)
    #get n2 raman extinction
    Cxx.beta_r_387 = lu.Rayleigh_cross_section(
         wavelength_n2_raman,sounding.pressures,sounding.temps,sounding.altitudes)
    return Cxx

def process_raman(consts,mean,process_control,rs_cal,Cxx,corr_adjusts):
    """
       process_raman(consts,mean,process_control,sounding)

       operates on the output of range process after range to altitude conversion
       
      
    """

    
    print 'processing raman slice===================='                 
    inv = raman_inversion(mean,consts,Cxx,corr_adjusts,process_control)
    
    if 0:
        import matplotlib.pylab as plt
    	plt.figure(3001)
        plt.plot(np.nansum(rs.nitrogen_counts,0),raw.altitudes)
	ax = plt.gca()
        ax.set_xscale('log')
        plt.grid(True)
        plt.show()
    
    inv.optical_depth,inv.optical_depth_aerosol,inv.mol_norm_index ,inv.mol_ref_aod = lu.compute_optical_depth(
        mean.nitrogen_counts, Cxx.beta_r_355*3.0/(8*np.pi), mean.altitudes, process_control,consts)
    if 0:
        import matplotlib.pylab as plt
        plt.figure(4000)
        plt.plot(np.nanmean(inv.optical_depth,0),inv.altitudes)
	ax = plt.gca()
        ax.set_xscale('log')
        plt.grid(True)
        plt.show()
    return inv

from datetime import timedelta

def getRollingDescription(processing_defaults):
    import lg_dpl_toolbox.filters.dpl_rolling_window_filter as drwf
    window=processing_defaults.get_value('extinction_processing','time_window_length')
    order=processing_defaults.get_value('extinction_processing','polynomial_order')
    return drwf.WindowedFilterDescription(lsge.polyFilter,width=3,edgemode='fullduplicate',
        varlist=['filtered_nitrogen'],time_width=timedelta(seconds=window),cargs=(order,),
        windowinfo_perframe="ext_rollingwindow_debug")



def raman_inversion(mean,consts,Cxx,corr_adjusts,process_control):
    """raman_inversion(mean,consts,Cxx,corr_adjusts,process_control)
    """
    
    inv = hau.Time_Z_Group(like=mean)     
    inv.beta_r_355 = Cxx.beta_r_355 #FIXME
    inv.beta_r_387 = Cxx.beta_r_387
    beta_r_n2 = Cxx.beta_r_387

  
    inv.times=mean.times.copy()
    inv.delta_t=mean.delta_t.copy()
    inv.altitudes=mean.altitudes.copy()


    #add backscatter ratios to inv.
    adj = consts['elastic_to_n2_gain_ratio']
    wfov_adj = consts['wfov_elastic_to_n2_gain_ratio']
   
    inv.aerosol_backscatter_ratio = (mean.elastic_counts - mean.nitrogen_counts)\
           /mean.nitrogen_counts
    inv.aerosol_backscatter_ratio_low = (wfov_adj * mean.elastic_counts_low - mean.nitrogen_counts_low)\
          /mean.nitrogen_counts_low
    inv.aerosol_backscatter_ratio_high = (adj * mean.elastic_counts_high - mean.nitrogen_counts_high)\
          / mean.nitrogen_counts_high


    
    inv.beta_a_backscat_low =  3 * inv.aerosol_backscatter_ratio_low * Cxx.beta_r_355 / (8.0 * np.pi)
    inv.beta_a_backscat = 3 * inv.aerosol_backscatter_ratio * Cxx.beta_r_355 / (8.0 * np.pi)
    inv.integrated_backscatter = lu.integrated_backscatter(inv.beta_a_backscat, inv.altitudes)

    print
    print
    print 'check depol adjustment in raman_inversion***********************************************'
    print
    print
    
    inv.linear_depol = mean.depolarization_counts_high / (mean.elastic_counts_high - mean.nitrogen_counts_high)

    if 0:
      
        
        import matplotlib.pylab as plt
        bin_vec = np.arange(inv.altitudes.shape[0])
    
        plt.figure(5555)
        plt.plot(np.nanmean(inv.aerosol_backscatter_ratio_low,0),bin_vec,'r'
                 ,np.nanmean(inv.aerosol_backscatter_ratio_high,0),bin_vec,'m'
                 ,np.nanmean(inv.aerosol_backscatter_ratio,0),bin_vec,'k')
        plt.grid(True)
    
    #new code for raman extinction
 
    #filter_params needs to be called only when one of the calling parameters changes
    filter_params = lsge.filter_setup(inv.altitudes,process_control,consts,mean.delta_t[0])

    #compute range integrated backscatter cross section
    temp = inv.beta_a_backscat.copy()
    temp[np.isnan(temp)]= 0.0
    inv.integ_backscat = np.cumsum(temp,1) * filter_params['dz']
    
    #for key in  filter_params:
    #    print key ,' = ',filter_params[key]

    inv.extinction = np.NaN * np.zeros_like(mean.nitrogen_counts)
    ntimes = len(mean.nitrogen_counts[:,0])

    inv.extinction_aerosol = np.NaN * np.ones_like(mean.nitrogen_counts)
    inv.p180 = np.NaN * np.ones_like(mean.nitrogen_counts)
    if ntimes == 1: #dont care
        alreadyfiltered=hasattr(mean,'filtered_nitrogen')
        Nm=mean.filtered_nitrogen if alreadyfiltered else mean.nitrogen_counts
    elif hasattr(mean,'filtered_nitrogen') and not (mean.nitrogen_counts==mean.filtered_nitrogen).all():
        alreadyfiltered=True
        Nm=mean.filtered_nitrogen
        print "@@@@@@@@@ SKIPPING ROLLING FILTER. Already done"
    else:
        print '@@@@@ Trying bad ROLLING FILTER. Already done'
        alreadyfiltered=False
        Nm=mean.nitrogen_counts
    if not alreadyfiltered:
        try:
            half_slice = filter_params['t_window_pts']/2
        except:
            half_slice = 1
    
    sg_ext = lsge.sg_extinction(filter_params)

    if ntimes == 1:  #for profile call
      sl=np.arange(1)
      inv.extinction[0,:],inv.extinction_aerosol[0,:],inv.p180[0,:] \
              = sg_ext(mean.times[sl],mean.delta_t[sl],Nm[sl,:]
              ,inv.beta_a_backscat[0,:],inv.integ_backscat[0,:],beta_r=Cxx.beta_r_355)
    else:    
        #needs to add edge times with padded intervals
        state=[]
        if alreadyfiltered:
            fullrange=range(ntimes)
        else:
            fullrange=range(half_slice,ntimes-half_slice-1)
        for i in fullrange:
            if alreadyfiltered:
                sl=np.arange(i,i+1)#already filtered
            else:
                sl=np.arange(i-half_slice,i+half_slice+1)
            inv.extinction[i,:],inv.extinction_aerosol[i,:],inv.p180[i,:] \
              = sg_ext(mean.times[sl],mean.delta_t[sl],Nm[sl,:]
              ,inv.beta_a_backscat[i,:],inv.integ_backscat[i,:],beta_r=Cxx.beta_r_355,state=state)


    #end of new code for Raman exticntion           
  
    if 0:
        import matplotlib.pylab as plt
        plt.ion()
    	plt.figure(3002)
        plt.plot(np.nanmean(inv.extinction,0),inv.altitudes/1000.0)
	ax = plt.gca()
        ax.set_xscale('log')
        plt.xlabel('extinction (1/m)')
        plt.ylabel('altitude (km)')
        plt.grid(True)
    return inv


def process_ranged_raman(consts,rs_merge,process_control,rs_cal,corr_adjusts):

#def range_process(rs,consts,rs_cal,corr_adjusts,process_control):
    """
       range_process(rs,consts,rs_cal,corr_adjusts,process_control
       Do all processing that must be done on basis of  range bins
       rather than altitude.
    """
    print
    print
    print
    print
    print 'entering process_ranged_raman'

    rs = copy.deepcopy(rs_merge)

    if hasattr(rs,'ranged_runcount'):
        rs.ranged_runcount=rs.ranged_runcount+1
    else:
        rs.ranged_runcount=1
    print 'process on',rs.times[0],rs.times[-1],'run count is',rs.ranged_runcount
    import time
    time.sleep(3)
    
    binwidth = consts['binwidth']
    s_bin = np.int(consts['laser_pulse_timing'][1]/binwidth)
    s_bin = len(rs.heights[rs.heights <=0]) 
    nbins = len(rs.heights)
    
        
    #subtract dark counts determined prior to laser pulse for low channels
    #dictionary of channel_names and corresponding dark_count_names
    #dict(counts_name_in_structure_rs = 'dark_count_name_in_structure_rs')
    #                         
    select = dict(elastic_counts_low = 'elastic_low_dark_counts'
                  ,nitrogen_counts_low = 'nitrogen_low_dark_counts'
                  ,water_counts_low = 'water_low_dark_counts')          

    rs = lu.dark_count_correction_prefire(select,rs,corr_adjusts)


    if 0:
        bin_vec = np.arange(rs.nitrogen_counts_low[0,:].shape[0])
        import matplotlib.pylab as plt
        plt.figure(30)
        plt.plot(bin_vec,rs.nitrogen_counts_low[0,:],'b'
                 ,bin_vec,rs.elastic_counts_low[0,:],'c'
                 ,bin_vec,rs.elastic_counts_high[0,:],'r'
                 ,bin_vec,rs.depolarization_counts_high[0,:],'g'
                 ,bin_vec,rs.elastic_counts_high[0,:],'r')

    clear_first_bins(select,rs,process_control,consts)
    if 0:
        plt.plot(bin_vec,np.nanmean(rs.nitrogen_counts_low,0)+1)
        plt.show()
        print j
   

    #subtract dark counts for high channels
    #dict(counts_name_in_rs = 'dark_counts_name_to_added_to_rs'
   
    select = dict(nitrogen_counts_high = 'nitrogen_high_dark_counts'
              ,elastic_counts_high = 'elastic_high_dark_counts'
              ,water_counts_high = 'water_high_dark_counts'
              ,depolarization_counts_high = 'depolarization_high_dark_count')
    end_index = len(rs.nitrogen_counts_high[0,:])
    start_index = end_index - 10
    rs = lu.dark_count_correction_from_signal(select,rs,start_index,end_index)
    
  
    clear_first_bins(select,rs,process_control,consts)
      
    #select variables for geometry correction
    
    select = ['nitrogen_counts_high'
              ,'elastic_counts_high'
              ,'water_counts_high'
              ,'depolarization_counts_high']
    
    #do overlap correction for high channels using only r-sqrd correction
    rs = lu.geometry_correction(select,rs,rs_cal,nbins,s_bin,0)

    select = ['nitrogen_counts_low'
              ,'elastic_counts_low'
              ,'water_counts_low']
    
    #do only r-sqrd overlap correction for low channels
    rs = lu.geometry_correction(select,rs,rs_cal,nbins,s_bin,apply_geo_corr=1.0)
    
    #possible nitrogen channel gain adjustment
    rs.nitrogen_counts_high *= corr_adjusts['nitrogen_to_combined_gain']
    rs.nitrogen_counts_low  *= corr_adjusts['nitrogen_to_combined_gain']


    print dir(rs)

    depol_gain = consts['elastic_hi_to_depol_gain_ratio']
    rs.nitrogen_counts = rs.nitrogen_counts_high.copy()
    rs.elastic_counts = (rs.elastic_counts_high + rs.depolarization_counts_high / depol_gain).copy()  
    rs.water_counts = rs.water_counts_high.copy()
   
    if process_control.enabled('use_wfov'):
        #replace nfov data below splice range with wfov data
        splice_bin = np.float(process_control.get_value('use_wfov','splice_range'))*1000.0
        splice_bin =np.int( splice_bin/ (binwidth * 1.5e8)) +s_bin#1.5e8 = 1/2 speed of light)
        elastic_gain_ratio,n2_gain_ratio,h2o_gain_ratio = hi_lo_channel_gain_ratios(rs,consts)
        fm = np.ones_like(rs.filter_mode).astype('bool')
        #fm = rs.filter_mode == 1
        if 0:
            import matplotlib.pylab as plt
            bin_vec = np.arange(rs.nitrogen_counts.shape[1])
            plt.figure(1000)
            plt.plot(bin_vec,np.nanmean(rs.nitrogen_counts,0),'r',bin_vec,np.nanmean(rs.nitrogen_counts_low,0)*50.0,'g')
            plt.grid(True)

        rs.nitrogen_counts[fm,:splice_bin] = (rs.nitrogen_counts_low[fm,:splice_bin]\
                     * n2_gain_ratio[fm,:splice_bin])

        if 0:
            import matplotlib.pylab as plt
            bin_vec = np.arange(rs.nitrogen_counts.shape[1])
            plt.figure(1001)
            plt.plot(bin_vec,np.nanmean(rs.nitrogen_counts,0))
            plt.grid(True)

        #adjustments to elastic gain required to make zero scat ratio in clear layers
        adj = consts['elastic_to_n2_gain_ratio']
        wfov_adj = consts['wfov_elastic_to_n2_gain_ratio']
        rs.elastic_counts = adj * rs.elastic_counts 
        rs.elastic_counts[:,:splice_bin] = rs.elastic_counts_low[:,:splice_bin]\
                     * elastic_gain_ratio[:,:splice_bin] * wfov_adj 

        rs.water_counts[fm,:splice_bin] = (rs.water_counts_low[fm,:splice_bin]\
                     * h2o_gain_ratio[fm,:splice_bin])
  
    rs.filtered_nitrogen=rs.nitrogen_counts.copy() #this is copied so a rolling filter will filter it
    
    print 'leaving range_process'
    print

   
    return rs


def clear_first_bins(select,rs,process_control,consts):
    first_bin = process_control.get_value('first_bin_to_process','bin_number')
    first_bin = first_bin + consts['data_bin_containing_laser_pulse']
    for name in select:
         var = getattr(rs,name)
         var[:,:first_bin] = np.NaN
         setattr(rs,name,var)
    return

def hi_lo_channel_gain_ratios(rs,consts):
    """
       hi_lo_channel_gain_ratio(rs,consts)
       provides ratio between high and low gain channels based on filter mode
    """

    ones_array = np.ones_like(rs.nitrogen_counts_high)
    nbins_low = len(rs.nitrogen_counts_low[0,:])
    elastic_gain_ratio = consts['elastic_hi_to_lo_gain_ratio'][1] \
              *ones_array
    elastic_gain_ratio[rs.filter_mode == 2,:nbins_low] = consts['elastic_hi_to_lo_gain_ratio'][2]\
                         * ones_array[rs.filter_mode == 2,:nbins_low]
    n2_gain_ratio = consts['nitrogen_hi_to_lo_gain_ratio'][1] * ones_array
    n2_gain_ratio[rs.filter_mode == 2,:nbins_low] = consts['nitrogen_hi_to_lo_gain_ratio'][2] \
                         * ones_array[rs.filter_mode==2,:nbins_low]
    h2o_gain_ratio = consts['water_hi_to_lo_gain_ratio'][1] * ones_array
    h2o_gain_ratio[rs.filter_mode == 2,:nbins_low] = consts['water_hi_to_lo_gain_ratio'][2] \
                         * ones_array[rs.filter_mode == 2,:nbins_low]

    return elastic_gain_ratio,n2_gain_ratio,h2o_gain_ratio


