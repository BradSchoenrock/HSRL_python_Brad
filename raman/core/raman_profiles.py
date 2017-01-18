import lg_base.core.array_utils as hau
import numpy as np
from datetime import datetime,timedelta
from lg_base.core.accumulation import accumulate

energies= ('pulse_energy',)

channel_shorthand=dict(nitrogen_counts= 'n2',nitrogen_counts_high= 'n2_high',nitrogen_counts_low= 'n2_low'
                           ,elastic_counts= 'e',elastic_counts_high= 'e_high',elastic_counts_low= 'e_low'
                           ,water_counts= 'h20',water_counts_high='h20_high', water_counts_low= 'h20_low'
                           ,depolarization_counts_high= 'cpol')

def accumulate_raman_inverted_profiles(consts,rs_mean,qc_mask,old_profiles,process_control=None,corr_adjusts=None):
    indices=np.arange(rs_mean.times.shape[0])
    if len(indices)==0:
      return old_profiles
  

    if qc_mask is not None and processing_defaults.get_value('averaged_profiles','apply_mask'):    
        #make mask array with NaN's for array elements where bit[0] of qc_mask==0 
        #all other elements of mask = 1
        mask = (np.bitwise_and(qc_mask,1)).astype('float')        #mask is float to allow use of NaN values
        mask[mask == 0] = np.NaN
       
        print 'qc_mask applied to time averaged profiles vs altitude'
    else:
        #set mask == 1
        mask = None
        #for sh in channel_shorthand.keys():
        #
        # if hasattr(rs_mean,sh):
    #        mask = np.ones_like(getattr(rs_mean,sh))
    #        break
        #if mask is None:
        #    mask = np.ones((rs_mean.times.shape[0],0)) 
        print 'qc_mask has not been applied to time averaged profiles'

    #energies=('transmitted_1064_energy','transmitted_energy')
  
    profiles = hau.Time_Z_Group(can_append=False,altname='altitudes')
    

    profiles.hist=hau.Time_Z_Group()
    ft=None
    tc=0
    #tv=0
    lt=None
    if old_profiles is not None:
      ft=old_profiles.hist.ft
      tc=old_profiles.hist.tc
      #tv=old_profiles.hist.tv
      lt=old_profiles.hist.lt
    if len(indices)>0:
      if ft is None:
        ft=rs_mean.times[indices][0]          
      lt=rs_mean.times[indices][-1]+timedelta(seconds=rs_mean.delta_t[indices][-1] if not np.isnan(rs_mean.delta_t[indices][-1]) else 0)                  
      for x in indices:
          if rs_mean.times[x] is not None:
              tc=tc+1;
              #tv=tv+(rs_mean.times[x]-ft).total_seconds()
    if tc>0:
      profiles.times = hau.T_Array([ft])#+timedelta(seconds=tv/tc), ])
      profiles.start = ft
      profiles.width = lt-ft
      profiles.delta_t = hau.T_Array([profiles.width.total_seconds()])
    else:
      profiles.times=hau.T_Array([])
      profiles.start = ft
      profiles.width = timedelta(seconds=0)
      profiles.delta_t = hau.T_Array([])
    profiles.start_time = ft
    profiles.end_time = lt
    profiles.hist.ft=ft
    profiles.hist.tc=tc
    #profiles.hist.tv=tv
    profiles.hist.lt=lt
    if rs_mean is not None and hasattr(rs_mean,'altitudes'):
      profiles.altitudes = rs_mean.altitudes.copy()
    elif hasattr(old_profiles,'altitudes'):
        profiles.altitudes = old_profiles.altitudes

    #FIXME need to accumulate the inverted products here, so here's a hack
    interval=hau.Time_Z_Group()
    interval.intervals=hau.T_Array(np.ones(rs_mean.times.shape))

    accumulate(profiles,old_profiles,interval,indices,'intervals')
    interval_count=profiles.intervals
    print 'Total intervals for profile =',interval_count
    for k,v in vars(rs_mean).items():
      if k.startswith('_') or k in ('times','start','width','delta_t'):
        continue
      if not isinstance(v,hau.T_Array):
        continue
      if isinstance(v,hau.TZ_Array):
        continue
      accumulate(profiles,old_profiles,rs_mean,indices,k,interval_count)
    # create TZ_Array with time dimension of '1', so hsrl_inversion doesn't choke
    for k,v in vars(rs_mean).items():
      if k.startswith('_'):
        continue
      if not isinstance(v,hau.TZ_Array):
        continue
      if len(v.shape)!=2:
        continue
      accumulate(profiles,old_profiles,rs_mean,indices,k,interval_count,mask)


    return profiles
    
def accumulate_raman_profiles(consts,rs_mean,qc_mask,old_profiles,process_control=None,rs_cal=None,Cxx=None,corr_adjusts=None):
    indices=np.arange(rs_mean.times.shape[0])
    if len(indices)==0:
      return old_profiles
  

    if qc_mask is not None and processing_defaults.get_value('averaged_profiles','apply_mask'):    
        #make mask array with NaN's for array elements where bit[0] of qc_mask==0 
        #all other elements of mask = 1
        mask = (np.bitwise_and(qc_mask,1)).astype('float')        #mask is float to allow use of NaN values
        mask[mask == 0] = np.NaN
       
        print 'qc_mask applied to time averaged profiles vs altitude'
    else:
        #set mask == 1
        mask = None
        #for sh in channel_shorthand.keys():
        #
        #	if hasattr(rs_mean,sh):
		#        mask = np.ones_like(getattr(rs_mean,sh))
		#        break
        #if mask is None:
        #    mask = np.ones((rs_mean.times.shape[0],0)) 
        print 'qc_mask has not been applied to time averaged profiles'

    #energies=('transmitted_1064_energy','transmitted_energy')
  
    profiles = hau.Time_Z_Group(can_append=False,altname='altitudes')
    

    profiles.hist=hau.Time_Z_Group()
    ft=None
    tc=0
    #tv=0
    lt=None
    if old_profiles is not None:
      #total_seeded_shots=total_seeded_shots+profiles.hist.total_seeded_shots
      ft=old_profiles.hist.ft
      tc=old_profiles.hist.tc
      #tv=old_profiles.hist.tv
      lt=old_profiles.hist.lt
    if len(indices)>0:
      if ft is None:
        ft=rs_mean.times[indices][0]          
      lt=rs_mean.times[indices][-1]+timedelta(seconds=rs_mean.delta_t[indices][-1] if not np.isnan(rs_mean.delta_t[indices][-1]) else 0)                  
      for x in indices:
          if rs_mean.times[x] is not None:
              tc=tc+1;
              #tv=tv+(rs_mean.times[x]-ft).total_seconds()
    if tc>0:
      profiles.times = hau.T_Array([ft])#+timedelta(seconds=tv/tc), ])
      profiles.start = ft
      profiles.width = lt-ft
      profiles.delta_t = hau.T_Array([profiles.width.total_seconds()])
    else:
      profiles.times=hau.T_Array([])
      profiles.start = ft
      profiles.width = timedelta(seconds=0)
      profiles.delta_t = hau.T_Array([])
    profiles.start_time = ft
    profiles.end_time = lt
    profiles.hist.ft=ft
    profiles.hist.tc=tc
    #profiles.hist.tv=tv
    profiles.hist.lt=lt
    #profiles.hist.total_seeded_shots=total_seeded_shots
    if rs_mean is not None and hasattr(rs_mean,'altitudes'):
      profiles.altitudes = rs_mean.altitudes.copy()


    #elif hasattr(old_profiles,'altitudes'):
    #  profiles.altitudes=old_profiles.altitudes
    elif hasattr(old_profiles,'heights'):
        profiles.heights = old_profiles.heights

    accumulate(profiles,old_profiles,rs_mean,indices,'shots',pref='mean_',filler=hau.T_Array([0]))
    total_shots=profiles.mean_shots
    profiles.shots=total_shots.copy()
    print 'Total shots for profile =',total_shots
    for e in energies:
      accumulate(profiles,old_profiles,rs_mean,indices,e,total_shots)
    # create TZ_Array with time dimension of '1', so hsrl_inversion doesn't choke
    for chan in channel_shorthand.keys():  
      accumulate(profiles,old_profiles,rs_mean,indices,chan,total_shots,mask)
   
    #compute inverted profiles from mean count profiles
    if Cxx is not None and hasattr(profiles,'elastic_counts'):
        import raman.core.raman_processing_utilities as rpu
        profiles.inv = rpu.process_raman(consts,profiles,process_control,rs_cal,Cxx,corr_adjusts)
        import lidar.sg_extinction as lsge
        filter_params = lsge.filter_setup(profiles.inv.altitudes,process_control,consts)
        sg_ext = lsge.sg_extinction(filter_params)
        profiles.inv.extinction,profiles.inv.extinction_aerosol,profiles.inv.p180 \
              = sg_ext(profiles.times,profiles.delta_t,profiles.nitrogen_counts ,profiles.inv.beta_a_backscat,profiles.inv.integ_backscat,beta_r=Cxx.beta_r_355)

    elif hasattr(old_profiles,'inv'):
        profiles.inv=old_profiles.inv
    return profiles
    
