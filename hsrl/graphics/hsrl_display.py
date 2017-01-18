#usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import hsrl.filters.savitzky_golay as sg
# Try to use the much faster nanmean from bottleneck, otherwise fall back
# to the scipy.stats version
try:
    from bottleneck import nanmean, nanmax
except ImportError:
    print
    print "No bottleneck.nanmean available! Falling back to SLOW scipy.stats.nanmean"
    print
    from scipy.stats import nanmean

import lg_base.core.read_utilities as hru
import lg_base.graphics.graphics_toolkit as gt
import lg_base.core.json_config as jc
import traceback
from operator import itemgetter
import datetime
import json
import string
import copy
from lg_base.core.open_config import open_config
from lg_base.core.locate_file import locate_file


def getframes(v):
  ret=[]
  for x in vars(v).keys():
    if x.startswith("_"):
      continue
    ret.append(x)
  return ret

def get_display_defaults(display_defaults_file,filetype="new"):
    """get display defaults from display_default.json.
       image selections"""
    
    if filetype == "old":
        raise RuntimeError('Update the code. "old" defaults file is deprecated')
        fd = open_config(display_defaults_file)
        dd = json.load(fd)
        display_defaults = dd['display_defaults']
        config = dd['config']
        fd.close()
        return (display_defaults, config)
    
    elif filetype == "new":
        print 'hsrl.graphics.hsrl_display.get_display_defaults() is deprecated. use lg_base.core.json_config.get_display_defaults()'
        return jc.get_display_defaults(display_defaults_file)


def define_multiple_alt_plots(var,altitudes,indices,lines,widths,colors,legends):
    """
       [lines,colors,widths,legends]=define_multiple_alt_plots(var,altitudes \
                        ,indices,lines,colors,legends)
       provides list of variables at altitudes[indices]
       var = variable to plot at altitudes[indices[:]]
       altitudes = altitude vector corresponding to var levels
       indices   = list of altitude indices to plot
       lines     = list of variables at requested levels
    """
  
    if len(lines) == 0:
        width = 2
    else:
        width = 1
    color_list = ['r','b','g','k','c','m','r','b','g','k','c','m']
    if not isinstance(indices,(list,tuple,np.ndarray)):
      indices=[indices]
   
    for i in range(len(indices)):
       
        lines.append(var[:,indices[i]])
        widths.append(width)
        colors.append(color_list[i])
        legends.append('%4.2f km'%(altitudes[indices[i]]/1000.0))
                   
    return lines,colors,widths,legends

def multiple_parameter_mask(rs,display_defaults):
    mask = np.ones(rs.beta_a_backscat.shape,dtype=bool)
    if not display_defaults.enabled('multiple_parameter_mask'):
        #no masking if multiple-parameter_mask is not found
        return mask
    else:
        #check for mask limits on the following:
             
        for name in ['beta_a_backscat','linear_depol'
                      ,'color_ratio','SN_beta_a_backscat'
                      ,'SN_mol']:
            if not display_defaults.get_value('multiple_parameter_mask',name)== None and hasattr(rs,name):
                limits = display_defaults.get_value('multiple_parameter_mask',name)
                mask[getattr(rs,name) < limits[0]] = 0 
                mask[getattr(rs,name) > limits[1]] = 0
                print "multiple parameter mask ",name,' set with limits = ',limits
               
    return mask

def show_images( instrument, rs, sounding, rs_constants,geo_corr,rs_Cxx,processing_defaults
                 ,display_defaults,last_sounding_time, max_alt, auto_loop, enable_masking=None, figlist=None):
   
    figs=figlist   #figure objects keyed to a unique name
    if figs is None:
         figs=gt.figurelist()

  
    if True:#display_defaults.get_value('select_plot_altitude','plot_altitude')!=None:
      getaltsfrom=None
      for f in ('rs_inv','rs_mean','profiles','rs_raw'):
        if hasattr(rs,f):
          print 'got frame ',f
          getaltsfrom=getattr(rs,f)
          break
      if hasattr(getaltsfrom,'msl_altitudes'):
        invalts=getaltsfrom.msl_altitudes
        shared_dz=invalts[2]-invalts[1]        
      else:
        invalts=np.array([0,max_alt])
        shared_dz=None
     
      #print 'PLOT ALT INDEX is ',getaltsfrom.plot_alt_index if hasattr(getaltsfrom,'plot_alt_index') else 'empty'
      #print 'LAYER INDICES is ',getaltsfrom.layer_indices if hasattr(getaltsfrom,'layer_indices') else 'empty'
      plot_alt_index = getaltsfrom.plot_alt_index if hasattr(getaltsfrom,'plot_alt_index') else 0
      layer_indices = getaltsfrom.layer_indices if hasattr(getaltsfrom,'layer_indices') else np.arange(len(invalts))
      if invalts is not None:
          if len(layer_indices)==0 and len(invalts) > 2:
              print
              print "WARNING--display layer average indices don't overlap data altitudes"
              print "        layer reset to include all available altitudes"
              print
              layer_indices = np.arange(len(invalts))

          if (not ('installation' in rs_constants) or rs_constants['installation'] == 'ground')\
               and display_defaults.enabled('image_altitude_above_ground_level'):
             invalts = invalts - rs_constants['lidar_altitude']
      if not isinstance(plot_alt_index,(list,tuple,np.ndarray)):
          plot_alt_index=[plot_alt_index]
      
      
      if hasattr(getaltsfrom,'start'):
        print 'for timerange, using start and width'
        timerange=[getattr(getaltsfrom,'start'),None]
        try:
            timerange[1]=timerange[0]+getattr(getaltsfrom,'width')
        except TypeError:
            try:
                print 'start is not a real time.'
                timerange[1]=getaltsfrom.times[0]+getaltsfrom.width
            except:
                print 'NO TIMES PFFF'
                timerange[1]=timerange[0]
        timerange=tuple(timerange)
      elif hasattr(getaltsfrom,'times') and getaltsfrom.times.shape[0]>0:
        if hasattr(getaltsfrom,'delta_t') and np.isfinite(getaltsfrom.delta_t[-1]):
          print 'for timerange, using times and delta_t'
          timerange=(getaltsfrom.times[0],getaltsfrom.times[-1]+datetime.timedelta(seconds=getaltsfrom.delta_t[-1]))
        else:
          print 'for timerange, using times'
          timerange=(getaltsfrom.times[0],getaltsfrom.times[-1])
      else:
        print "couldn't find a time range"
        timerange=None
    haveProfiles = hasattr(rs,'profiles')
    if display_defaults.enabled('sounding') and sounding is not None and hasattr(sounding,'times'):# and sounding.times != last_sounding_time:
       
        indices = np.arange(len(sounding.altitudes))
        top_index = np.max(indices[sounding.altitudes <= max_alt])
        top_raob_index = np.min([np.max(indices[sounding.altitudes <= sounding.top])
                                , top_index ])
        t_min = np.min(sounding.dew_points[:top_raob_index])-10
        t_max = np.max(sounding.temps[:top_index])+10
       
        if np.isnan(t_min):
            t_min =220
            print ' '
            print 'Display: error in t_min, setting plot min temp to 221 K'
        if np.isnan(t_max):
            t_max=299
            print ' '
            print 'Display: error in t_max, setting plot max temp to 299 K'
        t_min = int(t_min/10)*10
        t_max = int(t_max/10)*10
        if t_min < 233.15:
            t =[sounding.temps[:top_index]
                ,sounding.temps[:top_raob_index]
                ,sounding.dew_points[:top_raob_index]
               ,sounding.frost_points[:top_raob_index]  
                ,np.array([273.15, 273.15])
                ,np.array([233.15, 233.15])
                ,np.array([t_min, t_max])]
            a =[sounding.altitudes[:top_index]/1000
                ,sounding.altitudes[:top_raob_index]/1000
                ,sounding.altitudes[:top_raob_index]/1000
                ,sounding.altitudes[:top_raob_index]/1000
                ,np.array([0,sounding.altitudes[top_index]/1000])
                ,np.array([0,sounding.altitudes[top_index]/1000])
                ,np.array([sounding.bot/1000, sounding.bot/1000])]
        else:  #don't print -4OC line if no cold temps in sounding
            t =[sounding.temps[:top_index]
                ,sounding.temps[:top_raob_index]
                ,sounding.dew_points[:top_raob_index]
                ,sounding.frost_points[:top_raob_index]  
                ,np.array([273.15, 273.15])
                ,np.array([t_min, t_max])]
            a =[sounding.altitudes[:top_index]/1000
                ,sounding.altitudes[:top_raob_index]/1000
                ,sounding.altitudes[:top_raob_index]/1000
                ,sounding.altitudes[:top_raob_index]/1000
                ,np.array([0,sounding.altitudes[top_index]/1000])
                ,np.array([sounding.bot/1000, sounding.bot/1000])] 
        #compute potential temperature for highest point displayed in sounding
        t1000 = sounding.temps[top_index]*(1000.0/sounding.pressures[top_raob_index])**0.286 
        dt=(t_max-t_min)/10
        colors =[]
        linewidth=[]
        linetype =[]
        temps=[]
        alts=[]
      
        while np.isfinite(t1000) and t1000 > t_min:    
            adiabat =t1000*(1000./sounding.pressures[:top_index])**-0.286
            inx=indices[adiabat >t_min ]
            inx=inx[adiabat[inx] < t_max]
            inx=inx[sounding.altitudes[inx]>sounding.bot]
            alts.append(sounding.altitudes[inx]/1000.0)
            temps.append(adiabat[inx])
            colors.append('c')
            linewidth.append(1)
            linetype.append('-')
            t1000 = t1000 -dt
            
        temps=temps + t
        alts=alts + a
        colors=colors +['k','r','g','b','k','k','g']
        linetype=linetype +['-','-','-','-','-','-','-']
        linewidth=linewidth +[1,3,3,3,1,1,3]
        text_str=['temp(red)', 'dewpt(grn)','  frostpt(blue)',' adiabats(c)']
        text_position_x=[t_max-15
                        ,t_max-15
                        ,t_max-15
                        ,t_max-15]
        text_position_y=[0.95*max_alt/1000
                        ,0.9*max_alt/1000
                        ,0.85*max_alt/1000
                        ,0.8*max_alt/1000]
        text_angle=[0,0,0,0]
       
        gt.plot_xy('sounding'  #plot name
                ,sounding.station_id
                ,sounding.times
                ,temps
                ,alts
                ,colors
                ,None
                ,None
                ,linetype
                ,linewidth
                ,None
                ,'upper right'
                ,'Temperature'
                ,'deg-K'
                ,'MSL altitude'
                ,'km'
                ,'Sounding'
                ,text_str
                ,text_position_x
                ,text_position_y
                ,text_angle
                ,display_defaults
                ,figs)
       
    #plot image of qc_mask state
    #states plotted:
    #       0-no mask
    #       1-masked because I2 line locking is lost
    #       2-masked because mol signal below threshold counts
    #         or mol signal/noise below threshhold
    #       3-masked because cloud threshhold is exceeded
    #       4-masked because both mol low and I2 unlocked
    #       5-masked because both cloud threshold and mol lost
    #       6-masked because all of above occured    
    if display_defaults.enabled('qc_mask_field_image') and hasattr(rs,'rs_inv') and hasattr(rs.rs_inv,'qc_mask'): 
           title_str = 'qc_mask'
           qc_mask = rs.rs_inv.qc_mask
           if qc_mask.dtype!='uint32':
              nqc_mask=np.zeros_like(qc_mask,dtype='uint32')
              nqc_mask[np.isfinite(qc_mask)]=qc_mask[np.isfinite(qc_mask)]
              qc_mask=nqc_mask
              del nqc_mask
           temp = np.zeros_like(qc_mask,dtype='uint32')
           #note: both mol_lost and mol_sn_threshhold activate mol_lost in figure
           #if bit 0 is set none of the masks are active
           #temp[rs.rs_inv.qc_mask &1] = 0
           
           
           #mask=16 when only mol < mol_sn_noise_threshhold
           temp[np.bitwise_and(qc_mask,16)==0]=2
    
           #mask 64 when only mol_lost is active
           temp[np.bitwise_and(qc_mask,64)==0] = 2
           
           #mask=128 when only cloud mask is active
           temp[np.bitwise_and(qc_mask,128)==0]=3
           
           #mask=69 when mol_lost and i2_lock are active
           temp[np.bitwise_and(qc_mask,68)==0]=4

           #temp[mask == 20] = 4
           temp[np.bitwise_and(qc_mask,20)==0]=4

           #mask=192 when cloud_mask and mol_lost are active
           temp[np.bitwise_and(qc_mask,192)==0]=5
           temp[np.bitwise_and(qc_mask,144)==0]=5
           
           #mask=196 when i2_lock, cloud_mask and mol_lost are all active
           #temp[np.bitwise_and(qc_mask,196)==0]=6
           #temp[np.bitwise_and(qc_mask,146)==0]=6
           #temp[mask == 146] = 6

           #when only i2_lock warning is present (red)
           temp[np.bitwise_and(qc_mask,2**13+1)==1 ] = 1

           #mask=4 when  i2 lock mask is active (black)
           #temp[mask == 4] = 1
           temp[np.bitwise_and(qc_mask,2**2)==0]=6
           
           cb_labels=[gt.ClassColormapEntry(0,'no mask',(0,1,0)),   #green
                      gt.ClassColormapEntry(1,'I2 warn',(0,0.5,0)),  #dark green
                      gt.ClassColormapEntry(2,'m_lost','purple'),    #purple
                      gt.ClassColormapEntry(3,'cloud', 'blue'),
                      gt.ClassColormapEntry(4,'m_lost,i2',(.3,0,.3)),  #dark purple
                      gt.ClassColormapEntry(5,'cld,mol','darkcyan'),
                      gt.ClassColormapEntry(6,'i2_lock','black')]
           gt.rti_fig('qc_mask_field_image'
               ,instrument    
               ,temp
               ,rs.rs_inv.times
               ,invalts
               ,title_str               
               ,cb_labels
               ,None       
               ,display_defaults
               ,figs)
           del qc_mask #this is only used for the graph. masking is determined below

    if display_defaults.enabled('qa_mask_images'): 
        qa_flags=rs
        if hasattr(rs,'rs_inv') and hasattr(rs.rs_inv,'qa_flags'):
          qa_flags=rs.rs_inv.qa_flags
        flags=dict(qa_BS='Backscatter',qa_dep='Depolarization',qa_Ext='Extinction')
        for k,d in flags.items():
           if not hasattr(qa_flags,k):
              continue
           title_str = d + ' Quality Assurance Flag'
           qa_mask = getattr(qa_flags,k)
           
           cb_labels=[gt.ClassColormapEntry(0,'Not Checked',(0.5,0.5,0.0)), #dark yellow 
                      gt.ClassColormapEntry(2,'Bad',(0.5,0,0)),  #dark red
                      gt.ClassColormapEntry(1,'Good',(0,0.5,0)),  #dark green
                      gt.ClassColormapEntry(3,'Caution',(1.0,0.5,0))]  #dark orange
           gt.rti_fig('qa_mask_images'
               ,instrument    
               ,qa_mask
               ,qa_flags.times if hasattr(qa_flags,'times') else rs.rs_inv.times
               ,qa_flags.altitudes if hasattr (qa_flags,'altitudes') else invalts
               ,title_str               
               ,cb_labels
               ,None       
               ,display_defaults
               ,figs,fig_name=k+'_mask_image')

    qc_mask=None #none can be passed safely
    qc_mask_1064=qc_mask #dont need a copy. they're both dummies
    if ((enable_masking is None and display_defaults.enabled('mask_image')) or enable_masking):
         if hasattr(rs,'rs_inv') and hasattr(rs.rs_inv,'qc_mask'):
             tqcm=rs.rs_inv.qc_mask
         elif hasattr(rs,'rs_mean') and hasattr(rs.rs_mean,'qc_mask'): 
             tqcm=rs.rs_mean.qc_mask
         else:
             tqcm=None
             print 'Masking is enabled, but no mask. available subframes are ',getframes(rs)
         #select bits from qc_mask to mask images 
         if tqcm is not None:
               qc_mask=np.ones_like(rs.rs_inv.qc_mask)
               tqcm=rs.rs_inv.qc_mask
               qcbits={'mol_lost':64,'mol_sn_ratio':16,'particulate_backscat_SN_ratio':32,'cloud_mask':128,'I2_lock_lost':4}
               qcbits_1064={'1064_shutter':0x8000}

               if tqcm.dtype!='uint32':
                  ntqcm=np.zeros_like(tqcm,dtype='uint32')
                  ntqcm[np.isfinite(tqcm)]=tqcm[np.isfinite(tqcm)]
                  tqcm=ntqcm
                  del ntqcm

               for name,maskbit in qcbits.items():
                 print display_defaults.get_value('mask_image',name)
                 if display_defaults.get_value('mask_image',name):
                   qc_mask = np.logical_and(tqcm & maskbit > 0,qc_mask)
               qc_mask_1064=qc_mask.copy()
               for name,maskbit in qcbits_1064.items():
                 if display_defaults.get_value('mask_image',name):
                   qc_mask_1064 = np.logical_and(tqcm & maskbit > 0,qc_mask_1064)
    
         if hasattr(rs,'rs_inv'):  
             p_mask = multiple_parameter_mask(rs.rs_inv,display_defaults)
             qc_mask = qc_mask * p_mask
             qc_mask_1064 = qc_mask_1064 * p_mask



    if display_defaults.enabled('feature_mask') and hasattr(rs,'rs_inv') \
            and hasattr(rs.rs_inv,'color_ratio'):
        title_str = 'feature mask'
        mask_true =np.full( rs.rs_inv.beta_a_backscat.shape, True, dtype=bool )
        mask = mask_true.copy()
        temp = np.zeros_like(rs.rs_inv.beta_a_backscat,dtype='uint32')
        
        np.seterr(invalid = 'ignore')
        #dust
        mask = mask_true.copy()
        mask[rs.rs_inv.linear_depol < 0.15] = False
        mask[rs.rs_inv.color_ratio < 0.3] = False
        temp[mask] = 2 #set cirrus pixels
       

        #cirrus
        mask = mask_true.copy()
        mask[rs.rs_inv.beta_a_backscat < 3e-6] = False
        mask[rs.rs_inv.linear_depol < 0.3] = False
        mask[rs.rs_inv.color_ratio < 0.3] = False
        mask[rs.rs_inv.beta_a_backscat > 1e-4] = True #dense must be cloud
        temp[mask] = 1 #set cirrus pixels
        
     
        #water cloud
        mask = mask_true.copy()
        mask[rs.rs_inv.beta_a_backscat < 1e-4] = False
        mask[rs.rs_inv.linear_depol > 0.05] = False
        mask[rs.rs_inv.color_ratio < 0.6] = False
        temp[mask] = 3     
      

        #spherical coase mode
        mask = mask_true.copy()
        mask[rs.rs_inv.beta_a_backscat >= 1e-4] = False
        mask[rs.rs_inv.linear_depol > 0.05] = False
        mask[rs.rs_inv.color_ratio < 0.6] = False
        temp[mask] = 4
       

        #fine mode spherical aerosol
        mask = mask_true.copy()
        #mask[rs.rs_inv.beta_a_backscat < 1e-4] = False
        mask[rs.rs_inv.linear_depol > 0.05] = False
        mask[rs.rs_inv.color_ratio  >= 0.6] = False
        temp[mask] = 5 #set spherical aerosol pixels
     

        #coarse mode spherical + dust aerosol
        mask = mask_true.copy()
        mask[rs.rs_inv.linear_depol < 0.05] = False
        mask[rs.rs_inv.linear_depol > 0.15] = False
        temp[mask] = 6 #set spherical aerosol pixels
       

        #select fine mode portion of spherical + dust aerosol
        mask = mask_true.copy()
        mask[rs.rs_inv.linear_depol < 0.05] = False
        mask[rs.rs_inv.color_ratio >= 0.3] = False
        temp[mask] = 7 #set spherical aerosol pixels
       
        np.seterr(invalid = 'warn')


        # pixels with low aerosol S/N not classified
        SN_threshold = display_defaults.get_value('feature_mask','SN_beta_a_backscat_threshold')
	temp[rs.rs_inv.SN_beta_a_backscat < SN_threshold] = 0
      
        temp[np.isnan(temp)] = 0
       
        cb_labels=[gt.ClassColormapEntry(0,' ',(0,0,0)),   #undefined
                      gt.ClassColormapEntry(1,'ice',(1,0,0)),  #bright red
                      gt.ClassColormapEntry(2,'dust',(.7,0,0)),  #dark red
                      gt.ClassColormapEntry(3,'water',(1,1,1)), #white
                      gt.ClassColormapEntry(4,'cm s',(0,0,1)), #bright blue
                      gt.ClassColormapEntry(5,'fm s',(0,0,0.7)), #darker blue
                      gt.ClassColormapEntry(6,'cm s+d',(0,0.9,0)), #green
                      gt.ClassColormapEntry(7,'fm s+d',(0,0.5,0))] #dark green
        print 'rendering ',title_str
        if 1: #try:
            gt.rti_fig('feature_mask'
                    , instrument
                    , temp
                    , rs.rs_inv.times
                    , invalts
                    , title_str
                    , cb_labels
                    , qc_mask
                    , display_defaults
                    ,figs)
        if 0: #except AttributeError:
                    raise RuntimeError, "show_images - no data for feature id  plot"



    if display_defaults.enabled('atten_backscat_image') and hasattr(rs,'rs_inv') and hasattr(rs.rs_inv,'atten_beta_a_backscat'):

      # attenuated backscatter is displayed without mask

        title_str = instrument + '  atten backscatter cross section '
      
        gt.rti_fig('atten_backscat_image'
                ,instrument
                , rs.rs_inv.atten_beta_a_backscat
                , rs.rs_inv.times 
                , invalts
                , title_str
                , '1/(m sr)'
                , None             
                , display_defaults
                ,figs)
   
    if display_defaults.enabled('temperature_change') and hasattr(rs,'rs_inv') and hasattr(rs.rs_inv,'temps'):

        print 'rendering temperature change plot'

        title_str = instrument + '  sounding delta_t (deg C) '
       
        ones_array = np.ones_like(rs.rs_inv.temps)

        gt.rti_fig('temperature_change'
                ,instrument
                , rs.rs_inv.temps - rs.rs_inv.temps[0,:] * ones_array
                , rs.rs_inv.times 
                , invalts
                , title_str
                , ''
                , None             
                , display_defaults
                ,figs)
        
    if display_defaults.enabled('atten_IR_backscat_image') \
            and hasattr(rs,'rs_inv')\
            and hasattr(rs.rs_inv,'attenuated_1064_backscatter'):

        title_str = instrument + ' attenuated 1064 backscatter  '
        units_str= ''
        
        gt.rti_fig('atten_IR_backscat_image'
            ,instrument    
            ,rs.rs_inv.attenuated_1064_backscatter
            ,rs.rs_inv.times
            ,invalts
            ,title_str
            ,units_str
            ,qc_mask_1064
            ,display_defaults
            ,figs)
        
    if display_defaults.enabled('1064_backscat_image') and hasattr(rs,'rs_inv')\
            and hasattr(rs.rs_inv,'beta_a_1064_backscat'):     
        # particulate backscatter cross section   
        title_str = instrument + '  1064 backscatter, A='\
                    +str(processing_defaults.get_value('color_ratio','angstrom_coef'))

        gt.rti_fig('1064_backscat_image'
            ,instrument    
            ,rs.rs_inv.beta_a_1064_backscat
            ,rs.rs_inv.times
            ,invalts
            ,title_str
            ,'1/(m sr)'
            ,qc_mask_1064 
            ,display_defaults
            ,figs)
   
    if display_defaults.enabled('1064_backscat_image_klett') and hasattr(rs,'rs_inv')\
            and hasattr(rs.rs_inv,'beta_a_1064_backscat_klett'):     
        # particulate backscatter cross section   
        title_str = instrument + '  1064 backscatter-klett, LR=' \
                  +str(processing_defaults.get_value('klett','lidar_ratio_1064')) \
                  + ' ref_alt='\
                  +str(processing_defaults.get_value('klett','ref_altitude')) +'km'    
        gt.rti_fig('1064_backscat_image_klett'
            ,instrument    
            ,rs.rs_inv.beta_a_1064_backscat_klett
            ,rs.rs_inv.times
            ,invalts
            ,title_str
            ,'1/(m sr)'
            ,None 
            ,display_defaults
            ,figs)
        
    if display_defaults.enabled('532_backscat_image_klett') and hasattr(rs,'rs_inv')\
            and hasattr(rs.rs_inv,'beta_a_backscat_klett'):     
        # particulate backscatter cross section   
        title_str = instrument + '  532 backscatter-klett, LR=' \
                  +str(processing_defaults.get_value('klett','lidar_ratio_532')) \
                  + ' ref_alt='\
                  +str(processing_defaults.get_value('klett','ref_altitude')) +'km'    
        gt.rti_fig('532_backscat_image_klett'
            ,instrument    
            ,rs.rs_inv.beta_a_backscat_klett
            ,rs.rs_inv.times
            ,invalts
            ,title_str
            ,'1/(m sr)'
            ,None 
            ,display_defaults
            ,figs)
       
    if display_defaults.enabled('backscat_image') and hasattr(rs,'rs_inv'):
        print 'rendering backscat_image'
        # particulate backscatter cross section   
        title_str = instrument + '  532 backscatter '
        gt.rti_fig('backscat_image'
            ,instrument    
            ,rs.rs_inv.beta_a_backscat if hasattr(rs.rs_inv,'beta_a_backscat')\
                   else (rs.rs_inv.beta_a_backscat_par + rs.rs_inv.beta_a_backscat_perp)
            ,rs.rs_inv.times
            ,invalts 
            ,title_str
            ,'1/(m sr)'
            ,qc_mask 
            ,display_defaults
            ,figs)
            
    if display_defaults.enabled('backscat_SN_image') and hasattr(rs,'rs_inv') and hasattr(rs.rs_inv,'std_beta_a_backscat'):
    
        # particulate backscatter cross section signal-to-noise ratio
        title_str = instrument + '  particulate backscatter S/N '

        gt.rti_fig('backscat_SN_image'
            ,instrument    
            ,rs.rs_inv.SN_beta_a_backscat
            ,rs.rs_inv.times
            ,invalts
            ,title_str
            ,'Signal / Noise'
            ,None
            ,display_defaults
            ,figs)


    if display_defaults.enabled('i2_backscat_image') and hasattr(rs,'rs_inv') and hasattr(rs.rs_inv,'Nm_i2a'):
       
      # particulate backscatter cross section   
        title_str = instrument + '  i2 backscatter cross section '

        gt.rti_fig('i2_backscat_image'
            ,instrument    
            ,(rs.rs_inv.beta_a_backscat_par+rs.rs_inv.beta_a_backscat_perp)
            ,rs.rs_inv.times
            ,invalts
            ,title_str
            ,'1/(m sr)'
            ,qc_mask 
            ,display_defaults
            ,figs)
            
    if display_defaults.enabled('i2a_backscat_image') and hasattr(rs,'rs_inv') and hasattr(rs.rs_inv,'Nm_i2a'):
      # particulate backscatter cross section   
        title_str = instrument + '  i2a backscatter cross section '

        gt.rti_fig('i2a_backscat_image'
            ,instrument    
            ,(rs.rs_inv.beta_a_backscat_par_i2a+rs.rs_inv.beta_a_backscat_perp_i2a)
            ,rs.rs_inv.times
            ,invalts
            ,title_str
            ,'1/(m sr)'
            ,qc_mask 
            ,display_defaults
            ,figs)
           
    

    if display_defaults.enabled('gray_masked_image')\
              and (hasattr(rs,'rs_inv') or hasattr(rs,'rs_particle')):
        var_name = display_defaults.get_value('gray_masked_image','image_variable')
        mask_name = display_defaults.get_value('gray_masked_image','mask_variable')
        var =[]
        if hasattr(rs,'rs_inv') \
               and hasattr(rs.rs_inv,var_name) and hasattr(rs.rs_inv,mask_name):
            var = getattr(rs.rs_inv,var_name)
            mask_var =getattr(rs.rs_inv,mask_name)
        elif hasattr(rs,'rs_particle') and hasattr(rs.rs_particle,var_name) \
                 and hasattr(rs.rs_particle,mask_name):
            var = getattr(rs.rs_particle,var_name)
            mask_var =getattr(rs.rs_particle,mask_name)
        if len(var):    
            lo_mask_level =  display_defaults.get_value('gray_masked_image','lo_mask_level')
            hi_mask_level =  display_defaults.get_value('gray_masked_image','hi_mask_level')
            title_str = var_name + ' with '+ str(lo_mask_level)+ ' < '+ mask_name + ' < '+str(hi_mask_level)        
            mask = qc_mask.copy()
            mask[mask_var < lo_mask_level]= 0
            mask[mask_var > hi_mask_level] = 0
            gt.rti_fig('gray_masked_image'
                ,instrument    
                ,var
                ,rs.rs_inv.times
                ,invalts
                ,title_str
                ,'1/(m sr)'
                ,mask 
                ,display_defaults
                ,figs)      
        
    if display_defaults.enabled('second_backscat_image') and hasattr(rs,'rs_inv'):

      # alternate altitude range particulate backscatter cross section image
        temp = rs.rs_inv.beta_a_backscat_par \
                              +rs.rs_inv.beta_a_backscat_perp
    
        title_str = instrument + '  backscatter cross section 2 '
        max_index = np.int(1000.
                           * float(display_defaults.get_value('second_backscat_image',"max alt(km)")) / 
                           (shared_dz))
        

        #if installed on the ground, set below ground beta_a_backscat =1e-6 
        if not ('installation' in rs_constants) or rs_constants['installation'] == 'ground' :
           ground_index = np.int(rs_constants['lidar_altitude'] \
                         /(shared_dz))
           temp[:,:ground_index+1]=1e-8
       
        

        gt.rti_fig('second_backscat_image'
            ,instrument    
            ,temp[:,:max_index]
            ,rs.rs_inv.times
            ,invalts[ :max_index]
            ,title_str
            ,'1/(m sr)'
            ,None   #qc_mask[:, :max_index] if qc_mask is not None else None
            ,display_defaults
            ,figs
            )
  
    if display_defaults.enabled('circular_depol_image') \
                and hasattr(rs,'rs_inv') and hasattr(rs.rs_inv,'circular_depol'):

        title_str = instrument + '  circular depolarization ' 
        depol = 100 * rs.rs_inv.circular_depol
        units_str = '%'
   
        gt.rti_fig('circular_depol_image'
            ,instrument    
            ,depol
            ,rs.rs_inv.times
            ,invalts
            ,title_str
            ,units_str
            ,qc_mask
            ,display_defaults
            ,figs    
            )
       
    if display_defaults.enabled('linear_depol_image') \
            and hasattr(rs,'rs_inv') and hasattr(rs.rs_inv,'linear_depol'): 
        title_str = instrument + '  linear depolarization  '
        if not hasattr(rs.rs_inv,'circular_depol'):
            depol=100*rs.rs_inv.linear_depol
        else:    
            depol = 100 * rs.rs_inv.circular_depol / (2
                + rs.rs_inv.circular_depol)    
        units_str= '%'               
    
        gt.rti_fig('linear_depol_image'
               ,instrument 
               ,depol
               ,rs.rs_inv.times
               ,invalts
               ,title_str
               ,units_str
               ,qc_mask
               ,display_defaults
               ,figs 
                )
   
    if display_defaults.enabled('second_linear_depol_image') \
            and hasattr(rs,'rs_inv') and hasattr(rs.rs_inv,'circular_depol'):
        title_str = instrument + '  linear depolarization2  '
        depol=100*rs.rs_inv.linear_depol
        units_str= '%'
        max_index = np.int(1000.
                           * float(display_defaults.get_value('second_linear_depol_image',
                           "max alt(km)")) / (shared_dz))
        
        gt.rti_fig('second_linear_depol_image'
            ,instrument    
            ,depol[:, :max_index]
            ,rs.rs_inv.times
            ,invalts[:max_index] 
            ,title_str
            ,units_str
            ,None  #qc_mask[:, :max_index] if qc_mask is not None else None
            ,display_defaults
            ,figs)
    
    if display_defaults.enabled('inverted_mol_image') and hasattr(rs,'rs_inv') and hasattr(rs.rs_inv,'Nm'):
        title_str = instrument + '  inverted molecular, Nm '
        seeded_shots_array = rs.rs_inv.seeded_shots[:, np.newaxis]
        if rs.rs_inv.Nm.shape[0]-rs.rs_inv.seeded_shots.shape[0]!=0:
          seeded_shots_array=copy.deepcopy(seeded_shots_array)
          seeded_shots_array.extend(rs.rs_inv.Nm.shape[0]-rs.rs_inv.seeded_shots.shape[0])
        gt.rti_fig('inverted_mol_image'
            ,instrument    
            ,rs.rs_inv.Nm/seeded_shots_array
            ,rs.rs_inv.times
            ,invalts
            ,title_str
            ,'1/shot/bin'
            ,None
            ,display_defaults
            ,figs)
        
    if display_defaults.enabled('mol_signal_to_noise')  and hasattr(rs,'rs_inv') and hasattr(rs.rs_inv,'SN_mol'):
        title_str = instrument + '  molecular S/N '
        gt.rti_fig('mol_signal_to_noise'
            ,instrument    
            ,rs.rs_inv.SN_mol
            ,rs.rs_inv.times
            ,rs.rs_inv.msl_altitudes
            ,title_str
            ,' '
            ,None
            ,display_defaults
            ,figs)  
        
    if display_defaults.enabled('inverted_aerosol_image') and hasattr(rs,'rs_inv') and hasattr(rs.rs_inv,'Na'):
        title_str = instrument + '  inverted aerosol, Na '
        gt.rti_fig('inverted_aerosol_image'
            ,instrument    
            ,rs.rs_inv.Na
            ,rs.rs_inv.times
            ,invalts
            ,title_str
            ,' '
            ,None
            ,display_defaults
            ,figs)

    if display_defaults.enabled('inverted_cpol_image') and hasattr(rs,'rs_inv') and hasattr(rs.rs_inv,'Ncp'):
        title_str = instrument + '  inverted cpol, Ncp '
        #f=figs.figure(title_str)        
        gt.rti_fig('inverted_cpol_image'
            ,instrument    
            ,rs.rs_inv.Ncp
            ,rs.rs_inv.times
            ,invalts
            ,title_str
            ,' '
            ,None
            ,display_defaults
            ,figs)
   
    if display_defaults.enabled('raw_color_ratio_image') and hasattr(rs,'rs_inv')\
           and hasattr(rs.rs_inv,'raw_color_ratio'):
        title_str = instrument + '  1064/532 combined count ratio '
        gt.rti_fig('raw_color_ratio_image'
            ,instrument    
            ,rs.rs_inv.raw_color_ratio
            ,rs.rs_inv.times
            ,invalts
            ,title_str
            ,'1064/ 532 counts'
            ,qc_mask_1064
            ,display_defaults
            ,figs)
        
    if display_defaults.enabled('color_ratio_image') and hasattr(rs,'rs_inv')\
           and hasattr(rs.rs_inv,'color_ratio'):
        title_str = instrument + '  1064/532 aerosol backscat ratio, A= '\
                    +str(processing_defaults.get_value('color_ratio','angstrom_coef'))

        gt.rti_fig('color_ratio_image'
            ,instrument    
            ,rs.rs_inv.color_ratio
            ,rs.rs_inv.times
            ,invalts
            ,title_str
            ,'1064/532 backscatter'
            ,qc_mask_1064
            ,display_defaults
            ,figs)
       
    if display_defaults.enabled('second_color_ratio_image') and hasattr(rs,'rs_inv')\
           and hasattr(rs.rs_inv,'color_ratio'):
        title_str = instrument + '  1064/532 aerosol backscat ratio, A= '\
                    +str(processing_defaults.get_value('color_ratio','angstrom_coef'))

        gt.rti_fig('second_color_ratio_image'
            ,instrument
            ,rs.rs_inv.color_ratio
            ,rs.rs_inv.times
            ,invalts
            ,title_str
            ,'1064/532 backscatter'
            ,None #qc_mask[:, :max_index] if qc_mask is not None else None
            ,display_defaults
            ,figs)
    if display_defaults.enabled('aerosol_optical_depth_image') and hasattr(rs,'rs_inv') \
           and hasattr(rs.rs_inv,'optical_depth_aerosol'):
    
        title_str = instrument + '  532nm aerosol optical depth'
        print 'rendering '+title_str    
        gt.rti_fig('aerosol_optical_depth_image'
            ,instrument    
            ,rs.rs_inv.optical_depth_aerosol
            ,rs.rs_inv.times
            ,invalts
            ,title_str
            ,None
            ,qc_mask  if qc_mask is not None else None
            ,display_defaults
            ,figs)


    
    if display_defaults.enabled('second_extinction_image') and hasattr(rs,'rs_inv') and hasattr(rs.rs_inv,'extinction'):
        if 1: # processing_defaults.get_value('extinction_processing'
              #    ,'filter_type') == 'savitzky_golay':
            ddz =\
               processing_defaults.get_value('extinction_processing'
                    ,'alt_window_length')
            dt = processing_defaults.get_value('extinction_processing'
                    ,'time_window_length')                           
            title_str = instrument + '  extinction--Savitsky_Golay, dz= ' \
               + str(ddz) + ', dt= ' + str(dt) 
         
        
        gt.rti_fig('second_extinction_image'
            ,instrument    
            ,rs.rs_inv.extinction[:, :max_index]
            ,rs.rs_inv.times
            ,invalts[:max_index]
            ,title_str
            ,'1/m'
            ,None  #qc_mask[:, :max_index]  if qc_mask is not None else None
            ,display_defaults
            ,figs)
   
    if display_defaults.enabled('wfov_ext_corr_image') and hasattr(rs,'rs_inv') \
           and hasattr(rs.rs_inv,'wfov_extinction_corr'):
        print 'rendering wfov_ext_corr_image'
        if 1: 
            dt = processing_defaults.get_value('wfov_corr'
                  ,'window_durration')                             
               
            title_str = instrument + '  wfov_ext_corr '
        
        gt.rti_fig('wfov_ext_corr_image'
            ,instrument    
            ,rs.rs_inv.wfov_extinction_corr *1e6
            ,rs.rs_inv.times
            ,invalts
            ,title_str
            ,'1/m *1e6'
            ,qc_mask
            ,display_defaults
            ,figs)
        
    if display_defaults.enabled('extinction_image') and hasattr(rs,'rs_inv') and hasattr(rs.rs_inv,'extinction'):
       
        if 1: #processing_defaults.get_value('extinction_processing'
              #    ,'filter_type') == 'savitzky_golay':
        
            ddz =2*shared_dz * int(np.ceil(
                  processing_defaults.get_value('extinction_processing'
                  ,'alt_window_length')/(2*shared_dz)))
            dt = processing_defaults.get_value('extinction_processing'
                  ,'time_window_length')                             
               
            title_str = instrument + '  extinction--Savitsky_Golay, dz=' \
               + str(ddz)+ ',dt='+str(dt)
        
        gt.rti_fig('extinction_image'
            ,instrument    
            ,rs.rs_inv.extinction
            ,rs.rs_inv.times
            ,invalts
            ,title_str
            ,'1/m'
            ,qc_mask
            ,display_defaults
            ,figs)
        
    if display_defaults.enabled('phase_function_image') and hasattr(rs,'rs_inv') and hasattr(rs.rs_inv,'p180'):
        if 1: #processing_defaults.get_value('extinction_processing'
              #   ,'filter_type') == 'savitzky_golay':
            ddz =2*shared_dz * int(np.ceil( \
               processing_defaults.get_value('extinction_processing'
               ,'alt_window_length')/(2*shared_dz)))
            dt = processing_defaults.get_value('extinction_processing'
               ,'time_window_length')
            title_str = instrument + '  phase function--Savitsky_Golay, dz=' \
               + str(ddz) + ',dt=' + str(dt)

        gt.rti_fig('phase_function_image'
         ,instrument
         ,rs.rs_inv.p180
         ,rs.rs_inv.times
         ,invalts
         ,title_str
         ,'1/str'
         ,qc_mask
         ,display_defaults
         ,figs)
        
    if display_defaults.enabled('lidarRatio_image') and hasattr(rs,'rs_inv') and hasattr(rs.rs_inv,'p180'):
        if processing_defaults.get_value('extinction_processing'
                 ,'filter_type') == 'savitzky_golay' :
            ddz =2*shared_dz * int(np.ceil( \
               processing_defaults.get_value('extinction_processing'
               ,'alt_window_length')/(2*shared_dz)))
            dt = processing_defaults.get_value('extinction_processing'
               ,'time_window_length')                                  
            title_str = instrument + '  lidar ratio--Savitsky_Golay, dz= ' \
               + str(ddz)+', dt= '+str(dt)
        else:
            processing_defaults.get_value('extinction_processing','bin_delta'
                     )  
            title_str = instrument + ' lidar ratio--block_ave , dz= ' \
                 + str(ddz)
       
        lidarRatio = 1.0/rs.rs_inv.p180        
        print 'rendering lidar ratio image'
        gt.rti_fig('lidarRatio_image'
             ,instrument
             ,lidarRatio
             ,rs.rs_inv.times
             ,invalts
             ,title_str
             ,'str'
             ,qc_mask
             ,display_defaults
             ,figs)
           
  
    if display_defaults.enabled('raqms_total_extinction_image') and sounding is not None and hasattr(sounding,'soundings_used') and hasattr(sounding.soundings_used,'ext_total'):
           title_str = 'raqms extinction cross section'
           indices = np.arange(len(sounding.soundings_used.altitudes))
           top_index = np.max(indices[sounding.soundings_used.altitudes <= max_alt])
       
           gt.rti_fig('raqms_total_extinction_image'
               ,instrument    
               ,sounding.soundings_used.ext_total[:,:top_index]*1e-3
               ,sounding.soundings_used.times
               ,sounding.soundings_used.altitudes[:top_index]
               ,title_str
               ,'1/m'
               ,None
               ,display_defaults
               ,figs)


    if False:     
      start_time_str = timerange[0].strftime('%d-%b-%y %H:%M')
      if timerange[0].year == timerange[-1].year\
               and timerange[0].month == timerange[-1].month \
               and timerange[0].day == timerange[-1].day:
          end_time_str = timerange[-1].strftime('%H:%M')
      else:
          end_time_str = timerange[-1].strftime('%d-%b-%y %H:%M')

   # vertical profiles averaged over requested time interval

    if display_defaults.enabled('backscat_profile') and haveProfiles and hasattr(rs.profiles,'inv'):
        
        ref_index = rs.profiles.inv.mol_norm_index
        #beta_a_backscat = rs.profiles.inv.beta_a_backscat_par \
        #                          +rs.profiles.inv.beta_a_backscat_perp
        tmpbrb=rs.profiles.inv.beta_r_backscat.copy()
        tmpbrb[np.isnan(tmpbrb)]=0.0
        mol_od = (8*np.pi/3.0)*np.cumsum(tmpbrb)*(shared_dz)
       
        if hasattr(rs.profiles.inv,'Nm_i2a'):
            lines =  [rs.profiles.inv.beta_a_backscat[0, :]
                    ,(rs.profiles.inv.Nm_i2[0, :] + rs.profiles.inv.Nm_i2a[0,:])
                    * rs.profiles.inv.beta_r_backscat[ref_index]
                    *np.exp(2*mol_od)
                    / (rs.profiles.inv.Nm_i2[0,ref_index] + rs.profiles.inv.Nm_i2a[0,ref_index])
                    ,rs.profiles.inv.beta_r_backscat]
            colors = ['r','b','k']
            legend = ['beta_a','(mol+mol_i2a)*exp(-2*AOD)','beta_R']
        else:
            lines =  [rs.profiles.inv.beta_a_backscat[0, :]
                    ,rs.profiles.inv.Nm_i2[0, :] * rs.profiles.inv.beta_r_backscat[ref_index]
                    *np.exp(2*mol_od) /(rs.profiles.inv.Nm_i2[0,ref_index] * np.exp(2*mol_od[ref_index]))
                    ,rs.profiles.inv.beta_r_backscat]
            colors = ['r','b','k','c','c','c']
            legend = ['beta_a','mol*exp(-2*AOD)','beta_R']
            """
            lines.append(rs.profiles.inv.beta_a_backscat[0,:]*20)
            lines.append(rs.profiles.inv.beta_a_backscat[0,:]*100)
            lines.append(rs.profiles.inv.beta_a_backscat[0,:]/0.035)
            """
        gt.plot_vs_altitude('backscat_profile'                     
                 ,instrument                   
                 ,timerange                        
                 ,invalts# rs.rs_mean.msl_altitudes                    
                 ,lines
                 ,colors               
                 ,[]
                 ,legend                 
                 ,'upper right'          
                 ,'532nm ackscatter cross section'                  
                 ,'1/(m sr)'                         
                 ,'532 backscatter cross section'            
                 ,auto_loop                
                 ,display_defaults                    
                 ,figs)
    if display_defaults.enabled('532_backscat_profile_klett')\
            and haveProfiles and hasattr(rs.profiles,'inv')\
            and hasattr(rs.profiles.inv,'beta_a_532_backscat_klett')\
            and hasattr(rs.profiles.inv,'beta_a_backscat'):
        print 'Rendering 532 backscat profile klett'

        
        gt.plot_vs_altitude('532_backscat_profile_klett'                     
                 ,instrument                   
                 ,timerange                        
                 ,invalts                
                 ,[rs.profiles.inv.beta_a_532_backscat_klett[0,:],rs.profiles.inv.beta_a_backscat[0,:]]           
                 ,['c','r']
                 ,[2,2]           
                 ,['klett','std']                 
                 ,'upper right'          
                 ,'532 nm Backscatter cross section-klett'                  
                 ,'1/(m sr)'                         
                 ,'532 nm Backscatter-klett'            
                 ,auto_loop                
                 ,display_defaults                    
                 ,figs)      
    if display_defaults.enabled('1064_backscat_profile_klett')\
            and haveProfiles and hasattr(rs.profiles,'inv')\
            and hasattr(rs.profiles.inv,'beta_a_1064_backscat_klett')\
            and hasattr(rs.profiles.inv,'beta_a_1064_backscat'):
        print 'Rendering 1064 backscat profile klett'
        ref_bin = np.int(processing_defaults.get_value('klett','ref_altitude')*1000.0\
                         /(rs.profiles.inv.msl_altitudes[2]-rs.profiles.inv.msl_altitudes[1]))
        if ref_bin > len(invalts):
            ref_bin = len(invalts)-1
        backscat_no_atten = rs.profiles.combined_1064_counts[0,:] \
           * rs.profiles.inv.beta_r_backscat[ref_bin] \
           /(16.0 * rs.profiles.combined_1064_counts[0,ref_bin])\
           - rs.profiles.inv.beta_r_backscat/16.0
        
        gt.plot_vs_altitude('1064_backscat_profile_klett'                     
                 ,instrument                   
                 ,timerange                        
                 ,invalts                
                 ,[rs.profiles.inv.beta_a_1064_backscat_klett[0,:],rs.profiles.inv.beta_a_1064_backscat[0,:] \
                   ,backscat_no_atten]           
                 ,['r','k','c']
                 ,[2,2,1]           
                 ,['klett','std','no_att']                 
                 ,'upper right'          
                 ,'1064 nm Backscatter cross section-klett'                  
                 ,'1/(m sr)'                         
                 ,'1064 nm Backscatter-klett'            
                 ,auto_loop                
                 ,display_defaults                    
                 ,figs)
        
    if display_defaults.enabled('1064_backscat_profile') and haveProfiles and hasattr(rs.profiles,'inv') and hasattr(rs.profiles.inv,'beta_a_1064_backscat'):
        #ref_index = rs.rs_inv.mol_norm_index
        ref_index = rs.profiles.inv.mol_norm_index
        tmpbrb=rs.profiles.inv.beta_r_backscat.copy()
        tmpbrb[np.isnan(tmpbrb)]=0.0
        mol_od = (8*np.pi/3.0)*np.cumsum(tmpbrb)*(shared_dz)


       
        if hasattr(rs.profiles.inv,'Nm_i2a'):
            mol_profile = (rs.profiles.inv.Nm_i2[0, :] + rs.profiles.inv.Nm_i2a[0,:])\
                    * rs.profiles.inv.beta_r_backscat[ref_index]\
                    *np.exp(2*mol_od)\
                    / (rs.profiles.inv.Nm_i2[0,ref_index] + rs.profiles.inv.Nm_i2a[0,ref_index])
            lines =  [rs.profiles.inv.beta_a_backscat[0, :]
                    ,mol_profile  
                    ,rs.profiles.inv.beta_r_backscat
                    ,rs.profiles.inv.beta_r_backscat/16.0  
                    ,rs.profiles.inv.beta_a_1064_backscat[0,:]]
            colors = ['r','b','k','k','c']
            legend = ['beta_a','(mol+mol_i2a)*exp(-2*AOD)','beta_R','beta_R_1064','beta_a_1064']
        else:
            lines =  [rs.profiles.inv.beta_a_backscat[0, :]
                    ,rs.profiles.inv.Nm_i2[0, :] * rs.profiles.inv.beta_r_backscat[ref_index]
                    *np.exp(2*mol_od) /(rs.profiles.inv.Nm_i2[0,ref_index] * np.exp(2*mol_od[ref_index]))
                    ,rs.profiles.inv.beta_r_backscat
                    ,rs.profiles.inv.beta_r_backscat/16.0
                    ,rs.profiles.inv.beta_a_1064_backscat[0,:]]
                   
            colors = ['r','b','k','m','c']
            legend = ['beta_a','mol*exp(-2*AOD)','beta_R','beta_R_1064','beta_a_1064']
            
        gt.plot_vs_altitude('1064_backscat_profile'                     
                 ,instrument                   
                 ,timerange                        
                 ,invalts# rs.rs_mean.msl_altitudes                    
                 ,lines
                 ,colors               
                 ,None
                 ,legend                 
                 ,'upper right'          
                 ,'1064 nm Backscatter cross section'                  
                 ,'1/(m sr)'                         
                 ,'1064 nm Backscatter'            
                 ,auto_loop                
                 ,display_defaults                    
                 ,figs)
        
    if hasattr(rs,'profiles') and hasattr(rs.profiles.inv,'beta_a_backscat_par_i2a') \
             and display_defaults.enabled('i2_and_i2a_backscat_profiles'):
      
        ref_index = rs.profiles.inv.mol_norm_index
      
        mol_od = (8*np.pi/3.0)*np.cumsum(rs.profiles.inv.beta_r_backscat)*shared_dz
       
        lines =  [rs.profiles.inv.beta_a_backscat_par[0,:] \
                            + rs.profiles.inv.beta_a_backscat_perp[0,:]
                  ,rs.profiles.inv.beta_a_backscat_par_i2a[0,:] \
                            +rs.profiles.inv.beta_a_backscat_perp_i2a[0,:]
                  ,rs.profiles.inv.Nm_i2[0, :] * rs.profiles.inv.beta_r_backscat[ref_index]
                         *np.exp(2*mol_od) / rs.profiles.inv.Nm_i2[0,ref_index]
                  ,rs.profiles.inv.Nm_i2a[0, :] * rs.profiles.inv.beta_r_backscat[ref_index]
                         *np.exp(2*mol_od) / rs.profiles.inv.Nm_i2a[0,ref_index]
                  ,rs.profiles.inv.beta_r_backscat]
        colors = ['r','m','b','g','k']
        legend = ['beta_i2','beta_i2a','mol_i2*exp(-2*A)D)','mol_i2a*exp(-2*A)D)','beta_R']
        
            
        gt.plot_vs_altitude('i2_and_i2a_backscat_profiles'                     
                 ,instrument                   
                 ,timerange                       
                 ,invalts# rs.rs_mean.msl_altitudes                    
                 ,lines
                 ,colors               
                 ,[2,2,2,2,2]
                 ,legend                 
                 ,'upper right'          
                 ,'I2 and I2A Backscatter cross section'                  
                 ,'1/(m sr)'                         
                 ,'I2 and I2A backscatter cross section'            
                 ,auto_loop                
                 ,display_defaults                    
                 ,figs) 
    
    if display_defaults.enabled('extinction_and_p180_profile') and haveProfiles and hasattr(rs.profiles,'inv') and hasattr(rs.profiles.inv,'extinction_aerosol'):
        
        beta_a_backscat = rs.profiles.inv.beta_a_backscat_par \
                                  +rs.profiles.inv.beta_a_backscat_perp         
        temp_ae = rs.profiles.inv.extinction_aerosol.copy()
        temp_p180 = rs.profiles.inv.p180.copy()
        temp_p180[np.isinf(temp_p180)] = 1.0e-10
        # avoid overflow errors
        temp_ae[np.isnan(temp_ae)] = 1.0e-11
        temp_ae[np.less(temp_ae,1.0e-10)]= 1.0e-10
        beta_a_backscat[np.isnan(beta_a_backscat)] =0.0
        beta_a_backscat[beta_a_backscat<=0.0] =1e-15
        temp_ae[np.isinf(temp_ae)] = 1.0e-10
        
        gt.plot_vs_altitude('extinction_and_p180_profile'                     
                 ,instrument                   
                 ,timerange                
                 ,invalts              
                 ,[rs.profiles.inv.beta_r_backscat*8*np.pi/3.0
                   ,temp_p180[0,:]
                   ,temp_ae[0,:]] #,beta_a_backscat[0,:] / temp_ae[0,:],temp_ae[0,:]]
                 ,['b','c','r']               
                 ,None
                 ,['Ray_ext','P(180)/4pi','Aerosol_ext']                 
                 ,'upper right'          
                 ,'extinction cross section (1/m), p180/4pi (1/sr)'                  
                 ,None                         
                 ,'Extinction, P(180)/4pi'            
                 ,auto_loop                
                 ,display_defaults                    
                 ,figs)
   
    if display_defaults.enabled('extinction_profile') and hasattr(rs,'profiles') and hasattr(rs.profiles,'inv')\
             and hasattr(rs.profiles.inv,'extinction_aerosol'):

        lines =[rs.profiles.inv.beta_r_backscat*8*np.pi/3.0]
        legend = ['Rayleigh','aerosol','total']

        if hasattr(rs.profiles,'wfov_geo_adjust'):
            xlabel = 'wfov corrected extinction'
        else:
            xlabel = 'extinction cross section'

        lines.append(rs.profiles.inv.extinction_aerosol[0,:])
        lines.append(rs.profiles.inv.extinction_aerosol[0,:]
                   +rs.profiles.inv.beta_r_backscat*8*np.pi/3.0)
            
        gt.plot_vs_altitude('extinction_profile'                     
                 ,instrument                   
                 ,timerange                
                 ,invalts              
                 ,lines
                 ,['b','r','c']               
                 ,None
                 ,legend                
                 ,'upper right'          
                 , xlabel                  
                 ,'1/m'                         
                 ,xlabel            
                 ,auto_loop                
                 ,display_defaults                    
                 ,figs)    

    if display_defaults.enabled('integrated_backscatter_profile') \
                 and haveProfiles and hasattr(rs.profiles,'inv'):
        
        beta_a_backscat_profile = rs.profiles.inv.beta_a_backscat_par \
                              +rs.profiles.inv.beta_a_backscat_perp
        ddz=shared_dz

        if hasattr(rs.profiles,'telescope_pointing')\
               and rs.profiles.telescope_pointing[0]<.1:
            ddz = -ddz
        
        beta_a_backscat_profile[np.isnan(beta_a_backscat_profile)] = 0.0
        int_backscat = ddz * np.cumsum(beta_a_backscat_profile)
        int_backscat = int_backscat \
            - int_backscat[rs.profiles.inv.mol_norm_index]
        p180_4pi = display_defaults.get_value('integrated_backscatter_profile','p180/4pi')
       
        od_at_norm_alt = rs.profiles.inv.optical_depth_aerosol[0,rs.profiles.inv.mol_norm_index]
        aerosol_od = rs.profiles.inv.optical_depth_aerosol[0,:]
        gt.plot_vs_altitude('integrated_backscatter_profile'       
                 ,instrument                   
                 ,timerange                    
                 ,invalts                
                 ,[int_backscat / p180_4pi
                     ,aerosol_od
                     ,int_backscat * 70
                     ,int_backscat * 60
                     ,int_backscat * 50
                     ,int_backscat * 40
                     ,int_backscat * 30
                     ,int_backscat * 20]
                 ,['b','r','b','c','m','k','g','m']               
                 ,[3,3,1,1,1,1,1,1]
                 ,[('intBS*%i'%(1.0/p180_4pi))
                     ,'aerosol_od'
                     ,'intBS*70'
                     ,'intBS*60'
                     ,'intBS*50'
                     ,'intBS*40'
                     ,'intBS*30'
                     ,'intBS*20']                 
                 ,'lower right'          
                 ,'Optical depth '                  
                 ,None                         
                 ,'Integrated backscatter'            
                 ,auto_loop                
                 ,display_defaults                    
                 ,figs)
  
    if display_defaults.enabled('i2a_temperature') and haveProfiles \
           and hasattr(rs.profiles,'i2a_temperatures'):
       
        raob_str=  sounding.station_id +' '\
                + sounding.times.strftime("%HZ")
        max_index = len(invalts)
        temp_offset =-22
        legend3='i2a+' + str(temp_offset)
        gt.plot_vs_altitude('i2a_temperatures'                     
                 ,instrument                   
                 ,timerange                    
                 ,invalts                 
                 ,[rs.profiles.i2a_temperatures,sounding.temps[:max_index]]
                 ,['r','k','c']               
                 ,[2,1,1]
                 ,['i2a_temp',raob_str,legend3]                 
                 ,'upper right'          
                 ,'i2a temperatures'                  
                 ,'Deg K'                         
                 ,'i2a temperature'          
                 ,auto_loop                
                 ,display_defaults                    
                 ,figs)
     
    if display_defaults.enabled('depol_profile') and haveProfiles:
    
        beta_a_backscat_profile = rs.profiles.inv.beta_a_backscat_par \
                              +rs.profiles.inv.beta_a_backscat_perp
        if not hasattr(rs.profiles.inv,'linear_depol'):
             linear_depol = (100*rs.profiles.inv.circular_depol[0, :] / (2
                    + rs.profiles.inv.circular_depol[0, :])).transpose()
        else:
            linear_depol = 100*rs.profiles.inv.linear_depol[0,:]
        mask = np.ones_like(linear_depol)*np.NaN
        np.seterr(invalid='ignore')
        mask[beta_a_backscat_profile[0,:] > 1e-7] = 1
        np.seterr(invalid='warn')
        gt.plot_vs_altitude('depol_profile'                     
                 ,instrument                   
                 ,timerange                       
                 ,invalts                 
                 ,[linear_depol,linear_depol*mask]
                 ,['c','g']               
                 ,[1,3]
                 ,['b_a<1e-7','b_a>1e-7']                 
                 ,'lower right'          
                 ,'Linear depolarization '                  
                 ,'%'                         
                 ,'Linear depolarization'            
                 ,auto_loop                
                 ,display_defaults                    
                 ,figs)
  
    #od plot provided only if telescope is always pointed up or always down
    
    if display_defaults.enabled('od_profile') and haveProfiles \
             and (not ('installation' in rs_constants) \
                  or rs_constants['installation'] == 'ground' \
                  or rs_constants['installation'] == 'shipborne'\
                  or (hasattr(rs.profiles,'telescope_pointing') and
                    (rs.profiles.telescope_pointing[0] > 0.9 \
                  or rs.profiles.telescope_pointing[0] < 0.1))): 
        
        #Rayliegh backscatter phase function time altitude step size
        #pdr=(8*np.pi/3)*(rs.rs_mean.msl_altitudes[3]-rs.rs_mean.msl_altitudes[2])
        #pdr = pdr/np.cos(np.pi*rs_constants['telescope_zenith_angle']/180.0)
        #if no installation in calvals assume ground based instrument
        #or if airborne based the telescope always pointed up
        ref_alt = invalts[rs.profiles.inv.mol_norm_index]/1000.0
        xlabel = 'Optical depth (ref alt = '+ str(ref_alt)+' km)'
       
        if not('installation' in rs_constants) \
                or (rs_constants['installation'] == 'ground')\
                or (rs_constants['installation'] == 'shipborne')\
                or rs.profiles.telescope_pointing[0] > 0.9:
            print 'ground based with zenith pointing telescope'
            #Rayliegh backscatter phase function time altitude step size
            pdr=(8*np.pi/3)*(shared_dz)
            pdr = pdr/np.cos(np.pi*rs_constants['telescope_roll_angle_offset']/180.0)
           
            od_at_norm_alt= \
               pdr*np.sum(rs.profiles.inv.beta_r_backscat[0:rs.profiles.inv.mol_norm_index])
            
            if hasattr(rs.profiles.inv,'wfov_corr_optical_depth'):

                od_total=rs.profiles.inv.wfov_corr_optical_depth[0,:] 
                od_aerosol = rs.profiles.inv.wfov_corr_optical_depth_aerosol[0,:]
                xlabel = 'WFOV corrected OD (ref alt = '+ str(ref_alt)+' km)'
            else:
                od_total=rs.profiles.inv.optical_depth[0,:] 
                od_aerosol = rs.profiles.inv.optical_depth_aerosol[0,:]
                title = 'Optical depth'
                              
        #if airborne with telescope always pointed at nadir
        elif (rs_constants['installation']=='airborne')\
            and hasattr(rs.profiles,'telescope_pointing') and rs.profiles.telescope_pointing[0]<0.1:
                print 'airborne with telescope always pointed at nadir'
                mean_alt=rs.profiles.mean_GPS_MSL_Alt
                aircraft_alt_index = len(invalts[invalts<=mean_alt])
                #Rayliegh backscatter phase function time altitude step size
                pdr=-(8*np.pi/3)*shared_dz
                pdr = pdr/np.cos(np.pi*rs_constants['telescope_roll_angle_offset']/180.0)
                mol_optical_depth = np.zeros_like(rs.profiles.inv.beta_r_backscat)
        
                #make reverse list of indices for integration looking downward
                indices = range(len(mol_optical_depth[:aircraft_alt_index])-1,-1,-1)
                mol_optical_depth[:aircraft_alt_index]= pdr*np.cumsum(rs.profiles.inv.beta_r_backscat[indices])
                mol_optical_depth=mol_optical_depth \
                      -mol_optical_depth[rs.profiles.inv.mol_norm_index]
                 
                od_aerosol=rs.profiles.inv.optical_depth[0,:] - mol_optical_depth
                od_total = rs.profiles.inv.optical_depth[0,:]

        #plot if ground based or if nearly all telescope up or down pointing
        if not rs_constants['installation']=='airborne'\
            or  rs.profiles.telescope_pointing[0] <0.1 \
            or  rs.profiles.telescope_pointing[0] >0.9:

            gt.plot_vs_altitude('od_profile'       #display defaults plot name
                 ,instrument                       #instrument name
                 ,timerange                        #python datetimes vector
                 ,invalts                          #altitude vector (meters)
                 ,[od_aerosol,od_total]            #variables
                 ,['r','b']                        #colors, [] default colors
                 ,None                             #widths, [] sets widths = 2
                 ,['aerosol','total']              #legend list, [] = no legend
                 ,'lower right'                    #legend position, [] ok if list []
                 ,xlabel                           #xlabel
                 ,None                             #x units
                 ,'Optical depth'                  #plot title
                 ,auto_loop                        # =1, clear figure before new plot
                 ,display_defaults                    
                 ,figs)
        else:
            print ' '
            print 'mixed up and down pointing no OD plot'
            print ' '
    
               
    if display_defaults.enabled('raw_color_ratio_profile') and haveProfiles \
             and hasattr(rs,'profiles') and hasattr(rs.profiles,'inv') \
             and hasattr(rs.profiles.inv,'raw_color_ratio')\
             and hasattr(rs.profiles.inv,'color_ratio'):
            v_line = np.ones(len(invalts))
            v_line[invalts < 0] =np.NaN
            title = 'IR/532 color ratios corrected for energy and OD'
            traces = []
            legend = []
            angstrom_coef = [0.0, 0.5, 1.0, 1.5, 2.0]
                       
            for i in range(len(angstrom_coef)):
               traces.append(rs.profiles.inv.raw_color_ratio[0,:] \
                          *np.exp(-2*rs.profiles.inv.optical_depth_aerosol[0,:]*(1-0.5**angstrom_coef[i])))
               legend.append('A='+str(angstrom_coef[i]))
            traces.append(v_line)
            traces.append(v_line/16.0)
            legend.append('mol')
            legend.append('cloud')
            gt.plot_vs_altitude('raw_color_ratio_profile'     #display defaults plot name
                 ,instrument                                  #instrument name
                 ,timerange                                   #python datetimes vector
                 ,invalts                                     #altitude vector (meters)
                 ,traces                                      #variables
                 ,['c','r','m','b','k','k','k']                #colors, [] default colors
                 ,None                                        #widths, [] sets widths = 2
                 ,legend                                      #legend list, [] = no legend
                 ,'upper right'                               #legend position, [] ok if list []
                 ,'1064/532 backscatter--mol atten corrected'                             #xlabel
                 ,None                                        #x units
                 ,title                                       #plot title
                 ,auto_loop                                   # =1, clear figure before new plot
                 ,display_defaults                    
                 ,figs)    
                   
    if display_defaults.enabled('color_ratio_profile')\
            and haveProfiles \
            and hasattr(rs,'profiles') and hasattr(rs.profiles,'inv')\
            and hasattr(rs.profiles.inv,'color_ratio'):
       title = 'color ratio, Angst coef = '\
                    +str(processing_defaults.get_value('color_ratio','angstrom_coef'))                           
       gt.plot_vs_altitude('color_ratio_profile'              #display defaults plot name
                 ,instrument                                  #instrument name
                 ,timerange                                   #python datetimes vector
                 ,rs.profiles.inv.msl_altitudes               #altitude vector (meters)
                 ,[rs.profiles.inv.color_ratio[0,:]]          #variables
                 ,['r','k','b']                               #colors, [] default colors
                 ,None                                        #widths, [] sets widths = 2
                 ,None                                        #legend list, [] = no legend
                 ,''                                          #legend position, [] ok if list []
                 ,'1064/532 aerosol backscatter'              #xlabel
                 ,[]                                          #x units
                 ,title                                       #plot title
                 ,auto_loop                                   # =1, clear figure before new plot
                 ,display_defaults                    
                 ,figs)    
      
                 
    if display_defaults.enabled('color_ratio_profile_by_angstrom_coef')\
            and haveProfiles \
            and hasattr(rs,'profiles') and hasattr(rs.profiles,'inv')\
            and hasattr(rs.profiles.inv,'color_ratio'):
       assumed_angstrom_coef = processing_defaults.get_value('color_ratio','angstrom_coef')
       mol_optical_depth_532 = (8.0*np.pi/3.0) * np.cumsum(rs.profiles.inv.beta_r_backscat) \
                               * (rs.profiles.inv.msl_altitudes[2]-rs.profiles.inv.msl_altitudes[1])
       angstrom_coef =[0.0, 0.5, 1.0, 1.5, 2.0]
       #exp1 =np.exp(-2.0*0.5**assumed_angstrom_coef*rs.profiles.inv.optical_depth_aerosol[0,:])
       lines = []
       legend =[]
       scat_ratio_532 = rs.profiles.inv.beta_a_backscat[0,:]/rs.profiles.inv.beta_r_backscat
       for i in range(len(angstrom_coef)):
           lines.append(rs.profiles.inv.raw_color_ratio[0,:] * (1 + 1.0/scat_ratio_532) \
               *np.exp(-2 * rs.profiles.inv.optical_depth_aerosol[0,:] \
               * (1.0 - 0.5**angstrom_coef[i]) - mol_optical_depth_532 * 30.0 / 16.0)\
               - 1.0/(16.0 * scat_ratio_532))
           """
           lines.append(
              exp1
              *np.exp(2*0.5**angstrom_coef[i]*rs.profiles.inv.optical_depth_aerosol[0,:])\
              *(rs.profiles.inv.color_ratio[0,:]-1.0/(16.0 *scat_ratio_532))\
              +1.0/(16.0*scat_ratio_532))
           """            
           legend.append(str(angstrom_coef[i]))
       gt.plot_vs_altitude('color_ratio_profile_by_angstrom_coef'  #display defaults plot name
                 ,instrument                                  #instrument name
                 ,timerange                                   #python datetimes vector
                 ,rs.profiles.inv.msl_altitudes               #altitude vector (meters)
                 ,lines                                       #variables
                 ,None                                        #colors, [] default colors
                 ,None                                        #widths, [] sets widths = 2
                 ,legend                                      #legend list, [] = no legend
                 ,'lower right'                               #legend position, [] ok if list []
                 ,'1064/532 aerosol backscatter'              #xlabel
                 ,None                                        #x units
                 ,'1064/532 vs angstrom coef'              #plot title
                 ,auto_loop                                   # =1, clear figure before new plot
                 ,display_defaults                    
                 ,figs)          

    if display_defaults.enabled('sc_ratio_profile') and haveProfiles and hasattr(rs.profiles,'inv'):       
        sc_ratio_par =  rs.profiles.inv.Na[0, :] / rs.profiles.inv.Nm_i2[0, :]
        sc_ratio_perp = rs.profiles.inv.Ncp[0, :] / rs.profiles.inv.Nm_i2[0, :]
        lines = [sc_ratio_par,sc_ratio_perp]
        legend = ['i2_par','i2_perp']
        if hasattr(rs.profiles.inv,'Na_i2a'):
            lines.append(rs.profiles.inv.Na_i2a[0,:]\
                                 /rs.profiles.inv.Nm_i2a[0,:])
            lines.append(rs.profiles.inv.Ncp[0,:]\
                                 /rs.profiles.inv.Nm_i2a[0,:])
            legend.append('i2a_par')
            legend.append('i2a_perp')
        if (not 'installation' in rs_constants \
               or rs_constants['installation'] == 'ground'\
               or rs_constants['installation'] == 'shipborne') \
               and 'lidar_altitude' in rs_constants:
            sc_ratio_par[invalts < rs_constants['lidar_altitude']+100] = 0 
            sc_ratio_perp[invalts < rs_constants['lidar_altitude']+100] = 0
       

        gt.plot_vs_altitude('sc_ratio_profile'        #display defaults plot name
                 ,instrument                       #instrument name
                 ,timerange                        #python datetimes vector
                 ,invalts                          #altitude vector (meters)
                 ,lines                            #variables
                 ,['r','g','m','b']                #colors, [] default colors
                 ,[2,1,2,1]                               #widths, [] sets widths = 2
                 ,legend                           #legend list, [] = no legend
                 ,'upper right'                    #legend position, [] ok if list []
                 ,'Scattering ratio'               #xlabel
                 ,None                             #x units
                 ,'Scattering ratio'               #plot title
                 ,auto_loop                        # =1, clear figure before new plot
                 ,display_defaults                    
                 ,figs)  



    if display_defaults.enabled('sc_ratio_errors') and haveProfiles and hasattr(rs.profiles,'inv'):
          
        sc_ratio  =  (rs.profiles.inv.Na[0, :] + rs.profiles.inv.Ncp[0,:])/ rs.profiles.inv.Nm[0, :]

        if (not 'installation' in rs_constants \
               or rs_constants['installation'] == 'ground'\
               or rs_constants['installation'] == 'shipborne') \
               and 'lidar_altitude' in rs_constants:
            sc_ratio[rs.profiles.msl_altitudes < rs_constants['lidar_altitude']+100] = 0 
           
        gt.plot_vs_altitude('sc_ratio_errors'       #display defaults plot name
                 ,instrument                        #instrument name
                 ,timerange                         #python datetimes vector
                 ,rs.profiles.msl_altitudes         #altitude vector (meters)
                 ,[sc_ratio-rs.profiles.inv.SR_std[0,:]
                   ,sc_ratio + rs.profiles.inv.SR_std[0,:]
                   ,sc_ratio]                      #variables
                 ,['c','c','r']                    #colors, [] default colors
                 ,None                             #widths, [] sets widths = 2
                 ,['-std','SR_par','+std']         #legend list, [] = no legend
                 ,'upper right'                    #legend position, [] ok if list []
                 ,'Scattering ratio'               #xlabel
                 ,None                             #x units
                 ,'SR count errors'                #plot title
                 ,auto_loop                        # =1, clear figure before new plot
                 ,display_defaults                    
                 ,figs)
 
    
    # plot ratios of corrected profiles for evaluation of differential geometry

    if display_defaults.enabled('dif_geo_profiles') and haveProfiles \
           and hasattr(rs.profiles,'molecular_counts') :

        inx = rs.profiles.inv.mol_norm_index
        lines = [rs.profiles.combined_hi_counts[0, :]
                         / rs.profiles.molecular_counts[0, :]
                         * rs.profiles.molecular_counts[0, inx]
                         / rs.profiles.combined_hi_counts[0, inx]]
        legend=['chi/mol']
        if hasattr(rs.profiles,'combined_lo_counts'):
            lines = [rs.profiles.combined_lo_counts[0, :]
                         / rs.profiles.molecular_counts[0, :]
                         * rs.profiles.molecular_counts[0, inx]
                         / rs.profiles.combined_lo_counts[0, inx]]
            legend=['clo/mol']
            lines.append(rs.profiles.combined_hi_counts[0, :]
                         / rs.profiles.molecular_counts[0, :]
                         * rs.profiles.molecular_counts[0, inx]
                         / rs.profiles.combined_hi_counts[0, inx])
            legend.append('chi/mol')
        else:    
            lines = [rs.profiles.combined_hi_counts[0, :]
                         / rs.profiles.molecular_counts[0, :]
                         * rs.profiles.molecular_counts[0, inx]
                         / rs.profiles.combined_hi_counts[0, inx]]
            legend=['chi/mol']
            
        if hasattr(rs.profiles,'molecular_i2a_counts'):
            lines.append(rs.profiles.molecular_i2a_counts[0, :]
                         / rs.profiles.molecular_counts[0, :]
                         * rs.profiles.molecular_counts[0, inx]
                         / rs.profiles.molecular_i2a_counts[0, inx])
            legend.append('mol_i2a/mol')             
        gt.plot_vs_altitude('dif_geo_profiles'     #display defaults plot name
                 ,instrument                       #instrument name
                 ,timerange                        #python datetimes vector
                 ,rs.profiles.msl_altitudes        #altitude vector (meters)
                 ,lines  
                 ,['c','r','g']                    #colors, [] default colors
                 ,None                             #widths, [] sets widths = 2
                 ,legend                           #legend list, [] = no legend
                 ,'upper right'                    #legend position, [] ok if list []
                 ,'ratio of chanels'               #xlabel
                 ,None                             #x units
                 ,'Differential Geometry'          #plot title
                 ,auto_loop                        # =1, clear figure before new plot
                 ,display_defaults
                 ,figs)
   
    if display_defaults.enabled('raw_profiles') and hasattr(rs,'raw_profiles')\
           and hasattr(rs.raw_profiles,'sum_molecular_counts'):
         print 'plotting raw profiles'
         N = rs.raw_profiles.sum_mean_seeded_shots 
         bin_vec = np.arange(rs.raw_profiles.sum_molecular_counts.shape[1])
         if hasattr(rs.raw_profiles,'sum_combined_wfov_counts'):
            lines = [rs.raw_profiles.sum_combined_wfov_counts[0,:]/N
                    ,rs.raw_profiles.sum_molecular_counts[0,:]/N
                    ,rs.raw_profiles.sum_combined_hi_counts[0,:]/N]
            legend = ['cwfov','mol','chi']
            colors = ['k','b','r']
            y_vecs = [bin_vec,bin_vec,bin_vec]
         elif hasattr(rs.raw_profiles,'sum_molecular_wfov_counts'):
            lines = [rs.raw_profiles.sum_molecular_wfov_counts[0,:]/N
                    ,rs.raw_profiles.sum_molecular_counts[0,:]/N
                    ,rs.raw_profiles.sum_combined_hi_counts[0,:]/N]
                    
            legend = ['mwfov','mol','chi']
            colors = ['k','b','r']
            y_vecs = [bin_vec,bin_vec,bin_vec]
         else:   
            lines = [rs.raw_profiles.sum_molecular_counts[0,:]/N
                    ,rs.raw_profiles.sum_combined_hi_counts[0,:]/N]
            legend = ['mol','chi']
            colors =['b','r']
            y_vecs = [bin_vec,bin_vec]
         if hasattr(rs.raw_profiles,'sum_combined_lo_counts'):
            lines.append(rs.raw_profiles.sum_combined_lo_counts[0,:]/N)
            legend.append('clo')
            y_vecs.append(bin_vec)
            colors.append('c')
        
         if hasattr(rs.raw_profiles,'sum_molecular_i2a_counts'):
            lines.append(rs.raw_profiles.sum_molecular_i2a_counts[0,:]/N)
            legend.append('mol_I2A')
            y_vecs.append(bin_vec)
            colors.append('m')
         
             
         if hasattr(rs.raw_profiles,'sum_combined_1064_counts'):
            lines.append(rs.raw_profiles.sum_combined_1064_counts[0,:]/N)
            legend.append('comb_1064')
            y_vecs.append(bin_vec)
            colors.append('k')
         
         if hasattr(rs.raw_profiles,'sum_cross_pol_counts'):
             lines.append(rs.raw_profiles.sum_cross_pol_counts[0,:]\
                          /rs.raw_profiles.sum_mean_seeded_shots)
             legend.append('cpol')
             y_vecs.append(bin_vec)
             colors.append('g')

         gt.plot_xy('raw_profiles'  #plot name
                ,instrument
                ,tuple([rs.raw_profiles.start_time,rs.raw_profiles.end_time]) #timerange
                ,lines
                ,y_vecs
                ,colors
                ,['None','None','None','None','None','None','None']
                ,[2,2,2,2,2,2,2]
                ,['-','-','-','-','-','-','-']
                ,[2,2,2,2,2,2,2]
                ,legend
                ,'upper right'
                ,'pileup corr counts'
                ,'counts/bin/pulse'
                ,'Bin number'
                ,None
                ,'Pileup corrected '
                ,None           #list of strings to place on plot
                ,None           #list of x-positions for text_str entries
                ,None           #list of y-positions for text_str entries
                ,None           #list of text angles for text_str entries
                ,display_defaults
                ,figs)
   
    if display_defaults.enabled('raw_dark_corrected_profiles') and hasattr(rs,'raw_profiles')\
           and hasattr(rs.raw_profiles,'sum_molecular_counts'):

         N = rs.raw_profiles.sum_mean_seeded_shots 
         bin_vec = np.arange(rs.raw_profiles.sum_molecular_counts.shape[1])
         if hasattr(rs.raw_profiles,'sum_molecular_wfov_counts'):
            lines = [(rs.raw_profiles.sum_molecular_wfov_counts[0,:]\
                           -rs.raw_profiles.sum_m_wfov_dark_counts[0,0])/N
                    ,(rs.raw_profiles.sum_molecular_counts[0,:]\
                           -rs.raw_profiles.sum_mol_dark_counts[0,0])/N
                    ,(rs.raw_profiles.sum_combined_hi_counts[0,:]\
                            -rs.raw_profiles.sum_c_hi_dark_counts[0,0])/N]
            legend = ['mwfov','mol','chi']
            colors = ['k','b','r']
            y_vecs = [bin_vec,bin_vec,bin_vec]
         else:   
            lines = [(rs.raw_profiles.sum_molecular_counts[0,:]\
                        -rs.raw_profiles.sum_mol_dark_counts[0,0])/N
                    ,(rs.raw_profiles.sum_combined_hi_counts[0,:]\
                        -rs.raw_profiles.sum_c_hi_dark_counts[0,0])/N]
            legend = ['mol','chi']
            colors =['b','r']
            y_vecs = [bin_vec,bin_vec]
         
         if hasattr(rs.raw_profiles,'sum_combined_lo_counts'):
            lines.append((rs.raw_profiles.sum_combined_lo_counts[0,:]\
                     -rs.raw_profiles.sum_c_lo_dark_counts[0,0])/N)
            legend.append('clo')
            y_vecs.append(bin_vec)
            colors.append('c')
        
         if hasattr(rs.raw_profiles,'sum_molecular_i2a_counts'):
            lines.append((rs.raw_profiles.sum_molecular_i2a_counts[0,:]\
                          -rs.raw_profiles.sum_mol_i2a_dark_counts[0,0])/N)
            legend.append('mol_I2A')
            y_vecs.append(bin_vec)
            colors.append('m')
         
             
         if hasattr(rs.raw_profiles,'sum_combined_1064_counts'):
            lines.append((rs.raw_profiles.sum_combined_1064_counts[0,:]\
                          -rs.raw_profiles.sum_combined_1064_dark_counts[0,0])/N)
            legend.append('comb_1064')
            y_vecs.append(bin_vec)
            colors.append('k')
         
         if hasattr(rs.raw_profiles,'sum_cross_pol_counts'):
             lines.append((rs.raw_profiles.sum_cross_pol_counts[0,:]\
                           -rs.raw_profiles.sum_c_pol_dark_counts[0,0])/N)
             legend.append('cpol')
             y_vecs.append(bin_vec)
             colors.append('g')
    
         gt.plot_xy('raw_dark_corrected_profiles'  #plot name
                ,instrument
                ,tuple([rs.raw_profiles.start_time,rs.raw_profiles.end_time]) #timerange
                ,lines
                ,y_vecs
                ,colors
                ,['None','None','None','None','None','None','None']
                ,[2,2,2,2,2,2,2]
                ,['-','-','-','-','-','-','-']
                ,[2,2,2,2,2,2,2]
                ,legend
                ,'upper right'
                ,'pileup corr counts'
                ,'counts/pulse/'+str(rs_constants['binwidth']*1e9)+'ns'
                ,'Bin number'
                ,None
                ,'dark corrected '
                ,None           #list of strings to place on plot
                ,None           #list of x-positions for text_str entries
                ,None           #list of y-positions for text_str entries
                ,None           #list of text angles for text_str entries
                ,display_defaults
                ,figs)
         
 



   
    if display_defaults.enabled('count_profiles') and haveProfiles \
            and hasattr(rs.profiles,'raw_molecular_counts')  and np.isfinite(rs.profiles.raw_molecular_counts).any():
        
         bin_vec = np.arange(rs.profiles.raw_molecular_counts.shape[1])
         lines = [rs.profiles.raw_molecular_counts[0,:]
                ,rs.profiles.raw_combined_hi_counts[0,:]]
         legend = ['mol','chi']
         colors =['b','r']
         y_vecs = [bin_vec,bin_vec]

         if hasattr(rs.profiles,'raw_combined_wfov_counts'):
            lines.append(rs.profiles.raw_combined_wfov_counts[0,:])
            legend.append('cwfov')
            colors.append('k')
            y_vecs.append(bin_vec)
         
         if hasattr(rs.profiles,'raw_molecular_wfov_counts'):
            lines.append(rs.profiles.raw_molecular_wfov_counts[0,:])
            legend.append('mwfov')
            colors.append('k')
            y_vecs.append(bin_vec)
         
         if hasattr(rs.profiles,'raw_combined_lo_counts'):
            lines.append(rs.profiles.raw_combined_lo_counts[0,:])
            legend.append('clo')
            y_vecs.append(bin_vec)
            colors.append('c')
        
         if hasattr(rs.profiles,'raw_molecular_i2a_counts'):
            lines.append(rs.profiles.raw_molecular_i2a_counts[0,:])
            legend.append('mol_I2A')
            y_vecs.append(bin_vec)
            colors.append('m')
         
             
         if hasattr(rs.profiles,'raw_combined_1064_counts'):
            lines.append(rs.profiles.raw_combined_1064_counts[0,:])
            legend.append('comb_1064')
            y_vecs.append(bin_vec)
            colors.append('k')
         
         
         lines.append(rs.profiles.raw_cross_pol_counts[0,:])
         legend.append('cpol')
         y_vecs.append(bin_vec)
         colors.append('g')
        
         gt.plot_xy('count_profiles'  #plot name
                ,instrument
                ,timerange
                ,lines
                ,y_vecs
                ,colors
                ,['None','None','None','None','None','None','None']
                ,[2,2,2,2,2,2,2]
                ,['-','-','-','-','-','-','-']
                ,[2,2,2,2,2,2,2]
                ,legend
                ,'upper right'
                ,'pileup corr counts'
                ,'counts/bin/pulse'
                ,'Bin number'
                ,None
                ,'Pileup corrected '
                ,None           #list of strings to place on plot
                ,None           #list of x-positions for text_str entries
                ,None           #list of y-positions for text_str entries
                ,None           #list of text angles for text_str entries
                ,display_defaults
                ,figs)
   
    if display_defaults.enabled('dark_corrected_profiles') and haveProfiles \
           and np.isfinite(rs.profiles.dc_molecular_counts).any():
        bin_vec = np.arange(rs.profiles.dc_molecular_counts.shape[1])
        rs.profiles.dc_molecular_counts[0,
                rs.profiles.dc_molecular_counts[0, :] <= 0] = np.NaN
        rs.profiles.dc_combined_hi_counts[0,
                rs.profiles.dc_combined_hi_counts[0, :] <= 0] = np.NaN
        rs.profiles.dc_cross_pol_counts[0,
                rs.profiles.dc_cross_pol_counts[0, :] <= 0] = np.NaN
        lines= [rs.profiles.dc_cross_pol_counts[0,:]
               ,rs.profiles.dc_molecular_counts[0,:]
               ,rs.profiles.dc_combined_hi_counts[0,:]]
        legend = ['c_pol','mol','chi']
        colors = ['g','b','r']
        y_vecs = [bin_vec,bin_vec,bin_vec]
        if hasattr(rs.profiles,'dc_combined_wfov_counts'):
            rs.profiles.dc_combined_wfov_counts[0,
                rs.profiles.dc_combined_wfov_counts[0, :] <= 0] = np.NaN
            lines.insert(0,rs.profiles.dc_combined_wfov_counts[0,:])
            legend.insert(0,'cwfov')
            colors.insert(0,'k')
            y_vecs.insert(0,bin_vec)
        if hasattr(rs.profiles,'dc_molecular_wfov_counts'):
            rs.profiles.dc_molecular_wfov_counts[0,
                rs.profiles.dc_molecular_wfov_counts[0, :] <= 0] = np.NaN
            lines.insert(0,rs.profiles.dc_molecular_wfov_counts[0,:])
            legend.insert(0,'mwfov')
            colors.insert(0,'k')
            y_vecs.insert(0,bin_vec)
        if hasattr(rs.profiles,'dc_combined_lo_counts'):
            rs.profiles.dc_combined_lo_counts[0,
                rs.profiles.dc_combined_lo_counts[0, :] <= 0] = np.NaN
            lines.append(rs.profiles.dc_combined_lo_counts[0,:])
            legend.append('clo')
            y_vecs.append(bin_vec)
            colors.append('c')     
        if hasattr(rs.profiles,'dc_molecular_i2a_counts'):
            rs.profiles.dc_molecular_i2a_counts[0,
                rs.profiles.dc_molecular_i2a_counts[0, :] <= 0] = np.NaN
            lines.append(rs.profiles.dc_molecular_i2a_counts[0,:])
            legend.append('mol_I2A')
            y_vecs.append(bin_vec)
            colors.append('m')
        if hasattr(rs.profiles,'dc_combined_1064_counts'):
            rs.profiles.dc_combined_1064_counts[0,
                rs.profiles.dc_combined_1064_counts[0, :] <= 0] = np.NaN
            lines.append(rs.profiles.dc_combined_1064_counts[0,:])
            legend.append('comb_1064')
            y_vecs.append(bin_vec)
            colors.append('k')
        #current_binwidth_ns = (rs.profiles.msl_altitudes[2]-rs.profiles.msl_altitudes[1])/0.15    
        current_binwidth_ns = rs_constants['binwidth'] * 1e9
        print 'Rendering dark_corrected_profiles'    
        gt.plot_xy('dark_corrected_profiles'
                ,instrument
                ,timerange
                ,lines 
                ,y_vecs
                ,colors
                ,['None','None','None','None','None']
                ,[2,2,2,2,2,2,2]
                ,['-','-','-','-','-','-']
                ,[2,2,2,2,2,2,2]
                ,legend
                ,'upper right'
                ,'Pileup, dark corr counts'
                ,'counts/pulse/'+str(current_binwidth_ns)+'ns'
                ,'Bin number'
                ,None
                ,'pileup+dark corrected counts'
                ,None           #list of strings to place on plot
                ,None           #list of x-positions for text_str entries
                ,None           #list of y-positions for text_str entries
                ,None           #list of text angles for text_str entries
                ,display_defaults
                ,figs)
        
    if display_defaults.enabled('corrected_profiles') and haveProfiles \
           and hasattr(rs.profiles,'molecular_counts') and np.isfinite(rs.profiles.molecular_counts).any():           
        bin_vec = np.arange(rs.profiles.molecular_counts.shape[1])
        np.seterr(invalid='ignore')
        rs.profiles.molecular_counts[0,
                rs.profiles.molecular_counts[0, :] <= 0] = np.NaN
        rs.profiles.combined_hi_counts[0,
                rs.profiles.combined_hi_counts[0, :] <= 0] = np.NaN
        np.seterr(invalid='warn')
        lines = [rs.profiles.molecular_counts[0,:],rs.profiles.combined_hi_counts[0,:]]
        legend = ['mol','chi']
        y_vecs = [bin_vec,bin_vec]
        if hasattr(rs.profiles,'combined_lo_counts'):
            np.seterr(invalid = 'ignore')
            rs.profiles.combined_lo_counts[0,rs.profiles.combined_lo_counts[0,:]<=0] = np.NaN
            np.seterr(invalid = 'warn')
            gain_cor_comb_lo = rs.profiles.combined_lo_counts[0,:] \
                          *rs_constants['hi_to_low_combined_channel_gain_ratio']  
            #lines.extend([temp , gain_cor_comb_lo])
            lines.extend([rs.profiles.combined_lo_counts[0,:],gain_cor_comb_lo])
            legend.extend(['clo','gain*clo'])
            y_vecs.extend([bin_vec,bin_vec])
        np.seterr(invalid = 'ignore')    
        rs.profiles.dc_cross_pol_counts[0,
                rs.profiles.cross_pol_counts[0, :] <= 0] = np.NaN
        np.seterr(invalid = 'warn')
        lines.append(rs.profiles.cross_pol_counts[0,:])
        legend.append('cpol')
        y_vecs.append(bin_vec)

        if hasattr(rs.profiles,'molecular_i2a_counts'):
            rs.profiles.molecular_i2a_counts[0,
                rs.profiles.molecular_i2a_counts[0,:] <=0] =np.NaN                              
            lines.append(rs.profiles.molecular_i2a_counts[0,:])
            legend.append('mol_I2A')
            y_vecs.append(bin_vec)

        if hasattr(rs.profiles.inv,'combined_1064_counts'):
            rs.profiles.combined_1064_counts[0,
                rs.profiles.combined_1064_counts[0,:] <=0] =np.NaN               
            lines.append(rs.profiles.combined_1064_counts[0,:])
            legend.append('1064')
            y_vecs.append(bin_vec)
            
        gt.plot_xy('corrected_profiles'
                ,instrument
                ,timerange
                ,lines
                ,y_vecs
                ,['b','r','c','k','g','m']
                ,['None','None','None','None','None']
                ,[2,2,2,2,2]
                ,['-','-','-','-','-']
                ,[2,2,2,2,2]
                ,legend
                ,'upper right'
                ,'corrected counts per laser pulse'
                ,None
                ,'Bin number'
                ,None
                ,'Fully corrected counts'
                ,None           #list of strings to place on plot
                ,None           #list of x-positions for text_str entries
                ,None           #list of y-positions for text_str entries
                ,None           #list of text angles for text_str entries
                ,display_defaults
                ,figs)

    if display_defaults.enabled('geo_correction_profile') \
             and not geo_corr == None:
         s_bin=int(rs_constants['apd_pulse_timing'][1]
                          /rs_constants['binwidth'])-1
         alt_vec = np.arange(len(geo_corr[:-s_bin]))*rs_constants['binwidth']*3e8/2
       
         gt.plot_xy('geo_correction_profile'  #plot name
                ,instrument
                ,timerange
                ,[geo_corr[:-s_bin,1]/alt_vec**2]
                ,[alt_vec/1000.0]
                ,['r','k']
                ,['None','None','None']
                ,[2,2,2]
                ,['-','-','-']
                ,[1,2,2]
                ,None  #legend
                ,'upper right'
                ,'geo correction'
                ,None
                ,'Range '
                ,'km'
                ,'geo correction'
                ,None           #list of strings to place on plot
                ,None           #list of x-positions for text_str entries
                ,None           #list of y-positions for text_str entries
                ,None           #list of text angles for text_str entries
                ,display_defaults
                ,figs)    

    if display_defaults.enabled('wfov_geo_correction') and haveProfiles\
             and hasattr(rs.profiles,'dc_molecular_wfov_counts')\
             and hasattr(rs.profiles,'dc_molecular_wfov_counts'):    
        m_wfov_counts = rs.profiles.dc_molecular_wfov_counts.copy()
        m_counts = rs.profiles.dc_molecular_counts.copy()
        m_wfov_counts[m_wfov_counts <=0] = np.NaN
        m_counts[m_counts <= 0] = np.NaN
        s_bin=int(rs_constants['apd_pulse_timing'][1]
                          /rs_constants['binwidth'])-1
        wfov_ratio = nanmean(m_counts[:,s_bin:1000],0)/nanmean(m_wfov_counts[:,s_bin:1000],0)
        alt_vec = np.arange(len(wfov_ratio))*rs_constants['binwidth']*3e8/2
        geo = geo_corr[:1000-s_bin,1]/alt_vec**2
        wfov_ratio = geo * wfov_ratio
        wfov_ratio = wfov_ratio/nanmean(wfov_ratio[600:800])
        gt.plot_xy('wfov_geo_correction'  #plot name
                ,instrument
                ,timerange
                ,[wfov_ratio,geo/geo[900]]
                ,[alt_vec/1000.0,geo_corr[:1000-s_bin,0]/1000.0]
                ,['r','k']
                ,['None','None','None']
                ,[2,2,2]
                ,['-','-','-']
                ,[1,2,2]
                ,['(mol/wfov)*geo','geo']
                ,'upper right'
                ,'geo correction'
                ,None
                ,'Range '
                ,'km'
                ,'wfov geo correction'
                ,None           #list of strings to place on plot
                ,None           #list of x-positions for text_str entries
                ,None           #list of y-positions for text_str entries
                ,None           #list of text angles for text_str entries
                ,display_defaults
                ,figs)
        
    if display_defaults.enabled('wfov_profile')\
             and haveProfiles\
             and (hasattr(rs.profiles,'dc_combined_wfov_counts')\
             or hasattr(rs.profiles,'dc_molecular_wfov_counts')):
         
         s_bin=int(rs_constants['apd_pulse_timing'][1]
                          /rs_constants['binwidth'])-1
         bin_vec = np.arange(1500-s_bin)
         if hasattr(rs.profiles,'dc_combined_wfov_counts'):
              title_str ='combined wfov profile'
              legend0= 'wfov * '+str(rs_constants['wfov_to_combined_gain_ratio'])
              legend2='chi'
              lines = [rs.profiles.dc_combined_wfov_counts[0,s_bin:1500]\
                  * rs_constants['wfov_to_combined_gain_ratio']\
                  ,rs.profiles.dc_combined_wfov_counts[0,s_bin:1500]\
                 ,rs.profiles.dc_combined_hi_counts[0,s_bin:1500]]
              colors = ['c','k','r']
         else:  #molecular wfov 
             lines = [rs.profiles.dc_molecular_wfov_counts[0,s_bin:1500]
                  * rs_constants['molecular_to_wfov_gain_ratio']
                  ,rs.profiles.dc_molecular_wfov_counts[0,s_bin:1500]
                 ,rs.profiles.dc_molecular_counts[0,s_bin:1500]]
             title_str ='wfov molecular  profile'
             legend0= 'wfov * '+str(rs_constants['molecular_to_wfov_gain_ratio'])
             legend2= 'mol'
             colors = ['c','k','b']
         alt_vec = rs.profiles.msl_altitudes[:lines[0].size]/1000.0# bin_vec * rs_constants['binwidth']*1.5e8/1000.0
         gt.plot_xy('wfov_profile'  #plot name
                ,instrument
                ,timerange
                ,lines
                ,[alt_vec, alt_vec, alt_vec]
                ,colors
                ,['None','None','None']
                ,[2,2,2]
                ,['-','-','-']
                ,[1,2,2]
                ,[legend0,'wfov',legend2]
                ,'upper right'
                ,'wfov counts'
                ,None
                ,'Range '
                ,'km'
                ,title_str
                ,None           #list of strings to place on plot
                ,None           #list of x-positions for text_str entries
                ,None           #list of y-positions for text_str entries
                ,None           #list of text angles for text_str entries 
                ,display_defaults
                ,figs)

    if display_defaults.enabled('wfov_ratios')\
             and haveProfiles\
             and (hasattr(rs.profiles,'dc_combined_wfov_counts')\
                  or hasattr(rs.profiles,'dc_molecular_wfov_counts')) and hasattr(rs.profiles,'geo_corr'):
         s_bin=int(rs_constants['apd_pulse_timing'][1]
                          /rs_constants['binwidth'])-1
         bin_vec = np.arange(1500-s_bin)
         
         if hasattr(rs.profiles,'dc_combined_wfov_counts'):
               rs.profiles.dc_molecular_wfov_counts[
                   rs.profiles.dc_molecular_wfov_counts <= 0] = np.NaN
               rs.profiles.dc_molecular_counts[
                   rs.profiles.dc_molecular_counts <=0] = np.NaN
               rs.profiles.dc_combined_wfov_counts[
                       rs.profiles.dc_combined_wfov_counts <= 0] = np.NaN
               lines = [rs.profiles.dc_combined_wfov_counts[0,s_bin:1500]
                      /rs.profiles.dc_combined_hi_counts[0,s_bin:1500]
                      ,rs.profiles.dc_combined_wfov_counts[0,s_bin:1500]
                      /rs.profiles.dc_molecular_counts[0,s_bin:1500]]   
               title_str ='combined wfov ratio profiles'
               legend_str = ['cwfov/chi','cwfov/mol']
               y_vars =[bin_vec,bin_vec] 
         else: #system with wfov molecular channel
               #conversion to NaN's to allow log plot

               np.seterr(invalid = 'ignore')
               rs.profiles.dc_molecular_wfov_counts[
                   rs.profiles.dc_molecular_wfov_counts <= 0] = np.NaN
               rs.profiles.dc_molecular_counts[
                   rs.profiles.dc_molecular_counts <=0] = np.NaN
               np.seterr(invalid = 'warn')
               
               wfov_ratio = rs.profiles.dc_molecular_wfov_counts[0,:1500] * rs_constants['molecular_to_wfov_gain_ratio'] / rs.profiles.dc_molecular_counts[0,:1500]
               pgeo=rs.profiles.geo_corr[:wfov_ratio.size]
               pgeo_ranges = rs.profiles.msl_altitudes[:wfov_ratio.size]-rs_constants['lidar_altitude']
               r2 = 1e-6 * pgeo_ranges**2
               pgeo_altitudes = rs.profiles.msl_altitudes[:wfov_ratio.size]
               lines = [wfov_ratio 
                       ,pgeo / r2
                       ,wfov_ratio * r2
                       ,pgeo
                       ,wfov_ratio / (pgeo / r2)]
               y_vars = [pgeo_altitudes/1000.0
                        ,pgeo_altitudes/1000.0
                        ,pgeo_altitudes/1000.0
                        ,pgeo_altitudes/1000.0
                        ,pgeo_altitudes/1000.0]
         title_str ='mol wfov ratio profiles'
         legend_str = ['mwfov/mol' , 'geo_corr', 'mwfov/mol * r2', 'geo_corr * r2', '(mwfov/mol) / geo_corr']
         gt.plot_xy('wfov_ratios'  #plot name
                ,instrument
                ,timerange
                ,lines
                ,y_vars
                ,['c','b','g','k','m']
                ,['None','None','None','None', 'None']
                ,[1,1,1,1,1]
                ,['-','-','-','-','-']
                ,[1,1,1,1,1]
                ,legend_str
                ,'upper right'
                ,'wfov ratios'
                ,None
                ,'Altitude (km)'
                ,None
                ,'wfov ratios'
                ,None           #list of strings to place on plot
                ,None           #list of x-positions for text_str entries
                ,None           #list of y-positions for text_str entries
                ,None           #list of text angles for text_str entries 
                ,display_defaults
                ,figs)
    if display_defaults.enabled('wfov_geo_adjust') and hasattr(rs,'profiles') \
              and hasattr(rs.profiles,'wfov_geo_adjust'):
        gt.plot_vs_altitude('wfov_geo_adjust'      #display defaults plot name
                 ,instrument                       #instrument name
                 ,timerange                        #python datetimes vector
                 ,rs.profiles.msl_altitudes[rs.profiles.msl_altitudes <= 7000.]  #altitude vector (meters)
                 ,[rs.profiles.geo_corr[rs.profiles.msl_altitudes <= 7000.]
                   ,rs.profiles.wfov_geo_adjust[rs.profiles.msl_altitudes <= 7000.]]
                 ,['b','r']                        #colors, [] default colors
                 ,[]                               #widths, [] sets widths = 2
                 ,['geo*r^2','wfov_geo*r^2']                               #legend list, [] = no legend
                 ,'upper left'                    #legend position, [] ok if list []
                 ,'geometry correction'            #xlabel
                 ,None                               #x units
                 ,'wfov geo adjust '               #plot title
                 ,auto_loop                        # =1, clear figure before new plot
                 ,display_defaults
                 ,figs)

    if display_defaults.enabled('wfov_ratio_to_geo_corr_ratio')\
             and hasattr(rs,'rs_mean')\
             and hasattr(rs.rs_mean,'molecular_wfov_counts') \
             and hasattr(rs.rs_mean,'molecular_wfov_counts') and hasattr(rs.rs_mean,'geo_corr_array'):
        range_sq = (rs.rs_mean.msl_altitudes-rs_constants['lidar_altitude'])**2
        ratio = rs_constants['molecular_to_wfov_gain_ratio'] * range_sq/1e6 \
                * (rs.rs_mean.molecular_wfov_counts/rs.rs_mean.molecular_counts) \
                / rs.rs_mean.geo_corr_array
        title_str = instrument + ' wfov_ratio to geo_corr  '
        units_str= ''
        print 'rendering wfov_ratio_to_geo_corr image'
        gt.rti_fig('wfov_ratio_to_geo_corr_ratio'
            ,instrument    
            ,ratio
            ,rs.rs_mean.times
            ,invalts
            ,title_str
            ,units_str
            ,qc_mask_1064
            ,display_defaults
            ,figs)
                                 
    if display_defaults.enabled('raw_i2a_mol_ratio')\
             and haveProfiles\
             and hasattr(rs.profiles,'dc_molecular_i2a_counts'):
         title_str ='raw i2a/mol ratio '
         xlabel = 'i2a/i2, dark corr only--includes dark interval'
         ratio = rs.profiles.dc_molecular_i2a_counts[0,:] / rs.profiles.dc_molecular_counts[0,:]
         gt.plot_xy('raw_i2a_mol_ratio'  #plot name
                ,instrument
                ,timerange
                ,[ratio]
                ,[np.arange(len(ratio))]
                ,['c','b']
                ,['None','None']
                ,[2,2]
                ,['-','-']
                ,[2,2]
                ,None
                ,'upper right'
                ,'i2a/mol---dark corr only'
                ,None
                ,'bin number'
                ,None
                ,'raw i2a/mol ratios'
                ,None           #list of strings to place on plot
                ,None           #list of x-positions for text_str entries
                ,None           #list of y-positions for text_str entries
                ,None           #list of text angles for text_str entries
                ,display_defaults
                ,figs)
        
         
    if display_defaults.enabled('i2a_mol_ratio')\
             and haveProfiles\
             and hasattr(rs.profiles,'i2a_mol_ratio'):
         filter_window =np.int(processing_defaults.get_value('i2a_mol_ratio','filter_window'
              ))
         title_str ='i2a/mol ratio '
         xlabel = 'i2a/i2 ratio '+ str(filter_window) + 'm filter'
         gt.plot_vs_altitude('i2a_mol_ratio'       #display defaults plot name
                 ,instrument                       #instrument name
                 ,timerange                        #python datetimes vector
                 ,rs.profiles.msl_altitudes        #altitude vector (meters)
                 ,[rs.profiles.i2a_mol_ratio[0,:]] #variables
                 ,['r' ]                           #colors, [] default colors
                 ,None                             #widths, [] sets widths = 2
                 ,None                             #legend list, [] = no legend
                 ,'upper right'                    #legend position, [] ok if list []
                 ,xlabel                           #xlabel
                 ,None                             #x units
                 ,title_str                        #plot title
                 ,auto_loop                        # =1, clear figure before new plot
                 ,display_defaults                    
                 ,figs)
         
    if display_defaults.enabled('calibration_coefficients') and rs_Cxx is not None: #FIXME this is a cal profile
         title_str ='Calibration Coefficents '
         xlabel = 'Cal Coef'
         nalts = len(rs_Cxx.msl_altitudes)
         lines = [rs_Cxx.Cmm[:nalts], rs_Cxx.Cmc[:nalts]]
         legend= ['Cmm','Cmc']
         if hasattr(rs_Cxx,'Cmm_i2a'):
             lines.append(rs_Cxx.Cmm_i2a[:nalts])
             legend.append('Cmm_i2a')
         gt.plot_vs_altitude('calibration_coefficients'        #display defaults plot name
                 ,instrument                       #instrument name
                 ,rs_Cxx.sounding_time          #python datetimes vector
                 ,rs_Cxx.msl_altitudes          #altitude vector (meters)
                 ,lines                            #variables
                 ,['b','r','k']                    #colors, [] default colors
                 ,None                             #widths, [] sets widths = 2
                 ,legend                           #legend list, [] = no legend
                 ,'upper right'                    #legend position, [] ok if list []
                 ,xlabel                           #xlabel
                 ,None                               #x units
                 ,title_str                        #plot title
                 ,auto_loop                        # =1, clear figure before new plot
                 ,display_defaults                    
                 ,figs)
  
    if display_defaults.enabled('photon_counting_error_extinction')\
            and hasattr(rs,'rs_mean') and  hasattr(rs.rs_mean,'var_raw_molecular_counts'):
         title_str ='photon count extinction error '
         xlabel = 'Extinction error'
         delta_r = rs.rs_mean.msl_altitudes[2]-rs.rs_mean.msl_altitudes[1]
         dt =rs.rs_mean.delta_t[1]
         nalts = len(rs.rs_mean.msl_altitudes)
         ext_error = np.sqrt(2.0)*np.sqrt(nanmean(rs.rs_mean.var_raw_molecular_counts,0))\
                     /(nanmean(rs.rs_mean.molecular_counts,0) * delta_r)
         lines = [ext_error,ext_error/10.0,ext_error/100.0]
         legend =['dt='+str(dt)+'s dr='+str(delta_r)+'m', 'dt='+str(dt)+'s dr='+str(10*delta_r)+'m'\
                    ,'dt='+str(dt*10)+'s dr='+str(10*delta_r)+'m']
       
         gt.plot_vs_altitude('photon_counting_error_extinction'        #display defaults plot name
                 ,instrument                       #instrument name
                 ,timerange               #python datetimes vector
                 ,rs.rs_mean.msl_altitudes         #altitude vector (meters)
                 ,lines                            #variables
                 ,['b','r','g']                    #colors, [] default colors
                 ,None                             #widths, [] sets widths = 2
                 ,legend                           #legend list, [] = no legend
                 ,'upper left'                    #legend position, [] ok if list []
                 ,xlabel                           #xlabel
                 ,'1/m'                               #x units
                 ,title_str                        #plot title
                 ,auto_loop                        # =1, clear figure before new plot
                 ,display_defaults                    
                 ,figs)
         
    if display_defaults.enabled('geometry_correction') and hasattr(rs,'profiles') \
           and hasattr(rs.profiles,'geo_corr') and rs_Cxx is not None:
      try:
         title_str ='Geometry correction '
         xlabel = 'Geometry correction'
         
         #altitudes = geo_corr[:,0]/1000.0
         geo = rs.profiles.geo_corr.copy()
         geo[np.isnan(geo)]=0.0
         geo[geo <=0] = np.NaN
         #geo[rs.profiles.msl_altitudes < rs_constants['lidar_altitude']] =np.NaN

         r = rs.profiles.msl_altitudes - rs_constants['lidar_altitude']
         
         geo_shift = geo.copy()
         half_shift = 2
         shift = np.int(half_shift * 2.0)
         geo_shift[shift:] = geo_shift[:-shift]
         d_log_geo = np.ones_like(geo)*np.NaN
         d_log_geo[half_shift:-half_shift] = \
              -(geo_shift[half_shift:-half_shift]-geo[half_shift:-half_shift]) \
              /(geo[half_shift:-half_shift] * shift*(rs.profiles.msl_altitudes[shift+1]-rs.profiles.msl_altitudes[1]))
         d_log_geo[d_log_geo<=0]=np.NaN
         d_log_geo[rs.profiles.msl_altitudes <= rs_constants['lidar_altitude']] = np.NaN
         lines = [geo,geo / (r**2 *1e-6),d_log_geo]
         legend= ['geo * r2','geo','d_log_geo * r2']
         
         if hasattr(rs.profiles,'wfov_geo_adjust'):
             temp = rs.profiles.wfov_geo_adjust.copy()
             temp[rs.profiles.msl_altitudes <= rs_constants['lidar_altitude']]=np.NaN
             temp[temp <= 0] =np.NaN
             lines.append(temp)
             
             legend.append('wfov_corr')
         gt.plot_vs_altitude('geometry_correction'        #display defaults plot name
                 ,instrument                       #instrument name
                 ,rs.profiles.times                #python datetimes vector
                 ,rs.profiles.msl_altitudes        #altitude vector (meters)
                 ,lines                            #variables
                 ,['b','g','r','k']                #colors, [] default colors
                 ,None                             #widths, [] sets widths = 2
                 ,legend                           #legend list, [] = no legend
                 ,'upper right'                    #legend position, [] ok if list []
                 ,xlabel                           #xlabel
                 ,None                             #x units
                 ,title_str                        #plot title
                 ,auto_loop                        # =1, clear figure before new plot
                 ,display_defaults                    
                 ,figs) 
      except:
        print 'GEOCORR PLOT Broke'
        traceback.print_exc()
    """     
    if display_defaults.enabled('geometry_correction') and timerange!=None and geo_corr!=None:#FIXME    
         title_str ='Geometry correction '
         xlabel = 'Geometry correction'
         
         altitudes = geo_corr[:,0]/1000.0
         geo = geo_corr[:,1].copy()
         geo[geo <=0] = np.NaN
       
         geo_shift = geo.copy()
         half_shift = 10
         shift = np.int(half_shift * 2.0)
         geo_shift[shift:] = geo_shift[:-shift]
         d_log_geo = np.ones_like(geo)*np.NaN
         d_log_geo[half_shift:-half_shift] = -(geo_shift[half_shift:-half_shift]-geo[half_shift:-half_shift])/(geo[half_shift:-half_shift] * shift*(altitudes[shift+1]-altitudes[1]))
         d_log_geo[d_log_geo<=0]=np.NaN
        
         #r_corr_d_log_geo = -d_log_geo+1.0/altitudes
         #r_corr_d_log_geo[r_corr_d_log_geo<=0] = np.NaN
         lines = [geo,d_log_geo]
         legend= ['geo * r2','d_log_geo * r2']
         alt_list = [altitudes,altitudes]
         if hasattr(rs,'rs_mean') and hasattr(rs.mean,'wfov_geo_adjust'):
             lines.append(rs.rs_mean.wfov_geo_adjust*rs.profiles.geo_corr)
             alt_list.append[rs.profiles.msl_altitudes]
         gt.plot_xy('geometry_correction'  #plot name
                ,instrument
                ,timerange
                ,lines
                ,alt_list
                ,['c','b']
                ,['None','None']
                ,[2,2]
                ,['-','-','-']
                ,[2,2]
                ,legend
                ,'upper right'
                ,'Geometry Correction'
                ,None
                ,'altitude'
                ,'km'
                ,'Geometry correction'
                ,[]
                ,[]
                ,[]
                ,[]
                ,display_defaults
                ,figs)
       """ 
    #wfov ratio
    if display_defaults.enabled('filtered_wfov_mol_ratio_vs_time') and hasattr(rs,'rs_mean') \
          and hasattr(rs.rs_mean,'filtered_wfov_mol_ratio'):
       print 'rendering filtered_wfov_mol_ratio_vs_time'

       
       lines = []
       widths = []
       colors = []
       legends = []
       [lines,colors,widths,legends]=define_multiple_alt_plots(
                 rs.rs_mean.filtered_wfov_mol_ratio
                ,rs.rs_mean.msl_altitudes
                ,plot_alt_index
                ,lines
                ,widths
                ,colors
                ,legends)
      
       if 1:
           gt.plot_vs_time('filtered_wfov_mol_ratio_vs_time'
                 ,instrument
                 ,rs.rs_mean.times
                 ,lines
                 ,colors
                 ,widths
                 ,legends
                 ,'upper left'    
                 ,'(wfov/mol)/calib(wfov/mol)'
                 ,None
                 ,'wfov/mol ratio vs time'   
                 ,auto_loop
                 ,display_defaults
                 ,figs) 
       if 0: #except:
           print 'error ploting wfov_mol_ratio'
        
   # show other figures
    if display_defaults.enabled('depol_backscat_hist') and hasattr(rs,'rs_inv'):
       if 1:
           depol = 100*rs.rs_inv.linear_depol
           beta_a_backscat = rs.rs_inv.beta_a_backscat_par \
                             + rs.rs_inv.beta_a_backscat_perp
           if rs.rs_inv.qc_mask is not None:
             mask = np.bitwise_and(rs.rs_inv.qc_mask,1)
             print 'qc_mask applied to depolarization-backscatter histogram plot' 
             beta_a_backscat[mask==0]=np.NaN
             depol[mask==0] = np.NaN
           title= 'depol vs backscat'
           xlabel='Backscatter cross section'
           x_units = '1/(m sr)'
           ylabel='Depol'
           y_units = '%'
           gt.plot_2d_histogram('depol_backscat_hist'
                          ,instrument
                          ,timerange#.times
                          ,[beta_a_backscat]
                          ,[depol]   #100*rs.rs_inv.linear_depol]
                          ,None                          #no optional second plot
                          ,None                          #no legend for second plot
                          ,None                            #no position for legend      
                          ,xlabel
                          ,x_units
                          ,ylabel
                          ,y_units
                          ,title
                          ,display_defaults
                          ,figs)
       else:
           print 'no data for depol vs backscatter histogram'
    if display_defaults.enabled('color_ratio_backscat_hist') and hasattr(rs,'rs_inv')\
         and hasattr(rs.rs_inv,'color_ratio'):
       if 1:
           beta_a_backscat = rs.rs_inv.beta_a_backscat_par \
                             + rs.rs_inv.beta_a_backscat_perp
           color_ratio = rs.rs_inv.color_ratio.copy()
           if qc_mask is not None:
             print 'qc_mask applied to color ratio - 532nm backscatter histogram plot' 
             beta_a_backscat[qc_mask == 0]=np.NaN
             color_ratio[qc_mask == 0] = np.NaN
           title= 'color_ratio vs backscat, %4.2f-->%4.2f km'\
                  %(invalts[layer_indices[0]]/1000.0,invalts[layer_indices[-1]]/1000.0)
           xlabel='532nm Backscatter cross section'
           x_units = '1/(m sr)'
           ylabel='1064/532 backscat ratio'
           y_units = '%'
           gt.plot_2d_histogram('color_ratio_backscat_hist'
                          ,instrument
                          ,timerange#.times
                          ,[beta_a_backscat[:,layer_indices]]
                          ,[color_ratio[:,layer_indices]]
                          ,None                          #no optional second plot
                          ,None                            #no legend for second plot
                          ,None                            #no position for legend      
                          ,xlabel
                          ,x_units
                          ,ylabel
                          ,y_units
                          ,title
                          ,display_defaults
                          ,figs)
       else:
           print 'no data for color ratio vs 532nm backscatter histogram'
           
    if display_defaults.enabled('color_ratio_depol_hist') and hasattr(rs,'rs_inv')\
         and hasattr(rs.rs_inv,'color_ratio'):
       if 1:
           depol = rs.rs_inv.linear_depol[:,layer_indices].copy()
           color_ratio = rs.rs_inv.color_ratio[:,layer_indices].copy()
           if qc_mask is not None:
             print 'qc_mask applied to color_ratio depol histogram plot' 
             color_ratio[qc_mask[:,layer_indices] == 0]=np.NaN
             depol[qc_mask[:,layer_indices] == 0] = np.NaN
           title= 'color_ratio vs depol, %4.2f-->%4.2f km'\
                  %(invalts[layer_indices[0]]/1000.0,invalts[layer_indices[-1]]/1000.0)
           xlabel='532nm depolarization ratio'
           x_units = '%'
           ylabel='1064/532 backscatter ratio'
           y_units = '%'
         
           gt.plot_2d_histogram('color_ratio_depol_hist'
                          ,instrument
                          ,timerange
                          ,[100*depol]
                          ,[color_ratio]
                          ,None                          #no optional second plot
                          ,None                            #no legend for second plot
                          ,None                            #no position for legend      
                          ,xlabel
                          ,x_units
                          ,ylabel
                          ,y_units
                          ,title
                          ,display_defaults
                          ,figs)
       else:
           print 'no data for color ratio vs depol histogram'
           
    if display_defaults.enabled('p180_backscat_hist') and hasattr(rs,'rs_inv') and hasattr(rs.rs_inv,'p180'):

            beta_a_backscat = rs.rs_inv.beta_a_backscat[:,layer_indices].copy()         
            p180 = rs.rs_inv.p180[:,layer_indices].copy()
            if qc_mask is not None:
               print
               print 'qc_mask applied to p180-backscatter histogram plot'
               p180[qc_mask[:,layer_indices] == 0] = np.NaN
               beta_a_backscat[qc_mask[:,layer_indices] == 0] = np.NaN
            title= 'p180 vs backscat %4.2f-->%4.2f km '\
                   %(invalts[layer_indices[0]]/1000.0,invalts[layer_indices[-1]]/1000.0)
            xlabel='Backscatter cross section'
            x_units= '1/(m sr)'
            ylabel='P180/4pi'
            y_units = '1/sr'
            gt.plot_2d_histogram('p180_backscat_hist'
                          ,instrument
                          ,timerange
                          ,[beta_a_backscat]
                          ,[p180]
                          ,None                   #no optional second plot       
                          ,None                   #no legend for second plot
                          ,None                     #no legend position
                          ,xlabel
                          ,x_units
                          ,ylabel
                          ,y_units
                          ,title
                          ,display_defaults
                          ,figs)

    if display_defaults.enabled('lidarRatio_backscat_hist') and hasattr(rs,'rs_inv') and hasattr(rs.rs_inv,'extinction') and hasattr(rs.rs_inv,'beta_r_backscat'):
        
            #inv.extinction is created in filtered_extinction.py
            beta_a_backscat = rs.rs_inv.beta_a_backscat_par[:,layer_indices] \
                            + rs.rs_inv.beta_a_backscat_perp[:,layer_indices]
     
            lidarRatio=(rs.rs_inv.extinction[:,layer_indices]
                   -rs.rs_inv.beta_r_backscat[layer_indices]*8*np.pi/3)\
                   /beta_a_backscat
            #apply qc_mask
            if qc_mask is not None:
               print 'qc_mask applied to lidarRatio-backscatter histogram plot'
               lidarRatio[qc_mask[:,layer_indices] == 0] = np.NaN
               beta_a_backscat[qc_mask[:,layer_indices] == 0] = np.NaN
            title= 'lidarRatio vs backscat %4.2f-->%4.2f km '\
                   %(invalts[layer_indices[0]]/1000.0,invalts[layer_indices[-1]]/1000.0)
            xlabel='Backscatter cross section'
            x_units= '1/(m sr)'
            ylabel='lidar_ratio'
            y_units = 'sr'
            gt.plot_2d_histogram('lidarRatio_backscat_hist'
                          ,instrument
                          ,timerange
                          ,[beta_a_backscat]
                          ,[lidarRatio]
                          ,None                    #no optional second plot
                          ,None                  #no legend for second plot
                          ,None                      #no legend position
                          ,xlabel
                          ,x_units
                          ,ylabel
                          ,y_units
                          ,title
                          ,display_defaults
                          ,figs)


    if display_defaults.enabled('interferometer_snapshot') and \
       hasattr(rs,'rs_raw') and hasattr(rs.rs_raw,'interf_snapshot'):
         
        gt.plot_image('interferometer_snapshot'
                   ,instrument
                   ,rs.rs_raw.interferometer_snapshot_time if hasattr(rs.rs_raw,'interferometer_snapshot_time') else timerange
                   ,rs.rs_raw.interf_snapshot
                   ,None
                   ,'interferometer'
                   ,display_defaults
                   ,figs)
        
    if display_defaults.enabled('interferometer_vs_time') and \
        hasattr(rs,'rs_raw') and hasattr(rs.rs_raw,'interferometer_intensity'):
           tmp=np.array(rs.rs_raw.interferometer_intensity,dtype='float64')
           tmp-=tmp.min()

           gt.rti_fig('interferometer_vs_time'
               ,instrument    
               ,tmp/tmp.max()
               ,rs.rs_raw.times
               ,np.arange(rs.rs_raw.interferometer_intensity.shape[1])
               ,'Interferometer Image'               
               ,None
               ,None       
               ,display_defaults
               ,figs
               ,alt_units='unitless',alt_label='')

    if  display_defaults.enabled('interferometer_spectrum') and \
        hasattr(rs,'rs_raw') and hasattr(rs.rs_raw,'interferometer_intensity'):
        npixels=rs_constants['interferometer_fft_npixels']
        xform = np.fft.rfft(rs.rs_raw.interferometer_intensity[0,:npixels], axis=0)
        interf_spectrum = np.real(xform*xform.conj())
        spectral_bins =np.arange(xform.shape[0])
        
        gt.plot_vs_x('interferometer_spectrum'
                ,instrument
                ,timerange
                ,spectral_bins
                ,[interf_spectrum]
                ,['r']
                ,[2]
                ,None
                ,'lower left'
                ,'spectral component'
                ,None
                ,'power'
                ,None
                ,'Interferometer spectrum'
                ,auto_loop
                ,display_defaults
                ,figs)
                
            
    if display_defaults.enabled('snowscope_snapshot') and \
       hasattr(rs,'rs_raw') and hasattr(rs.rs_raw,'overhead_snapshot'):
        gt.plot_image('overhead_snapshot'
                   ,instrument
                   ,rs.rs_raw.overhead_snapshot_time if hasattr(rs.rs_raw,'overhead_snapshot_time') else timerange
                   ,rs.rs_raw.overhead_snapshot
                   ,'Greys'
                   ,'overhead'
                   ,display_defaults
                   ,figs,invert=True)

    if display_defaults.enabled('snowscope_snapshot') and \
       hasattr(rs,'rs_raw') and hasattr(rs.rs_raw,'snowscope_snapshot'):

        gt.plot_image('snowscope_snapshot'
                   ,instrument
                   ,rs.rs_raw.snowscope_snapshot_time if hasattr(rs.rs_raw,'snowscope_snapshot_time') else timerange
                   ,rs.rs_raw.snowscope_snapshot
                   ,'Greys'
                   ,'window'
                   ,display_defaults
                   ,figs,invert=True)
   
    if display_defaults.enabled('int_backscat_vs_time') and hasattr(rs,'rs_inv'):
        #plots optical depth from integrated backscatter between top and bottom of requested layer
        ylabel = 'Aerosol Optical Depth  %5.2f -->%5.2f km' %(
                 invalts[layer_indices[0]]/1000.0
                ,invalts[layer_indices[-1]]/1000.0)  
        lines = []
        legend = []
        colors =[]
        linewidth = []
       
        if hasattr(rs.rs_inv,'optical_depth_aerosol'):
             lines.append((rs.rs_inv.optical_depth_aerosol[:, layer_indices[-1]]
                          - rs.rs_inv.optical_depth_aerosol[:, layer_indices[0]]))
             legend.append('layer AOD')
             colors.append('b')
             linewidth.append(3)

        clist = ['c','k','r','g','y','m']
        lidar_ratios = display_defaults.get_value('int_backscat_vs_time','lidar_ratio_list')

        if lidar_ratios!=None:
          for i in range(len(lidar_ratios)):
            lines.append((rs.rs_inv.integrated_backscatter[:, layer_indices[-1]]\
                          - rs.rs_inv.integrated_backscatter[:, layer_indices[0]])*lidar_ratios[i])
            legend.append('LR='+str(lidar_ratios[i]))
            colors.append(clist[i])
            linewidth.append(1)
       
        gt.plot_vs_time('int_backscat_vs_time'
                 ,instrument
                 ,rs.rs_inv.times
                 ,lines
                 ,colors
                 ,linewidth
                 ,legend
                 ,'upper left'    
                 ,ylabel
                 ,None
                 ,'integrated backscatter vs time'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)

  

    if display_defaults.enabled('od_from_int_backscat_vs_time') and hasattr(rs,'rs_inv'):
        #compute aerosol optical depth between bottom layer and other selected levels
        LR = display_defaults.get_value('od_from_int_backscat_vs_time','lidar_ratio')
        ods = rs.rs_inv.integrated_backscatter.copy()
        for i in range(len(plot_alt_index)):
            ods[:,plot_alt_index[i]] = (ods[:,plot_alt_index[i]]
                          - rs.rs_inv.integrated_backscatter[:,plot_alt_index[0]]) * LR
        lines = []
        widths = []
        colors = []
        legends = []
        [lines,colors,widths,legends]=define_multiple_alt_plots(
                ods
                ,rs.rs_inv.msl_altitudes
                ,plot_alt_index
                ,lines
                ,widths
                ,colors
                ,legends)
        ylabel = 'Aerosol Optical Depth  %5.2f -->%5.2f km' %(
                 invalts[plot_alt_index[0]]/1000.0
                ,invalts[plot_alt_index[-1]]/1000.0)  
        #add layer OD computed from slope of molecular return
        if hasattr(rs.rs_inv,'optical_depth_aerosol'):
             lines.append((rs.rs_inv.optical_depth_aerosol[:, plot_alt_index[-1]]
                          - rs.rs_inv.optical_depth_aerosol[:, plot_alt_index[0]]))
             legends.append('OD '+str(rs.rs_inv.msl_altitudes[plot_alt_index[-1]]/1000.0)+'km')
             colors.append('c')
             widths.append(3)
             
        gt.plot_vs_time('od_from_int_backscat_vs_time'
                 ,instrument
                 ,rs.rs_inv.times
                 ,lines
                 ,colors
                 ,widths
                 ,legends
                 ,'upper left'    
                 ,ylabel
                 ,None
                 ,'OD from integ backscat, LR=' + str(LR)   
                 ,auto_loop
                 ,display_defaults
                 ,figs)

    if display_defaults.enabled('mol/comb_counts_vs_time') and hasattr(rs,'rs_mean') \
                and hasattr(rs.rs_mean,'molecular_i2a_counts') and shared_dz is not None:
    
         title ='mol/comb counts,averaged bin %5.3f -->%5.3f' \
                 %(rs.rs_mean.msl_altitudes[layer_indices[0]]/1000.0,rs.rs_mean.msl_altitudes[layer_indices[-1]]/1000.0)
    
         legend =['mol/comb','i2a/comb','i2a/mol']
         comb=rs.rs_mean.combined_counts if hasattr(rs.rs_mean,'combined_counts') else rs.rs_mean.combined_hi_counts
        
         gt.plot_vs_time('mol/comb_counts_vs_time'
                 ,instrument
                 ,rs.rs_mean.times
                 ,[nanmean(rs.rs_mean.molecular_counts[:,layer_indices],1) 
                       /nanmean(comb[:,layer_indices],1)
                       ,nanmean(rs.rs_mean.molecular_i2a_counts[:,layer_indices],1) 
                       /nanmean(comb[:,layer_indices],1)
                       ,nanmean(rs.rs_mean.molecular_i2a_counts[:,layer_indices],1)
                       /nanmean(rs.rs_mean.molecular_counts[:,layer_indices],1)]
                 ,['r','g','b']
                 ,[2,2,2]
                 ,legend
                 ,'upper left'    
                 ,'ratio to combined counts'
                 ,None
                 ,title  
                 ,auto_loop
                 ,display_defaults
                 ,figs)
       
    if display_defaults.enabled('lapse_rate') and sounding is not None \
                   and hasattr(sounding,'times'):

        
        npts = sounding.altitudes.shape[0]-1
        lapse_rate = np.zeros_like(sounding.altitudes)    
        lapse_rate[1:npts-1] = sounding.temps[2:npts]\
                       -sounding.temps[:npts-2]
        lapse_rate[1:npts-1] = 1000.0 \
               * lapse_rate[1:npts-1]/(sounding.altitudes[2:npts]\
                                -sounding.altitudes[0:npts-2])  
         
        max_alt_km = display_defaults.get_value('lapse_rate','max_alt_km')
        max_alt_km = np.float(max_alt_km)
        max_alt_index = sounding.altitudes[sounding.altitudes 
                        <= max_alt_km *1000.0].shape[0]
              
        gt.plot_vs_altitude('lapse rate'
                 ,''
                 ,sounding.times
                 ,sounding.altitudes[:max_alt_index]                    
                 ,[lapse_rate[:max_alt_index]]
                 ,'r'
                 ,[2]
                 ,None
                 ,'upper right'
                 ,'Lapse rate'
                 ,'deg K / km'
                 ,'Lapse rate'
                 ,auto_loop
                 ,display_defaults
                 ,figs)

                   

    if display_defaults.enabled('wfov_ratio_vs_time')\
             and hasattr(rs,'rs_mean') and hasattr(rs.rs_mean,'filtered_wfov_mol_counts'):
        print 'rendering wfov_ratio_vs_time plot'
        try:    
            [lines,colors,widths,legends]=define_multiple_alt_plots(
                   rs.rs_mean.filtered_wfov_mol_ratio
                  ,rs.rs_mean.msl_altitudes
                  ,plot_alt_index
                  ,lines
                  ,widths
                  ,colors,legends)

           
            gt.plot_vs_time('wfov_ratio_vs_time'
                 ,instrument
                 ,rs.rs_mean.times
                 ,lines  
                 ,colors 
                 ,widths 
                 ,legends
                 ,'upper left'
                 ,'meas(wfov / mol) / (cal wfov/mol)'
                 ,None
                 ,'wfov ratio vs time'
                 ,auto_loop
                 ,display_defaults
                 ,figs)
        except:
            print
            print 'no  wfov_ratio_vs_time plot--likely no data at requested altitudes'
            print

    if display_defaults.enabled('wfov_delta_extinction_vs_time')\
             and hasattr(rs,'rs_mean')\
             and hasattr(rs.rs_mean,'wfov_extinction_corr') and shared_dz is not None:

       
        print 'rendering delta extinction vs time plot'
        title ='mol/wfov as delta extinction'

        lines = []
        widths = []
        colors = []
        legends = []
        [lines,colors,widths,legends]=define_multiple_alt_plots(
               rs.rs_mean.wfov_extinction_corr * 1e6
               ,rs.rs_mean.msl_altitudes
               ,plot_alt_index
               ,lines
               ,widths
               ,colors
               ,legends)
              
        gt.plot_vs_time('wfov_delta_extinction_vs_time'
                 ,instrument
                 ,rs.rs_mean.times
                 ,lines
                 ,colors     #['r','k']
                 ,widths     #[2,2]
                 ,legends
                 ,'upper left'
                 ,'extinction correction'
                 ,'1/m *1e6'
                 ,title
                 ,auto_loop
                 ,display_defaults
                 ,figs)
        
    if display_defaults.enabled('scattering_ratio_vs_time') \
           and hasattr(rs,'rs_inv'):
        print 'Rendering scattering_ratio_vs_time'

        if display_defaults.get_value('scattering_ratio_vs_time','ave_over_layer'):
           #plot layer average scattering ratio' 
           title ='Scattering ratio, %4.2f-->%4.2f km,' \
                %(rs.rs_inv.msl_altitudes[layer_indices[0]]/1000.0\
                   ,rs.rs_inv.msl_altitudes[layer_indices[-1]]/1000.0)
           sc_ratio_par =  np.sum(rs.rs_inv.Na[:, layer_indices[0]:layer_indices[-1]],1) \
                       / np.sum(rs.rs_inv.Nm_i2[:,layer_indices[0]:layer_indices[-1]],1)
        
           alt_index = [np.int((layer_indices[0] + layer_indices[-1])/2)]       
          
                        
           sc_ratio_perp =  np.sum(rs.rs_inv.Ncp[:, layer_indices[0]:layer_indices[-1]],1) \
                       / np.sum(rs.rs_inv.Nm_i2[:,layer_indices[0]:layer_indices[-1]],1)
           lines = [sc_ratio_par,sc_ratio_perp]
           legends = ['par','perp']
           colors = ['r','g']
           widths = [2,2]
                 
        else:          
            title ='Scattering ratio vs time'
            sc_ratio_par =  rs.rs_inv.Na / rs.rs_inv.Nm_i2
            lines = []
            widths = []
            colors = []
            legends = []
            [lines,colors,widths,legends]=define_multiple_alt_plots(
               sc_ratio_par
               ,rs.rs_inv.msl_altitudes
               ,plot_alt_index
               ,lines
               ,widths
               ,colors
               ,legends)
            for i in range(len(plot_alt_index)):
                legends[i]= legends[i]+' par'
                
            sc_ratio_perp =  rs.rs_inv.Ncp / rs.rs_inv.Nm_i2        
            [lines,colors,widths,legends]=define_multiple_alt_plots(
               sc_ratio_perp
               ,rs.rs_inv.msl_altitudes
               ,plot_alt_index
               ,lines
               ,widths
               ,colors
               ,legends)
            
            for i in range(len(plot_alt_index)):
                legends[i+len(plot_alt_index)] = legends[i+len(plot_alt_index)]+' perp'
               
        gt.plot_vs_time('scattering_ratio_vs_time'
                 ,instrument
                 ,rs.rs_inv.times
                 ,lines
                 ,colors        #['r','g','k','b']
                 ,widths        #[2,1,2,1]
                 ,legends
                 ,'upper left'    
                 ,'Scattering ratio'
                 ,None
                 ,title  
                 ,auto_loop
                 ,display_defaults
                 ,figs)
        
    if display_defaults.enabled('od_vs_time') and hasattr(rs,'rs_inv'):

        lines = []
        widths = []
        colors = []
        legends = []
       
        [lines,colors,widths,legends]=define_multiple_alt_plots(
            rs.rs_inv.optical_depth_aerosol[:,:]
            ,rs.rs_inv.msl_altitudes
            ,plot_alt_index
            ,lines
            ,widths
            ,colors
            ,legends)                                                        
              
        if display_defaults.get_value('od_vs_time','show_od_difference'):
            od = rs.rs_inv.optical_depth_aerosol[:, plot_alt_index[-1]]\
                          - rs.rs_inv.optical_depth_aerosol[:, plot_alt_index[0]]
            if qc_mask is not None:
              od[qc_mask[:,plot_alt_index[0]]==0] =np.NaN
            lines.append(od)
            colors.append('c')
            widths.append(3)
            legends.append('delta_od') 
        gt.plot_vs_time('od_vs_time'
                 ,instrument
                 ,rs.rs_inv.times
                 ,lines
                 ,colors          
                 ,widths
                 ,legends
                 ,'upper left'    
                 ,'Aerosol optical depth'
                 ,None
                 ,'Aerosol optical depth vs time'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)   
         
    if display_defaults.enabled('mol_ref_AOD_vs_time') and hasattr(rs,'rs_inv'):
        print 'rendering mol_ref_AOD_vs_time'
        gt.plot_vs_time('mol_ref_AOD_vs_time'
                 ,instrument
                 ,rs.rs_inv.times
                 ,[rs.rs_inv.mol_ref_aod]
                 ,['r']          
                 ,[2]
                 ,None
                 ,'upper left'    
                 ,'AOD at mol ref alt'
                 ,None
                 ,'AOD at ref alt vs time'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
        
    if display_defaults.enabled('extinction_vs_time') and hasattr(rs,'rs_inv') and hasattr(rs.rs_inv,'extinction_aerosol'):
        lines = []
        widths = []
        colors = []
        legends = []
        [lines,colors,widths,legends]=define_multiple_alt_plots(
            rs.rs_inv.extinction_aerosol
            ,rs.rs_inv.msl_altitudes
            ,plot_alt_index
            ,lines
            ,widths
            ,colors
            ,legends)
        
        gt.plot_vs_time('extinction_vs_time'
                 ,instrument
                 ,rs.rs_inv.times
                 ,lines
                 ,colors
                 ,widths
                 ,legends
                 ,'upper left'    
                 ,'Aerosol extinction (1/m)'
                 ,None
                 ,'Aerosol extinction vs time'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
        
    if display_defaults.enabled('backscatter_vs_time') and hasattr(rs,'rs_inv') and hasattr(rs.rs_inv,'beta_a_backscat'):
        lines = []
        widths = []
        colors = []
        legends = []
        [lines,colors,widths,legends]=define_multiple_alt_plots(
            rs.rs_inv.beta_a_backscat
            ,rs.rs_inv.msl_altitudes
            ,plot_alt_index
            ,lines
            ,widths
            ,colors
            ,legends)
        
        gt.plot_vs_time('backscatter_vs_time'
                 ,instrument
                 ,rs.rs_inv.times
                 ,lines
                 ,colors
                 ,widths
                 ,legends
                 ,'upper left'    
                 ,'Aerosol backscatter'
                 ,'1/(m sr)'
                 ,'Aerosol backscatter vs time'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
        
    if display_defaults.enabled('counts_vs_time') and hasattr(rs,'rs_raw')\
           and hasattr(rs.rs_raw,'molecular_counts'):
        title= 'raw counts vs time'
        counts=[]
        legend=[]
        colors=[]
        widths=[]
        wid=3
        requested_altitudes = display_defaults.get_value('counts_vs_time','altitudes') 
        ddz = rs_constants['binwidth']*1.5e8
        raw_alts = ddz * np.arange(len(rs.rs_raw.molecular_counts[0,:]))
       
        for i in range(len(requested_altitudes)):            
            req_altitude = np.float(requested_altitudes[i])
            index = len(raw_alts[raw_alts <= req_altitude *1000.0 \
                        + rs_constants['lidar_altitude']])
            if index < rs.rs_raw.combined_hi_counts.shape[1]:
                #if energy normalization is requested
                if display_defaults.get_value('counts_vs_time','energy_normalized'):
                    counts.append((rs.rs_raw.combined_hi_counts[:,index]
                                 / rs.rs_raw.transmitted_power))
                    counts.append((rs.rs_raw.molecular_counts[:, index]
                                 / rs.rs_raw.transmitted_power))
                    counts.append((rs.rs_raw.cross_pol_counts[:,index]
                               / rs.rs_raw.transmitted_power))
                    ylabel = 'Counts/bin/mW'
                else: #shot normalized   
                    counts.append((rs.rs_raw.combined_hi_counts[:,index]
                                 / rs.rs_raw.seeded_shots))
                    counts.append((rs.rs_raw.molecular_counts[:, index]
                                 / rs.rs_raw.seeded_shots))
                    counts.append((rs.rs_raw.cross_pol_counts[:,index]
                               / rs.rs_raw.seeded_shots))
                    ylabel = 'Counts/bin/shot'
                legend.append('chi %3.1f km' %(req_altitude))
                legend.append('mol %3.1f km' %(req_altitude))
                legend.append('cpol %3.1f km' %(req_altitude))
                colors.extend(['r','g','b'])
                widths.extend([wid,wid,wid])
                wid=1
                                
       
        if len(counts)>0:         
            gt.plot_vs_time('counts_vs_time'
                 ,instrument
                 ,rs.rs_raw.times
                 ,counts
                 ,colors
                 ,widths
                 ,legend
                 ,'upper left'    
                 ,ylabel
                 ,None
                 ,title   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
            
    if display_defaults.enabled('linear_depolarization_vs_time') and hasattr(rs,'rs_inv')\
           and hasattr(rs.rs_inv,'linear_depol'):
        print 'rendering depolarization vs time'
        title= 'linear depol vs time'
        lines = []
        widths = []
        colors = []
        legends = []
        [lines,colors,widths,legends]=define_multiple_alt_plots(
            rs.rs_inv.linear_depol *100.0
            ,rs.rs_inv.msl_altitudes
            ,plot_alt_index
            ,lines
            ,widths
            ,colors
            ,legends)
      
        gt.plot_vs_time('linear_depolarization_vs_time'
                 ,instrument
                 ,rs.rs_inv.times
                 ,lines
                 ,colors     
                 ,widths     
                 ,legends
                 ,'upper left'    
                 ,'%'         
                 ,None
                 ,title   
                 ,auto_loop
                 ,display_defaults
                 ,figs)

    #plot sky counts after any time averaging is applied
    if display_defaults.enabled('sky_noise_counts') and hasattr(rs,'rs_mean'):
        nshots = 10000
        lines = []
        legend = []
        colors = []
        if rs_constants.has_key('comb_hi_detector_dark_count') and hasattr(rs.rs_mean,'c_hi_dark_counts'):
              chi = nshots * (rs.rs_mean.c_hi_dark_counts[:,0]/rs.rs_mean.seeded_shots\
                         - rs_constants['comb_hi_detector_dark_count'])
              mean_chi = nanmean(chi)
              lines.append(chi)
              legend.append('chi <%4.2f>' %(mean_chi))
              colors.append('r')
        if rs_constants.has_key('mol_detector_dark_count') and hasattr(rs.rs_mean,'mol_dark_counts'):  
               mol = nshots * (rs.rs_mean.mol_dark_counts[:,0]/rs.rs_mean.seeded_shots\
                          - rs_constants['mol_detector_dark_count'])
               mean_mol = nanmean(mol)
               lines.append(mol)
               legend.append('mol <%4.2f>' %(mean_mol))
               colors.append('b')      
        if rs_constants.has_key('cpol_detector_dark_count'):
             cpol = nshots * (rs.rs_mean.c_pol_dark_counts[:,0]/rs.rs_mean.seeded_shots\
                         - rs_constants['cpol_detector_dark_count'])
             mean_cpol = nanmean(cpol)
             lines.append(cpol)
             legend.append('cpol <%4.2f>' %(mean_cpol))
             colors.append('g')                                          
        if rs_constants.has_key('comb_lo_detector_dark_count') and hasattr(rs.rs_mean,'c_lo_dark_counts'):
             clo = nshots * (rs.rs_mean.c_lo_dark_counts[:,0]/rs.rs_mean.seeded_shots\
                         - rs_constants['comb_lo_detector_dark_count'])
             mean_clo = nanmean(clo)
             lines.append(clo)
             legend.append('clo <%4.2f>' %(mean_clo))
             colors.append('c')            
        if rs_constants.has_key('IR_detector_dark_count') and hasattr(rs.rs_mean,'combined_1064_dark_counts'):  
             c1064 = nshots * (rs.rs_mean.combined_1064_dark_counts[:,0]/rs.rs_mean.seeded_shots\
                          - rs_constants['IR_detector_dark_count'])
             mean_1064 = nanmean(c1064)
             lines.append(c1064)
             legend.append('clo <%4.2f>' %(mean_1064))
             colors.append('k')  
        if rs_constants.has_key('mol_i2a_detector_dark_count') and hasattr(rs.rs_mean,'mol_i2a_dark_counts'):  
             I2a_mol = nshots * (rs.rs_mean.mol_i2a_dark_counts[:,0]/rs.rs_mean.seeded_shots\
                          - rs_constants['mol_i2a_detector_dark_count'])                           
             mean_I2a = nanmean(I2a_mol)
             lines.append(I2a_mol)
             legend.append('I2a <%4.1f>' %(mean_I2a))
             colors.append('m')  
        if rs_constants.has_key('mol_wfov_detector_dark_count') and hasattr(rs.rs_mean,'m_wfov_dark_counts'):  
             wfov_mol = nshots * (rs.rs_mean.m_wfov_dark_counts[:,0]/rs.rs_mean.seeded_shots\
                          - rs_constants['mol_wfov_detector_dark_count'])                                  
             mean_wfov = nanmean(wfov_mol)
             lines.append(wfov_mol)
             legend.append('wfov <%4.1f>' %(mean_wfov))
             colors.append('m')  
        title = 'sky noise counts'

        print 'rendering sky noise counts'
        
        gt.plot_vs_time('sky_noise_counts'
                 ,instrument
                 ,rs.rs_mean.times
                 ,lines
                 ,colors
                 ,[1,1,1,1,1,1,1]
                 ,legend
                 ,'upper left'    
                 ,'sky noise counts per '+str(nshots)+ ' shots'
                 ,None
                 ,title   
                 ,auto_loop
                 ,display_defaults
                 ,figs)

    #plot dark counts after any time averaging is applied        
    if display_defaults.enabled('dark_counts') and hasattr(rs,'rs_mean'):   
        nshots = 10000
        mean_cpol = nshots * nanmean(rs.rs_mean.c_pol_dark_counts[:,0]/rs.rs_mean.seeded_shots)
        mean_chi = nshots * nanmean(rs.rs_mean.c_hi_dark_counts[:,0]/rs.rs_mean.seeded_shots)
        mean_mol = nshots * nanmean(rs.rs_mean.mol_dark_counts[:,0]/rs.rs_mean.seeded_shots)
        

        title = 'dark counts, chi=%4.1f,m=%4.1f,cp=%4.1f' \
                %(mean_chi,mean_mol,mean_cpol)

        lines = [nshots * rs.rs_mean.c_hi_dark_counts[:,0]/rs.rs_mean.seeded_shots]
        legend = ['chi']
        colors = ['r']
        
        if hasattr(rs.rs_mean,'c_lo_dark_counts'):
            mean_clo = nshots * nanmean(rs.rs_mean.c_lo_dark_counts[:,0]/rs.rs_mean.seeded_shots)
            title = title +',clo=%4.1f' %(mean_clo)
            lines.append(nshots * rs.rs_mean.c_lo_dark_counts[:,0]/rs.rs_mean.seeded_shots)
            legend.append('clo')
            colors.append('c')
            
        if hasattr(rs.rs_mean,'mol_i2a_dark_counts'):
            mean_i2a = nshots * nanmean(rs.rs_mean.mol_i2a_dark_counts[:,0]/rs.rs_mean.seeded_shots)
            title = title +',m_i2a=%4.1f' %(mean_i2a)
            lines.append(nshots * rs.rs_mean.mol_i2a_dark_counts[:,0]/rs.rs_mean.seeded_shots)
            legend.append(',mol_i2a')
            colors.append('m')
        if hasattr(rs.rs_mean,'combined_1064_dark_counts'):
            mean_ir = nshots * nanmean(rs.rs_mean.combined_1064_dark_counts[:,0]/rs.rs_mean.seeded_shots)
            title = title +',IR=%4.1f' %(mean_ir)
            lines.append(nshots * rs.rs_mean.combined_1064_dark_counts[:,0]/rs.rs_mean.seeded_shots)
            legend.append('1064')
            colors.append('k')
       
        lines.append(nshots * rs.rs_mean.mol_dark_counts[:,0]/rs.rs_mean.seeded_shots)
        legend.append('mol')
        colors.append('b')
       
        lines.append(nshots * rs.rs_mean.c_pol_dark_counts[:,0]/rs.rs_mean.seeded_shots)
        legend.append('cpol')
        colors.append('g')
        
        gt.plot_vs_time('dark_counts'
                 ,instrument
                 ,rs.rs_mean.times
                 ,lines
                 ,colors
                 ,[1,1,1,1,1,1]
                 ,legend
                 ,'upper left'    
                 ,'Dark counts per '+str(nshots)+ ' shots'
                 ,None
                 ,title   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
    #plot raw dark counts before any time averaging is applied        
    if display_defaults.enabled('raw_dark_counts') and hasattr(rs,'rs_raw'):   
        nshots = 10000
        mean_cpol = nshots * nanmean(rs.rs_raw.c_pol_dark_counts[:,0]/rs.rs_raw.seeded_shots)
        mean_chi = nshots * nanmean(rs.rs_raw.c_hi_dark_counts[:,0]/rs.rs_raw.seeded_shots)
        mean_mol = nshots * nanmean(rs.rs_raw.mol_dark_counts[:,0]/rs.rs_raw.seeded_shots)
        

        title = 'raw dark counts, chi=%4.1f,m=%4.1f,cp=%4.1f' \
                %(mean_chi,mean_mol,mean_cpol)

        lines = [nshots * rs.rs_raw.c_hi_dark_counts[:,0]/rs.rs_raw.seeded_shots[:]]
        legend = ['chi']
        colors = ['r']
        
        if hasattr(rs.rs_raw,'c_lo_dark_counts'):
            mean_clo = nshots * nanmean(rs.rs_raw.c_lo_dark_counts[:,0]/rs.rs_raw.seeded_shots)
            title = title +',clo=%4.1f' %(mean_clo)
            lines.append(nshots * rs.rs_raw.c_lo_dark_counts[:,0]/rs.rs_raw.seeded_shots)
            legend.append('clo')
            colors.append('c')
            
        if hasattr(rs.rs_raw,'mol_i2a_dark_counts'):
            mean_i2a = nshots * nanmean(rs.rs_raw.mol_i2a_dark_counts[:,0]/rs.rs_raw.seeded_shots)
            title = title +',m_i2a=%4.1f' %(mean_i2a)
            lines.append(nshots * rs.rs_raw.mol_i2a_dark_counts[:,0]/rs.rs_raw.seeded_shots)
            legend.append(',mol_i2a')
            colors.append('m')
        if hasattr(rs.rs_raw,'combined_1064_dark_counts'):
            mean_ir = nshots * nanmean(rs.rs_raw.combined_1064_dark_counts[:,0]/rs.rs_raw.seeded_shots)
            title = title +',IR=%4.1f' %(mean_ir)
            lines.append(nshots * rs.rs_raw.combined_1064_dark_counts[:,0]/rs.rs_raw.seeded_shots)
            legend.append('1064')
            colors.append('k')
       
        lines.append(nshots * rs.rs_raw.mol_dark_counts[:,0]/rs.rs_raw.seeded_shots)
        legend.append('mol')
        colors.append('b')
       
        lines.append(nshots * rs.rs_raw.c_pol_dark_counts[:,0]/rs.rs_raw.seeded_shots)
        legend.append('cpol')
        colors.append('g')
        
        gt.plot_vs_time('raw_dark_counts'
                 ,instrument
                 ,rs.rs_raw.times
                 ,lines
                 ,colors
                 ,[1,1,1,1,1,1]
                 ,legend
                 ,'upper left'    
                 ,'Dark counts per '+str(nshots)+ ' shots'
                 ,None
                 ,title   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
        
    if display_defaults.enabled('transmitted_energy'):
      fr=None
      pref=''
      for f,p in (('rs_raw','Raw '),('rs_mean','Mean ')):
        if hasattr(rs,f):
         fr=getattr(rs,f)
         pref=p
         break
      
      if fr is not None:
         lines = [fr.transmitted_power]
         legend =['532']

         if hasattr(fr,'transmitted_1064_power'):
              lines.append(fr.transmitted_1064_power)
              legend.append('1064')              
         gt.plot_vs_time('transmitted_energy'
                 ,instrument
                 ,fr.times
                 ,lines
                 ,['r','k']
                 ,[1,1]
                 ,legend
                 ,'upper left'    
                 ,pref+'Transmitted power'
                 ,'mW'
                 ,pref+'Transmitted power'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)

    if display_defaults.enabled('gv_qwp_rotation') \
            and hasattr(rs,'rs_mean') and hasattr(rs.rs_mean,'qwp_theta'):
         gt.plot_vs_time('gv_qwp_rotation'
                 ,instrument
                 ,rs.rs_mean.times
                 ,[rs.rs_mean.qwp_theta]
                 ,None
                 ,None
                 ,None
                 ,None
                 ,'quarter wave plate angle '
                 ,'deg'
                 ,'Quarter wave plate angle'
                 ,auto_loop
                 ,display_defaults
                 ,figs)


    
  
    if display_defaults.enabled('short_cell_ratio') and 'shortcell_locked_ratio' in rs_constants:
      
      fr=None
      if fr is None and hasattr(rs,'rs_mean'):
        fr = getattr(rs,'rs_mean')
      if fr is None and hasattr(rs,'rs_raw'):
        fr = getattr(rs,'rs_raw')
    
      if fr is not None and hasattr(fr,'filtered_energy'):
        #older system did not include min,mean, and max values in filtered and unfilted energies
        #always plot raw ratios for these systems
    
        if display_defaults.get_value('short_cell_ratio','type') == 'raw' or fr.filtered_energy.shape[1] < 2:
             shortcell_ratio =  fr.filtered_energy[:,0] / fr.nonfiltered_energy[:,0]
             if fr.filtered_energy.shape[1] > 1:
                 max_shortcell_ratio = fr.filtered_energy[:,2] / fr.nonfiltered_energy[:,0]
                 #if not rs_constants['shortcell_locked_ratio'][0] == -9999 :
                 lock_pt = np.zeros(len(fr.times))          
                 warn_pt = lock_pt + rs_constants['shortcell_locked_ratio'][0]
                 lost_pt = lock_pt + rs_constants['shortcell_locked_ratio'][1]
             else:
                 max_shortcell_ratio = np.NaN *np.zeros(len(fr.times))
                 warn_pt = np.NaN *np.zeros(len(fr.times))
                 lost_pt = np.NaN *np.zeros(len(fr.times))
             title = 'raw short cell ratio'
             legend = ['mean_ratio']
        else:   #plot shotcell_ratio scaled by locked and unlocked values of the raw ratio

            #the following is the ratio of the mean shortcell energies
            shortcell_ratio = (fr.filtered_energy[:,0] / fr.nonfiltered_energy[:,0] \
                           -rs_constants['shortcell_locked_ratio'][0])\
               / (rs_constants['shortcell_locked_ratio'][1]-rs_constants['shortcell_locked_ratio'][0])

            #this is the ratio mean of the max filtered energies
            if len(fr.filtered_energy[0,:])>1:
               max_shortcell_ratio = (fr.filtered_energy[:,2] / fr.nonfiltered_energy[:,0] \
                               -rs_constants['shortcell_locked_ratio'][0])\
               / (rs_constants['shortcell_locked_ratio'][1]-rs_constants['shortcell_locked_ratio'][0])
            else:
                max_shortcell_ratio = np.NaN * np.zeros(len(fr.times))
            title = 'scaled short cell ratio'
            legend = ['mean_ratio','max_ratio','warning','lock lost']
            if not rs_constants['shortcell_locked_ratio'][0] == -9999 :
                lock_pt = np.zeros(len(fr.times))          
                warn_pt = lock_pt + processing_defaults.get_value('I2_lock_mask','lock_warning')
                lost_pt = lock_pt + processing_defaults.get_value('I2_lock_mask','lock_lost')
            else:
                 max_shortcell_ratio = np.NaN *np.zeros(len(fr.times))
                 warn_pt = np.NaN *np.zeros(len(fr.times))
                 lost_pt = np.NaN *np.zeros(len(fr.times))
        
          
        gt.plot_vs_time('short_cell_ratio'
                ,instrument
                 ,fr.times
                 ,[shortcell_ratio,max_shortcell_ratio,warn_pt,lost_pt]
                 ,['r','c','k','g','b']
                 ,[2,1,1,1]
                 ,legend
                 ,'upper left'    
                 ,'Short cell I2 ratio '
                 ,None
                 ,title   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
      else:
            print 'no short cell ratio plot--short cell lock ratio not provided in calvals.xxxx.txt'
  
    if display_defaults.enabled('altitude') and hasattr(rs,'rs_mean') and \
       'installation' in rs_constants  and hasattr(rs.rs_mean,'telescope_pointing') and  rs_constants['installation']=='airborne' :
         mask=np.copy(rs.rs_mean.telescope_pointing.astype('float64'))
         mask[mask<.999]=np.NaN
         gt.plot_vs_time('altitude'
                 ,instrument
                 ,rs.rs_mean.times
                 ,[rs.rs_mean.GPS_MSL_Alt / 1000.
                   ,rs.rs_mean.GPS_MSL_Alt * mask / 1000]
                 ,['k','r']
                 ,None
                 ,['tel-down','tel-up']
                 ,'lower left'    
                 ,'MSL Altitude'
                 ,'km'
                 ,'Altitude'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
   

    if display_defaults.enabled('pitch_roll_angles')  and hasattr(rs,'rs_raw') \
       and 'installation' in rs_constants and hasattr(rs.rs_raw,'telescope_pointing') and  rs_constants['installation']=='airborne':
        mask = np.copy(rs.rs_raw.telescope_pointing).astype('float64')
        mask[mask<.999] = np.NaN
        gt.plot_vs_time('pitch_roll_angles'
                 ,instrument
                 ,rs.rs_raw.times
                 ,[rs.rs_raw.pitch_angle
                 ,rs.rs_raw.pitch_angle * mask
                 ,rs.rs_raw.roll_angle
                 , rs.rs_raw.roll_angle * mask]
                 ,['b','r','k','g']
                 ,[1,3,1,3]
                 ,['pitch dn','pitch up','roll dn','roll up']
                 ,'lower left'    
                 ,'Pitch, Roll '
                 ,'deg'
                 ,'Pitch, Roll, Tel pointing'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)

    if 'installation' in rs_constants and (hasattr(rs,"rs_mean") or hasattr(rs,"rs_raw")) \
          and (rs_constants['installation']=='airborne' \
               or rs_constants['installation'] == 'shipborne') \
                    and display_defaults.enabled('lat_long'):

           
        if hasattr(rs,"rs_mean") and not (hasattr(rs.rs_mean,'latitude')) \
               or hasattr(rs,"rs_raw") and not (hasattr(rs.rs_raw,'latitude')):
            print
            print 'WARNING:*****skipping lat_long plot -- no lat_long data'
            return
        
    
        else:
            #add time and altitude text along lat-long line
            angle=[]
            text_position_x=[]
            text_position_y=[]
            text_str=[]
            label_indices = []
            if hasattr(rs,'rs_mean'):
                times = rs.rs_mean.times.copy()
                latitude = rs.rs_mean.latitude.copy()
                longitude = rs.rs_mean.longitude.copy()
                nshots = rs.rs_mean.latitude.shape[0]
                if rs_constants['installation'] == 'airborne':
                   telescope_pointing = rs.rs_mean.telescope_pointing
                   GPS_MSL_Alt = rs.rs_mean.GPS_MSL_Alt
            else:
                times = rs.rs_raw.times.copy()
                latitude = rs.rs_raw.latitude.copy()
                longitude = rs.rs_raw.longitude.copy()
                nshots = rs.rs_raw.latitude.shape[0]
                if rs_constants['installation'] == 'airborne':
                    telescope_pointing = rs.rs_raw.telescope_pointing
                    GPS_MSL_Alt = rs.rs_raw.GPS_MSL_Alt
        
            dd = display_defaults.get_value("lat_long","label_interval_days")
            if dd == None:
                dd = 0  
            hh = display_defaults.get_value("lat_long","label_interval_hours")
            if hh == None:
                hh = 0
            mm = display_defaults.get_value("lat_long","label_interval_minutes")
            if mm == None:
                mm = 0
                
            if dd + hh + mm > 0:      
            
                dt = datetime.timedelta(days = dd, hours = hh, minutes = mm)
                if (times[-1] - times[0]) < 3 * dt:
                    print
                    print
                    print
                    print 'requested label time spacing too small in display defaults'
                    print 'no lat-lon map plotted'
                    print
                    print
                    return
                start_date = times[0]
        
                for i in range(nshots-10):
                    if times[i]> start_date:
                        start_date = start_date + dt
                        if not i < 10 and not i > nshots-10:
                           label_indices.append(i)
                                  
            else:
                 n_locs = display_defaults.get_value('lat_long','n plot points')
                 if np.int(len(latitude)/n_locs) ==0: 
                     print
                     print 'WARNING:***skipping lat_long plot - not enough data (is vehicle parked?)'
                     return
                 nshots = latitude.shape[0]
                 step = np.int(nshots/n_locs) 
                 label_indices = range(0,nshots-2-n_locs,np.int(nshots/n_locs))
            if rs_constants['installation'] == 'shipborne':
                t_up = np.ones_like(latitude)
            else:
                #airborne
                t_up = np.copy(telescope_pointing) * 1.0
                t_up[t_up != 1.0] = np.nan
                longitude_t_up = longitude * t_up
                latitude_t_up = latitude * t_up

                    
            for i in range(len(label_indices)):
                if not np.isnan(latitude[label_indices[i]]) and not np.isnan(longitude[label_indices[i]]):
                     # compute slope of lat/lon lin
                     dy = latitude[label_indices[i]+9] \
                              - latitude[label_indices[i]-9]
                     dx = longitude[label_indices[i]+9] \
                                 - longitude[label_indices[i]-9]

                     # avoid ValueError caused by trying to convert Nan to integer

                     if np.isnan(dy) or np.isnan(dx):
                         continue
                     if dx == 0:
                         ang = 0.0
                     else:
                         ang = 180 / np.pi * np.arctan(dy / dx) + 90
                     if ang > 90:
                         ang = ang - 180
                     angle.append(ang)
                     if rs_constants['installation'] == 'shipborne':
                          alt_km = rs_constants['lidar_altitude']/1000.0
                     else:
                          #airborne
                          alt_km = GPS_MSL_Alt[label_indices[i]] / 1000.
                     text_position_x.append(longitude[label_indices[i]])
                     text_position_y.append(latitude[label_indices[i]])
                     if rs_constants['installation'] == 'airborne':
                          text_str.append(times[label_indices[i]].strftime('%H:%M ') + '    ' \
                          + '%3.1fkm' % alt_km)
                     else:
                         print 'dd = ',dd,'  hh= ',hh
                         if dd > 0 or hh > 4 :

                             print 'dd > 0 or hh >4'
                             text_str.append(times[label_indices[i]].strftime('%m/%d ') + '    ' \
                                  + times[label_indices[i]].strftime('%H:%M'))
                         else:    
                             text_str.append(times[label_indices[i]].strftime('%H  %M '))
          
            #if not np.isnan(longitude_t_up).any():
            if rs_constants['installation'] == 'airborne':
                longi = [longitude
                    ,longitude_t_up,longitude[label_indices]]
                lat = [latitude
                    ,latitude_t_up,latitude[label_indices]]
                colors = ['k','r','k']
            else:
                longi = [longitude,longitude,longitude[label_indices]]
                lat = [latitude,latitude,latitude[label_indices]]
                colors = ['r','r','r']
            
    
            ftest=figs.figure('lat_long')
            doGraph=True
            if ftest.isNewFigure:
                figs.close(ftest)#was new, keep it new
            elif hasattr(rs,'rs_raw'):
                if ftest.sourceframe!='rs_raw':
                  figs.close(ftest)#was made, but we have raw now, recreate it
            elif ftest.sourceframe!='rs_mean':
                doGraph=False
            if doGraph:
              print 'rendering lat-lon plot'  
              gt.plot_map_latlon('lat_long'
                ,instrument
                ,times
                ,longi
                ,lat 
                ,colors   
                ,['None','None','o']
                ,[0,0,8]
                ,['-','-','None']
                ,[2,2,2]
                ,'Location'
                ,text_str
                ,text_position_x
                ,text_position_y
                ,angle
                ,display_defaults
                ,figs)
              
              setattr(figs.figure('lat_long'),'sourceframe','rs_raw' if hasattr(rs,'rs_raw') else 'rs_mean')
        
       
    if display_defaults.enabled('beam_position') and hasattr(rs,'rs_raw') \
             and hasattr(rs.rs_raw,'cg_xs')\
             and hasattr(rs.rs_raw,'cg_ys'):
        if hasattr(rs.rs_raw,'cg_xs2'):
            lines =[rs.rs_raw.cg_xs, rs.rs_raw.cg_ys
                    ,rs.rs_raw.cg_xs2,(rs.rs_raw.cg_ys2-250)*3.0 +220]
            legend = ['x','y','x2','y2']
        else:
            lines=[rs.rs_raw.cg_xs, rs.rs_raw.cg_ys]
            legend = ['x','y']
        gt.plot_vs_time('beam_position'
                 ,instrument
                 ,rs.rs_raw.times
                 ,lines
                 ,['b','r','g','k']
                 ,[3,3,1,1]
                 ,legend
                 ,'lower left'    
                 ,'Beam position '
                 ,'pixels'
                 ,'Beam position'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
            
    if display_defaults.enabled('ktp_temperature') \
            and hasattr(rs,'rs_raw') and hasattr(rs.rs_raw,'ktp_temp'):
        gt.plot_vs_time('ktp_temperature'
                 ,instrument
                 ,rs.rs_raw.times
                 ,[rs.rs_raw.ktp_temp,rs.rs_raw.ktp_temp_setpoint]
                 ,['b','r']
                 ,[2,2]
                 ,['temp','setpt']
                 ,'lower left'    
                 ,'KTP temp '
                 ,'deg_C'
                 ,'SHG temperature'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
        
    if display_defaults.enabled('l3_piezo_voltage') \
            and hasattr(rs,'rs_raw') and hasattr(rs.rs_raw,'l3cavityvoltage'):
        rs.rs_raw.piezo_voltage_ave[rs.rs_raw.piezo_voltage_ave >100] =np.NaN
        rs.rs_raw.piezo_voltage_min[rs.rs_raw.piezo_voltage_min >100] =np.NaN
        rs.rs_raw.piezo_voltage_max[rs.rs_raw.piezo_voltage_max >100] =np.NaN
        gt.plot_vs_time('l3cavityvoltage'
                 ,instrument
                 ,rs.rs_raw.times
                 ,[rs.rs_raw.piezo_voltage_ave
                    ,rs.rs_raw.piezo_voltage_min
                    ,rs.rs_raw.piezo_voltage_max]
                 ,['r','k','c']
                 ,[3,1,1]
                 ,['ave','min','max']
                 ,'lower left'    
                 ,'Cavity piezo'
                 ,'Volts'
                 ,'Laser cavity piezo'   
                 ,auto_loop
                 ,display_defaults
                 ,figs) 
    if display_defaults.enabled('l3_lockslope') \
            and hasattr(rs,'rs_raw') and hasattr(rs.rs_raw,'l3locking_stats'):
        if hasattr(rs.rs_raw,'l3frequency_offset'):
          metr=rs.rs_raw.l3frequency_offset
          units='Hz'
          label = 'L3 Frequency Offset'
          metr[metr >1e15] =np.NaN
        else:
          metr=rs.rs_raw.l3locking_stats.copy()
          units='slope'
          label = 'L3 Lock Slope Metric'
          metr[metr >100000] =np.NaN
        gt.plot_vs_time('l3_lockslope'
                 ,instrument
                 ,rs.rs_raw.times
                 ,[metr[:,0]
                    ,metr[:,1]
                    ,metr[:,2]]
                 ,['r','k','c']
                 ,[3,1,1]
                 ,['ave','min','max']
                 ,'lower left'    
                 ,'L3 Lock Slope'
                 ,units
                 ,label 
                 ,auto_loop
                 ,display_defaults
                 ,figs) 
    if display_defaults.enabled('chiller_temperatures') \
                and hasattr(rs,'rs_raw') and hasattr(rs.rs_raw,'chiller_temp'):
        gt.plot_vs_time('chiller_temperatures'
                 ,instrument
                 ,rs.rs_raw.times
                 ,[rs.rs_raw.chiller_temp, rs.rs_raw.chiller_setpt]
                 ,['b','r']
                 ,[2,1]
                 ,['temp','setpt']
                 ,'lower left'    
                 ,'Temperature '
                 ,'deg-C'
                 ,'Chiller temperature'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
    
    if display_defaults.enabled('laser_diode_temp') \
         and hasattr(rs,'rs_raw') and hasattr(rs.rs_raw,'laser_diode_temp'):
    
        
         gt.plot_vs_time('laser_diode_temp'
                 ,instrument
                 ,rs.rs_raw.times
                 ,[rs.rs_raw.laser_diode_temp, rs.rs_raw.laser_diode_temp_setpoint]
                 ,['b','r']
                 ,[2,1]
                 ,['temp','setpt']
                 ,'lower left'    
                 ,'Temperature '
                 ,'deg-C'
                 ,'laser diode temperature'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
         
    if display_defaults.enabled('laser_current') \
         and hasattr(rs,'rs_raw') and hasattr(rs.rs_raw,'laser_current'):
         lines = [rs.rs_raw.laser_current]
         legend = ['current']
         if hasattr(rs.rs_raw,'laser_current_setpoint'):
            lines.append(rs.rs_raw.laser_current_setpoint)
            legend.append('setpt')
            
         gt.plot_vs_time('laser_current'
                 ,instrument
                 ,rs.rs_raw.times
                 ,lines
                 ,['b','r']
                 ,[2,1]
                 ,legend
                 ,'lower left'    
                 ,'Current '
                 ,'A'
                 ,'Laser diode current'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
    if display_defaults.enabled('laser_voltage') \
         and hasattr(rs,'rs_raw') and hasattr(rs.rs_raw,'laser_voltage'):           
         
         gt.plot_vs_time('laser_voltage'
                 ,instrument
                 ,rs.rs_raw.times
                 ,[rs.rs_raw.laser_voltage]
                 ,['b']
                 ,[2]
                 ,None
                 ,None  
                 ,'Voltage '
                 ,'V'
                 ,'Laser diode fwd voltage'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
         
    if display_defaults.enabled('seed_percent') \
          and hasattr(rs,'rs_raw')   and hasattr(rs.rs_raw,'seeded_shots'):
        gt.plot_vs_time('seed_percent'
                 ,instrument
                 ,rs.rs_raw.times
                 ,[100.0*rs.rs_raw.seeded_shots/rs.rs_raw.shot_count]
                 ,['b']
                 ,[2]
                 ,None
                 ,None  
                 ,'Seed percent '
                 ,'%'
                 ,'seeding percentage'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
   
#fixme these should use a delta time value, with 0 as start time, not epoch
    
    if display_defaults.enabled('tcomp_interf_freq')\
           and hasattr(rs,'rs_raw')  and hasattr(rs.rs_raw,'tcomp_interf_freq'):
        #pc = np.polyfit(date2num(rs.rs_raw.times), rs.rs_raw.interf_freq,
        #        display_defaults.get_value('interferometer_freq_deviations',
        #        'polyfit degree'))
        #freq_fit = np.polyval(pc, date2num(rs.rs_raw.times))                
        gt.plot_vs_time('tc_comp_interf_freq'
                 ,instrument
                 ,rs.rs_raw.times
                 ,[rs.rs_raw.interf_freq/1e6,rs.rs_raw.tcomp_interf_freq/1e6]
                 ,['c','r']
                 ,[1,2]
                 ,None
                 ,['raw','T comp']    
                 ,'Frequency'
                 ,'MHz'
                 ,'Temp comp Laser frequency'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
   

    if display_defaults.enabled('interferometer_freq')\
           and hasattr(rs,'rs_raw')  and hasattr(rs.rs_raw,'interf_freq'):
        #pc = np.polyfit(date2num(rs.rs_raw.times), rs.rs_raw.interf_freq,
        #        display_defaults.get_value('interferometer_freq_deviations',
        #        'polyfit degree'))
        #freq_fit = np.polyval(pc, date2num(rs.rs_raw.times))     
        lines = [(rs.rs_raw.interf_freq - rs.rs_raw.interf_freq[0])/1e6]
        legend= ['interf']
        if hasattr(rs.rs_raw,'superseedlasercontrollog'):
            if rs.rs_raw.superseedlasercontrollog.shape[1]==11:
              index=7
            elif rs.rs_raw.superseedlasercontrollog.shape[1]==10:
              index=2
            else:
              raise NotImplemented("Can't get frequency from superseed board. %i element not supported" % (rs.rs_raw.superseedlasercontrollog.shape[1]))
            seed_freq = (rs.rs_raw.superseedlasercontrollog[:,index] \
                 - rs.rs_raw.superseedlasercontrollog[0,index])
            if 'seedlaser_temp_to_freq' in rs_constants:
                seed_freq *= rs_constants['seedlaser_temp_to_freq']*1e3
                legend.append('calibratied seeder')
            else:
                seed_freq *= lines[0][-1]/seed_freq[-1]
                legend.append('normalized seeder')
            lines.append(seed_freq)
        gt.plot_vs_time('interferometer_freq'
                 ,instrument
                 ,rs.rs_raw.times
                 ,lines
                 ,['b','c']
                 ,[2,2]
                 ,legend
                 ,'upper right'       
                 ,'Frequency'
                 ,'MHz'
                 ,'Laser frequency deviation'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
        
    if display_defaults.enabled('interferometer_temp')\
           and hasattr(rs,'rs_raw')  and hasattr(rs.rs_raw,'interferometer_temp'):     
        gt.plot_vs_time('interferometer_temp'
                 ,instrument
                 ,rs.rs_raw.times
                 ,[rs.rs_raw.interferometer_temp]
                 ,['b']
                 ,[2]
                 ,None
                 ,'upper right'       
                 ,'temperature'
                 ,'deg C'
                 ,'Interferometer temperature'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
   
    if display_defaults.enabled('superseed_controller_temps')\
           and hasattr(rs,'rs_raw') and hasattr(rs.rs_raw,'superseedlasercontrollog'):
         if rs.rs_raw.superseedlasercontrollog.shape[1]==11:
          indexes=(9,7,1)
          tempnames=['ambient','crystal','PCB']
         elif rs.rs_raw.superseedlasercontrollog.shape[1]==10:
          indexes=(1,2)
          tempnames=['ambient','crystal']
         else:
          raise NotImplemented("Can't get frequency from superseed board. %i element not supported" % (rs.rs_raw.superseedlasercontrollog.shape[1]))

         gt.plot_vs_time('superseed_controller_temps'
                 ,instrument
                 ,rs.rs_raw.times
                 ,[rs.rs_raw.superseedlasercontrollog[:,x] for x in indexes]
                 ,['b','g','r']
                 ,[2,2,2]
                 ,tempnames
                 ,'upper right'        
                 ,'Super seed laser temperatures'
                 ,'deg-C'
                 ,'Super seed laser temperatures'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
    if display_defaults.enabled('superseed_controller_voltages')\
           and hasattr(rs,'rs_raw') and hasattr(rs.rs_raw,'superseedlasercontrollog'):
         if rs.rs_raw.superseedlasercontrollog.shape[1]==11:
            index=6
         elif rs.rs_raw.superseedlasercontrollog.shape[1]==10:
            index=0
         else:
          raise NotImplemented("Can't get voltage from superseed board. %i element not supported" % (rs.rs_raw.superseedlasercontrollog.shape[1]))
         gt.plot_vs_time('superseed_controller_voltages'
                 ,instrument
                 ,rs.rs_raw.times
                 ,[rs.rs_raw.superseedlasercontrollog[:,index]]
                 ,['b','g']
                 ,[2,2]
                 ,None
                 ,None        
                 ,'Control voltage'
                 ,'Volts'
                 ,'Super seed laser voltage'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
    if display_defaults.enabled('superseed_peltier_power')\
           and hasattr(rs,'rs_raw') and hasattr(rs.rs_raw,'superseedlasercontrollog'):
         if rs.rs_raw.superseedlasercontrollog.shape[1]==11:
            indexes=(3,4)
            labels=['photo','diode']
            units='Amps'
         elif rs.rs_raw.superseedlasercontrollog.shape[1]==10:
            indexes=(4,5)
            labels=['crystal','diode']
            units='Watts'
         else:
          raise NotImplemented("Can't get voltage from superseed board. %i element not supported" % (rs.rs_raw.superseedlasercontrollog.shape[1]))
         gt.plot_vs_time('superseed_peltier_power'
                 ,instrument
                 ,rs.rs_raw.times
                 ,[rs.rs_raw.superseedlasercontrollog[:,x] for x in indexes]
                 ,['b','g']
                 ,[2,2]
                 ,labels
                 ,'upper right'       
                 ,'Power to Peltier'
                 ,units
                 ,'Super seed peltier power'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)   
    if display_defaults.enabled('seed_voltage')\
           and hasattr(rs,'rs_raw')  and hasattr(rs.rs_raw,'seedvoltage'):
         gt.plot_vs_time('seed_voltage'
                 ,instrument
                 ,rs.rs_raw.times
                 ,[rs.rs_raw.seedvoltage]
                 ,['b']
                 ,[2,]
                 ,None
                 ,None  
                 ,'Seed voltage'
                 ,'V'
                 ,'Seed Voltage'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
   
  
    if display_defaults.enabled('cal_pulse_time_averaged') and hasattr(rs,'rs_mean'):
      legendList=[]
      y_vars=[]
      colors=[]
      energy_norm = rs.rs_mean.transmitted_energy/nanmean(rs.rs_mean.transmitted_energy)
     
      if hasattr(rs.rs_mean,"transmitted_1064_energy"):
        energy_1064_norm = rs.rs_mean.transmitted_1064_energy/nanmean(rs.rs_mean.transmitted_1064_energy)
      else:
        energy_1064_norm = energy_norm
      if hasattr(rs.rs_mean,'molecular_cal_pulse')\
                   and display_defaults.get_value('cal_pulse_time_averaged','mol/chi/clo/m_i2a/wfov/IR',missing='mol').find('mol') >= 0:
            y_vars.append(rs.rs_mean.molecular_cal_pulse/energy_norm)
            colors.append('b')
            legendList.append('mol')
            
      if hasattr(rs.rs_mean,'combined_lo_cal_pulse')\
                   and display_defaults.get_value('cal_pulse_time_averaged','mol/chi/clo/m_i2a/wfov/IR',missing='clo').find('clo') >= 0:
            y_vars.append(rs.rs_mean.combined_lo_cal_pulse/energy_norm)
            colors.append('c')
            legendList.append('comb lo')
           
      if hasattr(rs.rs_mean,'combined_hi_cal_pulse')\
                  and display_defaults.get_value('cal_pulse_time_averaged','mol/chi/clo/m_i2a/wfov/IR',missing='chi').find('chi') >= 0:
            y_vars.append(rs.rs_mean.combined_hi_cal_pulse/energy_norm)
            colors.append('r')
            legendList.append('comb hi')
      if hasattr(rs.rs_mean,'molecular_i2a_cal_pulse')\
                   and display_defaults.get_value('cal_pulse_time_averaged','mol/chi/clo/m_i2a/wfov/IR',missing='m_i2a').find('m_i2a') >= 0:
                y_vars.append(rs.rs_mean.molecular_i2a_cal_pulse/energy_norm)
                colors.append('m')
                legendList.append('mol i2a')
      if hasattr(rs.rs_mean,'molecular_wfov_cal_pulse')\
                   and display_defaults.get_value(
                   'cal_pulse_time_averaged','mol/chi/clo/m_i2a/wfov/IR',missing='wfov').find('wfov') >= 0:
                y_vars.append(rs.rs_mean.molecular_wfov_cal_pulse/energy_norm)
                colors.append('k')
                legendList.append('mol wfov')
      if hasattr(rs.rs_mean,'combined_1064_cal_pulse')\
                   and display_defaults.get_value(
                   'cal_pulse_time_averaged','mol/chi/clo/m_i2a/wfov/IR',missing='IR').find('IR') >= 0:
                y_vars.append(rs.rs_mean.combined_1064_cal_pulse/energy_1064_norm)
                colors.append('k')
                legendList.append('comb 1064')              

      if hasattr(rs.rs_mean,'molecular_cal_pulse') and hasattr(rs.rs_mean,'combined_hi_cal_pulse'):
            i2_trans=rs.rs_mean.molecular_cal_pulse/rs.rs_mean.combined_hi_cal_pulse
            y_vars.append(100.0*i2_trans/nanmax(i2_trans))
            colors.append('g')                
            legendList.append('I2 trans')
      
      gt.plot_vs_time('cal_pulse_time_averaged'
                 ,instrument
                 ,rs.rs_mean.times
                 ,y_vars
                 ,colors
                 ,None
                 ,legendList
                 ,'lower left'    
                 ,'Counts/shot, trans'
                 ,'counts,%'
                 ,'scattered light pulse'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
      
    if display_defaults.enabled('cal_pulse_raw') and hasattr(rs,'rs_raw'):
      legendList=[]
      y_vars=[]
      colors=[]
      energy_norm = rs.rs_raw.transmitted_energy/nanmean(rs.rs_raw.transmitted_energy)

      if hasattr(rs.rs_raw,'transmitted_1064_energy'):
        energy_1064_norm = rs.rs_raw.transmitted_1064_energy/nanmean(rs.rs_raw.transmitted_1064_energy)
      else:
        energy_1064_norm = energy_norm

      if hasattr(rs.rs_raw,'molecular_cal_pulse')\
                   and display_defaults.get_value('cal_pulse_raw','mol/chi/clo/m_i2a/wfov/IR',missing='mol').find('mol') >= 0:
            y_vars.append(rs.rs_raw.molecular_cal_pulse/energy_norm)
            colors.append('b')
            legendList.append('mol')
            
      if hasattr(rs.rs_raw,'combined_lo_cal_pulse')\
                   and display_defaults.get_value('cal_pulse_raw','mol/chi/clo/m_i2a/wfov/IR',missing='clo').find('clo') >= 0:
            y_vars.append(rs.rs_raw.combined_lo_cal_pulse/energy_norm)
            colors.append('c')
            legendList.append('comb lo')
      
      if hasattr(rs.rs_raw,'combined_hi_cal_pulse')\
               and display_defaults.get_value(
               'cal_pulse_raw','mol/chi/clo/m_i2a/wfov/IR',missing='chi').find('chi') >= 0:
            y_vars.append(rs.rs_raw.combined_hi_cal_pulse/energy_norm)
            colors.append('r')
            legendList.append('comb hi')

      if hasattr(rs.rs_raw,'molecular_i2a_cal_pulse')\
                   and display_defaults.get_value('cal_pulse_raw','mol/chi/clo/m_i2a/wfov/IR',missing='m_i2a').find('m_i2a') >= 0:
                y_vars.append(rs.rs_raw.molecular_i2a_cal_pulse/energy_norm)
                colors.append('m')
                legendList.append('mol i2a')
      if hasattr(rs.rs_raw,'molecular_wfov_cal_pulse')\
                   and display_defaults.get_value(
                   'cal_pulse_raw','mol/chi/clo/m_i2a/wfov/IR',missing='wfov').find('wfov') >= 0:
                y_vars.append(rs.rs_raw.molecular_wfov_cal_pulse/energy_norm)
                colors.append('k')
                legendList.append('mol wfov')       
      if hasattr(rs.rs_raw,'molecular_cal_pulse') and hasattr(rs.rs_raw,'combined_hi_cal_pulse'):
            i2_trans=rs.rs_raw.molecular_cal_pulse/rs.rs_raw.combined_hi_cal_pulse
            y_vars.append(100.0*i2_trans/nanmax(i2_trans))
            colors.append('g')                
            legendList.append('I2 trans')
      if hasattr(rs.rs_raw,'combined_1064_cal_pulse') \
                 and display_defaults.get_value(
                 'cal_pulse_raw','mol/chi/clo/m_i2a/wfov/IR',missing='IR').find('IR') >= 0 :
            y_vars.append(rs.rs_raw.combined_1064_cal_pulse/energy_1064_norm)
            colors.append('k')                
            legendList.append('comb 1064')
            
      gt.plot_vs_time('cal_pulse_raw'
                 ,instrument
                 ,rs.rs_raw.times
                 ,y_vars
                 ,colors
                 ,None
                 ,legendList
                 ,'lower left'    
                 ,'Counts/shot, trans'
                 ,'counts,%'
                 ,'raw scattered light pulse'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
            
   
       
    if display_defaults.enabled('i2_spectrum'):
      
      fr=None
      if fr is None and hasattr(rs,'rs_mean'):
        fr=rs.rs_mean
      if fr is None and hasattr(rs,'rs_raw'):
        fr=rs.rs_raw
      if fr is not None and hasattr(fr,'interf_freq'):
              if display_defaults.get_value('i2_spectrum','show_theory'):
                    line_type = display_defaults.get_value('i2_spectrum','show_theory')
              else:
                    line_type = 'o'


              #plot with x_y_plot and discrete symbols
              #if fr!=None and line_type == 'o':
              if line_type == 'o':
                mol_temp =fr.molecular_cal_pulse.copy()
                var_x=[]
                var_y=[]
                mol_temp[mol_temp<= 0] = 1e-8
                mol_temp = mol_temp / fr.combined_hi_cal_pulse
                mol_temp[mol_temp == np.inf]= 1e-8
                #normalize to form i2 transmission
                freq_offset = np.float(display_defaults.get_value('i2_spectrum','freq_offset'))
                var_x.append(fr.interf_freq/1e9 + freq_offset)
                var_y.append(mol_temp/np.max(mol_temp))
                legend=['measured']
                if display_defaults.get_value('i2_spectrum','show_theory'):
                    filename = locate_file('I2cell_272_31_extended_freq.txt',systemOnly=True)
                    print 'reading file ',filename
                    [jnk, theory] = hru.readascii(filename )
                    var_x.append(theory[:,0])
                    var_y.append(theory[:,2])
                    legend.append('theory')
                if hasattr(fr,'molecular_i2a_cal_pulse'):
                    i2a_temp=fr.molecular_i2a_cal_pulse
                    i2a_temp=i2a_temp/fr.combined_hi_cal_pulse
                    i2a_temp[i2a_temp<=0] = 1e-8
                    i2a_temp[i2a_temp == np.inf] = 1e-8
                    var_x.append(fr.interf_freq/1e9 + freq_offset)
                    var_y.append(i2a_temp/np.max(i2a_temp))
                    legend.append('mol_I2A')
                    if 0:
                         import os
                         filename = 'i2_scan.txt'
                         fileid=open(filename,'w')
                         os.chmod(filename,0664)
                         print >>fileid, '#profile data from ',fr.times[0] \
                              ,' --> ', fr.times[-1], ' UTC'
                         print >>fileid, '#freq(GHz)     mol         mol_i2'
                         for i in range(np.size(var_y[0])):
                             print >>fileid, '%10.4e   %10.3e  %10.3e'\
                                   %(np.float(var_x[0][i]),np.float(var_y[0][i])
                                          ,np.float(var_y[2][i]))    
                         fileid.close()
              
                gt.plot_xy('i2_spectrum'
                        ,instrument
                        ,fr.times
                        ,var_x
                        ,var_y
                        ,['r','k','g']
                        ,['o','None','o']
                        ,[8,0,8,0]
                        ,['None','-','None']
                        ,[2,2,2]
                        ,legend
                        ,'lower left'
                        ,'Frequency'
                        ,'GHz'
                        ,'I2 transmission'
                        ,None
                        ,'I2 spectrum'
                        ,[]
                        ,[]
                        ,[]
                        ,[]
                        ,display_defaults
                        ,figs)


              #elif fr!=None:
              else:
                  trans = fr.molecular_cal_pulse.copy()
                  trans[trans <= 0] = 1e-8
                  #normalize to for i2 transmission
                  trans = trans / fr.combined_hi_cal_pulse
                  trans[trans == np.inf]= 1e-8
                  if display_defaults.get_value('i2_spectrum','trans_scaling'):
                     gain = np.float(display_defaults.get_value('i2_spectrum','trans_scaling'))
                     trans = trans / gain
                  
                  if display_defaults.get_value('i2_spectrum','freq_scaling'):
                      freq_gain = np.float(display_defaults.get_value('i2_spectrum','freq_scaling'))
                      freq_offset = np.float(display_defaults.get_value('i2_spectrum','freq_offset'))
                  else:
                      freq_gain = 1.0
                      freq_offset = 0.0
                  #rescale frequency    
                  freq = fr.interf_freq * freq_gain/1e9 + freq_offset
                
                  legend=['measured']
                  if display_defaults.get_value('i2_spectrum','show_theory'):
                      filename = locate_file('I2cell_272_31_extended_freq.txt',systemOnly=True)
                      print 'reading file ',filename
                      [jnk, theory] = hru.readascii(filename )
                      model = theory[:,2]
                      model_freq = theory[:,0]
                      legend.append('theory')
                  
                  if hasattr(fr,'filtered_energy'):
                      short_cell_ratio = fr.filtered_energy[:,0] \
                         / fr.nonfiltered_energy[:,0]
                      legend.append('short_cell')
                  else:
                      short_cell_ratio = np.zeros_like(freq)
                      
                  gt.plot_xy('i2_spectrum'
                      ,instrument
                      ,fr.times
                      ,[freq,model_freq,freq]
                      ,[trans,model,short_cell_ratio]
                      ,['r','k','g']
                      ,['None','None','None']
                      ,[1,1,1]       
                      ,['-','-','-']
                      ,[2,2,2]
                      ,['mol/comb','theory']
                      ,'lower left'
                      ,'Frequency'
                      ,'GHz'
                      ,'I2 transmission'
                      ,None
                      ,'I2 spectrum'
                      ,None
                      ,None
                      ,None
                      ,None
                      ,display_defaults
                      ,figs)
          
    if display_defaults.enabled('short_cell_energies')\
            and hasattr(rs,'rs_raw') and hasattr(rs.rs_raw,'nonfiltered_energy')\
           and hasattr(rs.rs_raw,'filtered_energy'):
      
        gt.plot_vs_time('short_cell_energies'
                 ,instrument
                 ,rs.rs_raw.times
                 ,[rs.rs_raw.nonfiltered_energy[:,0]
                   ,rs.rs_raw.filtered_energy[:,0]]  
                 ,['r','b']
                 ,None
                 ,['non-filt','filtered']
                 ,'lower left'    
                ,'Counts'
                 ,''
                 ,'short cell energies'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
   


    
    
    if display_defaults.enabled('mol_cal_vs_short_cell_ratio'):
        fr=None
        if fr is None and hasattr(rs,'rs_mean'):
          fr=rs.rs_mean
        if fr is None and hasattr(rs,'rs_raw'):
          fr=rs.rs_raw
        if fr is None and hasattr(fr,'filtered_energy'):
           #and hasattr(rs,'rs_raw'):
           #and hasattr(rs.rs_raw,'nonfiltered_lockpoint')\
           #and rs.rs_raw.nonfiltered_lockpoint[-1] > 0:
       
          #try:       
            short_cell_ratio = fr.filtered_energy[:,0] \
               / fr.nonfiltered_energy[:,0]
            mcal_pulse = fr.molecular_cal_pulse \
               / fr.transmitted_energy
            mcal_pulse[np.isnan(mcal_pulse)] = 0.0
            
            if hasattr(fr,'filtered_lockpoint')\
                   and hasattr(fr,'nonfiltered_lockpoint'):
                lock_pt = fr.filtered_lockpoint[-1] \
                      / fr.nonfiltered_lockpoint[-1]
            else:
              lock_pt=None
            ymin =np.min(mcal_pulse)
            ymax =np.max(mcal_pulse)
            rv=short_cell_ratio[np.isfinite(short_cell_ratio)]
            if len(rv)>0:
                xmin=np.min(rv)
                xmax=np.max(rv)
            else:
                xmin = 0
                xmax =1.0
                print 'error in lock pt graph --len(short_cell_ratio)==0'
            npoints = len(short_cell_ratio)
            indices = range(npoints)
           
            data_list = zip(indices, short_cell_ratio, mcal_pulse)
            data_list = sorted(data_list, key=itemgetter(1))
            delta_ratio = 0.0025
            i = 0
            mean_ratio = np.NaN * np.zeros((xmax - xmin) / delta_ratio)
            x = np.zeros(len(mean_ratio))
            k = 0
            while i < npoints and k < len(mean_ratio):
                x[k] = xmin + delta_ratio / 2.0
                sum = 0
                n = 0

                while i < npoints and data_list[i][1] < xmin + delta_ratio:
                   sum = sum + mcal_pulse[i]
                   n = n + 1
                   i = i + 1
                   if n > 0:
                        mean_ratio[k] = sum / n
                   else:
                        mean_ratio[k] = np.NaN
                xmin = xmin + delta_ratio
                k = k + 1
            if lock_pt is not None:
                x_vars = [short_cell_ratio, [lock_pt, lock_pt],x]
                y_vars = [mcal_pulse,[ymin,ymax],mean_ratio]
            else:
                x_vars = [short_cell_ratio,x]
                y_vars = [mcal_pulse,mean_ratio] 
            gt.plot_xy('mol_cal_vs_short_cell_ratio'
                ,instrument
                ,fr.times
                ,x_vars
                ,y_vars
                ,['g','r','r']
                ,['o','None','None']
                ,[5,0,0]
                ,['None','-','-']
                ,[0,1,2]
                ,None
                ,None
                ,'short_cell ratio'
                ,None
                ,'molecular counts/energy'
                ,'counts/mJ'
                ,'molecular vs short cell ratio'
                ,None
                ,None
                ,None
                ,None
                ,display_defaults
                ,figs)

        #except Exception , err:
        #    print 'WARNING:  erorr in attempt to plot mol vs short cell ratio'
    
    if display_defaults.enabled('coolant_temperature')\
           and hasattr(rs,'rs_raw')\
           and hasattr(rs.rs_raw,'coolant_temperature'):   
        gt.plot_vs_time('coolant_temperature'
                 ,instrument
                 ,rs.rs_raw.times
                 ,[rs.rs_raw.coolant_temperature]
                 ,['b']
                 ,[2]
                 ,None
                 ,None  
                 ,'Temperature'
                 ,'C-deg'
                 ,'coolant temperature'    
                 ,auto_loop
                 ,display_defaults
                 ,figs)
    
    if display_defaults.enabled('telescope_temperature')\
           and hasattr(rs,'rs_raw')\
           and hasattr(rs.rs_raw,'telescope_temperature'):   
        gt.plot_vs_time('telescope_temperature'
                 ,instrument
                 ,rs.rs_raw.times
                 ,[rs.rs_raw.telescope_temperature]
                 ,['b']
                 ,[2]
                 ,None
                 ,None  
                 ,'Temperature'
                 ,'C-deg'
                 ,'telescope temperature'    
                 ,auto_loop
                 ,display_defaults
                 ,figs)
        
    if display_defaults.enabled('etalon_temperature')\
           and hasattr(rs,'rs_raw')\
           and hasattr(rs.rs_raw,'etalon_temp'):   
        gt.plot_vs_time('etalon_temperature'
                 ,instrument
                 ,rs.rs_raw.times
                 ,[rs.rs_raw.etalon_temp]
                 ,['b']
                 ,[2]
                 ,None
                 ,None  
                 ,'Temperature'
                 ,'C-deg'
                 ,'etalon temperature'    
                 ,auto_loop
                 ,display_defaults
                 ,figs)
        
    if display_defaults.enabled('etalon_pressure')\
           and hasattr(rs,'rs_raw') and hasattr(rs.rs_raw,'etalon_pressure'):   
        gt.plot_vs_time('etalon_pressure'
                 ,instrument
                 ,rs.rs_raw.times
                 ,[rs.rs_raw.etalon_pressure]
                 ,['b']
                 ,[2,]
                 ,None
                 ,None  
                 ,'Pressure'
                 ,'mb'
                 ,'etalon pressure'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
    if display_defaults.enabled('qswitch_buildup_time')\
           and hasattr(rs,'rs_raw')\
             and hasattr(rs.rs_raw,'qswitch_buildup_time'):
        rs.rs_raw.qswitch_buildup_time[
                    rs.rs_raw.qswitch_buildup_time>1.0] =np.NaN
        rs.rs_raw.min_qswitch_buildup_time[
                   rs.rs_raw.min_qswitch_buildup_time>1.0]=np.NaN
        rs.rs_raw.max_qswitch_buildup_time[
                   rs.rs_raw.max_qswitch_buildup_time>1.0]=np.NaN
        rs.rs_raw.builduptime_threshhold[
                   rs.rs_raw.builduptime_threshhold>1.0]=np.NaN
       
        gt.plot_vs_time('qswitch_buildup_time'
                 ,instrument
                 ,rs.rs_raw.times
                 ,[rs.rs_raw.qswitch_buildup_time * 1e6 
                     ,rs.rs_raw.min_qswitch_buildup_time*1e6
                     ,rs.rs_raw.max_qswitch_buildup_time*1e6
                     ,rs.rs_raw.builduptime_threshhold*1e6]
                 ,['r','b','g','k']
                 ,[3,1,1,3]
                 ,['QBT','QBT_min','QBT_max','Threshold']
                 ,'upper left'   
                 ,'QBT'                 ,'us'
                 ,'Q-switch build up time'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
   

    if  display_defaults.enabled('optical_bench_air_pressure') \
           and hasattr(rs,'rs_raw')\
           and hasattr(rs.rs_raw,'opticalbenchairpressure'):      
        gt.plot_vs_time('optical_bench_air_pressure'
                 ,instrument
                 ,rs.rs_raw.times
                 ,[rs.rs_raw.opticalbenchairpressure]
                 ,['b']
                 ,[2]
                 ,None
                 ,None  
                 ,'Pressure'
                 ,'mb'
                 ,'optical_bench_air_pressure'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
    
    if display_defaults.enabled('tcs_currents')\
           and hasattr(rs,'rs_raw'):
        currents = []
        legends = []
        for tcs in ('tcsoptics','tcsopticstop','tcstelescope','thermal1','thermal2','tcsaft','tcsfore'):
          if hasattr(rs.rs_raw,tcs+'_main_current'):
            currents.append(getattr(rs.rs_raw,tcs+'_main_current'))
            legends.append(tcs.replace('tcs',' main'))
          if hasattr(rs.rs_raw,tcs+'_fan1_current'):
            currents.append(getattr(rs.rs_raw,tcs+'_fan1_current'))
            legends.append(tcs.replace('tcs','')+' fan1')
          if hasattr(rs.rs_raw,tcs+'_fan2_current'):
            currents.append(getattr(rs.rs_raw,tcs+'_fan2_current'))
            legends.append(tcs.replace('tcs','')+' fan2')
        if len(currents):    
            gt.plot_vs_time('tcs_currents'
                 ,instrument
                 ,rs.rs_raw.times
                 ,currents
                 ,['r','b','g']
                 ,[2,2,2]
                 ,legends
                 ,'upper left'    
                 ,'Current'
                 ,'Amps'
                 ,'TCS currents'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
        
    if display_defaults.enabled('tcs_temps') \
           and hasattr(rs,'rs_raw'):
        #and hasattr(rs.rs_raw,'tcsoptics_control_temp'):
      for tcs in ('tcsoptics','tcsopticstop','tcstelescope','thermal1','thermal2','tcsaft','tcsfore'):
        temps=[]
        legends=[]
        if hasattr(rs.rs_raw,tcs+'_control_temp'):
            for tm,leg in (('control_temp','control'),('in_temp','in box'),('out_temp','out'),('elect','elect')):
              temps.append(getattr(rs.rs_raw,tcs+'_' + tm))
              legends.append(tcs.replace('tcs','')+' '+leg)
        if len(legends)>0:
          gt.plot_vs_time(tcs+'_temps'
                 ,instrument
                 ,rs.rs_raw.times
                 ,temps
                 ,['r','b','g','c']
                 ,[3,2,2,1]
                 ,legends
                 ,'upper left'    
                 ,'Temperature'
                 ,'deg-C'
                 ,tcs+' temperatures'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
  
    if display_defaults.enabled('humidity') \
           and hasattr(rs,'rs_raw')\
                 and hasattr(rs.rs_raw,'humidity') :
        try:
            lines = [rs.rs_raw.humidity[:,0],rs.rs_raw.humidity[:,1]]
            legend = ['optics','ambient']
        except:
            lines = [rs.rs_raw.humidity]
            legend = []
            
        gt.plot_vs_time('humidity'
                 ,instrument
                 ,rs.rs_raw.times
                 ,lines
                 ,['b','r']
                 ,[2,2]
                 ,legend
                 ,'upper left'
                 ,'Humidity'
                 ,'%'
                 ,'Humidity'
                 ,auto_loop
                 ,display_defaults
                 ,figs)
 
   
    if display_defaults.enabled('one_wire_temps') \
           and hasattr(rs,'rs_raw')\
           and hasattr(rs.rs_raw,'one_wire_temperatures'):    
        temps=[]
        legendList = []
        for i in range(len(rs.rs_raw.one_wire_temperatures[0,:])):
          if len(rs.rs_raw.one_wire_attrib)>i and rs.rs_raw.one_wire_attrib[i] is not None:
            temps.append(rs.rs_raw.one_wire_temperatures[:,i])
            legendList.append(rs.rs_raw.one_wire_attrib[i])
        gt.plot_vs_time('one_wire_temps'
                 ,instrument
                 ,rs.rs_raw.times
                 ,temps
                 ,None
                 ,None
                 ,legendList
                 ,'upper left'    
                 ,'Temperature'
                 ,'deg-C'
                 ,'one wire temperatures'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
   
    if display_defaults.enabled('select_one_wire_temps') \
           and hasattr(rs,'rs_raw')\
             and hasattr(rs.rs_raw,'one_wire_temperatures'):
        temps=[]
        legendList = []
        select_string = display_defaults.get_value(
             'select_one_wire_temps','value')
         
        
        for i in range(len(rs.rs_raw.one_wire_temperatures[0,:])):
            if len(rs.rs_raw.one_wire_attrib)>i and rs.rs_raw.one_wire_attrib[i] is not None and select_string.find(rs.rs_raw.one_wire_attrib[i])>=0:
                temps.append(rs.rs_raw.one_wire_temperatures[:,i])
                legendList.append(rs.rs_raw.one_wire_attrib[i])
                
        if len(temps)>0:
            gt.plot_vs_time('selected_one_wire_temps'
                 ,instrument
                 ,rs.rs_raw.times
                 ,temps
                 ,None
                 ,[2,2,2,2,2,2,2,2]
                 ,legendList
                 ,'upper left'    
                 ,'Temperature'
                 ,'deg-C'
                 ,'one wire temperatures'   
                 ,auto_loop
                 ,display_defaults
                 ,figs)
        else:
            print 'did not find requested one_temperature--no plot'
           
    if display_defaults.enabled('cpol_vs_comb_hi_gain') and hasattr(rs,'rs_mean')\
           or  (display_defaults.enabled('cpol_vs_comb_hi_gain') and hasattr(rs,'rs_raw') \
           and not display_defaults.get_value('cpol_vs_comb_hi_gain','no_raw_plot')):
      fr=None
      pref=''
      for f,p in (('rs_raw','Raw '),('rs_mean','Mean ')):
        if hasattr(rs,f):
          fr=getattr(rs,f)
          pref=p
          for fiel in ('c_pol_dark_counts','c_hi_dark_counts'):
            if not hasattr(fr,fiel):
              fr=None
              break
          if fr is not None:
            break
      if fr is not None:
        #correct for detector dark count if available
        dc_cpol = fr.c_pol_dark_counts[:,0]
        dc_chi  = fr.c_hi_dark_counts[:,0]   
        if rs_constants.has_key('cpol_detector_dark_count'):
             cpol_d_dc = rs_constants['cpol_detector_dark_count']\
                    * fr.seeded_shots
             chi_d_dc = rs_constants['comb_hi_detector_dark_count']\
                    * fr.seeded_shots
        else:
            cpol_d_dc = 0.0 * fr.seeded_shots
            chi_d_dc = cpol_d_dc
            
        #strip nan's
        indices  = np.arange(len(dc_chi))
        mask = np.logical_and(np.logical_and(np.isfinite(dc_chi),np.isfinite(dc_cpol)),
                              np.logical_and(np.isfinite(cpol_d_dc),np.isfinite(chi_d_dc)))
        indices = indices[mask]
        fit = None
        cp_temp = (dc_cpol[indices] - cpol_d_dc[indices]) / fr.seeded_shots[indices]
        chi_temp = (dc_chi[indices]-chi_d_dc[indices]) / fr.seeded_shots[indices]
        
        try:
          if indices.size>0:
            #indices = np.nonzero(indices)             
            max_c_pol_dark_count = np.max(cp_temp)
            x = np.linspace(0, max_c_pol_dark_count, 100)
            pc = np.polyfit(cp_temp
               , chi_temp, 1)
            fit = np.polyval(pc, x)
        except ValueError:
          pass
        if fit is None:
          pc=np.array([np.NaN])
          x=np.array([np.NaN])
          fit=np.array([np.NaN])

      
        if hasattr(rs,'rs_raw'):
            title= 'chi vs cpol(raw), slope= ' + '%4.3f' %(pc[0])
            plot_id = 'cpol_vs_comb_hi_gain(raw)'
        else:
            plot_id = 'cpol_vs_comb_hi_gain(mean)'
            title= 'chi vs cpol(mean), slope= ' + '%4.3f' %(pc[0])    
        gt.plot_xy(plot_id                    #'cpol_vs_comb_hi_gain'
                ,instrument
                ,fr.times
                ,[cp_temp*10000,x*10000.0]     #[fr.c_pol_dark_counts-cpol_d_dc,x]
                ,[chi_temp*10000,fit*10000.0]  #[fr.c_hi_dark_counts-chi_d_dc,fit]
                ,['k','r']
                ,['*','None']
                ,[5,0]
                ,['None','-']
                ,[0,2]
                ,None
                ,None
                ,'Cpol dark counts * 10000'
                ,'/shot/bin'
                ,'Combined hi dark counts *10000'
                ,'/shot/bin'
                ,title
                ,None
                ,None
                ,None
                ,None
                ,display_defaults
                ,figs)
   
   
    if display_defaults.enabled('mode_bits') \
           and hasattr(rs,'rs_raw')\
                and hasattr(rs.rs_raw,'op_mode'):

        bit_names =['i2_scan','i2_lock','cal_scan','filt_1K','filt_1K'
                      ,'cpol_block']
        bv=2048
        bit_value =[1, 4, 16, 32, 128,2**11]
        if hasattr(rs.rs_raw,'i2_cell_out'):
            bit_names.append('no_i2Cell')
            bv=bv*2
            bit_value.append(bv)
        if hasattr(rs.rs_raw,'i2a_cell_out'):
            bit_names.append('no_i2aCell')
            bv=bv*2
            bit_value.append(bv)
        if hasattr(rs.rs_raw,'transmitted_1064_energy'):
            bit_names.append('IR gating off')
            bit_value.append(2**17)
            bit_names.append('IR_filt')
            bit_value.append(2**18)

        bits = np.zeros((len(rs.rs_raw.times),10))
        opmode=rs.rs_raw.op_mode.astype('uint32')
        for i in range(len(bit_value)):
            bits[(opmode[:] & bit_value[i])!=0,i]=1
       
        if hasattr(rs.rs_raw,'i2_cell_out'):
            i=bit_names.index('no_i2Cell')
            bits[:,i] = 0
            bits[rs.rs_raw.i2_cell_out[:]!=0,i] = 1
        if hasattr(rs.rs_raw,'i2a_cell_out'):
            i=bit_names.index('no_i2aCell')           
            bits[:,i] = 0
            bits[rs.rs_raw.i2a_cell_out[:]!=0,i] = 1
        bits[bits == 0] = np.NaN
        
        gt.plot_mode_bits('mode_bits'
                         ,instrument
                         ,rs.rs_raw.times
                         ,bits
                         ,bit_names
                         ,'Mode bits'
                         ,display_defaults
                         ,figs)

    if display_defaults.enabled('selected_sample_profile'):
         if display_defaults.get_value('selected_sample_profile','requested_variable'):
             requested_variable =display_defaults.get_value('selected_sample_profile','requested_variable')
             index = requested_variable.find('.')
             structure_name = requested_variable[:index]
             if hasattr(rs,structure_name):
                 substructure = getattr(rs,structure_name)
                 times = getattr(substructure, 'times')
                 year   = display_defaults.get_value('selected_sample_profile','year')
                 month  = display_defaults.get_value('selected_sample_profile','month')
                 day    = display_defaults.get_value('selected_sample_profile','day')
                 hour   = display_defaults.get_value('selected_sample_profile','hour')
                 minute = display_defaults.get_value('selected_sample_profile','min')
                 second = display_defaults.get_value('selected_sample_profile','sec')
                 requested_time = datetime.datetime(year,month,day,hour,minute,second)
                 index =(times -requested_time)  < datetime.timedelta(seconds =0.0)
                 index = len(index[index == True])-1
                 s = requested_variable.split('.')
                 v= rs
                 for f in s:
                    if hasattr(v,'msl_altitudes'):
                        msl_altitudes = getattr(v,'msl_altitudes')
                    if hasattr(v,f):    
                        v = getattr(v,f)
                    else:
                        print
                        print
                        v=None
                        print RuntimeError('select_sample plot error, ' + requested_variable +' not found') 
                 if v is not None and index>=0:       
                     title_str =s[-1] 
                     xlabel = s[-1]
                     units = display_defaults.get_value('selected_sample_profile','units')
                     gt.plot_vs_altitude('selected_sample_profile'        #display defaults plot name
                         ,instrument                       #instrument name
                         ,times[index]                  #python datetimes vector
                         ,msl_altitudes         #altitude vector (meters)
                         ,[v[index,:]]      #variables
                         ,['r' ]                        #colors, [] default colors
                         ,None                               #widths, [] sets widths = 2
                         ,None                             #legend list, [] = no legend
                         ,'upper right'                    #legend position, [] ok if list []
                         ,xlabel                           #xlabel
                         ,units                            #x units
                         ,title_str                        #plot title
                         ,auto_loop                        # =1, clear figure before new plot
                         ,display_defaults                    
                         ,figs)
    return figs
