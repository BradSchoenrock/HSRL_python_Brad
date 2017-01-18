
def show_corr_adjusts(corr_adjusts,instrument):
 corr_adjusts=corr_adjusts.copy()

 print 'variable      extended name                    scale_factor'
 print 'mol pileup    at high count rates              = ', corr_adjusts.pop('mol_pileup',1.0) 
 print 'c hi pileup                                    = ', corr_adjusts.pop('comb_hi_pileup',1.0)
 print 'c lo pileup                                    = ', corr_adjusts.pop('comb_lo_pileup',1.0)
 print 'mol wfov pileup                                = ', corr_adjusts.pop('mol_wfov_pileup',1.0)
 print 'mol dark corr background subtraction           = ', corr_adjusts.pop('mol_dark_count',1.0)
 print 'comb_hi dark corr background subtraction       = ', corr_adjusts.pop('comb_hi_dark_count',1.0)
 print 'comb_lo dark corr background subtraction       = ', corr_adjusts.pop('comb_lo_dark_count',1.0)
 print 'comb_wfov dark corr background subtraction     = ', corr_adjusts.pop('comb_wfov_dark_count',1.0)
 print 'mol_wfov dark corr background subtraction      = ', corr_adjusts.pop('mol_wfov_dark_count',1.0)
 print 'c_pol dark corr background subtraction         = ', corr_adjusts.pop('c_pol_dark_count',1.0)
 if instrument == 'bagohsrl' or instrument == 'ahsrl':
     print '1064 dark corr background subtraction          = ', corr_adjusts.pop('comb_1064_dark_count',1.0)
 print 'sig_in_dc     signal in dark count cor         = ', corr_adjusts.pop('signal_in_dark',1.0)
 print 'baseline_mol   mol internal scattering         = ', corr_adjusts.pop('mol_baseline',1.0)
 if instrument == 'bagohsrl':
     print 'baseline_i2a_mol                               = ',corr_adjusts.pop('mol_i2a_baseline',1.0)
 if instrument == 'bagohsrl' or instrument == 'ahsrl': 
     print 'baseline_comb_1064                             = ',corr_adjusts.pop('comb_1064_baseline',1.0)
 print 'baseline_chi combined hi internal scattering   = ', corr_adjusts.pop('comb_hi_baseline',1.0)
 print 'baseline_clo combined lo internal scattering   = ', corr_adjusts.pop('comb_lo_baseline',1.0)
 print 'baseline_cpol c pol internal scattering        = ', corr_adjusts.pop('c_pol_baseline',1.0)
 print 'geo corr      mol chan overlap correction      = ', corr_adjusts.pop('geo_corr',1.0)
 print 'i2            adj mol/comb cal ratio           = ', corr_adjusts.pop('i2_corr',1.0)
 if instrument == 'bagohsrl':
     print 'i2a           adj mol_i2a/comb cal ratio       = ', corr_adjusts.pop('i2a_corr',1.0)
 print 'Cam           adjust Cam                       = ', corr_adjusts.pop('Cam_corr',1.0)
 print 'pol_xtalk     polarization cross talk          = ', corr_adjusts.pop('pol_x_talk',1.0)
 print 'dif geo corr, range variation in channel gains = ', corr_adjusts.pop('dif_geo_corr',1.0)
 print 'dif_cpol_geo, range variation in cpol gain     = ', corr_adjusts.pop('cpol_dif_geo',1.0) 
 if instrument == 'bagohsrl' or instrument == 'ahsrl':
    print 'dif 1064/532 geo correction                 = ', corr_adjusts.pop('diff_1064_532_geo_corr',1.0) 
 if instrument == 'bagohsrl':
    print 'i2a dif geo corr, range variation in i2a_mol   = ', corr_adjusts.pop('i2a_dif_geo_corr',1.0)
 if instrument == 'bagohsrl' or instrument == 'ahsrl':
    print '1064 gain                                      = ', corr_adjusts.pop('1064_gain',1.0)   
 print ' '
 print 'instrument: ',instrument,'  ignored parameters:',(','.join(corr_adjusts.keys()))
 print ' '

#input scale factor that can be applied to data corrections
def get_scale_corrections(corr_adjusts,instrument):
#def get_scale_corrections(corr_adjusts):
 show_corr_adjusts(corr_adjusts,instrument)
 #enter scalings for data corrections
 if instrument == 'bagohsrl':
    print 'adjust calibrations: \n chi_pileup=%g, clo_pileup=%g,mol_wfov_pileup=%g \n m_dark=%g, chi_dark=%g, clo_dark=%g, mwfov_dark=%g, cwfov_dark=%g, cp_dark=%g, IR_dark=%g, sig_in_dc=%g \n bl_mol=%g, bl_i2a=%g, bl_chi=%g, bl_clo=%g, bl_cp=%g bl_1064=%g \n geo_corr=%g, i2=%g, i2a=%g, Cam=%g, dif_geo=%g, i2a_dgeo=%g, dif_cpol_geo=%g,dif_1064_532=%g,pol_xtalk=%g, i2a_ratio=%g, 1064_gain=%g  ?' %(
            corr_adjusts.get('comb_hi_pileup',1.0)
           ,corr_adjusts.get('comb_lo_pileup',1.0)
           ,corr_adjusts.get('mol_wfov_pileup',1.0)
           ,corr_adjusts.get('mol_dark_count',1.0)
           ,corr_adjusts.get('comb_hi_dark_count',1.0)
           ,corr_adjusts.get('comb_lo_dark_count',1.0)
           ,corr_adjusts.get('mol_wfov_dark_count',1.0)
           ,corr_adjusts.get('comb_wfov_dark_count',1.0)
           ,corr_adjusts.get('c_pol_dark_count',1.0)
           ,corr_adjusts.get('comb_1064_dark_count',1.0)
           ,corr_adjusts.get('signal_in_dark',1.0)
           ,corr_adjusts.get('mol_baseline',1.0)
           ,corr_adjusts.get('mol_i2a_baseline',1.0)
           ,corr_adjusts.get('comb_hi_baseline',1.0)
           ,corr_adjusts.get('comb_lo_baseline',1.0)
           ,corr_adjusts.get('c_pol_baseline',1.0)
           ,corr_adjusts.get('comb_1064_baseline',1.0)
           ,corr_adjusts.get('geo_corr',1.0)
           ,corr_adjusts.get('i2_corr',1.0)
           ,corr_adjusts.get('i2a_corr',1.0)
           ,corr_adjusts.get('Cam_corr',1.0)
           ,corr_adjusts.get('dif_geo_corr',1.0)
           ,corr_adjusts.get('i2a_dif_geo_corr',1.0)
           ,corr_adjusts.get('cpol_dif_geo',1.0)
           ,corr_adjusts.get('diff_1064_532_geo_corr',1.0)
           ,corr_adjusts.get('pol_x_talk',1.0)
           ,corr_adjusts.get('i2a_ratio',1.0)
           ,corr_adjusts.get('1064_gain',1.0))
    
 elif instrument == 'ahsrl':
    print 'adjust calibrations: \n chi_pileup=%g, clo_pileup=%g,mol_wfov_pileup=%g \n m_dark=%g, chi_dark=%g, clo_dark=%g, mwfov_dark=%g, cwfov_dark=%g, cp_dark=%g, IR_dark=%g, sig_in_dc=%g \n bl_mol=%g, bl_chi=%g, bl_clo=%g, bl_cp=%g, bl_1064=%g, \n geo_corr=%g, i2=%g, i2a=%g, Cam=%g, dif_geo=%g,dif_cpol_geo=%g,dif_1064_532=%g,pol_xtalk=%g, 1064_gain=%g ?' %(
            corr_adjusts.get('comb_hi_pileup',1.0)
           ,corr_adjusts.get('comb_lo_pileup',1.0)
           ,corr_adjusts.get('mol_wfov_pileup',1.0)
           ,corr_adjusts.get('mol_dark_count',1.0)
           ,corr_adjusts.get('comb_hi_dark_count',1.0)
           ,corr_adjusts.get('comb_lo_dark_count',1.0)
           ,corr_adjusts.get('mol_wfov_dark_count',1.0)
           ,corr_adjusts.get('comb_wfov_dark_count',1.0)
           ,corr_adjusts.get('c_pol_dark_count',1.0)
           ,corr_adjusts.get('comb_1064_dark_count',1.0)
           ,corr_adjusts.get('signal_in_dark',1.0)
           ,corr_adjusts.get('mol_baseline',1.0)
           ,corr_adjusts.get('comb_hi_baseline',1.0)
           ,corr_adjusts.get('comb_lo_baseline',1.0)
           ,corr_adjusts.get('c_pol_baseline',1.0)
           ,corr_adjusts.get('comb_1064_baseline',1.0)
           ,corr_adjusts.get('geo_corr',1.0)
           ,corr_adjusts.get('i2_corr',1.0)
           ,corr_adjusts.get('i2a_corr',1.0)
           ,corr_adjusts.get('Cam_corr',1.0)
           ,corr_adjusts.get('dif_geo_corr',1.0)
           ,corr_adjusts.get('cpol_dif_geo',1.0)
           ,corr_adjusts.get('diff_1064_532_geo_corr',1.0)
           ,corr_adjusts.get('pol_x_talk',1.0)
           ,corr_adjusts.get('1064_gain',1.0))
 else:
     print 'adjust calibrations:\n  chi_pileup=%g, clo_pileup=%g, m_dark=%g, chi_dark=%g, clo_dark=%g, mwfov_dark=%g, cwfov_dark=%g, cp_dark=%g \n sig_in_dc=%g, bl_mol=%g, bl_chi=%g, bl_clo=%g, bl_cp=%g \n  geo_corr=%g, i2=%g, Cam=%g, dif_geo=%g, dif_cpol_geo=%g,pol_xtalk=%g' %(
           corr_adjusts.get('comb_hi_pileup',1.0)
           ,corr_adjusts.get('comb_lo_pileup',1.0)
           ,corr_adjusts.get('mol_dark_count',1.0)
           ,corr_adjusts.get('comb_hi_dark_count',1.0)
           ,corr_adjusts.get('comb_lo_dark_count',1.0)
           ,corr_adjusts.get('mol_wfov_dark_count',1.0)
           ,corr_adjusts.get('comb_wfov_dark_count',1.0)
           ,corr_adjusts.get('c_pol_dark_count',1.0)
           ,corr_adjusts.get('signal_in_dark',1.0)
           ,corr_adjusts.get('mol_baseline',1.0)
           ,corr_adjusts.get('comb_hi_baseline',1.0)
           ,corr_adjusts.get('comb_lo_baseline',1.0)
           ,corr_adjusts.get('c_pol_baseline',1.0)
           ,corr_adjusts.get('geo_corr',1.0)
           ,corr_adjusts.get('i2_corr',1.0)
           ,corr_adjusts.get('Cam_corr',1.0)
           ,corr_adjusts.get('dif_geo_corr',1.0)
           ,corr_adjusts.get('cpol_dif_geo',1.0)
           ,corr_adjusts.get('pol_x_talk',1.0))    
 adj_str_full = raw_input('?')                                 
 directives=adj_str_full.strip().split()
 for directive in directives:
     splitdir=directive.split('=')
     if len(splitdir)!=2:
        print "Malformed directive ",directive
        continue
     cmd=splitdir[0].strip()
     value=float(splitdir[1].strip())
     #adj_str=raw_input(prompt_str)
     if cmd=='m_dark':
         corr_adjusts['mol_dark_count']=value
        
     elif cmd=='chi_dark':
         corr_adjusts['comb_hi_dark_count']=value
            
             
     elif cmd=='clo_dark':
         corr_adjusts['comb_lo_dark_count']=value
             
     elif cmd=='cwfov_dark':
         corr_adjusts['comb_wfov_dark_count']=value

     elif cmd=='mwfov_dark':
         corr_adjusts['mol_wfov_dark_count']=value
             
     elif cmd=='cp_dark':
         corr_adjusts['c_pol_dark_count']=value
     elif cmd=='IR_dark':
         corr_adjusts['comb_1064_dark_count']=value
     elif cmd=='bl_mol':
          corr_adjusts['mol_baseline']=value
     elif cmd=='bl_i2a':
          corr_adjusts['mol_i2a_baseline']=value                                
     elif cmd=='bl_chi':
          corr_adjusts['comb_hi_baseline']=value
     elif cmd=='bl_clo':
          corr_adjusts['comb_lo_baseline']=value
     elif cmd=='bl_cp':
          corr_adjusts['c_pol_baseline']=value
     elif (instrument == 'bagohsrl' or instrument == 'ahsrl') and cmd=='bl_1064':
          corr_adjusts['comb_1064_baseline']=value   
     elif cmd=='chi_pileup':
          corr_adjusts['comb_hi_pileup']=value
     elif cmd=='clo_pileup':
          corr_adjusts['comb_lo_pileup']=value
     elif instrument == 'bagohsrl' and cmd=='mol_wfov_pileup':
          corr_adjusts['mol_wfov_pileup']=value
                
     elif cmd=='geo_corr':
          corr_adjusts['geo_corr']=value      
     elif cmd=='dif_geo':
          corr_adjusts['dif_geo_corr']=value
     elif cmd=='dif_cpol_geo':
          corr_adjusts['cpol_dif_geo']=value
     elif (instrument == 'bagohsrl' or instrument == 'ahsrl') and cmd=='dif_1064_532':
            corr_adjusts['diff_1064_532_geo_corr']=value
     elif instrument == 'bagohsrl' and cmd=='i2a_dgeo':
            corr_adjusts['i2a_dif_geo_corr']=value
     elif cmd=='pol_xtalk':
            corr_adjusts['pol_x_talk']=value
            print corr_adjusts['pol_x_talk']    
     elif cmd=='i2':
            corr_adjusts['i2_corr']=value
     elif instrument == 'bagohsrl' and cmd=='i2a':
            corr_adjusts['i2a_corr']=value
     elif cmd=='Cam':
            corr_adjusts['Cam_corr']=value
     elif cmd=='sig_in_dc':
         corr_adjusts['signal_in_dark']=value
     elif instrument == 'bagohsrl' and cmd=='i2a_ratio':
         corr_adjusts['i2a_ratio']=value
     elif (instrument == 'bagohsrl' or instrument == 'ahsrl') and cmd=='1064_gain':
         corr_adjusts['1064_gain']=value
     else:
        print "Unknown directive ",directive
        continue

 show_corr_adjusts(corr_adjusts,instrument)

 return corr_adjusts

