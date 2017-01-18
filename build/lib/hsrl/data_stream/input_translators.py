import numpy as np
import lg_base.core.array_utils as hau 
import hsrl.data_stream.processing_utilities as pu
import lg_base.core.read_utilities as hru
from datetime import timedelta
import copy
try:
    from bottleneck import anynan
except ImportError:
    def anynan(x):
        return np.any(np.isnan(x))       

"""
 filterobj=raw_translator('bagohsrl',bagoconstants)
 # filterobj is an instance of the callable class 'raw_translator' with parameters as above
....
 rs_raw #contains data from netcdf
....
 filterobj(rs_raw)
 # now rs_raw has been manipulated by the __call__ function below
"""

class raw_translator:
    def __init__(self,instrument,constants,corr_adjusts):  
        self.instrument   = instrument
        self.constants    = constants
        self.corr_adjusts = corr_adjusts
        #self.ntime_ave    = ntime_ave
        self.unwrap_firstangle = 0.0
        self.unwrap_firstangle_atmagnitude = 0.0

    def update_constants(self,newvals):
        self.constants=newvals       
        
    def __call__(self,raw,cdf_attr):
        """input_translator convert raw netcdf variables into form
           used by the hsrl processing code and preforms pileup
           correction on photon counts""" 

        if hasattr(raw,'wfov_counts') and self.constants['wfov_type'] == 'molecular':
            raw.molecular_wfov_counts = raw.wfov_counts.copy()
        elif hasattr(raw,'wfov_counts'):
            raw.combined_wfov_hi_counts =raw.wfov_counts.copy()

        if hasattr(raw,'op_mode'):
          #extract i2 lock bit from operating mode
          #this will allow testing of bit even after averaging
          raw.i2_locked = (raw.op_mode[:].astype(int)&4)/4

        if hasattr(raw,'seeded_shots'):
          setattr(raw,'delta_t',raw.seeded_shots[:,0]/float(self.constants['laser_rep_rate']))
        else:
          setattr(raw,'delta_t',np.zeros([0]))
        #for i in np.arange(raw.times.size):
        #    raw.times[i]-=timedelta(seconds=raw.delta_t[i])

        if hasattr(raw,'transmitted_energy'):
            # convert to mJ per preaveraged accumulation interval
            raw.transmitted_energy[:] = raw.transmitted_energy \
                *self.constants['transmitted_energy_monitor'][0]\
                +self.constants['transmitted_energy_monitor'][1]\
                *raw.seeded_shots[:,0]
            #compute tranmitted power
            setattr(raw,'transmitted_power',raw.transmitted_energy/raw.delta_t)
            
        if hasattr(raw,'transmitted_1064_energy'):
            # convert to mJ per preaveraged accumulation interval
            raw.transmitted_1064_energy[:] = raw.transmitted_1064_energy \
                *self.constants['transmitted_1064_energy_monitor'][0]\
                +self.constants['transmitted_1064_energy_monitor'][1]\
                *raw.seeded_shots[:,0]
            #compute tranmitted 1064 power
            setattr(raw,'transmitted_1064_power',raw.transmitted_1064_energy/raw.delta_t)
             
        if hasattr(raw,'filtered_energy'):
            if raw.filtered_energy.dtype == 'int32':
                raw.nonfiltered_energy = raw.nonfiltered_energy.astype('float64')
                raw.filtered_energy = raw.filtered_energy.astype('float64')
            if len(raw.filtered_energy.shape) == 1:        
                raw.filtered_energy = raw.filtered_energy[:,np.newaxis]
                raw.nonfiltered_energy = raw.nonfiltered_energy[:,np.newaxis]
            raw.filtered_energy[raw.filtered_energy > 1e10] = np.NaN
            raw.nonfiltered_energy[raw.nonfiltered_energy > 1e10] = np.NaN
          
        if hasattr(raw,'builduptime') and raw.builduptime.size > 0:
            raw.qswitch_buildup_time=raw.builduptime[:,0]
            raw.min_qswitch_buildup_time = raw.builduptime[:,1]
            raw.max_qswitch_buildup_time = raw.builduptime[:,2]
       
        if hasattr(raw,'superseedlasercontrollog'):
            raw.superseedlasercontrollog[raw.superseedlasercontrollog > 1e10] = np.NaN
            
        if hasattr(raw,'energyRatioLockPoint') \
               and raw.energyRatioLockPoint.size>0:
            if len(raw.energyRatioLockPoint.shape)==2:
              raw.filtered_lockpoint = raw.energyRatioLockPoint[:,0]
              raw.nonfiltered_lockpoint = raw.energyRatioLockPoint[:,1]
            else:
              clarray=hau.T_Array(np.ones([raw.filtered_energy.shape[0]]))
              raw.filtered_lockpoint = clarray*raw.energyRatioLockPoint[0]
              raw.nonfiltered_lockpoint = clarray*raw.energyRatioLockPoint[1]
        
        if hasattr(raw,'raw_analog_interferometertemperature'):
             thermistor_cal = self.constants['interferometer_temp_cal']
             R = np.abs(raw.raw_analog_interferometertemperature / 0.000250)
             raw.interferometer_temp = 1/(thermistor_cal[0] + thermistor_cal[1] \
                     * np.log(R) + thermistor_cal[2] * np.log(R)** 3) - 273.15
           
             ntemps = len(raw.interferometer_temp)
             if 0: #ntemps > 5:
                 #do eleven point median filter
                 temps=np.zeros((ntemps,11))
                 temps[0:ntemps-5,0] = raw.interferometer_temp[5:]
                 temps[0:ntemps-4,1] = raw.interferometer_temp[4:]
                 temps[0:ntemps-3,2] = raw.interferometer_temp[3:]
                 temps[0:ntemps-2,3] = raw.interferometer_temp[2:]
                 temps[0:ntemps-1,4] = raw.interferometer_temp[1:]
                 temps[:,5] = raw.interferometer_temp
                 temps[1:,6] = raw.interferometer_temp[:ntemps-1]
                 temps[2:,7] = raw.interferometer_temp[:ntemps-2]
                 temps[3:,8] = raw.interferometer_temp[:ntemps-3]
                 temps[4:,9] = raw.interferometer_temp[:ntemps-4]
                 temps[5:,10]= raw.interferometer_temp[:ntemps-5]
                 raw.interferometer_temp = hau.T_Array(np.median(temps,1))
             else:
                 raw.interferometer_temp = hau.T_Array(raw.interferometer_temp)
             
        if hasattr(raw,'raw_analog_etalontemperature'):
                # convert etalon thermistor voltage to themistor resistance
                # T(degC) =1/( a + b(Ln R) + cLn R)^3)-273.15
                #(Steinhart  Hart equation)
                # Where:
                # a = 0.000862448
                # b = 0.000258456
                # c = 0.000000142
                # and
                # R = (Volts ADC Reading)/(0.000250 amps)
            thermistor_cal = self.constants['interferometer_temp_cal']
            R = np.abs(raw.raw_analog_etalontemperature / 0.000250)
            raw.etalon_temp = (hau.T_Array(np.array((1.0 / (thermistor_cal[0] + thermistor_cal[1]
                * np.log(R) + thermistor_cal[2] * np.log(R)** 3) - 273.15),dtype=np.float32,ndmin=1)))  
        if hasattr(raw,'raw_analog_coolanttemperature'):
                # convert coolant thermistor voltage to themistor resistance
                # T(degC) =1/( a + b(Ln R) + cLn R)^3)-273.15
                #(Steinhart  Hart equation)
                # Where:
                # a = 0.000862448
                # b = 0.000258456
                # c = 0.000000142
                # and
                # R = (Volts ADC Reading)/(0.000250 amps)
            thermistor_cal = self.constants['interferometer_temp_cal']
    
            R = np.abs(raw.raw_analog_coolanttemperature / 0.000250)
           

            raw.coolant_temperature = \
                  (hau.T_Array(np.array((1.0 / (thermistor_cal[0] + thermistor_cal[1]
                  * np.log(R) + thermistor_cal[2] * np.log(R)** 3)
                  - 273.15),dtype=np.float32,ndmin=1)))

        if hasattr(raw,'telescope_pointing'):
            if not hasattr(raw,'telescope_locked'):
                setattr(raw,'telescope_locked',np.ones_like(raw.telescope_pointing))
            raw.telescope_pointing=raw.telescope_pointing.astype('float64')
            raw.telescope_pointing[raw.telescope_locked==0]=.5
            #roll component of telescope mounting angle in degrees measured relative
            #to platform (zero degrees = vertical)
            #roll angle is + in clockwise direction
            if not hasattr(raw,'telescope_roll_angle_offset'):
                setattr(raw,'telescope_roll_angle_offset',np.ones_like(raw.telescope_pointing))
            raw.telescope_roll_angle_offset[:] = self.constants['telescope_roll_angle_offset']
            raw.telescope_roll_angle_offset[raw.telescope_pointing == 0] = \
                                                180.0 - self.constants['telescope_roll_angle_offset']

        if hasattr(raw,'raw_analog_telescope_temperature'):
                # convert coolant thermistor voltage to themistor resistance
                # T(degC) =1/( a + b(Ln R) + cLn R)^3)-273.15
                #(Steinhart  Hart equation)
                # Where:
                # a = 0.000862448
                # b = 0.000258456
                # c = 0.000000142
                # and
                # R = (Volts ADC Reading)/(0.000250 amps)
            thermistor_cal = self.constants['interferometer_temp_cal']
    
            R = np.abs(raw.raw_analog_telescope_temperature / 0.000250)
           

            raw.telescope_temperature = \
                  (hau.T_Array(np.array((1.0 / (thermistor_cal[0] + thermistor_cal[1]
                  * np.log(R) + thermistor_cal[2] * np.log(R)** 3)
                  - 273.15),dtype=np.float32,ndmin=1)))
        if hasattr(raw,'OutgoingBeamPosition_centermass')\
               and raw.OutgoingBeamPosition_centermass.size > 0 :
            raw.cg_xs = raw.OutgoingBeamPosition_centermass[:,0]
            raw.cg_ys = raw.OutgoingBeamPosition_centermass[:,1]

        if hasattr(raw,'OutgoingBeamPosition2_centermass')\
                  and raw.OutgoingBeamPosition2_centermass.size > 0 :
            raw.cg_xs2 = raw.OutgoingBeamPosition2_centermass[:,0]
            raw.cg_ys2 = raw.OutgoingBeamPosition2_centermass[:,1]
            
        if hasattr(raw,'interferometer_intensity') \
               and raw.interferometer_intensity.size > 0:
            interf_peak = \
                self.constants['interferometer_spectral_peak']
            phase_to_freq = \
                self.constants['interferometer_phase_to_freq']
            npixels=self.constants['interferometer_fft_npixels']
            xform = np.fft.rfft(raw.interferometer_intensity[:,:npixels], axis=1)
            tmp = np.concatenate(([self.unwrap_firstangle],np.angle(xform[:, interf_peak])))
            newlast=tmp[-1]
            tmp = np.unwrap(tmp)
            tmp = (self.unwrap_firstangle_atmagnitude-tmp[0])+tmp
            if np.isfinite(tmp[-1]):
                self.unwrap_firstangle_atmagnitude=tmp[-1]
                self.unwrap_firstangle=newlast
            raw.interf_freq = tmp[1:]
            raw.interf_freq = hau.T_Array(-raw.interf_freq * phase_to_freq[0])
          
        #compute temperature compensated interferometer freq 
        if 0 and hasattr(raw,'interferometer_temp') \
               and hasattr(raw,'interf_freq')\
               and self.constants.has_key('interf_temp_coef'):
            raw.tcomp_interf_freq = raw.interf_freq \
                        - (raw.interferometer_temp-raw.interferometer_temp[0])\
                         * self.constants['interf_temp_coef']*1e9
        for imagetime in ('interferometer_snapshot_time','outgoingbeamalignment_snapshot_time','overhead_snapshot_time','snowscope_snapshot_time'):
            if hasattr(raw,imagetime):
                setattr(raw,imagetime,hru.convert_to_python_times(getattr(raw,imagetime)[np.newaxis,:]))
            
        #replace missing values witn NaN's
        if hasattr(raw,'seedvoltage'):    
            raw.seedvoltage[raw.seedvoltage > 100] = np.NaN
        if hasattr(raw,'latitude'):    
            raw.latitude[raw.latitude > 100] = np.NaN
        if hasattr(raw,'longitude'):    
            raw.longitude[raw.longitude > 200] = np.NaN
            
        if hasattr(raw, 'laserpowervalues') and raw.laserpowervalues.size > 0:
            raw.laser_current = raw.laserpowervalues[:,0]
            raw.laser_voltage = raw.laserpowervalues[:,1]
            if raw.laserpowervalues.shape[1] > 2:
                raw.laser_current_setpoint = raw.laserpowervalues[:,2]
                raw.laser_diode_temp = raw.laserpowervalues[:,3]
                raw.laser_diode_temp_setpoint = raw.laserpowervalues[:,4]
            if raw.laserpowervalues.shape[1] > 6:
                raw.ktp_temp = raw.laserpowervalues[:,5]
                raw.ktp_temp_setpoint = raw.laserpowervalues[:,6]

        #remove spikes from tcs records
        for fiel in ('tcsopticstop_','tcsoptics_','tcstelescope_','thermal1_','thermal2_','tcsaft_','tcsfore_'):
            for f in vars(raw).keys():
                if f.startswith(fiel):
                    v=getattr(raw,f)
                    v[v>1000]=np.NaN
            
        if hasattr(raw,'one_wire_temperatures') \
                and raw.one_wire_temperatures.size >0 :

            #raw.one_wire_attrib = cdf_attr['one_wire_temperatures']
            raw.one_wire_attrib = []
            [ntime,ntemps]=raw.one_wire_temperatures.shape
            for i in range(ntemps):
                string = 'field'+ str(i)+'_name'
                try:
                  raw.one_wire_attrib.append(
                    cdf_attr['one_wire_temperatures'][string])
                except KeyError:
                  print "Couldn't find attribute for ",string
                  raw.one_wire_attrib.append(None)
            #remove spikes of 1e37 that appear in temperatues
            raw.one_wire_temperatures[raw.one_wire_temperatures>1000.]\
                         =np.NaN
        if hasattr(raw,'RemoveLongI2Cell'):
            servo_range = cdf_attr['RemoveLongI2Cell']['range']            
            raw.i2_cell_out = np.abs(raw.RemoveLongI2Cell-servo_range[1]) \
                     > np.abs(raw.RemoveLongI2Cell-servo_range[0])
            
        if hasattr(raw,'RemoveLongI2ArCell'):
            servo_range = cdf_attr['RemoveLongI2ArCell']['range']            
            raw.i2a_cell_out = np.abs(raw.RemoveLongI2ArCell-servo_range[1]) \
                     > np.abs(raw.RemoveLongI2ArCell-servo_range[0])    
           
        if hasattr(raw,'shot_count'):
          if raw.shot_count.size > 0:
            raw.shot_count= raw.shot_count[:,0]
          else:
            raw.shot_count= raw.shot_count.reshape([0])
            
        if hasattr(raw,'seeded_shots'):
          if raw.seeded_shots.size > 0 :
            raw.seeded_shots=raw.seeded_shots[:,0]
          else:
            raw.seeded_shots=raw.seeded_shots.reshape([0])


        #extract average dark counts from profiles and add dark counts to raw
        #dark count extracted from 'first_bins' or 'last_bins' as specified in constants
        #pu.extract_dark_count(raw,self.constants) #moved to after PILEUP 20140805

        #extract cal pulse from light scattered as laser pulse exits system
        #and place in raw
        #pu.extract_cal_pulse(raw,self.constants)


        if hasattr(raw,'l3cavityvoltage') and raw.l3cavityvoltage.size >0:
           raw.piezo_voltage_ave = raw.l3cavityvoltage[:,0]
           raw.piezo_voltage_min = raw.l3cavityvoltage[:,1]
           raw.piezo_voltage_max = raw.l3cavityvoltage[:,2]
        if hasattr(raw,'l3locking_stats') and 'l3slope_to_frequency' in self.constants:
           raw.l3frequency_offset = raw.l3locking_stats.copy()
           for x in range(0,3):
             raw.l3frequency_offset[:,x] = np.polyval(self.constants['l3slope_to_frequency'], raw.l3locking_stats[:,x])

        if hasattr(raw,'GPS_MSL_Alt'):
           #replace and spikes in a altitude with base altitude
           #this allows the code to run but produces garbage data
           if np.any(raw.GPS_MSL_Alt > 20000.0):
               raw.GPS_MSL_Alt[raw.GPS_MSL_Alt > 20000.0] = \
                          self.constants['lidar_altitude']
        if hasattr(raw,'roll_angle'):
            if anynan(raw.roll_angle):
                raw.roll_angle[np.isnan(raw.roll_angle)]= 0.0        
            if anynan(raw.pitch_angle):
                raw.pitch_angle[np.isnan(raw.pitch_angle)]=0.0
        if hasattr(raw,'opticalbenchairpressure'):
            #convert psi to mb
            #print 'pre--opticalbenchairpressure ',raw.opticalbenchairpressure.shape
            raw.opticalbenchairpressure=hau.T_Array(np.array((raw.opticalbenchairpressure\
                    *self.constants['optical_bench_air_pressure_cal']),ndmin=1,dtype=np.float32))
            #print 'optical_bench_air_pressure.size',raw.opticalbenchairpressure.size,raw.times.size
        if hasattr(raw,'chillertemperature') and raw.chillertemperature.size > 0:
                raw.chiller_temp = raw.chillertemperature[:,0]
                raw.chiller_setpt =raw.chillertemperature[:,1]
              
        if hasattr(raw,'etalon_pressure'):
            raw.etalon_pressure = raw.etalon_pressure * self.constants['etalon_pressure']
           
        if hasattr(raw,'qw_rotation_angle'):
            #convert gv quarter wave plate rotation angle from radians to deg
            raw.qw_rotation_angle=raw.qw_rotation_angle*180.0/np.pi

        if hasattr(raw,'GPS_MSL_Alt') or self.constants['installation'] in ('airborne','shipborne') :
            #do quality check on aircraft GPS and INS data
            pu.gps_quality_check(raw,self.constants)

        if hasattr(raw,'molecular_counts'):
          for k,v in vars(raw).items():
            if '_counts' in k:
              if raw.molecular_counts.shape[1]!=v.shape[1]:
                print 'raw field ',k,' is messed up. size difference'
                tmp=copy.deepcopy(raw.molecular_counts)
                minidx=min(tmp.shape[1],v.shape[1])
                tmp[:,:]=0
                tmp[:,:minidx]=v[:,:minidx]
                setattr(raw,k,tmp)

        #do pileup correction on signals before any averaging is done
        pu.pileup_correction(raw,self.constants,self.corr_adjusts)
        
        #extract average dark counts from profiles and add dark counts to raw
        #dark count extracted from 'first_bins' or 'last_bins' as specified in constants
        if hasattr(raw,'molecular_counts'):
          pu.extract_dark_count(raw,self.constants) #relocated from above 20140805

          #extract cal pulse from light scattered as laser pulse exits system
          #and place in raw
      
          pu.extract_cal_pulse(raw,self.constants)
          
          if 0:
            import hsrl.simulation.rb_simulation as sim
            #rescale for new energy  
            sim.rb_simulation(raw,self.constants)
            #redo dark count
            pu.extract_dark_count(raw,self.constants)
