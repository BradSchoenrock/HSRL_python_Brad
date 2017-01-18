A Sphinx-based documentation suite is available in the doc folder. There, you can type 'make html' to construct html files in doc/build/html which are usable on the local hard disk.





"typedarray_hsrl_python"  is the hsrl_python repository.

to launch the lidar data processing code on lidar:

cd ~/typedarray_hsrl_python
from maestro.rti_maestro import Rti

minimum cmd string--with optional parameters in hsrl_config/process_control.json
r = Rti('ahsrl','8-aug-12 13:00', '14:00', 0.00, 12.00)
r = Rti('instrument','start_time','end_time','lo_alt(km),'hi_alt(km'))

optional paramenters that may be specified.
r = Rti('ahsrl','8-aug-12 13:00', '14:00', 0.00, 12.00 , ext_pts_to_ave=1 \
, mol_norm_alt=0.5 , m_smooth=[0, 3, 200] , display='flight_plots.json' \
, make_qc_mask=[1, 1] , ext_bin_delta=1 )

ext_pts_to_ave -> allows multiple point averaging for ext, values = 1, 3, or 9
mol_norm_alt -> alt (km) where molecular signal is normalized to density profile
m_smooth -> provides smoothing of mol signal and inverted molecular return
           = [enable,polynomial_order,smoothing_window_width(m)]
           e.g.  [1, 3 , 500] 500 m window with 3rd order Savitzky-Golay filter
           or    [0, 1,200] no smoothing
make_qc_mask -> generate a data quality mask based on number of molecular photons
           when the number of molecular photons in a bin falls below this level
           mask all subsequent bins in this profile
           = [enable,mol_lost_level]
           eg [1,5]     
default display parameters are provided in a .json file:
    hsrl_config/xxxx_plots.json  
    these are selected via display='xxxx_plots.json'
    you may generate your own files with custom display selections

