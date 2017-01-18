import numpy as np
import datetime

def read_processed_hsrl_data(filename,var_name = None ,list_vars = None):
    """read_process_hsrl_data(netcdf_file,variable,list_vars = False)
       read a netcdf variable from a file created by write_netcdf
       if list_variables = True, provide a list of file variables"""

    from pycdf import CDF
    nc=CDF(filename)
    if not list_vars == None:
       nc_vars = nc.variables().keys()
       for name in nc_vars:
           print ' %30s  %40s' %(name.ljust(30), nc.var(name).dimensions())
       return    
    if not var_name == None:
        if var_name == 'time' or var_name == 'times':
            # inits a datetime object from the base time
            base_time = nc.var('base_time').get()
            time_offset = nc.var('time').get()
           #base = datetime.datetime.fromtimestamp(base_time,tz=datetime.tzinfo("UTC"))
            base = datetime.datetime(1970,1,1,0,0,0)+datetime.timedelta(seconds=int(base_time)) 
           # adds a time offset (seconds) to the base time to make a datetime object
            var = np.array([(base + datetime.timedelta(seconds=x)) for x in time_offset])
        else:
            var =nc.var(var_name).get()
        pass
    return var
    
if __name__ == '__main__':
    read_processed_hsrl_data('/home/eloranta/firefox_downloads/mf2hsrl_20130705T1200_20130705T1300_30s_7.5m(1).nc','profile_combined_counts_hi',list_vars= True)
   # read_processed_hsrl_data('/home/eloranta/processed_lidar_data/mf2hsrl_20120928T1200_201209281209.nc','molecular_counts',list_vars= True)

   
