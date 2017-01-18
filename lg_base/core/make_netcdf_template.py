import bz2
from os import system
import os.path


def make_netcdf_template(instrument,in_dir,netcdf_filename,outfilename):
    """reads a netcdf file and creates a outfilename.json
       with a template for reading this netcdf file--this file will
       must be edited to provide name translations
       and variable type information"""

    from pycdf import CDF, CDFError
    if netcdf_filename.find('.bz2') >= 0:
        index = netcdf_filename.find('zraw/')
        index2 = netcdf_filename.find('_')
        fileid = open(in_dir+netcdf_filename, 'r')
        file = fileid.read()
        fileid.close()
        temp = bz2.decompress(file)
        # write temp file to memory rather than disk
        netcdf_filename = '/dev/shm/' + instrument + '_temp_file.nc'
        fileid = open(netcdf_filename, 'w')
        fileid.write(temp)

    nc = CDF(netcdf_filename)
    nc_vars = nc.variables().keys()
    dims = nc.dimensions()
    print dims
    
    f = open(outfilename,'w')
    header  =('#\n'
            +'#edit this file and change extention to *.json\n'
            +'#entry format:"\n'
            +'#"desired_var_name":["name_in_netcdf","select_flag"]\n'
            +'#var_type   == "t" ,time var or time by other than z\n'
            +'#select_flag   when set = 0 and read_mode=="min",variable skipped\n' 
            +'#remove these comment lines before use\n'
            +'#\n'
            +'{"config":\n'
            +'    {"#read_mode":"all / select / min",\n'
            +'    "read_mode": "select"},\n'
            +'"selected_vars":\n'
            +'    {\n')
    f.write(header)
    print header
    first=True
    for name in nc_vars:
        var = nc.var(name)
       
        #print help(nc)
        if  first:
            print  '    "'+name+'":["'+name+'", 1]'
            f.write('    "'+name+'":["'+name+'",1]')
        else:
            print '    "'+name+'":["'+name+'",1]'
            f.write(',\n    "'+name+'":["'+name+'",1]')
        first=False
    f.write('    }\n')
    f.write('}\n')
if __name__ == '__main__':
    make_netcdf_template('gvhsrl','/lidar64/raid/gvhsrldata/2012/02/03/zraw/','gvhsrl_20120203T190000_data_jpg_February03.nc.bz2','../../hsrl_config/bagohsrl_netcdf_defaults.template')
