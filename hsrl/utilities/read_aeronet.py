import os
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, num2date
from datetime import datetime
import numpy as np
#from subprocess import call
import zipfile

def read_aeronet(filename,times):

    start_time = times[0] 
    start_time_num = date2num(start_time)
    print 'start_time ',start_time_num,num2date(start_time_num)
    end_time = times[-1] 
    end_time_num = date2num(end_time)
    print 'end_time ',end_time_num, num2date(end_time_num)
    
    if filename[-4:] == '.zip':    
          if not zipfile.is_zipfile(filename):
              print 'filename was not found or is not a zipfile'
          else:
              print filename, ' is a valid zipfile'
          zf = zipfile.ZipFile(filename)
          names = zf.namelist()
          file = zf.getinfo(names[0])
          print dir(file)
          print '%s is %d bytes' % (file.filename, file.file_size)
          f = zf.read(file.filename)
          print f
    
            
  
    for i in range(5):
        header = f.split(',')
        print header
        if len(header)>=13:
             print 'header[7] ',header[7],' header[12]',header[12]
    #skip to first requested time
    print
    print
    time_num = 0
    print 'skip aeronet times prior to requested start time'
    value = f.split('\n')
    i=0
    print len(f.split('\n'))

    while time_num < start_time_num:
        this_line = value[i].split(',')
        if i == len(value)-2:
            print
            print 'No aeronet data found in requested period'
            print
            break

        if len(this_line)>1:
            datestring = this_line[0]+' '+this_line[1]
            if len(datestring.split(':'))>4 :    
                print datestring                                   
                print len(datestring.split(':'))           
                time = datetime.strptime(datestring,"%d:%m:%Y %H:%M:%S")
                time_num = date2num(time)
                print 'time_num = ',time_num ,'time = ',time
        i=i+1    
    
    aeronet_time_nums = []
    aeronet_times = []
    AOT667 = []
    AOT500 = []
    AOT532 = []
    AOT380 = []
    AOT340 = []
    AOT370 = []
    print
    print 'should be at start of requested times'
    if end_time_num < time_num:
            print
            print 'Aeronet data ends prior to end of requested period'
            print
            print
            return
    while time_num <= end_time_num:
        next_data = value[i]
        if len(next_data) == 0:
           break
        this_line = next_data.split(',')    
        datestring = this_line[0]+' '+this_line[1]
        #print datestring, value[12]
        time = datetime.strptime(datestring,"%d:%m:%Y %H:%M:%S")
        time_num = date2num(time)
        #print 'time_num = ',time_num ,'time = ',time
        aeronet_times = np.append(aeronet_times,time)
        aeronet_time_nums = np.append(aeronet_time_nums,time_num)
        ao500  = np.float(this_line[12])
        AOT500 = np.append(AOT500,ao500)
        ao667 = np.float(this_line[6])
        AOT667 = np.append(AOT667,ao667)
        AOT340 = np.append(AOT340,np.float(this_line[18]))
        AOT380 = np.append(AOT380,np.float(this_line[17]))
        ao370 = (370.0-340.0)/(380.0-340.0) * (AOT380[-1]-AOT340[-1]) + AOT340[-1]
        AOT370 = np.append(AOT370,ao370)                                       
        ao532 = (532.0-500.0)/(667.0-500.0) * (AOT667[-1]-AOT500[-1]) + AOT500[-1]
        angstrom_coef =-np.log(ao500 / ao667) / np.log(500.0/667.0)
        ao532_ang = ao500 * (532.0/500.0)**(-angstrom_coef)
        #print 'time = ',time,' Angstrom Coef = ',angstrom_coef,' AOD500= ',ao500,' AOD532 = ',ao532_ang
        AOT532 = np.append(AOT532,ao532_ang)
        #print 'aeronet_times ',times[-1],time_num, AOT500[-1], AOT532[-1],AOT667[-1]
        i=i+1
    AOT500 = np.array(AOT500)
    AOT532 = np.array(AOT532)
    AOT667 = np.array(AOT667)
    AOT370 = np.array(AOT370)                                           
    aeronet_time_nums = np.array(aeronet_time_nums)
    
    return AOT500,AOT532,AOT667,AOT370,aeronet_times

    
if __name__ == '__main__':
    import matplotlib as mpl
    import datetime as dt
    start_time_str = '2013-07-14 05:00:00'
    start_time_dt = dt.datetime.strptime(start_time_str, '%Y-%m-%d %H:%M:%S')
    start_time_num = mpl.dates.date2num(start_time_dt) # 734

    end_time_num = start_time_num +0.2
    read_aeronet('/home/eloranta/firefox_downloads/130713_130713_UAHuntsville.zip',start_time_num,end_time_num)
