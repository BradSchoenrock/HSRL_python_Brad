def writeascii(filename,header,data,format):
   #y=writeascii(filename,header,data,format);
   #filename=name of file to be written.
   #header=matlab string containing ascii header 
   #data=matlab array containing data to be written.
   #format=( 1=5.2f, 2=7.3f, 3=10.2e,4=10.3e)

   formats=(None,'%5.2f ','%7.3f ','%10.2e ','%10.3e ')

   #open new file
   fid=open(filename,'w')
   fid.write(header)
   [nr,nc]=data.shape
   fmt=formats[format]
   for i in range(nr):
      for x in range(nc):
        fid.write(fmt % (data[i,x]))
      fid.write('\n');
   fid.close()
   #eval(['!cat<',filename]);
   return
