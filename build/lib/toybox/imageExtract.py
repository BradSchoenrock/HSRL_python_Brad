import os,sys
from datetime import datetime,timedelta
from PIL import Image
import numpy as np
from netCDF4 import Dataset

def saveAsImage(camname,time,array):
	img=Image.fromarray(array,mode="L")
	fn=camname+time.strftime('_%Y%m%dT%H%M%S.jpg')
	img.save(fn)
	print 'wrote '+fn

def readImages(src,imagelist):
	for nc,fn in src:
		for imname in imagelist:
			if imname+'_snapshot' not in nc.variables:
				print imname+' missing from file '+fn
			imarray=np.array(nc.variables[imname+'_snapshot'])
			imtime=nc.variables[imname+'_snapshot_time']
			imtime=datetime(*imtime[0:6])
			saveAsImage(imname,imtime,imarray)
		yield nc

def fileAsNC(src):
	for fr in src:
		#print fr
		try:
			nc=Dataset(fr['path'],'r')
			yield nc,fr['path']
		except:
			pass

def main(inst,st,en,*imagelist):
    from hsrl.dpl.HSRLLibrarian import HSRLLibrarian
    lib=HSRLLibrarian(instrument=inst)
    src=lib(st,en)
    src=fileAsNC(src)
    src=readImages(src,imagelist)
    for x in src:
    	pass

if __name__ == '__main__':
    p=os.path.realpath(os.path.join(os.path.dirname(sys.argv[0]),'..'))
    print(p)
    sys.path.append(p)
    inst=sys.argv[1]
    startdate=datetime.strptime(sys.argv[2],'%Y%m%dT%H%M%S')#(2015,5,1,0,0,0)
    enddate=datetime.strptime(sys.argv[3],'%Y%m%dT%H%M%S')#datetime(2015,5,18,0,0,0)
    main(inst,startdate,enddate,*sys.argv[4:])
