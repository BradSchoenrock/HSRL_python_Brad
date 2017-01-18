import numpy as np
import lg_base.core.decoratortools as nt

try:
    from bottleneck import anynan
except ImportError:
    def anynan(x):
        return np.any(np.isnan(x))       

def  polynomial_smoothing(array,delta_z,smoothing):
    """If smoothing[0]>0, smooth with running polynomial fit of order smoothing[1]
       width of smoothing width in meters increases linearly from
       smoothing[2] at lowest altitude  too smoothing[3] highest altitude.

       array        = data array to smooth with running polynomial 
       smoothing[0] = enable smoothing if True
       smoothing[1] = order of polynominal to use for local fit
       smoothing[2] = width of smoothing at lowest range (m)
       smoothing[3] = width of smoothing at highest range (m)
       smoothing[4] = first range to smooth
       delta_z      = bin width (m)"""

 
    if smoothing[0] == False:
        print 'no smoothing'
        return array
    
    #check to see if array is 2-d
    try:
       [ntimes,nbins]=array.shape
    except: 
       nbins = len(array)
       ntimes = 1
       
    #compute how much to increment smoothing half-width per altitude index
    delta_w=(smoothing[3]-smoothing[2])/(2.0*delta_z*nbins)
 
    #initial half-width in bins--note this is float
    w0 = smoothing[2]/(2.0*delta_z) 

    #start at larger of polynomial order, or specified start range
    first_bin = np.int(np.float(smoothing[4])/delta_z)
    #first_bin = np.max(first_bin,np.int(smoothing[1]))



    for i in range(ntimes):
       #smooth profiles with local 2nd order polynomial fit
       if ntimes >1:
          temp=array[i,:].copy()
       else:
          temp = array.copy()
          #set NaNs to 0.0
          np.nan_to_num(temp)
          
       #loop over ranges limited by number of points needed to fit polynomial 
       #where the polynomial order = smoothing[1]
   
    
       for bin in range(first_bin+1,nbins-np.int(smoothing[1])-1):  
          w=int(w0+bin*delta_w)
          
          if bin >= w and bin <= nbins-w-2:
              start = bin - w
              end  = bin + w
              #print 'w1 ',start,j,end
             
          elif bin < w:  #i is less than half_width
              start = 0
              end   = 2 * bin
              
          else:   #bin + half_width bumping against nbins
              #"top of profile not smoothed to prevent introduction of extra NaN's"
              a=1

          x=range(start,end+1)
          p=np.polyfit(x,temp[x],np.int(smoothing[1]))

          if not anynan(p):    
              if ntimes>1:
                  array[i,bin]=np.polyval(p,bin)
              else:
                  array[bin]=np.polyval(p,bin)
             
    return array        
