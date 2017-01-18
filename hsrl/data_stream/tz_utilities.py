#tz_utilities.py


""" The time_z_data structure must contain a vector "time_z_data.times"
    containing time entries express in matplotlib time numbers with
    dimensions [ntimes,0]. It may also contian an altitude vector
    "time_z_data.altitudes" with dimensions [0,nalts], Other arrays
    must have individual entries corresponding to these times and
    altitudes, with columuns refering to time and rows to altitudes.
    Entries with time and altitude dimensions must be of subclass "tz_array".
    Entries with time but no altitude dimension must be of subclass "t_array".
    Entries with altitude but no time dimension must be of subclass "z_array"
    Items that are not members of these subclasses may be included--they
    are passed through routines without change."""

import numpy as np
# Try to use the much faster nanmean from bottleneck, otherwise fall back
# to the scipy.stats version
try:
    from bottleneck import nanmean
except ImportError:
    print
    print "No bottleneck.nanmean available! Falling back to SLOW scipy.stats.nanmean"
    print
    from scipy.stats import nanmean

class time_z_data(object):
    pass

#definitions of lidar data array classes
#test_data = array(arange(1104, dtype=float))
#aa = altitude_array.from_array( test_data )
#assert( isinstance(aa, altitude_array) )
#from pprint import pprint
#pprint(aa)

class t_array( np.ndarray ): 
    """a 2-d array which contains time information
    in the first dimension"""
    @staticmethod
    def from_array(data):
        "copy numpy data into a new altitude_array"
        aa = t_array(data.shape, dtype=data.dtype)
        aa[:] = data[:]
        return aa
class z_array( np.ndarray ): 
    """a 2-d array which contains altitude information
    in the second dimension"""
    @staticmethod
    def from_array(data):
        "copy numpy data into a new altitude_array"
        aa = z_array(data.shape, dtype=data.dtype)
        aa[:] = data[:]
        return aa


class tz_array( np.ndarray ): 
    """a 2-d array which contains time information in the
    first dimension and altitude information in the
    second dimension"""
    @staticmethod
    def from_array(data):
        "copy numpy data into a new altitude_array"
        aa = tz_array(data.shape, dtype=data.dtype)
        aa[:] = data[:]
        return aa




#def merge_time_vectors( D,D_source):
def  merge_time_vectors(D,D_source):
    """Copies all of the 1-d time vectors and non-array items from D_source
       and merges them into D"""

    # create a blank dictionary to return
    R = {}
     
    # copy D into R
    for (name, value) in vars(D).items():
      R[name]=value 

    # copy time vectors and non-ndarray types to D
    # from D_source
    for (name, value) in vars(D_source).items():  
      if      isinstance(value,t_array) or \
              isinstance(value,z_array) or \
              isinstance(value,tz_array):
         if isinstance(value,t_array):
                 R[name]=value
      else:  
          R[name]= value

    # create a new object of the same type as D
    dtype = type(D)
    Robj = dtype()
    # transplant the attributes to the new object, which is technically sneaky
    vars(Robj).update(R)

    return Robj

def nanmean_altitudes(D,new_alts):
    """def nanmean_altitudes(D,new_alts):   
       create altitude means of time_z object centered on new altitude vector
       D        = time_z_object
       new_alts = altitude vector, center nameans of D on these altitudes
       returns smaller array with average values that properly treat NaN's
       It also checks for D.nave and multiplies entries according to the
       number of alt bins averaged together at that alititude"""

    # create a blank dictionary to return
    R = { }
    #create altitude vector centered on supplied new_alts vector
    alts=np.zeros_like(new_alts)
    alts[0,1:]=(new_alts[0,1:]+new_alts[0,:-1])/2
    ind=np.arange((D.msl_altitudes.shape[1]))
    #make a list of indices to indicate which elements to average for each entry
    #in new alts vector
        
    k=0
    temp=D.msl_altitudes[0,:]
    indices=[]
    while k< alts.shape[1]-1:  
      inxa = ind[temp >alts[0,k]]
      temp2=D.msl_altitudes[0,inxa]
      inx=inxa[temp2 <= alts[0,k+1]]
      indices.append(inx)
      k=k+1
     
    # iterate through all the items in a dictionary
    for (name, value) in vars(D).items():
        if name == 'var_mol':  #check if photon count stats requested
            n_ave=np.zeros((1,k+1))
            #store number of points averaged in altitude
            for i in range(k):
               n_ave[0,i]=len(indices[i])   
        if isinstance(value,t_array) or\
           isinstance(value,z_array) or\
           isinstance(value,tz_array):
        
          #if this is the altitude vector replace with new alts  
          if name == 'msl_altitudes':
            R[name]=new_alts
          elif name == 'times':
            R[name]=value
          else:  #not the msl_altitude vector  
            #if altitude vector 
            if isinstance(value,z_array):
               k=0
               temp=z_array.from_array(np.zeros((1,new_alts.shape[1])))
               while k <= (new_alts.shape[1]-2):
                 temp[0:1,k]=nanmean(value[0,indices[k]])
                 k=k+1
               R[name]=temp
            #elif value.shape[1]==1:    
            elif isinstance(value,t_array): #if time array just pass through
               R[name] = D.__dict__[name]  

            #if time-altitude array
            
            elif isinstance(value,tz_array):
               # do for each profile
               nshots=value.shape[0]
               temp=tz_array.from_array(np.zeros((nshots,new_alts.shape[1])))
               k=0 
               #do the alt means
               while k <= (new_alts.shape[1]-2):
                   temp[:,k]=nanmean(value[:,indices[k]],1)
                   k=k+1  
               R[name]=temp    
        else:
            
            #if not numpy array, just pass through
            R[name] = D.__dict__[name]

            
    # create a new object of the same type as D
    dtype = type(D)
    Robj = dtype()
    # transplant the attributes to the new object, which is technically sneaky
    vars(Robj).update(R)
    
    if 'n_ave' in locals():
        Robj.n_ave=n_ave
        
    return Robj

def nanmean_time_z_data_object(D,kt,ka):
    """def nanmean_time_z_data_object( D, kt,ka ):   
       D = time_z_object
       kt= number of times to ave together
       ka= number of altitudes to ave together
       warning--no treatment of remainder if ka or kt does not evenly divide buffer
       returns smaller array with average values that properly treat NaN's"""
  
    # create a blank dictionary to return
    R = { }
    
    # iterate through all the items in a dictionary
    for (name, value) in vars(D).items():
        #if this is a numpy array
        #if isinstance(value,np.ndarray):
        if isinstance(value,t_array) or\
                isinstance(value,z_array) or\
                isinstance(value,tz_array):           
            #if altitude vector
            if isinstance(value,z_array):
               k=0
               temp=np.zeros((1,int(value.shape[1]/ka)))
               while ka*(k+1) <= value.shape[1]:
                 temp[0:1,k]=nanmean(value[0,ka*k:ka*(k+1)])
                 k=k+1
               R[name]=z_array.from_array(temp)  
               
            #if time vector  
            elif isinstance(value,t_array): 
               k=0
               temp=np.zeros((int(value.shape[0]/kt),1))
               while kt*(k+1) <= value.shape[0]:
                 temp[k,0:1]=nanmean(value[kt*k:kt*(k+1),0:1])
                 k=k+1
               R[name]=t_array.from_array(temp)
               
            #if time-altitude array
            else:
               [nshots,nalts]=value.shape
               #start with time mean 
               k=0
               temp=np.zeros((int(value.shape[0]/kt),nalts))               

               if nshots >=1:  
                 while kt*(k+1) <= nshots: 
                   temp[k,0:nalts]=nanmean(value[kt*k:kt*(k+1),0:nalts],0)
                   k=k+1
               else:
                   temp[1,0:nalts]=value
               # proceed with altitude mean
               nshots=k
               k=0
               temp2=np.zeros((int(value.shape[0]/kt),int(value.shape[1]/ka)))
               while ka*(k+1) < nalts:
                   temp2[0:nshots,k]=nanmean(temp[0:nshots,ka*k:ka*(k+1)],1)
                   k=k+1  
               R[name]=tz_array.from_array(temp2)
        else:
            #if not numpy array, just pass through
            R[name] = D.__dict__[name]
            
    # create a new object of the same type as D
    dtype = type(D)
    Robj = dtype()
    # transplant the attributes to the new object, which is technically sneaky
    vars(Robj).update(R)

   
    return Robj



def append_time_z_data_object(D,T):
    """def append_time_z_data_object( D, T ):
       Append time_z_data_object D to time_z_data_object T"""
    #create a blank dictionary to return
    R = { }
    D_ok=True
    try:
       D.times.shape
    except AttributeError:
      D_ok=False
      
    if D_ok:
       # T.times.shape produces an error when T is empty,
       #this allows appending onto an empty array.
       T_ok=True
       try:
           T.times.shape
       except AttributeError:
           T_ok=False  
       if T_ok:
         # iterate through all the items in a dictionary
         for (name, value) in vars(D).items():
            if isinstance(value,tz_array):
                if T.__dict__[name].shape[1] == value.shape[1]:
                    R[name] = tz_array.from_array( \
                         np.vstack((T.__dict__[name],D.__dict__[name])))
                else: #different numbers of altitude bins
                    T_t_bins=T.__dict__[name].shape[0]
                    T_z_bins=T.__dict__[name].shape[1]
                    D_t_bins=value.shape[0]
                    D_z_bins=value.shape[1]
                    if  T_z_bins > D_z_bins: 
                       #make temp with time dim of D and alt dim of T
                       temp=np.NaN*np.ones((D_t_bins,T_z_bins))
                       temp[:,:D_z_bins]=value[:,:D_z_bins]
                       R[name] = tz_array.from_array( \
                         np.vstack((T.__dict__[name],temp)))
                    else:
                       
                       #make temp with time dim of T and alt dim of D
                       temp=-1e10*np.ones((T_t_bins,D_z_bins))
                       temp[:,:T_z_bins]=T.__dict__[name][:,:T_z_bins]
                       R[name] = tz_array.from_array( \
                            np.vstack((temp,value)))
            elif isinstance(value,t_array):
                R[name] = t_array.from_array( \
                    np.vstack((T.__dict__[name],D.__dict__[name])))
            else:
                R[name] = D.__dict__[name]
         # create a new object of the same type as D
         dtype = type(D)
         Robj = dtype()
         # transplant the attributes to the new object, which is technically sneaky
         vars(Robj).update(R)
       #except:
       else:
         Robj=D  #if T is empty output D
    else:     #if D is empty output T
         Robj=T
    return Robj

#FIXME this looks like it still uses matplotlib datenums
#def select_time_z_data_object( D,time ):
def select_time_z_data_object(D,time):
    """select object members nearest to a given time
       time_z objects include 2-d numpy arrays that must include D.times
       with times in matplotlib elapsed days format. First dim entries must
       correspond to times and 2nd dim entries to altitude.
       Elements that are not of type numpy.ndarray are copied to output 
       without change. A field bounds is added to the object to indicate
       if requested time is between the first and last entry in D.times.
       bounds is -1 before the first entry, 0 within time period and +1
       after last entry. A field expire_time is added to indicate a point
       10% beyound the half-way point between this entry and the next entry"""
    
    
    dt=abs(D.times-time)
    inx=dt.argmin()
    ntimes=D.times.shape[0]
    expire_time=np.empty([1,1])
    if D.times[0,0] > time:
      bounds=-1     #time is prior to all elements in object
      if ntimes>1:
         expire_time[0,0]=(D.times[0,0]+D.times[1,0])/2.0 \
                           +(D.times[1,0]-D.times[0,0])/10.0
      else:
         expire_time[0,0]=D.times[0.0]+0.01
    elif D.times[-1,0] < time:
      bounds=+1     #time is after all elements in object
      expire_time[0,0]=time+1/8.0
    else:  
      bounds=0      #time is within bounds
      if inx < ntimes-2:
          if time <= (D.times[inx,0]+D.times[inx+1,0])/2.0:
              #expire at half-way point + a little
              expire_time[0,0]=(D.times[inx,0]+D.times[inx+1,0])/2.0 \
                        +(D.times[inx+1,0]-D.times[inx,0])/10.0   
          else:  #beyound half-way point               
              expire_time[0,0]=(D.times[inx+1,0]+D.times[inx+2,0])/2.0 \
                        +(D.times[inx+1,0]-D.times[inx,0])/10.0
      elif inx == ntimes-2:
             if time < (D.times[inx,0]+D.times[inx+1,0])/2.0:
                  expire_time[0,0]=(D.times[inx,0]+D.times[inx+1,0])/2.0 \
                        +(D.times[inx+1,0]-D.times[inx,0])/10.0
             else:
                 expire_time[0,0]=D.times[-1,0]+1/1.9
      else: #inx=ntimes-1
         expire_time[0,0]=D.times[-1,0]+1/1.9
    # create a blank dictionary to return
    R = { }
    
    # iterate through all the items in a dictionary
    for (name, value) in vars(D).items():
        if isinstance(value,t_array) or \
           isinstance(value,tz_array):
          R[name] = value[np.newaxis,inx,:]
        else:  
          R[name] = value 
    # create a new object of the same type as D
    dtype = type(D)
    Robj = dtype()
    # transplant the attributes to the new object, which is technically sneaky
    vars(Robj).update(R)
    Robj.bounds_flag=bounds
    Robj.expire_time=expire_time   
    return Robj

#def trim_time_interval( D, start_time, end_time ):
def  trim_time_interval(D,start_time,end_time):
    """trim all the fields in a dictionary to [start_time>=D.times<=end_time]"""

    #create a blank dictionary to return
    R = { } 

    #find start_index and end_index
    [nshots,nalts]=D.times.shape
    indices=np.arange(nshots)
    time_mask=D.times[:,0]>=start_time
    if all(time_mask==0):
      print 'tz_utilities--no shots in requested time interval'  
      return R
    start=min(indices[time_mask])
    time_mask=D.times[:,0]<=end_time
    #test for empty array
    if len(time_mask)==0:
      print 'here'  
      return R
    check_ind = indices[time_mask]
    # does indices contain anything? 
    if len(check_ind) ==  0:
      return R

    end=max(indices[time_mask])+1
    
    # iterate through all the items in a dictionary
    for (name, value) in vars(D).items():  
      if    isinstance(value,t_array) or \
            isinstance(value,z_array) or \
            isinstance(value,tz_array):  
        R[name] = value[start:end,:]
      else:
        R[name] = value

    # create a new object of the same type as D
    
    dtype = type(D)
    Robj = dtype()
    # transplant the attributes to the new object, which is technically sneaky
    vars(Robj).update(R)        
    return Robj

#def trim_altitudes(D, min_alt, max_alt):
def  trim_altitudes(D,min_alt,max_alt):
    """trim all the fields in a dictionary to [min_alt>=D.altitudes<=min_alt]"""
    
    # create a blank dictionary to return
    R = { } 

    #find bot_index and top_index
    [nshots,nalts]=D.msl_altitudes.shape
    indices=np.arange(nalts)
    alt_mask=D.msl_altitudes[0,:]>=min_alt
    if all(alt_mask==0):
      print 'Trim_altitudes----No altitudes in range'  
      Robj=[]
      return
    bottom=min(indices[alt_mask])
    alt_mask=D.msl_altitudes[0,:]<=max_alt
    #test for empty array
    if all(alt_mask==0):
      print 'Trim_altitudes---No altitudes in range'  
      Robj=[]
      return  
    top=max(indices[alt_mask])+1
    
    # iterate through all the items in a dictionary
    for (name, value) in vars(D).items():
      #only trim arrays with altitude components  
      if isinstance(value,z_array) or isinstance(value,tz_array):
        R[name] = value[:,bottom:top]
      else:
        R[name] = value

    # create a new object of the same type as D
    dtype = type(D)
    Robj = dtype()
    # transplant the attributes to the new object, which is technically sneaky
    vars(Robj).update(R)
            
    return Robj

