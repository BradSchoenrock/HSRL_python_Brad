
from lg_dpl_toolbox.filters.dpl_rolling_window_filter import WindowedFilterDescription
from lg_dpl_toolbox.filters.FIR_streaming_filter import FIR_streaming_filter
import numpy as np
from datetime import timedelta

#from scipy.signal import savgol_coeffs


#def myarrayfilter(array_value,some_unnamed_value,instrument,processing_defaults):
#    print array_value.shape
#    return array_value[array_value.shape[0]/2,:]




def filter_setup(instrument,processing_defaults): 
    if (instrument == 'magkazrge' or instrument == 'magmwacr') \
               and processing_defaults.enabled('ship_motion_correction'):
        from lg_dpl_toolbox.filters.FIR_streaming_filter import FIR_streaming_filter
        filter_type = 'hi_pass_rectangular'
        filter_width_sec = processing_defaults.get_value('ship_motion_correction','filter_length')
        order = 3
        filter_width = timedelta(seconds=filter_width_sec)
   
        return WindowedFilterDescription(
            filter=FIR_streaming_filter(filter_type,order)
           ,width=3, time_width=filter_width,edgemode='fullduplicate'
           ,varlist=['vertically_averaged_doppler'])
    else:
       return None

  
"""

#requires a full frame
class myframefilterobj(object):
    def __init__(self,myname,processing_defaults):
        self.framesSeen=0
        self.processing_defaults=processing_defaults
        self.name=myname

    def __del__(self):
        print 'Filter',self.name,'saw',self.framesSeen,'frames'

    def __call__(self,frames):
        print len(frames),[x.times[0] for x in frames]
        self.framesSeen+=1
        return frames[len(frames)/2] #just return the middle one

class SG_filter_object(object):
    def __init__(self,width,order):
        self.coef=someinitthing(window,order,use='dot')
        self.matr=None

    def __call__(self,array):
        if self.matr==None:
            self.matr=self.coef*ones(array.shape[1]) #FIXME

        #apply the correction


        return arr#the resulting 1xR vector


def select_filter(instrument,processing_defaults,timeres):
    #returns are, in order:
    # callable_obj: callable filter object or function, or None for no-op
    # width: an integer number of the preferred window width
    # edgemode: 'fullduplicate' or 'short'
    #       - short         : the edge cases will  have as few as half the number of samples. 
    #       - fullduplicate : duplicate end points to pad the window to fullness
    # varlist: list of variable names to run the filter on. each variable will be passed to the filter in its own call
    #     If None, will pass the entire list of frames to the filter in one single call
    # cargs: list of additional parameters, in order, to the filter
    # kwcargs: keyed dictionary of additional paramters to the filter
    #
    # if the return is:
    #  callable_obj,5,'fullduplicate',['var1','var2','var3'],[12345],dict(processing_defaults=processing_defaults,instrument=instrument)
    #resulting calls will be equivalent to:
    # f.var1[2,:]=callable_obj(f.var1[0:5,:],12345,processing_defaults=processing_defaults,instrument=instrument)
    # f.var2[2,:]=callable_obj(f.var2[0:5,:],12345,processing_defaults=processing_defaults,instrument=instrument)
    # f.var3[2,:]=callable_obj(f.var3[0:5,:],12345,processing_defaults=processing_defaults,instrument=instrument)

    #if False and processing_defaults.get_value('time_filter','filter_type')=='filterx':
        #from somewhere.core import myfilterx
        # returning the function defined above with a stored width value, to run on var1 and backscatter
        # with the extra in-order parameter 3.14159, and the extra named parameters instrument and processing_defaults set.

   
     
    if processing_defaults.enabled('ship_motion_correction'):
        #import hsrl.filters.FIR_rolling as FIR_rolling
        filter_length = processing_defaults.get_value('ship_motion_correction','filter_length')
        #convert filter_length from seconds to filter_width_bins
        delta_t=timeres.total_seconds()
        print 'delta_t', delta_t
        if delta_t >0:
            filter_width_bins = int(filter_length/delta_t)
        else:
            print
            print 'WARNING--delta_t = zero,   ship_motion_correction()'
            print
            return None
        
        #filter_bins must be odd
        if (filter_width_bins/2)*2 == filter_width_bins:
            filter_width_bins = filter_width_bins +1
        #if window is long enough compute savitzky-golay coeffs    
        if filter_width_bins >3:
           filter_order =3
           #fir_coeffs=savgol_coeffs(filter_wdith_bins,filter_order,use='dot')
           fir_coeffs=np.ones(5)
        #otherwise define coeffs to return the mean over the window width
        elif filter_width_bins == 3:
           fir_coeffs = np.ones(3)/3.0
        elif filter_width_bins == 2:
           fir_coeffs = np.ones(2)/2.0
        else:
           return None
       
        print 'fir_coeffs',fir_coeffs
        print j

        return WindowedFilterDescription(filter=FIR_rolling,width=filter_width_bins
                ,edgemode='fullduplicate',varlist=['var1','Backscatter'])
    
        ======      
        return WindowedFilterDescription(filter=myarrayfilter,width=processing_defaults.get_value('time_filter','width'),
                edgemode='fullduplicate',varlist=['var1','Backscatter'],
                cargs=[3.14159],kwcargs=dict(processing_defaults=processing_defaults,instrument=instrument))
        =====
        
    if False:# and processing_defaults.get_value('time_filter','filter_type')=='filterx':
        # returning an initialized callable object defined above with a width value 5, to run on the entire frame
        # with no addtional args, and the extra named parameter processing_defaults set.
        return WindowedFilterDescription(filter=myframefilterobj('dummyfilter',processing_defaults=processing_defaults),
            width=5,edgemode='fullduplicate',varlist=None,cargs=None,kwcargs=None)
 
    if False:#SG function
        width=xxxxxx #processing_defaults.get_value('time_filter','width')
        order=xxxxx #processing_defaults.get_value('time_filter','order'),
        return WindowedFilterDescription(filter=SG_filter_object(width,order),
            width=width,edgemode='fullduplicate',varlist=['Backscatter'],cargs=None,kwcargs=None)


    #fallback no-op
    return None

    """

