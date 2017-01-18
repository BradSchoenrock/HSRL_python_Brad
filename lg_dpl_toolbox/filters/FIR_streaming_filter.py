import numpy as np
import lg_base.core.array_utils as hau

class FIR_streaming_filter(object):
    "Finite Impulse Filter applied to a rolling window"
    def __init__(self,filter_type,order = None):
        """__init__(self,filter_type,width,time_res,order = None)
           filter_type = 'savitzky_golay | 'mean'
           order       = order of filter (if required)"""
        
        self._coeffs_matrix1d = {}
        self._coeffs_matrix2d = {}
        self._fir_coeffs = {}
        self.filter_type=filter_type
        self.order=order

    @property
    def coeffs_matrix1d(self):
        if self.usewidth in self._coeffs_matrix1d:
          return self._coeffs_matrix1d[self.usewidth]
        return None
        
    @coeffs_matrix1d.setter
    def coeffs_matrix1d(self,val):
        self._coeffs_matrix1d[self.usewidth]=val
        
    @property
    def coeffs_matrix2d(self):
        if self.usewidth in self._coeffs_matrix2d:
          return self._coeffs_matrix2d[self.usewidth]
        return None
        
    @coeffs_matrix2d.setter
    def coeffs_matrix2d(self,val):
        self._coeffs_matrix2d[self.usewidth]=val

    @property
    def fir_coeffs(self):
        if self.usewidth in self._fir_coeffs:
          return self._fir_coeffs[self.usewidth]
        return None
        
    @fir_coeffs.setter
    def fir_coeffs(self,val):
        self._fir_coeffs[self.usewidth]=val

    def loadcoefsFor(self,width):
        self.usewidth=width
        if self.fir_coeffs!=None:
          return
        print 'Making coefficients for width ',width,' with type ',self.filter_type

        if self.filter_type == 'mean':
            self.fir_coeffs = 1.0/width * np.ones(width)
            

        elif self.filter_type == 'hi_pass_rectangular':
           #rectangular window high pass filter 
           self.fir_coeffs = -1.0/(width-1)*np.ones(width) 
           center_bin = width/2
           self.fir_coeffs[center_bin] = 1.0
            
        elif  self.filter_type == 'savitzy_golay':
           #if window is long enough compute savitzky-golay coeffs    
           if width >3:
               filter_order =3
               self.fir_coeffs=savgol_coeffs(width,filter_order,use='dot')
               self.fir_coeffs=np.ones(5)
           #otherwise define coeffs to return the mean over the window width
           elif width == 3:
              self.fir_coeffs = np.ones(3)/3.0
           elif width == 2:
              self.fir_coeffs = np.ones(2)/2.0
           else:
              return None
           import scipy.signal.savgol_coeffs as savgol_coeffs
           self.fir_coeffs = savgol_coeffs(width,self.order,use='dot')
        else:
            print "WARNING--FIR_rolling---" + type + " is an unknown filter type"
            return None
        
    def __call__(self,array):
        """Implelments a general Finite Impulse Filter, providing
           filtered values for the time mid-point of the supplied
           [time,alt] data array. This is the convolution of FIR
           coeffs and data window----it multiplies coef_matrix[time,alt]
           *  data_window[times,alt],sums in the time dimension and
           returns an altitude vector representing the filtered value
           at the mid-point of the time window"""
        self.loadcoefsFor(array.shape[0])
         
        if isinstance(array,hau.TZ_Array):
            if self.coeffs_matrix2d == None:
                self.coeffs_matrix2d = self.fir_coeffs[:,np.newaxis]**np.ones(array.shape[1])
            
        else:  #for T_array, which is a superclass of TZ_Array
            if self.coeffs_matrix1d == None:
                self.coeffs_matrix1d = self.fir_coeffs
        
        if isinstance(array,hau.TZ_Array):
            ret=np.sum(self.coeffs_matrix2d * array,0) 
            #ret=ret.reshape([1]+list(ret.shape))#FIXME
        else:
            ret=np.sum(self.coeffs_matrix1d * array,0)
        ret=ret.reshape([1]+list(ret.shape))
        
        return ret
