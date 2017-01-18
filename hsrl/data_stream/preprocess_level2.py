import numpy as np
import lg_base.core.array_utils as hau

class preprocess_level2:
    def __init__(self,instrument,post_operator=None):
        self.instrument = instrument
        self.post_operator = post_operator
    
    def cleanup(self,raw):
        for (name,value) in vars(raw).items():
            if hau.safeIsInstance(value,hau.T_Array) \
                   and value.size >0 \
                   and  (type(value[0]) == 'numpy.float32' \
                         or type(value[0]) == 'numpy.float64'):             
                value[value>1e10]=np.NaN 

    def flush(self,raw):
        if raw is not None and raw.times.size>0:
            self.cleanup(raw)
        if self.post_operator:
            return self.post_operator.flush(raw)
        return raw

    def __call__(self,raw):
        """Use this to do any pre processing that is best done after
           the initial block ntimes_ave"""

        #replace any "out of range" values in T_Arrays by NaN's
        self.cleanup(raw)

        #for those variables that are best represented as sums,
        #multiply by ntimes_ave to compensate for pre averaging
        
        #if hasattr(raw,'seeded_shots'):
        #    raw.seeded_shots*=self.ntime_ave
        #if hasattr(raw,'shot_count'):
        #    raw.shot_count*=self.ntime_ave

        if self.post_operator:
            self.post_operator(raw)
        return raw

            
