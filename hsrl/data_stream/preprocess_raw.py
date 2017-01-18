import numpy as np
import lg_base.core.array_utils as hau
import copy

"""
 filterobj = preprocess_raw('ahsrl',ahsrlconstants,ntime_ave)"""

class time_frame:
    def __init__(self,post_operator=None):
        self.post_operator=post_operator

    def cleanup(self,raw):
        if raw.times.shape[0] ==0:
            raw.start_time = None
        else:    
            raw.start_time = raw.times[0]
            
        if raw.times.shape[0]>0:
            raw.end_time = raw.times[-1]
        else:    
            raw.end_time = None        

    def flush(self,raw=None):
        if raw is not None:
            self.cleanup(raw)
        if self.post_operator!=None:
            return self.post_operator.flush(raw)
        return raw



    def __call__(self,raw):
        #if these lines ever crash, no data was appended to no data. why would this run?     
        self.cleanup(raw)

        if self.post_operator!=None:
            return self.post_operator(raw)
        return raw



class time_select:
    def __init__(self,ntime_ave,post_operator=None):
        self.ntime_ave = ntime_ave
        self.post_operator=post_operator
        self.Rem=hau.Time_Z_Group()

    def flush(self,raw=None):
        if raw is not None:
            self.Rem.append(raw)
        raw=self.Rem
        self.Rem=hau.Time_Z_Group(like=raw)
        if self.post_operator!=None:
            return self.post_operator.flush(raw)
        return raw
 
    def selectTimesMod_kt(self,raw):
        """ selectTimesMod(self,kt):     
            Iterate thought all items in dictionary, add any remainder from
            last call then return vars(self) with number of times evenly divided
            by kt and a new remainder.
        """
        kt=self.ntime_ave
        #entries in self.Rem.times but not in raw.times, copy Rem to Raw
        if hasattr(self.Rem,'times') and( not hasattr(raw,'times') or raw.times.shape[0] ==0):
            #for (name,value) in vars(self.Rem).items():
            #    raw[name]=self.Rem[name]
            #    self.Rem[name]=[]
             pass   
        #note--if entries in raw.times but not in self.Rem.times
        else:
          if hasattr(raw,'times') and (not hasattr(self.Rem,'times') or self.Rem.times.shape[0] == 0):
            arrlen=raw.times.shape[0]
            self.Rem=copy.copy(raw)#copy.deepcopy(raw)
        #both Rem.times and raw.times have entries, concatenate Rem and raw
        #and define new remainder if it exists
          else:
            arrlen=raw.times.shape[0] + vars(self.Rem)['times'].shape[0]
            self.Rem.append(raw)

          last_index = arrlen-(arrlen%kt)
          mask =np.arange(arrlen)
          #mzero=mask<0
          #mrem= mask>=last_index
          #mret= mask<last_index
                   
          for (name,value) in vars(self.Rem).items():
                #print 'empty Rem--mask.size,last_index,raw.times.size',mask.size,last_index,raw.times.size
                if hau.safeIsInstance(value,hau.T_Array) or hau.safeIsInstance(value,hau.TZ_Array):
                  if last_index<arrlen:
                    setattr(self.Rem,name,value[last_index:])
                  else:
                    l=list(value.shape)
                    l[0]=0
                    setattr(self.Rem,name,type(value)(np.ones(l),summode=value.summode,dtype=value.dtype))
                  if last_index>0:
                    setattr(raw,name,value[:last_index])                                             
                  else:
                    l=list(value.shape)
                    l[0]=0
                    setattr(raw,name,type(value)(np.ones(l),summode=value.summode,dtype=value.dtype))
                else:
                  setattr(raw,name,copy.deepcopy(value))
        #if these lines ever crash, no data was appended to no data. why would this run?     
        if raw.times.shape[0] ==0:
            raw.start_time = None
        else:    
            raw.start_time = raw.times[0]
            
        if raw.times.shape[0]>0:
            raw.end_time = raw.times[-1]
        else:    
            raw.end_time = None
        if self.Rem.times.shape[0] ==0:
            self.Rem.start_time = None
        else:    
            self.Rem.start_time = self.Rem.times[0]
            
        if self.Rem.times.shape[0]>0:
            self.Rem.end_time = self.Rem.times[-1]
        else:    
            self.Rem.end_time = None
        
         
        
        
      
       
       
 
    def __call__(self,raw):
        """ appends the contents of raw to any profiles stored in rem.
        It then selects first int(ns)hots/ntime_ave) profiles in the resulting
        buffer and does any processing that must be completed before time
        averaging on these profiles. It then time averages in ntime_ave
        blocks. The remaining profiles are held in rem to be appended on the
        next call."""
       
        self.selectTimesMod_kt(raw)
        
        if self.post_operator!=None:
            return self.post_operator(raw)
        return raw


class time_ave:
    def __init__(self,ntime_ave,post_operator=None):
        self.ntime_ave = ntime_ave
        self.post_operator=post_operator

    def work(self,raw,ntime_ave):
        raw.binMean(ntime_ave,1,complete_only=False)
        raw.preprocess_ntime_ave = ntime_ave

    def flush(self,raw=None):
        if raw is not None and raw.times.size>0:
            self.work(raw,raw.times.size)

        if self.post_operator!=None:
            return self.post_operator.flush(raw)
        return raw

    def __call__(self,raw):

        self.work(raw,self.ntime_ave)

        if self.post_operator!=None:
            return self.post_operator(raw)
        return raw
