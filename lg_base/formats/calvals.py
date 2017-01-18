#read a calvals text file


#import string
import string as st
from datetime import datetime
import numpy as np
import collections
from lg_base.core.open_config import open_config
from lg_base.core.locate_file import locate_file

#create a calvals object
class calvals_class(object):
    def __init__(calvals,filename=None,instrument=None,systemOnly=True,varlist=None,level=0):
        if filename is None:
            filename='calvals_'+instrument+'.txt'
            filename=locate_file(filename,systemOnly=systemOnly,level=1+level)#level 1 means find based on caller
        print 'cal_file_reader: opening ', filename
        fileid = open(filename)
    # fileid = open_config(filename)
        import lg_base.core.read_utilities as hru
        name=None
        date_value=[]
        cal_vec=np.ndarray(3)
        realunits=None
        for line in fileid:
            if '#' in line:
                line=line[:line.find('#')]#strip any comment suffix
                if len(line)==0:
                    continue
            if line[0] != ' ':#header line
                #if this is not the first time through
                if name: 
                    if vars(calvals).has_key(name):
                        raise RuntimeError, 'Error - duplicate entry for %s in %s' % (name, filename)
                    setattr(calvals,name,date_value)
                date_value=[]
                index=line.find('(')
                index2=line.find(')') 
                if index >=0:
                    name=st.rstrip(st.lstrip(line[:index]))
                    units=st.rstrip(st.lstrip(line[index+1:index2]))
                    if len(st.rstrip(st.lstrip(line[index2+1:])))>0:
                        raise RuntimeError, "Format error. Variable has excess content after units. %s in %s" % (name,filename)
                else:
                    if index2>=0:
                        raise RuntimeError, "Format error. Variable names can't have parentheses. %s in %s" % (name,filename)
                    name=st.rstrip(st.lstrip(line)) 
                    units=None
                if varlist is not None and name not in varlist:
                    name=None
                realunits=units
                k=0                      
            elif len(st.lstrip(line))>0:#record line
                #look for required comma after effective date string
                units=realunits
                try:
                    index=line.index(',')
                except ValueError:
                    print ' '
                    print 'ERROR ***************missing comma in calvals file **********'
                    print line
                    print ' '
                    raise RuntimeError('Missing comma in calvals file on line %i' %(line))  

                e_date=hru.convert_date_str(line[:index],twodigityearok_fubar=True)['datetime']

                line=line[(index+1):]
                #check for string value as return
                if line.find("'")>=0 or line.find('"')>=0:
                    string_var=st.rstrip(st.lstrip(line))
                    string_var=st.replace(string_var,"'","")
                    string_var=st.replace(string_var,'"','')
                    date_value.append(list())
                    date_value[k].append((e_date,string_var,None))
                #check for calval vector   
                elif line.find('[')>=0 and line.find(']')>=0:
                    index=line.index('[')
                    index2=line.index(']')
                    tmp_vec=eval(line[(index+1):index2])               
                    line=line[(index2+1):]
                    if line[0]=='/':
                        divisor=float(eval(line[1:]))
                        cal_vec=np.array(tmp_vec)/divisor  
                    elif line[0]=='*':
                        multiplier=float(eval(line[1:]))
                        cal_vec=np.array(tmp_vec)*multiplier
                    elif isinstance(tmp_vec,tuple):
                        cal_vec=np.array(tmp_vec)
                    else: # special case for scalars - don't save as array
                        cal_vec = tmp_vec
                    date_value.append(list())
                    if units is not None and units == "deg W":
                        cal_vec=-cal_vec
                        units="deg E"
                    date_value[k].append((e_date,cal_vec,units))
                else:
                    print ' '
                    print 'Error*************calvals syntax error****************'
                    print line
                    print ' '
                    raise RuntimeError('Calvals syntax error on line %i' %(line))
                k=k+1
        if name is not None:
            if hasattr(calvals,name):
                raise RuntimeError, 'Error - duplicate entry for %s in %s' % (name, filename)
            setattr(calvals,name,date_value)
        # validate all entries contain descending times
        for (name, value) in vars(calvals).items():
            cur_time = datetime(2200, 1,1, 0,0,0)
            for i in range(len(value)):
                if cur_time < value[i][0][0]:
                    print 'cal_file_reader for %s: %s < %s' % (name, cur_time, value[i][0][0])
                    raise RuntimeError ,'Error - non-descending time values in %s for %s' % (filename, name)
                cur_time = value[i][0][0]

    def units(calvals,keyname):
        if not hasattr(self,keyname):
            return None
        return getattr(self.keyname)[0][2]

    def __call__(self, time):
        """select calibration constants at a given time
       from 'calvals_xxhsrl.txt'. Return current values
       on call and append any values that have changed to
       the archive.
       """

        return self.select_time(time)

    def select_time(calvals,time):
        rs_constants={}
        #set expiration date far in the future if none is specified   
        next_cal_time=datetime(2100,1,1,0,0,0)
        prior_cal_time=None
        rs_constants['next_cal_time']=next_cal_time
        rs_constants['prior_cal_time']=prior_cal_time
        for (name, value) in vars(calvals).items():
            if name==None or value==None:
                continue

            #loop through all of the calibration constants
            #find values for given time
            # if this calibration is later than current time skip to next
            cal_index=len(value)
            for i in range(len(value)):
                if time >= value[i][0][0] :
                    cal_index = i
                    break
            # if the specified time is less than the smallest calibration time, flag as an error,
            # since we need a calibration time less then the specified time
            if cal_index == len(value):
                print 'error on calvals',name,time,cal_index,len(value)
                raise RuntimeError, """ select_sys_constants: no calibration value for %s at time %s ;
                          please try a time after %s. If this is an error, update calvals
                          """ %  (name, str(time), str(vars(calvals)[name][-1][0][0]))          
            # we never found a smaller calibration time.  Use the 0th value (for the latest calibration time)
            # NOTE - there is no expiration time for this calibration
            elif cal_index == 0:
                exp_index = -1

            else:
                exp_index = cal_index-1
                assert(exp_index >= 0)
                # exp_index now contains the calibration time just after this time , when the current calibration time expires
                # always choose the earliest expiration time
                if value[exp_index][0][0]<next_cal_time:
                    next_cal_time=value[exp_index][0][0]
                    #print 'updating next_cal_time'
                    rs_constants['next_cal_time']=next_cal_time
            if (prior_cal_time is None or value[cal_index][0][0]>prior_cal_time) and value[cal_index][0][0]<time:
                    prior_cal_time=value[cal_index][0][0]
                    rs_constants['prior_cal_time']=prior_cal_time

            rs_constants[name]=(value[cal_index][0][1])
            rs_constants[name+'_units']=value[cal_index][0][2]

        return rs_constants

def cal_file_reader(instrument,level=1,**kwargs):
    """   
    cal_file_reader(instrument)
    reads calvals_xxxxx.txt file
    """

    return calvals_class(instrument=instrument,level=level,**kwargs)

def matching(a,b):
    return set(a).intersection(set(b))==set(a) and set(a).union(set(b))==set(a)

def lint(filename):
    print 'linting',filename
    c=calvals_class(filename=filename)
    t=datetime.utcnow()
    keys=None
    try:
        while True:
            x=c(t)
            if keys is not None:
                assert(matching(keys,x.keys()))
            keys=[y for y in x.keys()]
            t=x['prior_cal_time']
            if t is None:
                break
    except RuntimeError:
        print 'breaks at',t
    #print 'keys are ',x.keys()
    print 'Done'

#select the system constants for the current time
#rs_constants={}
def select_sys_constants(calvals,time):    
    """returns a dictionary, rs_constants, that provides the constant value

       calvals is object with members like 'calvals.cross_pol_dead_time', with an sequence of datetimes/value pairs
       	e.g:
       polarization_cross_talk
         20-apr-12, [0.0]
         24-feb-12, [0.015]
         10-jan-12, [0.002]
         12-oct-11, [0.0095]
         25-aug-11, [0.01]
         19-aug-11, [0.014]

       NOTE - datetimes are stored in DECREASING order
       (calvals could have been a dictionary, rather than a class with dynamically added members)


       as we add values, determine the calibration constant that will expire soonest

       """
    return calvals(time)

def main():
    import sys
    for f in sys.argv[1:]:
        lint(f)

if __name__ == '__main__':
    main()