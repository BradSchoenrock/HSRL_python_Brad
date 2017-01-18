from collections import OrderedDict
import re
import lg_base.core.json_config as jc
import lg_base.core.locate_file as lf

def filterDictionary(*args,**kwargs):
    """call filterDictionaryGenerator to create a dictionary"""
    return OrderedDict([(k,v) for k,v in filterDictionaryGenerator(*args,**kwargs)])

def filterDictionaryGenerator(d,matching=None,keylist=None,keylist2=None):
    """ create a generator that iterates key/value pairs of a dictionary based on parameters
        
        :param d: source dictionary to iterate and filter
        :param matching: regular expression to match keys to
        :param keylist: list of acceptable keys or None
        :param keylist2: list of acceptable keys or None

        any key must pass each individual test before it can be yielded.
        If there is no overlap between the lists or the regex, nothing will be returned.
    """

    assert(isinstance(d,dict))
    match=None
    if matching is not None:
        match=re.compile(matching)
    
    for k,v in d.items():
        if match is not None and not match.match(k):
            continue
        if keylist is not None and k not in keylist:
            continue
        if keylist2 is not None and k not in keylist2:
            continue
        yield k,v

def mergedDictionaryCopy(*args):
    """ merge a series of dictionaries using the update() function.
        the first parameter will be shallow-copied, and then update() with each progressive parameter
    """ 
    ret=args[0].copy()
    for x in args[1:]:
        ret.update(x)
    return ret

def rip_json_config_parameters(kwargs):
    """ remove json_config init parameters from dictionary (modifies parameter) and return them as a separate dictionary
    """
    ret=dict()
    for arg in ('allow_missing_values','allow_missing_enables'):
        if arg in kwargs:
            ret[arg]=kwargs.pop(arg)
    return ret 

def only_locate_file_parameters(kwargs):
    ret=kwargs.copy()
    ret.pop('enable_masking',None)
    return ret

def maybeMakeJC(displays,majorkey,*args,**kwargs):
    """ Attempt to make a json_config object from a display and majorkey (required subsection used as default)
        if the first parameter is a string, an attempt to locate it will be made.
           if this fails, the string will be returned as-is
           if this succeeds, the resulting json_config is loaded
        if its not a string, it will be treated like a dictionary for a json_config initialization
        if the loading above fails because the major key is missing, the located string or dictionary is returned as is
        if it fails for another reason, the exception as to why is raised
    """
    if majorkey is not None and not isinstance(displays,jc.json_config):
        jsd=rip_json_config_parameters(kwargs)
        try:
            if isinstance(displays,basestring):
                try:
                    displays=lf.locate_file(displays,**only_locate_file_parameters(kwargs))
                except IOError:
                    return displays
            #if its located or is not a string, try loading as a json_config
            return jc.json_config(displays,majorkey,**jsd)
        except KeyError:
            pass
        except ValueError:
            if displays.endswith('.json'):
                raise
    return displays


def jsonParameterListGenerator(displays,majorkey=None,*args,**kwargs):
    """ given a displays parameter, generate artist parameter set parameters

        if is a tiered json, will yield multiple entries, possibly with parameters for artist
        if is a list, same
        if is just a name, will yield a single entry

        returns a subset key (default if non-tiered), displays parameter (json_config or string if couldn't load), and artist parameter kwargs dictionary
    """
    jsd=rip_json_config_parameters(kwargs)
    if isinstance(displays,dict):
        if majorkey is not None and majorkey in displays:
            yield 'default',maybeMakeJC(displays,majorkey,**jsd),{}
            return
        for k, disp in filterDictionaryGenerator(displays,*args,**kwargs):
            if isinstance(disp,(list,tuple)):
                disp=OrderedDict([(d,{}) for d in disp])
            for filename,params in disp.items():
                params=params.copy()
                if not params.pop('enable',True):
                    continue
                if majorkey is not None and majorkey in params:
                    yield k,maybeMakeJC(params,majorkey,**jsd),{}
                else:
                    jsd.update(params)
                    yield k,maybeMakeJC(filename,majorkey,**jsd),params
    elif isinstance(displays,(tuple,list)):
        for filename in displays:
            yield 'default',maybeMakeJC(filename,majorkey,**jsd),{}
    else:
        yield 'default',maybeMakeJC(displays,majorkey,**jsd),{}

class jsonListAccumulator(object):
    """ container for display defaults, with funtion to store runtime value (for reuse or modification)
    """
    def __init__(self,jsonconfig,majorkey=None,**kwargs):
        self.conf=jsonconfig
        if isinstance(self.conf,basestring):
            import json
            import lg_base.core.open_config as oc
            self.conf=json.load(oc.open_config(self.conf), object_pairs_hook=OrderedDict)
        self.majorkey=majorkey
        self.caching=OrderedDict()
        self.sourcename={}
        self.parms={}
        self.extraargs=kwargs

    def __call__(self,*args,**kwargs):
        tmp=self.extraargs.copy()
        tmp.update(kwargs)
        for i,x in enumerate(jsonParameterListGenerator(self.conf,self.majorkey,*args,**tmp)):
            k,v,p=x
            if k not in self.caching:
                self.caching[k]=[]
                self.sourcename[k]=[]
                self.parms[k]=[]
            if i>=len(self.caching[k]):
                self.caching[k].append(v)
                if isinstance(v,jc.json_config):
                    self.sourcename[k].append(v.confg_file)
                else:
                    self.sourcename[k].append(v)
                self.parms[k].append(p)
                print 'new entry for '+k,i
            else:
                print 'reusing '+k,i
            self.sk=k
            self.si=i
            yield k,self.caching[k][i],self.parms[k][i]

    def storeStructure(self,stru):
        if isinstance(self.caching[self.sk][self.si],(dict,basestring)):
            self.caching[self.sk][self.si]=stru
        elif stru is not self.caching[self.sk][self.si]:
            print 'WARNING structure changed for',self.sk,self.si
            self.caching[self.sk][self.si]=stru

    def __iter__(self):
        for k,v in self.caching.items():
            for i,d in enumerate(v):
                if isinstance(self.sourcename[k][i],basestring):
                    yield k,self.sourcename[k][i],d
                else:
                    yield k,'',d
