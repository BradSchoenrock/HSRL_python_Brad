import os
import json
import plistlib
from datetime import datetime,timedelta
import lg_base.core.open_config as oc

# def get_path_to_data(instrument,time, pref_file='/usr/local/etc/hsrl_python.json'):
global Data_path_dict
Data_path_dict = {}

def get_paths_to_data(instrument=None):
    global Data_path_dict
    if len(Data_path_dict)==0 or \
        (instrument is not None and not Data_path_dict.has_key(instrument)):

        fd = oc.open_config('hsrl_python.json')
        Data_path_dict = json.load(fd)['data_dir']
        fd.close()
    for k in Data_path_dict.keys():
        if instrument is not None and instrument.lower()!=k.lower():
            continue
        inst={'Name':k,'Path':Data_path_dict[k],'Instruments':[k]}
        yield inst
    return
 
def get_path_to_data(instrument, time=None):
    if isinstance(instrument,basestring) and instrument[0]=='/' and os.path.exists(instrument):
        return instrument
    for d in get_paths_to_data(instrument):
        return d['Path']
    raise RuntimeError('Instrument '+instrument+' not in configuration')

def timeintersection(startt,endt,windows):
    ret=[]
    if startt>endt:
        return ret
    if windows==None or len(windows)==0:
        tmp={'Start':startt,'End':endt}
        ret.append(tmp)
        return ret
    for i in windows:
        if i['Start']>=endt or ('End' in i and i['End']<=startt):
            continue
        tmp={'Start':startt,'End':endt}
        if startt<i['Start']:
            tmp['Start']=i['Start']
        if 'End' in i and endt>i['End']:
            tmp['End']=i['End']
        ret.append(tmp)
    return ret

class DataArchive(object):
    def __init__(self,dataarchivepath=None):
        self.filename=dataarchivepath or os.getenv('HSRL_DATA_ARCHIVE_CONFIG',"/etc/dataarchive.plist")
        self.archive=plistlib.readPlist(self.filename)


    def instrument(self,name):
        """instrument info
            instrument      - get info about individual instrument streams (always from dataarchive)
                metaframe takes shape of:
                    type
                    datasets
                    imagesets (prefix)
                    thumbsets
                        prefix
                        name
        """

        ins=self.archive['Instruments']
        k=ins.keys()
        kl=[x.lower() for x in k]
        nl=name.lower()
        return ins[k[kl.index(nl)]]

    def filterTimes(self,gen,start=None,end=None):
        if start==None:
            start=datetime(1990,1,1,0,0,0)
        if end==None:
            end=datetime(2100,1,1,0,0,0)
        for site in gen:
            assert('Windows' in site)
            w=timeintersection(start,end,site['Windows'])
            if len(w)==0:
                continue
            site['Windows']=w
            yield site

    def filterActive(self,gen,isactive):#always use before filterTimes
        for site in gen:
            assert('Windows' in site)
            hadActive=False
            for w in site['Windows']:
                if 'End' not in w:
                    hadActive=True
            if isactive!=hadActive:
                continue
            yield site

    def filterHidden(self,gen,hidden):
        for site in gen:
            ishid=False if 'hidden' not in site else site['hidden']
            if hidden!=ishid:
                continue
            yield site

    def filterReverse(self,gen):
        g=[x for x in gen]
        g.reverse()
        for site in g:
            yield site

    def latestTime(self,site):
        assert('Windows' in site)
        newest=datetime(1970,1,1,0,0,0)
        for window in site['Windows']:
            if 'End' in window:
                if window['End']>newest:
                    newest=window['End']
            elif 'Start' in window:
                if window['Start']>newest:
                    newest=window['Start']
        return newest

    def filterOrderTimes(self,gen):
        g=[x for x in gen]
        g.sort(key=self.latestTime)
        for site in g:
            yield site

    def _prepDataset(self,dsetidx,dset=None):
        dset=(dset or self.archive['Datasets'][dsetidx]).copy()
        dset['DatasetID']=dsetidx
        if "Instrument" in dset:
            dset['Instruments']=[dset['Instrument']]
        else:
            dset['Instruments']=[dset['Name']]
        return dset

    def dataset(self,dataset=None):
        if dataset is not None:
            try:
                dsi=int(dataset)
                dset=self._prepDataset(dsi)
                yield dset
                return
            except ValueError:
                pass
        for dsetidx,dset in enumerate(self.archive['Datasets']):
            if dataset!=None and dataset.lower()!=dset['Name'].lower():
                    continue
            dset=self._prepDataset(dsetidx,dset.copy())
            yield dset

    def _prepSite(self,siteid,site=None):
        site=(site or self.archive['Sites'][siteid]).copy()
        site['SiteID']=siteid
        return site

    def site(self,site=None):
        if site is not None:
            try:
                siteid=int(site)
                site=self._prepSite(siteid)
                yield site
                return
            except ValueError:
                pass
        sitename=site
        for siteidx,site in enumerate(self.archive['Sites']):
            if sitename is not None and site['Name']!=sitename:
                continue
            site=self._prepSite(siteidx,site)
            yield site

def exhaustIfNotNone(ret,parm):
    if parm is None:
        return ret
    newret=None
    for x in ret:
        if newret is not None:
            raise RuntimeError('too many inputs matching '+parm)
        newret=x
    if newret is None:
        raise RuntimeError('no inputs matching '+parm)
    return newret

def selectSource(*args,**kwargs):
    if len(args)>1:
        raise RuntimeError('Too many unnamed arguments: %i' % len(args))
    if len(args)==1:
        assert('instrument' not in kwargs)
        kwargs['instrument']=args[0]
    ret=None
    if 'instrument' in kwargs:
        assert(ret is None)
        inst=kwargs.pop('instrument')
        ret=get_paths_to_data(inst)
        ret=exhaustIfNotNone(ret,inst)
    if 'site' in kwargs:
        assert(ret is None)
        site=kwargs.pop('site')
        ret=DataArchive(kwargs.pop('dataarchive_path',None)).site(site)
        ret=exhaustIfNotNone(ret,site)
    if 'dataset' in kwargs:
        assert(ret is None)
        dataset=kwargs.pop('dataset')
        ret=DataArchive(kwargs.pop('dataarchive_path',None)).dataset(dataset)
        ret=exhaustIfNotNone(ret,dataset)

    if len(kwargs)>0:
        raise RuntimeError('Unused parameters exist for selectSource: '+(','.join([x for x in kwargs.keys()])))
    if ret is None:
        raise RuntimeError('Failed to determine source')

    return ret
