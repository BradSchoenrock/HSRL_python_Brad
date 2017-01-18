import lg_base.core.array_utils as hau
import os
from datetime import datetime,timedelta
import lg_dpl_toolbox.formats.VectorTableLibrarian

def getBasedir(siteid):
        import lg_dpl_toolbox.core.archival as hru
        if isinstance(siteid,basestring):
            if os.access(siteid,os.F_OK):
                return siteid
            else:
                try:
                    tmp = hru.selectSource(instrument=siteid)
                except:
                    tmp = hru.selectSource(site=siteid)
        else:
            try:
                tmp = hru.selectSource(instrument=siteid)
            except:
                tmp = hru.selectSource(site=siteid)
        return tmp['Path']

_callist=None

def callist():
    global _callist
    if _callist is None:
        import hsrl.data_stream.hsrl_read_utilities as hru #OH HOW I HATE THIS
        _callist=dict(baseline=dict(prefix='baseline_correction',suffix='.blc',reader=hru.read_baseline),
                    geo=dict(prefix='geofile_default_',suffix='.geo',reader=hru.read_geo_corr),
                    d_geo=dict(prefix='diff_default_geofile',suffix='.geo',reader=hru.read_diff_geo),
                    d_geo_1064_532=dict(prefix='diff_1064_532_geofile',suffix='.geo',reader=hru.read_diff_1064_532_geo),       
                    i2a_d_geo=dict(prefix='i2a_mol_diff_geo',suffix='.geo',reader=hru.read_i2a_diff_geo),
                    n_geo=dict(prefix='nadir_default_geofile',suffix='.geo',reader=hru.read_nadir_geo_corr),
                    i2=dict(prefix='i2-default-scan',suffix='.cal',reader=hru.read_i2_scan),
                    rb=dict(prefix='rb-default-scan',suffix='.cal',reader=hru.read_i2_scan),#FIXME
                    qw_baseline=dict(prefix='qwave_baseline',suffix='.blc',reader=hru.read_qw_baseline),
                    pol=dict(prefix='pol_cal_default',suffix='.json',reader=hru.read_pol_cal),
                    i2a_temp_table=dict(prefix='i2a_temp_table',suffix='.cal',reader=hru.read_i2a_temp_table),
                    #cp_d_geo=dict(prefix='cross_pol_diff_geofile',suffix='.geo',reader=hru.read_cross_poll_diff_geo),
                    cpol_diff_geo=dict(prefix='cpol_diff_geo',suffix='.geo',reader=hru.read_cross_poll_diff_geo),
                    i2offset=dict(prefix='i2offsettable',suffix='.i2off',hasDate=False,reader=hru.read_i2_offset))
    return _callist

def calDataInfo(name):
    if name not in callist():
        r='Unknown cal file type '+name+'. Add it to callist in TableLibrarian.py'
        print r
        raise RuntimeError(r)
    return callist()[name].copy()

class TableLibrarian(lg_dpl_toolbox.formats.VectorTableLibrarian.VectorTableLibrarian):
    """ Librarian Initialzation for Calibration Tables
  
            :param siteid: source site id for the data source. typically an hsrl instrument or base directory
            :param datatype: file type to list
    """
    def __init__(self, siteid,datatype,**kwargs):
        super(TableLibrarian,self).__init__(instrumentname=siteid,datatype=calDataInfo(datatype),basedir=getBasedir(siteid),**kwargs)

def findFile(instrument,caltype,moment,*args,**kwargs):
    lib=TableLibrarian(instrument,caltype,*args,**kwargs)#,completeList=True)
    nar=lib(moment,moment+timedelta(seconds=1))#TableNarrator(getBasedir(datatype),datatype,caltype,caltype,st,et,completeList=True)
    for f in nar:
        return f
    return None

def findCalFile(*args,**kwargs):
    expire_time=kwargs.pop('expire_time',None)
    v=findFile(*args,**kwargs)
    if v is None:
        altp=kwargs.pop('alternativePath',None)
        r="Couldn't find cal file "+args[1]+' for '+args[0]+' time '+args[2].strftime('%Y.%m.%dT%H:%M:%S'+('' if altp is None else (' using alt path '+altp)))
        print r
        return (None, datetime(2200,1,1,0,0,0))
        #raise RuntimeError(r)
    return (v['path'],expire_time or (v['start']+v['width']))

if __name__=='__main__':
    import sys
    datatype='ahsrl' if len(sys.argv)<2 else sys.argv[1]
    caltype='geo'  if len(sys.argv)<3 else sys.argv[2]
    et=datetime.utcnow() if len(sys.argv)<5 else datetime.strptime(sys.argv[4],'%Y%m%dT%H%M%S')
    st=(et-timedelta(days=.5)) if len(sys.argv)<4 else datetime.strptime(sys.argv[3],'%Y%m%dT%H%M%S')
    lib=TableLibrarian(datatype,caltype)#,completeList=True)
    nar=lib(st,et)#TableNarrator(getBasedir(datatype),datatype,caltype,caltype,st,et,completeList=True)
    print 'first is ',findCalFile(datatype,caltype,st)
    for i,f in enumerate(nar):
        print i,f
        
