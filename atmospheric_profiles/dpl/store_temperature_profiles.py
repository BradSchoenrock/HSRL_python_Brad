#parameters: instrument, filename, start, optional end (implies month)

def main():
    from dpl_temperature_profiles import dpl_virtualradiosonde
    import sys,os
    from datetime import datetime,timedelta
    import numpy
    from lg_dpl_toolbox.dpl.NetCDFZookeeper import GenericTemplateRemapNetCDFZookeeper
    import hsrl.dpl.HSRLLibrarian as hsrllib
    from lg_dpl_toolbox.filters.time_frame import TimeGinsu
    import lg_dpl_toolbox.dpl.dpl_artists as dpl_artists
    datatype=sys.argv[1]
    filename=sys.argv[2]
    deltat=float(sys.argv[3])
    st=datetime.strptime(sys.argv[4],'%Y%m%dT%H%M%S')
    if len(sys.argv)<6:
        import calendar
        monthrange=calendar.monthrange(st.year,st.month)
        monthdur=timedelta(days=monthrange[1])
        tmp=st+monthdur
        et=datetime(tmp.year,tmp.month,1,0,0,0)
        if et>datetime.utcnow():
            et=datetime.utcnow()
    else:   
        et=datetime.strptime(sys.argv[5],'%Y%m%dT%H%M%S')
    #fields=['times','telescope_position','telescope_rotation','telescope_rotation_measured','telescope_elevation','telescope_accelerometer_raw']#,'superseedlasercontrollog','laserpowervalues']
    zoo=GenericTemplateRemapNetCDFZookeeper(datatype,user_read_mode='position',forModule=hsrllib)#,keepfields=fields)
    lib=hsrllib.HSRLLibrarian(instrument=datatype,zoo=zoo)#site=16)#None,datatype)
    m=lib(start=st,end=et)#,filetype='data')
    m=TimeGinsu(m,'times',None)
    vr_lib=dpl_virtualradiosonde('name',os.getenv('GRIB_CACHE','/arcueid/data/grib_cache'),timedelta(minutes=deltat),numpy.arange(0,30000+.1,15),do_interpolate=False)
    m=vr_lib(m,expire_duration=timedelta(minutes=deltat))
    art=dpl_artists.dpl_netcdf_artist(m,'NWS_Profile_Archive.cdl',filename,format='NETCDF4',usecfradial=False)
    #art()
    for f in art:
        print vars(f)
        print f.temps.shape,f.temps



if __name__ == '__main__':
    main()
