from collections import OrderedDict
import traceback

class removeGraphParams(object):
    def __init__(self,callnext=None):
        self.callnext=callnext

    def __call__(self,stream,args,kwargs):

        graphmode=False
        try:
            import lg_dpl_toolbox.dpl.dpl_dispatch_graph as ddg
            graphmode=isinstance(stream,ddg.DPLPrimativeDispatchGraph)
        except ImportError:
            pass
        if not graphmode:
            for x in ['templatesubscopes']:
                if x in kwargs:
                    del kwargs[x]

        if self.callnext!=None:
            self.callnext(stream,args,kwargs)

def get_hsrl_parameters(consts):
    attrs=OrderedDict()
    attrs['hsrl_wavelength_nm']=consts['wavelength']
    if not 'installation' in consts or consts['installation'] == 'ground':
        attrs['hsrl_altitude_m']=consts['lidar_altitude']
        attrs['hsrl_latitude_degN']=consts['latitude']
        attrs['hsrl_longitude_degE']=consts['longitude']
    return attrs

def getdictionary(obj):
    return obj.get_dict()

class addConstantsToParms(object):

    def __init__(self,callnext=None):
        self.callnext=callnext
        self.graphmode=None
        self.omitkeys=('doc','docs','documentation','parameters','Parameters')
        self.additional_attributes=OrderedDict()
        self.additional_attributes['hsrl_instrument']=None#('hsrl_instrument',None)
        self.additional_attributes['hsrl_constants_first']=('',get_hsrl_parameters)
        self.additional_attributes['radarType']=('radar_type',None)
        self.additional_attributes['radarLambda']=('radar_wavelength_m',None)
        self.additional_attributes['hsrl_process_control']=('hsrl_processing_parameter',getdictionary)
        self.additional_attributes['radar_parameters']=('radar_parameters',getdictionary)
        self.additional_attributes['spheroid_particle_parameters']=None
        self.additional_attributes['mass_dimension_particle_parameters']=None
        self.additional_attributes['multiple_scattering_parameters']=None


    def assembleAttributeDictionary(self,prefix,content):
        if isinstance(content,dict):
            ret=OrderedDict()
            keys=[k for k in content.keys()]
            keys.sort()
            for k in keys:
                if k in self.omitkeys or k.startswith('#'):
                    continue
                v=content[k]
                ret.update(self.assembleAttributeDictionary(((prefix+'__') if len(prefix)>0 else '')+k,v))
            return ret
        else:
            return {prefix:content}

    def streamWith(self,stream,attrname):
        if self.graphmode:
            for sname in stream.endNames:
                n=stream[sname]
                if hasattr(n,attrname):
                    return n
            return None
        if hasattr(stream,attrname):
            return stream
        return None

    def getAnAttribute(self,stream,attrname):
        s=self.streamWith(stream,attrname)
        if s==None:
            return None
        return getattr(s,attrname)

    def __call__(self,stream,args,kwargs):
        attrs=OrderedDict()
        self.graphmode=False
        try:
            import lg_dpl_toolbox.dpl.dpl_dispatch_graph as ddg
            self.graphmode=isinstance(stream,ddg.DPLPrimativeDispatchGraph)
        except ImportError:
            pass

        for a,k in self.additional_attributes.items():
            attr=self.getAnAttribute(stream,a)
            if attr!=None:
                attrs.update(self.assembleAttributeDictionary(k[0] if k is not None and k[0] is not None else a,k[1](attr) if k is not None and k[1] is not None else attr))

        if 'addAttributes' in kwargs:
            kwargs['addAttributes'].update(attrs)
        else:
            kwargs['addAttributes']=attrs
        if self.callnext!=None:
            self.callnext(stream,args,kwargs)

class addCalibrationsToNetCDF(object):
    def __init__(self,callnext=None):
        self.callnext=callnext
        self.output=None

    def replaceOutput(self,output):
        if self.output!=None:
            self.output.close()
        self.output=output

    def __call__(self,stream,args,kwargs):
        import lg_dpl_toolbox.dpl.dpl_artists as artists
        graphmode=False
        # begin specific hack to get separate cal stream into file if available
        if 'template' in kwargs:
            fmt,cfradial=artists.datasetParametersFor(kwargs['template'])
        elif 'templateparameters' in kwargs:
            fmt,cfradial=artists.datasetParametersFor(kwargs['templateparameters'])
        cals=[] 
        sn=None
        provides=None
        if hasattr(stream,'hsrl_cal_stream'):
            provides=stream.hsrl_cal_stream.provides
            for cal in stream.hsrl_cal_stream:#only applies to non-graph mode
                cals.append(cal)
            sn='hsrl_calibrations'
        else:
            graphmode=False
            try:
                import lg_dpl_toolbox.dpl.dpl_dispatch_graph as ddg
                graphmode=isinstance(stream,ddg.DPLPrimativeDispatchGraph)
            except ImportError:
                pass
            if graphmode and '3' in fmt:
                    cstream=None
                    print 'looking for stream'
                    sn=None
                    print 'possible cals',stream.endNames
                    for s in ('hsrl_calibrations','hsrl_constants'):
                        if cstream is None and s in stream:
                            print 'FOUND ',s
                            cstream=stream[s]
                            sn=s
                            break
                    if cstream is not None and True:#FIXME this assumes we can run it here
                        provides=cstream.provides
                        if 'force_dimsize' not in kwargs:
                            kwargs['force_dimsize']={}
                        kwargs['force_dimsize']['calibration']=0
                        count=0
                        for f in cstream:
                            count+=1
                        print 'THERE ARE ',count,'CALIBRATIONS from ',sn
                        kwargs['force_dimsize']['calibration']=count
        if len(cals)>0:
            if 'output' in kwargs:
                output=kwargs['output']
            else:
                from netCDF4 import Dataset
                try:
                    output=Dataset(kwargs['outputfilename'],'w',clobber=True,format=fmt)
                except IOError:
                    print 'Error on filename '+kwargs['outputfilename']
                    traceback.print_exc()
                    raise
                del kwargs['outputfilename']
                kwargs['output']=output
                self.replaceOutput(output)
            mykwargs=kwargs.copy()
            if 'template' in mykwargs and mykwargs['template'].endswith('.json') and sn is not None:
                import json
                d=json.load(file(mykwargs['template']), object_pairs_hook=OrderedDict)
                templates=[]
                from lg_base.core.locate_file import locate_file
                for streamconf in d['templates']:
                    if streamconf['substream']==sn:
                        tl=streamconf['template']
                        if not isinstance(tl,(list,tuple)):
                            tl=[tl]
                        if 'forModule' in streamconf:
                            mykwargs['forModule']=streamconf['forModule']
                        for t in tl:
                            templates.append(t)
                mykwargs['template']=templates
            import lg_dpl_toolbox.filters.substruct as frame_substruct
            #print stream.hsrl_cal_stream.provides
            #sb=kwargs['selected_bindings']
            #del kwargs['selected_bindings']
            try:
                nar=frame_substruct.TupleNarrator(cals,stream.hsrl_cal_stream.provides)
            except TypeError:
                nar=frame_substruct.TupleNarrator(cals)
            #formodule=[type(stream.hsrl_cal_stream).__module__,artists]
            #tmpkwargs=kwargs
            #if 'forModule' in kwargs:
            #    tmpkwargs=tmpkwargs.copy()
            #    formodule.extend(tmpkwargs.pop('forModule'))
            if 'basetime' not in mykwargs:
                mykwargs['basetime']=cals[0]['chunk_start_time']
                kwargs['basetime']=cals[0]['chunk_start_time']
            if len(mykwargs['template'])>0:
                for f in artists.dpl_netcdf_artist(nar,usecfradial=cfradial, withUnlimited=len(cals) if '3' in fmt else None,*args,**mykwargs):
                    kwargs.pop('basetime',None)
            #kwargs['selected_bindings']=sb
            output.sync()
        if self.callnext!=None:
            self.callnext(stream,args,kwargs)

    def __del__(self):
        self.replaceOutput(None)

