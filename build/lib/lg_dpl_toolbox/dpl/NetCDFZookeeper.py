import dplkit.role.zookeeper
import json
import lg_base.core.open_config as oc
import lg_base.core.read_utilities as hru
import os
import struct
import bz2
from collections import OrderedDict

class FileSystemZookeeper(dplkit.role.zookeeper.aZookeeper):
    def __init__(self):
        super(FileSystemZookeeper,self).__init__()

    def obtain(self,uri, *args, **kwargs):
        return uri

class GenericTemplateRemapNetCDFZookeeper(dplkit.role.zookeeper.aZookeeper):
    """ GenericTemplateRemapNetCDFZookeeper - Zookeeper for reading raw files

        :param template_prefix: netcdf defaults json file prefix
        :param template: optional template json dictionary. loaded from file named above if not provided
        :param maxbin: maximum bin to read in
        :param preprocess: preprocessor object filter
        :param keepfields: if provided, lists fields to be kept explicitly
        :param verbose: if true, will be more verbose in its reading preprocess
        :param user_read_mode: if provided, overrides the template default read mode with a user provided read_mode value
    """
    def __init__(self,template_prefix,template=None,maxbin=50000,preprocess=None,keepfields=None,verbose=False,user_read_mode=None,forModule=None):
        super(GenericTemplateRemapNetCDFZookeeper,self).__init__()
        self.template_prefix=template_prefix
        self.read_raw_parms={}
        if verbose!=None:
            self.read_raw_parms['verbose']=verbose
        if user_read_mode!=None:
            self.read_raw_parms['user_read_mode']=user_read_mode
        self.template=template if template else json.load(oc.open_config(template_prefix+'_netcdf_defaults.json',level=1,forModule=forModule), object_pairs_hook=OrderedDict)
        if keepfields:
            tmp=self.template.copy()
            tmp['selected_vars']={}
            for f in keepfields:
                if f in self.template['selected_vars']:
                    tmp['selected_vars'][f]=self.template['selected_vars'][f]
                else:
                    tmp['selected_vars'][f]={'aliases':(f,),'read_sel':self.template['config']['read_mode']}
            self.template=tmp
        self.foundvars=None
        self.lastfoundvars=None
        self.globalattributes=None
        self.lastglobalattributes=None
        self.foundvars_uri=None
        self.maxbin=maxbin
        self.preprocess=preprocess
        self.decompFileSource=None
        self.decompFileToDelete=None

    def __del__(self):
        self.cleanDecomp()

    def cleanDecomp(self):
        if hasattr(self,'decompFileToDelete') and self.decompFileToDelete!=None:
            os.unlink(self.decompFileToDelete)
            self.decompFileToDelete=None
            self.decompFileSource=None

    def doDecompress(self,fn):
        if fn==self.decompFileSource:
            return self.decompFileToDelete
        if fn.endswith('.bz2'):
            self.cleanDecomp()
            randomformat='>Q'
            rnd=open('/dev/urandom','r')
            rndn=rnd.read(struct.calcsize(randomformat))
            rnd.close()
            randomtag='%08x' % struct.unpack(randomformat,rndn)
            print 'opening ',fn
            bzf=bz2.BZ2File(fn,'r')
            if os.access(os.path.join('/dev','shm'),os.W_OK):
                filename=os.path.join('/dev','shm','readraw_' + randomtag + '.nc')
            else:
                filename='tmp_readraw_' + randomtag + '.nc'
            fileid=open(filename,'w')
            while True:
                tmp=bzf.read(1024*1024*64)
                if len(tmp)==0:
                    break
                fileid.write(tmp)
            fileid.close()
            bzf.close()
            self.decompFileToDelete=filename
            self.decompFileSource=fn
            return filename
        return fn

    def obtain(self,uri,*args,**kwargs):
        return uri['path']

    def doReadraw(self,uri,*args,**kwargs):
        realuri=uri
        real_kwargs=self.read_raw_parms.copy()
        real_kwargs.update(kwargs)
        try:
            return hru.read_raw(self.doDecompress(realuri),self.template_prefix,self.template,self.maxbin,self.preprocess,*args,**real_kwargs)
        except IOError:#bad file
            return None,None,None

    def open(self,uri,*args,**kwargs):
        """ open the named file (decompressing if needed), and return content

            :param uri: path
        """
        #return netcdf structure with metaframe
        a,b,c=self.doReadraw(uri,*args,**kwargs)
        self.foundvars=b
        self.globalattributes=c
        if b!=None:
            self.lastfoundvars=b
        if c!=None:
            self.lastglobalattributes=c
        self.foundvars_uri=uri
        #print 'Zookeeper read from',realuri,':', a
        return a

    def getFoundVars(self,uri):#FIXME this should go into being provides...
        if self.foundvars_uri==uri and self.foundvars!=None:
            return self.foundvars
        a,b,c=self.doReadraw(uri,doread=False)
        if b!=None:
            self.lastfoundvars=b
        if c!=None:
            self.lastglobalattributes=c
        return b

    def getAttributes(self,uri):#FIXME this should go into being provides...
        if self.foundvars_uri==uri and self.globalattributes!=None:
            return self.globalattributes
        a,b,c=self.doReadraw(uri,doread=False)
        if b!=None:
            self.lastfoundvars=b
        if c!=None:
            self.lastglobalattributes=c
        return c
