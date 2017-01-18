
import dplkit.role.narrator
import numpy
import copy
import lg_base.core.array_utils as hau

@dplkit.role.decorator.exposes_attrs_in_chain(['requested_altitudes'])
@dplkit.role.decorator.exposes_attrs_of_field('framestream')
class ResampleXd(dplkit.role.narrator.aNarrator):
    """ interpolate the altitude axis of a framestream filter

    :param framestream: input framestream
    :param sourceaxname: source altitude axis field name in the frame of the stream
    :param destax: array object of the intended result altitude axis
    :param left: parameter to numpy.interp
    :param right: parameter to numpy.interp
    """
    def __init__(self,framestream,sourceaxname,destax,order=1,left=None,right=None):
        self.framestream=framestream
        self.sourceaxname=sourceaxname
        self.destax=destax
        self.order=order
        self.interpparams=dict(left=left,right=right)

    def requested_altitudes(self):
        return self.destax

    def read(self):
        daxi=-1
        destlen=0
        for shpi in range(len(self.destax.shape)):
            if self.destax.shape[shpi]>1:
                daxi=shpi
                destlen=self.destax.shape[shpi]
        for f in self.framestream:
            if hasattr(f,self.sourceaxname):
                f=copy.copy(f)
                sourceax=getattr(f,self.sourceaxname).copy()
                saxi=-1
                sourcelen=0
                for shpi in range(len(sourceax.shape)):
                    if sourceax.shape[shpi]>1:
                        saxi=shpi
                        sourcelen=sourceax.shape[shpi]
                for k,v in vars(f).items():
                    if not hasattr(v,'shape'):
                        continue
                    shp=list(v.shape)
                    shpi=-1
                    if isinstance(f,hau.Time_Z_Group):
                        if isinstance(v,hau.TZ_Array):
                            shpi=1
                        elif isinstance(v,hau.T_Array):
                            continue
                        elif isinstance(v,hau.Z_Array):
                            shpi=0
                    else:
                        for i in range(len(shp)):
                            if shp[i]==sourcelen:
                                shpi=i
                    if shpi<0 or (len(shp)==1 and k==self.sourceaxname):
                        continue
                    shp[shpi]=destlen
                    print 'resizing',k,'to',shp,'type',type(v)
                    #time.sleep(5)
                    #print hasattr(sourceax.data,'ravel'),hasattr(v.data,'ravel'),hasattr(self.destax.data,'ravel')
                    tmp=numpy.zeros(shp)
                    if type(tmp)!=type(v):
                        tmp=type(v)(tmp)
                    setattr(f,k,tmp)
                    if len(shp)==1:
                        tmp[:]=numpy.interp(self.destax.ravel(),sourceax.ravel(),v.ravel(),**self.interpparams)
                    elif len(shp)==2:
                        tmp[:,:]=numpy.NaN
                        if shpi==0:
                            for x in range(shp[1]):
                                tmp[:,x]=numpy.interp(self.destax.ravel(),sourceax.ravel(),v[:,x].ravel(),**self.interpparams)
                        elif shpi==1:
                            for x in range(shp[0]):
                                tmp[x,:]=numpy.interp(self.destax.ravel(),sourceax.ravel(),v[x,:].ravel(),**self.interpparams)
                setattr(f,self.sourceaxname,copy.deepcopy(self.destax))
            else:
                print "frame doesn't have",self.sourceaxname
                #time.sleep(5)
            yield f

if __name__ == '__main__':
    main()
