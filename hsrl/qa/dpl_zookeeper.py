import dplkit.role.zookeeper
from parsing import fileParser

class QualityAssuranceZookeeper(dplkit.role.zookeeper.aZookeeper):
    def __init__(self):
        super(QualityAssuranceZookeeper,self).__init__()
        self.parser=fileParser()

    def obtain(self,uri, *args, **kwargs):
        if isinstance(uri,dict):
            return self.parser.parseFile(uri['path'],uri['start'],*args,**kwargs)
        try:
            ret=None
            for f in iter(uri):
                ret=self.obtain(f,destination=ret)
            return ret
        except:
            raise#FIXME
        raise RuntimeError('Unknown URI type',type(uri),':',uri)
