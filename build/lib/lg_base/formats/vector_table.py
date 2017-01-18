import json
import lg_base.core.read_utilities as ru
import numpy
from collections import namedtuple

ColumnHeader=namedtuple("ColumnHeader","name unit")

def breakUpLine(l,colcount):
    ret=[]
    for c in l.strip().split():
        unit=None
        if '(' in c and c[-1]==')':
            unit=c[(c.find('(')+1):-1]
            c=c[:c.find('(')]
            if len(c)==0:
                if len(ret)>0 and ret[-1].unit is None:
                    ret[-1]=ColumnHeader(ret[-1].name,unit)
                else:
                    ret.append(ColumnHeader(None,unit))
                continue
        c=c.replace('/','_divby_').replace('-','_')
        ret.append(ColumnHeader(c,unit))
    if len(ret)!=colcount:
        return None
    return tuple(ret)

def mergeheaders(first,second):
    if first is None:
        return second
    if second is None:
        return first
    ret=[]
    for i in range(len(first)):
        assert(first[i].name is not None or second[i].name is not None)
        if first[i].name is None or second[i].name is None:
            n=second[i].name or first[i].name
        else:
            n=first[i].name+'_'+second[i].name
        assert(first[i].unit is None or second[i].unit is None)
        u=first[i].unit or second[i].unit
        ret.append(ColumnHeader(n,u))
    return tuple(ret)

def printColumnNames(cols):
    if cols is None:
        return
    for i,x in enumerate(cols):
        print i,':',x.name+("" if x.unit is None else ('('+x.unit+')'))

def findColumnNames(filename,header,colcount):#THIS IS WAY TO COMPLICATED
    lines=header.strip().split('\n')
    colrow=len(lines)-1
    colheaders=breakUpLine(lines[colrow][1:],colcount)
    while colrow>0:
        colrow-=1
        for x in ('=','-----',':','>','<'):
            if x in lines[colrow]:
                lines[colrow]=None
                break
        if lines[colrow] is None:
            continue
        ncolheaders=breakUpLine(lines[colrow][1:],colcount)
        if ncolheaders is not None:
            colheaders=mergeheaders(ncolheaders,colheaders)
            print "WARNING: LAST LINE IN HEADER DOESN'T HAVE COLUMN NAMES IN "+filename
        elif colheaders is not None:
            break
    if colheaders is not None:
        chk={}
        for x in colheaders:
            if x.name in chk:
                print 'Duplicate column name in column headers in '+filename
                colheaders=None
                break
            if len(x.name)<3:
                print 'WARNING insist column headers are longer than 2 characters in '+filename
                colheaders=None
                break
            chk[x.name]=x
    if colheaders is not None and len(colheaders)!=colcount:
        print 'Header names mismatched in %s. ignoring' % (filename)
    #else:
    #      printColumnNames(colheaders)
    return colheaders


class calibration_vector_table(object):

    def __init__(self,filename,expire_time,isJson=None,max_range_bin=None,requireColumns=False,noMatrix=False,fallbackColumns=[]):
        self._filename = filename
        self._expire_time = expire_time
        if isJson is None and filename is not None:
            isJson=filename.endswith('.json')
        if filename is None:
            self.header = None
            self.data = None
        elif not isJson:
            self.header, self.data = ru.readascii(filename)
        else:
            f = open(filename, 'r')
            obj = json.load(f)
            f.close()
            self.header = obj['header']
            self.data = obj['calibration']
        if max_range_bin and self.data is not None:
            self.data = self.data[:max_range_bin, :]
        setcolumns=False
        if self.header is not None:
            dat=self.data
            if not hasattr(dat,'shape'):
                dat=numpy.array(dat)
            colheaders=None
            if len(dat.shape)>=2:
                colheaders=findColumnNames(filename,self.header,dat.shape[1])
                if colheaders is None and fallbackColumns is not None and len(fallbackColumns)>=dat.shape[1]:
                    print "WARNING: Using fallback column names for file "+filename
                    colheaders=[ColumnHeader(c,None) for c in fallbackColumns[:dat.shape[1]]]
            if colheaders is not None:
                setcolumns=True
                fields=[]
                for i,c in enumerate(colheaders):
                    if i==0:
                        setattr(self,'axis',c.name)
                    else:
                        fields.append(c.name)
                    setattr(self,c.name,dat[:,i])
                    if c.unit is not None:
                        setattr(self,c.name+'_unit',c.unit)
                setattr(self,'fields',tuple(fields))
        if (noMatrix or requireColumns) and not setcolumns:
            raise RuntimeError(filename+" is required to have column names in the header.")
        if noMatrix:
            delattr(self,'data')

    @property
    def filename(self):
        return self._filename

    @property
    def expire_time(self):
        return self._expire_time


def main():
    import sys,os
    for x in sys.argv[1:]:
        v=vars(calibration_vector_table(x,None))
        #print v.keys()
        for k,vec in v.items():
            if k in ('header','data') or k.startswith('_'):
                continue
            print k,':',vec

if __name__ == '__main__':
    main()