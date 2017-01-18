import dplkit.role.librarian
from datetime import datetime,timedelta
import calendar
import os

#this will contain the librarian, which locates files for a given window
class QualityAssuranceLibrarian(dplkit.role.librarian.aLibrarian):
    def __init__(self,instrument,siteid=None,basedirectory=None,fileprefix=None):
        super(QualityAssuranceLibrarian,self).__init__()
        from hsrl.dpl import HSRLLibrarian #FIXME
        if basedirectory!=None:
            dirname=basedirectory
        elif siteid!=None:
            dirname=HSRLLibrarian.HSRLLibrarian(site=siteid).basedir
        else:
            dirname=HSRLLibrarian.HSRLLibrarian(instrument=instrument).basedir
        self.basedirectory=dirname
        self.instrumentname=instrument
        self.fileprefix=fileprefix or (instrument+'_qa')
        #self.fileprefix='testing_qa'

    def search(self,start_time=None,end_time=None):
        atime=start_time.replace(day=1,hour=0,minute=0,second=0,microsecond=0)
        while atime<(end_time or datetime.utcnow()):
            basedir=os.path.join(self.basedirectory,atime.strftime('%Y'),atime.strftime('%m'))
            pref=atime.strftime(self.fileprefix+'_%Y%m')
            suff='.txt'
            lastfile=os.path.join(basedir,pref+suff)
            try:
                files=os.listdir(basedir)
            except OSError:
                files=[]
            files.sort()
            for f in files:
                if f.startswith(pref) and f.endswith(suff):
                    lastfile=os.path.join(basedir,f)#gets the last
            dur=timedelta(days=calendar.monthrange(atime.year,atime.month)[1])
            #if os.access(lastfile,os.R_OK):
            yield dict(start=atime,width=dur,path=lastfile,filename=os.path.basename(lastfile))
            atime+=dur

def main():
    import sys
    l=QualityAssuranceLibrarian(sys.argv[1])
    s=datetime.strptime(sys.argv[2],"%Y%m%d")
    e=datetime.strptime(sys.argv[3],"%Y%m%d")
    for x in l(s,e):
        print x

if __name__ == '__main__':
    main()