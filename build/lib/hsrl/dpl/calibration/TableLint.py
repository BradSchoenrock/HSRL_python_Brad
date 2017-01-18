
import sys
from datetime import datetime,timedelta
import os
import calendar
import functools
import dplkit.role.narrator
import dplkit.role.decorator
import traceback

class redirect:
    def doredirect(self,orig,newish):
        assert(isinstance(orig,file))
        tmp=None
        if isinstance(newish,basestring):
           newish=file(newish,"a")
        if isinstance(newish,file):
           tmp=newish
           newish=newish.fileno()
        assert(type(newish)==int)
        os.dup2(newish,orig.fileno())

    def unredirect(self,orig,newish):
        assert(isinstance(orig,file))
        assert(type(newish)==int)
        os.dup2(newish,orig.fileno())
        os.close(newish) 

    def __init__(self,stderr=None,stdout=None):
        self.outfd=None
        self.errfd=None
        if stdout is not None:
            self.outfd=os.dup(sys.stdout.fileno())
            self.doredirect(sys.stdout,stdout)
        if stderr is not None:
            self.errfd=os.dup(sys.stderr.fileno())
            self.doredirect(sys.stderr,stderr)

    def __del__(self):
        self.close()

    def close(self):
        if self.outfd is not None:
            self.unredirect(sys.stdout,self.outfd);
            self.outfd=None
        if self.errfd is not None:
            self.unredirect(sys.stderr,self.errfd);
            self.errfd=None


class TableTestNarrator(dplkit.role.narrator.aNarrator):
    """ Librarian Initialzation for Calibration Tables
  
            :param source: table librarian object
    """
    def __init__(self, source, datatype):
        super(TableTestNarrator,self).__init__()
        self.source=source
        self.datatype=datatype

    def read(self):
        reader=self.datatype
        if isinstance(reader,basestring):
            import TableLibrarian as tl
            reader=tl.calDataInfo(reader)
        if isinstance(reader,dict):
            reader=reader['reader']
        for v in self.source:
            f=v['path']
            e=(v['start']+v['width'])
            r=None
            swapOut=redirect(stderr="/dev/null",stdout="logout")
            try:
                r=reader(None,None,filename=f,expire_time=e)
            except Exception as ex:
                swapOut=None
                print traceback.format_exc()
                yield ex,f,e
            else:
                swapOut=None
                yield r,f,e

def filesmatch(*files):
    oldcontent=None
    oldfile=None
    ret=True
    for f in files:
        content=file(f,"r").read()
        if oldcontent is not None:
            if oldcontent!=content:
                print f,"doesn't match",oldfile
                ret=False
        oldcontent=content
        oldfile=f
    return ret

class CommandLog(object):
    def __init__(self,filename):
        self.script=file(filename,"w")
        os.chmod(filename,0640)

    def __call__(self,*parms):
        for p in parms:
            self.script.write(str(p))
        self.script.write('\n')
        self.script.flush()
   

class CommandScript(object):
    def __init__(self,filename,nogit=False):
        self.script=file(filename,"w")
        self.script.write("#!/bin/bash\n")
        os.chmod(filename,0750)
        self.nogit=nogit

    def __call__(self,*parms):
        if self.nogit and parms[0]=='git':
            if parms[1]=='rm':
                parms=list(parms)
                parms[0]='rm'
                parms[1]='-f'
            else:
                return
        isfirst=True
        for p in parms:
            doquote=False
            for ws in (' ','\t'):
                if ws in p:
                    doquote=True
            if p[0]=='#':
                doquote=False
            if doquote:
                p='"'+p+'"'
            if not isfirst:
                p=' '+p
            self.script.write(p)
            isfirst=False
        self.script.write('\n')
        self.script.flush()

def nextmonth(dt):
    _,md=calendar.monthrange(dt.year,dt.month)
    return dt+timedelta(days=md-dt.day+1)

def monthsteps(start,end):#noninclusive on both ends
    x=start
    while True:
        #_,md=calendar.monthrange(x.year,x.month)
        x=nextmonth(x)#x+timedelta(days=md)
        if x>=end:
            break
        yield x

def monthpath(base,d):
    return os.path.join(base,'%04i' % d.year,'%02i' % d.month)


def startp():
    print 'Tables missing for months:',

def runp(monthpath,date):
    print (date.strftime('%Y-%m')),

def endp():
    print 

defFuncSet=(startp,runp,endp)

def noop():
    pass

def makeCopyMissingMonth(**kwargs):
    return (noop,functools.partial(copyToMissingMonth,**kwargs),noop)

def hasData(path):
    days=os.listdir(path)
    for day in days:
        try:
            d=int(days)
        except:
            continue
        datfolders = os.listdir(os.path.join(path,day))
        for f in datfolders:
            if 'raw' not in datfolders:
                continue
            datafiles=os.listdir(os.path.join(path,day,f))
            for df in datafiles:
                if 'hsrl' in df and '.nc' in df:
                    return True
    return False

def copyMissing(dest,prior,pushcommand,caltype,lint):
    if not os.path.exists(dest) or not hasData(dest):
        return
    lint('should copy ',prior['path'],' forward to ',dest)
    print 'should copy',prior['path'],'to',dest
    pushcommand('cd',dest)#change directories
    pushcommand('cp',prior['path'],prior['filename'])#move the newer one to the other month
    pushcommand('git','add',prior['filename']) #add the new
    pushcommand('git','commit','-m','copying prior for coverage of beginning of month '+caltype)


def copyToMissingMonth(monthpath,date,prior,**kwargs):
    dest=monthpath+os.path.dirname(prior['path'])[len(monthpath):]
    copyMissing(dest,prior,**kwargs)

def checkPaths(base,start,end,funcset=None):
    hadPrint=False
    if funcset is None:
        funcset=defFuncSet
    for m in monthsteps(start,end):
        p=monthpath(base,m)
        if os.path.exists(p):
            if not hadPrint:
                funcset[0]()
                hadPrint=True
            funcset[1](p,m)
    if hadPrint:
        funcset[2]()

def storeError(fi,err):
    fi.write('%s %s\n' % (err[1],err[0]))

def main():
    import TableLibrarian as tl
    parms=dict(completeList=True)
    if "-h"  in sys.argv:
        print 'TableLint - HSRL Calibration Table Lint - verify time order accessibility for consistency, and verify format optionally'
        print '\t-h          this help'
        print "\t-v          read each file, verifying it is properly formatted"
        print "\t -a path    instead of default paths, verify an alternative path"
        return
    verify='-v' in sys.argv
    alternative='-a' in sys.argv
    if verify:
        sys.argv.remove('-v')
    if alternative:
        i=sys.argv.index('-a')
        if len(sys.argv[i])>2:
            alternative=sys.argv[i][2:]
        else:
            alternative=sys.argv[i+1]
            sys.argv.remove(sys.argv[i+1])
        sys.argv.remove('-a')
        parms['alternativePath']=alternative
    inst=[sys.argv[1]] if len(sys.argv)>1 else ('ahsrl','gvhsrl','mf2hsrl','nshsrl','bagohsrl')
    caltypes=[sys.argv[2]] if len(sys.argv)>2 else tl.callist().keys()
    verbose=len(sys.argv)>1
    proximity=0.0 #this has a gotcha of needing to check things across directories. makes audit impossible
    errorlist=file('errfiles.txt','w')
    lint=CommandLog("lint.txt")
    pushcommand=CommandScript('fixit_cals.sh',nogit=alternative)
    for i in inst:
        for c in caltypes:
            #if verbose:
            print i,c
            lib=tl.TableLibrarian(i,c,**parms)
            nar=lib(datetime(2000,1,1,0,0,0),None)
            prior=None
            base=None
            priorgood=None
            priormonthgaptime=None
            monthgaptime=None #used to check there is a copy of the last file from prior month or one from time 0 of the month
            for f in nar:
                if verbose:
                    print f
                if base is None:
                    yr='%04i' % f['pathtime'].year
                    base=f['path']
                    while yr in base:
                        base=os.path.dirname(base)
                    #print 'base is ',base
                if monthgaptime is None or monthgaptime!=f['pathtime']:
                    priormonthgaptime=monthgaptime
                    monthgaptime=f['pathtime']
                    if f['start']>f['pathtime']:
                        if prior is None:
                            #print "its the first, doesnt matter"
                            pass
                        elif prior['start']>monthgaptime:
                            print 'premature order '+c+' at '+prior['path']
                            lint('premature order '+c+' at '+prior['path']+".")
                            lint('should delete (',prior['path'],' has time ',prior['start'],' is after month directory end time ', monthgaptime,')')
                            pushcommand('cd',os.path.dirname(prior['path']))#change directories
                            pushcommand('git','rm',prior['filename'],"# "+str(monthgaptime)+" is before files "+str(prior['start']))
                            if priorgood:
                                pushcommand('cp',priorgood['path'],priorgood['filename'],"# "+str(monthgaptime)+" will include "+str(priorgood['start']))#move the newer one to the other month
                                lint("  and use ",priorgood['path'])
                                pushcommand('git','add',priorgood['filename']) #add the new
                            pushcommand('git','commit','-m','fix premature file '+c+' with immediately prior one')
                        else:
                            print 'first file in folder',f['path'],'starts after the path time',monthgaptime,'(',f['start'],')'
                            lint('first file in folder ',f['path'],' starts after the path time ',monthgaptime,' (',f['start'],')')
                            copyMissing(os.path.dirname(f['path']),prior,pushcommand,c,lint)
                    else:
                        priorgood=f
                    if priormonthgaptime is not None:
                        if (monthgaptime-priormonthgaptime).total_seconds()>(32*24*60*60):#32 day gap
                            #print 'Missing at least one month between ',priormonthgaptime,'and',monthgaptime
                            checkPaths(base,priormonthgaptime,monthgaptime,makeCopyMissingMonth(prior=prior,pushcommand=pushcommand,caltype=c,lint=lint))

                timeeq=False if prior is None else (f['start']==prior['start'])
                difftime = (proximity+1.0) if prior is None else (f['start']-prior['start']).total_seconds()
                if timeeq or difftime>proximity:
                    if timeeq and not filesmatch(f['path'],prior['path']):
                        lint(f['path']," has same name as ",prior['path'],". This is crap.")
                        print 'Files',f['path'],'and',prior['path'],'have same time',prior['start'],'but differ in content! FAIL!'
                        if f['filename']!=prior['filename']:
                            print 'NAMES ALSO MISMATCH'
                        if False:
                            pushcommand("diff","-u",prior['path'],f['path'])
                        else:
                            lint(f['path']," needs to be copied to ",os.path.dirname(prior['path']))
                            pushcommand('cd',os.path.dirname(prior['path']))#change directories
                            pushcommand('cp',f['path'],f['filename'])#move the newer one to the other month
                            if f['filename']!=prior['filename']:
                                lint("  and ",prior['path']," must be deleted")
                                pushcommand('git','rm',prior['filename'])
                            pushcommand('git','add',f['filename']) #add the new
                            pushcommand('git','commit','-m','fix matching file content '+c)
                    prior=f
                elif difftime>0.0 and difftime<=proximity:#too close
                    lint(f['path']," comes too close after ",prior['path'])
                    print 'Close proximity override',f['path'],'comes seconds after',prior['start']
                    #df=os.path.join(os.path.dirname(fnr['path']),pfnr['filename'])
                    pushcommand('cd',os.path.dirname(prior['path']))#change directories
                    pushcommand('cp',f['path'],'.')#copy the newer one over the prior
                    if os.path.dirname(prior['path'])==os.path.dirname(f['path']):
                        pushcommand('git','rm',f['filename'])
                    else:
                        pushcommand('cd',os.path.dirname(f['path']))#change directories
                        pushcommand('git','rm',f['filename'])
                        pushcommand('git','commit','-m','fix proximity '+c)
                        pushcommand('cd',os.path.dirname(prior['path']))#change directories
                    pushcommand('git','add',prior['filename'])
                    pushcommand('git','commit','-m','fix proximity '+c)
                    #print 'Recommended commands:'
                    for x in cmdlist:
                        script.write('"'+('" "'.join(x))+'"\n')                   
                else:
                    lint(f['path']," comes before file ",prior['path'],". copy newer backward")
                    print 'Out of order',f['path'],'comes before',prior['start']
                    #df=os.path.join(os.path.dirname(fnr['path']),pfnr['filename'])
                    pushcommand('cd',os.path.dirname(f['path']))#change directories
                    pushcommand('cp',prior['path'],prior['filename'])#move the newer one to the other month
                    pushcommand('git','rm',f['filename']) #remove the old one that is out of order
                    pushcommand('git','add',prior['filename']) #add the new
                    pushcommand('git','commit','-m','fix order '+c)
            if verify:
                sys.stdout.write('Verify Parsing %s ' % c)
                sys.stdout.flush()
                for f in TableTestNarrator(nar,c):
                    if isinstance(f[0],Exception):
                        storeError(errorlist,f)
                        sys.stdout.write('X')
                    else:
                        sys.stdout.write(".")# f
                    sys.stdout.flush()
                sys.stdout.write('done!\n')
                sys.stdout.flush()
        if monthgaptime is None:
            print 'Instrument has no such tables'
        elif (datetime.utcnow()-monthgaptime).total_seconds()>(32*24*60*60):
            #print 'Missing any new tables. Newest month is',monthgaptime
            checkPaths(base,monthgaptime,datetime.utcnow(),makeCopyMissingMonth(prior=prior,pushcommand=pushcommand,caltype=c,lint=lint))
    del pushcommand

if __name__ == '__main__':
    sys.path.append(os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),'..','..','..'))
    main()
