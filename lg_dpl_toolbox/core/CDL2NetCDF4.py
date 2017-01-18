import sys
import os

def CDL2NetCDF4(cdl,nc):
    """ Creates a mostly empty NetCDF4 file from a CDL file. Contents are uninitialized

    :param cdl: CDL template file name
    :param nc: NetCDF4 template file name for Output
    
    Will raise a RuntimeError if the command fails
    """
    (readp,writep)=os.pipe()
    pid=os.fork()
    if pid==0:
        os.close(readp)
        os.dup2(writep,sys.stdout.fileno())
        os.dup2(writep,sys.stderr.fileno())
        os.execvp('ncgen',('ncgen','-b','-k','hdf5','-o',nc,'-x',cdl))
        print 'failed to exec. ncgen executable in path?'
        exit(-1)
    os.close(writep)
    if pid<0:
        raise RuntimeError, 'Failed to fork'
    (p,res)=os.waitpid(pid,0)
    if res!=0:
        raise RuntimeError, 'Error %i.  Output:%s' % (res,os.read(readp,4096))
    os.close(readp)
