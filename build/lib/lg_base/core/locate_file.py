import os, os.path, sys, types
import logging
LOG = logging.getLogger(__name__)
import copy
DBG = logging.info
from time import sleep

def codebase_dir():
    ret=os.path.dirname(os.path.dirname(os.path.dirname(sys._getframe(0).f_code.co_filename)))
    ret=os.path.realpath(ret)
    return ret

def execution_dir(level=0):
    """ return the directory this code is executing from """
    ret=os.path.dirname(sys._getframe(1+level).f_code.co_filename)
    ret=os.path.realpath(ret)
    return ret
    
def getModuleNamed(name):
    if name not in sys.modules:
        raise ImportError('Module '+name+" hasn't been loaded. Will not hot-load for a locate call")
        #import importlib
        #return importlib.import_module(modul)
    return sys.modules[name]

def module_path(modul):
    assert(hasattr(modul,'__file__'))
    return os.path.realpath(os.path.dirname(modul.__file__))

def pathAscend(depth,topdir):
    adir=copy.copy(depth)
    if not adir.startswith(topdir):
        yield adir
        return
    while adir.startswith(topdir):
        yield adir
        adir=os.path.dirname(adir)

def moduledirs(forModuleParm,default_val,top_dir):
    if forModuleParm==None:
        yield default_val
    elif isinstance(forModuleParm,types.ModuleType):
        yield module_path(forModuleParm)
    elif len(forModuleParm)==0:
        yield top_dir
    else:
        for m in forModuleParm:
            yield module_path(m)

def locate_file(filename,forModule=None,level=0,systemOnly=False):
    """ locate a specified file in various standard locations
    
    if filename exists relative to the current directory, use it
    if the file exists relative to $HSRL_CONFIG use it
    if the file exists relative to the hsrl_config directory - use it
    
    relative locations are default of the calling module. this can be overridden by passing a module to forModule
    if forModule is iterable, it will instead iterate thru each module in the iterable
    if the iterable is of zero length, it will only to the top of this module

    """    

    pathschecked=[]
    if os.environ.has_key('OVERRIDE_SYSTEM_DEFAULTS'):
        systemOnly=False
    if not systemOnly:
        pathschecked.append('./')
        DBG('checking for %s' % filename)
    if os.path.isfile(filename):
        if systemOnly:
            LOG.warning('Found '+filename+' in current directory. Ignoring. Define environment OVERRIDE_SYSTEM_DEFAULTS to override')
        else:
            DBG('returning for '+filename)
            return filename
    top_dir = codebase_dir()
    user_tops = []
    if os.environ.has_key('HSRL_CONFIG'):
        _user_tops = os.environ['HSRL_CONFIG'].split(':')
        for user_top in _user_tops:
            err=None
            if len(user_top)==0:
                continue
            if user_top.startswith(top_dir):
                newerr = 'Potential Inconsistency problem: HSRL_CONFIG user path "%s" is within the system path "%s". Misconfigured. Should be moved to another private and separate location' % (user_top,top_dir)
                LOG.error(newerr)
                err = ((err + '\n') if err else '') + newerr
            if not os.path.isdir(user_top):
                newerr = 'Potential Inconsistency problem: HSRL_CONFIG user path "%s" isn\'t a directory. Misconfigured. Shouldn\'t be set if it\'s not used' % (user_top)
                LOG.error(newerr)
                err = ((err + '\n') if err else '') + newerr
            if err:
                raise RuntimeError(err)
            user_tops.append(user_top)
    #caller_dir = execution_dir(level=1+level) if forModule==None else module_path(forModule)
    for caller_dir in moduledirs(forModule,execution_dir(level=1+level),top_dir):
        
        DBG('SEARCH path for modules: '+caller_dir + ' '+top_dir)
        for user_top in user_tops:
            for check_dir in pathAscend(caller_dir,top_dir):
            #while check_dir.startswith(top_dir):
                config_file=os.path.join(check_dir.replace(top_dir,user_top),filename)
                if not systemOnly:
                    DBG('checking for HSRL_CONFIG '+config_file)
                    pathschecked.append(config_file)
                if os.path.isfile(config_file):
                    if systemOnly:
                        LOG.warning('Found '+filename+' at '+config_file+'. Ignoring. Define environment OVERRIDE_SYSTEM_DEFAULTS to override')
                    else:
                        DBG('found '+config_file)
                        return config_file

        for check_dir in pathAscend(caller_dir,top_dir):
            config_dir=os.path.join(check_dir,'config')
            config_file=os.path.join(config_dir,filename)
            DBG('checking for system '+config_file)
            pathschecked.append(config_file)
            if os.path.isfile(config_file):
                DBG('found '+config_file)
                return config_file

    config_dir = os.path.join(  top_dir ,"hsrl_config" ) #this path is not expected to be 'deployed' with the code, so it shouldn't be used. Recommended explicitly pathed with HSRL_CONFIG if is a do-not-care out of your own home directory
    config_file = os.path.join( config_dir, filename )
    DBG('checking for %s' % config_file)
    pathschecked.append(config_file)
    if os.path.isfile(config_file):
        LOG.error("FALLING BACK TO USING DEPLOYED hsrl_config FOLDER. This is not supposed to be needed.\n" \
            "Any files that exist should be in your HSRL_CONFIG environment-configured folder, or\n" \
            "deployed with the code in the appropriate modules 'config' directory.\n" \
            "Returning value "+config_file)
        sleep(10)
        return config_file
    
    LOG.error("preference file, '"+filename+"' not found")

    raise IOError('could not locate %s in %s' % (filename, ', '.join(pathschecked))) 
    
