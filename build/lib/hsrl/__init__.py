# hsrl python modules
import os,os.path
import sys

def local_path():
    return os.path.dirname(__file__)

# to allow a smooth transition, add the subdirectories to sys.path

#for d in __all__:
#    sys.path.append(os.path.join(local_path(), d))

