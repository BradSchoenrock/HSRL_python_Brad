#!/usr/bin/python
# -*- coding: utf-8 -*-
from locate_file import locate_file


def open_config(filename, verbose=True,level=0,*args,**kwargs):
    """ open 'filename' from specified path, DEFAULT_DIR, or
    os.environ['HSRL_CONFIG'], which allows the user to locate configuration
    files anywhere
    """
    return open(locate_file(filename,level=1+level,*args,**kwargs)) #set searchpath to path of caller
