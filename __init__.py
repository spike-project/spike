#!/usr/bin/env python 
 # encoding: utf-8
"""
The Spike Package

"""
# from __future__ import absolute_import not needed anymore in 2.7
import os, sys

# Import exceptions
from NPKError import NPKError

# version.py defines static version names
from version import version as __version__
from version import VersionName as __version_info__
from version import ProgramName as __program_name__
#from version import revision as __revision__
from version import rev_date as __date__

__author__ = "Marc A. Delsuc <delsuc@igbmc.fr>, Marie-Aude Coutouly, Lionel Chiron"
SPIKE_version = __version__


#### Header to set-up the whole SPIKE environment
import NPKData

### plugins to the spike.NPKData class.
# simply put a xxx.py in the plugins folder - and define the interface as described in the doc
from plugins import load

load()

#####################################################
# This code allows to pickle methods, and thus to use multiprocessing on methods
# This very nice trick comes from :
# Steven Bethard
# http://bytes.com/topic/python/answers/552476-why-cant-you-pickle-instancemethods
#
# check also
# http://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-pythons-multiprocessing-pool-ma
#
def _pickle_method(method):
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    return _unpickle_method, (func_name, obj, cls)

def _unpickle_method(func_name, obj, cls):
    for cls in cls.mro():
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
    return func.__get__(obj, cls)
import copy_reg
import types
copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)
# end of the trick
#####################################################

