#!/usr/bin/env python
# encoding: utf-8
"""
SetConfig.py

Created by mac on 2012-02-27.
Copyright (c) 2012 __NMRTEC__. All rights reserved.
"""

from __future__ import print_function
import os
from os import path
import sys
#from tempfile import mkstemp
from configobj import ConfigObj     # standard dans Enthought !!!

def set_config(paramdic, configfilename="process.mscf"):
    """
    update the configfile for the parameters given in argument
    parameters are given in a dictionary
    """
    config = ConfigObj(configfilename)  # load it as a dictionnary
    params = set(paramdic.keys())            # get all parameters to set
    done = set()
    for section in config.keys():
        for item in config[section].keys():
            if item in params:
                done.add(item)        # keep a copy for final check
                config[section][item] = paramdic[item]
    config.write()  # c'est tout !
    if params != done:      # je laisse quand meme le test
        print("WARNING, problem in %s"%configfilename)
        print("the following parameter(s) were not found, and thus not updated:")
        for i in sparams.difference(done):      # contains element in sparams but not in done
            print(" - "+i)



if __name__ == '__main__':
    set_config({"GLOBAL_PATH":"/home/mac"})
        
   