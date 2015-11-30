#!/usr/bin/env python 
 # encoding: utf-8
"""
Plug-ins for the Spike package

All the plugin files located in spike/plugins folder the loaded automatically 
when import importing spike the first time.

the variable spike.plugins.plugins contains the list of the loaded plugin modules.

It is allways possible to load a plugin afterward by import a plugin definition at a later time during run-time.

Each plugin file should define the needed functions :

def myfunc(npkdata, args):
    "myfunc doc"
    ...do whatever, assuming npkdata is a NPKData
    return npkdata     # THIS is important, that is the standard NPKData mechanism

and register them into NPKData as follows :
NPKData_plugin("myname", myfunc)

then, any NPKData will inherit the myname() method

For the moment, only NPKData plugins are handled.

"""
from __future__ import print_function
import imp
import traceback
import unittest

plugins = []  # contains plugin modules loaded so far
import glob, os
from ..NPKData import NPKData_plugin

def load():
    "the load() function is called at initialization, and loads all files found in the plugins folder"
# import all python code found here ( __path__[0] ) except me !
    for pgfile in sorted( glob.glob( os.path.join(__path__[0],"*.py") ) ):
        b = os.path.basename(pgfile)    # code.py
        pgmod = os.path.splitext(b)[0]  # code
        if not pgfile.endswith('__init__.py'):
            loadone(pgmod, pgfile)
                
def loadone(pluginname, pgfile=None):
    """
    loads the plugin "pluginname" found in file pgfile
    if pgfile is None, the file is found in the plugins folder (without the .py extension)
    """
    global plugins
    direc = __path__[0] # plugins folder
    if pgfile is None:
        pgfile = os.path.join(direc, pluginname + ".py")
    print("Importing plugin << %s >>"%pluginname)
    try:
        plugins.remove(pluginname)
        print("WARNING existing plugin %s is overwritten"%pluginname)
    except    ValueError:
        pass
    try:
        m = imp.load_source(pluginname, pgfile)
        plugins.append(pluginname)
    except:
        print("*** Failed ***")
        traceback.print_exc()
        print("*** Continuing ***")

class PluginTests(unittest.TestCase):
    def test_plugin(self):
        '''Test of plugin mechanism'''
        from ..NPKData import NPKData
        d1 = NPKData()
        d1.fake("this is a fake title").rfft()
        self.assertTrue( d1.test_title == "this is a fake title")

