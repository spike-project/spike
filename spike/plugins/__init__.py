
#!/usr/bin/env python 
# encoding: utf-8
"""
Plug-ins for the Spike package

All the plugin files located in the spike/plugins folder are loaded automatically 
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
from __future__ import print_function #
import imp
import traceback
import unittest
from collections import defaultdict
import glob, os
from pathlib import Path

#from ..NPKData import NPKData_plugin

plugins = {}  # contains plugin modules loaded so far
codes = defaultdict(list)  # loaded at plugin injection with all codes

def load(debug=True):
    """
    the load() function is called at initialization, and loads all files found in the plugins folders
    typically
        - plugins folder in distribution
        - $HOME/Spike/plugins
    """
    for folder in [Path(__path__[0]), Path.home()/'Spike'/'plugins']:
        loadfolder(folder, debug=debug)
    print( "\nspike.plugins.report() for a short description of each plugins")
    print( "spike.plugins.report('module_name') for complete documentation on one plugin") 


def loadfolder(folder, debug=True):
    """the loads plugins from a given folder
    import all python code found here  except the ones with names starting with a _
    """
    done = []
    try:
        folder.name
    except AttributeError:
        folder = Path(folder)  # assume it was a string
    for pgfile in sorted( folder.glob("*.py") ):
        b = pgfile.name    # code.py
        pgmod = pgfile.stem  # code
        if not b.startswith('_'):
            loadone(pgmod, str(pgfile), debug=debug)
            done.append(pgmod)
    if done:
        print( "plugins loaded:\n" + 
           " ".join(["{}, ".format(k) for k in done if k in plugins.keys()]) )


def report(module=None, mode=None):
    """ print a nice report
    mode=choose the verbosity
        "short": only one line (default) 
        "detailed": with list of all added methods
        "full": prints the whole plugin doc
    """
    if module is None:
        todo = sorted(plugins.keys())
    else:
        todo = [module]
        if mode is None:
            mode = "full"
    for k in todo:
        doc = plugins[k]
        if mode != 'full':
            print("%s : %s"%(k, doc.split('\n')[0]))  # only first line
        else:
            print("="*60)
            print("%s : %s"%(k, doc))
        if  mode != "short" or mode is None:
            print("    implements: "+" ".join(["{}(), ".format(c) for c in codes[k]]))
    if module is None and mode is None:    # add doc if default values
        print("spike.plugins.report('module_name') for complete documentation on one plugin") 

def loadone(pluginname, pgfile=None, debug=True):
    """
    loads the plugin "pluginname" found in file pgfile
    if pgfile is None, the file is found in the plugins folder (adding the .py extension)
    """
    global plugins
    direc = __path__[0] # plugins folder
    if pgfile is None:
        pgfile = os.path.join(direc, pluginname + ".py")
    if debug: print("---- Importing  << %s >>"%pluginname)
    try:
        plugins[pluginname]
        if debug: print("WARNING existing plugin %s is overwritten"%pluginname)
    except    KeyError:
        pass
    try:
        m = imp.load_source(pluginname, pgfile)
        doc = m.__doc__.split('\n')[0]    # print first line of documentation
        if debug: print("     "+doc)
        plugins[pluginname] = m.__doc__
    except:
        print("*** %s not loaded ***"%pluginname)
        if debug:
            traceback.print_exc()
            print("*** Continuing ***")

class PluginTests(unittest.TestCase):
    def test_plugin(self):
        '''Test of plugin mechanism'''
        from ..NPKData import _NPKData
        d1 = _NPKData()
        d1.fake("this is a fake title").rfft()
        self.assertTrue( d1.test_title == "this is a fake title")

