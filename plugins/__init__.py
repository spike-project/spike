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
import imp
import traceback
import unittest

plugins = []  # contains plugin modules loaded so far

def load():
    "the load() function is called at initialization, and loads all files found in the plugins folder"
    import glob, os
    from ..NPKData import NPKData_plugin
    global plugins
# import all python code found here ( __path__[0] ) except me !
    for pgfile in glob.glob( os.path.join(__path__[0],"*.py") ):
        b = os.path.basename(pgfile)    # code.py
        pgmod = os.path.splitext(b)[0]  # code
        if not pgfile.endswith('__init__.py'):
            print("Importing plugin << %s >>"%pgmod)
            try:
                m = imp.load_source(pgmod, pgfile)
                plugins.append(pgmod)
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

