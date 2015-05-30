"""
Plug-ins for the Spike package

Each plugin should define the needed functions :

def myfunc(npkdata, args):
    "myfunc doc"
    ...do whatever, assuming npkdata is a NPKData
    return npkdata     # THIS is important, that is the standard NPKData mechanism

and register them into NPKData as follows :
NPKData_plugin("myname", myfunc)

then, any NPKData will inherit the myname() method
"""
import imp
import traceback
import unittest

def load():
    "the load() function is called at initialization, and loads all files found in the plugins folder"
    import glob, os
    from spike.NPKData import NPKData_plugin

# import all python code found here ( __path__[0] ) except me !
    for pgfile in glob.glob( os.path.join(__path__[0],"*.py") ):
        b = os.path.basename(pgfile)    # code.py
        pgmod = os.path.splitext(b)[0]  # code
        if not pgfile.endswith('__init__.py'):
            print "Importing plugin << %s >>"%pgmod
            try:
                m = imp.load_source(pgmod, pgfile)
            except:
                print "*** Failed ***"
                print traceback.print_exc()
                print "*** Continuing ***"

class PluginTests(unittest.TestCase):
    def test_plugin(self):
        '''Test of plugin mechanism'''
        from spike.NPKData import NPKData
        d1 = NPKData()
        d1.fake("this is a fake title").rfft()
        self.assertTrue( d1.test_title == "this is a fake title")

