#!/usr/bin/env python 
# encoding: utf-8
"""
Tests.py

Created by Marc-André  on 2010-07-20.

Runs tests on selected modules using the integrated unittests in the different SPIKE modules.

most default values can be overloaded with run time arguments

Example on a module : python -m spike.Tests -D DATA_test -t File.Apex

"""

#added arguments on 2015-09-03

from __future__ import  print_function, division
import unittest
import os
import os.path as op
import sys
from .Display import testplot
# Beware  !! Awful Hack !!
# testplot is used to switch display off 
# BUT ALSO to monkey patch at Tests run time to store in testplot.config values used by the tests (DATA_dir)


__version__ = "2.0"

#############################################################################
# configuration - overwriten by arguments at run time

# when CLEAN is set to true, unwanted files in test directory will be removed.
CLEAN = True

# when MAIL is set to true, a mail with results is sent to e-mail add defined in list_of_mails
# defined with
MAIL = False
list_of_mails = []

#DATA_dir defines where the DATA for tests are located
DATA_dir = "/Volume/DATA_test"

RUN = True

# Add your module here
mod_util = ("plugins", "plugins.Peaks", 'util.dynsubplot', 'util.debug_tools') #'util.read_msh5', 
mod_algo = ('Algo.Cadzow', 'Algo.Linpredic', 'Algo.urQRd', 'Algo.SL0', 'Algo.maxent', 'Algo.BC') 
 
mod_file = ("File.BrukerNMR", "File.GifaFile", 'File.HDF5File', 'File.Apex', 'File.csv', 'File.Solarix')
mod_basicproc = ("NPKData", "FTICR", "Orbitrap")
mod_user = ('processing', )

list_of_modules = mod_basicproc + mod_file  + mod_util + mod_algo # + mod_user

# end of configuration
#############################################################################


# utilities to be called by tests using files in DATA_dir
def directory():
    "returns the location of the directory containing dataset for tests"
    return testplot.config["DATA_dir"]
def filename(name):
    "returns the full name of a test dataset located in the test directory"
    return op.join(testplot.config["DATA_dir"], name)

def msg(st, sep = '='):
    '''
    Message in Tests.py
    '''
    s = sep*(len(st) + 4)+"\n"
    s = s+ '| '+ st+ ' |'+"\n"
    s = s + sep*(len(st) + 4)+"\n"
    print(s)
    return s

def cleanspike():
    '''
    Removes the .pyc in spike
    '''
    for root, dirs, files in os.walk('.'):
        #print root, dirs, files
        for f in files:
            r,ext = os.path.splitext(f)
            if ext == '.pyc':
                addr = os.path.join(root,f)
                print(addr)
                os.remove(addr)

def cleandir():
    "checking files in DATA_dir directory and removes files created by previous tests"
    import glob
    files_to_keep = ('ubiquitin_5_scan_res_30000_1.dat','cytoC_ms_1scan_000001.d', 'cytoC_2D_000001.d',
                'dosy-cluster2-corr.gs2', 'dosy-cluster2.gs2',
                'proj.gs1', 'ubiquitine_2D_000002.d','file_fticr',
                'file_npar','npar_fticr','Lasalocid-Tocsy', 'Sampling_file.list', 
                'ubiquitine_2D_000002_Sampling_2k.list',
                'Sampling_file_aposteriori_cytoCpnas.list','angio_ms_000005.d')
    for i in glob.glob(filename("*")):
        print(i, end=' ')
        if os.path.basename(i) in files_to_keep:
            print(" Ok")
        else:
            if CLEAN:
                try:
                    os.remove(i)
                    print(" removed")
                except OSError:
                    print(" **** could not be removed ****")
            else:
                print(" should be removed")

class NPKTest(unittest.TestCase):
    """overload unittest.TestCase for default verbosity - Not Used - """
    def setUp(self):
            self.verbose = 1    # verbose > 0 switches messages on
    def announce(self):
        if self.verbose >0:
            print("\n========",self.shortDescription(),'===============')

def do_Test():
    '''
    Performs all tests then indicates if successfull.
    Gives total time elapsed.
    '''
    import time
    global list_of_modules
    subject = "SPIKE tests perfomed on {0} {2}  running on host {1}".format(*os.uname())
    to_mail = [msg(subject), "Test program version %s"%__version__]

    # add spike prefix
    list_of_modules = ["spike."+mod for mod in list_of_modules]
    msg("modules to be tested are:")
    for mod in list_of_modules:
        print(mod)

    stdata = "\nDatasets for tests are located in : %s"%DATA_dir
    to_mail.append(stdata)
    print(stdata)

    if not RUN:
        print ("\nDry run - stopping now, without executing the tests")
        sys.exit(0)

    msg("First removing leftover files")
    cleandir()
    msg("removing .pyc in spike")
    cleanspike()
    msg("Running automatic Tests")
    t0 = time.time()
    suite = unittest.defaultTestLoader.loadTestsFromNames( list_of_modules )
    nbtests = suite.countTestCases()
    results = unittest.TextTestRunner(verbosity = 2).run(suite)
    elaps = time.time()-t0
    to_mail.append("Ran %d tests in %.3fs"%(nbtests, elaps))
    if results.wasSuccessful():
        to_mail.append( msg("CONGRATULATIONS - all the {} SPIKE tests performed succesfully ".format(len(list_of_modules))) )
        subject = "SUCCESS :)  " + subject
        to_mail.append( msg("modules tested were:"))
        for mod in list_of_modules:
            print(mod)
            to_mail.append(mod)
        to_mail.append( msg("test performed in %.2f sec"%elaps) )
    else:
        subject = "Failed :(  " + subject
        to_mail.append( msg("Tests Failed, Please revise error codes", sep = '!') )
        to_mail.append( "%d test(s) failed :"%len(results.errors) )
        for err in results.errors:
            to_mail.append( msg("%s"%(err[0]) ) )
            to_mail.append( err[1] )
    print("\n".join(to_mail))
    if MAIL:
        from util.sendgmail import mail
        for address in list_of_mails:
            mail(address, subject, "\n".join(to_mail) )
    # finally clean dir
    cleandir()

def main():
#    import Display.testplot as testplot
    import argparse
    global MAIL, list_of_mails, RUN, list_of_modules, DATA_dir, CLEAN
    global config
    
    # Parse and interpret options.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-m', '--mail', action='append', help="a mail address to which the result is sent - multiple entries possible -")
    parser.add_argument('-D', '--Data', help="the location of the data_test folder")
    parser.add_argument('-t', '--test_modules', action='append', help="overwrite the list of modules to test - multiple entries possible -")
    parser.add_argument('-n', '--dry',  action='store_true', help="list parameters and do not run the tests")
    parser.add_argument('-d', '--dirty', action='store_true', help="do not remove temporary files")
    
    args = parser.parse_args()

    print( "mail", args.mail)
    print( "data", args.Data)
    print( "module", args.test_modules)
    print("dry", args.dry)
    print('dirty', args.dirty)
    if args.dry:
        RUN = False
    if args.mail is not None:
        MAIL = True
        list_of_mails =  args.mail  # is a list because of action='append'
    if args.test_modules is not None:
        list_of_modules = args.test_modules
    if args.Data is not None:
        DATA_dir = args.Data
    if args.dirty:
        CLEAN = False
    testplot.PLOT = False   # switches off the display for automatic tests
    testplot.config = {}
    testplot.config["DATA_dir"] = DATA_dir
#    print(sys.modules)
    do_Test()
    
if __name__ == '__main__':
    main()