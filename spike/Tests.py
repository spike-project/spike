#!/usr/bin/env python 
# encoding: utf-8
"""
Tests.py

Created by Marc-André  on 2010-07-20.

Runs tests on selected modules using the integrated unittests in the different SPIKE modules.

most default values can be overloaded with run time arguments

Example on a module :
    python -m spike.Tests -D DATA_test -t File.Apex

"""

#added arguments on 2015-09-03

from __future__ import  print_function, division
import unittest
import os
import os.path as op
import shutil
import sys
import warnings


from .Display import testplot
# Beware  !! Awful Hack !!
# testplot is used to switch display off 
# BUT ALSO to monkey patch at Tests run time to store in testplot.config values used by the tests (DATA_dir)

__version__ = "2.1"

#############################################################################
# configuration - overwriten by arguments at run time

# when CLEAN is set to true, unwanted files in test directory will be removed.
CLEAN = True

# when MAIL is set to true, a mail with results is sent to e-mail add defined in list_of_mails
# defined with
MAIL = False
list_of_mails = []

#DATA_dir defines where the DATA for tests are located
DATA_dir = "/DATA/DATA_test"

RUN = True

# Add your module here
mod_util = ("plugins", 'util.dynsubplot', 'util.debug_tools')  #'util.read_msh5', 
mod_algo = ('Algo.Cadzow', 'Algo.Linpredic', 'Algo.urQRd', 'Algo.sane', 'Algo.SL0', 'Algo.maxent', 'Algo.BC') 
mod_plugins = ("plugins.Peaks", "plugins.Fitter", "plugins.NMR.Bruker_NMR_FT", "plugins.Linear_prediction", "plugins.NMR._PALMA-Test")
mod_file = ("File.BrukerNMR", "File.GifaFile", 'File.HDF5File', 'File.Apex', 'File.csv', 'File.Solarix', 'File.mzXML')
mod_basicproc = ("NPKData", "NMR", "FTICR", "Orbitrap", 'NPKConfigParser')
mod_user = ('processing', )

list_of_modules = mod_basicproc + mod_file  + mod_util + mod_algo + mod_plugins + mod_user

# end of configuration
#############################################################################


# utilities to be called by tests using files in DATA_dir
def directory():
    "returns the location of the directory containing dataset for tests"
    try:
        dd = testplot.config["DATA_dir"]
    except AttributeError:  # may happen 
        dd = DATA_dir
    return dd
def filename(name):
    "returns the full name of a test dataset located in the test directory"
    return op.join(directory(), name)

def msg(st, sep = '='):
    '''
    Message in Tests.py
    '''
    s = sep*(len(st) + 4)+"\n"
    s = s+ '| '+ st+ ' |'+"\n"
    s = s + sep*(len(st) + 4)+"\n"
#    print(s)
    return s

def cleanspike():
    '''
    Removes the .pyc in spike
    '''
    for root, dirs, files in os.walk('spike'):
        #print root, dirs, files
        for f in files:
            r,ext = os.path.splitext(f)
            if ext in ('pyc','pyo'):
                addr = os.path.join(root,f)
                print(addr)
                os.remove(addr)
        for d in dirs:
            if d == '__pycache__':
                shutil.rmtree(os.path.join(root,d))

def cleandir(verbose=True):
    "checking files in DATA_dir directory and removes files created by previous tests"
    import glob
    files_to_keep = ('ubiquitin_5_scan_res_30000_1.dat','cytoC_ms_1scan_000001.d', 'cytoC_2D_000001.d',
                'dosy-cluster2-corr.gs2', 'dosy-cluster2.gs2',
                'proj.gs1', 'ubiquitine_2D_000002.d','Lasalocid-Tocsy', 'Sampling_file.list', 
                'ubiquitine_2D_000002_Sampling_2k.list','test.mscf','test_.mscf',
                'ubiquitine_2D_000002.msh5','ubiquitine_2D_000002_mr.msh5',   # these two for testing 2D FT-ICR
                'Sampling_file_aposteriori_cytoCpnas.list','angio_ms_000005.d',
                'SubsP_220615_2DFT_2k_128k_000001.d', 'HEKmedia',
                'testsmzXML')
    for i in glob.glob(filename("*")):
        if verbose: print(i, end=' ')
        if os.path.basename(i) in files_to_keep:
            if verbose: print(" Ok")
        else:
            if CLEAN:
                try:
                    os.remove(i)
                    if verbose: print(" removed")
                except OSError:
                    if verbose: print(" **** could not be removed ****")
            else:
                if verbose: 
                    print(" should be removed")
                else:
                    print(i, " should be removed")

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
    import numpy
    import platform
    global list_of_modules
    python = "{0}.{1}.{2}".format(*sys.version_info)
    npv = numpy.version.version
    subject = "SPIKE tests performed on {2} {4} running python {0} / numpy {1} on host {3}".format(python, npv, *platform.uname())
    to_mail = [msg(subject), "Test program version %s"%__version__]

    # add spike prefix
    list_of_modules = ["spike."+mod for mod in list_of_modules]
    print(msg("modules to be tested are:"))
    for mod in list_of_modules:
        print(mod)

    stdata = "\nDatasets for tests are located in : %s"%DATA_dir
    to_mail.append(stdata)
    print(stdata)

    if not RUN:
        print ("\nDry run - stopping now, without executing the tests")
        sys.exit(0)
    print(msg("First removing leftover files"))
    cleandir()   # (CLEAN is handeld internally )
    if CLEAN:
        print(msg("removing .pyc in spike"))
        cleanspike()
    print(msg("Running automatic Tests"))
    Alltests = 0
    t0 = time.time()
    suite = unittest.defaultTestLoader.loadTestsFromNames( list_of_modules )
    nbtests = suite.countTestCases()
    Alltests += nbtests
    results = unittest.TextTestRunner(verbosity = 2).run(suite)
    elaps = time.time()-t0
    to_mail.append("Ran %d tests in %.3fs"%(nbtests, elaps))
    if results.wasSuccessful():
        if len(list_of_modules)==1:
            to_mail.append( msg("CONGRATULATIONS - The SPIKE tested module performed succesfully ") )
        else:
            to_mail.append( msg("CONGRATULATIONS - all the {} SPIKE tested modules performed succesfully ".format(len(list_of_modules))) )
        subject = "SUCCESS :)  " + subject
        to_mail.append( msg("modules tested were:"))
        for mod in list_of_modules:
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
    cleandir(verbose=False)
    # and report
    print("\n" + "="*80 + "\n")
    print(len(list_of_modules), 'SPIKE modules tested for a total of', Alltests, 'tests')
    if results.wasSuccessful():
        print("All tests successfull")
        print("performed on {2} {4} running python {0} / numpy {1} on host {3}".format(python, npv, *platform.uname()))
        print("test performed in %.2f sec"%elaps) 
    else:
        print(msg("Some tests Failed\nPlease revise results above"))
    print("\n" + "="*80 + "\n")

def main():
#    import Display.testplot as testplot
    import argparse
    global MAIL, list_of_mails, RUN, list_of_modules, DATA_dir, CLEAN
    global config
    
    # Parse and interpret options.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-m', '--mail', action='append', help="a mail address to which the result is sent at the end of the test - multiple entries possible -")
    parser.add_argument('-D', '--Data', help="the location of the folder containing files for the tests")
    parser.add_argument('-t', '--test_modules', action='append', help="overwrite the list of modules to test - multiple entries possible -")
    parser.add_argument('-n', '--dry',  action='store_true', help="list parameters and do not run the tests")
    parser.add_argument('-d', '--dirty', action='store_true', help="do not remove temporary files")
    parser.add_argument('-c', '--clean', action='store_true', help="just remove left-over temporary files")
    parser.add_argument('-g', '--graphic',  action='store_true', help="restore graphic output (off by default for background testing)")
    parser.add_argument('-W', '--Warnings',  action='store_true', help="Activate all warnings")
    
    args = parser.parse_args()

    print( "mail", args.mail)
    print( "data", args.Data)
    print( "module", args.test_modules)
    print("dry", args.dry)
    print('dirty', args.dirty)
    print('clean', args.clean)
    print('graphic', args.graphic)
    print('Warnings', args.Warnings)

    if args.Warnings:
# disable eventual warning blocking
        if not sys.warnoptions:
            warnings.simplefilter("error") # Change the filter in this process
            os.environ["PYTHONWARNINGS"] = "error" # Also affect subprocesses

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
    if args.graphic:
        testplot.PLOT = True
    else:
        testplot.PLOT = False   # switches off the display for automatic tests
    # Beware  !! Awful Hack !!
    # testplot is used to switch display off 
    # BUT ALSO to monkey patch at Tests run time to store in testplot.config values used by the tests (DATA_dir)
    testplot.config = {}
    testplot.config["DATA_dir"] = DATA_dir

#    print(sys.modules)
    if args.clean:
        CLEAN = True
        cleandir(verbose=True)
    else:
        do_Test()
    
if __name__ == '__main__':
    main()
