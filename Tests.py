# encoding: utf-8
"""
Tests.py

Created by Marc-André on 2010-07-20.
Copyright (c) 2010 IGBMC. All rights reserved.

Runs tests on selected modules using the integrated unittests. 
"""
import unittest
import os
from util.sendgmail import mail

# when CLEAN is set to true, unwanted files in test directory will be removed.
CLEAN = True #True

# Add your module here
mod_util = ('util.dynsubplot', 'util.debug_tools') #'util.read_msh5', 
mod_algo = ('Algo.Cadzow','Algo.Linpredic','Algo.urQRd',
                'Algo.SL0','Algo.maxent') 
 
mod_file = ("File.GifaFile",'File.HDF5File', 'File.Apex', 'File.csv')
mod_basicproc = ("NPKData", "FTICR", "Orbitrap")
mod_user = ('processing',)

list_of_modules =  mod_algo + mod_user + mod_basicproc + mod_file  + mod_util

list_of_mails =  ["madelsuc@unistra.fr", "lionel.chiron@gmail.com"]  # []

def msg(st, sep = '='):
    '''
    Message in Tests.py
    '''
    s = sep*(len(st) + 4)+"\n"
    s = s+ '| '+ st+ ' |'+"\n"
    s = s + sep*(len(st) + 4)+"\n"
    print s
    return s

def cleandraft():
    '''
    Removes the .pyc
    '''
    for root, dirs, files in os.walk('.'):
        #print root, dirs, files
        for f in files:
            r,ext = os.path.splitext(f)
            if ext == '.pyc':
                addr = os.path.join(root,f)
                print addr
                os.remove(addr)

def cleandir():
    "checking files in DATA_test directory and removes files created by previous tests"
    import glob
    DATA_dir = os.path.join('..', 'DATA_test')
    files_to_keep = ('ubiquitin_5_scan_res_30000_1.dat','cytoC_ms_1scan_000001.d', 'cytoC_2D_000001.d',
                'dosy-cluster2-corr.gs2', 'dosy-cluster2.gs2',
                'proj.gs1', 'ubiquitine_2D_000002.d','file_fticr',
                'file_npar','npar_fticr','Lasalocid-Tocsy', 'Sampling_file.list', 
                'ubiquitine_2D_000002_Sampling_2k.list',
                'Sampling_file_aposteriori_cytoCpnas.list','angio_ms_000005.d')
    for i in glob.glob(os.path.join(DATA_dir,"*")):
        print i,
        if os.path.basename(i) in files_to_keep:
            print " Ok"
        else:
            if CLEAN:
                try:
                    os.remove(i)
                    print " removed"
                except OSError:
                    print " **** could not be removed ****"
            else:
                print " should be removed"

class NPKTest(unittest.TestCase):
    """overload unittest.TestCase for default verbosity - Not Used - """
    def setUp(self):
            self.verbose = 1    # verbose > 0 switches messages on
    def announce(self):
        if self.verbose >0:
            print "\n========",self.shortDescription(),'==============='

def do_Test():
    '''
    Performs all tests then indicates if successfull.
    Gives total time elapsed.
    '''
    import time
    msg("modules to be tested are:")
    for mod in list_of_modules:
        print mod
    msg("First removing leftover files")
    cleandir()
    msg("removing .pyc in draft")
    cleandraft()
    msg("Running automatic Tests")
    t0 = time.time()
    suite = unittest.defaultTestLoader.loadTestsFromNames( list_of_modules )
    results = unittest.TextTestRunner(verbosity = 2).run(suite)
    elaps = time.time()-t0
    to_mail = []
    if results.wasSuccessful():
        to_mail.append( msg("CONGRATULATIONS - all the {} SPIKE tests performed succesfully ".format(len(list_of_modules))) )
        to_mail.append( msg("modules tested were:"))
        for mod in list_of_modules:
            print mod
            to_mail.append(mod)
        to_mail.append( msg("test performed in %.2f sec"%elaps) )
    else:
        to_mail.append( msg("Tests Failed, Please revise error codes", sep = '!') )
    for address in list_of_mails:
        mail(address, "msg from python Test.py in SPIKE", "\n".join(to_mail) )
        
if __name__ == '__main__':
    import Display.testplot as testplot
    testplot.PLOT = False   # switches off the display for automatic tests
    do_Test()
