'''
Created by Lionel Chiron  03/12/2013
Copyright (c) 2013 __NMRTEC__. All rights reserved.
'''
from NPKConfigParser import NPKConfigParser
import re, os, sys

class Proc_Parameters(object):
    """this class is a container for applying Fista"""
    def __init__(self, configfile = None):
        "initialisation"
        # 
        self.rank = 50                              # urQRd rank
        self.iterations = 1                         # number of urQRd iterations
        self.zerofill = 2                           # zerofilling for Fourier processing.
        self.save_ref_spec = False   #
        self.save_urqrd_spec = False
        self.save_urqrd_fid = True
        self.save_date = True
        self.mp = True
           
        if configfile:
            self.load(configfile)
            
    def load(self, cp):
        "load from cp config file - should have been opened with ConfigParser() first"
        ###
        self.npk_dir =  cp.get( "config", "npk_dir")        # directory of NPKv2
        ###
        self.file_to_treat =  cp.get( "data", "file_to_treat")        # input file
        self.filename, self.ext = os.path.splitext(self.file_to_treat)
        ###
        self.rank = cp.getint( "urqrd", "rank", self.rank) #
        self.iterations = cp.getint( "urqrd", "iterations", self.iterations)
        self.zerofill = cp.getint( "urqrd", "zerofill", self.zerofill)
        self.mp = cp.getboolean( "urqrd", "mp", str(self.mp))   # multiprocessing
        ###
        self.save_ref_spec = cp.getboolean( "store", "save_ref_spec", str(self.save_ref_spec))   #
        self.save_urqrd_spec = cp.getboolean( "store", "save_urqrd_spec", str(self.save_urqrd_spec))
        self.save_urqrd_fid = cp.getboolean( "store", "save_urqrd_fid", str(self.save_urqrd_fid))#
        self.save_date = cp.getboolean( "store", "save_date", str(self.save_date))   #
        self.addr_data_saved =  cp.get( "store", "addr_data_saved")  # address for saving data.
     
    def report(self):
        "print a formatted report"
        print "------------ processing parameters ------------------"
        for i in dir(self):
            if not i.startswith('_'):
                v = getattr(self,i)
                if not callable(v):
                    print i, ' :', v
        print "-----------------------------------------------------"
        
def config_urqrd(configfile = "Applic_urqrd.mscf"):
    print "using %s as configuration file"%configfile

    #### get parameters from configuration file - store them in a parameter object
    cp = NPKConfigParser()
    cp.readfp(open(configfile))
    param = Proc_Parameters(cp)
    return param