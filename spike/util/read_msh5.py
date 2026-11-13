
# encoding: utf-8
"""
read_msh5.py
Small class that reads mr_MSH5 file and return a FTICRData
Has to be called with the full path to the mr_msh5 file
Created by mac on 2012-03-29.
Copyright (c) 2012 __NMRTEC__. All rights reserved.
"""

import sys
import os
import unittest
from .File.HDF5File import HDF5File

#import FTICR

class read_msh5:
    '''
    resmin method for recuperating the minimal resolution
    resi method for recuperating the resoluton nb i
    '''
    def __init__(self,filename):
        self.filename = HDF5File(filename)
        #print dir(self.filename)
        
    def resmin(self): # minimal resolution
        for gp in self.filename.hf.walk_groups("/"):  
            resolname = gp._v_name #
        print("d is the FTICRDATA for the minimal resolution")
        return self.filename.get_data(resolname) # smallest resolution

    def resi(self,i): # resolution nb i
        index = i
        for nbgp, gp in enumerate(self.filename.hf.walk_groups("/")):
            if gp._v_name == "attached":
                index = i+1
            if nbgp == index :  
                if "resol" in gp._v_name:
                    resolname = gp._v_name #
        return self.filename.get_data(resolname) # smallest resolution


class read_msh5Tests(unittest.TestCase):
    def setUp(self):
        pass


if __name__ == '__main__':
    RM = read_msh5(sys.argv[1])
    d = RM.resi(1)

# python -m spike.visu2D /Volumes/3ToBiak/MassSpec/FTICR_data/2D/cytoC_2D_000001_mr.msh5