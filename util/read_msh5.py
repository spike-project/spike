
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
import File.HDF5File as HDF5File
#import FTICR

class read_msh5:
    '''
    resmin  minimal resolution
    resi(i) resoluton nb i
    '''
    def __init__(self,filename):
        self.filename = HDF5File.HDF5File(filename)
        #print dir(self.filename)
    @property
    def resmin(self): # minimal resolution
        for gp in self.filename.hf.walkGroups("/"):  
            resolname = gp._v_name #
        print "d is the FTICRDATA for the minimal resolution"
        return self.filename.get_data(resolname) # smallest resolution
        
    def resi(self,i): # resolutiono nb i
        for nbgp, gp in enumerate(self.filename.hf.walkGroups("/")):
            if nbgp == i :   
                resolname = gp._v_name #
        return self.filename.get_data(resolname) # smallest resolution


class read_msh5Tests(unittest.TestCase):
    def test_msh5read(self):
        RM = read_msh5('../../DATA_test/1D_test.msh5')
        d = RM.resmin()
        
if __name__ == '__main__':
    unittest.main()
    