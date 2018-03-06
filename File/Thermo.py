#!/usr/bin/env python 
# encoding: utf-8

"""
    Utility to Handle Thermofisher files

    Marc-André from first draft by Lionel
"""

from __future__ import print_function

__author__ = "Marc André Delsuc"
__date__ = "april 2014"

import os
import unittest
import numpy as np
from ..Orbitrap import OrbiData
import re

def read_thermo(filename):
    """
    reads a thermofisher orbitrap file
    """
    with open(filename,'rb') as F:
        param = read_param(F)
        if param['Storage Type'] == 'float':
            typ = 'f4'
        elif param['Storage Type'] == 'double':
            typ = 'f8'
        elif param['Storage Type'] == 'int':
            typ = 'i4'
        else:
            raise Exception( "Unknown Storage type : " + param['Storage Type'] )
        # adapted spectralwidth
        data = read_data(F,typ)
        # swold = data.axis1.specwidth
        swnew = float(param['Bandwidth'])
        # m/z = A + B/f^2 + C/f^4
        data.axis1.calibA = float(param["Source Coeff1"])
        data.axis1.calibB = float(param["Source Coeff2"])*1E6     # Thermo works internally in kHz, we use Hz
        data.axis1.calibC = float(param["Source Coeff3"])*1E12
        data.axis1.specwidth = swnew
        # data.axis1.ref_freq *= swnew/swold
    return (param, data)
def read_param(F):
    """
        given F, an opend file , retrieve all parameters found in file header
        
        read_param returns  values in a plain dictionnary
    """
    dic = {}
    for l in F:
        if l.startswith("Data:"):
            break
        v = l.rstrip().split(':')       # remove trailing chars and split around :
        if len(v)<2:    # comment lines
            print(l)
        else:
            dic[v[0]] = v[1].lstrip()    # metadata lines
    return dic
def read_data(F, typ='float'):
    """
        given F, an opened file, reads the values and 
        read_param returns  values in a dictionnary
    """
    F.seek(0)
    pos = 0
    for l in F:
        pos += len(l)
        if l.startswith("Data Points"):
            print(re.findall(r'\d+', l)[0])
        if l.startswith("Data:"):
            break
    F.seek(pos)
    data_interm = np.array(np.fromfile(F, dtype=typ).tolist())  # ndarray to list to np.array
    data = OrbiData(buffer = data_interm )
    return data
#-----------------------------------------
def Import_1D(filename):
    """
    Entry point to import 1D spectra
    It returns an Orbitrap data
    """
    param, data = read_thermo(filename)
    data.params = param
    return data
#----------------------------------------------
class Thermo_Tests(unittest.TestCase):
    """ A FAIRE """
    def setUp(self):
        from ..Tests import filename, directory
        try:
            import ConfigParser
        except:
            import configparser as ConfigParser
        rootfiles = os.getcwd()        
        self.verbose = 1    # verbose > 0 switches messages on
    def announce(self):
        if self.verbose >0:
            print("\n========",self.shortDescription(),'===============')
    #-------------------------------------------------

if __name__ == '__main__':
    unittest.main()
    
