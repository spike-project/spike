# encoding: utf-8
"""
    Utility to Handle Thermofisher files

    Marc-André from first draft by Lionel
"""

__author__ = "Marc André Delsuc"
__date__ = "april 2014"

import os
import unittest
import numpy as np
from Orbitrap import OrbiData
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
        swold = data.axis1.specwidth
        swnew = float(param['Bandwidth'])
        data.axis1.specwidth = swnew
        data.axis1.ref_freq *= swnew/swold
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
            print l
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
            print re.findall(r'\d+', l)[0]
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
    return data
#----------------------------------------------
class Thermo_Tests(unittest.TestCase):
    """ A FAIRE"""
    def setUp(self):
        import ConfigParser
        rootfiles = os.getcwd()
        self.TestFolder = '../DATA_test'
        self.DataFolder = '../DATA_test/cytoC_2D_000001.d'
        self.serfile = '../DATA_test/cytoC_2D_000001.d/ser'
        self.outHDF = '../DATA_test/cytoC_2D_000001_direct.msh5'
        self.name_write = os.path.join(self.TestFolder,"file_write")
        self.name_fticr = os.path.join(self.TestFolder,"file_fticr")
        self.name_npar = os.path.join(self.TestFolder,"file_npar")
        self.npar_fticr = os.path.join(self.TestFolder,"npar_fticr")
        self.name_chunk = os.path.join(self.TestFolder,"Chunk.hf")
        self.name_get = os.path.join(self.TestFolder,"file_fticr")
        
        self.verbose = 1    # verbose > 0 switches messages on
    def announce(self):
        if self.verbose >0:
            print "\n========",self.shortDescription(),'==============='
    #-------------------------------------------------

if __name__ == '__main__':
    unittest.main()
