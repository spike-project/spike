#!/usr/bin/env python 
# encoding: utf-8

"""
    Utility to Handle Apex files
"""

from __future__ import print_function, division

__author__ = "Marc André Delsuc, Marie-Aude Coutouly <mac@nmrtec.com>"
__date__ = "July 2011"

import sys
import os
import unittest
import glob
import os.path as op
import numpy as np
import tables
from ..FTICR import FTICRData
from . import HDF5File as hf
# tables.parameters.NODE_CACHE_SLOTS = 0
#tables.parameters.CHUNK_CACHE_SIZE = 0*1024*1024
# tables.parameters.METADATA_CACHE_SIZE  = 10*1024*1024

if sys.version_info[0] < 3:
    pass
else:
    xrange = range

def read_param(filename):
    """
        Open the given file and retrieve all parameters written initially for apexAcquisition.method
        NC is written when no value for value is found
        
        structure : <param><name>C_MsmsE</name><value>0.0</value></param>
        
        read_param returns  values in a dictionnary
    """
    from xml.dom import minidom
    
    xmldoc = minidom.parse(filename)
    
    x = xmldoc.documentElement
    pp = {}
    children = x.childNodes
    for child in children:
        if (child.nodeName == 'paramlist'):
            params = child.childNodes
            for param in params:
                if (param.nodeName == 'param'):
                    if 'name' in param.attributes.keys():  # some versions
                        k = param.attributes['name'].value
                    for element in param.childNodes:   
                        if element.nodeName == "name":     # other versions
                            k = element.firstChild.toxml()
                        elif element.nodeName == "value":
                            try:
                                v = element.firstChild.toxml()
                            except: 
                                v = "NC"
                    pp[k] = v
    return pp
#-----------------------------------------
def read_scan(filename):
    """
    Function that returns the number of scan that have been recorded
    It is used to see wether the number of recorded points correspond to the L_20 parameter
    """
    from xml.dom import minidom

    xmldoc = minidom.parse(filename)
    
    x = xmldoc.documentElement
    pp = {}
    children = x.childNodes
    count_scan = 0
    for child in children:
        if (child.nodeName == 'scan'):
            count_scan += 1
    return count_scan
#-----------------------------------------
def get_param(param,names, values):
    """
    From params, this function returns the  value of the given param
    """
    for i in xrange(len(names)):
        if names[i] == param:
            return values[i]
#-----------------------------------------
def locate_acquisition(folder):
    """
        From the given folder this function return the absolute path to the apexAcquisition.method file
        It should always be in a subfolder 
    """
    L = glob.glob(op.join(folder,"*","apexAcquisition.method"))
    if len(L)>1:
        raise Exception( "You have more than 1 apexAcquisition.method file in the %s folder, using the first one"%folder )
    elif len(L) == 0: 
        raise Exception( "You don't have any apexAcquisition.method file in the  %s folder, please double check the path"%folder )
        
    return L[0]
#-----------------------------------------
def Import_1D(inifolder,outfile=""):
    """
    Entry point to import 1D spectra
    It returns a FTICRData
    It writes a HDF5 file if an outfile is mentionned
    """
    import array

    if sys.maxsize  == 2**31-1:   # the flag used by array depends on architecture - here on 32biy
        flag = 'l'              # Apex files are in int32
    else:                       # here in 64bit
        flag = 'i'              # strange, but works here.
    if op.isfile(inifolder):
        folder = op.dirname(inifolder)
    elif op.isdir(inifolder):
        folder = inifolder
    elif not op.exists(inifolder):
        raise Exception("File does not exist: "+inifolder)
    else:
        raise Exception("File is undecipherable: "+inifolder)
    try:
        parfilename = locate_acquisition(folder)
    except:
        raise Exception('%s does not seem to be a valid Apex spectrum'%(inifolder,))
    parfilename = locate_acquisition(folder)
    params = read_param(parfilename)
    
    
    # Import parameters : size in F1 
    sizeF1 = int(params["TD"])
    
    if os.path.isfile(os.path.join(folder, "fid")):
        fname = os.path.join(folder, "fid")
    else:
        raise Exception("You are dealing with 2D data, you should use Import_2D")
    data = FTICRData( dim=1 )   # create dummy 1D
    data.axis1.size = sizeF1    # then set parameters
    data.axis1.calibA = float(params["ML1"])
    data.axis1.calibB = float(params["ML2"])
    data.axis1.calibC = float(params["ML3"])
    data.axis1.specwidth = float(params["SW_h"])
    data.axis1.highfreq = data.axis1.calibA/float(params["EXC_low"])  # these two are in m/z !
    data.axis1.lowfreq = data.axis1.calibA/float(params["EXC_hi"])
    data.axis1.highmass = float(params["MW_high"])
    data.axis1.left_point = 0
    data.axis1.offset = 0.0
    
    data.params = params   # add the parameters to the data-set
        
    if (outfile): # Creates the fticrdata with no array, just infos for the table
        HF = hf.HDF5File(outfile,"w")
        HF.create_from_template(data)
        HF.store_internal_object(params, h5name='params')    # store params in the file
        # then store files xx.methods and scan.xml
        HF.store_internal_file(parfilename)
        HF.store_internal_file( os.path.join(folder,"scan.xml") )
#        HF.store_internal_file( os.path.join(folder,"ExciteSweep") )
        data.hdf5file = HF
        # I need a link back to the file in order to close it, however this creates a loop - experimental ! 
    else:
        data.buffer = np.zeros((sizeF1))
    data.adapt_size()
    with open(fname,"rb") as f:
        tbuf = f.read(4*sizeF1)
        abuf = np.array(array.array(flag,tbuf))
        data.buffer[:] = abuf[:]
    if (outfile):
        HF.flush()
    return data
#-----------------------------------------

def Import_2D(folder,outfile = "",F1specwidth = None):
    """
    Entry point to import 2D spectra
    It returns a FTICRData
    It writes a HDF5 file if an outfile is mentionned
    """
    import array

    if sys.maxsize  == 2**31-1:   # the flag used by array depends on architecture - here on 32bit
        flag = 'l'              # Apex files are in int32
    else:                       # here in 64bit
        flag = 'i'              # strange, but works here.
    parfilename = locate_acquisition(folder)
    params = read_param(parfilename)
    #count_scan = read_scan(os.path.join(folder,"scan.xml"))
    # Import parameters : size in F1 and F2
    sizeF1 = int(params["L_20"])
    # if (count_scan != sizeF1):
    #     print "Mismatch between Apex file and scan.xml size F1 set to %i"%count_scan
    #     sizeF1 = count_scan
    sizeF2 = int(params["TD"])

    if os.path.isfile(os.path.join(folder,"ser")):
        fname = os.path.join(folder, "ser")
    else:
        raise Exception("You are dealing with 1D data, you should use Import_1D")

    #size, specwidth,  offset, left_point, highmass, calibA, calibB, calibC, lowfreq, highfreq
    data = FTICRData( dim=2 )   # create dummy 2D

    data.axis2.size = sizeF2    # then set parameters along F2 - classical axis -
    data.axis2.calibA = float(params["ML1"])
    data.axis2.calibB = float(params["ML2"])
    data.axis2.calibC = float(params["ML3"])
    data.axis2.specwidth = float(params["SW_h"]) 
    data.axis2.highfreq = data.axis2.calibA/float(params["EXC_low"])  # these two are in m/z !
    data.axis2.lowfreq = data.axis2.calibA/float(params["EXC_hi"])
    data.axis2.highmass = float(params["MW_high"])
    data.axis2.left_point = 0
    data.axis2.offset = 0.0

    data.axis1.size = sizeF1    # then set parameters along F1- non classical axis -  
    # assumes most parameter are equivalent - except specwidth
    if F1specwidth is not None:     # if given in arguments
        data.axis1.specwidth = F1specwidth 
    else:
        f1 =  float(params["IN_26"])    # IN_26 is used in 2D sequence as incremental time
        if f1 < 1E-4:   # seems legit
            data.axis1.specwidth = 1.0/(2*f1)
        else:
            data.axis1.specwidth = data.axis2.specwidth     # else assume square...
    data.axis1.calibA = float(params["ML1"])
    data.axis1.calibB = float(params["ML2"])
    data.axis1.calibC = float(params["ML3"])
    data.axis1.highfreq = data.axis1.calibA/float(params["EXC_low"])  # these two are in m/z !
    data.axis1.lowfreq = data.axis1.calibA/float(params["EXC_hi"])
    data.axis1.highmass = float(params["MW_high"])
    data.axis1.left_point = 0
    data.axis1.offset = 0.0

    data.params = params   # add the parameters to the data-set

    c1 = int(sizeF1//8. +1)
    c2 = int(sizeF2//8. +1)
    #print "chunkshape",c1,c2
    #tables.parameters.CHUNK_CACHE_PREEMPT = 1
    #tables.parameters.CHUNK_CACHE_SIZE = c1*c2*8
    if (outfile): # Creates the fticrdata with no array, just infos for the table
        HF = hf.HDF5File(outfile,"w")
        HF.create_from_template(data)
        HF.store_internal_object(params, h5name='params')    # store params in the file
        # then store files xx.methods and scan.xml
        HF.store_internal_file(parfilename)
        HF.store_internal_file( os.path.join(folder,"scan.xml") )
#        HF.store_internal_file( os.path.join(folder,"ExciteSweep") )
        data.hdf5file = HF
        # I need a link back to the file in order to close it, however this creates a loop - experimental ! 
    else:
        data.buffer = np.zeros((sizeF1, sizeF2))
#    data.adapt_size()
    with open(fname,"rb") as f:
        for i1 in xrange(sizeF1-1):
            tbuf = f.read(4*sizeF2)
            abuf = np.array(array.array(flag,tbuf))
            data.buffer[i1,:] = abuf[:]
    # if (outfile):
    #     HF.close()
    if (outfile):
        HF.flush()
    return data
#-----------------------------------------
def read_2D(sizeF1, sizeF2, filename="ser"):
    """
    Reads in a Apex 2D fid

    sizeF1 is the number of fid
    sizeF2 is the number of data-points in the fid
    uses array
    """
    import array
#    import platform # platform seems to be buggy on MacOs, see http://stackoverflow.com/questions/1842544
    if sys.maxsize  == 2**31-1:   # the flag used by array depends on architecture - here on 32bit
        flag = 'l'              # Apex files are in int32
    else:                       # here in 64bit
        flag = 'i'              # strange, but works here.
    sz1, sz2 = int(sizeF1), int(sizeF2)
    # read binary
    fbuf = np.empty( (sz1, sz2) )
    with open(filename, "rb") as f:
        for i1 in xrange(sz1):
            tbuf = f.read(4*sz2)
            abuf = array.array(flag, tbuf)
            fbuf[i1,0:sz2] = abuf[0:sz2]
    return FTICRData(buffer = fbuf)
#-----------------------------------------
def read_3D(sizeF1, sizeF2, sizeF3, filename="ser"):  
    """
    Ebauche de fonction
    
    Reads in a Apex 3D fid


    uses array
    """
    import array
#    import platform # platform seems to be buggy on MacOs, see http://stackoverflow.com/questions/1842544
    if sys.maxsize  == 2**31-1:   # the flag used by array depends on architecture - here on 32biy
        flag = 'l'              # Apex files are in int32
    else:                       # here in 64bit
        flag = 'i'              # strange, but works here.
    sz1, sz2, sz3 = int(sizeF1), int(sizeF2), int(sizeF3)
    # read binary
    fbuf = np.empty( (sz1,sz2,sz3) )
    with open(filename,"rb") as f:
        for i1 in xrange(sz1):
            tbuf = f.read(4*sz2)
            abuf = array.array(flag, tbuf)
            fbuf[i1,0:sz2] = abuf[0:sz2]
    return FTICRData(buffer=fbuf)
#-----------------------------------------
def write_ser(bufferdata,filename="ser"):
    """
    Write a ser file from FTICRData
    """
    with open(filename, "wb") as f:
        for i in range (len(bufferdata)):
            for j in range (len(bufferdata[0])):
                
                f.write(bufferdata[i][j].astype("int32").tostring() )
#----------------------------------------------
class Apex_Tests(unittest.TestCase):
    def setUp(self):
        from ..Tests import filename, directory
        
        rootfiles = os.getcwd()
        self.TestFolder = directory()
        self.DataFolder = filename('ubiquitine_2D_000002.d')#'cytoC_2D_000001.d')
        # self.serfile = filename('cytoC_2D_000001.d/ser')
        # self.outHDF = filename('cytoC_2D_000001_direct.msh5')
        # self.name_write = filename("file_write")
        # self.name_fticr = filename("file_fticr")
        # self.name_npar = filename("file_npar")
        # self.npar_fticr = filename("npar_fticr")
        # self.name_chunk = filename("Chunk.hf")
        self.name_get = filename("file_fticr_2D.msh5")
        
        self.verbose = 1    # verbose > 0 switches messages on
    def announce(self):
        if self.verbose >0:
            print("\n========",self.shortDescription(),'===============')
    #-------------------------------------------------
    def test_Import_2D(self):
        "Test and time routine that import 2D from the MS-FTICR folder to the file given as second argument "
        from time import time
        self.announce()
        t0 = time()
        d = Import_2D(self.DataFolder, self.name_get)
        #d = Import_2D("/Volumes/XeonData/Developement/MS-FTICR/cytoC_2D_000006.d")
        print("import",time()-t0,"secondes")
        """
        sur mon ordi
        d = Import_2D("/DATA/FT-ICR-Cyto/cytoC_2D_000006.d")                prend 1080Mo de mémoire en 18 secondes
        d = Import_2D("/DATA/FT-ICR-Cyto/cytoC_2D_000006.d","test.hdf5")    prend 80Mo de mémoire en 21.6 secondes et créé un fichier de 1Go
        """
    #-------------------------------------------------
    def test_Import_2D_keep_mem(self):
        "Test and time routine that import 2D from the MS-FTICR folder in memory"
        from time import time
        self.announce()
        t0 = time()
        d = Import_2D( self.DataFolder )
        #d = Import_2D("/Volumes/XeonData/Developement/MS-FTICR/cytoC_2D_000006.d")
        print("import",time()-t0,"secondes")
        """
        sur mon ordi
        d = Import_2D("/DATA/FT-ICR-Cyto/cytoC_2D_000006.d")                prend 1080Mo de mémoire en 18 secondes
        d = Import_2D("/DATA/FT-ICR-Cyto/cytoC_2D_000006.d","test.hdf5")    prend 80Mo de mémoire en 21.6 secondes et créé un fichier de 1Go
        """
    #-------------------------------------------------
    def _test_2(self):
        " Test and time direct processing"
        self.announce()
        from time import time
        d = Import_2D(self.DataFolder,self.name2)
        t0 = time()
        t00 = t0
        d.rfft(axis=2)
        print("rfft2",time()-t0,"secondes")
        t0 = time()
        d.rfft(axis=1)
        print("rfft1",time()-t0,"secondes")
        t0 = time()
        d.modulus()
        print("modulus",time()-t0,"secondes")
        t0 = time()
        print("modulus",time()-t0,"secondes")
        print("calcul",time()-t00,"secondes")
        d.display(scale=5, show=True)
        d.hdf5file.close()
        """
        sur mon ordi - en mémoire
        rfft2 2.57822608948 secondes
        rfft1 9.06142997742 secondes
        modulus 0.851052999496 secondes
        modulus 2.92809915543 secondes
        calcul 15.4189431667 secondes

        - sur fichier -
        rfft2 7.82480788231 secondes
        rfft1 est incroyablement long > 5 minutes - et fait des I/O de folie

        Sur Xeon - en mémoire
        rfft2 4.27431201935 secondes
        rfft1 19.972933054 secondes
        modulus 1.25446605682 secondes
        modulus 5.96046447754e-06 secondes
        calcul 25.5018758774 secondes

        - sur fichier -
        """
    #-------------------------------------------------
    def _test_3(self):
        "Another strategy close the file after fft on F2 and reopen everything for F1 fft"
        self.announce()
        from time import time
        d = Import_2D(self.DataFolder,self.name3)
        t0 = time()
        t00 = t0
        d.rfft(axis=2)
        d.hdf5file.close()  # je ferme
        print("rfft2",time()-t0,"secondes")
        t0 = time()
        # je réouvre
        H = hf.HDF5File(self.name3,"rw")
        H.load()
        d2 = H.data      # B is a FTICRdata
        d2.rfft(axis=1)
        print("rfft1",time()-t0,"secondes")  # toujours aussi lent !

        t0 = time()
        d2.modulus()
        print("modulus",time()-t0,"secondes")
        t0 = time()
        print("modulus",time()-t0,"secondes")
        print("calcul",time()-t00,"secondes")
        print(type(d))
        d2.display(scale=5, show=True)
        H.close()
        H = hf.HDF5File(self.name3,"r")
        H.load()
        B = H.data      # B is a FTICRdata
        H.close()

        """
        Sur Xeon
        rfft2 19.9385159016 secondes
        rfft1 3348.65607595 secondes
        modulus 5.01055693626 secondes
        modulus 5.96046447754e-06 secondes
        calcul 3373.60540509 secondes

        en rfft(axis=2) 0.004 s par ligne --> 8.192 s
        en rfft(axis=1) 0.05 s par ligne --> 3276.8 s

        """
    # def test_1D(self):
    #     "Test the import_1D  routine"
    #     self.announce()
    #     d = Import_1D("/Volumes/XeonData/Developement/MS-FTICR/glyco_ms_000002.d","/Users/mac/Desktop/test1D_comp.hdf5")
if __name__ == '__main__':
    unittest.main()
