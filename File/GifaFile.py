#!/usr/bin/env python 
# encoding: utf-8

"""
GifaFile.py

Created by Marc-André on 2010-03-17.
Copyright (c) 2010 IGBMC. All rights reserved.


This module provides a simple access to NMR files in the Gifa format.


"""

from __future__ import print_function, division
import re
import numpy as np
from .. import NPKData as npkd
from ..NPKError import NPKError
import unittest
import os

import sys #
if sys.version_info[0] < 3:
    pass
else:
    xrange = range

HEADERSIZE = 4096   # the size of the ascii header
BLOCKIO = 4096      # the size of the chunk written to disk

__version__ = "0.3"
# ported to python 3 jan 2017
# MAD&MAC 12 12 2011 :  nasty bug in reading 2D which let one point to zero - introduced in previous modif
# MAD 13 07 2010 :   2D read in windows
#                   speed-up of read thanks to Numpy
# CR Adapted for Windows CR Easter 2010 
#    All modifications start  by a comment CR

########################################################################
class GifaFile(object):
    """
    defines the interface to simply (read/write) acces Gifa v4 files
    standard methods are load() and save()
    
    standard sequence to read is
    F = GifaFile(filename,"r")
    B = F.get_data()      # B is a NPKdata
    F.close()

    or
    F = GifaFile(filename,"r")
    F.load()
    B = F.data      # B is a NPKdata
    F.close()
    

    and to write
    F = GifaFile(filename,"w")
    F.set_data(B)         # where B is a NPKdata; do not use    F.data = B
    F.save()
    F.close()

    The file consists of a header (of size headersize) and data
    The header is handled as a dictionnary  self.header
    data is handled as a NPKdata    self.data
        so numpy ndarray are in                self.data.buffer
    """
    def __init__(self, fname, access="r", debug = 0):
        self.debug = debug
#        self.buffer = []
        if isinstance(fname, str):  # check if we have a string
            self.fname = fname
            if access not in ("w","r"):
                raise Exception(access + " : acces key not valid")
            self.file = open(fname, access) # open file at once  - CR change file to open Python 2.6
            if access == "r":
                self.fileB = open(fname, "rb") # CR "b" for Windows  - MAD and python3
        else:       # if not string, assume a File object
            try:
                self.fname = fname.name
            except AttributeError:
                self.fname = "internal_buffer"
            self.file = fname
            self.fileB = fname
        if access == "r":     # check if it is a real one
            l = self.fileB.readline(32)
            hsz = re.match('HeaderSize\s*=\s*(\d+)', l.decode())
            if not hsz:
                self.file.close()
                raise Exception("file %s not valid"%fname)
            self.headersize = int(hsz.group(1))
            self.fileB.seek(0)   # and rewind.
        if access == "w":
            self.headersize = HEADERSIZE
        self.header = None
        self.header_valid = False
        self.data = None
        self.data_valid = False
        
    #----------------------------------------------
    def report(self):
        "prints a little debugging report"
        print("Dim", self.dim)
        print("header_valid", self.header_valid)
        print("data_valid", self.data_valid)
        print("header", self.header)
        
    #----------------------------------------------
    def get_data(self):
        """returns the NPKdata attached to the (read) file"""
        if not self.data_valid:
            self.load()
            self.data_valid = True
        return self.data
    #----------------------------------------------
    def set_data(self, buff):
        """sets the NPKdata attached to the (to be written) file"""
        self.data = buff
        self.setup_header()
        self.data_valid = True
        self.header_valid = True
    #----------------------------------------------
    def copyaxesfromheader(self, n_axis):
        """
        get values from axis "n_axis" from header, and creates and returns a new (NMRAxis) axis with this values
        itype is not handled (not coded per axis in header)
        used internally
        """
        axis = npkd.NMRAxis()
        axis.size = int(self.header["Dim%d"%n_axis])
        axis.specwidth = float(self.header["Specw%d"%n_axis])
        axis.offset = float(self.header["Offset%d"%n_axis])
        axis.frequency = float(self.header["Freq%d"%n_axis])
        return axis
    #----------------------------------------------
    def copydiffaxesfromheader(self):
        """
        get values from axis "n" from header, and creates and returns a new (LaplaceAxis) axis with this values
        used internally
        """
        axis = npkd.LaplaceAxis()
        axis.size = int(self.header["Dim1"])
        if (self.header["Dmin"] != "NaN"):
            axis.dmin = float(self.header["Dmin"])
            axis.dmax = float(self.header["Dmax"])            
            axis.dfactor = float(self.header["Dfactor"])
        else:
        #except KeyError:    # if one is missing, ressort to default values
            axis.dmin = 1.0
            axis.dmax = 10.0
            axis.dfactor = 1.0
        return axis
    #----------------------------------------------
    def load(self):
        """creates a NPKdata loaded with the file content"""
        if not self.header_valid:
            self.load_header()
        if not self.data_valid:
#            ndata = NPKdata(dim =self.dim)
#            ndata.data = self.readc().data               # modified buffer to data
            ndata = npkd.NPKData(buffer=self.readc(), dim =self.dim)
            #ndata.display(scale = 10.0,show=True)
            # setup axes
            ndata.axis1 = self.copyaxesfromheader(1)
            if ndata.dim >= 2:
                ndata.axis2 = self.copyaxesfromheader(2)
            if ndata.dim >= 3:
                ndata.axis3 = self.copyaxesfromheader(3)
            # setp itype
            if ndata.dim == 1:
                ndata.axis1.itype = self.itype
            elif ndata.dim == 2:
                ndata.axis1.itype = self.itype//2
                ndata.axis2.itype = self.itype%2
            if ndata.dim == 3:
                ndata.axis1.itype = self.itype//4
                ndata.axis2.itype = (self.itype//2)%2
                ndata.axis3.itype = self.itype%4
            ndata.diffaxis = self.copydiffaxesfromheader()
            self.data = ndata
            self.data_valid = True
    #----------------------------------------------
    def save(self):
        """save the NPKdata to the file"""
        if not self.header_valid:
            raise Exception("header not set")
        if not self.data_valid:
            raise Exception("buffer data not set")
#        self.file = sys.stdout
        self.file = open(self.fname,'w')
        self.fileB = open(self.fname,'wb') # CR for Windows
        self.writec()
        self.close()
    #----------------------------------------------
    def read_header(self):
        """
        return a dictionnary of the file header
        internal use
        """
        self.fileB.seek(0)   # rewind.
        buf = self.fileB.read(self.headersize).decode()
        self.fileB.seek(0)   # rewind.
        dic = {}
        for line in buf.split("\n"): # go through lines in buf
            lsp = re.split(r'(?<!\\)=', line, 1)   #matches = but not
            if len(lsp) == 2:
                dkey = re.sub(r'\\=', '=', lsp[0])
                fval = re.sub(r'\\=', '=', lsp[1])      # replaces \= by =
                dic[dkey.strip()] = fval.strip()
        return dic
    #----------------------------------------------
    def write_header_line(self, key):
        """
        write into the header the entry key
        returns the number of byte written
        internal use
        """
        l = ("%-12s = %s\n")%(key, self.header[key])
        self.fileB.write(l.encode())    #  CR write in binary mode to presserve the UNIX EOL character - MAD now general in py3
        if self.debug > 0 : print(l, end=' ')
        return len(l)
    #----------------------------------------------
    def setup_header(self):
        """setup file header, from self.data"""
        h = {}
        h["Cacheversion"] = "2"
        h["Cacherelease"] = "0"
        h["Byteorder"] = "big_endian"   # default for Intel like endianness
        h["Dim"] = "%d"%self.data.dim
        h["Freq"] = "%f"%self.data.frequency
        it = 0  # will accumulate itype value
        for ax in range(self.data.dim):
            try:
                axis = self.data.axes(ax+1)
                """
                Here we struggled a bit on how to compute the general it. 
                We have :
                        1
                1D    |___|
                        \  \    axis 1 become 2 in 2D
                        1    2   
                2D    |___||___|  
                        \  \\   \   axis 2 becomes 3 in 3D, axis 1 becomes 2 
                        1    2    3
                3D    |___||___||___|
                
                We need itype between 0-1 in 1D
                                      0-2 in 2D (axis1)
                                      0-1 in 2D (axis2)
                                      0-4 in 3D (axis1)
                                      0-2 in 3D (axis2)
                                      0-1 in 3D (axis3)
                                      
                we came out with this line:
                it = it + axis.itype*(2**(self.data.dim-ax-1))
                """
                it = it + axis.itype*(2**(self.data.dim-ax-1))
                if self.debug>0: print("ICI",ax,axis.itype*(2**(self.data.dim-ax-1)),it)
            except:
                print("we don't have data axis")
                axis = None
            if axis:
                for (key, param, def_val) in (  ("Dim","size",64),
                                                ("Offset","offset",0.0),
                                                ("Freq","frequency",1.0),
                                                ("Specw","specwidth",1.0) ):
                    val = getattr(axis, param, def_val)   # get axis values, if exists, otherwise def_value
                    h["%s%d"%(key,ax+1)] = val
        
        try:
            h["Dmin"] = "%f"%self.data.diffaxis.dmin
            h["Dmax"] = "%f"%self.data.diffaxis.dmax
            h["Dfactor"] = "%f"%self.data.diffaxis.dfactor
            
        except:
            h["Dmin"] = 1.0
            h["Dmax"] = 10.0
            h["Dfactor"] = 1.0
        # try:
        #     h["Offset1"] = "%f"%self.data.axis1.offset
        # except:
        #     h["Offset1"] = 1.0
        # try:
        #     h["Freq1"] = "%f"%self.data.axis1.frequency
        # except:
        #     h["Freq1"] = 1.0
        # # set-up axes
        # if self.data.dim > 1:
        #     h["Dim2"] = "%d"%self.data.axis2.size
        #     h["Specw2"] = "%f"%self.data.axis2.specwidth
        #     try:
        #         h["Offset2"] = "%f"%self.data.axis2.offset
        #     except:
        #         h["Offset2"] = 1.0
        #     try:
        #         h["Freq2"] = "%f"%self.data.axis2.frequency
        #     except:
        #         h["Freq2"] =  1.0
        #     it = it + 2*self.data.axis2.itype
        # if self.data.dim > 2:
        #     h["Dim3"] = "%d"%self.data.axis3.size
        #     h["Specw3"] = "%f"%self.data.axis3.specwidth
        #     h["Offset3"] = "%f"%self.data.axis3.offset
        #     h["Freq3"] = "%f"%self.data.axis3.frequency
        #     it = it + 4*self.data.axis3.itype
        h["Type"] = "%d"%it
        # compute submatrix blocks
        if self.data.dim == 1:
            h["Szblk1"] = "1024"
            h["Nbblock1"] = "%d"%(self.data.axis1.size//1024)
        elif self.data.dim == 2:
            sz12 = float(self.data.axis2.size) / self.data.axis1.size
            n2 = 1
            n1 = BLOCKIO // 4
            while  ( (float(n2)// n1) < sz12) and (n1 > 1 ) :
                n1 = n1//2
                n2 = n2*2
            if self.debug > 0:
                print("si1 x si2 : %d %d   n1 x n2 : %d %d"%(self.data.axis1.size,self.data.axis2.size,n1,n2))
            h["Szblk1"] = "%d"%n1
            h["Szblk2"] = "%d"%n2
            h["Nbblock1"] = "%d"%(1+(self.data.axis1.size-1)//n1)
            h["Nbblock2"] = "%d"%(1+(self.data.axis2.size-1)//n2)
            
        elif self.data.dim == 3:
            sz12 = float(self.data.axis2.size*self.data.axis3.size) // (self.data.axis1.size*self.data.axis1.size)
            n1 = BLOCKIO // 4
            n2 = 1
            n3 = 1
            while  ((float(n2*n3)/(n1*n1)) < sz12 ) and (n1 > 1):
                n1 = n1//2
                n2 = n2*2
            sz12 = float(self.data.axis3.size) / (self.data.axis1.size)
            while  ((float(n3)/ n2) < sz12 ) and (n2 > 1):
                n2 = n2//2
                n3 = n3*2
            if self.debug > 0:
                print("si1 x si2 x si3: %d %d %d   n1 x n2 x n3 : %d %d %d"%(self.data.axis1.size, self.data.axis2.size, self.data.axis3.size, n1, n2, n3))
            h["Szblk1"] = "%d"%n1
            h["Szblk2"] = "%d"%n2
            h["Szblk3"] = "%d"%n3
            h["Nbblock1"] = "%d"%(1+(self.data.axis1.size-1)//n1)
            h["Nbblock2"] = "%d"%(1+(self.data.axis2.size-1)//n2)
            h["Nbblock3"] = "%d"%(1+(self.data.axis3.size-1)//n3)
        self.header = h
        
    #----------------------------------------------
    def write_header(self):
        """
        write file header
        setup_header() should have been called first
        """
        # write first line
        # CR write in binary mode to presserve the UNIX EOL character 
        self.fileB.seek(0) #CR not neccesary, better to be carefull
        len_so_far = 0
        l = "HeaderSize   = %d\n"%HEADERSIZE
        self.fileB.write(l.encode())
        len_so_far = len_so_far+len(l)
        # then the other
        for k in self.header.keys():
            len_so_far = len_so_far + self.write_header_line(k)
        # then flush buffer up to Headersize
        self.fileB.write( b"0"*(HEADERSIZE-len_so_far) )
        
    #----------------------------------------------
    def load_header(self):
        """load the header from file and set-up every thing"""
        if not self.header_valid:
            self.header = self.read_header()
            self.header_valid = True
    #----------------------------------------------
    def close(self):
        """ closes the associated file"""
        self.file.close()
        self.fileB.close()
    #----------------------------------------------
    # properties
    @property
    def dim(self):
        """dimensionality of the dataset 1 2 or 3"""
        return int(self.header["Dim"])
    @property
    def size1(self):
        """size along the F1 axis (either 1D, or slowest varyong axis in nD)"""
        return int(self.header["Dim1"])
    @property
    def size2(self):
        """size along the F2 axis (fastest varying in 2D)"""
        return int(self.header["Dim2"])
    @property
    def size3(self):
        """size along the F3 axis (fastest varying in 3D)"""
        return int(self.header["Dim3"])
    @property
    def szblock1(self):
        """size of data block on disk along F1 axis"""
        return int(self.header["Szblk1"])
    @property
    def szblock2(self):
        """size of data block on disk along F2 axis"""
        return int(self.header["Szblk2"])
    @property
    def szblock3(self):
        """size of data block on disk along F3 axis"""
        return int(self.header["Szblk3"])
    @property
    def nblock1(self):
        """number of data block on disk along F1 axis"""
        return int(self.header["Nbblock1"])
    @property
    def nblock2(self):
        """number of data block on disk along F2 axis"""
        return int(self.header["Nbblock2"])
    @property
    def nblock3(self):
        """number of data block on disk along F3 axis"""
        return int(self.header["Nbblock3"])
    @property
    def itype(self):
        """
        Real/complex type of the dataset
        in 1D :     0 : real  1: complex
        in 2D :     0 : real on both; 
                    1 : complex on F2
                    2 : complex on F1
                    3 : complex on both
        in 3D :     0 : real on all; 
                    1 : complex on F3
                    2 : complex on F2
                    3 : complex on F3-F2
                    4 : complex on F1
                    5 : complex on F1-F3
                    6 : complex on F1-F2
                    7 : complex on all
        """
        return int(self.header["Type"])
    @property
    def byte_order(self):
        """pour intel"""
        try:
            if self.header["Byteorder"] == "big_endian":
                r = False
            else:
                r = True
        except KeyError:
            r = True
        return r
    #----------------------------------------------
    def readc(self):
        """
        read a file in Gifa format, and returns the binary buffer as a numpy array
        internal use - use load()
        """
        import array
        self.load_header()
        self.fileB.seek(self.headersize) #CR binary handler for data
        if self.dim == 1:
            print("loading 1D")
            sz = self.size1
            fbuf = self.fileB.read(4*sz)   # Gifa data are in 4 byte floats CR binary handler for data
            abuf = array.array('f', fbuf)    # simple precision 4 bytes
            if self.byte_order : abuf.byteswap()
    #        self.fbuf = np.ndarray( (sz,), dtype='f', data = abuf)   # copy without casting
            fbuf = np.empty( (sz,), dtype='float_')  # double precision
            fbuf[:] = abuf[:]
            # for i in xrange(sz):      # this is 3-10 times slower
            #     fbuf[i] = abuf[i]
        elif self.dim == 2:
            print("loading 2D")
            sz1 = self.size1
            sz2 = self.size2
            if self.debug > 0: print("2D", sz1, sz2)
            fbuf = np.empty( (sz1, sz2))
#            fbuf = -1000.0*np.ones( (sz1, sz2))
            i1 = 0
            i2 = 0
            if self.debug > 0: 
                print("sz", self.szblock1, self.szblock2)
                print("nb", self.nblock1, self.nblock2)
            for b1 in xrange(self.nblock1):
                for b2 in xrange(self.nblock2):
                    #print b1,b2,i1,i2
                    tbuf = self.fileB.read(4*self.szblock1*self.szblock2)   # Gifa buffer are in 4 byte floats
                    abuf = array.array('f', tbuf)
                    if self.byte_order: abuf.byteswap()
                    imax = min(i1+self.szblock1, sz1)-i1
                    for i in xrange(imax):
                        jmax = min(i2+self.szblock2, sz2)-i2
                        fbuf[i1+i, i2:i2+jmax] = abuf[i*self.szblock2:i*self.szblock2+jmax]
#                         jmax = min(i2+self.szblock2, sz2)-i2        # this is 3-10 times slower
#                         for j in xrange(jmax):
# #                            print i1+i,i2+j
#                             fbuf[i1+i,i2+j] = abuf[i*self.szblock2+j]
                    i2 = i2+self.szblock2
                i2 = 0
                i1 = i1+self.szblock1
        elif self.dim == 3:
            print("reading 3D")
            print(" A VERIFIER")
            sz1 = self.size1
            sz2 = self.size2
            sz3 = self.size3
            fbuf = np.empty( (sz1, sz2, sz3))
#            fbuf = -1000.0*np.ones( (sz1, sz2))
            print("3D:", sz1, sz2, sz3)
            i1 = 0
            i2 = 0
            i3 = 0
            if self.debug > 0:
                print("3D:", sz1, sz2, sz3)
            for b1 in xrange(self.nblock1):
                for b2 in xrange(self.nblock2):
                    for b3 in xrange(self.nblock3):
                        tbuf = self.fileB.read(4*self.szblock1*self.szblock2*self.szblock3)   # Gifa buffer are in 4 byte floats
                        abuf = array.array('f', tbuf)
                        if self.byte_order: abuf.byteswap()
                        imax = min(i1+self.szblock1, sz1)-i1
                        for i in xrange(imax):
                            jmax = min(i2+self.szblock2, sz2)-i2
                            for j in xrange(jmax):
                                kmax = min(i3+self.szblock3, sz3)-i3
                                fbuf[i1+i,i2+i,i3:i3+kmax] = abuf[i*self.szblock3:i*self.szblock3+kmax]
                        i3 = i3+self.szblock3
                    i3 = 0
                    i1 = i1+self.szblock1
                    i2 = i2+self.szblock2
        return fbuf
    #----------------------------------------------
    def writec(self):
        """
        write a file in Gifa format
        internal use - use save()
        """
#        import array
        self.write_header()
        self.fileB.seek(self.headersize) #CR binary handler for data
        if self.dim == 1:
            self.fileB.write( self.data.buffer.astype("float32").tostring() ) #CR binary handler for data
        elif self.dim == 2:
            print("writing 2D")
            sz1 = self.size1
            sz2 = self.size2
            if self.debug > 0:
                print("2D:", sz1, sz2)
            i1 = 0
            i2 = 0
            fbuf = np.zeros((BLOCKIO//4,), dtype='float32')
            for b1 in xrange(self.nblock1):
                for b2 in xrange(self.nblock2):
    #                print b1,b2,i1,i2
                    imax = min(i1+self.szblock1, sz1)-i1
                    for i in xrange(imax):
                        jmax = min(i2+self.szblock2, sz2)-i2
                        #print fbuf[i*self.szblock2:i*self.szblock2+jmax] 
                        fbuf[i*self.szblock2:i*self.szblock2+jmax] = self.data.buffer[i1+i,i2:i2+jmax]
#                        for j in xrange(self.szblock2):
#                            fbuf[i*self.szblock2+j] = self.data.buffer[i1+i,i2+j]
                    self.fileB.write( fbuf.tostring() ) #CR binary handler for data
                    i2 = i2+self.szblock2
                i2 = 0
                i1 = i1+self.szblock1
        elif self.dim == 3:
            print("writing 3D")
            print(" A VERIFIER")
            sz1 = self.size1
            sz2 = self.size2
            sz3 = self.size3
            if self.debug > 0:
                print("3D:", sz1, sz2, sz3)
            i1 = 0
            i2 = 0
            i3 = 0
            fbuf = np.zeros((BLOCKIO//4,), dtype='float32')
            print(self.nblock1, self.nblock2, self.nblock3)
            for b1 in xrange(self.nblock1):
                for b2 in xrange(self.nblock2):
                    for b3 in xrange(self.nblock3):
                        imax = min(i1+self.szblock1, sz1)-i1
                        for i in xrange(imax):
                            jmax = min(i2+self.szblock2, sz2)-i2
                            for j in xrange(jmax):
                                kmax = min(i3+self.szblock3, sz3)-i3
                                fbuf[i*self.szblock2+j*self.szblock3:i*self.szblock2+j*self.szblock3+kmax] = self.data.buffer[i1+i, i2+1, i3:i3+kmax]
        #                        for j in xrange(self.szblock2):
        #                            fbuf[i*self.szblock2+j] = self.data.buffer[i1+i,i2+j]
                        self.fileB.write( fbuf.tostring() ) #CR binary handler for data
                        i3 = i3+self.szblock3
                    i3 = 0
                    i1 = i1+self.szblock1
                    i2 = i2+self.szblock2
class GifaFileTests(unittest.TestCase):
    "  - Testing GifaFile on various 1D and 2D files - "
    # basic files

    verbose = 1    # verbose > 0 switches messages on
    def announce(self):
        if self.verbose >0:
            print("\n========",self.shortDescription(),'===============')
    def test_read(self):
        """ - testing read capacities - """
        from ..Tests import filename, directory
        self.announce()
        # 1D
        name1D = filename("proj.gs1")
        name2D = filename("dosy-cluster2.gs2")       # Byteorder = big_endian
        name2D_little_endian = filename("dosy-cluster2-corr.gs2")   # Byteorder = little_endian
        G = GifaFile(name1D,"r")
        # show header
        G.load_header() # load header
        self.assertEqual(G.dim,1)
        self.assertEqual(G.header["Spec_date"],"2006-05-06")
        G.load()
        B = G.get_data()
        self.assertAlmostEqual(B.buffer[0], 1869.4309082)
        self.assertAlmostEqual(B.buffer.max(), 603306.75)
        G.close()
        # 2D
        G = GifaFile(name2D,"r")
        # show header
        G.load_header() # load header
        self.assertEqual(G.dim,2)
        G.load()
        B = G.get_data()
        # print B.buffer[0,0]
        # print B.buffer[133,1101]
        # print B.buffer.max()
        # B.display(show=True)
        G.close()
        self.assertEqual(B.buffer[0,0], 0.0)
        self.assertAlmostEqual(B.buffer[133,1101], 5164615.5)
        self.assertAlmostEqual(B.buffer.max(), 6831767.0)
        # 2D little_endian
        G = GifaFile(name2D_little_endian,"r")
        # show header
        G.load_header() # load header
        self.assertEqual(G.dim,2)
        G.load()
        B = G.get_data()
        # print B.buffer[0,0]
        # print B.buffer[133,1101]
        # print B.buffer.max()
        # B.display(show=True)
        G.close()
        self.assertEqual(B.buffer[0,0], 0.0)
        self.assertAlmostEqual(B.buffer[133,1101], 5164615.5)
        self.assertAlmostEqual(B.buffer.max(), 6831767.0)

    def test_write1D(self):
        """ - test 1D write capacities -"""
        import os
        from ..Tests import filename
        self.announce()
        nameout = filename("test_write.gs1")
        G = GifaFile(nameout,"w")
        # generate data
        x = np.arange(1024)/1000.   # 1000 Hz on 1024 points
        fid = np.zeros_like(x)
        LB = 5.0
        for i in range(1, 6):
            fid = fid + i*20*np.exp(2*i*432.1j*x)*np.exp(-LB*x)   # that's 5 lines altogether
        # put it into a NPKBuffer
        print("***", fid[10])
        B = npkd.NPKData(buffer=fid)
        B.axis1 = npkd.NMRAxis(size=2*1024, specwidth=1000, offset=0.0, frequency=400.0, itype=1)
        # and save it
        G.set_data(B)
        G.save()
        G.close()
        G2 = GifaFile(nameout,"r")
        B2 = G2.get_data()
        G2.close()
#        os.unlink(nameout)
        print("===============\n",B2.report())
        self.assertEqual(B2.axis1.itype, 1)
        self.assertAlmostEqual(B2.buffer[20], 18.7938525625, places=5)  # places=5 because GifaFile are single precision !
        self.assertAlmostEqual(B2.buffer[21], -51.1309912819, places=5)

    def test_write2D(self):
        """ - testing 2D read/write capacities - """
        from ..Tests import filename
        self.announce()
        # first read
        #G = GifaFile(self.name2D,"r")
        G = GifaFile(filename("dosy-cluster2.gs2"), "r")
        G.load()        # load dataset
        A = G.get_data()
        A.buffer *= 3   # modify trivially
        print(type(A))
        G.close()
        del(G)
        # then save it
        nameout = filename("test_write2D2.gs2")
        #f.close()
        H = GifaFile(nameout,"w")
        H.set_data(A)
        H.save()        # save dataset
        H.close()
        # finally read back
        GG = GifaFile(nameout,"r")
        GG.load()        # load dataset
        GG.close()
        B = GG.get_data()
        os.unlink(nameout)
        #self.assertEqual(B.axis1.itype, 0)
        # self.assertEqual(B.buffer[0,0], 0.0)
        # self.assertAlmostEqual(B.buffer[133,1101], 15493846.0)
        # self.assertAlmostEqual(B.buffer.max(), 20495300.0)
        # self.assertAlmostEqual(B.buffer.max()/A.buffer.max(), 1.0)
        # self.assertAlmostEqual(B.buffer.min(),A.buffer.min())
    def base(self):
        "test basic function"
        from ..Tests import filename
        nameout = filename("toto.gs2")
        try:
            os.unlink(nameout)
        except:
            pass
        dd = 2*np.ones((508,2*1000)) # + np.arange(2*1000)
        print(dd.shape)
        H = GifaFile(nameout,"w")
        A = npkd.NPKData(buffer=dd)
        H.set_data(A)
        H.save()        # save dataset
        H.close()
        # finally read back
        GG = GifaFile(nameout,"r")
        GG.debug=2
        GG.load()        # load dataset
        GG.close()
        B = GG.get_data()
        print(A.buffer.min(), B.buffer.min(), B.buffer.argmin())
        print(A.buffer.max(), B.buffer.max())
        print(A.buffer.shape, B.buffer.shape)
        os.unlink(nameout)
    
if __name__ == '__main__':
    unittest.main()
