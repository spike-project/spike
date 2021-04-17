#!/usr/bin/env python 
# encoding: utf-8

"""
NMR.py

This file defines generic classes for NMR Spectroscopy

Used to be inside NPKData
"""

from __future__ import print_function, division
import os
import math
import unittest
import numpy as np
from . import NPKData
from .NPKData import _NPKData, Axis, Unit, copyaxes
from .NPKError import NPKError


class NMRAxis(Axis):
    """
    hold information for one NMR axis
    used internally
    """
    def __init__(self,size = 64, specwidth = 2000.0*math.pi, offset = 0.0, frequency = 400.0, itype = 0, currentunit = "points"):
        """
        all parameters from Axis, plus
        specwidth   spectral width, in Hz
        offset      position in Hz of the rightmost point
        frequency   carrier frequency, in MHz
        zerotime    position (in points) on the time zero
        
        """
        super(NMRAxis, self).__init__(size = size, itype = itype)
        self.specwidth = specwidth  # spectral width, in Hz
        self.offset = offset        # position in Hz of the rightmost point
        self.frequency = frequency  # carrier frequency, in MHz
        self.zerotime = 0.0         # position (in points) on the time zero
        self.P0 = 0.0      # phase parameters
        self.P1 = 0.0
        self.NMR = "NMR"
        self.kind = "NMR"
        self.units["ppm"] = Unit(name="ppm", converter=self.itop, bconverter=self.ptoi, reverse=True)
        self.units["Hz"]= Unit(name="Hz", converter=self.itoh, bconverter=self.htoi, reverse=True)
        self.units["sec"]= Unit(name="Hz", converter=self.itos, bconverter=self.stoi)  # for FID
        for i in ("specwidth", "offset", "frequency", "NMR"):  # updates storable attributes
            self.attributes.insert(0, i)
        self.currentunit = currentunit
        
    def report(self):
        "high level reporting"
        rst = "NMR axis at %f MHz, "%(self.frequency,)
        if self.itype == 0:
            rst += "%d real points, "%(self.size,)
            if self.size>1:
                rst += "from %f ppm (%f Hz) to %f ppm  (%f Hz)"% \
                    (self.itop(self.size-1), self.itoh(self.size-1), self.itop(0), self.itoh(0))
        else:
            rst += "%d complex pairs,  "%(self.cpxsize,)
            if self.cpxsize>1:
                rst += "from %f ppm (%f Hz) to %f ppm  (%f Hz)"%  \
                    (self.itop(self.size-1), self.itoh(self.size-1), self.itop(0), self.itoh(0))
        return rst
    #-------------------------------------------------------------------------------
    def __extract(self, zoom):
        """
        redefines the axis parameters so that the new axis is extracted for the points [start:end] 
        
        zoom is given in current unit - does not modify the Data, only the axis definition
        """
        start, end = self.getslice(zoom)
        self.specwidth = (self.specwidth * (end - start)) /self.size
        self.offset = self.offset + self.specwidth * (self.size - end)/self.size
        self.size = end-start
        return (start, end)
    #-------------------------------------------------------------------------------
    def extract(self, zoom):
        """
        redefines the axis parameters so that the new axis is extracted for the points [start:end] 
        
        zoom is given in current unit - does not modify the Data, only the axis definition
        """
        start, end = self.getslice(zoom)
        new_specwidth = self.itoh(start+1)-self.itoh(end)  # check itoh() for the +1
        new_offset = self.itoh(end-1)
        self.specwidth = new_specwidth
        self.offset = new_offset
        #self.size = end-start
        return (start, end)
    #-------------------------------------------------------------------------------
    def itos(self,value):
        """
        returns time value (s) from point value
        """
        return 0.5*(value-2*self.zerotime)/self.specwidth 
    def stoi(self,value):
        """
        returns point value (i) from time value (s)
        """
        return 2.0*value*self.specwidth + 2*self.zerotime
    #-------------------------------------------------------------------------------
    def itop(self,value):
        """
        returns ppm value (p) from point value (i)
        """
        ppm_value = self.itoh(value) / self.frequency
        return ppm_value
    #-------------------------------------------------------------------------------
    def htop(self,value):
        """
        returns ppm value (p) from Hz value (h)
        """
        
        ppm_value = value / self.frequency 
        return ppm_value
    #-------------------------------------------------------------------------------
    def htoi(self,value):
        """
        returns point value (i) from Hz value (h)
        """
        pt_value = (self.size -1)*(self.offset - value)/self.specwidth + self.size-1
        return pt_value
    #-------------------------------------------------------------------------------
    def ptoh(self,value):
        """
        returns Hz value (h) from ppm value (p)
        """
        Hz_value = value * self.frequency
        return Hz_value
    #-------------------------------------------------------------------------------
    def ptoi(self,value):
        """
        returns point value (i) from ppm value (p)
        """
        pt_value = self.htoi((value*self.frequency))
        return pt_value
    #-------------------------------------------------------------------------------
    def itoh(self,value):
        """
        returns Hz value (h) from point value (i)
        """
        # N points define N-1 intervals ! hence the -1 in itoh() and htoi()
        hz_value =   (self.size-value-1)*self.specwidth / (self.size-1) + self.offset
        return hz_value
    def freq_axis(self):
        """return axis containing Hz values, can be used for display"""
        return self.itoh(self.points_axis())
    Hz_axis = freq_axis     # two names for this function
    def ppm_axis(self):
        """return axis containing ppm values, can be used for display"""
        return self.itop( self.points_axis() )
########################################################################
class NMRData(_NPKData):
    """
    a working data used by the NPK package
    
    The data is a numpy array, found in self.buffer     can also be accessed directly d[i], d[i,j], ...
    
    1D 2D and 3D are handled, 3 axes are defined : axis1 axis2 axis3
    axes are defined as in NMR
    in 1D, every is in axis1
    in 2D, the fastest varying dimension is in axis2, the slowest in axis1
    in 3D, the fastest varying dimension is in axis3, the slowest in axis1
    see axis_index
    typical properties and methods are :
    utilities:
        .display() 
        .check()
    properties
        .itype
        .dim .size1, .size2, .size3 ...
    moving data :
        .row(i) .col(i) .set_row(i)  .set_col(i)
        .copy()
        .load() .save()
    processing :
        .fft() .rfft() .modulus() .apod_xxx()  sg()  transpose() ...
    arithmetics :
        .fill() .mult .add()
        also direct arithmetics : f = 2*d+e
    all methods return self, so computation can be piped
    etc...
    """
    def __init__(self, dim = 1, shape = None, buffer = None, name = None, debug = 0):
        """
        data initialisation,
        four alternative posibilities :
        - name : file-name is given, the file is read and loaded
        - buffer : numpy buffer is given - used as is, not copied !
        - shape eg : (si1,si2) is given 
        - dim is given (default is dim=1)
        
        the first found takes over the others which are not used
        """
        from .File.GifaFile import GifaFile
        self.debug = debug
        self.frequency = 400.0
        self._absmax = 0.0
        self.noise = 0.0
        self.name = None
        self.level = None

        if buffer is not None:
            dim = len(buffer.shape)
        if name is not None:
            Go = GifaFile(name,"r")
            Go.load()
            Go.close()
            B = Go.get_data()
            del(Go)
            self.buffer = B.buffer  
            copyaxes(B,self)
            self.name = name
        elif shape is not None:
            self.buffer = np.zeros(shape)
            self.axis1 = NMRAxis()
            if (len(shape) == 2):
                self.axis1 = NMRAxis()
                self.axis2 = NMRAxis()
            elif (len(shape) == 3):
                self.axis1 = NMRAxis()
                self.axis2 = NMRAxis()
                self.axis3 = NMRAxis()
        else:
            if dim == 1:
                self.buffer = np.zeros((64,))
                self.axis1 = NMRAxis()
            elif dim == 2:
                self.buffer = np.zeros((64,64))
                self.axis1 = NMRAxis()
                self.axis2 = NMRAxis()
            elif dim == 3:
                self.buffer = np.zeros((64,64,64))
                self.axis1 = NMRAxis()
                self.axis2 = NMRAxis()
                self.axis3 = NMRAxis()
            else:
                raise NPKError("invalid dimension")
        if buffer is not None:
            self.set_buffer(buffer)
        self.adapt_size()
        self.check()
########################################################################
    #-------------------------------------------------------------------------------
    def copy(self):
        """return a copy of itself"""
        c = NPKData._NPKData.copy(self)
        c.frequency = self.frequency
        return c
    def itop(self,axis,value):
        todo = self.test_axis(axis)
        if self.dim == 1:
            return self.axis1.itop(value)
        elif self.dim == 2:
            return self.axes(todo).itop(value)
    def htop(self,axis,value):
        todo = self.test_axis(axis)
        if self.dim == 1:
            return self.axis1.htop(value)
        elif self.dim == 2:
            return self.axes(todo).htop(value)
    def htoi(self,axis,value):
        todo = self.test_axis(axis)
        if self.dim == 1:
            return self.axis1.htoi(value)
        elif self.dim == 2:
            return self.axes(todo).htoi(value)
    def ptoh(self,axis,value):
        todo = self.test_axis(axis)
        if self.dim == 1:
            return self.axis1.ptoh(value)
        elif self.dim == 2:
            return self.axes(todo).ptoh(value)
    def ptoi(self,axis,value):
        todo = self.test_axis(axis)
        if self.dim == 1:
            return self.axis1.ptoi(value)
        elif self.dim == 2:
            return self.axes(todo).ptoi(value)
    def itoh(self,axis,value):
        todo = self.test_axis(axis)
        if self.dim == 1:
            print("going for F1 , on ", value)
            return self.axis1.itoh(value)
        elif self.dim == 2:
            return self.axes(todo).itoh(value)

class NMRDataTests(unittest.TestCase):
    def test_load(self):
        " - Testing load methods"
        from .Tests import filename
        name1D = filename("proj.gs1")
        E = NMRData(name=name1D)
        self.assertAlmostEqual(E[0], 1869.4309082)
        self.assertAlmostEqual(E.get_buffer().max(), 603306.75)
    def test_unitval(self):
        "testing unit conversion functions"
        F = NMRAxis(size=1000, specwidth=12345.0, frequency = 400.0, offset=321)
        self.assertAlmostEqual(F.itoh(0), F.specwidth+F.offset)
        self.assertAlmostEqual(F.itoh(F.size-1), F.offset)   # last point is size-1 !!!
        for f,g in ((F.htoi,F.itoh), (F.ptoi,F.itop), (F.htop,F.ptoh)):
            self.assertAlmostEqual( g(f(4321)), 4321)
        for u in ("points", "Hz", "ppm"):
            F.currentunit = u
            self.assertAlmostEqual( F.ctoi( F.itoc(4321)), 4321)


from . import plugins
plugins.loadfolder(os.path.join(plugins.__path__[0],'NMR'), debug=False)

if __name__ == '__main__':
    # minitest
    d = NMRData(dim=1, debug=1)
    print(d.report())
    d = NMRData(shape=(33,33),debug=1)
    print(d.report())
    d = NMRData(buffer=np.ones((12,24,48)),debug=1)
    print(d.report())
    print('Hello from NMR')