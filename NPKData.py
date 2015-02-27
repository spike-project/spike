#!/usr/bin/env python
# encoding: utf-8
"""
NPKData.py

Implement the basic mechanisms for NMR data-sets

Created by Marc-André and Marie-Aude on 2010-03-17.
Copyright (c) 2010 IGBMC and NMRTEC. All rights reserved.
"""

import version
import numpy as np
import numpy.fft as npfft
import copy
import spike.File.GifaFile
import spike.File.csv
from NPKError import NPKError
import itertools as it
import unittest
import math
import re
import time

__version__ = "0.3.3"   # this is NPKData version
__date__ = "11/Jul/2011"

########################################################################
# series of utilities

def hypercomplex_modulus(arr, size1, size2):
    """
    Calculates the modulus of an array of hypercomplex numbers.
    input:
        arr : hypercomplex array
        size1 : size counting horizontally each half quadrant.
        size2 : siez counting vertically each half quadrant.
    eg:
        arr = np.array([[1, 4],[3, 7],[1, 9],[5, 7]])
        is an hypercomplex with size1 = 2 and size2 = 2
    """
    b = np.zeros((size1/2, size2/2))
    brr = arr[::2, ::2]
    bri = arr[::2, 1::2]
    bir = arr[1::2, ::2]
    bii = arr[1::2, 1::2]
    b = np.sqrt(brr**2 + bri**2 + bir**2 + bii**2 )
    return b
def _t(x):
    return _base_irfft(as_float(x))
def _tt(x):
    return as_cpx(_base_rfft(x))
def as_cpx(arr):
    """
    interpret arr as a complex array
    useful to move between complex and real arrays (see as_float)
    
    >>> print as_cpx(np.arange(4.0))
    [ 0.+1.j  2.+3.j]
    
    """
#    return np.frombuffer(arr,dtype = "complex
    return arr.view(dtype = "complex")
def as_float(arr):
    """
    interpret arr as a float array
    useful to move between complex and real arrays (see as_float)

    >>> print as_float(np.arange(4)*(1+1j))
    [ 0.  0.  1.  1.  2.  2.  3.  3.]
    """
#    return np.frombuffer(arr,dtype = "float")
    return arr.view(dtype = "float")
def conj_ip(a):
    """
    computes conjugate() in-place

    >>> conj_ip(np.arange(4)*(1+1j))
    [ 0.-0.j  1.-1.j  2.-2.j  3.-3.j]
    """
    if a.dtype == np.complex:
        b = as_float(a)[1::2]       # create a view of the imaginary part of a
        np.multiply(b,-1.0,b)       # and inverse it on-place
#        b *= -1.0                  # is equivalent
    return a
def _conj_ip_from_float(a): 
    return np.conjugate(a.view(dtype = "complex"))
def _conj_ip_to_float(a):
    return np.conjugate(a).view(dtype = "float")
def _base_fft(a):
    """
    should not be used for regular use - called by wrapper routines
    
    computes the complex Fourier transform, the NMR way
    built as the conjugate of the fft found in numpy 
    
    WARNING - destroy the buffer given as parameter
    
    test : 
    >> print _base_fft(np.arange(4.0))
    [ 2.  4. -2. -2.]
    
    """
    return _conj_ip_to_float( npfft.fft( _conj_ip_from_float(a) ) )   # reverse,
#    v = npfft.fft( conj_ip( as_cpx(a) ))     # creates a new buffer
#    return as_float( conj_ip(v) )   # reverse,
def _base_ifft(a):
    """
    should not be used for regular use - called by wrapper routines

    computes the inverse complex Fourier transform, the NMR way
    built as the conjugate of the ifft found in numpy 

    WARNING - destroy the buffer given as parameter

    test : 
    >> print _base_ifft(np.arange(4.0))
    array([ 1.,  2., -1., -1.])

    """
    v = npfft.ifft( conj_ip( as_cpx(a) ))     # creates a new buffer
    return as_float( conj_ip(v) )   # reverse,
       
def _base_rfft(a):
    """
    should not be used for regular use - called by wrapper routines
    
    imaginary parts of first and last freq is zero ( [1] and [-1]) so they are dropped and rfft is inplace.
    This works only if self.size2 is even !!!

    test : 
    >> print _base_rfft(np.arange(4.0))
    [ 3.  1. -2. -2.]
    """
    v = as_float(npfft.rfft(as_float(a)))
#            print "dropped:",v[1],v[-1]     
    v[1] = 0.5*v[-2]        # 0.5 ??? c'est peut-etre un bug dans NPK
    v[0] *= 0.5
    return as_float(conj_ip(as_cpx(v[:-2])))
def _base_irfft(a):
    """
    should not be used for regular use - called by wrapper routines
    
    inverse of _base_rfft
    This works only if self.size2 is even !!!

    test : 
    >> print _base_irfft(np.arange(4.0))
    [ 0.5,  2. , -1.5, -1. ]
    """
    v = np.zeros(len(a)+2)
    v[:-2] = as_float((conj_ip(as_cpx(a[:]))))
    v[0] = 2*v[0]
    v[-2] = 2*v[1]
    v[1] = 0.0
    return npfft.irfft(as_cpx(v))
def _base_fftr(a):
    """
    should not be used for regular use - called by wrapper routines

    complex to real direct FT
    This works only if self.size2 is even !!!

    test : 
    >> print _base_fftr(np.arange(4.0))
    [ 2.0, -3.0, -2.0, 3.0 ]
    >> print _base_fftr(np.arange(8.0)+1.0)
    [16.0, -16.31, 0.0, 1.34, -4.0, 6.313, -8.0, 12.65]

    """
    v = np.zeros(len(a)+2)
#    v[:-2] = as_float((conj_ip(as_cpx(a[:]))))
    v[:-2] = a[:]
    v[0] = 2*v[0]
    v[-2] = 2*v[1]
    v[-2] = 0.0
    v =  npfft.irfft(as_cpx(v))
    return v*(len(v)/2)
def _base_ifftr(a):
    """
    should not be used for regular use - called by wrapper routines

    inverse of fftr

    test : 
    >> print _base_ifftr(np.arange(4.0))
    [ 1.5,  0.5, -1.0, -1.0 ]
    >> print _base_ifftr(np.arange(8.0)+1.0)
    [4.5, 0.5, -1.0, 2.41, -1.0, 1.0, -1.0, 0.414]
    """
    v = as_float(npfft.rfft(as_float(a)))
    v[1] = 0.5*v[-2]        # 0.5 ??? c'est peut-être un bug dans NPK
    v[0] *= 0.5
    return v[:-2]*(2.0/(len(v)-2))
def flatten(*arg):
    """
    flatten recursively a list of lists

    >>>print flatten( ( (1,2), 3, (4, (5,), (6,7) ) ) )
    [1, 2, 3, 4, 5, 6, 7]
    """
    import collections
    r = []
    for i in arg:
        if isinstance(i, collections.Sequence):
            #print i
            for j in i:
                r += flatten(j)
        else:
            r += [i]
    return r
def warning(msg):
    """issue a warning message to the user"""
    print "WARNING"
    print msg
########################################################################
class Axis(object):
    """
    hold information for one spectral axis
    used internally
    """
    def __init__(self, size = 64, itype = 0, units = "point"):
        """
        size        number of points along axis
        itype       0 == real, 1 == complex
        units       string which hold the unit name (defaut is "points")
        
        """
        self.size = size            # number of ponts along axis
        self.itype = itype          # 0 == real, 1 == complex
        self.units = units          #
        self.unit_types = ["points"]    # holds possible unit. each unit is associated with a xxx_axis() where xxx is the unit, which returns a axis values
        self.sampling = None        # index list of sampled points if instantiated
        self.sampling_info = {}
        self.attributes = ["itype", "units", "sampling"]    # storable attributes
    def report(self):
        if self.sampling:
            return "size : %d   sampled from %d   itype %d   units %s"%(self.size, max(self.sampling), self.itype, self.units)
        else:
            return "size : %d   itype %d   units %s"%(self.size, self.itype, self.units)
    def typestr(self):
        " returns its type (real or complex) as a string"
        if self.itype == 1:
            tt = "complex"
        else:
            tt =  "real"
        return tt
    def copy(self):
        return copy.copy(self)
    def check_zoom(self, zoom):
        """
        check whether a zoom window, given as (low,high) is valid
        - check low<high and within axis size
        - check that it starts on a real index in itype is complex
        return a boolean
        """
        test =  zoom[0] >= 0 and zoom[0] <= (self.size)
        test = test and zoom[1] >= 0 and zoom[1] <= (self.size)
        test = test and zoom[0] <= zoom[1]
        if self.itype == 1:
            test = test and zoom[0]%2 == 0
        return test
    def load_sampling(self, filename):
        """
        loads the sampling scheme contained in an external file
        file should contain index values, one per line, comment lines start with a #
        complex axes should be sampled by complex pairs, and indices go up to self.size1/2
        
        sampling is loaded into self.sampling  and self.sampling_info is a dictionnary with information
        """
        import Algo.CS_transformations as cstr
        S = cstr.sampling_load(filename)
        self.sampling = S[0]
        #print 'self.sampling ', self.sampling
        self.sampling_info = S[1]
    def get_sampling(self):
        """returns the sampling scheme contained in current axis"""
        return self.sampling
    def set_sampling(self, sampling):
        """sets the sampling scheme contained in current axis"""
        self.sampling = sampling
        return self.sampling
    @property
    def sampled(self):
        """true is sampled axis"""
        return self.sampling is not None
    def points_axis(self):
        """return axis in points units, actually 0..size-1"""
        return np.arange(self.size)
    def unit_axis(self):
        """returns an axis in the unit defined in self.units"""
        u = self.units
        uu = "".join(re.findall("[\w]*",u))
        return getattr(self, uu+"_axis")()
        
class NMRAxis(Axis):
    """
    hold information for one NMR axis
    used internally
    """
    def __init__(self,size = 64, specwidth = 2000.0*math.pi, offset = 0.0, frequency = 400.0, itype = 0, units = "points"):
        """
        all parameters from Axis, plus
        specwidth   spectral width, in Hz
        offset      position in Hz of the rightmost point
        frequency   carrier frequency, in MHz
        zerotime    position (in points) on the time zero
        
        """
        super(NMRAxis, self).__init__(size = size, itype = itype, units = units)
        self.specwidth = specwidth  # spectral width, in Hz
        self.offset = offset        # position in Hz of the rightmost point
        self.frequency = frequency  # carrier frequency, in MHz
        self.zerotime = 0.0         # position (in points) on the time zero
        self.NMR = "NMR"
        self.unit_types.append("ppm")
        self.unit_types.append("Hz")
        for i in ("specwidth", "offset", "frequency", "NMR"):  # updates storable attributes
            self.attributes.insert(0, i)
    def _report(self):
        "low level reporting"
        return "size : %d   freq %f    sw %f   off %f   itype %d units %s"%(self.size, self.frequency, self.specwidth, self.offset, self.itype, self.units)
    def report(self):
        "high level reporting"
        if self.itype == 0:
            return "NMR axis at %f MHz,  %d real points,  from %f ppm (%f Hz) to %f ppm  (%f Hz)"%  \
            (self.frequency, self.size, self.itop(self.size-1), self.itoh(self.size-1), self.itop(0), self.itoh(0))
        else:
            return "NMR axis at %f MHz,  %d complex pairs,  from %f ppm (%f Hz) to %f ppm  (%f Hz)"%  \
            (self.frequency, self.size/2, self.itop(self.size-1), self.itoh(self.size-1), self.itop(0), self.itoh(0))

    #-------------------------------------------------------------------------------
    def extract(self, (start, end)):
        """redefines the axis parameters so that the new axe is extracted for the points [start:end]"""
        if end <0: end = self.size+end+1
        if start<0 or end>self.size or start>=end:
            raise NPKError("The current axis contains %d points, it cannot be extracted from %d to %d"%(self.size, start, end))
        self.specwidth = (self.specwidth * (end - start)) /self.size
        self.offset = self.offset + self.specwidth * (self.size - end)/self.size
        self.size = end-start
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
        pt_value = (self.size -1)*(self.offset - value)/self.specwidth + self.size
        return pt_value
    #-------------------------------------------------------------------------------
    def ptoh(self,value):
        """
        returns Hz value (h) from ppm value (p)
        """
        Hz_value = (value * self.frequency)/1e6 + self.frequency
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
        hz_value =   (self.size-value)*self.specwidth / (self.size-1) + self.offset
        return hz_value
    def freq_axis(self):
        """return axis containing Hz values, can be used for display"""
        return self.itoh(self.points_axis())
    Hz_axis = freq_axis     # two names for this function
    def ppm_axis(self):
        """return axis containing ppm values, can be used for display"""
        return self.itop( self.points_axis() )
########################################################################
class LaplaceAxis(Axis):
    """
    hold information for one Laplace axis (DOSY)
    used internally
    """
    
    def __init__(self, size = 64, dmin = 1.0, dmax = 10.0, dfactor = 1.0, units = "points"):
        super(LaplaceAxis, self).__init__(size = size, itype = 0, units = units)
        self.dmin = dmin
        self.dmax = dmax
        self.dfactor = dfactor
        self.unit_types.append("damping")
        for i in ("dmin", "dmax", "dfactor", "Laplace"):  # updates storable attributes
            self.attributes.append(i)
        
    def itod(self, value):
        """
        returns damping value (d) from point value (i)
        """
        print "itod might have to be checked"
        cst = (math.log(self.dmax)-math.log(self.dmin)) / (float(self.size)-1)
        d = (self.dmin )* math.exp(cst*(value -1.0))
        return d
    def dtoi(self, value):
        """
        returns point value (i) from damping value (d)
        """
        print "dtoi might have to be checked"
        cst = (math.log(self.dmax)-math.log(self.dmin)) / (float(self.size)-1)

        i = 1.0 + (math.log(value) - math.log(self.dmin)) / cst
        return i
    def _report(self):
        "low level report"
        return "Dmin %f   Dmax %f  Dfactor %f"%(self.dmin,self.dmax,self.dfactor)
    def report(self):
        "hight level report"
        return "Laplace axis of %d points,  from %f to %f  using a scaling factor of %f"%  \
            (self.size, self.itod(0.0), self.itod(self.size-1), self.dfactor)
        
########################################################################
def copyaxes(inp,out):
    """
    copy axes values from NPKDAta in to out.
   
    internal use
    """
    for ii in range(inp.dim):
        i = ii +1       # 1 2 3
        setattr(out, "axis%1d"%(i), copy.deepcopy(inp.axes(i)) )
########################################################################
class NPKData(object):
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
        self.debug = debug
        self.frequency = 400.0
        self.absmax = 0.0
        self.noise = 0.0
        self.name = None
        self.level = None

        if name is not None:
            Go = File.GifaFile.GifaFile(name,"r")
            Go.load()
            Go.close()
            B = Go.get_data()
            del(Go)
            self.buffer = B.buffer  
            copyaxes(B,self)
            self.name = name            
        elif buffer is not None:
            self.buffer = buffer
            self.axis1 = NMRAxis()
            t = buffer.shape    # will raise exception if not an array
            try:                # pbs may appear with pytables buffer
                dt = buffer.dtype
            except:
                dt = np.float
            if self.dim>1:
                self.axis2 = NMRAxis()
            if self.dim>2:
                self.axis3 = NMRAxis()
            if dt == np.float:
                self.axes(dim).itype = 0
            elif dt == np.complex:
                self.buffer.dtype = np.float
                self.axes(self.dim).itype = 1
            else:
                raise NPKError("Your data are neither float nor complex", data=self)
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
        self.adapt_size()
        self.check()
    #----------------------------------------------
    def get_buffer(self,copy = False):
        """
        returns a view or a copy of the numpy buffer containing the NPKData values
        dtype is either real or complex if axis is complex.
        remarks :
         - default is a view, if you want a copy, simply do d.get_buffer(copy=True)
         - if you use a view, do not modify the size, nor the dtype
         - see set_buffer()
        WARNING
        - In nD with n>1 and if NPKData is hypercomplex, only the fastest (n) axis is considered, all other imaginary parts are left as real.
        """
        if self.axes(self.dim).itype == 1:
            buf = as_cpx(self.buffer)
        else:
            buf = self.buffer
        if copy:
            return buf.copy()
        else:
            return buf
            
    def set_buffer(self, buff):
        """
        modify the internal buffer of the NPKData.
        allows real or complex arrays to be used
        remarks
         - see get_buffer()
        """
        t = buff.shape    # will raise exception if not an array
        #print "buff.shape ",buff.shape
        #print "self.dim ", self.dim
        if len(t) != self.dim:
            raise NPKError("set_buffer() cannot change the data dimension", data=self)
        try:                # pbs may appear with pytables buffer
            dt = buff.dtype
        except:
            dt = np.float
        if dt == np.complex:
            buff = as_float(buff)
            self.axes(self.dim).itype = 1
        elif dt == 'float':
            self.axes(self.dim).itype = 0
        self.buffer = buff
        self.absmax = 0.0
        self.adapt_size()
        self.check()
        return self
    #----------------------------------------------
    def axes(self,axis):
        """
        returns the required axis : 1, 2 or 3
        """
        return getattr(self,"axis%1d"%(axis))
    @property
    def dim(self):
        "returns the dimension of data : 1 2 or 3 (for 1D 2D or 3D)"
        return len(self.buffer.shape)
    @property
    def size1(self):
        """
        returns the size of the F1 spectral axis in 1D 2D and 3D
        i.e. the unique axis in 1D, the slowest axis in 2D and 3D
        warning, if data along axis is complex, the size is twice the number of complex pairs
        """
        return self.axis1.size
    @property
    def size2(self):
        """
        returns the size of the F2 spectral axis in 2D and 3D
        i.e. the slowest axis in 2D and the intermediate in 3D
        warning, if data along axis is complex, the size is twice the number of complex pairs
        """
        return self.axis2.size
    @property
    def size3(self):
        """
        returns the size of the F3 spectral axis in 3D
        i.e. the slowest axis in 3D
        warning, if data along axis is complex, the size is twice the number of complex pairs
        """
        return self.axis3.size
    @property
    def itype(self):
        "returns complex type of each axes coded as single number, using NPKv1 code"
        if self.dim == 1:
            t =  self.axis1.itype
        elif self.dim == 2:
            t =  self.axis2.itype + 2*self.axis1.itype
        elif self.dim == 3:
            t =  self.axis3.itype + 2*self.axis2.itype + 4*self.axis1.itype
        return t
    def __getitem__(self, key):
        """
        allows d[i] where d is an NPKData
        will always return as if data is real, independently of itype
        """ 
        return self.buffer.__getitem__(key)
    def __setitem__(self, key, value):
        """
        allows d[i] where d is an NPKData
        will always set as if data is real, independently of itype
        """ 
        return self.buffer.__setitem__(key, value)
    #----------------------------------------------
    def _gunits(self):
        "copy units to all the axes"
        return (self.axes(i+1).units for i in range(self.dim) )
    def _sunits(self, units):
        for i in range(self.dim):
            ax = self.axes(i+1)
            ax.units = units
    units = property(_gunits, _sunits)
    #------------------------------------------------
    def check(self, warn = False):
        """
        check basic internal validity
        raises exceptions unless warn is set to True - in which case, only warnings are issued
        can be used in pipes as it returns self if everything is ok
        """
        #----------------------------------------------
        def check_msg(string):
            if warn:
                warning( "WARNING in NPKData.check() : "+string)
            else:
                print self.report()
                raise Exception(string)
        #----------------------------------------------
        try:                            # self.buffer might come from HDF5File and doesn't have flags
            if not self.buffer.flags['OWNDATA']:
                if self.debug >0:
                    print "WARNING in NPKData.check() : NPKData does not own its buffer"
#            print self.buffer.base.flags['UPDATEIFCOPY']
            # I am not sure this is a concern...
        except:
            pass
        
        try:                            # self.buffer might come from HDF5File and doesn't have flags
            if not self.buffer.flags['C_CONTIGUOUS']:
                if self.debug >0:
                    print "WARNING in NPKData.check() : NPKData does not own its buffer"
#            print self.buffer.base.flags['UPDATEIFCOPY']
            # I am not sure this is a concern...
        except:
            pass
        try:
            dt = self.buffer.dtype
        except:
            dt = np.float
        if dt != np.float:
            check_msg("wrong buffer type : %s"%str(dt))
        if len(self.buffer.shape) != self.dim:
            check_msg("wrong dim value : %d while buffer is %d"%(self.dim,len(self.buffer.shape)))
        if self.dim == 1:
            if self.buffer.shape[0] != self.axis1.size:
                check_msg("wrong size value : %d while buffer is %d"%(self.axis1.size, self.buffer.size))
        elif self.dim == 2:
            if self.buffer.shape[0] != self.axis1.size or self.buffer.shape[1] != self.axis2.size:
                check_msg("wrong size value : %d x %d while buffer is %d x %d" % \
                    (self.axis1.size, self.axis2.size, self.buffer.shape[0],self.buffer.shape[1]))
        elif self.dim == 3:
            if self.buffer.shape[0] != self.axis1.size or self.buffer.shape[1] != self.axis2.size  or self.buffer.shape[2] != self.axis3.size:
                check_msg("wrong size value : %d x %d x %d while buffer is %d x %d x %d" % \
                    (self.axis1.size, self.axis2.size, self.axis3.size, self.buffer.shape[0], self.buffer.shape[1], self.buffer.shape[2]))
        for i in range(self.dim):
            if (self.axes(i+1).itype == 1) and (self.axes(i+1).size%2 == 1):
                check_msg("axis %d as size and type mismatch : %d - %d"%(i,self.axes(i+1).itype, self.axes(i+1).size))
        return self
    #----------------------------------------------
    def checknD(self,n):
        if self.dim != n:
            raise NPKError("The dataset is not a %1dD experiment, as required"%n, data=self)
        else:
            return True
    def check1D(self):
        "true for a 1D"
        self.checknD(1)
    def check2D(self):
        "true for a 2D"
        self.checknD(2)
    def check3D(self):
        "true for a 3D"
        self.checknD(3)
    #========================================= 
    # begin of processing functions
    def copy(self):
        """return a copy of itself"""
        Data = type(self)   # NPKData get subclassed, so subclass creator is to be used
        c = Data(buffer = self.buffer.copy())
        copyaxes(self,c)
        return c
    #---------------------------------------------------------------------------
    def adapt_size(self):
        """
        adapt the sizes held in the axis objects to the size of the buffer
        TO BE CALLED each time the buffer size is modified
        otherwise strange things will happen
        """
        sizes = self.buffer.shape
        for i in range (self.dim):
            self.axes(i+1).size = sizes[i]
        self.check()
    #---------------------------------------------------------------------------
    def _chsize1d(self,sz1=-1):
        """
        Change size of data, zero-fill or truncate. 
        DO NOT change the value of OFFSET and SPECW, so EXTRACT should 
        always be preferred on spectra (unless you know exactly what your are doing).
        """
        self.check1D()
        if sz1 == -1:
            sz1 = self.axis1.size
        if sz1<= self.size1:
            self.buffer = self.buffer[:sz1]
        else:
            b = np.zeros(sz1)
            b[:self.size1] = self.buffer
            self.buffer = b
        self.adapt_size()
        return self
    #---------------------------------------------------------------------------
    def _chsize2d(self,sz1=-1,sz2=-1):
        """
        Change size of data, zero-fill or truncate. 
        DO NOT change the value of OFFSET and SPECW, so EXTRACT should 
        always be preferred on spectra (unless you know exactly what your are doing).
        """
        self.check2D()
        if sz1 == -1:
            sz1 = self.axis1.size
        if sz2 == -1:
            sz2 = self.axis2.size
        b = np.zeros((sz1,sz2))
        s1 = min(sz1,self.size1)
        s2 = min(sz2,self.size2)
        b[:s1,:s2] = self.buffer[:s1,:s2]
        self.buffer = b
        self.adapt_size()
        return self
    #---------------------------------------------------------------------------
    def _chsize3d(self,sz1=-1,sz2=-1,sz3=-1):
        """
        Change size of data, zero-fill or truncate. 
        DO NOT change the value of OFFSET and SPECW, so EXTRACT should 
        always be preferred on spectra (unless you know exactly what your are doing).
        """
        self.check3D()
        if sz1 == -1:
            sz1 = self.axis1.size
        if sz2 == -1:
            sz2 = self.axis2.size
        if sz3 == -1:
            sz3 = self.axis3.size
        b = np.zeros((sz1,sz2,sz3))
        s1 = min(sz1,self.size1)
        s2 = min(sz2,self.size2)
        s3 = min(sz3,self.size3)
        b[:s1,:s2,:s3] = self.buffer[:s1,:s2,:s3]
        self.buffer = b
        self.adapt_size()
        return self
    #---------------------------------------------------------------------------
    def chsize(self,sz1=-1,sz2=-1,sz3=-1):
        """
        Change size of data, zero-fill or truncate. 
        DO NOT change the value of OFFSET and SPECW, so EXTRACT should 
        always be preferred on spectra (unless you know exactly what your are doing).
        """   
        if self.dim == 1 :
            self._chsize1d(sz1)
        elif self.dim == 2 :
            self._chsize2d(sz1,sz2)
        else:
            self._chsize3d(sz1,sz2,sz3)
        self.absmax = 0.0
        return self
    #---------------------------------------------------------------------------
    def zf(self, zf1=None, zf2=None, zf3=None, ):
        """
        Zerofill data by adding zeros.
        for a dataset of length size, will add zeros up to zf*size

        do nothing by default unless axis is sampled,
        in which case, missing unsampled points are replaced by 0.0
        """
        if self.dim == 1 :
            if self.axis1.sampled:
                if self.axis1.itype == 0:
                    data = np.zeros(max(self.axis1.sampling) + 1)   # create empty data
                    data[self.axis1.sampling] = self.buffer         # fill them with values
                    self.buffer = data                              # and update
                else:
                    data = np.zeros(max(self.axis1.sampling) + 1)*1j   # create empty data
                    data[self.axis1.sampling] = as_cpx(self.buffer)         # fill them with values
                    self.buffer = as_float(data)                              # and update
                self.axis1.sampling = None
                self.adapt_size()
                    
            if zf1:
                self._chsize1d(self.size1*zf1)
        elif self.dim == 2 :
            if self.axis2.sampled: raise("This is to be done")
            if self.axis1.sampled:
                if self.axis1.itype == 0:
                    data = np.zeros((max(self.axis1.sampling) + 1, self.size2))   # create empty data
                    data[self.axis1.sampling,:] = self.buffer[:,:]              # fill them with values
                    self.buffer = data                                          # and update
                else:
                    data = np.zeros((2*max(self.axis1.sampling) + 2, self.size2))    # create empty data
                    data[2*self.axis1.sampling,:] = self.buffer[::2,:]               # fill  with real values
                    data[2*self.axis1.sampling+1,:] = self.buffer[1::2,:]            # fill  with imag values
                    self.buffer = data                                    # and update
                self.axis1.sampling = None
                self.adapt_size()
            if zf2: self._chsize2d(self.size1, self.size2*zf2)
            if zf1: self._chsize2d(self.size1*zf1, self.size2)
        else:
            if self.axis3.sampled: raise("This is to be done")
            if zf3: self._chsize3d(self.size1, self.size2, self.size3*zf3)
            if zf2: self._chsize3d(self.size1, self.size2*zf2, self.size3)
            if zf1: self._chsize3d(self.size1*zf1, self.size2, self.size3)
        self.absmax = 0.0
        return self
    #---------------------------------------------------------------------------
    def extract(self, *args):
        """
        extract([[x1, y1]])
        extract([x1, y1], [x2, y2]) or extract([x1, y1, x2, y2])
        etc...

        Permits to extract a portion of the data.
        Data can then be processed as a regular data-set.
        EXTRACT changes the value of OFFSET and SPECW accordingly.

            * extract(x1,y1) for 1D datasets.
            * extract(x1, y1, x2, y2) for 2D datasets.

        see also : chsize
        """
        limits = flatten(args)
        if len(limits) != 2*self.dim:
            raise NPKError(msg="wrong arguments for extract :"+str(args), data=self)
        if self.dim == 1:
            self._extract1d(limits)
        elif self.dim == 2:
            self._extract2d(limits)
        elif self.dim == 3:
            self._extract3d(limits)
        self.absmax = 0.0
        return self
    #---------------------------------------------------------------------------
    def _extract1d(self, (x1, y1) ):
        self.check1D()
        self.axis1.extract([x1, y1])
        self.buffer = self.buffer[x1:y1]
        self.adapt_size()
        return self
    #---------------------------------------------------------------------------
    def _extract2d(self, (x1, y1, x2, y2)):
        self.check2D()
        self.axis1.extract([x1, y1])
        self.axis2.extract([x2, y2])
        self.buffer = self.buffer[x1:y1, x2:y2]
        self.adapt_size()
        return self
    #---------------------------------------------------------------------------
    def _extract3d(self, (x1, y1, x2, y2, x3, y3)):
        self.check3D()
        self.axis1.extract([x1, y1])
        self.axis2.extract([x2, y2])
        self.axis3.extract([x3, y3])
        self.buffer = self.buffer[x1:y1,x2:y2,x3:y3]
        self.adapt_size()
        return self
    #---------------------------------------------------------------------------
    def zeroing(self, threshold):
        """
        Sets to zero points below threshold (in absolute value)
        see also :  plus, minus
        """
        self.buffer[abs(self.buffer)<threshold] = 0.0
        return self
    #---------------------------------------------------------------------------
    def plus(self):
        """
        Sets to zero the negative part of the data set
        see also :  minus, zeroing
        """
        self.buffer[self.buffer<0] = 0.0
        self.absmax = 0.0
        return self
    #---------------------------------------------------------------------------
    def minus(self):
        """
        Sets to zero the positive part of the data set
        see also :  minus, zeroing
        """
        self.buffer[self.buffer>0] = 0.0
        return self
        
    #--------------------------------------------------------------------------
    def diag(self,direc = "F12"):
        """
        In 2D, extracts the diagonal of the 2D and put into the 1D buffer.

        In 3D, extracts one diagonal plane of the 3D cube, chosen with the direc parameter 
        and put it into the 2D buffer
        direct values are :

        "F12" is the F1=F2 diagonal
        "F23" is the F2=F3 diagonal
        "F13" is the F1=F3 diagonal

        """
        if self.itype != 0 : 
            print "You must have real data"
        else:
            if self.dim == 1:
                print "Can not extract diagonal from a 1D buffer"
            elif self.dim == 2:
                c =self.diag2D()
            elif self.dim == 3:
                c = self.diag3D(direc)
            return c
    #--------------------------------------------------------------------------
    def diag2D(self):
        self.check2D()
        Data = type(self)   # NPKData get subclassed, so subclass creator is to be used
        if self.size1 > self.size2:
            c = Data(shape =(self.size1,))
            z = float(self.size1)/self.size1
            for i in range(1,self.size2):
                for j in range(0,int(z)):
                    c.buffer [(i-1)*z+j] = self.buffer[i][float(i)/z]
            c.axis1 = copy.copy(self.axis1)
        else:
            z = float(self.size2)/self.size1
            c = Data(shape = (self.size2,))
            for i in range (1,self.size1):
                for j in range(0,int(z)):
                    c.buffer [(i-1)*z+j] = self.buffer[i-1][(i-1)*z+j]
            c.axis1 = copy.copy(self.axis2)
        return c
    #--------------------------------------------------------------------------
    def diag3D(self,direc):
        print "RESTE A FAIRE"
        # self.check3D()
        # print int(direc[1])-1
        # print int(direc[2])-1
#        self._plane2d.buffer = self._image.buffer.diagonal(axis1 = int(direc[1])-1,axis2 = int(direc[2])-1)
    #---------------------------------------------------------------------------
    def zoom(self,dim,*args):
        """
        The basic command for defining region of interest window


        * if n=0, zoom mode is off.
        * if n=1, zoom mode is on,
        """  
        print " ZOOM EST A REVOIR en PROFONDEUR"
        if dim == 1:
            self.zoom1d(*args)
        else:
            if dim == 2 :
                self.zoom2d(*args)
            elif dim == 3 :
                self.zoom3d(*args)
            else:
                print "This should never happen"
        return self
    #---------------------------------------------------------------------------
    def zoom1d(self,*args):     # parameter order : n, f1_low, f1_max
        if args[0] == 0:
            print "Zoom Off"
            self.zo_1_1l = 1
            self.zo_1_1m = self.size1
        elif args[0] == 1:
            if self.itype == 1:
                if args[1]%2 == 0:
                    if self.debug>0: print "We change left value"
                    args[1] = args[1]+1
                if args[2]%2 == 1:
                    if self.debug>0: print "We change right value"
                    args[2] = args[2]+1
            self.zo_1_1l = min(args[1],self.size1)
            self.zo_1_1m = min(args[2],self.size1)
        return self
    #---------------------------------------------------------------------------
    def zoom2d(self,*args): # parameter order : n, f1_low, f2_low, f1_max, f2_max

        if args[0] == 0:
            print "Zoom Off"
            self.zo_2_2l = 1
            self.zo_2_1l = 1
            self.zo_2_2m = self.size2
            self.zo_2_1m = self.size1
        elif args[0] == 1:
            if self.itype == 1:
                if args[2]%2 == 0:
                    args[2] = args[2]+1
                if args[4]%2 == 1:
                    args[2] = args[2]+1
            elif self.itype == 2:
                if args[1]%2 == 0:
                    args[1] = args[1]+1
                if args[3]%2 == 1:
                    args[1] = args[1]+1
            elif self.itype == 3:
                if args[1]%2 == 0:
                    args[1] = args[1]+1
                if args[2]%2 == 0:
                    args[2] = args[2]+1
                if args[3]%2 == 1:
                    args[3] = args[3]+1
                if args[4]%2 == 1:
                    args[4] = args[4]+1
            self.zo_2_2l = min(args[2],self.size2)
            self.zo_2_1l = min(args[1],self.size1)
            self.zo_2_2m = min(args[4],self.size2)
            self.zo_2_1m = min(args[3],self.size1)
        return self
    #---------------------------------------------------------------------------
    def zoom3d(self,*args): # parameter order : n, f1_low, f2_low, f3_low, f1_max, f2_max, f3_max
        if args[0] == 0:
            print "Zoom Off"
            self.zo_3_1l = 1
            self.zo_3_2l = 1
            self.zo_3_3l = 1
            self.zo_3_1m = self.size1
            self.zo_3_2m = self.size2
            self.zo_3_3m = self.size3    # MAD
        elif args[0] == 1:
            if self.itype == 1:
                if args[3]%2 == 0:
                    args[3] = args[3]+1
                if args[6]%2 == 1:
                    args[6] = args[6]+1
            elif self.itype == 2:
                if args[2]%2 == 0:
                    args[2] = args[2]+1
                if args[5]%2 == 1:
                    args[5] = args[5]+1
            elif self.itype == 3:
                if args[2]%2 == 0:
                    args[2] = args[2]+1
                if args[3]%2 == 0:
                    args[3] = args[3]+1
                if args[5]%2 == 1:
                    args[5] = args[5]+1
                if args[6]%2 == 1:
                    args[6] = args[6]+1
            elif self.itype == 4:
                if args[1]%2 == 0:
                    args[1] = args[1]+1
                if args[4]%2 == 1:
                    args[4] = args[4]+1 
            elif self.itype == 5:
                if args[3]%2 == 0:
                    args[3] = args[3]+1
                if args[6]%2 == 1:
                    args[6] = args[6]+1
                if args[1]%2 == 0:
                    args[1] = args[1]+1
                if args[4]%2 == 1:
                    args[4] = args[4]+1
            elif self.itype == 6:
                if args[2]%2 == 0:
                    args[2] = args[2]+1
                if args[5]%2 == 1:
                    args[5] = args[5]+1
                if args[1]%2 == 0:
                    args[1] = args[1]+1
                if args[4]%2 == 1:
                    args[4] = args[4]+1
            elif self.itype == 7:
                if args[3]%2 == 0:
                    args[3] = args[3]+1
                if args[6]%2 == 1:
                    args[6] = args[6]+1
                if args[2]%2 == 0:
                    args[2] = args[2]+1
                if args[5]%2 == 1:
                    args[5] = args[5]+1
                if args[1]%2 == 0:
                    args[1] = args[1]+1
                if args[4]%2 == 1:
                    args[4] = args[4]+1
            self.zo_3_1l = min(args[1],self.size1)
            self.zo_3_2l = min(args[2],self.size2)
            self.zo_3_3l = min(args[3],self.size3)
            self.zo_3_1m = min(args[4],self.size1)
            self.zo_3_2m = min(args[5],self.size2)
            self.zo_3_3m = min(args[6],self.size3)
        
    #----------------------------------------------
    def col(self, i):
        """returns a 1D extracted from the current 2D at position 0<=i<=size2-1 """
        self.check2D()
        Data = type(self)   # NPKData get subclassed, so subclass creator is to be used
        c = Data(buffer = self.buffer[:,i].copy())
        c.axis1 = copy.copy(self.axis1)
        return c
    #----------------------------------------------
    def set_col(self, i, d1D):
        """set into the current 2D the given 1D, as the column at position 0<=i<=size2-1 """
        self.check2D()
        d1D.check1D()
        self.buffer[:,i] = d1D.buffer
        return self
    #----------------------------------------------
    def xcol(self, start=0, stop=None, step=1):
        """
        an iterator over columns of a 2D
        so 
        for c in matrix.xcol():
            do something with c...

        is equivalent to
        for i in range(matrix.size2):     # i.e. all cols
            c = matrix.col(i)
            do something with c...

        you can limit the range by giving start, stop and step arguments - using the same syntax as xrange()
        
        on hypercomplex data
        matrix.xcol( step=matrix.axis2.itype+1 )
        will step only on cols associated to the real point 
        """
        if not stop:
            stop = self.size2
        for index in xrange(start, stop, step):
            yield self.col(index)
        
    @property
    def iter_col(self):
        warning( "iter_col is OBSOLETE, use xcol() instead")
        self.check2D()
        for index in xrange(0, self.size2, self.axis2.itype+1):
            yield self.col(index)
    #----------------------------------------------
    @property
    def iter_all_col(self):
        warning( "iter_col is OBSOLETE, use xcol() instead ")
        self.check2D()
        for index in xrange(self.size2):
            yield self.col(index)
    #----------------------------------------------
    def row(self, i):
        """returns a 1D extracted from the current 2D at position 0<=i<=size1-1 """
        self.check2D()
        Data = type(self)   # NPKData get subclassed, so subclass creator is to be used
        r = Data(buffer = self.buffer[i,:].copy())
        r.axis1 = copy.copy(self.axis2)
        return r
    #----------------------------------------------
    def plane(self, axis, i):
        """returns a 2D extracted from the current 3D at position 0<=i<=size1-1 """
        todo = self.test_axis(axis)
        self.check3D()
        Data = type(self)   # NPKData get subclassed, so subclass creator is to be used        
        if todo == 1 :
            r = Data(buffer = self.buffer[i,:,:].copy())
            r.axis1 = copy.copy(self.axis2)
            r.axis2 = copy.copy(self.axis3)
            return r
        elif todo == 2:
            r = Data(buffer = self.buffer[:,i,:].copy())
            r.axis1 = copy.copy(self.axis1)
            r.axis2 = copy.copy(self.axis3)
            return r
        elif todo == 3:
            r = Data(buffer = self.buffer[:,:,i].copy())
            r.axis1 = copy.copy(self.axis1)
            r.axis2 = copy.copy(self.axis2)
            return r
        else:
            print "problem with plane"
    #----------------------------------------------
    def set_row(self, i, d1D):
        """set into the current 2D the given 1D, as the row at position 0<=i<=size1-1 """
        self.check2D()
        d1D.check1D()
        self.buffer[i,:] = d1D.buffer[:]
        return self
    #----------------------------------------------
    def xrow(self, start=0, stop=None, step=1):
        """
        an iterator over rows of a 2D
        so 
        for r in matrix.xrow():
            do something with r...

        is equivalent to
        for i in range(matrix.size1):     # i.e. all rows 
            r = matrix.row(i)
            do something with r...

        you can limit the range by giving start, stop and step arguments - using the same syntax as xrange()
        
        on hypercomplex data
        matrix.xrow( step=matrix.axis1.itype+1 )
        will step only on rows associated to the real point 
        """
        if not stop:
            stop = self.size1
        for index in xrange(start, stop, step):
            yield self.row(index)
    #----------------------------------------------
    @property
    def iter_row(self):
        warning("iter_row is OBSOLETE, use xrow() instead")
        self.check2D()
        for index in xrange(0, self.size1, self.axis1.itype+1):
            yield self.row(index)
    #----------------------------------------------
    @property
    def iter_all_row(self):
        warning("iter_row is OBSOLETE, use xrow() instead")
        self.check2D()
        for index in xrange(self.size1):
            yield self.row(index)

    #----------------------------------------------
    def check_zoom(self, zoom):
        """
        check whether a zoom window, given as (low,high) or ((low1,high1),(low2,high2))  is valid
        - check low<high and within axis size
        - check that it starts on a real index in itype is complex
        return a boolean
        """
        if not zoom:
            return True
        if self.dim == 1:
            test = self.axis1.check_zoom(zoom)
        elif self.dim == 2:
            test = self.axis1.check_zoom(zoom[0]) and self.axis2.check_zoom(zoom[1])        
        return test
    def display(self, scale = 1.0, absmax = 0.0, show = False, label = None, new_fig = True, axis = None,
                mode3D = False, zoom = None, xlabel="_def_", ylabel = "_def_", figure = None ):

        """
        not so quick and dirty display using matplotlib or mlab - still a first try
        
        scale   allows to increase the vertical scale of display
        absmax  overwrite the value for the largest point, which will not be computed
            display is scaled so that the largest point is first computed (and stored in absmax),
            and then the value at absmax/scale is set full screen 
        show    will call plot.show() at the end, allowing every declared display to be shown on-screen
                useless in ipython
        label   add a label text to plot
        xlabel, ylabel : axes label (default is self.units - use None to remove)
        axis    used as axis if present, axis length should match experiment length
                in 2D, should be a pair (xaxis,yaxis)
        new_fig will create a new window if set to True (default) (active only is figure==None)
        mode3D  use malb 3D display instead of matplotlib contour for 2D display
        zoom    is a tuple defining the zoom window (left,right) or   ((F1_limits),(F2_limits))
        figure  if not None, will be used directly to display instead of using its own
        
        can actually be called without harm, even if no graphic is available, it will just do nothing.
        
        """
        if not self.check_zoom(zoom):
            raise NPKError("wrong zoom window : %s "%(str(zoom)), data=self)
        if not figure:
            import Display.testplot as testplot
            plot = testplot.plot()
            if new_fig:
                plot.figure()
            fig = plot.subplot(111)
        else:
            fig = figure
        if not absmax:  # absmax is the largest point on secptrum, either given from call, or handled internally
            if not self.absmax:     # compute it if absent
                #print "computing absmax...",
                self.absmax = np.nanmax( np.abs(self.buffer) )
        else:
            self.absmax = absmax
        mmin = -self.absmax/scale
        mmax = self.absmax/scale
        if self.dim == 1:
            step = self.axis1.itype+1
            if zoom:
                z1 = zoom[0]
                z2 = zoom[1]
            else:
                z1 = 0
                z2 = self.size1
            if axis is None:
                if self.axis1.units == "Hz":
                    ax = self.axis1.itoh(np.arange(z1,z2,step))
                else:
                    ax = np.arange(z1,z2,step)
            else:
                ax = axis[z1:z2:step]
            fig.plot(ax, self.buffer[z1:z2:step].clip(mmin,mmax), label=label)
            if xlabel == "_def_":
                xlabel = self.axis1.units
            if ylabel == "_def_":
                ylabel = "a.u."
        if self.dim == 2:
            step2 = self.axis2.itype+1
            step1 = self.axis1.itype+1
            if zoom:        # should ((F1_limits),(F2_limits))
                z1lo=zoom[0][0]
                z1up=zoom[0][1]
                z2lo=zoom[1][0]
                z2up=zoom[1][1]
            else:
                z1lo=0
                z1up=self.size1-1
                z2lo=0
                z2up=self.size2-1
            if mode3D:
                print "3D not implemented"
                # from enthought.tvtk.tools import mlab
                # from enthought.pyface.api import GUI
                # gui = GUI()
                # fig = mlab.figure()
                # x = np.arange(self.size1)
                # y = np.arange(self.size2)
                # surf =  mlab.SurfRegularC(x,y,self.f)
                # if self.debug>0: print self.axis1.itype,self.axis2.itype
                # #surf =  mlab.ImShow(self.buffer[::self.axis1.itype+1,::self.axis2.itype+1],scale=[1,1,1])
                # fig.add(surf)
                # gui.start_event_loop()
            else:
                if self.level:
                    level = self.level
                else:
                    m = self.absmax/scale
                    level = (m*0.5, m*0.25, m*0.1, m*0.05)
                    
                    if xlabel == "" and ylabel == "":
                        fig.set_xticklabels('')
                        fig.set_yticklabels('')
                    
                if axis is None:
                    fig.contour(np.arange(z2lo,z2up,step2),np.arange(z1lo,z1up,step1),self.buffer[z1lo:z1up:step1,z2lo:z2up:step2], level)
                else:
#                    for jj in (axis[1][z2lo:z2up:step2],axis[0][z1lo:z1up:step1],self.buffer[z1lo:z1up:step1,z2lo:z2up:step2]):
#                        print jj.shape
                    fig.contour(axis[1][z2lo:z2up:step2],axis[0][z1lo:z1up:step1],self.buffer[z1lo:z1up:step1,z2lo:z2up:step2], level )
            if xlabel == "_def_":
                xlabel = self.axis2.units
            if ylabel == "_def_":
                ylabel = self.axis1.units
        if xlabel is not None:
            fig.set_xlabel(xlabel)
        if ylabel is not None:
            fig.set_ylabel(ylabel)
        if label: fig.legend()
        if show and not figure: plot.show()
        return self
    def f(self,x,y):
        """used by 3D display"""
        return (self.buffer[x,y]/self.absmax)*100
    #----------------------------------------------
    def load(self, name):
        """load data from a file"""
        Go = File.GifaFile.GifaFile(name,"r")
        Go.load()
        Go.close()
        B = Go.get_data()
        del(Go)
        self.buffer = B.buffer  
        copyaxes(B,self)
        self.name = name
        del(B)
        return self
    #----------------------------------------------
    def save(self, name):
        """save data to a file"""
        Go = File.GifaFile.GifaFile(name,"w")
        Go.set_data(self)
        Go.save()
        Go.close()
        del(Go)
        self.name = name
        return self
        
    #----------------------------------------------
    def save_txt(self, name):
        "save 1D data in texte, single column, no unit - with attributes as pseudo comments "
        if self.dim>1:
            raise NPKError("text only possible on 1D", data=self)
        File.csv.save(self, name)
        return self
    #----------------------------------------------
    def load_txt(self, name):
        "load 1D data in texte, single column, no unit - with attributes as pseudo comments "
        buf, att = File.csv.load(name)
        self.buffer = buf
        for k,v in att.items():
            if k in self.axis1.attributes:
                setattr(self.axis1, k, v)
            else:
                print "Warning - wrong attributes : ",k,v
        self.adapt_size()
        return self
    #----------------------------------------------
    def save_csv(self, name):
        """save 1D data in csv,
        in 2 columns : 
        x, y   x values are conditions by the .units attribute
        data attributes are stored as pseudo comments
        
        data can be read back with File.csv.Import_1D()
        """
        if self.dim>1:
            raise NPKError("csv only possible on 1D", data=self)
        File.csv.save_unit(self, name)
        return self
    #----------------------------------------------
    def report(self):
        """reports itself"""
        self.check(warn=True)
        s = "Dim %d"%self.dim
        s =  s+"\nAxis F1 : " + self.axis1.report()
        if self.dim >1:
            s =  s+"\nAxis F2 : " + self.axis2.report()
        if self.dim >2:
            s =  s+"\nAxis F3 : " + self.axis3.report()
        return s
    def __str__(self,*args,**kw):
        """
        express itself as a string,
        add checking here, as this is used during interactive use in ipython
        """
        self.check(warn=True)
        return self.report(*args,**kw)
    def __repr__(self,*args,**kw):
        """
        express itself as a string,
        add checking here, as this is used during interactive use in ipython
        """
        return self.report(*args,**kw)
    #---------------------------------------------------------------------------
    def fill(self, value):
        self.buffer.fill(value)
        self.absmax = value
        return self
    #---------------------------------------------------------------------------
    # Basic arithmetics
    def mult(self,multiplier):
        """
        Multiply data-set by a scalar
        eg : d.mult(alpha) multiplies d buffer by alpha
        """
        import numbers
        # several cases
        # self      real / complex
        # multiplier       real / complex
        
        if not isinstance(multiplier, numbers.Number):
            raise NotImplementedError
        if multiplier.imag == 0:  # if real
            self.buffer *= multiplier
            self.absmax *= multiplier
        else:
            if self.itype == 1:     # means complex
                bb = as_cpx(self.buffer)
                bb *= multiplier
                self.buffer = as_float(bb)
                self.absmax *= multiplier
            else:                           # if complex
                raise NPKError("Multiplication of a real buffer by a complex scalar is not implemented yet",data=self)
        return self
    def __mul__(self,multiplier):
        # as mult but create a new object
        return self.copy().mult(multiplier)
    def __rmul__(self,multiplier):
        # as mult but create a new object
        return self.copy().mult(multiplier)
    def __imul__(self,multiplier):
        # as mult 
        return self.mult(multiplier)
    def __neg__(self):
        return self.copy().mult(-1)
    def __pos__(self):
        return self.copy()
    #-------------------------------------------------------------------------------
    def add(self, otherdata):
        """
        add the provided data : otherdata to the current one
        eg : data.add(otherdata) add content of otherdata to data buffer
        
        can add NPKData and numbers
        """
        import numbers
        if isinstance(otherdata,NPKData):
            if self.itype != otherdata.itype:
                raise NPKError("addition of dataset with different complex states is not implemented yet", data=self)
            self.buffer += otherdata.buffer
            self.absmax = 0.0
        elif isinstance(otherdata,complex):
            if self.itype != 1:
                raise NPKError("cannot add a complex value to this data-set", data=self)
            self.buffer += otherdata
            self.absmax += otherdata
        elif isinstance(otherdata, numbers.Number):
            self.buffer += otherdata
            self.absmax += otherdata
        else:       # then its a number or it fails
            raise NotImplementedError
        return self
    def __iadd__(self, otherdata):
        # as add()
        return self.add(otherdata)
    def __add__(self, otherdata):
        import numbers
        # as add() but creates a new object
        print 'radd'
        if isinstance(otherdata, NPKData):
            return self.copy().add(otherdata)
        elif isinstance(otherdata,numbers.Number):
            return self.copy().add(otherdata)
        else:
            return NotImplemented
    def __radd__(self, otherdata):
        # as add()) but creates a new object
        import numbers
        if isinstance(otherdata, NPKData):
            return self.copy().add(otherdata)
        elif isinstance(otherdata,numbers.Number):
            return self.copy().add(otherdata)
        else:
            return NotImplemented
    def __sub__(self, otherdata):
        # as -add but creates a new object
        import numbers
        if isinstance(otherdata, NPKData):
            return self.copy().add(otherdata.copy().mult(-1))
        elif isinstance(otherdata,numbers.Number):
            return self.copy().add(-otherdata)
        else:
            return NotImplemented
    def __isub__(self, otherdata):
        # as -add 
        import numbers
        if isinstance(otherdata, NPKData):
            return self.add(otherdata.copy().mult(-1))
        elif isinstance(otherdata,numbers.Number):
            return self.add(-otherdata)
        else:
            return NotImplemented
    #-------------------------------------------------------------------------------
    def addbase(self, constant):
        """
        add a constant to the data
        """
        self.buffer += constant
        self.absmax += constant
        return self
    #--------------------------------------------------------------------------
    def addnoise(self, noise, seed=None):
        """
        add to the current data-set (1D, 2D, 3D) a white-gaussian, 
        characterized by its level noise, and the random generator seed.
        """
        print "A VERIFIER!!!"
        if seed is not None:
            np.random.seed(seed)
        self.absmax = 0.0
        self.buffer += noise*np.random.standard_normal(self.buffer.shape)
        return self
    #--------------------------------------------------------------------------
    def addfreq(self, freq, amp=1.0):
        """
        add to the current data-set (1D, 2D, 3D) a single frequency sinusoid
        characterized by its frequency (from axis.specwidth) and amplitude
        """
        print "A VERIFIER!!!"
        if self.dim == 1:
            if self.axis1.itype == 0:
                t = np.arange(self.size1)*freq/self.axis1.specwidth
                self.buffer += amp*np.cos(np.pi*t)
            else:
                t = np.arange(self.size1/2)*freq/self.axis1.specwidth
                self.buffer += as_float( amp*np.exp(1j*np.pi*t) )
        else:
            raise "to be done"
        self.absmax = 0.0
        return self
    #-----------------
    def mean(self, zone):  # ((F1_limits),(F2_limits))
        """
        computes mean value  in the designed spectral zone
        Consider array as real even if itype is 1
        """
        if self.dim == 1:
            ll = zone[0]
            ur = zone[1]
            shift = self.buffer[ll:ur].mean()
        elif self.dim == 2:
            ll = zone[0][0]
            lr = zone[0][1]
            ul = zone[1][0]
            ur = zone[1][1]
            shift = self.buffer[ll:lr,ul:ur].mean()
        return shift
    #-----------------
    def std(self, zone):  # ((F1_limits),(F2_limits))
        """
        computes standard deviation in the designed spectral zone
        Computes value on the real part only ** CHANGED ON July 2012 **
        
        """
        if self.dim == 1:
            ll = zone[0]
            ur = zone[1]
            noise = self.buffer[ll:ur:self.axis1.itype+1].std()
        elif self.dim == 2:
            ll = zone[0][0]
            lr = zone[0][1]
            ul = zone[1][0]
            ur = zone[1][1]
            noise = self.buffer[ll:lr:self.axis1.itype+1,ul:ur:self.axis2.itype+1].std()
        return noise
    #-----------------
    #-----------------
    def test_axis(self, axis=0):
        """
        tests on axis

         in 1D,  axis is not used
                 axis has to be 1 or "F1"
         in 2D,  axis is either 2 for horizontal / faster incremented dimension  == "F2"
                 or 1 for the other dimension == "F1"
                 defaut is 2
         in 3D, axis is 3, 2 or 1 with 3 the faster incremented and 1 the slower == F3 F2 F1
                 defaut is 3
         alternativaly, you may use the strings "F1", "F2" or "F3"
         BUT not F12 F23 as 
         0 is rest to default

        """
        if self.dim == 1:
            if axis in (0,1,"F1","f1"):
                r = 1
            else:
                raise NPKError("Wrong axis parameter : %s in %dD"%(str(axis),self.dim))
        elif self.dim == 2:
            if axis in (0,2,"F2","f2"):
                r = 2
            elif axis in (1,"F1","f1"):
                r = 1
            else:
                raise NPKError("Wrong axis parameter : %s in %dD"%(str(axis),self.dim))
        elif self.dim == 3:
            if axis in  (0,3,"F3","f3"):
                r = 3
            elif axis in (2,"F2","f2"):
                r = 2
            elif axis in (1,"F1","f1"):
                r = 1
            else:
                raise NPKError("Wrong axis parameter : %s in %dD"%(str(axis),self.dim))
        else:
            raise NPKError("Wrong axis parameter : %s in %dD"%(str(axis),self.dim))
        return r
    #---------------------------------------------------------------------------
    def abs(self):
        """
        This command takes the absolute value of the current the data set 
        """
        self.buffer[:] = abs(self.buffer[:])
        return self
    #---------------------------------------------------------------------------
    def real(self, axis=0):
        """
        This command extract the real part of the current the data set 
        considered as complex.
        <ul>
        <li>axis is not needed in 1D, 
        <li>can be F1, F2 or F12 in 2D,
        <li>and can be F1, F2, F3, F12, F13, F23, or F123 in 3D.
        </ul>
        """
        todo = self.test_axis(axis)
        it = self.axes(todo).itype
        if it == 0:
            print "You already have real data - nothing done"
            return self
        if self.dim == 1:
            self.buffer = self.buffer[::2]
        elif self.dim ==2:
            if todo == 1:
                self.buffer = self.buffer[::2]
            elif todo == 2:
                self.buffer = self.buffer[:,::2]
            else:
                print "this should never happen"
        elif self.dim ==3:
            if todo == 1:
                self.buffer = self.buffer[:,::2]
            elif todo == 2:
                self.buffer = self.buffer[::2]
            elif todo == 3:
                self.buffer = self.buffer[:,:,::2]
            else:
                print "this should never happen"
        else:
            print " this should never happen"
        self.axes(todo).itype = 0
        self.adapt_size()
        return self
    #-------------------------------------------------------------------------------
    def swap(self, axis=0):
        """
        swap both parth to complex
        this is the opposite of swa()
        >>>aa=NPKData(buffer=arange(8.))
        >>>aa.axis1.itype=1
        >>>print aa.buffer
        array([ 0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.])
        >>>print aa.swa().buffer
        array([ 0.,  4.,  1.,  5.,  2.,  6.,  3.,  7.])
        """
        todo = self.test_axis(axis)
        it = self.axes(todo).itype
        if it == 1:
            raise NPKError("Dataset should be real along given axis", data=self)
        if self.dim ==1:
            self.buffer = as_float(self.buffer[:self.size1/2] + 1j*self.buffer[self.size1/2:]).copy()
        elif self.dim == 2:
            if todo == 1:
                for i in xrange(self.size2):
                    r = self.col(i).swap()
                    self.set_col(i,r)
            elif todo == 2:
                for i in xrange(self.size1):
                    r = self.row(i).swap()
                    self.set_row(i,r)
        elif self.dim == 3:
            raise NPKError("reste a faire", data=self)
        self.axes(todo).itype = 1
        return self
    #-------------------------------------------------------------------------------
    def unswap(self, axis=0):
        """
        this is the opposite of swap()
        >>>aa=NPKData(buffer=arange(8.))
        >>>aa.axis1.itype=1
        >>>print aa.buffer
        array([ 0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.])
        >>>print aa.unswa().buffer
        array([ 0.,  2.,  4.,  6.,  1.,  3.,  5.,  7.])
        """
        todo = self.test_axis(axis)
        it = self.axes(todo).itype
        if it == 0:
            raise NPKError("Dataset should be complex along given axis", data=self)
        if self.dim ==1:
            cop = as_cpx(self.buffer).copy()
            self.buffer[:self.size1/2] = cop.real
            self.buffer[self.size1/2:] = cop.imag
        elif self.dim == 2:
            if todo == 1:
                for i in xrange(self.size2):
                    r = self.col(i).unswap()
                    self.set_col(i,r)
            elif todo == 2:
                for i in xrange(self.size1):
                    r = self.row(i).unswap()
                    self.set_row(i,r)
        elif self.dim == 3:
            raise NPKError("reste a faire")
        self.axes(todo).itype = 0
        return self
    #-------------------------------------------------------------------------------
    def flip(self):
        """
        on a 2D with axis2.itype==1 and axis1.itype==0
        copies the imaginary from on axis to the other
        after this, we have
            axis2.itype==0 and axis1.itype==1
            size1 is doubled
            size2 is halved
        Useful for complex FT
        this is the opposite of flop()
        
        >>>bb=NPKData(buffer=array([[  0.,   1.,   2.,   3.],[  4.,   5.,   6.,   7.],[  8.,   9.,  10.,  11.],[ 12.,  13.,  14.,  15.]]))
        >>>print bb.buffer
        array([[  0.,   1.,   2.,   3.],
               [  4.,   5.,   6.,   7.],
               [  8.,   9.,  10.,  11.],
               [ 12.,  13.,  14.,  15.]])
        >>>bb.axis2.itype=1
        >>>bb.flip()
        >>>print bb.buffer
        array([[  0.,   2.],
               [  1.,   3.],
               [  4.,   6.],
               [  5.,   7.],
               [  8.,  10.],
               [  9.,  11.],
               [ 12.,  14.],
               [ 13.,  15.]])
        """
#        print "flip A REFAIRE SUR UNSWAP"
        if self.dim ==1:
            raise NPKError("Only in 2D or higher",data=self)
        elif self.dim == 2:
            if not(self.axis2.itype==1 and self.axis1.itype==0):
                raise NPKError("wrong axis itype", data=self)
            self.unswap(axis=2)       # separate real and imag
            self.buffer.resize((2*self.size1,self.size2/2)) # then change size
            self.axis2.itype = 0
            self.axis1.itype = 1
        elif self.dim == 3:
            raise NPKError("reste a faire")
        self.adapt_size()
        return self
    #-------------------------------------------------------------------------------
    def flop(self):
        """
        on a 2D with axis2.itype==0 and axis1.itype==1
        copies the imaginary from on axis to the other
        after this, we have
            axis2.itype==1 and axis1.itype==0
            size1 is halved
            size2 is doubled
        Useful for complex FT
        this is the opposite of flip()
        """
        #print "flop A FAIRE SUR SWAP"
        if self.dim ==1:
            raise NPKError("Only in 2D or higher", data=self)
        elif self.dim == 2:
            
            if not(self.axis2.itype==0 and self.axis1.itype==1):
                raise NPKError("wrong axis itype", data=self)
            #print "resize"
            self.buffer.resize((self.size1/2,self.size2*2)) # then change size
            #print "adapt_size"
            self.adapt_size()
            #print "swap"
            self.swap(axis=2)       # separate real and imag
            #print "initializing the axis"
            self.axis2.itype = 1
            self.axis1.itype = 0
        elif self.dim == 3:
            raise NPKError("reste a faire")
        
        return self
    #-------------------------------------------------------------------------------
    def proj(self, axis=0, projtype="s"):
        """
        returns a projection of the dataset on the given axis
        projtype determines the algorithm :
            "s" is for skyline projection (the highest point is retained)
            "m" is for mean,
        """
        ptype = projtype.lower()
        todo = self.test_axis(axis)
        Data = type(self)   # NPKData get subclassed, so subclass creator is to be used
        if self.dim == 2 :
            
            if todo == 1:
                if ptype == "s":
                    c = Data(buffer =self.buffer.max(1))
                elif ptype == "m":
                    c = Data(buffer =self.buffer.mean(axis = 1))      # mean of each line
            elif todo == 2:
                if ptype == "s":
                    c = Data(buffer =self.buffer.max(0))
                elif ptype == "m":
                    c = Data(buffer =self.buffer.mean(axis = 0))     # mean of each column
        elif self.dim == 3 :
            print "3D"
            # if todo == 1:
            #     if projtype == "s":
            #     elif projtype == "m":
            # elif todo == 2:
            #     if projtype == "s":
            #     elif projtype == "m":
            # elif todo == 3:
            #     if projtype == "s":
            #     elif projtype == "m":
        else :
            print "Dim should be at least 2"
        return c

    #-------------------------------------------------------------------------------
    def phase(self, ph0, ph1, axis=0):
        """
        apply a phase correction along given axis
        phase corrections are in degree
        """
        import math as m
        todo = self.test_axis(axis)
        # compute shape parameters
        it = self.axes(todo).itype
        if it == 0:
            raise NPKError(msg="no phase correction on real data-set", data=self)
        size = self.axes(todo).size/2       # compute correction in e
        if ph1==0:  # we can optimize a little
            e = np.exp(1J*m.radians(float(ph0))) * np.ones( size, dtype=complex)   # e = exp(j ph0)
        else:
            e = np.empty( size,dtype=complex)
            p0 = float(ph0) - 0.5*ph1
            p1 = float(ph1)/(size-1)
            for i in range(size):
                z = m.radians(p0 + i*p1)
                e[i] = m.cos(z) + 1J*m.sin(z)           # e_i = exp(j ph0 + i)
#            e = np.exp( 1J* e)
        # then apply
#        print e
        return self.mult_by_vector(axis,e,mode="complex")
    #-------------------------------------------------------------------------------
    def flipphase(self, ph0, ph1, axis=1):
        """
        equivalent to   flip(); phase();flop()   but much faster
        apply a phase correction along F1 axis of a 2D.
        on 2D where axis1.itype = 0   and   axis2.itype = 1
        using pairs of columns as real and imaginary pair
        phase corrections are in degree
        """
        import math as m
#        print "flipphase A TESTER"
        todo = self.test_axis(axis)
        if todo != 1 or self.dim != 2:
            raise NPKError(msg="works only along F1 axis of 2D", data=self)
        # compute shape parameters
        it = self.axis2.itype
        if it == 0:
            raise NPKError(msg="no phase correction on real data-set", data=self)
        size = self.axis1.size       # compute correction in e
        if ph1==0:  # we can optimize a little
            e = np.exp(1J*m.radians(float(ph0))) * np.ones( size, dtype=complex)   # e = exp(j ph0)
        else:
            e = np.empty( size,dtype=complex)
            p0 = float(ph0) - 0.5*ph1
            p1 = float(ph1)/(size-1)
            for i in range(size):
                z = m.radians(p0 + i*p1)
                e[i] = m.cos(z) + 1J*m.sin(z)           # e_i = exp(j ph0 + i)
        #c = np.empty( size,dtype=complex)
        for i in xrange(self.axis2.size/2):
            #c[:,2*i] = e * (self.buffer[:,2*i] + 1j*self.buffer[:,2*i+1])
            c = e * (self.buffer[:,2*i] + 1j*self.buffer[:,2*i+1])
            self.buffer[:,2*i] = c.real
            self.buffer[:,2*i+1] = c.imag
        return self

    
    #-------------------------------------------------------------------------------
    def bruker_corr(self):
        if self.dim == 1:
            delay = self.axis1.zerotime
        elif self.dim == 2:
            delay = self.axis2.zerotime
        elif self.dim == 3:
            delay = self.axis3.zerotime
        self.phase(0, -360.0*delay, axis=0) # apply phase correction
        return self
    #-------------------------------------------------------------------------------
    def conv_n_p(self):
        """
        realises the n+p to SH conversion
        """
        self.check2D()
        for i in xrange(0,self.size1,2):
            a = self.row(i)
            b = self.row(i+1)
            self.set_row(i,a.copy().add(b) )    # put sum
            a.buffer += -2*b.buffer             # compute diff
            a.buffer = as_float( 1j*as_cpx(a.buffer) )  # dephase by 90°
            self.set_row(i+1,a)                 # and put back
        return self

    #-------------------------------------------------------------------------------
    def urqrd(self, k, orda = None, iterations = 1, axis=0):
        """
        Apply urQRd denoising to data
        k is about 2 x number_of_expected_lines
        Manages real and complex cases.
        Handles the case of hypercomplex for denoising of 2D FTICR for example.
        """
        #from Algo.urQRd import urQRd
        from Algo.urQRd_trick import urQRd
        from util.signal_tools import filtering
        if self.dim == 1:
            if self.axis1.itype == 0:   # real
                buff = as_cpx(_base_ifft(_base_rfft(self.buffer)))       # real case, go to analytical signal
            else:   #complex
                buff = self.get_buffer()                       # complex case, makes complex
            urqrd_result = urQRd( buff, k, orda = orda, iterations = iterations) # performs denoising
            if self.axis1.itype == 0:   # real
                buff = _base_irfft(_base_fft(as_float(urqrd_result)))      # real case, comes back to real
                self.set_buffer(buff)
            else:
                self.buffer = as_float(urqrd_result)             # complex case, makes real
        elif self.dim == 2:
             todo = self.test_axis(axis)
             if todo == 2:
                 for i in xrange(self.size1):
                     r = self.row(i).urqrd(k=k, orda=orda, iterations=iterations)
                     self.set_row(i,r)
             elif todo == 1:
                 for i in xrange(self.size2):
                     r = self.col(i).urqrd(k=k, orda=orda, iterations=iterations)
                     self.set_col(i,r)
        elif self.dim == 3:
             raise Exception("not implemented yet")            
#             urqrd_result = np.zeros(self.buffer.shape)
#             separated = True # Separation of real and imaginary part for processing
#             remove_edge = False  # remove edges in the 2D spectrum
#             passband = [0.05, 0.2] # passband on which signal is conserved vertically in the 2D spectrum
#             if self.axis1.itype == 0:   # real
#                 if separated:
#                     for i in [0,1]:
# #                        self.put_col(self.col(i).urqrd(),i )
#                         if remove_edge: 
#                             self.buffer[:,i] = filtering(self.buffer[:,i], passband)
#                         if self.axis1.itype == 0:   # real
#                             buff = as_cpx(_base_ifft(_base_rfft(self.buffer[:,i])))       # real case, go to analytical signal
#                         else:   #complex
#                             buff = as_cpx(self.buffer[:,i])                       # complex case, makes complex
#                         result = urQRd( buff , k, orda = orda, iterations = iterations) # performs denoising
#                         if self.axis1.itype == 0:   # real
#                             urqrd_result[:,i] = _base_irfft(_base_fft(as_float(result)))      # real case, comes back to real
#                             self.set_buffer(urqrd_result)
#                 else:
#                     result = urQRd( self.buffer[:,0] +1j*self.buffer[:,1], k, orda = orda, iterations = iterations) # performs denoising
#                     urqrd_result[:,0] = np.real(result)
#                     urqrd_result[:,1] = np.imag(result)
#                     self.set_buffer(urqrd_result)
            

        return self
    #-------------------------------------------------------------------------------
    def fft(self, axis=0):
        """
        computes the complex Fourier transform,
        
        takes complex time domain data and returns complex frequency domain data

        see test_axis for information on axis
        """
        self._fft_nD(_base_fft, axis,1,1)
        return self
    #-------------------------------------------------------------------------------
    def fftr(self, axis=0):
        """
        computes the alternate Fourier transform,

        takes complex time domain data and returns real frequency domain data

        see test_axis for information on axis
        """
        self._fft_nD(_base_fftr,axis,1,0)
        return self
    #-------------------------------------------------------------------------------
    def rfft(self, axis=0):
        """
        computes the real Fourier transform,
        takes real time domain data and returns complex frequency domain data

        see test_axis for information on axis
        """
        self._fft_nD(_base_rfft,axis,0,1)
        return self
    #-------------------------------------------------------------------------------
    def ifft(self, axis=0):
        """
        computes the inverse of fft(),
        takes complex frequency domain data and returns complex time domain data

        see test_axis for information on axis
        """
        self._fft_nD(_base_ifft,axis,1,1)
        return self
    #-------------------------------------------------------------------------------
    def irfft(self, axis=0):
        """
        computes the inverse of rfft(),
        takes complex frequency domain data and returns real time domain data

        see test_axis for information on axis
        """
        self._fft_nD(_base_irfft,axis,1,0)
        return self
    #-------------------------------------------------------------------------------
    def ifftr(self, axis=0):
        """
        computes the inverse of fftr,
        takes real frequency domain data and returns complex time domain data

        see test_axis for information on axis
        """
        self._fft_nD(_base_ifftr,axis,0,1)
        return self
    #-------------------------------------------------------------------------------
    def _fft_nD(self, fft_base, axis, it_before, it_after):
        """
        should not be used for regular use - called by wrapper routines
        
        computes the nD Fourier transform, 
        uses fft_base as the method to do 1D fft on a given buffer
        it_before, it_after determines the state before and after transform, and should be either 0 (real) or 1 (cpx)
        FT is performed inplace
        """
        from time import time
        t0 = time()
        todo = self.test_axis(axis)
        if self.axes(todo).itype != it_before:      # check complex along axis to be transformed
            raise NPKError("wrong itype", data=self)     
        if self.dim == 1:                           # then select depending on dim and todo
            self.buffer = fft_base(self.buffer)
        elif self.dim == 2:
            if todo == 2:
                for i in xrange(self.size1):
#                    if i%16==0: print i,self.size1
                    self.buffer[i,:] = fft_base(self.buffer[i,:])
            elif todo == 1:
                a = np.zeros(self.size1)
                for i in xrange(self.size2):
#                    if i%16==0: print i,self.sizeZ
                    #a = np.zeros(self.size1)
                    a[:] = self.buffer[:,i]
                    self.buffer[:,i] = fft_base(a)
        elif self.dim == 3:
            print "A TESTER"
            if todo == 3:
                for i in xrange(self.size1):
                    for j in xrange(self.size2):
                        self.buffer[i,j,:] = fft_base(self.buffer[i,j,:])
            elif todo == 2:
                for i in xrange(self.size1):
                    for j in xrange(self.size3):
                        self.buffer[i,:,j] = fft_base(self.buffer[i,:,j])
            elif todo == 1:
                for i in xrange(self.size2):
                    for j in xrange(self.size3):
                        self.buffer[:,i,j] = fft_base(self.buffer[:,i,j])
        self.axes(todo).itype = it_after
        self.absmax = 0.0
        #print "********** FFT :",todo,time()-t0
    
    #---------------------------------------------------------------------------
    def apply_process(self, axis_it, process, axis = 0, mp = True, N_proc = None):
        """
        scans through given data, using axis_it which is an iterator,
        applying process method (by its name)
        store results into self, along to axis
        if axis_it iterates over self, then processing is in-place
            however it can iterate over an other data-set, thus importing the data

        if self.dim is equal to axis_it().dim, then data are 

        if mp,  does it in a multiprocessing fashion using multiprocessing.Pool()
        if N_proc is None, finds the optimum number itself.
        """
        def _do_it(arg):
            data, process = arg
            p = getattr(data, process)
            return p()
        todo = self.test_axis(axis)
        results = it.imap(_do_it, it.izip(axis_it, it.repeat(process)))  # generate results
        r0 = results.next() # store first one
        if r0.dim <= self.dim:      # dim should reduce !
            raise NPKError("wrong dimension", data=self)
        if self.dim == 1:
            self.buffer[:] = r0.buffer[:]   # initialize
            for r in results:             # then accumulate
                self.buffer[:] += r.buffer[:]
        elif self.dim == 2:
            if todo == 2:
                self.buffer[0,:] = r0.buffer[:]
                i = 0
                for r in results:
                    i += 1
                    self.buffer[i,:] = r.buffer[:]
            elif todo == 1:
                self.buffer[:,0] = r0.buffer[:]
                i = 0
                for r in results:
                    i += 1
                    self.buffer[:,i] = r.buffer[:]
        elif self.dim == 3:
            print "A FAIRE"
            
    #---------------------------------------------------------------------------
    def transpose(self, axis=0):
        """
        Transposes the 2D matrix or planes of the 3D cube. The sizes of 
        the matrix must be a power of two for this command to be used. After 
        transposition, the two dimensions are completely permuted
        
        axis is used in 3D to tell which submatrices should be transposed
        
        
        see also : sym chsize modifysize
        """
        todo = self.test_axis(axis)
        if self.dim == 2:
            if np.isInt(self.size1/2) :
                print "Ok"
            
    #---------------------------------------------------------------------------
    def reverse(self, axis=0):
        """
        reverse the order of the current data-set (i.e. first points are 
        last, last points are first).
        If dataset is complex, REVERSE will reverse the complex  vector (2  by 2).
        """
        
        todo = self.test_axis(axis)
        it = self.axes(todo).itype
        if self.dim == 1 :
            if it == 0:
                self.buffer = self.buffer[::-1].copy() # we want a genuine buffer, not a view
            elif it == 1:
                v = self.buffer[::-1].copy()
                self.buffer[::2]=v[1::2]    # reverse separatly real and imag
                self.buffer[1::2]=v[::2]
            else:
                raise NPKError("internal error")
        elif self.dim == 2:  # en 2D F1 est suivant y , F2 suivant x
            if todo == 1:
                if it == 0:
                    self.buffer = self.buffer[::-1].copy()
                elif it == 1:
                    v = self.buffer[::-1].copy()
                    self.buffer[::2]=v[1::2]    # reverse separatly real and imag
                    self.buffer[1::2]=v[::2]
            elif todo == 2:
                if it == 0:
                    self.buffer = self.buffer[:,::-1].copy()
                elif it == 1:
                    v = self.buffer[:,::-1].copy()
                    self.buffer[:,::2]=v[:,1::2]    # reverse separatly real and imag
                    self.buffer[:,1::2]=v[:,::2]
                else:
                    raise NPKError("internal error")
        elif self.dim == 3:  # en 3D F1 est suivant y , F2 est suivant z et F3 est suivant x
            if todo == 1:
                if it == 0:
                    self.buffer = self.buffer[:,::-1].copy()
                elif it == 1:
                    v = self.buffer[:,::-1].copy()
                    self.buffer[:,::2]=v[:,1::2]    # reverse separatly real and imag
                    self.buffer[:,1::2]=v[:,::2]
            elif todo == 2:
                if it == 0:
                    self.buffer = self.buffer[::-1].copy()
                elif it == 1:
                    v = self.buffer[::-1].copy()
                    self.buffer[::2]=v[1::2]    # reverse separatly real and imag
                    self.buffer[1::2]=v[::2]
            elif todo == 3:
                if it == 0:
                    self.buffer = self.buffer[:,:,::-1].copy()
                elif it == 1:
                    v = self.buffer[:,:,::-1].copy()
                    self.buffer[:,:,::2]=v[:,:,1::2]    # reverse separatly real and imag
                    self.buffer[:,:,1::2]=v[:,:,::2]
                else:
                    raise NPKError("internal error")
        else:
            print "This should never happen ",self.dim
        return self
    #-------------------------------------------------------------------------------
    def conjg(self, axis=0):
        """
        take the inverse conjugate of the buffer
        """
        todo = self.test_axis(axis)
        if self.axes(todo).itype != 1:      # check complex along axis to be transformed
            raise NPKError("wrong itype")
        if self.dim == 1:
            self.buffer[1::2] *= -1
        elif self.dim == 2:
            if todo == 2:
                for i in xrange(self.size1):
                    self.buffer[i,1::2] *= -1
            elif todo == 1:
                for i in xrange(1,self.size1,2):    # by line is faster
                    self.buffer[i,:] *= -1
        else:
            raise NPKError("a faire")
        return self
    #-------------------------------------------------------------------------------
    def revf(self, axis=0):
        """
         Processes FID data-sets by multiplying by -1 2 points out of 4. 
         Permits to preprocess Bruker FIDs in Dim 2 (Bruker trick) before 
         RFT, or permits to bring back zero frequency in the center for some 
         other data formats
        
        """
        todo = self.test_axis(axis)
        if self.dim == 1:
            self.buffer[2::4] *= -1
            self.buffer[3::4] *= -1
        elif self.dim == 2:
            if todo == 2:
                for i in xrange(self.size1):        # columns
                    self.buffer[i,2::4] *= -1
                    self.buffer[i,3::4] *= -1
            elif todo == 1:
                for i in xrange(2,self.size1,4):    # by line is faster
                    self.buffer[i,:] *= -1
                    self.buffer[i+1,:] *= -1
        else:
            raise NPKError("a faire")
        return self
    #-------------------------------------------------------------------------------
    def apod_gm(self, axis=0,gb=1.0):
        """
        apply an gaussian apodisation, gb is in Hz
        WARNING : different from common definition of apodisation
        """
        todo = self.test_axis(axis)
        it = self.axes(todo).itype
        sw = self.axes(todo).specwidth
        size = self.axes(todo).size
        if it == 1: # means complex
            size = size/2
        e = np.exp(-gb*(np.arange(size)/sw)**2)
        if it == 1:
            e = as_float((1 + 1.0j)*e)
        return self.apod_apply(axis,e)
    #-------------------------------------------------------------------------------
    def _spline_interpolate(self, xpoints, kind = 3):
        """compute and returns a spline function 
            we are using splrep and splev instead of interp1d because interp1d needs to have 0 and last point
            it doesn't extend.
        """
        from scipy import interpolate
        self.check1D()
        y=[]
        for point in xpoints:
            y.append(self.buffer[point])
        if len(xpoints) > 2:
            tck = interpolate.splrep(xpoints, y,k=kind)
        elif len(xpoints) == 2 :
            tck = interpolate.splrep(xpoints, y,k=1)
        else:                        # if only one points given, returns a constant, which is the value at that point.
            pass
        if len(xpoints) > 1:
            return interpolate.splev(np.arange(len(self.buffer)),tck,der=0)
        else:                        # if only one points given, returns a constant, which is the value at that point.
            return self.buffer(xpoints)
            
        #return interpolate.interp1d(xpoints, self.buffer[xpoints], kind=kind)
    #-------------------------------------------------------------------------------
    def _linear_interpolate(self,xpoints):
        """computes and returns a linear interpolation"""
        from scipy.optimize import leastsq
        def affine(x,P):
            """docstring for affine"""         
            return P[0]*x+P[1]
        def residual(P,x,y):
            """docstring for residual"""
            return y-affine(x,P)
        y=[]
        x= np.array(xpoints)       # If we keep xpoints as it is, we get:leastsq operands could not be broadcast together with shapes (2) (0)
        for point in xpoints:
            y.append(self.buffer[point])
        res = leastsq(residual,[0.0,0.0], args=(x,y))
        return lambda x: affine(x,res[0])
    #-------------------------------------------------------------------------------
    def linear_interpolate(self, xpoints, axis = 'F2'):
        """"compute and applies a linear function as a baseline correction"""
        
        if self.dim == 1:
            f = self._linear_interpolate(xpoints)
            for i in xrange(len(self.buffer)):
                self.buffer[i] -= f(i)
        elif self.dim == 2:
            if self.test_axis(axis) == 2:
                for i in xrange(self.size1):
                    r = self.row(i)
                    f = r._spline_interpolate(xpoints)
                    print f,self.buffer
                    self.buffer[i,:] -= f
                    print self.buffer
        else:
            raise NPKError("not implemented")
        return self
    #-------------------------------------------------------------------------------
    def spline_interpolate(self, xpoints, axis = 'F2', kind = 3):
        """compute and applies a spline function as a baseline correction"""
        if self.dim == 1:
            f = self._spline_interpolate(xpoints, kind = kind)
            self.buffer -= f
        elif self.dim == 2:
            if self.test_axis(axis) == 2:
                for i in xrange(self.size1):
                    r = self.row(i)
                    f = r._spline_interpolate(xpoints, kind=kind)
                    self.buffer[i,:] -= f
        else:
            raise NPKError("not implemented")
        return self
    #-------------------------------------------------------------------------------
    def apod_tm(self, axis=0, tm1=0, tm2=0):
        """
        apply a trapezoide apodisation, lb is in Hz
        WARNING : different from common definition of apodisation
        This commands applies a trapezoid filter function to the data-
        set. The function raises from 0.0 to 1.0 from the first point to 
        point n1. The function then stays to 1.0 until point n2, from which 
        it goes down to 0.0 at the last point.
        If in 2D or 3D then Fx tells on which axis to apply the filter.
        """
        todo = self.test_axis(axis)
        it = self.axes(todo).itype
        size = self.axes(todo).size
        if it == 1: # means complex
            #size = size/2
            tm1 = min(size,2*(tm1/2)+1)
            tm2 = 2*(tm2/2)+1
        ftm1 = tm1
        ftm2 = size-tm2+1
        e = np.zeros(size)
        if it==0:
            for i in range(1,ftm1):
                e[i] = float(i)/ftm1+1
            for i in range(ftm1,tm2):
                e[i] = 1.
            for i in range(tm2, size):
                e[i] =float(i)/size
            print "****** ",e[0]
        elif it ==1:
            for i in range(1,ftm1,2):
                e[i] = float(i)/ftm1
                e[i+1] = float(i)/ftm1
            for i in range(ftm1,tm2,2):
                e[i] = 1.
                e[i+1] = 1.
            for i in range(tm2, size-1,2):
                e[i] = float(size-i+1)/ftm2
                e[i+1] = float(size-i+1)/ftm2
        #if it == 1:
        #    e = as_float((1 + 1.0j)*e)
        print "APOD_TM still to be doublechecked",e
        #return self.apod_apply(axis,e)
    #-------------------------------------------------------------------------------
    def apod_em(self, axis=0,lb=1.0):
        """
        apply an exponential apodisation, lb is in Hz
        WARNING : different from common definition of apodisation
        """
        todo = self.test_axis(axis)
        it = self.axes(todo).itype
        sw = self.axes(todo).specwidth
        size = self.axes(todo).size
        if it == 1: # means complex
            size = size/2
        e = np.exp(-lb*np.arange(size)/sw)
        if it == 1:
            e = as_float((1 + 1.0j)*e)
        return self.apod_apply(axis,e)
    #-------------------------------------------------------------------------------
    def apod_sq_sin(self, axis=0, maxi=0):
        """
        apply a squared sinebell apodisation
        maxi ranges from 0 to 0.5
        """
        import math as m
        if maxi<0.0 or maxi>0.5:
            raise ValueError
        todo = self.test_axis(axis)
        # compute shape parameters
        it = self.axes(todo).itype
        size = self.axes(todo).size
        if it == 1:
            size = size/2
        s = 2*(1-maxi)
        zz = m.pi/((size-1)*s)
        yy = m.pi*(s-1)/s         #yy has the dimension of a phase
        # then draw buffer
        e = np.sin( zz*np.arange(size)+yy)**2
        if it == 1:
            e = as_float((1 + 1.0j)*e)
        # then apply
        return self.apod_apply(axis,e)
        
    #-------------------------------------------------------------------------------
    def apod_sin(self, axis=0, maxi=0):
        """
        apply a sinebell apodisation
        maxi ranges from 0 to 0.5
        """
        import math as m
        if maxi<0.0 or maxi>0.5:
            raise ValueError
        todo = self.test_axis(axis)
        # compute shape parameters
        it = self.axes(todo).itype
        size = self.axes(todo).size
        if it == 1:
            size = size/2
        s = 2*(1-maxi)
        zz = m.pi/((size-1)*s)
        yy = m.pi*(s-1)/s         #yy has the dimension of a phase
        # then draw buffer
        e = np.sin( zz*np.arange(size)+yy)
        if it == 1:
            e = as_float((1 + 1.0j)*e)
        # then apply
        return self.apod_apply(axis,e)
        
    #-------------------------------------------------------------------------------
    def apod_apply(self, axis, apod_buf):
        """
        apply an apodisation, held into the buffer apod_buf
        """
        todo = self.test_axis(axis)
        if self.dim == 1:
            self.buffer = self.buffer*apod_buf
        if self.dim == 2:
            if todo == 2:
                self.buffer = self.buffer*apod_buf         # broadcast does its work
            if todo == 1:
                for i in xrange(self.size2):
                    self.buffer[:,i] = self.buffer[:,i]*apod_buf
        if self.dim == 3:
            if todo == 3:
                self.buffer = self.buffer*apod_buf         # broadcast does its work
            if todo == 2:
                for i in xrange(self.size1):
                    for j in xrange(self.size3):
                        self.buffer[i,:,j] = self.buffer[i,:,j]*apod_buf
            if todo == 1:
                for i in xrange(self.size2):
                    for j in xrange(self.size3):
                        self.buffer[:,i,j] = self.buffer[:,i,j]*apod_buf
        return self

    #-------------------------------------------------------------------------------
    def mult_by_vector(self, axis, vector, mode="real"):
        """
        multiply the data-set by a vector, along a given axis
        if mode == "real", does it point by point regardles of itype
        if mode == "complex" uses axis.itype to determine how handle complex values
            in all cases vector can be real or complex
        """
        todo = self.test_axis(axis)
        if mode == "complex":
            if self.axes(todo).itype == 1:
                tf = as_cpx
            else:
                tf = as_float
        elif mode == "real":
            raise NPKError('reste a faire')
        else:
            raise NPKError("error with mode", data=self)
        if self.dim == 1:
            self.buffer = as_float(tf(self.buffer)*vector)
        if self.dim == 2:
            if todo == 2:
                self.buffer = as_float(tf(self.buffer)*vector)         # broadcast does its work
            if todo == 1:
                for i in xrange(self.size2):
                    self.buffer[:,i] = as_float(tf(self.buffer[:,i].copy())*vector)
        if self.dim == 3:
            print "A VERIFIER"
            if todo == 3:
                self.buffer = as_float(tf(self.buffer)*vector)         # broadcast does its work
            if todo == 2:
                for i in xrange(self.size1):
                    for j in xrange(self.size3):
                        self.buffer[i,:,j] = as_float(tf(self.buffer[i,:,j].copy())*vector)
            if todo == 1:
                for i in xrange(self.size2):
                    for j in xrange(self.size3):
                        self.buffer[:,i,j] = as_float(tf(self.buffer[:,i,j].copy())*vector)
        self.absmax = 0.0
        return self
        
    #-------------------------------------------------------------------------------
    def median(self):
        """
        Executes a median filter on the data-set (1D or 2D).a window of x 
        points (or y by x in 2D) is moved along the data set, the point are 
        ordered, and the indexth point is taken as the new point for the 
        data set.
        """
      
    #-------------------------------------------------------------------------------
    def modulus(self):
        """
        takes the modulus of the dataset
        depends on the value f axis(i).itype
        """
        if self.dim == 1:
            if self.axis1.itype != 1:      # has to be complex
                raise NPKError("wrong itype", data=self)
            d = as_cpx(self.buffer)
            self.buffer = np.real(np.sqrt(d*d.conj()))
            self.axis1.itype = 0
        elif self.dim == 2:
            if self.axis1.itype == 0 and self.axis2.itype == 0:      # has to be complex
                print ("real data, nothing to do")
            else:   # do something
                if self.axis1.itype == 1 and self.axis2.itype == 0: # along F1 only
                    b = np.zeros((self.size1/2, self.size2))    # create new smaller array
                    for i in xrange(0,self.size1,2):    # new version is about x3 faster
                        dr = self.buffer[i,:]
                        di = self.buffer[i+1,:]
                        b[i/2,:] = np.sqrt(dr**2 + di**2)        # and copy modulus in b 
                    self.axis1.itype = 0                        # finally update axis
                elif self.axis1.itype == 0 and self.axis2.itype == 1: # along F2 only
                    b = np.zeros((self.size1,self.size2/2))
                    for i in xrange(self.size1):
                        d = as_cpx(self.buffer[i,:])    # copy() not needed in F2
                        b[i,:] = np.sqrt(np.real(d*d.conj()))
                    self.axis2.itype = 0
                elif self.axis1.itype == 1 and self.axis2.itype == 1: # along both dim
                    self.axis1.itype = 0
                    self.axis2.itype = 0
                    b =  hypercomplex_modulus(self.get_buffer(), self.size1, self.size2)
                    
                self.buffer = b
        elif self.dim == 3:
            raise NPKError("reste a faire")
        self.absmax = 0.0
        self.adapt_size()
        return self
    #----------------------------------------------------------
    def peaks2d(self, threshold = 0.1, zoom = None, value = False):
        '''
        Extract peaks from 2d Array dataset
        if value is True, return the magnitude at position (x,y)
        '''
        self.check2D()
        #print threshold
        if zoom:        # should ((F1_limits),(F2_limits))
            z1lo=zoom[0][0]
            z1up=zoom[0][1]
            z2lo=zoom[1][0]
            z2up=zoom[1][1]
        else:
            z1lo=0
            z1up=self.size1-1
            z2lo=0
            z2up=self.size2-1
        buff = self.buffer[z1lo:z1up, z2lo:z2up]            # take the zoom window
        listpk=np.where(((buff > threshold*np.ones(buff.shape))&            # thresholding
                        (buff > np.roll(buff,  1, 0)) &         # laplacian - kind of
                        (buff > np.roll(buff, -1, 0)) &
                        (buff > np.roll(buff,  1, 1)) &
                        (buff > np.roll(buff, -1, 1))))
        listpk = [int(z1lo) + listpk[0], int(z2lo) + listpk[1]]         # absolute coordinates
        if value: 
            return listpk[0], listpk[1], self.buffer[listpk[0], listpk[1]]
        else:
            return listpk[0], listpk[1] # list f1, f2
    #-------------------------------------------------------
    # def peak(self, threshold=0.1, zoom=None):
    #     '''
    #     Extract peaks from 1d Array dataset
    #     '''
    #     print threshold
    #     if zoom:        # should ((F1_limits),(F2_limits))
    #         left = zoom[0]
    #         right = zoom[1]
    #        
    #     else:
    #         left = 0
    #         right = self.size1-1
    # 
    # 
    #     buff = self.buffer[left:right]# take the zoom window
    #     print "flags", buff.flags
    #     print "OWNDATA", buff.flags['OWNDATA']
    #     print "type(buff)",type(buff)
    # 
    # 
    #     listpk=np.where(((buff > threshold*np.ones(buff.shape))&# thresholding
    #                     (buff > np.roll(buff,  1, 0)) &     # laplacian - kind of
    #                     (buff > np.roll(buff, -1, 0))))
    #     print "from peak",len(listpk[0])
    #     listpk=int(left)+listpk[0]# absolute coordinates
    #     listamp = buff[listpk]
    #     return listpk,listamp
         
    def peak(self, pos_neg = 1, threshold = 0.1, offset = None):
        """
        first trial for peakpicker
        1D only
        pos_neg = 1 / -1 / 0   :  type of peaks positive / negative / 
        threshold = minimum level, as absolute value
        self.peaks : index of the peaks
        self.peaks_ordered : index of the ordered peaks from maximum to minimum.
        """
        self.check1D()
        a_buf = self.get_buffer()
        if self.itype == 1:
            a_buf = abs(a_buf)
        a_rest = a_buf[1:-1]
        #valmin = threshold*np.max(np.abs(a_rest))            # minimumvalue to be a peak
        valmin = threshold
        diff = (a_rest-a_buf[:-2])*(a_rest-a_buf[2:])       # positive at extremum
        peaks = diff > 0                              # true at extremum
        if pos_neg == 1:
            peaks = np.logical_and(peaks , (a_rest-a_buf[:-2])>0)       # true at positive extremum
            peaks = np.logical_and(peaks, (a_rest>=valmin))
        elif pos_neg == -1:
            peaks = np.logical_and(peaks , (a_rest-a_buf[:-2])<0)       # true at negative extremum
            peaks = np.logical_and(peaks, (a_rest<=valmin))
        elif pos_neg == 0:
            peaks = np.logical_and(peaks, np.abs(a_rest)>=valmin)
        (ipeaks,) = np.where(peaks)     # returns index list
        self.peaks = ipeaks+1   # +1 because to compensate initial shifts
        self.peaks_ordered = self.peaks[np.argsort(a_buf[self.peaks])[::-1]]        # list from maximum to minimum
        return self.peaks
    #-------------------------------------------------------
    def centroid1d(self, npoints=3):
        """
        from peak lists determined with peak()
        realize a centroid fit of the peak summit and width,
        computes Full width at half maximum
        creates lists self.centered_peaks and self.width_peaks
        
        Temporary  so far,
            only based on regular sampling, not unit axis.
            ah-hoc structure, waiting for a real PEAK object
        """
        from scipy import polyfit
        noff = (int(npoints)-1)/2
        self.centered_peaks = []
        self.width_peaks = []
        if (2*noff+1 != npoints) or (npoints<3):
            raise NPKError("npoints must odd and >2 ",data=self)
        for i, pk in enumerate(self.peaks):
            xdata = np.arange(pk-noff, pk+noff+1)
            ydata = self.get_buffer()[xdata]
            coeffs = polyfit(xdata, ydata, 2)
            pospeak = -coeffs[1]/(2*coeffs[0])#/factexp# position of maximum of parabola in m/z
            c = coeffs[2]/2 + coeffs[1]**2/(8*coeffs[0]) #half height 
            b = coeffs[1]
            a = coeffs[0]
            delt = b**2-4*a*c
            width =  np.sqrt(delt)/abs(a) # width at half height
            self.centered_peaks.append(pospeak)
            self.width_peaks.append(width)
        self.centered_peaks = np.array(self.centered_peaks)
        self.width_peaks = np.array(self.width_peaks)
    #-------------------------------------------------------
    def display_peaks(self, axis = None, peak_label = False, zoom = None, show = False):
        """displays peaks generated with peak()"""
        import Display.testplot as testplot
        plot = testplot.plot()
        
        if zoom:
            z0=zoom[0]
            z1=zoom[1]
            if z1 == -1 : z1=self.size1
            pk = self.peaks[np.where(np.logical_and(self.peaks>=z0,self.peaks<=z1))]
        else:
            z0 = 0
            z1 = self.size1
            pk = self.peaks
        if axis is None:
            plot.plot(pk-z0,self.buffer[pk], "o")
            if peak_label:
                for p in pk:
                    plot.text(1.05*self.buffer[p], "%.2f"%axis[p])
        else:
            plot.plot(axis[pk], self.buffer[pk], "o")
            if peak_label:
                for p in pk:
                    plot.text(axis[p], 1.05*self.buffer[p], "%.2f"%axis[p])
        if show: plot.show()
        return self

########################################################################
    def fastclean(self, nsigma=2.0, nbseg=20, axis=0):
        """
        set to zeros all points below nsigma times the noise level
        This allows the corresponding data-set, once stored to file, to be considerably more compressive.
        
        nsigma: float
            the ratio used, typically 1.0 to 3.0 (higher compression)
        nbseg: int
            the number of segment used for noise evaluation, see util.signal_tools.findnoiselevel
        axis: int
            the axis on which the noise is evaluated, default is fastest varying dimension
        """
        from util.signal_tools import findnoiselevel
        todo = self.test_axis(axis)
        if self.dim == 1:
            noise = findnoiselevel(self.get_buffer(), nbseg=nbseg)
            self.zeroing(nsigma*noise)
        elif self.dim == 2:
            if todo == 2:
                for i in xrange(self.size1):
                    self.set_row(i, self.row(i).fastclean(nsigma=nsigma, nbseg=nbseg))
            elif todo == 1:
                for i in xrange(self.size2):
                    self.set_col(i, self.col(i).fastclean(nsigma=nsigma, nbseg=nbseg))
        else:
            raise NPKError("a faire")
        return self
########################################################################
    def sg(self, window_size, order, deriv=0, axis=0):
        """applies saviski-golay of order filter to data
        window_size : int
            the length of the window. Must be an odd integer number.
        order : int
            the order of the polynomial used in the filtering.
            Must be less than `window_size` - 1.
        deriv: int
            the order of the derivative to compute (default = 0 means only smoothing)
        axis: int
            the axis on which the filter is to be applied, default is fastest varying dimension
        """
        import Algo.savitzky_golay as sgm
        todo = self.test_axis(axis)
        m = sgm.sgolay_coef(window_size, order, deriv=0)
        if self.dim == 1:
            self.buffer = sgm.sgolay_comp(self.buffer, m, window_size)
        elif self.dim == 2:
            if todo == 2:
                for i in xrange(self.size1):
                    self.buffer[i,:] = sgm.sgolay_comp(self.buffer[i,:], m, window_size)
            elif todo == 1:
                for i in xrange(1,self.size2):
                    self.buffer[:,i] = sgm.sgolay_comp(self.buffer[:,i], m, window_size)
        else:
            raise NPKError("a faire")
        return self
########################################################################
    def sg2D(self, window_size, order, deriv=None):
        """applies a 2D saviski-golay of order filter to data
        window_size : int
            the length of the square window. Must be an odd integer number.
        order : int
            the order of the polynomial used in the filtering.
            Must be less than `window_size` - 1.
        deriv: None, 'col', or 'row'.   'both' mode does not work.
            the direction of the derivative to compute (default = None means only smoothing)
        can be applied to a 2D only.
        """
        import Algo.savitzky_golay as sgm
        self.check2D()
        self.buffer[:] = sgm.savitzky_golay2D(self.buffer[:], window_size, order, derivative=None)
        return self
########################################################################
    #-------------------------------------------------------------------------------
    def itop(self,axis,value):
        todo = self.test_axis(axis)
        if self.dim == 1:
            return self.axis1.itop(value)
        elif self.dim == 2:
            if todo == 1:
                return self.axis1.itop(value)
            elif todo == 2:
                return self.axis2.itop(value)
    def htop(self,axis,value):
        todo = self.test_axis(axis)
        if self.dim == 1:
            return self.axis1.htop(value)
        elif self.dim == 2:
            if todo == 1:
                return self.axis1.htop(value)
            elif todo == 2:
                return self.axis2.htop(value)
    def htoi(self,axis,value):
        todo = self.test_axis(axis)
        if self.dim == 1:
            return self.axis1.htoi(value)
        elif self.dim == 2:
            if todo == 1:
                return self.axis1.htoi(value)
            elif todo == 2:
                return self.axis2.htoi(value)
    def ptoh(self,axis,value):
        todo = self.test_axis(axis)
        if self.dim == 1:
            return self.axis1.ptoh(value)
        elif self.dim == 2:
            if todo == 1:
                return self.axis1.ptoh(value)
            elif todo == 2:
                return self.axis2.ptoh(value)
    def ptoi(self,axis,value):
        todo = self.test_axis(axis)
        if self.dim == 1:
            return self.axis1.ptoi(value)
        elif self.dim == 2:
            if todo == 1:
                return self.axis1.ptoi(value)
            elif todo == 2:
                return self.axis2.ptoi(value)
    def itoh(self,axis,value):
        todo = self.test_axis(axis)
        if self.dim == 1:
            print "going for F1 , on ", value
            return self.axis1.itoh(value)
        elif self.dim == 2:
            if todo == 1:
                return self.axis1.itoh(value)
            elif todo == 2:
                return self.axis2.itoh(value)

class NPKDataTests(unittest.TestCase):
    """ - Testing NPKData basic behaviour - """
    name1D = "../DATA_test/proj.gs1"
    name2D = "../DATA_test//dosy-cluster2.gs2"       # Byteorder = big_endian
    def test_fft(self):
        " - Testing FFT methods - "
        print self.test_fft.__doc__
        x = np.arange(1024.0)
        E=NPKData(buffer=x)
        E.axis1.itype = 0
        E.ifftr()
        sim = np.exp(-x*(0.1+1j))
        D=NPKData(buffer=sim)
        p = 0 # 56.78                # surprisingly,  tests fail by 10-1/10-4 if p!=0.0 
        D1 = D.copy().phase(p,0).chsize(2*D.size1).fft().real()     #.display(label='fft')
        D2 = D.copy().phase(p,0).fftr()                             #.display(new_fig=False,label='fftr')
        D3 = D2.ifftr()                                             #.display(show=True,new_fig=False)
        D.phase(p,0).add(D3.mult(-1))                               #.display(show=True,label='difference sur les FID')
        print "ecart max :",   np.max(D.get_buffer())
        print "ecart moyen :", np.mean(D.get_buffer())
        self.assertTrue(np.max(D.get_buffer()) < 1E-14)
        self.assertTrue(np.mean(D.get_buffer()) < 1E-14)

    def test_load(self):
        " - Testing load methods"
        print self.test_load.__doc__
        E=NPKData(name=self.name1D)
        self.assertAlmostEqual(E[0], 1869.4309082)
        self.assertAlmostEqual(E.get_buffer().max(), 603306.75)

    def test_math(self):
        " - Testing math methods - "
        M=np.zeros((20,20))
        d = NPKData(buffer=M)
        d[5,7]=10
        d[10,12]=20
        d += 1
        self.assertAlmostEqual(d[10,12] , 21)
        e = d+2
        self.assertAlmostEqual(e[10,12] , 23)
        e = d-2
        self.assertAlmostEqual(e[10,12] , 19)
        f = 2*d+e
        self.assertAlmostEqual(f[10,12] , 2*21+19)
        f += d
        self.assertAlmostEqual(f[10,12] , 3*21+19)
        f *= 3.0
        self.assertAlmostEqual(f[10,12] , 3*(3*21+19))
        f = 2*d-e
        self.assertAlmostEqual(f[10,12] , 2*21-19)
        re = e.row(10)
        re.axis1.itype = 1
        re *= 2j
        self.assertAlmostEqual(re[13], 2*19.0)
    def test_flatten(self):
        " test the flatten utility "
        self.assertEqual( [1,2,3,4,5,6,7], flatten( ( (1,2), 3, (4, (5,), (6,7) ) ) ))
    
    def test_peaks2d(self):
        "test 2D peak picker"
        print self.test_peaks2d.__doc__
        M=np.zeros((30, 30))
        M[5,7] = 20
        M[10,12] = 20
        d = NPKData(buffer = M)
        thresh = 10
        x,y = d.peaks2d(threshold = thresh)
        print "hou 2D",list(x),list(y)
        self.assertEqual(list(x) , [ 5, 10])
        self.assertEqual(list(y) , [ 7, 12])
        
    def test_peaks1d(self):
        "test 1D peak picker"
        print self.test_peaks1d.__doc__
        M = np.zeros((30))
        M[5] = 20
        M[7] = 8
        M[15] = 11
        M[10] = 20
        d = NPKData(buffer = M)
        thresh = 10
        x = d.peak(threshold = thresh)
        self.assertEqual(list(x) , [ 5, 10, 15])
        #self.assertEqual(list(y) , [ 20, 20, 11])
        
    def test_dampingunit(self):
        "test itod and dtoi"
        print self.test_dampingunit.__doc__
        LA = LaplaceAxis( dmin = 10.0, dmax = 10000.0, dfactor = 35534.34765625)
        damping = LA.itod(50)
        point = LA.dtoi(damping)
        self.assertAlmostEqual(damping, 2154.43469003)
        self.assertAlmostEqual(point, 50)

    def test_zf(self):
        for ty in range(2):  # check both real and complex cases, and 1D and 2D (in F1)
            d1 = NPKData(buffer = np.arange(6.0))
            d2 = NPKData(buffer = np.zeros((6, 8)))
            for i in range(d2.size2):
                d2.set_col(i,  NPKData(buffer = 0.1*i + np.arange((d2.size1))) )
            # d2[i,j] = i + 0.1*j
            d1.axis1.itype = ty
            d2.axis1.itype = ty
            d2.axis2.itype = 1
            if ty == 0:
                samp = np.array([0,3,6,8,10,11])   # real
            else:
                samp = np.array([0,3,5])   # complex   - 
            d1.axis1.set_sampling(samp)
            d2.axis1.set_sampling(samp)
            self.assertTrue(d1.axis1.sampled)
            self.assertTrue(d2.axis1.sampled)
            d1.zf()
            d2.zf()
            if ty == 0:
                self.assertAlmostEqual(d1[6],2.0)
                self.assertAlmostEqual(d1[7],0.0)
                self.assertAlmostEqual(d2[6,4],2.4)
                self.assertAlmostEqual(d2[7,4],0.0)
            else:
                self.assertAlmostEqual(d1[6],2.0)
                self.assertAlmostEqual(d1[7],3.0)
                self.assertAlmostEqual(d2[6,4],2.4)
                self.assertAlmostEqual(d2[7,4],3.4)
    
    def test_hypercomplex_modulus(self):
        '''
        Test of hypercomplex modulus
        '''
        arr = np.array([[1, 4],[3, 7],[1, 9],[5, 7]]) # hypercomplex of size 2x2
        modulus = hypercomplex_modulus(arr, 2, 2) 
        np.testing.assert_almost_equal(modulus, np.array([[np.sqrt(75)],[np.sqrt(156)]])) # 

if __name__ == '__main__':
#     x = np.arange(1024.0)
#     sim = np.exp(-x*(0.1+1j))
#     D=NPKData(buffer=sim)
#     D.fft()
#     print D
# #    D.display(show=True,zoom=(100,400))
#     E=NPKData(buffer=np.ones((1000,2000)))
#     E.axis1.itype=1
#     E.axis2.itype=1
#     E.apod_sin(maxi=0.5,axis=1).apod_sin(maxi=0.5,axis=2)
#     print E
#     E.display(show=True,zoom=((100,500),(100,1000)))
    unittest.main()
