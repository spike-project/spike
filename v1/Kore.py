#!/usr/bin/env python 
# encoding: utf-8

"""
Kore.py

Created by Marie-Aude Coutouly on 2010-03-26.
Copyright (c) 2010 NMRTEC. All rights reserved.
"""
from __future__ import print_function

import numpy as np
from .. import NPKData as npkd
from ..File import GifaFile as gf
import array
import sys
import inspect

###########################################################################
class Kore(object):
    def __init__(self, debug=0):
        
        self.debug = debug
        self._column = npkd.NPKData(dim = 1)
        self._plane2d = npkd.NPKData(dim = 2)
        self._image = npkd.NPKData(dim = 3)
        self._datab = npkd.NPKData(dim = 1)
        self._window = npkd.NPKData(dim = 1)
        self._filter = npkd.NPKData(dim = 1)
        self._tab = npkd.NPKData(dim = 1)
        self._last_row = 0
        self._last_plane = 0
        self._last_col = 0
        self._last_ph0 = 0
        self._last_ph1 = 0
        self._shift = 0.0
        self._noise = 0.0
        self.peaks = []
        
        self.dim(1)
    def _getcurrent(self):
        """
        getter for _current
        """
        if self._dim == 1:
            return self._column
        elif self._dim == 2:
            return self._plane2d
        elif self._dim == 3:
            return self._image
    def _setcurrent(self,npkdata):
        """
        setter for _current
        """
        if self._dim == 1:
            npkdata.check1D()
            self._column = npkdata
        elif self._dim == 2:
            npkdata.check2D()
            self._plane2d = npkdata
        elif self._dim == 3:
            npkdata.check3D()
            self._image = npkdata
        else:
            raise Exception("Kore internal problem")
    _doc_current = "doc"
    _current = property(_getcurrent,_setcurrent,None,_doc_current)
    #--------------------------------------------------------------------------
    def report(self):
        """print a summary of the internal state of the kernel"""
        report = """
        buffer 1D : %s
        buffer 2D : %s
        buffer 3D : %s
        current   : %s
        """%(self._column.report(),self._plane2d.report(),self._image.report(),self._current.report())
        print(report)
    #-------------------------------------------------------------------------
    def status(self):
        """
        print a summary of the internal state of the kernel
        """
        d = self.get_dim()
        if d == 1:
            d1 = "- current working buffer"
        else:
            d1 = ""
        if d == 2:
            d2 = "- current working buffer"
        else:
            d2 = ""
        if d == 3:
            d3 = "- current working buffer"
        else:
            d3 = ""
        self.dim(1)
        self.com_max()
        
        report = """
    DIM 1 %s
    =====
       buffer size : %i - itype %i
       values from : %f to %f
       Spectral width : %f
       %i peak(s) in database

    """%(d1,self.get_si1_1d(),self.get_itype_1d(),self.geta_max(2),self.geta_max(1),self.get_specw_1d(),self.get_npk1d())

        self.dim(2)
        self.com_max()
        report = report+"""
    DIM 2 %s
    =====
        buffer sizes : %i x %i - itype %i
         values from : %f to %f
    Spectral widthes : %f x %f
    %i peak(s) in database

    """%(d2,self.get_si1_2d(),self.get_si2_2d(),self.get_itype_2d(),self.geta_max(2),self.geta_max(1),self.get_specw_1_2d(),self.get_specw_2_2d(),self.get_npk2d())
        self.dim(3)
        self.com_max()
        report = report+"""
    DIM 3 %s
    =====
        buffer sizes : %i x %i x %i - itype %i
         values from : %f to %f
    Spectral widthes : %f x %f x %f
    %i peak(s) in database

    """%(d3,self.get_si1_3d(),self.get_si2_3d(), self.get_si3_3d(),self.get_itype_3d(),self.geta_max(2),self.geta_max(1),self.get_specw_1_3d(),self.get_specw_2_3d(),self.get_specw_3_3d(),self.get_npk3d())
       
        self.dim(d)
        return  report
    #--------------------------------------------------------------------------
    def dim(self,d):
        """
        Declaration of the ._current buffer
        """
        
        if d not in (1,2,3):
            raise Exception("dim can only take the value : 1 2 3")
        self._dim = d
        
    #--------------------------------------------------------------------------
    def checknD(self,n):
        if self._dim != n:
            raise Exception("The buffer set is not a %1dD experiment, as required"%n)
    #--------------------------------------------------------------------------
    def check1D(self):
        "true for a 1D"
        self.checknD(1)
    def check2D(self):
        "true for a 2D"
        self.checknD(2)
    def check3D(self):
        "true for a 3D"
        self.checknD(3)
    #--------------------------------------------------------------------------
    def power2(self, i):
        """
        Compute the power of 2 that is under or equal to i 
        """
        import math as m
        return int(m.pow(2,m.floor(m.log(i)/m.log(2))))
    #--------------------------------------------------------------------------
    def _test_axis(self,axis):
        """
        takes axis as
            F1 f1 F2 f2 F12 or f12 in 2D
        or  F1 f1 F2 f2 F3 f3 F12 f12 F13 f13 F23 f23 F123 or f123 in 3D
        and return a list of axis to process :
        1 / 2 / 3
        used by axis relative commands
        """
        ret = []
        if axis in ("F3","F13","F23","F123", "f3","f13","f23","f123",3):
            ret.append(3)
        if axis in ("F2","F12","F23","F123", "f2","f12","f23","f123",2):
            ret.append(2)
        if axis in ("F1","F12","F13","F123","f1", "f12","f13","f123",1):
            ret.append(1)
        if (ret == []):
            raise Exception("Error with axis")
        return ret
    #--------------------------------------------------------------------------
    def addbase(self,constant):
        """
        Removes a constant to the data. The default value is the value of 
        SHIFT (computed by EVALN).
        
        see also : bcorr evaln shift
        """
        #print "A VERIFIER!!!"
        self._current.addbase(constant)
    #--------------------------------------------------------------------------
    def addnoise(self,noise,seed=0):
        """
        add to the current data-set (1D, 2D, 3D) a white-gaussian, 
        characterized by its level noise, and the random generator seed.
        """
        self._current.addnoise(noise,seed)
    #--------------------------------------------------------------------------
    def adddata(self,debug = False):
        """
        Add the contents of the DATA buffer to the current data-set. 
        Equivalent to ADD but in-memory.
        """
        #print "A VERIFIER!!!"
        if self._datab.buffer.shape == self._current.buffer.shape:
            self._current.buffer += self._datab.buffer
        else:
            raise Exception("diff sizes ", self._datab.buffer.shape, self._current.buffer.shape)
    
    #--------------------------------------------------------------------------
    def mult1d(self, axis = 0):
        """
        multiply the current 2D or 3D with the contents of the 1d buffer considered as a f1(i)f2(j) concatenated buffer
        
        
        see also : multdata add adddata filter
        """
        for ax in self._test_axis(axis):
            if ax == 2:
                c = npkd.NPKData(buffer = self._column.buffer)
                #c.display(label ="row")
                c.chsize(self.get_si1_2d() + self.get_si2_2d())
                c.buffer[self.get_si2_2d():self.get_si1_2d() + self.get_si2_2d()] = 1.0
                
                
                for j in range (self.get_si1_2d()):
                    a = c.buffer[self.get_si2_2d() +j]
                    for i in range(self.get_si2_2d()):
                        self._current.buffer[j,i] = self._current.buffer[j,i]* c.buffer[i]* a
                        #print self._current.buffer[i,j]* c.buffer[i]* a
            elif ax == 1:
                c = npkd.NPKData(buffer = self._column.buffer)
                
                c.chsize(self.get_si1_2d() + self.get_si2_2d())
                c.buffer[0:self.get_si2_2d() ] = 1.0
                c.buffer[self.get_si2_2d():self.get_si1_1d() + self.get_si2_2d()] = self._column.buffer[:]
                #c.display()
                for j in range (self.get_si1_2d()):
                    a = c.buffer[self.get_si2_2d() +j]
                    for i in range(self.get_si2_2d()):
                        #print self._current.buffer[i,j],c.buffer[i],a
                        self._current.buffer[j,i] = self._current.buffer[j,i]* c.buffer[i]* a
                        
    #--------------------------------------------------------------------------
    def mult(self,constant):
        #print "A VERIFIER!!!"
        self._current = self._current.mult(constant)
    #---------------------------------------------------------------------------
    def multdata(self):
        """
        Multiplies point by point, the content of the current working buffer 
        with the content of the DATA buffer. Permits to realize convolution 
        product. Works in 1D, 2D, in real, complex and hypercomplex modes.

        see also : ADDDATA MINDATA MAXDATA EXCHDATA MULT PUT
        """
        if self._current.dim != self._datab.dim:
            raise Exception("wrong buffer dim: %s"%str(self._current.dim()))
        if self._current.dim == 1:
            if self.get_si1_1d() != self._datab.size1:
                raise Exception("wrong buffer size: %s %s"%(str(self.get_si1_1d()),str(self._datab.size1)))
        elif self._current.dim == 2:
            if self.get_si1_2d() != self._datab.size1 or self.get_si2_2d() != self._datab.size2:
                raise Exception("wrong buffer size 2D : %s"%str(self.get_si1_2d()))
        elif self._current.dim == 3:
            if self.get_si1_3d() != self._datab.size1 or self.get_si2_3d() != self._datab.size2 or self.get_si3_3d() != self._datab.size3:
                raise Exception("wrong buffer size 3D : %s"%str(self.get_si1_3d()))
        self._current.buffer *= self._datab.buffer
        return self
    #-------------------------------------------------------------------------------
    def maxdata(self):
        """
        Compare the content of the current buffer with the content of the 
        DATA buffer, and leave in memory the largest of the 2 values. 
        Usefull for projections or symetrisation macros.

        see also : mindata exchdata adddata multdata sym put
        """
        if self._current.dim != self._datab.dim:
            raise Exception("wrong buffer dim: %i %i"%(self._current.dim,self._datab.dim))
        if self._current.dim == 1:
            if self.get_si1_1d() != self._datab.size1:
                raise Exception("wrong buffer size: %s %s"%(str(self.get_si1_1d()),str(self._datab.size1)))
        elif self._current.dim == 2:
            if self.get_si1_2d() != self._datab.size1 or self.get_si2_2d() != self._datab.size2:
                raise Exception("wrong buffer size 2D : %s"%str(self.get_si1_2d()))
        elif self._current.dim == 3:
            if self.get_si1_3d() != self._datab.size1 or self.get_si2_3d() != self._datab.size2 or self.get_si3_3d() != self._datab.size3:
                raise Exception("wrong buffer size 3D : %s"%str(self.get_si1_3d()))
        self._current.buffer = np.maximum(self._current.buffer,self._datab.buffer)
        return self
    #-------------------------------------------------------------------------------
    def mindata(self):
        """
        Compare the content of the current buffer with the content of the 
        DATA buffer, and leave in memory the smallest of the 2 values. 
        Usefull for projections or symetrisation macros.

        see also : maxdata exchdata adddata multdata sym put
        """
        if self._current.dim != self._datab.dim:
            raise Exception("wrong buffer dim: %i %i"%(self._current.dim,self._datab.dim))
        if self._current.dim == 1:
            if self.get_si1_1d() != self._datab.size1:
                raise Exception("wrong buffer size: %s %s"%(str(self.get_si1_1d()),str(self._datab.size1)))
        elif self._current.dim == 2:
            if self.get_si1_2d() != self._datab.size1 or self.get_si2_2d() != self._datab.size2:
                raise Exception("wrong buffer size 2D : %s"%str(self.get_si1_2d()))
        elif self._current.dim == 3:
            if self.get_si1_3d() != self._datab.size1 or self.get_si2_3d() != self._datab.size2 or self.get_si3_3d() != self._datab.size3:
                raise Exception("wrong buffer size 3D : %s"%str(self.get_si1_3d()))
        self._current.buffer = np.minimum(self._current.buffer,self._datab.buffer)
        return self
    #--------------------------------------------------------------------------
    def modulus(self):
        #print "A VERIFIER!!!"
        self._current = self._current.modulus()
    #--------------------------------------------------------------------------
    def itype(self,value):
        #print "A VERIFIER!!!"
        if self._dim == 1:
            if value not in (0,1):
                raise Exception("wrong value for itype")
            self._column.axis1.itype = value
        elif self._dim == 2:
            if value not in (0,1,2,3):
                raise Exception("wrong value for itype")
            self._plane2d.axis1.itype = value/2
            self._plane2d.axis2.itype = value%2
        elif self._dim == 3:
            if value not in (0,1,2,3,4,5,6,7):
                raise Exception("wrong value for itype")
            self._image.axis1.itype = value/4
            self._image.axis2.itype = (value/2)%2
            self._image.axis3.itype = value%2
    #--------------------------------------------------------------------------
    def col(self,i):
        self.check2D()
        self._column = self._plane2d.col(i-1)
        self._last_col = i
    #---------------------------------------------------------------------------
    def diag(self,direc = "F12"):
        if self._dim == 2:
            self._column = self._plane2d.diag()
        elif self._dim == 3:
            self._plane2d = self._image.diag(direc)
    #---------------------------------------------------------------------------
    def one(self):
        if self._dim == 1:
            if self._current.axis1.itype == 0:
                self._current.buffer = np.ones_like(self._current.buffer)
            else:
                self._current.buffer[::2] = 1.0
                self._current.buffer[1::2] = 0.0
        elif self._dim == 2:
            if self._current.axis1.itype == 0:
                 if self._current.axis2.itype == 0:
                     self._current.buffer = np.ones_like(self._current.buffer)
                 else :             # itype is 1 for axis2
                     self._current.buffer[:,::2] = 1.0
                     self._current.buffer[:,1::2] = 0.0
            else:                   # itype is 1 for axis1
                if self._current.axis2.itype == 0:
                     self._current.buffer[::2] = 1.0
                     self._current.buffer[1::2] = 0.0
                else :             # itype is 1 for axis2
                     self._current.buffer[::2] = 1.0
                     self._current.buffer[1::2] = 0.0
                     self._current.buffer[:,1::2] = 0.0
        elif self._dim == 3: 
            print("ONE 3D A FINIR")
            if self._current.axis1.itype == 0:
                 if self._current.axis3.itype == 0:
                     self._current.buffer = np.ones_like(self._current.buffer)
                 else :             # itype is 1 for axis3
                     self._current.buffer[:,:,::2] = 1.0
                     self._current.buffer[:,:,1::2] = 0.0
            elif self._current.axis2.itype == 0:                   # itype is 1 for axis2
                if self._current.axis1.itype == 0:
                     self._current.buffer[:,::2] = 1.0
                     self._current.buffer[:,1::2] = 0.0
                else :             # itype is 1 for axis2
                     self._current.buffer[:,::2] = 1.0
                     self._current.buffer[:,1::2] = 0.0
                     self._current.buffer[:,:,1::2] = 0.0
    #---------------------------------------------------------------------------
    def bcorr(self, mode, *arg): #i,radius,axis, points):
        """Apply a baseline correction
        Computes and applies a base-line correction to the current data set.
        mode   describe the algorithm used:
          *  1 is linear correction
          *  2 is cubic spline correction.
          *  3 is polynomial (and related) correction  NOT IMPLEMENTED YET !

          if mode == 1 or 2

          then in 1D *arg is radius, list_of_points
            or in 2D *arg is radius, axis, list_of_points

        axis in 2D is either f1 or f2 (dimension in which correction is  applied).
        radius   is the radius around which each pivot point is averaged.
        
         list_of_points is then the list of the pivot points used for the 
        base-line correction.
        Linear correction can use 1 or more pivot points. 1 point 
        corresponds to correction of a continuous level. Spline corrections 
        needs at least 3 points.
        In any case maximum is 100 pivot points.
        """

        if mode in (1,2):
            radius = arg[0]
            if self.get_dim() == 1:
                if self.get_itype_1d() != 0:
                    raise Exception("not available")
                xpoints = arg[1]
                if mode == 1:
                    #print "bcorr mode 1 a verifier"
                    self._current.linear_interpolate(xpoints)
                elif mode == 2:
                    self._current.spline_interpolate(xpoints, kind=3)
            elif self.get_dim() == 2:
                axis = arg[1]
                xpoints = arg[2]
                if mode == 1:
                    print("attention, c'est faux")
                    self._current.linear_interpolate(xpoints, axis=axis)
                elif mode == 2:
                    self._current.spline_interpolate(xpoints, axis=axis, kind=3)
            else:
                raise Exception("not implemented")
        elif mode == 3:
            print(" bcorr 3 a faire")
    def bcorrp1(self):
        print("bcorrp1 a faire")
    def bcorrp0(self):
        print("bcorrp0 a faire")
        # from scipy import interpolate
        # 
        # A = arg[0]
        # R = arg[1]
        # if len(args) == 4:
        #     axis = arg[2]
        #     points = arg[3]
        # else:
        #     points = arg[2]
        #     
        # 
        # if axis :
        #     for ax in self._test_axis(axis):
        #         print "from bcorr",ax
        #         # if ax == "f2":
        #         #     for row in range (self._current.axis1.size):
        #         #         values = []
        #         #         for point in points:
        #         #             value = 0
        #         #             for j in range(point-int(radius),point+int(radius)+1):
        #         #                 value = value + self._current.buffer[row,j]
        #         #             
        #         #             values.append(value/(2*int(radius) +1))
        #         #         tck = interpolate.splrep(points,values, k = 2,s = 0)
        #         #         xnew = xrange(self._current.axis2.size)
        #         #         ynew = interpolate.splev(xnew,tck,der=0)
        #         #         
        #         #         self._current.buffer[row,:] -= ynew[:]
        #         # if ax == "f1":
        #         #     for col in range (self._current.axis2.size):
        #         #         self._current.buffer[:,col] = self._current.buffer[:,col] - self._current.buffer[point,col]
        #     # self.write("/Users/mac/Desktop/After_bcorr.gs2")
        # else:
        #     for point in points:
                
#        avec scipy.interpolate !
   
    #---------------------------------------------------------------------------
    def chsize(self,*args):
        """
        Change size of data, zero-fill or truncate.
        DO NOT change the value of OFFSET and SPECW, so EXTRACT should
        always be preferred on spectra (unless you know exactly what your are doing).


        see also : extract modifysize
        """
        self._current.chsize(*args)
    #---------------------------------------------------------------------------
    def modifysize(self, si1, si2=-1, si3=-1):
        """
        modifysize( si1, si2 )
        modifysize( si1, si2, si3 )

        Permits to modify the leading sizes of a 2D or a 3D data-set, 
        provided the product of the sizes : si1*si2{*si3} is equal to the 
        product of the old ones.
        
        Does not actually modify the data.
        
        see also : chsize
        """
        print("modifisize A VALIDER")
        if self._dim == 2:
            self._plane2d.chsize(si1, si2)
        elif self._dim == 3:
            self._image.chsize(si1, si2, si3)
        self._current.adapt_size()
    #---------------------------------------------------------------------------
    def plus(self):
        self._current.plus()
    #---------------------------------------------------------------------------
    def minus(self):
        self._current.minus()
    #---------------------------------------------------------------------------
    def extract(self,*args):
        self._current.extract(*args)
    #---------------------------------------------------------------------------
    def noise(self,value):
        """
        Contains the level of noise in the data-set. When loading data (1 or 
        2D) the noise level is evaluated automatically from the last 10th of 
        the data. Can also be set with EVALN.
        Used by INTEG and by Maximum Entropy run.
        """
        self._noise = value
    #---------------------------------------------------------------------------
    def shift(self,value):
        """
        This context holds the systematic baseline shift of the current 
        data-set, computed automatically by EVALN.
        Used by INTEG.
        see also : evaln noise addbase
        """
        self._shift = value
    #---------------------------------------------------------------------------
    def evaln(self,a,b,c=-1,d=-1):
        """
        evaluates the noise level as well as the overall offset of the data,over a area of the data. The results are stored in the NOISE and SHIFT 
        contexts This command is called automatically whenever a data set is read. The command will prompt for the last selected region with the POINT command
        
        in 2D, a,b,c,d is llf1, llf2, ur1, ur2
        """
        if self._dim == 1:
            shift = self._current.mean((a,b))
            noise = self._current.std((a,b))
        elif self._dim == 2:
            shift = self._current.mean( ((a,c),(b,d)) )
            noise = self._current.std( ((a,c),(b,d)) )
        else:
            raise Exception('Not available in 3D')
        self._shift = shift
        self._noise = noise
    #---------------------------------------------------------------------------
    def fill(self,value):
        self._current.fill(value)
    #---------------------------------------------------------------------------
    def pkclear(self):
        self.peaks = []
    #---------------------------------------------------------------------------
    def peak(self, pkradius = 0):
        if self._dim == 1:
            self.peaks = self._current.peak(threshold = self.mini, offset = None)
            self.index = np.where(self.peaks <= self.maxi)[0]
            
            # for peak in peaks:
            #     if self._current.buffer[peak] <= self.maxi :
            #         self.peaks.append(peak)
        elif self._dim == 2:
            self.peaks2d = self._current.peaks2d(threshold = self.mini, zoom = ((self._current.zo_2_1l,self._current.zo_2_1m),(self._current.zo_2_2l,self._current.zo_2_2m)))
        else:
            print("3D Pick peaker is not yet written")
    #---------------------------------------------------------------------------
    def geta_pk1d_a(self,i):
        return self._current.buffer[self.peaks[self.index[i-1]]]
    def geta_pk1d_a_err(self,i):
        return 0
    #---------------------------------------------------------------------------
    def geta_pk1d_f(self,i):
        return self.peaks[self.index[i-1]]
    def geta_pk1d_f_err(self,i):
        return 0
    #---------------------------------------------------------------------------
    def geta_pk1d_p(self,i):
        return 0
    #---------------------------------------------------------------------------
    def geta_pk1d_t(self,i):
        return 0
    #---------------------------------------------------------------------------
    def geta_pk1d_w(self,i):
        return 0
    def geta_pk1d_w_err(self,i):
        return 0
    #---------------------------------------------------------------------------
    def geta_pk2d_a(self,i):
        return self._current.buffer[self.peaks2d[0][i-1],self.peaks2d[1][i-1]]
    def geta_pk2d_a_err(self,i):
        return 0
    #---------------------------------------------------------------------------
    def geta_pk2d_f1f(self,i):
        return self.peaks2d[0][i-1]
    def geta_pk2d_f1f_err(self,i):
        return 0
    #---------------------------------------------------------------------------
    def geta_pk2d_f1w(self,i):
        return 0
    def geta_pk2d_f1w_err(self,i):
        return 0
    #---------------------------------------------------------------------------
    def geta_pk2d_f2f(self,i):
        return  self.peaks2d[1][i-1]
    def geta_pk2d_f2f_err(self,i):
        return 0
    #---------------------------------------------------------------------------
    def geta_pk2d_f2w(self,i):
        return 0
    def geta_pk2d_f2w_err(self,i):
        return 0
    #---------------------------------------------------------------------------
    def geta_pk3d_a(self,i):
        return 0
    #---------------------------------------------------------------------------
    def geta_pk3d_f1f(self,i):
        return 0
    #---------------------------------------------------------------------------
    def geta_pk3d_f1w(self,i):
        return 0
    #---------------------------------------------------------------------------
    def geta_pk3d_f2f(self,i):
        return 0
    #---------------------------------------------------------------------------
    def geta_pk3d_f2w(self,i):
        return 0
    #---------------------------------------------------------------------------
    def geta_pk3d_f3f(self,i):
        return 0
    #---------------------------------------------------------------------------
    def geta_pk3d_f3w(self,i):
        return 0
    #---------------------------------------------------------------------------
    def freq(self,*args):
        """
        
        The context FREQ holds the basic frequency of the spectrometer (in MHz).
        freq_H1 is meant to be the basic frequency of the spectrometer (1H freq)
        and is not used in the program. freq2 (and freq1 in 2D) are the freq
        associated to each dimension (different if in heteronuclear mode).
        Values are in MHz.
        
        see also : specw offset
        """
        if self._dim  == 1 :
            self.freq1d(*args)
        elif self._dim == 2 :
            self.freq2d(*args)
        elif self._dim == 3 :
            self.freq3d(*args)
        else:
            print("This should never happen", self._dim)
    #---------------------------------------------------------------------------
    def freq1d(self,freq_h1, freq1):
        self._column.frequency = freq_h1
        self._column.axis1.frequency =  freq1
    #---------------------------------------------------------------------------
    def freq2d(self,freq_h1, freq1, freq2):
        self._plane2d.frequency = freq_h1
        self._plane2d.axis1.frequency = freq1
        self._plane2d.axis2.frequency = freq2
    #---------------------------------------------------------------------------
    def freq3d(self,freq_h1, freq1, freq2,freq3):
        self._image.frequency = freq_h1
        self._image.axis1.frequency = freq1
        self._image.axis2.frequency = freq2  
        self._image.axis3.frequency = freq3      
    #---------------------------------------------------------------------------
    def ft(self,axis="F1"):
        """
        Performs in-place complex Fourier Transform on the current data-set; 
        Data-set must be Complex.
         
         
        All FT commands work in 1D, 2D or 3D
         
        <ul>
        <li> in 1D axis, is not needed
        <li> in 2D axis, is F1, F2 or F12
        <li> in 3D axis, is F1, F2, F3, F12, F13, F23 or F123
        </ul>
         
         
        Here is a complete overview of FT routines : C stands for Complex, R  stands for Real
        <pre>
         	FIDs		Spectra
         	C	---FT--->	C
         	C	<--IFT---	C
         	R	--RFT-->	C
         	R	<--IRFT--	C
         	C	-FTBIS->	R
         	C	<-IFTBIS-	R
         	R	Does not exist 	R
        </pre>
        """
        #print "A TESTER"
        for ax in self._test_axis(axis):
            self._current.fft(ax)
    #---------------------------------------------------------------------------

    def ift(self,axis="F1"):
        """
        Performs in-place inverse complex Fourier Transform on the current data-set; 
        Data-set must be Complex.
        """
        #print "A TESTER"
        for ax in self._test_axis(axis):
            print(ax)
            self._current.ifft(ax)
    #---------------------------------------------------------------------------
    def ftbis(self,axis="F1"):
        """
        Data-set must be Complex.
        """
        #print "A TESTER"
        for ax in self._test_axis(axis):
            self._current.fftr(ax)
    #---------------------------------------------------------------------------
    def iftbis(self,axis="F1"):
        """
        Data-set must be Real.
        """
        #print "A TESTER"
        for ax in self._test_axis(axis):
            self._current.ifftr(ax)
    #---------------------------------------------------------------------------
    def com_max(self):
        self.maxi = self._current.buffer.max()          # built in function
        self.mini = self._current.buffer.min()          # built in function
        self.imaxi = self._current.buffer.argmax()       # built in function
        self.imini = self._current.buffer.argmin()       # built in function
        return self
    #---------------------------------------------------------------------------
    def minimax(self,mini,maxi):   
        self.mini = mini
        self.maxi = maxi
        return self
    #---------------------------------------------------------------------------
    def geta_max(self,index):
        if index == 1:
            return self.maxi
        elif index == 2:
            return self.mini
        else:
            print("we have a problem")
    #---------------------------------------------------------------------------
    def dmax(self,value):
        """
        Determines the fastest decaying component during Laplace analysis
        Given in arbitrary unit, use DFACTOR to relate to actual values.
        
        see also : dmin dfactor laplace tlaplace invlap invtlap
        
        sets the final value for the laplace transform
        """
        self._current.diffaxis.dmax = value
    #---------------------------------------------------------------------------
    def dmin(self,value):
        """
        Determines the fastest decaying component during Laplace analysis
        Given in arbitrary unit, use DFACTOR to relate to actual values.

        see also : dmin dfactor laplace tlaplace invlap invtlap

        sets the final value for the laplace transform
        """
        self._current.diffaxis.dmin = value
    #--------------------    -------------------------------------------------------
    def dfactor(self,value):
        self._current.diffaxis.dfactor = value
    #---------------------------------------------------------------------------
    def offset(self,*args):
        """
        Permits to specify the offset of the right-most (upper right most in
        2D) point of the data set. The value for offset are changed by @extract
        see also : specw
        """
        if self._dim  == 1 :
            self.offset1d(*args)
        elif self._dim == 2 :
            self.offset2d(*args)
        elif self._dim == 3 :
            self.offset3d(*args)
        else:
            print("This should never happen", self._dim)
    #---------------------------------------------------------------------------
    def offset1d(self, off1): 
        self._column.axis1.offset = off1
    #---------------------------------------------------------------------------
    def offset2d(self, off1, off2): 
        self._plane2d.axis1.offset = off1
        self._plane2d.axis2.offset = off2
    #---------------------------------------------------------------------------
    def offset3d(self, off1, off2, off3): 
        self._image.axis1.offset = off1
        self._image.axis2.offset = off2
        self._image.axis3.offset = off3
    #---------------------------------------------------------------------------
    def read(self,file_name):
        """
        read( file_name )

        Reads the file as the new data set in standard format . 
        Same as readc 

        see also : write 
        """

        F = gf.GifaFile(file_name,"r")
        F.load()
        self.dim(F.dim)
        self._current = F.get_data()    # A VERIFIER !!!!!!
        self._current.check()
        #self = F.get_data() 
    #---------------------------------------------------------------------------
    def real(self,axis='F1'):
        for ax in self._test_axis(axis):
            self._current.real(ax)
    #---------------------------------------------------------------------------
    def reverse(self,axis='F1'):
        for ax in self._test_axis(axis):
            self._current.reverse(ax)
    #---------------------------------------------------------------------------
    def revf(self,axis='F1'):
        """
        Processes FID data-sets by multiplying by -1 2 points out of 4. 
        Permits to preprocess Bruker FIDs in Dim 2 (Bruker trick) before 
        RFT, or permits to bring back zero frequency in the center for some 
        other data formats
        """
        # print "A TESTER"
        for ax in self._test_axis(axis):
            self._current.revf(ax)
    #---------------------------------------------------------------------------
    def invf(self,axis='F1'):
        """
        Process data-sets by multiplying by -1 1 point every 2 points. 
        Equivalent to taking the conjugated on complex data-sets, or 
        hyperconjugated on hypercomplex data-sets. If applied on a complex 
        FID, inverses the final spectrum obtained after Fourier transform. 

        see also : revf itype ft reverse
        """
        #print "A TESTER"
        for ax in self._test_axis(axis):
            self._current.revf(ax)
    #---------------------------------------------------------------------------
    def irft(self,axis="F1"):
        """
        Perform  real-to-complex Fourier Transform on data
        """
        #print "A TESTER"
        for ax in self._test_axis(axis):
            self._current.irfft(ax)
    #---------------------------------------------------------------------------
    def rft(self,axis="F1"):
        """
        Perform  real-to-complex Fourier Transform on data
        """
        for ax in self._test_axis(axis):
            self._current.rfft(ax)
    #---------------------------------------------------------------------------
    def row(self,i):
        """
        Extract the nth 1D row (along F2 axis) from the 2D data-set, and put 
        it in the 1D buffer. The row will be available as a 1D data set when 
        going from 2D to 1D
        """
        self.check2D()
        self._column = self._plane2d.row(i-1)
        self._last_row = i
    #---------------------------------------------------------------------------
    def plane(self,axis,i):
        """
        Extract the nth 1D row (along F2 axis) from the 2D data-set, and put 
        it in the 1D buffer. The row will be available as a 1D data set when 
        going from 2D to 1D
        """
        self.check3D()
        self._plane2d = self._image.plane(axis,i-1)
        self._last_plane = axis
    #---------------------------------------------------------------------------
    def proj(self,axis,projtype):
        
        if self._dim == 2:
            for axes in self._test_axis(axis):
                self._column = self._current.proj(axes,projtype)
            
        return self
    #---------------------------------------------------------------------------
    def vert(self,i,j):
        """
        In 3D mode, extract a column orthogonal to the last displayed plane.
        The column is taken at coordinates i and j in this plane.
        
        see also : plane col row dim
        """
        #print "A TESTER"
        self.check3D()
        if self._last_plane == 1:
            self._column.buffer = self._image.buffer[:,i,j].copy()
        elif  self._last_plane == 2:
            self._column.buffer = self._image.buffer[i,:,j].copy()
        elif  self._last_plane == 3:
            self._column.buffer = self._image.buffer[i,j,:].copy()
    #---------------------------------------------------------------------------
    def setval(self,*args):
        """
        Will set the value of the data point to x. The number of coordinates 
        of the point depends of dim. In dim 2 or 3, coordinates are F1 F2 or 
        F1 F2 F3. Can be usefully used when associated to the functions 
        valnd() to change data point value.
        """
        if self._dim  == 1 :
            self.setval1d(*args)
        elif self._dim == 2 :
            self.setval2d(*args)
        elif self._dim == 3 :
            self.setval3d(*args)
        else:
            print("This should never happen", self._dim)
    #---------------------------------------------------------------------------
    def setval1d(self,i,x):
        self._column.buffer[i-1] = x
            
    #---------------------------------------------------------------------------
    def setval2d(self,i,j,x):
        try:
            self._plane2d.buffer[i-1][j-1] = x
        except:
            print("Pb in setval 2D " , i,j ,x)
    #---------------------------------------------------------------------------
    def setval3d(self,i,j,k,x):
        try:
            self._image.buffer[i-1][j-1][k-1] = x
        except:
            print("Pb in setval 3D " , i,j ,k,x)
    #---------------------------------------------------------------------------
    def val1d(self,i):
        return self._column.buffer[i-1]
    #---------------------------------------------------------------------------
    def val2d(self,i,j):
        return self._plane2d.buffer[i-1][j-1]
    #---------------------------------------------------------------------------
    def val3d(self,i,j,k,x):
        
        return self._image.buffer[i-1][j-1][k-1]
    #---------------------------------------------------------------------------
    def specw(self,*args):
        """
        Permits to enter the value for the spectral width of the current
        data-set. One parameter will be needed for each dimension of the
        data-set.

        When reading a file the spectral width is set to 2000 * 3.1416 if no
        parameter block is available.

        The value for spectral width are changed by EXTRACT


        see also : offset extract
        """
        if self._dim  == 1 :
            self.specw1d(*args)
        elif self._dim == 2 :
            self.specw2d(*args)
        elif self._dim == 3 :
            self.specw3d(*args)
        else:
            print("This should never happen", self._dim)
    #---------------------------------------------------------------------------
    def specw1d(self,x):
        self._column.axis1.specwidth = x            
    #---------------------------------------------------------------------------
    def specw2d(self,x,y):
        self._plane2d.axis1.specwidth = x
        self._plane2d.axis2.specwidth = y
    #---------------------------------------------------------------------------
    def specw3d(self,x,y,z):
        self._image.axis1.specwidth = x
        self._image.axis2.specwidth = y
        self._image.axis3.specwidth  = z
    #---------------------------------------------------------------------------
    def write(self,file_name):
        """
        write( file_name )

        Writes the current data set to a file in standard format. 
        same as writec

        see also : read
        """
        F = gf.GifaFile(file_name,"w")

        F.set_data(self._current)

        F.save()
        F.close()
    writec = write
    #---------------------------------------------------------------------------
    def set_task(self, task):
        print(task)
    #---------------------------------------------------------------------------
    def window(self):
        """
        window( {axis}, x, y)

        Define the window (with the starting point and the ending
        point) on which data is actually used for the iteration. Data
        outside this window(displayed as 0 during the interactive input) are
        just ignored for the processing. Window can be entered several time,
        the result being cumulative.

        see also : window_reset window_mode put apply
        """
        
    #---------------------------------------------------------------------------
    def window_reset(self):
        """
        window_reset( {axis})

        Resets the window to 1.0

        see also : window window_mode
        """
        
    #---------------------------------------------------------------------------
    def zero(self):
        self._current.buffer = np.zeros_like(self._current.buffer)
    #---------------------------------------------------------------------------
    def phase(self,ph0,ph1,axis = 1):
       for ax in self._test_axis(axis):
            self._current.phase(ph0,ph1,ax)
    #---------------------------------------------------------------------------
    def zoom(self,*args):
        if self._dim == 1:
            self._current.zoom(1,*args)
        elif self._dim == 2:
            self._current.zoom(2,*args)
        elif self._dim == 3:
            self._current.zoom(3,*args)
        return self
    #---------------------------------------------------------------------------
    def get_debug(self):
        return self.debug
    def get_si1_1d(self):
        return self._column.size1
    def get_si1_2d(self):
        return self._plane2d.size1
    def get_si2_2d(self):
        return self._plane2d.size2
    def get_si1_3d(self):
        return self._image.size1
    def get_si2_3d(self):
        return self._image.size2
    def get_si3_3d(self):
        return self._image.size3
    #---------------------------------------------------------------------------
    def get_freq(self):
        return self._current.frequency
    def get_freq_1d(self):
        return self._column.axis1.frequency
    def get_freq_1_2d(self):
        return self._plane2d.axis1.frequency
    def get_freq_2_2d(self):
        return self._plane2d.axis2.frequency
    def get_freq_1_3d(self):
        return self._image.axis1.frequency
    def get_freq_2_3d(self):
        return self._image.axis2.frequency
    def get_freq_3_3d(self):
        return self._image.axis3.frequency
    #---------------------------------------------------------------------------
    def get_itype_1d(self):
        return self._column.axis1.itype
    def get_itype_2d(self):
        return 2*self._plane2d.axis1.itype + self._plane2d.axis2.itype
    def get_itype_3d(self):
        return 4*self._image.axis1.itype + 2*self._image.axis2.itype + self._image.axis3.itype
    #---------------------------------------------------------------------------
    def get_shift(self):
        return self._shift
    def get_noise(self):
        return self._noise
    def get_offset_1d(self):
        return self._column.axis1.offset
    def get_offset_1_2d(self):
        return self._plane2d.axis1.offset
    def get_offset_1_3d(self):
        return self._image.axis1.offset
    def get_offset_2_2d(self):
        return self._plane2d.axis2.offset
    def get_offset_2_3d(self):
        return self._image.axis2.offset
    def get_offset_3_3d(self):
        return self._image.axis3.offset
    #---------------------------------------------------------------------------
    def get_specw_1d(self):
        return self._column.axis1.specwidth
    def get_specw_1_2d(self):
        return self._plane2d.axis1.specwidth
    def get_specw_2_2d(self):
        return self._plane2d.axis2.specwidth
    def get_specw_1_3d(self):
        return self._image.axis1.specwidth
    def get_specw_2_3d(self):
        return self._image.axis2.specwidth
    def get_specw_3_3d(self):
        return self._image.axis3.specwidth
    #---------------------------------------------------------------------------
    def get_dmin(self):
        return self.diffaxis.dmin
    def get_dmax(self):
        return self.diffaxis.dmax
    def get_dfactor(self):
        return self.diffaxis.dfactor
    #---------------------------------------------------------------------------
    def get_dim(self):
        return self._dim
    def get_version(self):
        return "%s" %(npkd.__version__)
    def get_npk1d(self):
        return len(self.peaks)
    def get_npk2d(self):
        return len(self.peaks2d[0])
    def get_npk3d(self):
        return len(self.peaks)

    #---------------------------------------------------------------------------
    def get_ph0(self):
        return self._last_ph0
    def get_ph1(self):
        return self._last_ph1
    def get_row(self):
        return self._last_row
    def get_col(self):
        return self._last_col
    #---------------------------------------------------------------------------
    def get_si_tab(self):
        return self._tab.size1
    #---------------------------------------------------------------------------
    def itoh(self,index,dim,axis):
        if dim == 1:
            return self._column.itoh(axis,index)
        elif dim == 2:
            return self._plane2d.itoh(axis,index)
        elif dim ==3:
            return self._image.itoh(axis,index)
        else:
            print("problem in itoh")
    def itop(self,index,dim,axis):
        if dim == 1:
            return self._column.itop(axis,index)
        elif dim == 2:
            return self._plane2d.itop(axis,index)
        elif dim == 3:
            return self._image.itop(axis,index)
        else:
            print("problem in itop")
    def htop(self,index,dim,axis):
        if dim == 1:
            return self._column.htop(axis,index)
        elif dim == 2:
            return self._plane2d.htop(axis,index)
        elif dim ==3:
            return self._image.htop(axis,index)
        else:
            print("problem in htop")
    def htoi(self,index,dim,axis):
        if dim == 1:
            return self._column.htoi(axis,index)
        elif dim == 2:
            return self._plane2d.htoi(axis,index)
        elif dim == 3:
            return self._image.htoi(axis,index)
        else:
            print("problem in htoi")
    def ptoh(self,index,dim,axis):
        if dim == 1:
            return self._column.ptoh(axis,index)
        elif dim == 2:
            return self._plane2d.ptoh(axis,index)
        elif dim ==3:
            return self._image.ptoh(axis,index)
        else:
            print("problem in ptoh")
    def ptoi(self,index,dim,axis):
        if dim ==1:
            return self._column.ptoi(axis,index)
        elif dim == 2:
            return self._plane2d.ptoi(axis,index)
        elif dim == 3:
             return self._image.ptoi(axis,index)
        else:
            print("problem in ptoi")
    #---------------------------------------------------------------------------
    def bruker_corr(self,):
        self._current.bruker_corr()
    def lb(self,value):
        self.lb = value
    def em(self,axis=0,lb=1.0):
        self._current.apod_em(axis,lb)
        self.lb = lb
    def tm(self,tm1,tm2,axis=0):
        print("TM still to do")
        self._current.apod_tm(axis,tm1,tm2)
        self.lb = lb
    def sqsin(self,maxi,axis=1):

        for ax in self._test_axis(axis):
            self._current.apod_sq_sin(ax,maxi)
    def sin(self, maxi, axis=1):
        # NPK_v1 has this format : sin(0.0, "f12")
        # testing parameters to map NPK_v1
        for ax in self._test_axis(axis):
            self._current.apod_sin(ax,maxi)
           
    join = writec
    #---------------------------------------------------------------------------
    def put(self,parameter, n = 0):
        """
        put(parameter)
        put(parameter, n) 
        
        Moves the content of the current buffer to an other buffer
        With parameter equal to:
        
        xx* DATA
            load the data to be used for MaxEnt processing or
                 as a off-hand place for processing
        
        in 1D only
        FILTER   load the filter used for Deconvolution.  If NCHANNEL is 
          greater than 1, then which channel you want to put.  eg.  PUT FILTER 
          3.  PUT FILTER 0  will consider the current data set as the 
          multichannel filter, and will load the whole filter. Useful when 
          associated with GET FILTER to store filters as files.
        WINDOW   load the window to be used for MaxEnt processing
        TAB      load the TAB buffer, used for tabulated fit.
        
        in 2D only
        ROW n   load the 1D buffer in the ROW n
        COL n   load the 1D buffer in the COL n
        
        in 3D only
        PLANE Fx n   load the 2D buffer in the plane Fx n
         
        see also : GET SHOW APPLY
        
        """
        if parameter == "DATA" or parameter == "data" :
            self._datab = self._current.copy()
        elif self._dim ==1:
            if parameter == "WINDOW":
                self._window = self._current.copy()
            elif parameter == "TAB":
                self._tab = self._current.copy()
            else:
                print("********************* Nothing has been done yet")
        else:
            print("********************* Nothing has been done yet")
                 
    #---------------------------------------------------------------------------
    def get(self, buffer_name):
        """
        get(buffer_name)

        if parameter == "DATA":
                self._datab = self._current.copy()
            

        
        Moves the content of another buffer, back to the current buffer with 
        buffer_name equal to:
        "data":      get the content of the data buffer
        "linefit":   get the simulated spectrum obtained form the current peak table
        "window":    get actual window used to compute the chisquare
        "filter":    get filter used for deconvolution
        "residue":   get residue of the spectrum after a maxent run
        "tab":        get the tab buffer used for tabulated fit
        see also : put apply
        """    
        if buffer_name == "DATA" or buffer_name == "data":
            self._current = self._datab.copy()
        elif buffer_name == "WINDOW":
            self._current = self._window.copy()
        elif buffer_name == "TAB":
            self._current = self._tab.copy()
        elif buffer_name == "FILTER":
            self._current = self._filter.copy()
    #---------------------------------------------------------------------------
    def exchdata(self):
        """
        
        Exchange the contents of the DATA buffer with the current data-set. 
        
        
        see also : adddata multdata maxdata mindata add mult put  
        
        
        """
        if self._current.dim != self._datab.dim:
            raise Exception("wrong buffer dim: %s"%str(self._current.dim()))
        if self._current.dim == 1:
            if self.get_si1_1d() != self._datab.size1 or self.get_itype_1d() != self._datab.itype:  
                raise Exception("wrong buffer size: %s %s"%(str(self.get_si1_1d()),str(self._datab.size1)))
        elif self._current.dim == 2:
            if self.get_si1_2d() != self._datab.size1 or self.get_si2_2d() != self._datab.size2 or self.get_itype_2d() != self._datab.itype:
                raise Exception("wrong buffer size 2D : %s"%str(self.get_si1_2d()))
        elif self._current.dim == 3:
            if self.get_si1_3d() != self._datab.size1 or self.get_si2_3d() != self._datab.size2 or self.get_si3_3d() != self._datab.size3 or self.get_itype_3d() != self._datab.itype:
                raise Exception("wrong buffer size 3D : %s"%str(self.get_si1_3d()))
        c = npkd.NPKData(buffer = self._current.buffer.copy())
        self._current = self._datab.copy()
        self._datab = c
        return self
    #---------------------------------------------------------------------------
    def escale(self, value = 1.0):
        """
        The Entropy expression during Maximum Entropy run is computed as follow :
        
            A = Escale * Sum(F(i))
            P(i) = F(i)/A
            S = -Sum( log(P(i)) * P(i) )
        
        Escale should be set to 1.0 for normal operation
        
        see also : maxent
        
        """
        if value != 0 :
            self.escale = value
        else:
            self.escale = 1
    #---------------------------------------------------------------------------
    def get_Kore_1D(self):
        "return a working copy of the 1D Kore internal buffer"
        return self._column.copy()
    #---------------------------------------------------------------------------
    def get_Kore_2D(self):
        "return a working copy of the 2D Kore internal buffer"
        return self._plane2d.copy()
    #---------------------------------------------------------------------------
    def get_Kore_3D(self):
        "return a working copy of the 3D Kore internal buffer"
        return self._image.copy()
    #---------------------------------------------------------------------------
    def set_Kore_1D(self, npkdata):
        "uses npkdata as the 1D Kore buffer"
        if npkdata.dim != 1:
            NPKError("SHould be a 1D", data=npkdata)
        self._column = npkdata
    #---------------------------------------------------------------------------
    def set_Kore_2D(self, npkdata):
        "uses npkdata as the 1D Kore buffer"
        if npkdata.dim != 2:
            NPKError("SHould be a 2D", data=npkdata)
        self._plane2d = npkdata
    #---------------------------------------------------------------------------
    def set_Kore_3D(self, npkdata):
        "uses npkdata as the 3D Kore buffer"
        if npkdata.dim != 3:
            NPKError("SHould be a 3D", data=npkdata)
        self._image = npkdata
#----------------------------------------------------------------------------
def compatibility(context):
    """
    inject Kore definition into context given by the caller
    """
    global kore
    for i in dir(kore):
        f = getattr(kore,i)
        if callable(f):
            if not i.startswith("_"):
                context[i] = f
                context["com_"+i] = f

global kore
kore = Kore()
compatibility(globals())

    