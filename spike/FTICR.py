#!/usr/bin/env python 
# encoding: utf-8

"""
This file implements all the tools for handling FT-ICR data-sets

It allows to work with 1D and 2D 

To use it :

import FTICR
d = FTICR.FTICRData(...)    # There are several possible initialisation : empty, from file
play with d

d will allows all NPKData methods, plus a few specific ones.

alternatively, use an importer :
from File.(Importer_name) import Import_1D
d = Import_1D("filename)")


Created by Marc-André on 2014-08
Copyright (c) 2014 IGBMC. All rights reserved.
"""

from __future__ import print_function
import math
import unittest
import numpy as np
from . import NPKData
from .File import HDF5File
from . import FTMS
from .NPKError import NPKError

FREQ0 = 1E6
REF_FREQ = 419620.0     # Ad Hoc values for   REF_FREQ / REF_MASS
REF_MASS = 344.0974

class FTICRAxis(FTMS.FTMSAxis):
    """
    hold information for one FT-ICR axis
    used internally
    """
    def __init__(self, itype=0, currentunit="points", size=1024, specwidth=1E6,  offsetfreq=0.0, left_point = 0.0, highmass=10000.0, calibA=1E8, calibB=0.0, calibC=0.0, lowfreq=1E4, highfreq=1E6 ):
        """
        all parameters from Axis, plus
        specwidth   highest frequency,
        offsetfreq      carrier frequency in heterodyn or lowest frequency if acquisition does not contains 0.0,

        calibA, calibB, calibC : calibration constant, allowing 1 2 or 3 parameters calibration.
            set to zero if unused
            correspond to Bruker parameter ML1 ML2 ML3 for FTICR
            correspond to Thermo parameter  'Source Coeff1', 'Source Coeff2', 'Source Coeff3' for Orbitrap
        highmass    highest physical m/z of interest
        left_point  coordinates of first data point; usually 0.0 after Fourier Transform; may be different after extraction        
        currentunit default unit used for display and zoom,
            possible values for unit are "points" "m/z"

        lowfreq     lowest excitation pulse frequency
        highfreq     highest excitation pulse frequency
        
        """
        super(FTICRAxis, self).__init__(itype=itype, currentunit=currentunit, size=size,
            specwidth=specwidth, offsetfreq=offsetfreq, left_point=left_point, highmass=highmass,
            calibA=calibA, calibB=calibB, calibC=calibC)
        self.FTICR = "FTICR"

        self.lowfreq = lowfreq
        self.highfreq = highfreq

        self.attributes.insert(0,"FTICR") # updates storable attributes
        self.attributes.insert(0,"highfreq")
        self.attributes.insert(0,"lowfreq")
        
    #-------------------------------------------------------------------------------
    def report(self):
        "high level reporting"
        if self.itype == 0: # Real
            return "FT-ICR report axis at %f kHz,  %d real points,  from physical mz = %8.3f   to m/z = %8.3f  R max (M=400) = %.0f"%  \
            (self.specwidth/1000, self.size, self.lowmass, self.highmass, 400.0/self.deltamz(400.))
        else: # Complex
            return "FT-ICR report axis at %f kHz,  %d complex pairs,  from physical mz = %8.3f   to m/z = %8.3f  R max (M=400) = %.0f"%  \
            (self.specwidth/1000, self.size/2, self.lowmass, self.highmass, 400.0/self.deltamz(400.))

    #-------------------------------------------------------------------------------
# these are Bruker equations:
#
#    mzaxis = ml1/(ml2+faxis)
#    faxis = ml1/mzaxis - ml2
# or
# f = (ML1 / m) + (ML3 / m^2) - ML2
# m^2 F = -m^2*ML2 + m*ML1 +ML3
# m^2*(F+ML2) - m*ML1 - ML3 = 0
# delta = ML1^2 + 4*(ML2+F)ML3
# m = ML1 + sqrt(delta) / (2 * (ML2 + f))

    def htomz(self, value):
        """
        return m/z (mz) from hertz value (h)
        """
        h = np.maximum(value,0.1)       # protect from divide by 0
        if self.calibC == 0.0:
            m = self.calibA/(self.calibB + h)
        else:
#            delta = self.calibA**2 + 4*self.calibC*(self.calibB + h)
#            m = self.calibA + np.sqrt(delta) / (2 * (self.calibB + h))
            m = self.calibA/(2*(self.calibB + h)) + np.sqrt(self.calibA**2 + 4*self.calibB*self.calibC + 4*self.calibC*h)/(2*(self.calibB + h))
            #m = -1*(self.calibA+np.sqrt((self.calibA**2)-4*(self.calibB-h)*self.calibC))/(2*(self.calibB-h)) #WK Edit 13/11/2017 based on Solarix info from DPAK
        return m
    def mztoh(self, value):
        """
        return Hz value (h) from  m/z (mz) 
        """
        m = np.maximum(value,0.1)             # protect from divide by 0
        return self.calibA/m + self.calibC/(m**2)  - self.calibB  # MAD
        #return (self.calibC/(m**2)) + (self.calibA/m)  + self.calibB # WK

#-------------------------------------------------------------------------------

class FTICRData(FTMS.FTMSData):
    """
    subclass of FTMS.FTMSData, meant for handling FT-ICR data
    allows 1D and 2D data-sets
    
    """
    def __init__(self, dim=1, shape=None, mode="memory", group="resol1", buffer=None, name=None, debug=0):
        """
        dim : dimension of dataset
        shape : shape of the buffer (size1,size2)
        mode : memory : data-set is kept in-memory    /  onfile : data-set is read from file when needed
        buffer : if is not None; used as data
        name : if is not None, data is read from file
            group : when reading a hdf5 file, default group name used
        """
        self.axis1 = FTICRAxis()    # this creates an FTMSAxis so that pylint does not complain - will be overwritten
        if dim == 2:
            self.axis2 = FTICRAxis()
        if name:
            if name.endswith(".msh5"):  # try loading .msh5 file
                if debug>0: print("reading msh5")
                H = HDF5File.HDF5File(name,"r")
                H.load(mode=mode, group=group)      # load into memory by default !
                super(FTICRData, self).__init__(buffer=H.data.buffer, debug=debug)
                NPKData.copyaxes(H.data, self)  # and deep copy all axes from file
                try:
                    self.params = H.retrieve_object('params')
                except:
                    print('params block is missing in this file')
                self.name = name
                self.hdf5file = H
            else:
                raise Exception("Filename should have a .msh5 extension")
        else:
            if debug>0: print("calling super")
            super(FTICRData, self).__init__(dim=dim, shape=shape, buffer=buffer, name=name, debug=debug)
            for i in range(self.dim):
                axis = self.axes(i+1)
                setattr(self, "axis%d"%(i+1), FTICRAxis(size=axis.size, specwidth=axis.specwidth, itype=axis.itype) )
        if debug>1: print(self.report())


#-------------------------------------------------------------------------------
class FTICR_Tests(unittest.TestCase):
    def setUp(self):
        self.verbose = 1    # verbose > 0 switches messages on
    def announce(self):
        if self.verbose >0:
            print("\n========",self.shortDescription(),'===============')
    def test_atob(self):
        "testing unit conversion functions"
        self.announce()
        F = FTICRAxis(size=1000, specwidth=1667000, itype=0, currentunit="points", calibA=344.0974*419620.0, highmass=2000.0, offsetfreq=10000.0)
#        self.assertAlmostEqual(F.deltamz(F.itomz(1023))/F.itomz(1023), 1.0/1023)    # delta_m / m is 1/size at lowest mass - before extraction
#        self.assertAlmostEqual((F.lowmass/F.deltamz(F.lowmass)) /( F.size*(F.specwidth+F.offsetfreq)/F.specwidth-1), 1.0, 5)
        self.assertAlmostEqual(123*F.htomz(123), 321*F.htomz(321))    # verify that f*m/z is constant !
        self.assertAlmostEqual(123*F.mztoh(123), 321*F.mztoh(321))    # verify that f*m/z is constant !
        self.assertAlmostEqual(F.itoh(0), F.offsetfreq)
        self.assertAlmostEqual(F.itoh(F.size-1), F.specwidth+F.offsetfreq)   # last point is size-1 !!!
        for twice in range(2):
            m0 = 555.0
            f0 = 260162.434213
            F = FTICRAxis(size=1000, specwidth=1667000, itype=0, currentunit="points", calibA=344.0974*419620.0, highmass=2000.0)
            if twice == 1:  # second time tests 2nd order calibration
                print ("TWICE") 
                F.calibB = 55.0
                F.calibC = 0.0
                f0 = 260107.434213
            for i in (1,2):
                print(F.report())
                print(F._report())   # low level report
                for x in (301.0, 1.0, F.size-20.0, F.size-1.0):    # last point is size-1 !!!
                    print("point at index %d is at freq %f, m/z %f"%(x, F.itoh(x), F.itomz(x)))
                    self.assertAlmostEqual(F.itoh(F.htoi(x)), x)
                    self.assertAlmostEqual(F.mztoi(F.itomz(x)), x)
                # print("point at m/z %f is at freq %f"%(m0, F.mztoh(m0)))
                self.assertAlmostEqual(F.mztoh(m0), f0, 5)
                F.extract((300,F.size-20))
        self.assertAlmostEqual(F.mztoh(m0), f0, 5)

    def test_axis(self):
        "testing FTICRAxis object"
        self.announce()
        A = FTICRAxis(specwidth=1667000)
        ax1 = A.itomz(np.arange(0.0,A.size))[1:]
        ax2 = A.mass_axis()[1:]
        print (ax1)
        print (ax2)
        self.assertAlmostEqual(ax1.sum(), ax2.sum())
        self.assertAlmostEqual(ax1.min(), 5.99880024e+01)
        self.assertAlmostEqual(ax1.max(), 61367.72645470906)
        #        6.13677265e+04   3.06838632e+04   2.04559088e+04 ...,   6.01055107e+01, 6.00466991e+01   5.99880024e+01
    def test_trim(self):
        """
        Test trimz 
        """
        A = FTICRData(buffer=np.zeros((500, 10000)))
        A.specwidth = 1667000
        A.ref_mass = 344.0974
        A.ref_freq = 419620.0
        A.highmass = 1000.0
        print(A.report())
        l1 = int(A.axis1.mztoi(A.axis1.highmass))
        l2 = int(A.axis2.mztoi(A.axis2.highmass))
        print("l values:",l1,l2)
        A.trimz()
        self.assertEqual(l1, A.axis1.left_point)
        self.assertEqual(l2, A.axis2.left_point)
        self.assertEqual(A.size1, 500-l1)
        self.assertEqual(A.size2, 10000-l2)
        print("2D trimz gain is : %d %%" % (100*(1-(A.size1*A.size2/(500.*10000)))))
    def test_saving_1D(self):
        """
        Testing how save_msh5 works on 1D spectrum
        """
        from .Tests import filename
        self.announce()
        A = FTICRData(buffer=np.zeros(10000))
        A.specwidth = 1667000
        A.ref_mass = 344.0974
        A.ref_freq = 419620.0
        A.highmass = 1000.0
        print(A.report())
        A.save_msh5(filename("1D_test.msh5"))
#-------------------------------------------------------------------------------
if __name__ == '__main__':
    unittest.main()

