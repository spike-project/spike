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
import File.HDF5File
from spike import NPKData
from spike import FTMS
from .NPKError import NPKError

FREQ0 = 1E6
REF_FREQ = 419620.0     # Ad Hoc values for   REF_FREQ / REF_MASS
REF_MASS = 344.0974

class FTICRAxis(FTMS.FTMSAxis):
    """
    hold information for one FT-ICR axis
    used internally
    """
    def __init__(self, size=1024, specwidth=FREQ0, itype=0, currentunit="points", ref_mass=REF_MASS, ref_freq=REF_FREQ, highmass=10000.0, left_point = 0.0):
        """
        all parameters from Axis, plus
        specwidth   highest frequency,
        ref_mass    m/z used for calibration
        ref_freq    frequency of the calibration m/z
        highmass    highest m/z of interest - usually defined by the excitation pulse low frequency limit
        left_point  coordinates of first data point; usually 0.0 after Fourier Transform; may be different after extraction
        
        possible values for unit are "points" "m/z"
        
        conversion methods work on numpy arrays as well
        """
        super(FTICRAxis, self).__init__(size=size, specwidth=specwidth, itype=itype, currentunit=currentunit, ref_mass=ref_mass, ref_freq=ref_freq, highmass=highmass, left_point = left_point)
        self.FTICR = "FTICR"
        self.calibA = self.ref_mass*self.ref_freq
        self.calibB = 0.0
        self.calibC = 0.0
        
        self.attributes.insert(0,"FTICR") # updates storable attributes
        self.attributes.insert(0,"calibA") # ML1
        self.attributes.insert(0,"calibB") # ML2
        self.attributes.insert(0,"calibC") # ML3
        
    #-------------------------------------------------------------------------------
    def report(self):
        "high level reporting"
        
        if self.itype == 0: # Real
            return "FT-ICR axis at %f kHz,  %d real points,  from mz = %8.3f   to m/z = %8.3f  R max (M=400) = %.0f"%  \
            (self.specwidth/1000, self.size, self.highmass, self.lowmass, 400.0/self.deltamz(400.))
        else: # Complex
            return "FT-ICR axis at %f kHz,  %d complex pairs,  from mz = %8.3f   to m/z = %8.3f  R max (M=400) = %.0f"%  \
            (self.specwidth/1000, self.size/2, self.highmass, self.lowmass, 400.0/self.deltamz(400.))


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
            delta = self.calibA**2 + 4*self.calibC*(self.calibB + h)
            m = self.calibA + np.sqrt(delta) / (2 * (self.calibB + h))
        return m
    def mztoh(self, value):
        """
        return Hz value (h) from  m/z (mz) 
        """
        m = np.maximum(value,0.1)             # protect from divide by 0
        return self.calibA/m + self.calibC/(m**2)  - self.calibB

#-------------------------------------------------------------------------------
def fticr_mass_axis(length, spectral_width, ref_mass, ref_freq):
    """
    returns an array which will calibrate a FT-ICR experiment
    length : number of points in the axis
    spectral_width : of the ICR measure
    ref_mass : value of the m/z reference
    ref_freq =: frequence at which is is observed.
    """
    return ref_mass*ref_freq/(spectral_width * (1+np.arange(length))/length )

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
                H = File.HDF5File.HDF5File(name,"r")
                H.load(mode=mode, group=group)      # load into memory by default !
                super(FTICRData, self).__init__(buffer=H.data.buffer, debug=debug)
                NPKData.copyaxes(H.data, self)  # and deep copy all axes from file
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
        F = FTICRAxis(size=1000, specwidth=1667000, itype=0, currentunit="points", ref_mass=344.0974, ref_freq=419620.0, highmass=2000.0)
        self.assertAlmostEqual(F.deltamz(F.itomz(1023))/F.itomz(1023), 1.0/1023)    # delta_m / m is 1/size at highest mass - before extraction
        self.assertAlmostEqual(123*F.itomz(123), 321*F.itomz(321))    # verify that f*m/z is constant !
        self.assertAlmostEqual(123*F.mztoi(123), 321*F.mztoi(321))    # verify that f*m/z is constant !
        self.assertAlmostEqual(F.itoh(0), 0)
        self.assertAlmostEqual(F.itoh(F.size-1), F.specwidth)   # last point is size-1 !!!
        for twice in range(2):
            F = FTICRAxis(size=1000, specwidth=1667000, itype=0, currentunit="points", ref_mass=344.0974, ref_freq=419620.0, highmass=2000.0)
            if twice == 1:  # second time tests 2nd order calibration
                print ("TWICE") 
                F.calibB = 55.0
                F.calibC = 0.0
            for i in (1,2):
                print(F.report())
                print(F._report())   # low level report
                for x in (301.0, 1.0, F.size-20.0, F.size-1.0):    # last point is size-1 !!!
                    print("point at index %d is at freq %f, m/z %f"%(x, F.itoh(x), F.itomz(x)))
                    self.assertAlmostEqual(F.itoh(F.htoi(x)), x)
                    self.assertAlmostEqual(F.mztoi(F.itomz(x)), x)
                F.extract((300,F.size-20))
    def test_axis(self):
        "testing FTICRAxis object"
        self.announce()
        A = FTICRAxis(specwidth=1667000)
        print(fticr_mass_axis(A.size, A.specwidth, A.ref_mass, A.ref_freq))
        print(A.itomz(np.arange(1,A.size)))
        print(A.mass_axis())
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

