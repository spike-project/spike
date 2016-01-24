#!/usr/bin/env python 
# encoding: utf-8

"""
This file implements all the tools for handling Orbitrap data-sets

To use it :

import Orbitrap
d = Orbitrap.OrbiData(...)    # There are several possible initialisation : empty, from file
play with d

d will allows all NPKData methods, plus a few specific ones.

alternatively, use an importer :
from File.(Importer_name) import Import_1D
d = Import_1D("filename)")


Created by Marc-Andre' on 2014-09
Copyright (c) 2014 IGBMC. All rights reserved.
"""

from __future__ import print_function
import math
import unittest
import numpy as np
from .File.HDF5File import HDF5File
from spike import NPKData
from spike import FTMS
from .NPKError import NPKError


FREQ0 = 1E7
REF_FREQ = 1887533.975611561     # Ad Hoc values for   REF_FREQ / REF_MASS
REF_MASS = 715.3122

class OrbiAxis(FTMS.FTMSAxis):
    """
    hold information for one Orbitrap axis
    used internally
    """    
    def __init__(self, itype=0, currentunit="points", size=1024, specwidth=1E6,  offsetfreq=0.0, left_point = 0.0, highmass=10000.0, calibA=0.0, calibB=1E14, calibC=0.0):
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
        
        conversion methods work on numpy arrays as well
        """
        super(OrbiAxis, self).__init__(itype=itype, currentunit=currentunit, size=size,
            specwidth=specwidth, offsetfreq=offsetfreq, left_point=left_point, highmass=highmass,
            calibA=calibA, calibB=calibB, calibC=calibC)
        self.Orbitrap = "Orbitrap"
        self.calibA = calibA
        self.calibB = calibB
        self.calibC = calibC
        self.attributes.insert(0,"Orbitrap") # updates storable attributes
        self.attributes.insert(0,"calibA") # updates storable attributes
        self.attributes.insert(0,"calibB") # updates storable attributes
        self.attributes.insert(0,"calibC") # updates storable attributes

    #-------------------------------------------------------------------------------
    def report(self):
        "high level reporting"
        
        if self.itype == 0: # Real
            return "Orbitrap axis at %f kHz,  %d real points,  from mz = %8.3f   to m/z = %8.3f  M/DeltaM (M=400) = %.0f"%  \
            (self.specwidth/1000, self.size, self.lowmass, self.highmass, 400.0/self.deltamz(400.))
        else: # Complex
            return "Orbitrap axis at %f kHz,  %d complex pairs,  from mz = %8.3f   to m/z = %8.3f  M/DeltaM (M=400) = %.0f"%  \
            (self.specwidth/1000, self.size/2, self.lowmass, self.highmass, 400.0/self.deltamz(400.))

    #-------------------------------------------------------------------------------
    def htomz(self, value):
        """
        return m/z (mz) from Hertz value (h)
        """
        # m/z = A + B/f^2 + C/f^4
        v = np.maximum(value,0.1)   # protect from divide by 0
        return self.calibA + self.calibB/(v**2) + self.calibC/(v**4)
 
    def mztoh(self, value):
        """
        return Hz value (h) from  m/z (mz) 
        """
        v = np.maximum(value,0.1)   # protect from divide by 0
        Delta = self.calibB**2 - 4*self.calibC*(self.calibA - v)
        f2 = (-self.calibB - np.sqrt(Delta)) / (2*(self.calibA - v))
        return np.sqrt(f2)

#-------------------------------------------------------------------------------
class OrbiData(FTMS.FTMSData):
    """
    subclass of FTMS.FTMSData, meant for handling Orbitrap data
    doc to be written ...
    """
#    print "in Orbitrap"
    def __init__(self, dim=1, shape=None, mode="memory", buffer=None, name=None, debug=0):
        self.axis1 = OrbiAxis()    # this creates an OrbiAxis so that pylint does not complain - will be overwritten
        if dim == 2:
            raise Exception("2D Orbitrap is not physcally defined (yet ?)")
        if name:
            if name.endswith(".msh5"):  # try loading .msh5 file
                if debug>0: print("reading msh5")
                H = HDF5File(name,"r")
                H.load(mode=mode)      # load into memory by default !
                super(OrbiData, self).__init__(buffer=H.data.buffer, debug=debug)
                NPKData.copyaxes(H.data, self)  # and deep copy all axes from file
                self.name = name
                self.hdf5file = H
            else:
                raise Exception("Filename should have a .msh5 extension")
        else:
            if debug>0: print("calling super")
            super(OrbiData, self).__init__(dim=dim, shape=shape, buffer = buffer, name = name, debug=debug)
            for i in range(self.dim):
                axis = self.axes(i+1)
                setattr(self, "axis%d"%(i+1), OrbiAxis(size=axis.size, itype=0) )
        if debug>1: print(self.report())


#-------------------------------------------------------------------------------
class Orbi_Tests(unittest.TestCase):
    def setUp(self):
        self.verbose = 1    # verbose > 0 switches messages on
    def announce(self):
        if self.verbose >0:
            print("\n========",self.shortDescription(),'===============')
    def test_atob(self):
        "testing unit conversion functions"
        self.announce()
        Oaxis = OrbiAxis(size = 1000, specwidth = 1667000, itype = 0, currentunit = "points", calibB=344.0974*(419620.0**2), highmass = 2000.0)
        self.assertAlmostEqual(Oaxis.itoh(0), 0)
        self.assertAlmostEqual(Oaxis.itoh(Oaxis.size-1), Oaxis.specwidth)   # last point is size-1 !!!
        self.assertAlmostEqual(Oaxis.itomz(1023)/Oaxis.deltamz(Oaxis.itomz(1023)), 1023./2, places = 2)    # delta_m / m is 1/size at highest mass - before extraction
        self.assertAlmostEqual(123**2*Oaxis.itomz(123), 321**2*Oaxis.itomz(321))    # verify that f*m/z is constant !
        self.assertAlmostEqual(123*Oaxis.mztoi(123)**2, 321*Oaxis.mztoi(321)**2)    # verify that f*m/z is constant !
        for i in (1,2):
            print(Oaxis.report())
            print(Oaxis._report())   # low level report
            for x in (1.0, 301.0, Oaxis.size-20.0, Oaxis.size-1.0):    # last point is size-1 !!!
                print("point at index %d is at freq %f, m/z %f"%(x, Oaxis.itoh(x), Oaxis.itomz(x)))
                self.assertAlmostEqual(Oaxis.mztoi(Oaxis.itomz(x)), x)
                self.assertAlmostEqual(Oaxis.itoh(Oaxis.htoi(x)), x)
            Oaxis.extract([300,Oaxis.size-20])
            
    def test_trim(self):
        """
        Test trimz 
        """
        O = OrbiData(buffer=np.zeros((500, 10000)))
        O.specwidth = 1667000
        O.ref_mass = 344.0974
        O.ref_freq = 419620.0
        O.highmass = 1000.0
        print(O.report())
        l1 = int(O.axis1.mztoi(O.axis1.highmass))
        l2 = int(O.axis2.mztoi(O.axis2.highmass))
        O.trimz()
        self.assertEqual(l1, O.axis1.left_point)
        self.assertEqual(l2, O.axis2.left_point)
        self.assertEqual(O.size1, 500-l1)
        self.assertEqual(O.size2, 10000-l2)
        print("2D trimz gain is : %d %%" % (100*(1-(O.size1*O.size2/(500.*10000)))))

#-------------------------------------------------------------------------------
if __name__ == '__main__':
    unittest.main()

