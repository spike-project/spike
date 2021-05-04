#!/usr/bin/env python 
# encoding: utf-8

"""
Phase correction for MS

includes 2nd order correction

"""

from __future__ import print_function
from spike.NPKData import NPKData_plugin
from spike.NPKError import NPKError
import numpy as np
import math as m

def phase(self, ph0, ph1, ph2=0.0, pivot=0.0, axis=0):
    """
    apply a phase correction along given axis
    PH0 is in degrees
    PH1 PH2 corrections are in nb of turns over the whole spectral width, 
        warning, this is different from usual NMR phases
    pivot is the center of the 1st & 2nd order term goes 0.0 ... 1.0
        default value is 0 (left/slow frequency side)
        warning, this is different from usual NMR phases
    for a N complex spectrum, correction on ith point is 
        ph = ph0 + ph1*(i/N - pv) + ph2*(i/N - pv)^2 
    """
    todo = self.test_axis(axis)
    # compute shape parameters
    it = self.axes(todo).itype
    sw = self.axes(todo).specwidth
    if it == 0:
        raise NPKError(msg="no phase correction on real data-set", data=self)
    size = self.axes(todo).size//2       # compute correction in e
    coef2 = self.axes(todo).specwidth//2
    p = pivot
    le = m.radians(float(ph0)) + \
        2*m.pi*ph1*np.linspace(-p, 1-p, size) + \
        2*m.pi*ph2*(np.linspace(-p, 1-p, size))**2
    e = np.cos(le) + 1J*np.sin(le)
#            e = np.exp( 1J* e)
    # then apply
    self.axes(todo).P0 = ph0
    self.axes(todo).P1 = ph1
    self.axes(todo).P2 = ph2
    self.axes(todo).PV = pivot
    return self.mult_by_vector(axis, e, mode="complex")

def movepivot(p0, p1, p2, pvbef, pvaft):
    """computes a new set of parameters whil moving pivot from pvbef to pvaft"""
# p (pbef) goes from 0 to 1 so, for a buffer of size, pivot position is p*size
# in position x, phase correction is prop to (x/size - p) 
# ph = P0 + (x/size - p) P1 +  (x/size - p)² P2 
# ph = P0 + x/size P1 - p P1  + (x/size)² P2 +  p² P2 - 2 x p/size P2
# ph = P0 - p P1 +  p² P2   - 2 x p/size P2 + x/size P1     + (x/size)² P2
# and with a new pivot paft:
# ph = P0' - paft P1' +  paft² P2'   - 2 x paft/size P2' + x/size P1'     + (x/size)² P2'
#      P0' - paft P1' +  paft² P2'   + (x/size)(- 2 paft P2' +  P1')     + (x/size)² P2'
# for the phase to be equal with paft we need that cst, x and x² terms to be equal, which means
# so
# => P2' = P2
# - 2 paft P2' + P1' = - 2 p P2 + P1
# => P1'  =  P1 - 2 p P2 + 2 paft P2
#         =  P1 + 2(paft-p)P2
# P0' - paft P1' +  paft² P2'  = P0 - p P1 +  p² P2 
# => P0' =  P0 - p P1 +  p² P2 + paft P1' - paft² P2'
#        =  P0 - p P1 +  (p²-paft²) P2 + paft P1'    
    p2p = p2
    p1p = p1 + 2*(pvaft-pvbef)*p2
    p0p = 360*(p0/360 - pvbef*p1 + (pvbef**2 - pvaft**2)*p2 + pvaft*p1p)
    # then bring ph0 in [-180, 180]
    while p0p > 180:
        p0p -= 360
    while p0p < -180:
        p0p += 360
    
    return p0p, p1p, p2p, pvaft
# NPKData_plugin("phase", phase)  # supercede the regular phase correction code - not compatible with Test suite !
NPKData_plugin("phaseMS", phase)
