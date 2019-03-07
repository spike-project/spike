#!/usr/bin/env python 
# encoding: utf-8

"""set of function for Peak fitter

Very First functionnal - Not finished !


reqruies the Peaks plugin installed 

July 2016 M-A Delsuc
"""

from __future__ import print_function
import numpy as np
import unittest
import scipy
from scipy.optimize import minimize, curve_fit

import spike
from spike import NPKError
from spike.NPKData import NPKData_plugin, NPKData, flatten, parsezoom
from spike.util.counter import timeit
# This counter is used to count function evaluation
count = 0

# epsilon
eps = 1E-7

def Lor(Amp, Pos, Width, x):
    """
    One Lorentzian
    Param contains in sequence Amp_i, Pos_i, Width_i
    """
    L = Amp/(1+(2*(x-Pos)/(Width+eps))**2)
    return L
def cSpec(x, *Param):
    return Spec(Param, x)
def Spec(Param, x):
    """
    x is the spectral coordinates
    Param contains in sequence Amp_i, Pos_i, Width_i
    all coordinates are in index
    """
    global count
    count += 1
    y = np.zeros_like(x)
    for i in range(0,len(Param),3):
        a = Param[i]
        p = Param[i+1]
        w = Param[i+2]
        y += Lor(a, p , w, x)
    return y
def residu(Params, x, y):
    """
    The residue function, returns a vector Ycalc(Params) - y_experimental
    can be used by leastsq
    """
    Yc = Spec(Params, x)
    res = Yc-y
    return res
def tofit(Params, x, y):
    """
    calls the residue function, and return sum( residue**2 )
    can be used by minimize
    """
    res = residu(Params, x, y)
    val = (res**2).sum()
    return val/1E4

# the derivatives of the functions above
def dLor(Amp, Pos, Width, x):
    #L = Amp/(1+(2*(x-Pos)/Width)**2)
    a = 2* (x-Pos)/ Width
    b = 1/(1+a**2)
    c = Amp * b
    c = c * b * a * 2 / (Width+eps)
    return (b, 2*c, a*c)  # d/dA  d/dP  d/dW
def dSpec(Param, x, y=None):
    """
    Param contains in sequence Amp_i, Pos_i, Width_i
    """
    global count
    count += 1
    dy = np.zeros((len(x), len(Param)))
    for i in range(0,len(Param),3):
        a = Param[i]
        p = Param[i+1]
        w = Param[i+2]
        dA, dP, dW =  dLor(a, p , w, x)
        dy[:,i] = dA
        dy[:,i+1] = dP
        dy[:,i+2] = dW
    return dy
def cdSpec(x, *Param):
    return dSpec(Param, x)
def dToFit(Param, x, y):
    dS = 2* dSpec(Param,x)
    res = residu(Param, x, y)
    dT = np.dot(dS.T, res)
    return dT/1E4

def simulate(npkd, zoom=None):
    """
    Simulate the 1D npkd data-set using the content of the Peak List
    replace the current data-set
    """
    # 3 parameters per peaks : Amp, Pos, Width
    z1, z2 = parsezoom(npkd, zoom)
    PP = []
    for i,pk in enumerate(npkd.peaks):
        if pk.pos>=z1 and pk.pos<=z2:
            PP.append(pk.intens)     # Amp
            PP.append(pk.pos)   # Pos
            PP.append(max(1.0,pk.width))    # Width - mini is one pixel !
#    print (PP)
    x = np.arange(1.0*npkd.cpxsize1)
    npkd.set_buffer( Spec(PP, x) )
    return npkd

def fit(npkd, zoom=None):
    """
    fit the 1D npkd data-set for Lorentzian line-shape
    current peak list is used as an initial values for fitting
    Only peaks within the zoom windows are fitted
    
    fitting is contraint from the initial values
        - intensity will not allowed to change by more than x0.5 to x2
        - positions by more than 5 points
        - width by more than x5
        (constraints work only for scipy version >= 0.17 )
    It may help to use centroid() to pre-optimize the peak list before calling fit(), or calling fit() twice (slower)
    
    """
    # 3 parameters per peaks : Amp, Pos, Width
    z1, z2 = parsezoom(npkd, zoom)
    PP = []
    minbound = []
    maxbound = []
    # load initial values and constraints from peak list
    for i,pk in enumerate(npkd.peaks):
        if pk.pos>=z1 and pk.pos<=z2:
            PP.append(pk.intens)     # Amp
            minbound.append(0.5*pk.intens)
            maxbound.append(2*pk.intens)
            PP.append(pk.pos)   # Pos
            minbound.append(pk.pos-5)
            maxbound.append(pk.pos+5)
            PP.append(max(1.0,pk.width))    # Width - mini is one pixel !
            minbound.append(1E-3)
            maxbound.append(max(5.0,5*pk.width))
#    print (PP)
    x = np.arange(1.0*npkd.size1)[z1:z2]
    Y = npkd.get_buffer()[z1:z2]
    Y = Y.real
#    kwargs={"jac":cdSpec}
    if scipy.__version__ > '0.17.0':
        PP1 = curve_fit(cSpec, x, Y, PP, bounds=(minbound,maxbound), method='dogbox')
    else:
        PP1 = curve_fit(cSpec, x, Y, PP)
    results = PP1[0]
    errors = np.sqrt(np.diag(PP1[1]))
    chi2 = tofit(results,x,Y)   # computes error and store it
    npkd.peaks.chi2 = chi2
    # copy back
    for i,pk in enumerate(npkd.peaks):
        if pk.pos>=z1 and pk.pos<=z2:
            pk.intens = results[3*i]
            pk.pos = results[3*i+1]
            pk.width = results[3*i+2]
            pk.intens_err = errors[3*i]
            pk.pos_err = errors[3*i+1]
            pk.width_err = errors[3*i+2]
    return npkd

def display_fit(npkd, **kw):
    """
    displays the result of the fit
    accept the same arguments than display()
    """
    d = npkd.copy()
    d.peaks = npkd.peaks
    try:
        z = kw['zoom']
    except:
        z = None
    simulate(d, zoom=z)
    d.display(**kw)
    return npkd

class FitTests(unittest.TestCase):
    """
    Test for fitter, assumes Peaks plugin is loaded
    """
    def test_fit1d(self):
        # create 1D spectrum
        t = np.linspace(0,10,1000)
        y = np.zeros_like(t)
        A = (100,100,100)
        W = (100, 110, 115)
        TAU = (0.3, 1, 3)
        for a,w,tau in zip(A,W, TAU):
            y += a*np.cos(w*t)*np.exp(-t*tau)
        Y = np.fft.rfft(y).real
        Y -= Y[0]
        # load and peak pick
        d=spike.NPKData.NPKData(buffer=Y)
        d.pp(threshold=1000)
        # check
        self.assertEqual(list(d.peaks.pos) , [159.0, 175.0, 183.0])
        d.fit()
        if scipy.__version__ > '0.17.0':
            # first fit is not full because of constraints on widthes (third peak)
            self.assertAlmostEqual(d.peaks.chi2, 121.72613405, places=2)
        d.fit()
        self.assertAlmostEqual(d.peaks.chi2, 15.0445981291, places=2)    # second is complete
        # other possibility is centroid
        d.pp(threshold=1000)
        d.centroid()
        d.fit(zoom=(140,200))
        self.assertAlmostEqual(d.peaks.chi2, 12.4304236435, places=1)    # lower because of zoom.
        self.assertAlmostEqual( sum(list(d.peaks.pos)), 517.74817237246634, places=2)

NPKData_plugin("simulate", simulate)
NPKData_plugin("fit", fit)
NPKData_plugin("display_fit", display_fit)