#!/usr/bin/env python 
# encoding: utf-8

"""set of function for the baseline correction

First version - Not finished !

improved July 2016
"""

from __future__ import print_function, division
import numpy as np
from scipy import interpolate
from scipy.optimize import leastsq

from spike import NPKError
from spike.NPKData import NPKData_plugin

import sys
if sys.version_info[0] < 3:
    pass
else:
    xrange = range
#-------------------------------------------------------------------------------
def get_ypoints(buff, xpoints, nsmooth=0):
    """
    from  buff and xpoints, returns ypoints = buff[xpoints]
    eventually smoothed by moving average over 2*nsmooth+1 positions
    """
    nsmooth = abs(nsmooth)
    xp = np.array(xpoints).astype(int)
    y = np.zeros(len(xpoints))
    for i in range(2*nsmooth+1):
        xi = np.minimum( np.maximum(xp-i, 0), len(buff)-1)  # shift and truncate
        y += buff[xi]
    y /= 2*nsmooth+1
    return y
#-------------------------------------------------------------------------------
def _spline_interpolate(buff, xpoints, kind = 3, nsmooth=0):
    """compute and returns a spline function 
        we are using splrep and splev instead of interp1d because interp1d needs to have 0 and last point
        it doesn't extend.
    """
    if len(xpoints) == 2 :
        return _linear_interpolate(buff, xpoints)
    elif len(xpoints) > 2:
        xpoints.sort()
        y = get_ypoints(buff, xpoints, nsmooth=nsmooth)
        tck = interpolate.splrep(xpoints, y, k=kind)
        def f(x):
            return interpolate.splev(x, tck, der=0, ext=0)
        return f
    else:  # if only one points given, returns a constant, which is the value at that point.
        raise NPKError("too little points in spline interpolation")
#-------------------------------------------------------------------------------
def _linear_interpolate(buff, xpoints, nsmooth=0):
    """computes and returns a linear interpolation"""
    xdata = np.array(xpoints)
    ydata = get_ypoints(buff, xpoints, nsmooth=nsmooth)
    coeffs = np.polyfit(xdata, ydata, 1)
    return np.poly1d(coeffs)
#-------------------------------------------------------------------------------
def _interpolate(func, npkd, xpoints, axis = 'F2', nsmooth=0):
    """"
    compute and applies a linear function as a baseline correction
    xpoints are the location of pivot points
    """
    if npkd.dim == 1:
        f = func(npkd.buffer, xpoints, nsmooth=nsmooth)
        x = np.arange(npkd.size1)
        npkd.buffer -= f(x)
    elif npkd.dim == 2:
        if npkd.test_axis(axis) == 2:
            x = np.arange(npkd.size2)
            for i in xrange(npkd.size1):
                f = func(npkd.buffer[i,:], xpoints, nsmooth=nsmooth)
                npkd.buffer[i,:] -= f(x)
        elif npkd.test_axis(axis) == 1:
            x = np.arange(npkd.size1)
            for i in xrange(npkd.size2):
                f = func(npkd.buffer[:,i], xpoints, nsmooth=nsmooth)
                npkd.buffer[:,i] -= f(x)
    else:
        raise NPKError("not implemented")
    return npkd

def linear_interpolate(npkd, xpoints, axis='F2', nsmooth=0):
    """"
    compute and applies a linear function as a baseline correction
    xpoints are the location of pivot points
    """
    return _interpolate(_linear_interpolate, npkd, xpoints, axis=axis, nsmooth=nsmooth )

def spline_interpolate(npkd, xpoints, axis='F2', nsmooth=0):
    """"
    compute and applies a spline function as a baseline correction
    xpoints are the location of pivot points
    """
    return _interpolate(_spline_interpolate, npkd, xpoints, axis=axis, nsmooth=nsmooth)


########################################################################
import spike.Algo.savitzky_golay as sgm
import spike.Algo.BC as BC
########################################################################
def bcorr_auto(npkd, iterations=10, nbchunks=40, degree=1, nbcores=2, smooth=True):
    """applies an automatic baseline correction
    
    Find baseline by using low norm value and then high norm value to attract the baseline on the small values.
    Parameters : 
    iterations : number of iterations for convergence toward the small values. 
    nbchunks : number of chunks on which is done the minimization. Typically, each chunk must be larger than the peaks. 
    degree : degree of the polynome used for approaching each signal chunk. 
    nbcores : number of cores used for minimizing in parallel on many chunks (if not None)
    
    smooth i True, applies a final Savitsky-Golay smoothing
    """
    npkd.check1D()
    bl = BC.correctbaseline(npkd.get_buffer(), iterations=iterations, nbchunks=nbchunks, degree=degree, nbcores=nbcores)
    #baseline(npkd.get_buffer(), degree=degree)
    if smooth:
        bl = sgm.savitzky_golay( bl, 205, 7)
    npkd.set_buffer( npkd.get_buffer() - bl)
    return npkd

def autopoints(npkd, Npoints=8, modulus=True):
    """
    computes Npoints (defaut 8) positions for a spline baseline correction
    """
    if Npoints is None :
        N = 8
    else:
        N = Npoints
    if npkd.itype == 1 and modulus:
        bf = npkd.copy().modulus().get_buffer()
    else:  
        bf = npkd.copy().get_buffer()
    bf -= np.percentile(bf,20)  # assumes at least ~20% of data-set is baseline...
    bf = abs(bf)
    L = len(bf)
    chunksize = L//N
    #print (chunksize)
    xpoints = np.array([i+bf[i:i+chunksize-8].argmin() for i in range(4, L, chunksize)])
    if npkd.itype == 1:
        xpoints *= 2
    return xpoints

def bcorr(npkd, method='spline', xpoints=None, nsmooth=0, modulus=True):
    """
    recapitulate all baseline correction methods, only 1D so far
    
    method is either
        auto: 
            use bcorr_auto, uses an automatic determination of the baseline
            does not work with negative peaks.
        linear:
            simple 1D correction
        spline:
            a cubic spline correction
    both linear and spline use an additional list of pivot points 'xpoints' used to calculate the baseline
    if xpoints absent,  pivots are estimated automaticaly
    if xpoints is integer, it determines the number of computed pivots (defaut is 8 if xpoints is None)
    if xpoints is a list of integers, there will used as pivots

    if nsmooth >0, buffer is smoothed by moving average over 2*nsmooth+1 positions around pivots.
    if dataset is complex, the xpoints are computed on the modulus spectrum, unless modulus is False

    default is spline with automatic detection of 8 baseline points
    """
    if method=='auto':
        return bcorr_auto(npkd)
    else:
        if xpoints is None or isinstance(xpoints,int):
            xpoints = autopoints(npkd, xpoints, modulus=modulus)
    if method=='linear':
        return linear_interpolate(npkd, xpoints, nsmooth=nsmooth)
    elif method=='spline':
        return spline_interpolate(npkd, xpoints, nsmooth=nsmooth)
    else:
        raise Exception("Wrong method in bcorr plugin")

NPKData_plugin("bcorr_lin", linear_interpolate)
NPKData_plugin("bcorr_spline", spline_interpolate)
NPKData_plugin("bcorr_auto", bcorr_auto)
NPKData_plugin("bcorr", bcorr)
