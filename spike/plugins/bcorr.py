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
def _spline_interpolate(buff, xpoints, kind = 3):
    """compute and returns a spline function 
        we are using splrep and splev instead of interp1d because interp1d needs to have 0 and last point
        it doesn't extend.
    """
    if len(xpoints) == 2 :
        return _linear_interpolate(buff, xpoints)
    elif len(xpoints) > 2:
        xpoints.sort()
        y = buff[xpoints]
        tck = interpolate.splrep(xpoints, y, k=kind)
        def f(x):
            return interpolate.splev(x, tck, der=0, ext=0)
        return f
    else:  # if only one points given, returns a constant, which is the value at that point.
        raise NPKError("too little points in spline interpolation")
#-------------------------------------------------------------------------------
def _linear_interpolate(buff, xpoints):
    """computes and returns a linear interpolation"""
    xdata = np.array(xpoints)
    ydata = buff[xdata]
    coeffs = np.polyfit(xdata, ydata, 1)
    return np.poly1d(coeffs)
#-------------------------------------------------------------------------------
def _interpolate(func, npkd, xpoints, axis = 'F2'):
    """"
    compute and applies a linear function as a baseline correction
    xpoints are the location of pivot points
    """
    if npkd.dim == 1:
        f = func(npkd.buffer, xpoints)
        x = np.arange(npkd.size1)
        npkd.buffer -= f(x)
    elif npkd.dim == 2:
        if npkd.test_axis(axis) == 2:
            x = np.arange(npkd.size2)
            for i in xrange(npkd.size1):
                f = func(npkd.buffer[i,:], xpoints)
                npkd.buffer[i,:] -= f(x)
        elif npkd.test_axis(axis) == 1:
            x = np.arange(npkd.size1)
            for i in xrange(npkd.size2):
                f = func(npkd.buffer[:,i], xpoints)
                npkd.buffer[:,i] -= f(x)
    else:
        raise NPKError("not implemented")
    return npkd

def linear_interpolate(npkd, xpoints, axis='F2'):
    """"
    compute and applies a linear function as a baseline correction
    xpoints are the location of pivot points
    """
    return _interpolate(_linear_interpolate, npkd, xpoints, axis=axis)

def spline_interpolate(npkd, xpoints, axis='F2'):
    """"
    compute and applies a spline function as a baseline correction
    xpoints are the location of pivot points
    """
    return _interpolate(_spline_interpolate, npkd, xpoints, axis=axis)


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

def bcorr(npkd, method='spline', xpoints=None):
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

    default is spline with automatic detection of 8 baseline points
    """
    if method=='auto':
        return bcorr_auto(npkd)
    else:
        if xpoints is None or isinstance(xpoints,int):
            if xpoints is None :
                N = 8
            else:
                N = xpoints
            bf = abs(npkd.get_buffer())
            L = len(bf)
            chunksize = L//N
            #print (chunksize)
            xpoints = np.array([i+bf[i:i+chunksize-8].argmin() for i in range(4, L, chunksize)])
            if npkd.itype == 1:
                xpoints *= 2
            #print (xpoints)
    if method=='linear':
        return linear_interpolate(npkd, xpoints)
    elif method=='spline':
        return spline_interpolate(npkd, xpoints)
    else:
        raise Exception("Wrong method in bcorr plugin")

NPKData_plugin("bcorr_lin", linear_interpolate)
NPKData_plugin("bcorr_spline", spline_interpolate)
NPKData_plugin("bcorr_auto", bcorr_auto)
NPKData_plugin("bcorr", bcorr)
