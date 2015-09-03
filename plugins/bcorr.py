#!/usr/bin/env python 
# encoding: utf-8

"""
set of function for the baseline correction

Very Sloppy - Not finsihed !
"""

from __future__ import print_function
import numpy as np
from scipy import interpolate
from scipy.optimize import leastsq

from spike import NPKError
from spike.NPKData import NPKData_plugin

#-------------------------------------------------------------------------------
def _spline_interpolate(buff, xpoints, kind = 3):
    """compute and returns a spline function 
        we are using splrep and splev instead of interp1d because interp1d needs to have 0 and last point
        it doesn't extend.
    """
    y = buff[point]
    if len(xpoints) > 2:
        tck = interpolate.splrep(xpoints, y, k=kind)
    elif len(xpoints) == 2 :
        tck = interpolate.splrep(xpoints, y, k=1)
    else:  # if only one points given, returns a constant, which is the value at that point.
        raise NPKError("too little points in spline interpolation")
    return interpolate.splev(np.arange(len(buff)), tck, der=0)

#-------------------------------------------------------------------------------
def _linear_interpolate(buff, xpoints):
    """computes and returns a linear interpolation"""
    xdata = np.array(xpoints)
    ydata = buff[xdata]
    coeffs = polyfit(xdata, ydata, 1)
    return np.poly1d(coeffs)
#-------------------------------------------------------------------------------
def _interpolate(func, npkd, xpoints, axis = 'F2'):
    """"
    compute and applies a linear function as a baseline correction
    xpoints are the location of pivot points
    """
    if npkd.dim == 1:
        f = _linear_interpolate(npkd.buffer, xpoints)
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
def baseline(npkd, degree=4, smooth=True):
    """applies a polynomial baseline correction"""
    npkd.check1D()
    bl = BC.baseline(npkd.get_buffer(), degree=degree)
    if smooth:
        bl = sgm.savitzky_golay( bl, 205, 7)
    self.set_buffer( npkd.get_buffer() - bl)
    return npkd

NPKData_plugin("bcorr_lin", linear_interpolate)
NPKData_plugin("bcorr_spline", spline_interpolate)
NPKData_plugin("baseline", baseline)
