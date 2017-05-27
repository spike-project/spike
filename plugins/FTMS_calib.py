#!/usr/bin/env python 
# encoding: utf-8

"""A utility for calibration of MS experiments

based on a list of experimentally measured m/z and theoretical ones
The method will reduce the difference between lists.

adds ppm() ppm_error() and _display_calib() methods to the FTMSAxis object
"""

from __future__ import print_function
import numpy as np
from numpy.linalg import norm
from scipy.optimize import curve_fit, least_squares, leastsq, minimize
from spike import NPKError
from spike.NPKData import NPKData_plugin
from spike.util.signal_tools import findnoiselevel
import spike.Display.testplot as testplot
from spike.FTMS import FTMSAxis

########################################################################
def mzcalib(xind, ref, axis, method='l1'):
    """
    fit axis parameters so that points located at xind be closest to ref values
    fits on two parameters if the 3rd (ML3 ie axis.calibC) is zero
    method = l2      uses levenberg-marqtart on l2 norm : classical method
    method = l1      uses powell on l1 norm : robust method
    """
    axfit = axis.copy()                              # create copy 
    if axfit.calibC == 0:
        Po = [axfit.calibA, axfit.calibB]
    else:
        Po = [axfit.calibA, axfit.calibB, axfit.calibC]
    if method == 'l2':
        Pn, info = leastsq(distcalib, Po, args=(xind, ref, axfit))
    elif method == 'l1':
        res =  minimize(l1calib, Po, args=(xind, ref, axfit), method='Powell')
        Pn = res.x
    else:
        raise Exception("method is either l1 or l2")
    if axfit.calibC == 0:
        axfit.calibA, axfit.calibB = Pn
    else:
        axfit.calibA, axfit.calibB, axfit.calibC = Pn
    return axfit

def l1calib(param, xref, mzref, axis):
    """
    computes the sum of residual to minimize when calibrating
    basically a wrapper around axis.ppm
    """
    axis.calibA = param[0]
    axis.calibB = param[1]
    try:
        axis.calibC = param[2]
    except IndexError:
        pass
    return axis.ppm(xref, mzref)

def distcalib(param, xref, mzref, axis):
    """
    computes the residual to minimize when calibrating
    basically a wrapper around axis.ppm_error
    """
    axis.calibA = param[0]
    axis.calibB = param[1]
    try:
        axis.calibC = param[2]
    except IndexError:
        pass
    return axis.ppm_error(xref, mzref)

def computes_index(meas, axis):
    "locates in index measured values"
    # locates in index measured values
    return np.array([axis.mztoi(r) for r in meas])
    
def calib(npk, meas, ref, axis=1, method='l1', verbose=False):
    """
    given a list of measured m/z 'meas' and of theoretical values 'ref'
    the current FTMS experiment is recalibrated optimatly along its axis 'axis' (usefull only in 2D)
    uses the current (2 or 3 parameters) calibration

    method is either 'l1' (robust) or 'l2' (classic)

    The current calibration is copied to a new unused axis called RefAxis
    """
    # find axis
    todo = npk.test_axis(axis)
    # locates in index measured values
    xind = computes_index(meas, todo)
    return icalib(npk, xind, ref, axis=axis, method=method, verbose=verbose)
    
def icalib(npk, xind, ref, axis=1, method='l1', verbose=False):
    """
    given a list of location in index 'xind' and of theoretical values 'ref'
    the current FTMS experiment is recalibrated optimatly along its axis 'axis' (usefull only in 2D)
    uses the current (2 or 3 parameters) calibration
    
    method is either 'l1' (robust) or 'l2' (classic)
    
    The current calibration is copied to a new unused axis called RefAxis
    """
    # find axis
    todo = npk.test_axis(axis)
    npk.RefAxis = npk.axes(todo).copy()     # for latter comparison
    ax2fit = npk.axes(todo).copy()          # for action
    # fit
    if verbose:
        print( 'before - mean offset:', ax2fit.ppm(xind, ref), 'ppm')
    fitted = mzcalib(xind, ref, ax2fit, method=method)
    if verbose:
        print( ' after - mean offset:', fitted.ppm(xind, ref), 'ppm')
    # and assign
    npk.axes(todo).calibA = fitted.calibA
    npk.axes(todo).calibB = fitted.calibB
    npk.axes(todo).calibC = fitted.calibC
    if verbose:
        print( ' after - mean offset:', npk.axes(todo).ppm(xind, ref), 'ppm')
        print('**',npk.axis1.calibA, npk.axis1.calibB, npk.axis1.calibC)
    return npk

def display_calib(self, meas, mzref, axis=1, compare=False):
    """
    generates a plot of the current calibration
    is compare is True, will try to draw the previous calibration curve along with 
    """
    todo = self.test_axis(axis)
    axistodo = self.axes(todo)
    xref = computes_index(meas, axistodo)
    return display_icalib(self, xind, mzref, axis=axis, compare=compare)
def display_icalib(self, xind, mzref, axis=1, compare=False):
    """
    generates a plot of the current calibration
    is compare is True, will try to draw the previous calibration curve along with 
    """
    import spike.Display.testplot as testplot
    plt = testplot.plot()
    
    todo = self.test_axis(axis)
    axistodo = self.axes(todo)
    xref = computes_index(meas, axistodo)
    if compare:
        if self.RefAxis is None:
            print('Warning: No Reference Axis available')
        else:
            axistodo.display_calib(xref, mzref, self.RefAxis)
    else:
        axistodo.display_calib(xref, mzref)
    return self

######### The following definitions are meant for being inserted into FTMS.FTMSAxis objects
def ppm_error(self, xref, mzref):
    """
    computes the error from a list of positions (xref) and the theoretical m/z (mzref)
    returns an array with errors in ppm
    xref : list of point coordinates of the reference points
    mzref: list of reference m/z
    """
    return 1E6*(self.itomz(xref)-mzref)/mzref
setattr(FTMSAxis, "ppm_error", ppm_error)

def ppm(self, xref, mzref):
    """
    computes the mean error in ppm from a list of positions (xref) and the theoretical m/z (mzref)
    uses l1 norm !
    xref : list of point coordinates of the reference points
    mzref: list of reference m/z
    """
    shift = self.ppm_error(xref, mzref)
    return abs(shift).mean()
setattr(FTMSAxis, "ppm", ppm)

def _display_calib(self, xref, mzref, RefAxis = None):
    """
    generates a plot of the current calibration
    xref : list of point coordinates of the reference points
    mzref: list of reference m/z
    if RefAxis != None:
        it should be a FTMS.FTMSAxis, then RefAxis calibration is drawn in superposition
    """
    plt = testplot.plot()
    maxis = np.linspace(self.lowmass, self.highmass,1000)  # sample mass axis
    plt.plot(maxis, np.zeros_like(maxis),'k--')               # and draw a line at zero
    plt.plot(maxis, np.ones_like(maxis),'k:')
    plt.plot(maxis, -np.ones_like(maxis),'k:')
    plt.plot(mzref,1E6*(self.itomz(xref)-mzref)/mzref,'o') # plot position of ref points
    if RefAxis is not None:
        plt.plot(mzref,1E6*(RefAxis.itomz(xref)-mzref)/mzref,'rx') # plot position of ref points        
        plt.plot(maxis, -self.ppm_error( RefAxis.mztoi( maxis), maxis), 'r:' )
    plt.xlabel('$m/z$')
    plt.ylabel('$ppm$')
setattr(FTMSAxis, "display_calib", _display_calib)


# and plug the whole stuf into NPKData
NPKData_plugin("calib", calib)
NPKData_plugin("display_calib", display_calib)

