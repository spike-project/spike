#!/usr/bin/env python 
# encoding: utf-8

"""A utility for calibration of MS experiments

based on a list of experimentally measured m/z and theoretical ones
The method will reduce the difference between lists.

adds ppm() ppm_error() and display_icalib() methods to the FTMSAxis object
     and imzmeas and mzref attributes

Adds the following methods():
to FTICR datasets
    set_calib(mzmeas, mzref, axis=1)
    calib(axis=1, method='l1', verbose=False)
    display_calib(axis=1, compare=False)

to FTICR axes
    display_icalib(xref, mzref, symbol='bo')
    ppm_error(xref, mzref)
    ppm(xref, mzref)
    display_icalib(xref, mzref, symbol='bo')

and the following attributes:
    RefAxis: a backup FTICRAxis, used to store the previous calibration

to FTICR axes
    mzref : list of m/z of reference values
    imzmeas : list of pea indices of reference peaks (to be match with mzref)

"""

from __future__ import print_function
import numpy as np
from numpy.linalg import norm
from scipy.optimize import curve_fit, least_squares, leastsq, minimize
from spike import NPKError
from spike.NPKData import NPKData_plugin
from spike.util.signal_tools import findnoiselevel
import spike.Display.testplot as testplot
plt = testplot.plot()
from spike.FTMS import FTMSAxis

########################################################################
# mathematics
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
    return np.array([axis.mztoi(r) for r in meas])
    
   
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
        print( ' final - mean offset:', npk.axes(todo).ppm(xind, ref), 'ppm')
        print('**',npk.axis1.calibA, npk.axis1.calibB, npk.axis1.calibC)
    return npk


##########################################################################################
######### The following definitions are meant for being inserted into FTMS.FTMSAxis objects
##########################################################################################
def ppm_error(axis, xref, mzref):
    """
    computes the error from a array of positions (xref) and the theoretical m/z (mzref)
    returns an array with errors in ppm
    xref : array of point coordinates of the reference points
    mzref: array of reference m/z
    """
    return 1E6*(axis.itomz(xref)-mzref)/mzref
setattr(FTMSAxis, "ppm_error", ppm_error)

def ppm(axis, xref, mzref):
    """
    computes the mean error in ppm from a array of positions (xref) and the theoretical m/z (mzref)
    uses l1 norm !
    xref : array of point coordinates of the reference points
    mzref: array of reference m/z
    """
    shift = axis.ppm_error(xref, mzref)
    return abs(shift).mean()
setattr(FTMSAxis, "ppm", ppm)

def _display_icalib(axis, xref, mzref, symbol='bo'):
    """
    generates a plot of the current calibration
    xref : list of point coordinates of the reference points
    mzref: list of reference m/z
    symbol: matplotlib notation for points (default is blue rounds)
    """
    maxis = np.linspace(min(axis.mzref)-10, max(axis.mzref)+10,10)  # sample mass axis
    plt.plot(maxis, np.zeros_like(maxis),'k--')               # and draw a line at zero
    plt.plot(maxis, np.ones_like(maxis),'k:')
    plt.plot(maxis, -np.ones_like(maxis),'k:')
    plt.plot(mzref, axis.ppm_error(xref, mzref), symbol)
#    plt.plot(mzref,1E6*(axis.itomz(xref)-mzref)/mzref, symbol) # plot position of ref points
    plt.xlabel('$m/z$')
    plt.ylabel('$ppm$')
setattr(FTMSAxis, "display_icalib", _display_icalib)
setattr(FTMSAxis, "mzref", np.empty(10))  # template for reference m/z
setattr(FTMSAxis, "imzmeas", np.empty(10))  # template for index of measured m/z


##########################################################################################
######### high level tools
##########################################################################################
def display_calib(npkd,  axis=1, compare=False):
    """
    generates a plot of the current calibration
    if compare is True, will try to draw the previous calibration curve along with the current one
    """
    todo = npkd.test_axis(axis)
    axistodo = npkd.axes(todo)
#    xref = computes_index(meas, axistodo)
    plt.figure()
    axistodo.display_icalib(axistodo.imzmeas, axistodo.mzref)
    if compare:
        try:
            npkd.RefAxis.display_icalib(axistodo.imzmeas, axistodo.mzref, symbol='rx')
            maxis = np.linspace(min(axistodo.mzref)-10, max(axistodo.mzref)+10,1000)
            plt.plot(maxis, -axistodo.ppm_error( npkd.RefAxis.mztoi( maxis), maxis), 'r:' ) 
        except AttributeError:
            raise Exception('No Reference Axis available')
    plt.title("Mean error %.2f ppm"%(axistodo.ppm(axistodo.imzmeas, axistodo.mzref)))
    return npkd

def set_calib(npkd, mzmeas, mzref, axis=1):
    """
    preset parameters for mass calibration
    mzref : list of reference m/z of calibrants
    mzmeas : list of measured m/z of calibrants

    axis is the axis to be calibrated, defaut is 1
    """
    todo = npkd.test_axis(axis)
    npkd.axes(todo).mzref = np.array(mzref)
    npkd.axes(todo).imzmeas = npkd.axes(todo).mztoi( np.array(mzmeas) ) # convert to index (invariant !)
    return npkd

def calib(npk, axis=1, method='l1', verbose=False):
    """
    the current FTMS experiment is recalibrated optimatly 
    along its axis 'axis' (usefull only in 2D) using parameters provided with set_calib()
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
        print( 'before - mean offset:', ax2fit.ppm(ax2fit.imzmeas, ax2fit.mzref), 'ppm')
    fitted = mzcalib(ax2fit.imzmeas, ax2fit.mzref, ax2fit, method=method)
    if verbose:
        print( ' after - mean offset:', fitted.ppm(ax2fit.imzmeas, ax2fit.mzref), 'ppm')
    # and assign
    npk.axes(todo).calibA = fitted.calibA
    npk.axes(todo).calibB = fitted.calibB
    npk.axes(todo).calibC = fitted.calibC
    if verbose:
        print( ' final - mean offset:', npk.axes(todo).ppm(ax2fit.imzmeas, ax2fit.mzref), 'ppm')
        print('**',npk.axis1.calibA, npk.axis1.calibB, npk.axis1.calibC)
    return npk
 
# and plug the whole stuf into NPKData
NPKData_plugin("set_calib", set_calib)
NPKData_plugin("calib", calib)
NPKData_plugin("display_calib", display_calib)

