#!/usr/bin/env python 
# encoding: utf-8

"""set of function Savitsky-Golay smoothing

"""

from __future__ import print_function
from spike import NPKError
from spike.NPKData import NPKData_plugin
from spike.util.signal_tools import findnoiselevel

########################################################################
def sg(npkd, window_size, order, deriv=0, axis=0):
    """applies Savitzky-Golay of order filter to data
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less than `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    axis: int
        the axis on which the filter is to be applied, default is fastest varying dimension
    """
    import spike.Algo.savitzky_golay as sgm
    todo = npkd.test_axis(axis)
    m = sgm.sgolay_coef(window_size, order, deriv=0)
    if npkd.dim == 1:
        npkd.set_buffer( sgm.sgolay_comp(npkd.get_buffer(), m, window_size) )
    elif npkd.dim == 2:
        if todo == 2:
            for i in xrange(npkd.size1):
                npkd.buffer[i,:] = sgm.sgolay_comp(npkd.buffer[i,:], m, window_size)
        elif todo == 1:
            for i in xrange(1,npkd.size2):
                npkd.buffer[:,i] = sgm.sgolay_comp(npkd.buffer[:,i], m, window_size)
    else:
        raise NPKError("a faire")
    return npkd
    
########################################################################
def sg2D(npkd, window_size, order, deriv=None):
    """applies a 2D Savitzky-Golay of order filter to data
    window_size : int
        the length of the square window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less than `window_size` - 1.
    deriv: None, 'col', or 'row'.   'both' mode does not work.
        the direction of the derivative to compute (default = None means only smoothing)
    can be applied to a 2D only.
    """
    import spike.Algo.savitzky_golay as sgm
    npkd.check2D()
    npkd.set_buffer( sgm.savitzky_golay2D(npkd.get_buffer(), window_size, order, derivative=deriv) )
    return npkd

NPKData_plugin("sg", sg)
NPKData_plugin("sg2D", sg2D)
