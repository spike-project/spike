#!/usr/bin/env python 
# encoding: utf-8

"""A utility to set to zero all points below a ratio

"""

from __future__ import print_function
from spike import NPKError
from spike.NPKData import NPKData_plugin
from spike.util.signal_tools import findnoiselevel

########################################################################
def fastclean(npkd, nsigma=2.0, nbseg=20, axis=0):
    """
    set to zeros all points below nsigma times the noise level
    This allows the corresponding data-set, once stored to file, to be considerably more compressive.
    
    nsigma: float
        the ratio used, typically 1.0 to 3.0 (higher compression)
    nbseg: int
        the number of segments used for noise evaluation, see util.signal_tools.findnoiselevel
    axis: int
        the axis on which the noise is evaluated, default is fastest varying dimension
    """
    todo = npkd.test_axis(axis)
    if npkd.dim == 1:
        noise = findnoiselevel(npkd.get_buffer(), nbseg=nbseg)
        npkd.zeroing(nsigma*noise)
    elif npkd.dim == 2:
        if todo == 2:
            for i in xrange(npkd.size1):
                npkd.set_row(i, npkd.row(i).fastclean(nsigma=nsigma, nbseg=nbseg))
        elif todo == 1:
            for i in xrange(npkd.size2):
                npkd.set_col(i, npkd.col(i).fastclean(nsigma=nsigma, nbseg=nbseg))
    else:
        raise NPKError("a faire")
    return npkd

NPKData_plugin("fastclean", fastclean)
