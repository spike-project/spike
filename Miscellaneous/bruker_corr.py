#!/usr/bin/env python 
# encoding: utf-8

"""
bruker_corr.py

Created by Marc-André on 2010-06-21.
Copyright (c) 2010 IGBMC. All rights reserved.

This file contains bruker specific stuff, as digital filter correction
DSPFIRM is named form bruker document, however, it seems it is called dspfvs in the header bruker files (acqu)
It is thus called accordingly in the Gifa User Interface.
"""
from __future__ import print_function

#/* this table holds the delay values. -1 indicates an invalid entry */
delay = [
    [ 179,201,533,709,1097,1449,2225,2929,4481,5889,8993,11809,18017,23649,36065,47329,72161, 94689,144353,189409,288737 ],
    [ 184,219,384,602, 852,1668,2312,3368,4656,6768,9344,13568,18560,27392,36992,50040,73856,110336,147584,220928,295040 ],
    [ 184,219,384,602, 852,1668,2292,3369,4616,6768,9264,13568,18560,27392,36992,50040,73856,110336,147584,220928,295040 ],
    [  11, 17, 23, 35,  47,  71,  95, 143, 191, 287, 383,  575,   -1,   -1,   -1,   -1,   -1,    -1,    -1,    -1,    -1 ],
    [  60, 90,118,179, 244, 360, 492, 724, 980,1444,1958, 2886, 3912, 5768, 7820,11532,   -1,    -1,    -1,    -1,    -1 ],
    [  -1, -1, 58,152, 202, 318, 418, 642, 842,1290,1690, 2586, 3386,   -1,   -1,   -1,   -1,    -1,    -1,    -1,    -1 ] ]

#/* this one contains valid DECIM values */
decim_offset = [ 2,3,4,6,8,12,16,24,32,48,64,96,128,192,256,384,512,768,1024,1536,2048 ]

N_DECIM = len(decim_offset)
N_DSPFIRM = len(delay)


#**********************************************************************/
def brukerdelay(dspfvs, dspfirm, decim, sw):
    """
    The brukerphase function computes the correction for the Bruker figital filter.
    the phase correction to apply is computed given the 3 parameters :
    DSPFIRM DSPFVS DECIM
    as found in the acqus parameter file in XwinNMR

    dspfvs is not used so far
    version 1.0
    10 oct 2001 - M.A.D.
    """
#/* first special cases */
    if (dspfirm == 15 and decim == 3  and sw >= 104000.):
        z = (-110.0*180.0)/(decim)
        return(z)
    if (dspfirm == 0 or decim == 1):
        return(0.0)
#/* otherwise, First determines DECIM offset */
    j = -1
    for i in xrange(N_DECIM):
        if (decim_offset[i] == decim):
            j=i
            break
    if (j == -1):
        raise ("*** wrong value for decim %d\n"%decim)

#/* then get delay value */
    if (dspfirm < 10 or dspfirm > (N_DECIM+9)):
        raise ("*** wrong value for dspfirm %d\n"%dspfirm)
    d = -delay[dspfirm-10][j]
    if (d == 1):
        raise ("*** wrong parameter combination")
    z = (d*180.00)/(decim)
#    print ("decim = %d   delay = %d   phase = %f"%(j,d,z)
    return(z)

def brukerdelay_2(dspfirm, dspfvs, decim, sw):
  ph = brukerdelay( dspfirm, dspfvs, decim, sw)
  return ph*decim/180.0

