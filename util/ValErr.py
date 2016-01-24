#!/usr/bin/env python 
# encoding: utf-8

"""
untitled.py

Created by Marc-André on 2011-03-20.
Copyright (c) 2011 IGBMC. All rights reserved.
"""

from __future__ import print_function

def round(x):
    """takes the nearest integer of a float
    >>> round(0.6)
    1.0
    >>> round(0.4)
    0.0
    >>> round(-0.6)
    -1.0
    >>> round(-0.4)
    0.0
    """
    import math
    return math.floor(x+0.5)

def ValErr(v,e=0.0):
    """
    nice print of a value with error bars

    values are rounded to significative digits
    >>> print ValErr(1234.567,123.4)
    1230 +/- 123

    can be called with a pair of values or a tuple :
    >>> V=(1234.567,123.4); print ValErr(V)
    1230 +/- 123

    values are rounded to significative digits
    >>> print ValErr(7654.321,1234)
    7700 +/- 1230
    >>> print ValErr(7654.321,123.4)
    7650 +/- 123
    >>> print ValErr(7654.321,12.34)
    7654 +/- 12

    format depends on error
    >>> print ValErr(1234.567,12.34)
    1235 +/- 12
    >>> print ValErr(1234.567,1.234)
    1234.6 +/- 1.2
    >>> print ValErr(1234.567,0.1234)
    1234.57 +/- 0.12
    >>> print ValErr(1234.567,0.01234)
    1234.567 +/- 0.012
    """
    import math
    deb = 0 # set it to one for debugging info
    if deb: print((v,e))
    if isinstance(v,tuple): # if v is tuple, unpack value and error
        (val,err) = v
    else:   # assume v is value
        val = v
        err = e
    try:
        dig = math.log10(err)
    except:
        dig = 3     # if err = 0
    # the larger the value of tweek (within 0 .. 1), the more the number of digits on the error value
    # 0 : always 1      1 : always 2
    tweek=1     # auto tests are for tweek=1
    idig = round(dig-tweek+0.5)
    if deb: print(dig,idig)
    err = float(round(err*10**(2-idig))) / 10**(2-idig)
    if deb : print(err)
    # choose depending on the size of err
    if (idig > 1):
        format = "%.0f +/- %.0f"
    elif (idig<=1):
        format = "%%.%df +/- %%.%df"%((-idig+1),(-idig+1))
    if deb : print(format)
    if (idig>0):
        val = float(round(val*10**(1-idig))) / 10**(1-idig)
    return format%(val,err)
