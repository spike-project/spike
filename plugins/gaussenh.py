#!/usr/bin/env python
# encoding: utf-8
"""Gaussian enhancement apodisation 

d.gaussenh(width, enhancement=1.0, axis=0)

    apply an gaussian enhancement, width is in Hz
    enhancement is the strength of the effect
    axis is either F1, or F2 in 2D, 0 is default axis.
    multiplies by gauss(width) * exp(-enhancement*width)

Created by DELSUC Marc-André on February 2019
Copyright (c) 2019 IGBMC. All rights reserved.
"""
import numpy as np
from spike.NPKData import as_float, NPKData_plugin

#-------------------------------------------------------------------------------
def gaussenh(npkd, width, enhancement=2.0, axis=0):
    """
    apply an gaussian enhancement, width is in Hz
    enhancement is the strength of the effect
    multiplies by gauss(width) * exp(-enhancement*width)
    """
    todo = npkd.test_axis(axis)
    it = npkd.axes(todo).itype
    sw = npkd.axes(todo).specwidth
    size = npkd.axes(todo).size
    if it == 1: # means complex
        size = size//2
    baseax = width*np.arange(size)/sw
    e = np.exp(enhancement*baseax)
    e *= np.exp(-(baseax)**2)   
    if it == 1:
        e = as_float((1 + 1.0j)*e)
    # check NPKData.py to see how apodisations are handled.
    return npkd.apod_apply(axis,e)

NPKData_plugin("gaussenh", gaussenh)
