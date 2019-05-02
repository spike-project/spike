#!/usr/bin/env python
# encoding: utf-8

"""A set of tools for computing Integrals for 1D  NMR spectra

If present, it can guess integral zones from an existing peak-list
Add .integzones .integcurves and .integvalues into NPKDataset

First version by DELSUC Marc-André on May-2019.

"""

from __future__ import print_function
import numpy as np
import unittest
from spike.NPKData import NPKData_plugin

def compzones(data, separation = 3, wings = 5):
    """
    computes integral zones from peak list, 
    return a list [(star,end)...] in index unit
    separation : if two peaks are less than separation x width n they are aggregated, default = 3
    wings : integrals sides are extended by wings x width, default = 5
    """
    try:
        if len(data.peaks) == 0:
            return []                   # return [] if empty peaklist
    except AttributeError:
        data.pp().centroid()            # create one if missing
    # then build integral list
    integrals = []
    prev = data.peaks[0]    # initialize
    start = prev.pos - wings*prev.width
    for pk in data.peaks[1:]: # then through remaining
        # extending or creating a new zone
        if (pk.pos - separation*pk.width) > (prev.pos + separation*prev.width): # we're done
            end = prev.pos + wings*prev.width
            integrals.append([start,end])
            start = pk.pos - wings*pk.width
        prev = pk
    end = data.peaks[-1].pos + wings*data.peaks[-1].width
    integrals.append([start,end])
    data.integzones = integrals
    return data

def compsums(data, bias=0.0):
    "from integral lists, computes curves and values, sets integcurves and integvalues"
    curves = []
    sums = []
    buff = (data.get_buffer().real-bias)/data.cpxsize1
    for (a,b) in data.integzones:
        curves.append(buff[int(a):int(b)].cumsum())
    sums = np.array( [c[-1] for c in curves] )
    data.integcurves = curves
    data.integvalues = sums
    return data

def display(data, zoom=None, integoff=0.3, integscale=0.5, color='red', label=False):
    "displays integrals"
    import matplotlib.transforms as transforms
    from spike.Display import testplot
    plt = testplot.plot()
    sumax = max(data.integvalues)
    ax = plt.gca()
    trans = transforms.blended_transform_factory( ax.transData, ax.transAxes )
    for ((a,b),c) in zip(data.integzones, data.integcurves):
#        print(a,b,max(c)/sumax)
        xinteg = data.axis1.itoc( np.linspace(a,b,len(c)) )
        yinteg = integoff + integscale*c/sumax
        plt.plot(xinteg, yinteg, transform=trans, color=color)
        if label:
            ax.text(xinteg[-1],yinteg[-1],"%.1f"%(100*c[-1]/sumax,), transform=trans, color=color)

class IntegralTests(unittest.TestCase):
    def setUp(self):
        self.verbose = 1    # verbose >0 switches on messages
    def announce(self):
        if self.verbose >0:
            print (self.shortDescription())
    def _test_log(self):
        """testing log"""
        import math
        self.announce()
        x = 0.0
        y = math.log(1.0)
        self.assertAlmostEqual(x, y )

def integrate(data, separation=3, wings=5, bias=0.0):
    """
    computes integral zones and values from peak list, 

    separation : if two peaks are less than separation x width n they are aggregated, default = 3
    wings : integrals sides are extended by wings x width, default = 5
    bias: this value is substracted to data before integration
    """
    compzones(data, separation=separation, wings=wings)
    compsums(data, bias=bias)
    return data

NPKData_plugin("integrate", integrate)
NPKData_plugin("display_integrals", display)
