#!/usr/bin/env python
# encoding: utf-8

"""A set of tools for computing Integrals for 1D  NMR spectra

If present, it can guess integral zones from an existing peak-list
Add .integzones .integcurves and .integvalues  into NPKDataset

First version by DELSUC Marc-André on May-2019.

"""

from __future__ import print_function
import numpy as np
import unittest
from spike.NPKData import NPKData_plugin, parsezoom
from spike.NPKError import NPKError

def delintegrals(data):
    "horrible hack to remove integrals - everything should be rewritten OOP!"
    try:
        del data.integzones
        del data.integcurves
        del data.integvalues
    except:
        pass
    return data

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

def compsums(data, bias=0.0, calib=None):
    "from integral lists, computes curves and values, sets integcurves and integvalues"
    curves = []
    buff = (data.get_buffer().real-bias)/data.cpxsize1
    for (a,b) in data.integzones:
        curves.append(buff[int(a):int(b)].cumsum())
    data.integcurves = curves
    data.integvalues = np.array( [c[-1] for c in curves] )
    # then calibrate
    sumax = max(data.integvalues)
    if calib is None:
        data.integvalues *= 100/sumax   # max at 100
    else:
        data.integvalues *= calib
    return data

def display(data, integoff=0.3, integscale=0.5, color='red', label=False, zoom=None, figure=None):
    "displays integrals"
    import matplotlib.transforms as transforms
    from spike.Display import testplot
    plt = testplot.plot()
    if figure is None:
        ax = plt.subplot(111)
    else:
        ax = figure
    trans = transforms.blended_transform_factory( ax.transData, ax.transAxes )
    z1, z2 = parsezoom(data, zoom)
    sumax = max([c[-1] for c in data.integcurves])
    for ((a,b),c,v) in zip(data.integzones, data.integcurves, data.integvalues):
#        print(a,b,max(c)/sumax)
        if a>z2 or b<z1:
            continue   # we're outside
        xinteg = data.axis1.itoc( np.linspace(a,b,len(c)) )
        yinteg = integoff + integscale*c/sumax
        ax.plot(xinteg, yinteg, transform=trans, color=color)
        if label:
            ax.text(xinteg[-1],yinteg[-1],"%.2f"%v, transform=trans, color=color)

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

def calibrate(npkd, entry, value):
    """
    on a dataset alreat integrated, the integrals are adapted so that
    the given entry is set to the given value.
    """
    try:
        vals = npkd.integvalues
    except:
        raise NPKError('data set should be integrated with .integrate() first!',data=npkd)
    npkd.integvalues *= value/vals[entry]
    return npkd
def integrate(npkd, separation=3, wings=5, bias=0.0, calibration=None):
    """
    computes integral zones and values from peak list, 

    separation : if two peaks are less than separation x width n they are aggregated, default = 3
    wings : integrals sides are extended by wings x width, default = 5
    bias: this value is substracted to data before integration
    calibration: a coefficient to multiply all integrals / if None (default) largest is set at 100
    """
    compzones(npkd, separation=separation, wings=wings)
    compsums(npkd, bias=bias, calib=calibration)
    return npkd
#-------------------------------------------------------
def int2pandas(npkd):
    "export extract of current integrals list to pandas Dataframe"
    import pandas as pd
    I1 = pd.DataFrame({
        'Start': [npkd.axis1.itoc(ii[0]) for ii in npkd.integzones],
        'End': [npkd.axis1.itoc(ii[1]) for ii in npkd.integzones],
        'Value': [c[-1] for c in npkd.integcurves],
        'Calibration': npkd.integvalues
    })
    return I1

NPKData_plugin("integrate", integrate)
NPKData_plugin("calibrate", calibrate)
NPKData_plugin("display_integrals", display)
NPKData_plugin("int2pandas", int2pandas)
NPKData_plugin("delintegrals", delintegrals)

