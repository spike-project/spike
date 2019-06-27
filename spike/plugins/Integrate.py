#!/usr/bin/env python
# encoding: utf-8

"""A set of tools for computing Integrals for 1D  NMR spectra

If present, it can guess integral zones from an existing peak-list
Adds .integrals into NPKDataset  which is an object with its own methods.

First version by DELSUC Marc-André on May-2019.

"""

from __future__ import print_function
import numpy as np
import unittest
from spike.NPKData import NPKData_plugin, parsezoom
from spike.NPKError import NPKError


class Integralitem(object):
    def __init__(self, start, end, curve, value):
        """
        the elemental integral item - used by Integrals
        start, end : the delimited zone, in pixels
        curve : the cumsum over the zone (eventually modified)
        value : the calibrated value
        """
        self.start = start
        self.end = end
        self.curve = curve
        self.value = value
    def _report(self):
        return "%d - %d : %f  on %d points\n"%(self.start, self.end, self.value, len(self.curve))
    def __str__(self):
        return self._report()
    def __repr__(self):
        return "Integralitem %s"%self._report()
class Integrals(list):
    """
    the class to hold a list of Integral

    an item is [start, end, curve as np.array(), value]
    start and end are in index !

    """
    def __init__(self, data, *args, calibration=None, bias=0.0, separation=3, wings=5, compute=True, **kwds):
        """
        computes integral zones and values from peak list, 

        separation : if two peaks are less than separation x width n they are aggregated, default = 3
        wings : integrals sides are extended by wings x width, default = 5
        bias: this value is substracted to data before integration
        calibration: a coefficient to multiply all integrals / if None (default) largest is set at 100
        """
        # I can't figure out how to explictly specify a keyword arg with *args:
        #   def __init__(self, *arg, threshold=None, source=None): ...
        # so I use **kwds and sqauwk if something unexpected is passed in.
        # code taken from   lib/python2.7/pstats.py
        #
        # additional kw are source: the originating dataset, compute: initialize values
        self.source = data
        self.calibration = calibration  # global calibration
        self.bias = bias          # global bias
        self.separation = separation
        self.wings = wings
        self.compute = compute
        super(Integrals, self).__init__(*args, **kwds)
        if self.compute:
            self.do_compute()

    def do_compute(self):
        "realize the computation, using internal parameters"
        self.peakstozones()
        self.zonestocurves()
    def _report(self):
        ll = ['calibration: %f\n'%self.calibration]
        for i,ii in enumerate(self):
            ll.append( "%d: %s"%(i,ii._report ))
        return "\n".join(ll)
    def report(self):
        for ii in self:
            print(ii)
    def to_pandas(self):
        "export extract of current integrals list to pandas Dataframe"
        import pandas as pd
        I1 = pd.DataFrame({
            'Start': [self.source.axis1.itoc(ii.start) for ii in self],
            'End': [self.source.axis1.itoc(ii.end) for ii in self],
            'Value': [ii.curve[-1] for ii in self],
            'Calibration': self.integvalues
        })
        return I1
    def peakstozones(self):
        """
        computes integrals zones from peak list, 
        separation : if two peaks are less than separation x width n they are aggregated, default = 3
        wings : integrals sides are extended by wings x width, default = 5
        """
        data = self.source
        try:
            pk = data.peaks
        except AttributeError:
            data.pp().centroid()            # create one if missing
            pk = data.peaks
        # then build integral list
        if len(pk) == 0:
                return []                   # return [] if empty peaklist
        prev = data.peaks[0]    # initialize
        start = prev.pos - self.wings*prev.width
        for pk in data.peaks[1:]: # then through remaining
            # extending or creating a new zone
            if (pk.pos - self.separation*pk.width) > (prev.pos + self.separation*prev.width): # we're done
                end = prev.pos + self.wings*prev.width
                self.append(   Integralitem(start, end, [], 0.0)  )
                start = pk.pos - self.wings*pk.width
            prev = pk
        end = data.peaks[-1].pos + self.wings*data.peaks[-1].width
        self.append(   Integralitem(start, end, [], 0.0)  )
    def zonestocurves(self):
        "from integral lists, computes curves and values"
        curves = []
        buff = (self.source.get_buffer().real-self.bias)/self.source.cpxsize1
        for iint in self:
            curves = buff[int(iint.start):int(iint.end)].cumsum()
            iint.curve = curves
            iint.value = curves[-1]
        # then calibrate
        self.calibrate(calibration=self.calibration)
    def calibrate(self, calibration=None):
        """computes integration values from curves
        either use calibration value as a scale,
        if calibration is None put the largest to 100.0 
        """
        if not calibration:
            intmax = 0.0
            for iint in self:
                iint.value = iint.curve[-1]
                intmax = max(intmax,iint.value)
            calibration = 100/intmax
        for iint in self:
            iint.value = iint.curve[-1]*calibration
        self.calibration = calibration
    @property
    def integzones(self):
        "the list of (start, end) of integral zones"
        return [(iint.start, iint.end) for iint in self]
    @property
    def integvalues(self):
        "the list of calibrated values"
        return [iint.value for iint in self]
    def recalibrate(self, entry, calib_value):
        """
        on a dataset already integrated, the integrals are adapted so that
        the given entry is set to the given value.
        """
        self.calibrate(calibration = calib_value/self[entry].curve[-1])
    def display(self, integoff=0.3, integscale=0.5, color='red', label=False, 
        labelyposition=None, regions=False, zoom=None, figure=None):
        "displays integrals"
        import matplotlib.transforms as transforms
        from spike.Display import testplot
        plt = testplot.plot()
        if figure is None:
            ax = plt.subplot(111)
        else:
            ax = figure
        trans = transforms.blended_transform_factory( ax.transData, ax.transAxes )
        z1, z2 = parsezoom(self.source, zoom)
        sumax = max([c.curve[-1] for c in self])
        for iint in self:
    #        print(a,b,max(c)/sumax)
            if iint.start>z2 or iint.end<z1:
                continue   # we're outside
            xinteg = self.source.axis1.itoc( np.linspace(iint.start, iint.end, len(iint.curve)) )
            yinteg = integoff + integscale*iint.curve/sumax
            ax.plot(xinteg, yinteg, transform=trans, color=color)
            if label:
                if labelyposition:
                    xl = xinteg[0] + 0.3*(xinteg[-1]- xinteg[0])
                    yl = labelyposition
                else:
                    xl = xinteg[-1]
                    yl = yinteg[-1]
                ax.text(xl,yl,"%.2f"%iint.value, transform=trans, color=color, fontsize=7)
            if regions:
                ax.plot([xinteg[0],xinteg[0]], [0,1], transform=trans, color=color, alpha=0.1)
                ax.plot([xinteg[-1],xinteg[-1]], [0,1], transform=trans, color=color, alpha=0.1 )

def integrate(npkd, **kw):
    """
    computes integral zones and values from peak list, 

    separation : if two peaks are less than separation x width n they are aggregated, default = 3
    wings : integrals sides are extended by wings x width, default = 5
    bias: this value is substracted to data before integration
    calibration: a coefficient to multiply all integrals / if None (default) largest is set at 100
    """
    I1 = Integrals(npkd, **kw)
    npkd.integrals = I1
    return npkd

def calibrate(npkd, entry, calib_value):
    npkd.integrals.recalibrate(entry, calib_value)
    return npkd
calibrate.__doc__ = Integrals.recalibrate.__doc__

def display(npkd, integoff=0.3, integscale=0.5, color='red', label=False, 
        labelyposition=None, regions=False, zoom=None, figure=None):
    npkd.integrals.display(integoff=integoff, integscale=integscale, color=color, label=label, 
        labelyposition=labelyposition, regions=regions, zoom=zoom, figure=figure)
    return npkd

NPKData_plugin("integrate", integrate)
NPKData_plugin("integral_calibrate", calibrate)
NPKData_plugin("display_integral", display)

