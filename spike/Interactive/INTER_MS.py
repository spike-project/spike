#!/usr/bin/env python 
# encoding: utf-8

"""
A set of utilities to use spike in NMR or FTMS within jupyter


First version MAD june 2017
definitive ? version MAD october 2019

"""

from __future__ import print_function, division
import unittest
import sys
import os
import os.path as op
import glob

import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import MultiCursor
from ipywidgets import fixed, Layout, HBox, VBox, Label, Output, Button, Tab, HTML
import ipywidgets as widgets
from IPython.display import display, Javascript, Markdown
import numpy as np

from spike.File.BrukerNMR import Import_1D
from spike import NPKData
from spike.NMR import NMRData
from spike.Interactive.INTER import Show1D, jsalert
try:
    import spike.plugins.bcorr as bcorr
except:
    print('Baseline correction plugin not installed !')

try:
    from spike.plugins.MS.PhaseMS import movepivot
    PHASEQUAD = True
except:
    print('Quadratic phase correction plugin not installed !')
    PHASEQUAD = False

# REACTIVE modify callback behaviour
# True is good for inline mode / False is better for notebook mode
REACTIVE = True
HEAVY = False

class FloatButt(HBox):
    """
    a float text with buttons
    rounds is a list on power of 10 to set
    """
    def __init__(self, rounds, value=0, description="", **kw):
        super(FloatButt, self).__init__(**kw)
        self.field = widgets.FloatText(value=value, layout=Layout(width='20%'), step=0.1)
        self.incbuts = []
        self.decbuts = []
        self.description = description
        self.wdesc = Label(self.description, layout=Layout(width='20%'))
        for i in range(rounds):
            v = 10**i
            if i>2:
                des = "10^%d"%i
            else:
                des = "%d"%(v)
            b = widgets.Button(description="+"+des)
            b.val = v
            b.on_click( self.increm )
            self.incbuts.append(  b )
            b = widgets.Button(description="-"+des)
            b.val = v
            b.on_click( self.decrem )
            self.decbuts.append(  b )
        self.children = [self.wdesc] + list(reversed(self.decbuts)) + [self.field] + self.incbuts
    def increm(self, e):
        self.field.value += e.val
    def decrem(self, e):
        self.field.value -= e.val
    def _gvalue(self):
        return self.field.value
    def _svalue(self, value):
        self.field.value = value
    value = property(_gvalue, _svalue)

def firstguess(data):
    """
    compute phase initial guess - still quite imprecise
    """
    import math
    if data.params['SwpDir'] != '0':
        print('*** WARNING pulse with Increasing frequencies ***')
    sw = float(data.params['SW_h'])
    try:
        f0 = data.axis1.highfreq[0]
        f1 = data.axis1.lowfreq[0]
    except:
        f0 = data.axis1.highfreq
        f1 = data.axis1.lowfreq
    deltaf = f1-f0
    te = 3.5E-3
    p3 = float(data.params['P_3'])*1E-6
    l31 = int(data.params['L_31'])
    tau = (l31+1)*p3
    ph1 = tau*f1/deltaf + te 
    ph2 = -tau/deltaf/2
    return (ph1*sw, ph2*sw*sw)

def firstguess0(data):
    """
    as firstguess but for old Bruker systems (with acqus file)
    """
    import math
    acqu = data.params['acqu']
    sw = float(acqu['$SW_h'])
    f0 = sw
    f1 = float(acqu['$FR_low'])
    deltaf = f1-f0
    te = 3.5E-3
    p3 = float(acqu['$P'][3])*1E-6
    l31 = int(acqu['$L'][31])
    tau = (l31+1)*p3
    ph1 = tau*f1/deltaf + te 
    ph2 = -tau/deltaf/2
    return (ph1*sw, ph2*sw*sw)
class Phaser1D(Show1D):
    """
    An interactive phaser in 1D MS

        Phaser1D(spectrum)

    requires %matplotlib widget

    """
    def __init__(self, data, figsize=None, title=None, reverse_scroll=False, show=True):
        data.check1D()
        self.list = {}          # list is actually a dictionnary - sorry about that
                                # hold {pivot_in_points: (ph0, ph1)}
        if data.itype == 0:
            jsalert('Data is Real - Please redo Fourier Transform')
            return
        super().__init__( data, figsize=figsize, title=title, reverse_scroll=reverse_scroll, show=False)
        self.p0 = widgets.FloatSlider(description='P0:',min=-200, max=200, step=1,
                            layout=Layout(width='50%'), continuous_update=HEAVY)
        
        self.p1 = widgets.FloatSlider(description='P1:',min=-10000000, max=10000000, step=1000.0,
                            layout=Layout(width='100%'), continuous_update=HEAVY)
        self.p1 = FloatButt(4, description='P1  ', layout=Layout(width='100%'))
        self.p2 = FloatButt(4, description='P2  ', layout=Layout(width='100%'))
        P1, P2 = firstguess(data)
        self.p1.value, self.p2.value = round(P1), round(P2)
        if self.data.axis1.currentunit == 'm/z':
            pvmin = self.data.axis1.lowmass
            pvmax = self.data.axis1.highmass
        else:
            pvmin = 0
            pvmax = self.data.axis1.itoc(self.data.size1)
        self.pivot = widgets.BoundedFloatText(description='Pivot',
                        value=0, 
                        min=self.data.axis1.itoc(0),
                        max=self.data.axis1.itoc(self.data.size1),
                        step=0.1, layout=Layout(width='20%'))
        self.cancel = widgets.Button(description="Cancel", button_style='warning')
        self.cancel.on_click(self.on_cancel)
        # remove done button and create an Apply one
        self.done.close()
        self.apply = widgets.Button(description="Done", button_style='success')
        self.apply.on_click(self.on_Apply)
        # list managment
        self.clearlist = widgets.Button(description="Clear")
        self.clearlist.on_click(self.on_clearlist)
        self.addlist = widgets.Button(description="Add Entry")
        self.addlist.on_click(self.on_addlist)
        self.listlabel = Label(value="0 entry")
        self.printlist = widgets.Button(description="Print")
        self.printlist.on_click(self.on_printlist)

        # draw HBox
        if PHASEQUAD:
            phbutton = [self.p1, self.p2]
        else:
            phbutton = [self.p1]
        self.children = [VBox([
                            HBox([self.apply, self.cancel, ]),
                            HBox([HTML('<i>Phase List Managment</i>'), self.clearlist, self.addlist, self.listlabel, self.printlist]),
                            HBox([self.p0, self.pivot, HTML('<i>set with right-click on spectrum</i>')]),
                            *phbutton,
                            HBox([VBox([self.blank, self.reset, self.scale]), self.fig.canvas]) ])]
        # add interaction
        for w in [self.p0, self.p1.field, self.p2.field, self.scale]:
            w.observe(self.ob)
        self.pivot.observe(self.on_movepivot)
        # add click event on spectral window
        def on_press(event):
            if event.button == 3:
                self.pivot.value = event.xdata
        cids = self.fig.canvas.mpl_connect('button_press_event', on_press)
        self.lp0, self.lp1, self.lp2, self.lpv = self.ppivot()
        if show: self.show()
    def on_clearlist(self, b):
        self.list = {}
        self.listlabel.value = '0 entry'
    def on_printlist(self, b):
        if self.list != {}:
            for pv in sorted(self.list.keys()):
                print("%.1f:\t(%.0f %.0f)"%(self.data.axis1.itoc(pv), *self.list[pv]))
    def on_addlist(self, b):
        pv = self.data.axis1.ctoi(self.pivot.value)
        self.list[pv] = ( self.p0.value, self.p1.value )
        if len(self.list)>1:
            self.listlabel.value = '%d entries'%len(self.list)
            P = np.array(list(self.list.keys()))
            L = np.array(list(self.list.values()))
            dPH = np.polyfit(P, L[:,1], 1)
            self.p2.value = 0.5*dPH[0]*self.data.size1
        else:
            self.listlabel.value = '%d entry'%len(self.list)
        self.xb = self.ax.get_xbound()   # get current zoom
        self.ax.plot(self.data.axis1.itoc(pv), self.data[2*(int(pv)//2)], 'ro')
        self.ax.set_xbound( self.xb )
    def show(self):
#        self.data.display(figure=self.ax)
        if PHASEQUAD:
            pv = self.data.axis1.ctoi(self.lpv)/self.data.size1  #  0...1
            self.data.copy().phaseMS(self.lp0, self.lp1, self.lp2, pv).display(scale=self.scale.value, new_fig=False, figure=self.ax)
        else:
            self.data.copy().phase(self.lp0, self.lp1).display(scale=self.scale.value, new_fig=False, figure=self.ax)

        self.xb = self.ax.get_xbound()  # initialize zoom
        ppos = self.pivot.value
        self.ax.plot([ppos,ppos], self.ax.get_ybound())
        for pv in self.list.keys():
            self.ax.plot(self.data.axis1.itoc(pv), 0, 'ro')
        display(self)
    def on_cancel(self, b):
        # self.p0.value = 0  # because widget remains active...
        # self.p1.value = 0
        self.close()
        print("no applied phase")
    def on_Apply(self, b):
        self.close()
        lp0, lp1, lp2, lpv = self.ppivot() # get centered values
        if PHASEQUAD:
            pv = self.data.axis1.ctoi(self.lpv)/self.data.size1
            lp00, lp10, lp20, _ = movepivot(lp0, lp1, lp2, pv, 0.0)
            self.data.phaseMS(lp00, lp10, lp20, 0)
            msg = "Applied: phaseMS(%.1f,  %.1f, %.1f, %.2f)"%(lp00, lp10, lp20, 0)
        else:
            self.data.phase(lp0, lp1)
            msg = "Applied: phase(%.1f,  %.1f)"%(lp0, lp1)
        self.disp()
        for pv in self.list.keys():
            self.ax.plot(self.data.axis1.itoc(pv), self.data[2*(int(pv)//2)], 'ro')
        self.on_done(b)
        print(msg)

    def ppivot(self):
        "converts current pivot values to centered ones - not used in PHASEQUAD"
        if PHASEQUAD:
                return self.p0.value, self.p1.value, self.p2.value, self.pivot.value
        else:
                pp = 1.0-(self.pivot.value/self.data.size1)
                return (self.p0.value + (pp-0.5)*self.p1.value, self.p1.value, self.p2.value, self.pivot.value)
    def ctopivot(self, p0, p1, p2, pv):
        "convert from centered to pivot values - not used in PHASEQUAD"
        if PHASEQUAD:
                return p0, p1, p2, pv
        else:
                pp = 1.0-(self.pivot.value/self.data.size1)
                return p0- (pp-0.5)*p1, p1, 0, 0
    def on_movepivot(self, event):
        if event['name']=='value':
            if PHASEQUAD:
                pvaft = self.data.axis1.ctoi(self.pivot.value)/self.data.size1 
                pvbef = self.data.axis1.ctoi(self.lpv)/self.data.size1
                self.p0.value, self.p1.value, self.p2.value, _ = movepivot(self.lp0, self.lp1, self.lp2, pvbef, pvaft)
            else:
                self.p0.value, self.p1.value = self.ctopivot(self.lp0, self.lp1)
            self.phase()
    def ob(self, event):
        "observe changes and start phasing"
        if event['name']=='value':
            self.phase()
    def phase(self):
        "apply phase and display"
        self.xb = self.ax.get_xbound()   # get current zoom
        self.ax.clear()
        self.lp0, self.lp1, self.lp2, self.lpv = self.ppivot()         # get centered values
        if PHASEQUAD:
            pv = self.data.axis1.ctoi(self.lpv)/self.data.size1  #  0...1
            self.data.copy().phaseMS(self.lp0, self.lp1, self.lp2, pv).display(scale=self.scale.value, new_fig=False, figure=self.ax)
        else:
            self.data.copy().phase(self.lp0, self.lp1).display(scale=self.scale.value, new_fig=False, figure=self.ax)
        ppos = self.pivot.value
        self.ax.plot([ppos,ppos], self.ax.get_ybound())
        for pv in self.list.keys():
            self.ax.plot(self.data.axis1.itoc(pv), 0, 'ro')
        self.ax.set_xbound( self.xb )

# if __name__ == '__main__':
#    unittest.main()    
