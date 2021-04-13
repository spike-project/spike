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

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import MultiCursor
from ipywidgets import fixed, Layout, HBox, VBox, Label, Output, Button, Tab
import ipywidgets as widgets
from IPython.display import display, HTML, Javascript
import numpy as np

from ..File.BrukerNMR import Import_1D
from .. import NPKData
from ..NMR import NMRData
from ..Interactive.INTER import Show1D, Colors, baseline1D

try:
    import spike.plugins.bcorr as bcorr
except:
    print('Baseline correction plugins not installed !')

# REACTIVE modify callback behaviour
# True is good for inline mode / False is better for notebook mode
REACTIVE = True
HEAVY = False



class baseline2D_F2(baseline1D):
    def __init__(self, data, figsize=None):
        print('WARNING this tool is not functional/tested yet')
        self.data2D = data
        super().__init__( self.data2D.projF2, figsize=figsize)
    def on_done(self, e):
        super().on_done(e)
        ibsl_points = [int(self.data2D.axis2.ptoi(x)) for x in self.bsl_points]
        self.data2D.bcorr(method='spline', xpoints=ibsl_points)

class Show2D(object):
    """
    A display for 2D NMR with a scale cursor
    Show2D(spectrum) where spectrum is a NPKData object
    - special display for DOSY.
    """
    def __init__(self, data, title=None, figsize=None):
        self.data = data
        self.isDOSY =  isinstance(data.axis1, NPKData.LaplaceAxis)
        try:
            self.proj2 = data.projF2
        except:
            self.proj2 = data.proj(axis=2).real()
        try:
            self.proj1 = data.projF1
        except:
            self.proj1 = data.proj(axis=1).real()
        self.title = title
        self.scale = widgets.FloatLogSlider(description='scale:', value=1.0, min=-1, max=3,  base=10, step=0.01,
                            layout=Layout(width='80%'), continuous_update=HEAVY)
        self.posview = widgets.Checkbox(value=True,description='Positive', tooltip='Display Positive levels')
        self.negview = widgets.Checkbox(value=False,description='Negative', tooltip='Display Negative levels')
        self.cursors = widgets.Checkbox(value=False,description='Cursors', tooltip='show cursors (cpu intensive !)')
        for w in (self.scale, self.posview, self.negview, self.cursors):
            w.observe(self.ob)
        grid = {'height_ratios':[1,4],'hspace':0,'wspace':0}
        if self.isDOSY:
            fsize = (10,5)
            grid['width_ratios']=[7,1]
        else:
            fsize = (8,8)
            grid['width_ratios']=[4,1]
#        fig, self.axarr = plt.subplots(2, 1, sharex=True, figsize=fsize, gridspec_kw=grid)
        self.fig = plt.figure(figsize=fsize, constrained_layout=False)
        spec2 = gridspec.GridSpec(ncols=2, nrows=2, figure=self.fig, **grid)
        axarr = np.empty((2,2), dtype=object)
        axarr[0,0] = self.fig.add_subplot(spec2[0, 0])
        axarr[1,0] = self.fig.add_subplot(spec2[1, 0],sharex=axarr[0, 0])
        axarr[1,1] = self.fig.add_subplot(spec2[1, 1],sharey=axarr[1, 0])
        self.top_ax = axarr[0,0]
        self.spec_ax = axarr[1,0]
        self.side_ax = axarr[1,1]
        self.multitop = None
        self.multiside = None
        self.Box = HBox( [self.scale, self.posview, self.negview, self.cursors])
        display( self.Box )
        self.disp(new=True)
    def on_done(self, b):
        self.scale.close()
    def ob(self, event):
        "observe events and display"
        if event['name'] != 'value':
            return
        self.disp()
    def disp(self,new=False):
        if new:
            self.proj2.display(figure=self.top_ax, title=self.title)
            xb = self.top_ax.get_xbound()
            dataxis = self.proj1.axis1.itoc( self.proj1.axis1.points_axis() )
            self.side_ax.plot(self.proj1.get_buffer(),dataxis)
            yb = self.side_ax.get_ybound()
        else:
            yb = self.side_ax.get_ybound()
            xb = self.top_ax.get_xbound()
            self.spec_ax.clear()
        if self.cursors.value:
            self.multitop = MultiCursor(self.fig.canvas, (self.spec_ax, self.top_ax), color='r', lw=1, horizOn=False, vertOn=True)
            self.multiside = MultiCursor(self.fig.canvas, (self.spec_ax, self.side_ax), color='r', lw=1, horizOn=True, vertOn=False)
        else:
            self.multitop = None
            self.multiside = None
        if self.posview.value:
            self.data.display(scale=self.scale.value, new_fig=False, figure=self.spec_ax)
        if self.negview.value:
            self.data.display(scale=-self.scale.value, new_fig=False,
                figure=self.spec_ax, mpldic={'cmap':'Wistia'})
        self.spec_ax.set_xbound(xb)
        self.spec_ax.set_ybound(yb)



class Phaser2D(object):
    """
    An interactive phaser in 2D NMR

        Phaser2D(spec)

    best when in %matplotlib inline

    """
    def __init__(self, data):
        if data.itype != 3:
            print('Dataset should be complex along both axes, Phasing is not possible')
            return
        print('WARNING this tool is not functional/tested yet')
        self.data = data
        self.scale = widgets.FloatLogSlider(description='scale:', value=1.0, min=-1, max=2,  base=10, step=0.01,
                            layout=Layout(width='80%'), continuous_update=HEAVY)
        self.F1p0 = widgets.FloatSlider(min=-180, max=180, step=0.1, description='F1 p0',continuous_update=HEAVY)
        self.F1p1 = widgets.FloatSlider(min=-250, max=250, step=1.0, description='F1 p1',continuous_update=HEAVY)
        self.F2p0 = widgets.FloatSlider(min=-180, max=180, step=0.1, description='F2 p0',continuous_update=HEAVY)
        self.F2p1 = widgets.FloatSlider(min=-250, max=250, step=1.0, description='F2 p1',continuous_update=HEAVY)
        for w in [self.F1p0, self.F1p1, self.F2p0, self.F2p1, self.scale]:
            w.observe(self.ob)
        self.button = widgets.Button(description="Apply correction",button_style='success')
        self.button.on_click(self.on_Apply)
        self.cancel = widgets.Button(description="Cancel",button_style='warning')
        self.cancel.on_click(self.on_cancel)
#       interact(self.phase, scale=self.scale, F1p0=self.F1p0, F1p1=self.F1p1, F2p0=self.F2p0, F2p1=self.F2p1)
        display(VBox([self.scale,
            HBox([VBox([self.F1p0, self.F1p1], layout=Layout(width='40%')), VBox([self.F2p0, self.F2p1],  layout=Layout(width='40%'))], layout=Layout(width='80%'))
            ], layout=Layout(width='100%')))
        display(HBox([self.button, self.cancel]))
        fi,ax = plt.subplots()
        self.ax = ax
        self.display()
        #self.data.display(figure=self.ax)
    def ob(self, event):
        "observe changes and start phasing"
        if event['name']=='value':
            self.phase()
    def close(self):
        for w in [self.F1p0, self.F1p1, self.F2p0, self.F2p1, self.scale, self.button, self.cancel]:
            w.close()
    def on_cancel(self, b):
        print("No action")
        self.ax.clear()
        self.data.display(figure=self.ax,scale=self.scale.value)
        self.ax.set_xlim(xmin=self.data.axis2.itop(0), xmax=self.data.axis2.itop(self.data.size2))
        self.ax.set_ylim(ymin=self.data.axis1.itop(0), ymax=self.data.axis1.itop(self.data.size1))
        self.close()
    def on_Apply(self, b):
        print("Applied: phase(%.1f,%.1f,axis='F1').phase(%.1f,%.1f,axis='F')"%(self.F1p0.value, self.F1p1.value, self.F2p0.value, self.F2p1.value))
        self.data.phase(self.F2p0.value, self.F2p1.value, axis='F2').phase(self.F1p0.value, self.F1p1.value, axis='F1')
        self.data.display(figure=self.ax,scale=self.scale.value)
        self.ax.set_xlim(xmin=self.data.axis2.itop(0), xmax=self.data.axis2.itop(self.data.size2))
        self.ax.set_ylim(ymin=self.data.axis1.itop(0), ymax=self.data.axis1.itop(self.data.size1))
        self.close()
    def display(self,todisplay=None):
        "display either the current data or the one provided - red and blue"
        self.ax.clear()
        if not todisplay:
            todisplay = self.data
        todisplay.display(scale=self.scale.value, new_fig=False, figure=self.ax,color='blue')
        todisplay.display(scale=-self.scale.value, new_fig=False, figure=self.ax, color='red')
        self.ax.set_xlim(xmin=self.data.axis2.itop(0), xmax=self.data.axis2.itop(self.data.size2))
        self.ax.set_ylim(ymin=self.data.axis1.itop(0), ymax=self.data.axis1.itop(self.data.size1))
    def phase(self):
        "compute phase and display"
        dp = self.data.copy().phase(self.F2p0.value, self.F2p1.value, axis='F2').phase(self.F1p0.value, self.F1p1.value, axis='F1')
        self.display(dp)
    # def phase(self, scale, F1p0, F1p1, F2p0, F2p1):
    #     self.data.copy().phase(F1p0,F1p1,axis='F1').phase(F2p0,F2p1,axis='F2').display(scale=scale);



# if __name__ == '__main__':
#    unittest.main()    
