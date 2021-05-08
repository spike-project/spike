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

class Show2D(Show1D):
    """
    A display for 2D NMR with a scale cursor
    Show2D(spectrum) where spectrum is a NPKData object
    - special display for DOSY.
    """
    def __init__(self, data, title=None, figsize=None):
        super().__init__(data, title=title, create_children=False)
        self.isDOSY =  isinstance(data.axis1, NPKData.LaplaceAxis)
        try:
            self.proj2 = data.projF2
        except:
            self.proj2 = data.proj(axis=2).real()
        try:
            self.proj1 = data.projF1
        except:
            self.proj1 = data.proj(axis=1).real()
        # Controls
        self.scale.min = 0.2
        self.posview = widgets.Checkbox(value=True,description='Positive', tooltip='Display Positive levels', layout=Layout(width='20%'))
        self.negview = widgets.Checkbox(value=False,description='Negative', tooltip='Display Negative levels', layout=Layout(width='20%'))
        self.cursors = widgets.Checkbox(value=False,description='Cursors', tooltip='show cursors (cpu intensive !)', layout=Layout(width='20%'))
        self.showlogo = widgets.Checkbox(description="Logo", value = True)
        def switchlogo(e):
            if self.showlogo.value:
                self.axlogo.set_visible(True)
            else:
                self.axlogo.set_visible(False)
        self.showlogo.observe(switchlogo)
        for w in (self.scale, self.posview, self.negview, self.cursors):
            w.observe(self.ob)
        # Grid
        grid = {'height_ratios':[1,4],'hspace':0,'wspace':0}
        if self.isDOSY:
            if figsize is None:
                fsize = (10,5)
            else:
                fsize = figsize
            grid['width_ratios']=[7,1]
        else:
            if figsize is None:
                fsize = (8,8)
            else:
                fsize = figsize
            grid['width_ratios']=[4,1]
#        fig, self.axarr = plt.subplots(2, 1, sharex=True, figsize=fsize, gridspec_kw=grid)
        # Figure
        plt.ioff()
        self.fig = plt.figure(figsize=fsize, constrained_layout=False, tight_layout=True)
        plt.ion()
        self.fig.canvas.toolbar_position = 'left'
        spec2 = gridspec.GridSpec(ncols=2, nrows=2, figure=self.fig, **grid)
        axarr = np.empty((2,2), dtype=object)
        axarr[0,0] = self.fig.add_subplot(spec2[0, 0])
        axarr[1,0] = self.fig.add_subplot(spec2[1, 0],sharex=axarr[0, 0])
        axarr[1,1] = self.fig.add_subplot(spec2[1, 1],sharey=axarr[1, 0])
        axarr[0,1] = self.fig.add_subplot(spec2[0, 1])
        self.top_ax = axarr[0,0]
        self.spec_ax = axarr[1,0]
        self.side_ax = axarr[1,1]
        self.axlogo = axarr[0,1]
        self.axlogo.set_visible(False)
        self.multitop = None
        self.multiside = None
        self.ax = self.spec_ax
        # Children
        self.topbar = HBox( [self.posview, self.negview, self.cursors])
        self.controlbar = VBox([self.reset, self.scale, self.savepdf, self.done])
        self.middlebar = HBox( [  self.controlbar, self.fig.canvas ] )

        self.children = [  VBox([self.topbar, self.middlebar ]) ]
        self.set_on_redraw()
        self.disp(new=True)
    def on_reset(self, b):
        self.scale.value = 1.0
        self.ax.set_ybound( (self.data.axis1.itoc(0),self.data.axis1.itoc(self.data.size1)) )
        self.ax.set_xbound( (self.data.axis2.itoc(0),self.data.axis2.itoc(self.data.size2)) )
    def disp(self,new=False):
        if new:
            self.proj2.display(figure=self.top_ax, title=self.title)
            xb = self.top_ax.get_xbound()
            sidedisplay(self.proj1, self.side_ax)
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
            self.data.display(scale=self.scale.value, new_fig=False, figure=self.spec_ax, mpldic={'cmap':'winter'})
        if self.negview.value:
            self.data.display(scale=-self.scale.value, new_fig=False, figure=self.spec_ax, mpldic={'cmap':'YlOrRd'})
        self.spec_ax.set_xbound(xb)
        self.spec_ax.set_ybound(yb)
        self.fig.canvas.header_visible = False
        for s in ["left", "top", "right"]:
            self.top_ax.spines[s].set_visible(False)
        self.top_ax.yaxis.set_visible(False)
        for s in [ "top", "right", "bottom"]:
            self.side_ax.spines[s].set_visible(False)
        self.side_ax.xaxis.set_visible(False)


def sidedisplay(dt1d, ax):
    step = dt1d.axis1.itype+1
    dataxis = dt1d.axis1.itoc( dt1d.axis1.points_axis() )
    ax.plot(dt1d.buffer[::step],dataxis[::step])

class Phaser2D(Show2D):
    """
    An interactive phaser in 2D NMR

        Phaser2D(spec)

    """
    def __init__(self, data):
        if data.itype != 3:
            print('Dataset should be complex along both axes, Phasing is not possible')
            return
        super().__init__(data)
        self.data_ref = data
        # print('WARNING this tool is not functional/tested yet')
        # create additional widgets
        slidersize = Layout(width='500px')
        self.F1p0 = widgets.FloatSlider(min=-180, max=180, step=1.0, description='P0',continuous_update=HEAVY, layout=slidersize)
        self.F1p1 = widgets.FloatSlider(min=-250, max=250, step=1.0, description='P1',continuous_update=HEAVY, layout=slidersize)
        self.F2p0 = widgets.FloatSlider(min=-180, max=180, step=1.0, description='P0',continuous_update=HEAVY, layout=slidersize)
        self.F2p1 = widgets.FloatSlider(min=-250, max=250, step=1.0, description='P1',continuous_update=HEAVY, layout=slidersize)
        pivotsize = Layout(width='200px')
        self.pivotF1 = widgets.BoundedFloatText(description='Pivot',
                value=round(self.data.axis1.itoc(0.5*self.data.size1),2), 
                min=self.data.axis1.itoc(self.data.size1),
                max=self.data.axis1.itoc(0),
                format='%.2f',
                layout=pivotsize,
                step=0.1)
        self.pivotF2 = widgets.BoundedFloatText(description='Pivot',
                value=round(self.data.axis2.itoc(0.5*self.data.size2),2), 
                min=self.data.axis2.itoc(self.data.size2),
                max=self.data.axis2.itoc(0),
                format='%.2f',
                layout=pivotsize,
                step=0.1)
        # modify defaults
        self.negview.value = True
        self.done.description = 'Apply'
        for w in [self.F1p0, self.F1p1, self.F2p0, self.F2p1, self.pivotF1, self.pivotF2]:
            w.observe(self.ob)
        self.cancel = widgets.Button(description="Cancel",button_style='warning', layout=self.blay)
        self.cancel.on_click(self.on_cancel)

        stcenter = "<b>F%d</b>"
        box_layout = widgets.Layout(display='flex',
                flex_flow='column',
                align_items='center',
                grid_gap="100px",
                width='100%')
        grid_layout = widgets.Layout(grid_template_columns="40% 40%", 
                    justify_items='center')
        self.phasebar = \
            widgets.GridBox( [widgets.HTML(stcenter%1),     widgets.HTML(stcenter%2), 
                             self.F1p0,                             self.F2p0, 
                             self.F1p1,                             self.F2p1,
                             self.pivotF1,                          self.pivotF2],
                             layout=grid_layout)
        doc = widgets.HTML("""
            <p>Use the sliders to adjust the phase parameters, &nbsp; the pivot can be set with a right click on the spectrum<br>
            Top and Side spectra are taken at the pivot level.<br>
            </p>
            """)

        self.topbar = HBox( [self.posview, self.negview])
        self.controlbar = VBox([self.reset, self.scale, self.savepdf, self.done, self.cancel])
        self.middlebar = HBox( [  self.controlbar, self.fig.canvas ] )

        self.children = [  VBox([self.phasebar, doc, self.topbar, self.middlebar ]) ]
        self.pivotF1.observe(self.on_movepivot)
        self.pivotF2.observe(self.on_movepivot)
        # add right-click event on spectral window
        def on_press(event):
            if event.button == 3:
                self.pivotF1.value = round(event.ydata,4)
                self.pivotF2.value = round(event.xdata,4)
        cids = self.fig.canvas.mpl_connect('button_press_event', on_press)
    # pivot handling
    def ppivot(self):
        "converts from pivot values to centered ones"
        pp1 = 1.0-(self.data.axis1.ctoi(self.pivotF1.value)/self.data.size1)
        pp2 = 1.0-(self.data.axis2.ctoi(self.pivotF2.value)/self.data.size2)
        return (self.F1p0.value + (pp1-0.5)*self.F1p1.value, self.F1p1.value, \
                self.F2p0.value + (pp2-0.5)*self.F2p1.value, self.F2p1.value)
    def ctopivot(self, F1p0, F1p1, F2p0, F2p1):
        "convert from centered to pivot values"
        pp1 = 1.0-(self.data.axis1.ctoi(self.pivotF1.value)/self.data.size1)
        pp2 = 1.0-(self.data.axis2.ctoi(self.pivotF2.value)/self.data.size2)
        return F1p0- (pp1-0.5)*F1p1, F1p1, F2p0- (pp2-0.5)*F2p1, F2p1
    def on_movepivot(self, event):
        if event['name']=='value':
            self.F1p0.value, self.F1p1.value, self.F2p0.value, self.F2p1.value = \
                self.ctopivot(self.lF1p0, self.lF1p1, self.lF2p0, self.lF2p1)
            self.phase()
            #self.disp()

    def ob(self, event):
        "observe changes and start phasing"
        if event['name']=='value':
            self.phase()
    # def close(self):
    #     for w in [self.F1p0, self.F1p1, self.F2p0, self.F2p1, self.scale, self.done, self.cancel]:
    #         w.close()
    def on_cancel(self, b):
        print("No action")
        self.ax.clear()
        self.data.display(figure=self.ax,scale=self.scale.value)
        self.ax.set_xlim(xmin=self.data.axis2.itop(0), xmax=self.data.axis2.itop(self.data.size2))
        self.ax.set_ylim(ymin=self.data.axis1.itop(0), ymax=self.data.axis1.itop(self.data.size1))
        self.close()
    def on_Apply(self, b):
        self.lF1p0, self.lF1p1, self.lF2p0, self.lF2p1 = self.ppivot()         # get centered values
        self.data.phase(self.lF2p0, self.lF2p1, axis='F2').phase(self.lF1p0, self.lF1p1, axis='F1')
        self.data.display(figure=self.ax,scale=self.scale.value)
        self.ax.set_xlim(xmin=self.data.axis2.itop(0), xmax=self.data.axis2.itop(self.data.size2))
        self.ax.set_ylim(ymin=self.data.axis1.itop(0), ymax=self.data.axis1.itop(self.data.size1))
        self.close()
        print("Applied: phase(%.1f,%.1f,axis='F1').phase(%.1f,%.1f,axis='F2')"%(self.lF1p0, self.lF1p1, self.lF2p0, self.lF2p1))
        
    def disp(self,todisplay=None, new=False):
        "display either the current data or the one provided - red and blue"
        if not todisplay:
            todisplay = self.data
        if not new:
            yb = self.spec_ax.get_ybound()
            xb = self.spec_ax.get_xbound()
            for axx in (self.spec_ax, self.top_ax, self.side_ax):
                axx.clear()
            icol = todisplay.axis2.ctoi(self.pivotF2.value)
            irow = todisplay.axis1.ctoi(self.pivotF1.value)
            todisplay.row(irow).display(figure=self.top_ax, title=self.title)
            sidedisplay( todisplay.col(icol), self.side_ax) 
        if self.posview.value:
            todisplay.display(scale=self.scale.value, new_fig=False, figure=self.spec_ax,color='blue')
        if self.negview.value:
            todisplay.display(scale=-self.scale.value, new_fig=False, figure=self.spec_ax, color='red')
        if new:
            yb = self.spec_ax.get_ybound()
            xb = self.spec_ax.get_xbound()
        else:
            self.spec_ax.set_xbound(xb)
            self.spec_ax.set_ybound(yb)
            self.spec_ax.scatter(self.pivotF2.value, self.pivotF1.value, s=200, c='r', alpha=0.5)
            self.spec_ax.plot([self.pivotF2.value, self.pivotF2.value], self.spec_ax.get_ybound(), 'r--', alpha=0.5)
            self.spec_ax.plot(self.spec_ax.get_xbound(), [self.pivotF1.value, self.pivotF1.value], 'r--', alpha=0.5)
        self.fig.canvas.header_visible = False
        self.fig.canvas.header_visible = False
        for s in ["left", "top", "right"]:
            self.top_ax.spines[s].set_visible(False)
        self.top_ax.yaxis.set_visible(False)
        for s in [ "top", "right", "bottom"]:
            self.side_ax.spines[s].set_visible(False)
        self.side_ax.xaxis.set_visible(False)
        # self.ax.set_xlim(xmin=self.data.axis2.itop(0), xmax=self.data.axis2.itop(self.data.size2))
        # self.ax.set_ylim(ymin=self.data.axis1.itop(0), ymax=self.data.axis1.itop(self.data.size1))

    def phase(self):
        "compute phase and display"
        self.lF1p0, self.lF1p1, self.lF2p0, self.lF2p1 = self.ppivot()         # get centered values
        dp = self.data.copy().phase(self.lF2p0, self.lF2p1, axis='F2').phase(self.lF1p0, self.lF1p1, axis='F1')
        self.disp(dp)
    # def phase(self, scale, F1p0, F1p1, F2p0, F2p1):
    #     self.data.copy().phase(F1p0,F1p1,axis='F1').phase(F2p0,F2p1,axis='F2').display(scale=scale);



# if __name__ == '__main__':
#    unittest.main()    
