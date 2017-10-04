#!/usr/bin/env python 
# encoding: utf-8

"""
A set of utilities to use within jupyter


First version MAD june 2017
preliminary and not fully tested !

"""

from __future__ import print_function, division
import unittest
import sys

from ipywidgets import interact, interactive, fixed, interact_manual, Layout, HBox, VBox
import ipywidgets as widgets
from IPython.display import display


# REACTIVE modify callback behaviour
# True is good for inline mode / False is better for notebook mode
REACTIVE = True

class Show(object):
    """
    An interactive display, 1D or 2D
        Show(spec)

    best when in %matplotlib inline

    """
    def __init__(self, data):
        self.data = data
        self.scale = widgets.FloatSlider(min=1, max=100, step=0.5,
                            layout=Layout(width='100%'), continuous_update=REACTIVE)
        if data.dim == 1:
            a,b = data.axis1.itoc(0), data.axis1.itoc(data.axis1.size-1)
            self.z1 = widgets.FloatSlider(min=min(a,b),
                            max=max(a,b), # whole window,
                            value=min(a,b),
                            step=0.1,
                            layout=Layout(width='100%'),
                            continuous_update=REACTIVE)
            self.z2 = widgets.FloatSlider(min=min(a,b),
                            max=max(a,b), # whole window,
                            value=max(a,b),
                            step=0.1,
                            layout=Layout(width='100%'),
                            continuous_update=REACTIVE)
        else:
            self.z1 = widgets.Label('No zoom in 2D yet')
            self.z2 = self.z1
        self.done = widgets.Button(description="Done")
        self.done.on_click(self.on_done)
        interact(self.show, scale=self.scale, z1=self.z1, z2=self.z2)
        display(self.done)
    def on_done(self, b):
        for w in [self.scale, self.done, self.z1, self.z2]:
            w.close()
    def show(self, scale, z1=0, z2=10):
        if self.data.dim == 1:
            self.data.display(scale=scale, zoom=[z1,z2]);
        else:
            self.data.display(scale=scale);

class Phaser1D(object):
    """
    An interactive phaser in 1D

        Phaser1D(spec)

    best when in %matplotlib inline

    """
    def __init__(self, data):
        self.data = data
        self.scale = widgets.FloatSlider(min=1, max=100, step=0.5,
                            layout=Layout(width='100%'), continuous_update=REACTIVE)
        self.p0 = widgets.FloatSlider(min=-180, max=180, step=0.1,
                            layout=Layout(width='100%'), continuous_update=REACTIVE)
        self.p1 = widgets.FloatSlider(min=-250, max=250, step=1.0,
                            layout=Layout(width='100%'), continuous_update=REACTIVE)
        self.button = widgets.Button(description="Apply correction")
        self.button.on_click(self.on_Apply)
        self.cancel = widgets.Button(description="Cancel")
        self.cancel.on_click(self.on_cancel)
        interact(self.phase, scale=self.scale, p0=self.p0, p1=self.p1)
        display(HBox([self.button, self.cancel]))
    def on_cancel(self, b):
        print("No action")
#        self.data.display(scale=self.scale.value);
        for w in [self.p0, self.p1, self.scale, self.button, self.cancel]:
            w.close()
    def on_Apply(self, b):
        print("Applied: phase(%.1f,  %.1f)"%(self.p0.value, self.p1.value))
        for w in [self.p0, self.p1, self.scale, self.button, self.cancel]:
            w.close()
        self.data.phase(self.p0.value, self.p1.value)
    def phase(self, scale, p0, p1):
        self.data.copy().phase(p0,p1).display(scale=scale);

class Phaser2D(object):
    """
    An interactive phaser in 2D

        Phaser2D(spec)

    best when in %matplotlib inline

    """
    def __init__(self, data):
        self.data = data
        self.scale = widgets.FloatSlider(min=1, max=100, step=0.5,
                                        layout=Layout(width='50%'), continuous_update=REACTIVE)
        self.F1p0 = widgets.FloatSlider(min=-180, max=180, step=0.1, description='F1 p0',
                                        layout=Layout(width='50%'), continuous_update=REACTIVE)
        self.F1p1 = widgets.FloatSlider(min=-250, max=250, step=1.0, description='F1 p1',
                                        layout=Layout(width='50%'), continuous_update=REACTIVE)
        self.F2p0 = widgets.FloatSlider(min=-180, max=180, step=0.1, description='F2 p0',
                                        layout=Layout(width='50%'), continuous_update=REACTIVE)
        self.F2p1 = widgets.FloatSlider(min=-250, max=250, step=1.0, description='F2 p1',
                                        layout=Layout(width='50%'), continuous_update=REACTIVE)
        self.button = widgets.Button(description="Apply correction")
        self.button.on_click(self.on_Apply)
        self.cancel = widgets.Button(description="Cancel")
        self.cancel.on_click(self.on_cancel)
        interact(self.phase, scale=self.scale, F1p0=self.F1p0, F1p1=self.F1p1, F2p0=self.F2p0, F2p1=self.F2p1)
        display(HBox([self.button, self.cancel]))
    def on_cancel(self, b):
        print("No action")
#        self.data.display(scale=self.scale.value);
        for w in [self.F1p0, self.F1p1, self.F2p0, self.F2p1, self.scale, self.button, self.cancel]:
            w.close()
    def on_Apply(self, b):
        print("Applied: phase(%.1f,%.1f,axis='F1').phase(%.1f,%.1f,axis='F')"%(self.F1p0.value, self.F1p1.value, self.F2p0.value, self.F2p1.value))
        for w in [self.F1p0, self.F1p1, self.F2p0, self.F2p1, self.scale, self.button, self.cancel]:
            w.close()
        self.data.phase(self.F1p0.value, self.F1p1.value, axis='F1').phase(self.F2p0.value, self.F2p1.value, axis='F2')
    def phase(self, scale, F1p0, F1p1, F2p0, F2p1):
        self.data.copy().phase(F1p0,F1p1,axis='F1').phase(F2p0,F2p1,axis='F2').display(scale=scale);

#if __name__ == '__main__':
#    unittest.main()