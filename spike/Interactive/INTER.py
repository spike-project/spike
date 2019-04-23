#!/usr/bin/env python 
# encoding: utf-8

"""
A set of utilities to use spike in NMR or FTMS within jupyter


First version MAD june 2017
preliminary and not fully tested !

"""

from __future__ import print_function, division
import unittest
import sys
import os
import os.path as op
import glob

import matplotlib.pylab as plt
from ipywidgets import interact, interactive, fixed, interact_manual, Layout, HBox, VBox, Label, Output
import ipywidgets as widgets
from IPython.display import display

from ..File.BrukerNMR import Import_1D

# REACTIVE modify callback behaviour
# True is good for inline mode / False is better for notebook mode
REACTIVE = True

class FileChooser:
    """a simple file chooser for Jupyter"""
    def __init__(self, base=None, filetype=['fid','ser'], mode='r', show=True):
        if base is None:
            self.curdir = "/"
        else:
            self.curdir = base
        self.filetype = filetype
        self.mode = mode
        self.wfile = widgets.Text(layout=Layout(width='70%'),description='File to load')
        self.ldir = widgets.Label(value="Chosen dir:  "+self.curdir)
        self.wdir = widgets.Select(
            options=self.dirlist(),
            description='Choose Dir',
            layout=Layout(width='70%'))
        if mode=='r':
            self.wchooser = widgets.Select(
                options=self.filelist(),
                description='Choose File',
                layout=Layout(width='70%'))
            self.wchooser.observe(self.wob)
            self.wfile.disabled=True   # not mofiable in read mode
        elif mode=="w":
            self.wfile.description = 'File to create'
            self.wfile.disabled=False
            self.wchooser = widgets.Select(
                options=self.filelist(),
                description='Files',
                layout=Layout(width='70%'))
        else:
            raise Exception('Only "r" and  "w" modes supported')
        self.wsetdir = widgets.Button(description='❯',layout=Layout(width='20%'),
                button_style='success', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='descend in directory')
        self.wup = widgets.Button(description='❮',layout=Layout(width='20%'),
                button_style='success', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='up to previous directory')
        self.wsetdir.on_click(self.setdir)
        self.wup.on_click(self.updir)
        if show: self.show()
    def filelist(self):
        fl = []
        if self.mode == 'r':
            if type(self.filetype) is str:
                fl = glob.glob(op.join(self.curdir,self.filetype))
            elif type(self.filetype) in (tuple, list):
                for f in self.filetype:
                    fl += glob.glob(op.join(self.curdir,f))
            else:
                raise Exception('TypeError, filetype should be either a string or a list')
        else:   # 'w'
            fl = [f for f in glob.glob(op.join(self.curdir,'*')) if op.isfile(f)]
            self.wfile.value = op.join(self.curdir,self.filetype)
        if fl == []:
            fl = [" "]
        return fl
    def dirlist(self):
        base = self.curdir
        return [d for d in glob.glob(op.join(base,'*')) if (op.isdir(d) and not d.endswith('__pycache__'))]
    def wob(self, e):
        self.wfile.value = self.wchooser.value
    def updir(self, e):
        self.curdir = op.dirname(self.curdir)
        self.ldir.value = "Chosen dir:  "+self.curdir
        self.wdir.options = self.dirlist()
        self.wchooser.options = self.filelist()
    def setdir(self, e):
        self.curdir = self.wdir.value
        self.ldir.value = "Chosen dir:  "+self.curdir
        self.wdir.options = self.dirlist()
        self.wchooser.options = self.filelist()
    def show(self):
        display(VBox(
            [   self.ldir,
            HBox( [self.wdir, VBox([self.wup,self.wsetdir])] ),    
                self.wchooser,
                self.wfile
            ])
        )
        return self
    @property
    def file(self):
        "the chosen complete filename"
        return self.wfile.value
    @property
    def dirname(self):
        "the final dirname containing the chosen file"
        if op.isdir(self.wfile.value):
            return op.basename(self.wfile.value)
        else:
            return op.basename(op.dirname(self.wfile.value))
    @property
    def nmrname(self):
        "the final dirname containing the chosen file for TopSpin files"
        if op.isdir(self.wfile.value):
            return op.join(
                op.basename(op.dirname(self.wfile.value)), self.dirname)
        else:
            return op.join(
                op.basename(op.dirname(op.dirname(self.wfile.value))), self.dirname)
    @property
    def basename(self):
        "the basename of the chosen file"
        return op.basename(self.wfile.value)
    

class Show(object):
    """
    An interactive display, 1D or 2D NMR
        Show(spectrum)
    """
    def __init__(self, data):
        self.data = data
        self.done = widgets.Button(description="Done")
        self.scale = widgets.FloatSlider(description='scale:', min=1, max=100, step=0.5,
                            layout=Layout(width='30%'), continuous_update=REACTIVE)
        self.zoom = widgets.FloatRangeSlider(value=[0, 100],
            min=0, max=100.0, step=0.1, layout=Layout(width='60%'), description='zoom (%):',
            continuous_update=REACTIVE, readout=True, readout_format='.1f')
        self.scale.observe(self.ob)
        self.zoom.observe(self.ob)
        fi,ax = plt.subplots()
        self.ax = ax
        self.data = data
        self.data.display(figure=self.ax)
        self.done.on_click(self.on_done)
        if data.dim == 1:
            display( VBox([HBox([self.scale,self.done]),  self.zoom]) )
        else:
            display( VBox([self.scale, self.zoom]) )
    def on_done(self, b):
        for w in [self.scale, self.done, self.zoom]:
            w.close()
    def ob(self, event):
        "observe events and display"
        if event['name']=='value':
            self.ax.clear()
            self.data.display(scale=self.scale.value, new_fig=False, figure=self.ax)
            z = self.zoom.value
            self.ax.set_xlim(left=self.data.axis1.itop(z[0]*self.data.size1/100), right=self.data.axis1.itop(z[1]*self.data.size1/100) )

class Phaser1D(object):
    """
    An interactive phaser in 1D NMR

        Phaser1D(spectrum)

    best when in %matplotlib inline

    """
    def __init__(self, data):
        self.scale = widgets.FloatSlider(description='scale:', min=1, max=100, step=0.5,
                            layout=Layout(width='30%'), continuous_update=REACTIVE)
        self.p0 = widgets.FloatSlider(description='P0:',min=-180, max=180, step=0.1,
                            layout=Layout(width='100%'), continuous_update=REACTIVE)
        self.p1 = widgets.FloatSlider(description='P1:',min=-360, max=360, step=1.0,
                            layout=Layout(width='100%'), continuous_update=REACTIVE)
        self.pivot = widgets.FloatSlider(description='pivot:',max=data.axis1.itop(0), min=data.axis1.itop(data.size1),
                            layout=Layout(width='30%'),  continuous_update=REACTIVE)
        self.pivot.value = (data.axis1.itop(data.size1/2))
        self.zoom = widgets.FloatRangeSlider(value=[0, 100],
            min=0, max=100.0, step=0.1, layout=Layout(width='60%'), description='zoom (%):',
            continuous_update=REACTIVE, readout=True, readout_format='.1f',)
        self.scale.observe(self.ob)
        self.zoom.observe(self.ob)
        self.p0.observe(self.ob)
        self.p1.observe(self.ob)
        self.pivot.observe(self.ob)
        self.button = widgets.Button(description="Apply correction")
        self.button.on_click(self.on_Apply)
        self.cancel = widgets.Button(description="Cancel")
        self.cancel.on_click(self.on_cancel)
        display(HBox([self.button, self.cancel]))
#       interact( self.phase, scale=self.scale, p0=self.p0, p1=self.p1, pivot=self.pivot)
        display( VBox([HBox([self.scale, self.pivot]), self.zoom, self.p0, self.p1]) )
        fi,ax = plt.subplots()
        self.ax = ax
        self.data = data
        self.data.display(figure=self.ax)

    def on_cancel(self, b):
        print("No action")
#        self.data.display(scale=self.scale.value);
        for w in [self.p0, self.p1, self.scale, self.button, self.cancel, self.zoom, self.pivot]:
            w.close()
        self.ax.clear()
        self.data.display(new_fig=False, figure=self.ax)
    def on_Apply(self, b):
        lp0, lp1 = self.ppivot() # get centered values
        print("Applied: phase(%.1f,  %.1f)"%(lp0, lp1))
        for w in [self.p0, self.p1, self.scale, self.button, self.cancel, self.zoom, self.pivot]:
            w.close()
        self.ax.clear()
        self.data.phase(lp0, lp1).display(new_fig=False, figure=self.ax)
    def ppivot(self):
        "converts from pivot values to centered ones"
        pp = (self.pivot.value - self.data.axis1.itop(self.data.size1))/  \
             (self.data.axis1.itop(0)-self.data.axis1.itop(self.data.size1)) # this one is in 0.0 ... 1.0
        return (self.p0.value + (pp-0.5)*self.p1.value, self.p1.value)

    def ob(self, event):
        "observe changes and start phasing"
        if event['name']=='value':
            self.phase()

    def phase(self):
        "apply phase and display"
        self.ax.clear()
        lp0, lp1 = self.ppivot() # get centered values
        self.data.copy().phase(lp0, lp1).display(scale=self.scale.value, new_fig=False, figure=self.ax)
        plt.plot([self.pivot.value,self.pivot.value],[0,self.data.absmax/self.scale.value])
        z = self.zoom.value
        self.ax.set_xlim(left=self.data.axis1.itop(z[0]*self.data.size1/100), right=self.data.axis1.itop(z[1]*self.data.size1/100) )

class Phaser2D(object):
    """
    An interactive phaser in 2D NMR

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

class AvProc1D:
    "Detailed 1D NMR Processing"
    def __init__(self, filename=""):
        self.wfile = widgets.Text(description='File to process',layout=Layout(width='80%'), value=filename)
        self.wapod = widgets.Dropdown(
            options=['None', 'apod_sin (sine bell)', 'apod_em (Exponential)', 'apod_gm (Gaussian)', 'gaussenh (Gaussian Enhacement)', 'kaiser'],
            value='apod_sin (sine bell)',
            description='Apodisation')
        self.wpapod_Hz = widgets.FloatText(
            value=1.0,
            min=0, # max exponent of base
            max=30, # min exponent of base
            description='Width in Hz',
            layout=Layout(width='15%'),
            disabled = True)
        self.wpapod_enh = widgets.FloatText(
            value=2.0,
            min=0.0, # max exponent of base
            max=5.0, # min exponent of base
            description='strength',
            layout=Layout(width='15%'),
            step=1,
            disabled = True)
        self.wpapod_sin = widgets.FloatText(
            value=0.0,
            min=0, # max exponent of base
            max=0.5, # min exponent of base
            description='bell shape',
            layout=Layout(width='15%'),
            step=0.01,
            tooltip='value is the maximum of the bell, 0 is pure cosine, 0.5 is pure sine',
            disabled = False)
        self.wzf = widgets.Dropdown(
            options=[0, 1, 2, 4, 8],
            value=1,
            description='Zero-Filling')
        self.wphase0 = widgets.FloatText(
            value=0, description='Phase : P0', layout=Layout(width='20%'), disabled = True)
        self.wphase1 = widgets.FloatText(
            value=0, description='P1', layout=Layout(width='20%'), disabled = True)
        self.wapmin = widgets.Checkbox(
            value=True, description='AutoPhasing', tooltip='Perform AutoPhasing')
        self.wapmin.observe(self.apmin_select)
        self.wbcorr = widgets.Checkbox(
            value=False, description='Baseline Correction', tooltip='Perform AutoPhasing')
        self.wapod.observe(self.apod_select)
        self.bapod = widgets.Button(description='Show effect on FID')
        self.bapod.on_click(self.show_apod)
        self.bdoit = widgets.Button(description='Process')
        self.bdoit.on_click(self.process)
        self.show()
        fi,ax = plt.subplots()
        self.ax = ax
        if os.path.exists(filename):
            self.load()
            #self.data.set_unit('sec')
            self.display()
    def apod_select(self, e):
        test = self.wapod.value.split()[0]
        self.wpapod_sin.disabled = True
        self.wpapod_Hz.disabled = True
        self.wpapod_enh.disabled = True
        if test == "apod_sin":
            self.wpapod_sin.disabled = False
        if test in ('apod_em', 'apod_gm','gaussenh'):
            self.wpapod_Hz.disabled = False
        if test == 'gaussenh':
            self.wpapod_enh.disabled = False
    def apmin_select(self, e):
        for w in self.wphase0, self.wphase1:
            w.disabled = self.wapmin.value
    def load(self):
        self.data = Import_1D(self.wfile.value)
    def apod(self):
        func = self.wapod.value.split()[0]
        todo = None
        if func == 'apod_sin':
            todo = 'self.data.apod_sin(%f)'%(self.wpapod_sin.value,)
        elif func in ('apod_em', 'apod_gm'):
            todo = 'self.data.%s(%f)'%(func, self.wpapod_Hz.value)
        elif func == 'gaussenh':
            todo = 'self.data.gaussenh(%f,enhancement=%f)'%(self.wpapod_Hz.value, self.wpapod_enh.value)
        if todo is not None:
            eval(todo)
        return self.data
    def show_apod(self, e):
        self.load()
        self.apod()
        self.display()
    def process(self, e):
        self.load()
        self.apod().zf(self.wzf.value).ft_sim().bk_corr().set_unit('ppm')
        if self.wapmin.value:
            self.data.apmin()
            self.wphase0.value = round(self.data.axis1.P0,1)
            self.wphase1.value = self.data.axis1.P1
        else:
            self.data.phase(self.wphase0.value, self.wphase1.value)
        self.display()
    def display(self):
        self.ax.clear()
        self.data.display(figure=self.ax)
    def show(self):
        display(
            VBox([self.wfile,
                HBox([self.wapod, self.wpapod_sin, self.wpapod_Hz, self.wpapod_enh, self.bapod]),
                self.wzf,
                HBox([self.wapmin, self.wphase0, self.wphase1]),
#                self.wbcorr,
                self.bdoit]) )


#if __name__ == '__main__':
#    unittest.main()    