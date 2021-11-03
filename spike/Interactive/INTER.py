#!/usr/bin/env python 
# encoding: utf-8

"""
A set of utilities to use spike in NMR or FTMS within jupyter


First version MAD june 2017
Intermediate version MAD october 2019
Improvements MAD spring - summer 2021
"""

from __future__ import print_function, division
import unittest
import sys
import os
import os.path as op
from pathlib import Path 
import glob
import math as m
import time
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import MultiCursor
from ipywidgets import fixed, Layout, HBox, VBox, Label, Output, Button, Tab
import ipywidgets as widgets
from IPython.display import display, HTML, Javascript, Markdown, Image
import numpy as np

from .. import version
if version.revision < "533":
    warnmsg = """
There is version missmatch between the core program and the interactive tools
You may experiment some difficulties or halting with this notebook
"""
    warnings.warn(warnmsg)

from .. import NMR
from ..NMR import NMRData

from ..File.BrukerNMR import Import_1D

from .ipyfilechooser import FileChooser as FileChooser_code
try:
    import spike.plugins.bcorr as bcorr
except ModuleNotFoundError:
    print('Baseline correction plugins not installed !')

__version__ = "1.3.0"

# REACTIVE modify callback behaviour
# True is good for inline mode / False is better for notebook mode
REACTIVE = True
HEAVY = False
# Activate_Wheel: scale with wheel control in the graphic cells 
Activate_Wheel = True
# reverse_scroll : if True, reverses direction of mouse scroll
reverse_scroll = True

# default color set
Colors = ('black','red','blue','green','orange',
'blueviolet','crimson','turquoise','indigo',
'magenta','gold','pink','purple','salmon','darkblue','sienna')

##############################
# First General Utilities
##############################
def initialize(verb=1):
    "initialize notebook interface for Spike"
    if verb>0:
        print("\nInteractive module version,",__version__)
    if verb>1:
        hidecode()
        hidedoc()
    else:
        hidecode(message="")
        hidedoc(message="")
    initCSS()
    Logo()

def initCSS():
    "adapt some CSS"
    display(HTML('''
<style>
hr         {height: 2px; border: 0;border-top: 1px solid #ccc;margin: 1em 0;padding: 0; }
.container {width:100% !important; }
</style>
    '''))

def hidecode(initial='show', message="<i>useful to show/print a clean screen when processing is finished</i>"):
    """
    this func adds a button to hide/show the code on a jupyter page
    initial is either 'show' or 'hide'
    inspired from: https://stackoverflow.com/questions/27934885/how-to-hide-code-from-cells-in-ipython-notebook-visualized-with-nbviewer/28073228#28073228
    """
    if initial == 'show':
        init = 'false'
    elif initial == 'hide':
        init = 'true'
    display(HTML('''
<script>
code_show=%s; 
function code_toggle()
    { if (code_show)
        { $('div.input').hide(); $('#but').val("show python code");
        } else { $('div.input').show(); $('#but').val("hide python code");
    }
    code_show = !code_show } 
$(document).ready(code_toggle);
</script>
<form action="javascript:code_toggle()">
<input id="but" type="submit" style="border:1px solid black; background-color:#DDD">
%s
</form>'''%(init, message)))

def hidedoc(initial='show', message="<i>useful to show/print a clean screen when processing is finished</i>"):
    """
    this func adds a button to hide/show the doc on a jupyter page
    initial is either 'show' or 'hide'
    inspired from: https://stackoverflow.com/questions/27934885/how-to-hide-code-from-cells-in-ipython-notebook-visualized-with-nbviewer/28073228#28073228
    """
    if initial == 'show':
        init = 'false'
    elif initial == 'hide':
        init = 'true'
    display(HTML('''
<script>
doc_show=%s; 
function doc_toggle()
    { if (doc_show)
        { $('div.text_cell').hide(); $('#butdoc').val("show documentation");
        } else { $('div.text_cell').show(); $('#butdoc').val("hide documentation");
    }
    doc_show = !doc_show } 
$(document).ready(doc_toggle);
</script>
<form action="javascript:doc_toggle()">
<input id="butdoc" type="submit" style="border:1px solid black; background-color:#DDD">
%s
</form>'''%(init, message)))

from . import __path__
def UserLogofile():
    "return the address on disk of the User Logo file - located in $HOME/Spike/Logo.png"
    ufile = Path.home()/'Spike'/'Logo.png'
    if ufile.exists():
        return ufile
    else:
        return Logofile()
def Logofile():
    "return the address on disk of the Logo file"
    return op.join(__path__[0],'Logo.png')
def Logo(width=150):
    "display the Logo"
    display(Image(filename=Logofile(),width=width))
def jsalert(msg="Alert text"):
    "send a javascript alert"
    # use \\n if you want to add several lines
    display(Javascript("alert('%s')"%msg))
class _FileChooser:
    """a simple file chooser for Jupyter - obsolete - not used any more"""
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
        self.wsetdir = Button(description='❯',layout=Layout(width='20%'),
                button_style='success', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='descend in directory')
        self.wup = Button(description='❮',layout=Layout(width='20%'),
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
import asyncio

class Timer:
    def __init__(self, timeout, callback):
        self._timeout = timeout
        self._callback = callback

    async def _job(self):
        await asyncio.sleep(self._timeout)
        self._callback()

    def start(self):
        self._task = asyncio.ensure_future(self._job())

    def cancel(self):
        self._task.cancel()

def debounce(wait):
    """ Decorator that will postpone a function's
        execution until after `wait` seconds
        have elapsed since the last time it was invoked. """
    def decorator(fn):
        timer = None
        def debounced(*args, **kwargs):
            nonlocal timer
            def call_it():
                fn(*args, **kwargs)
            if timer is not None:
                timer.cancel()
            timer = Timer(wait, call_it)
            timer.start()
        return debounced
    return decorator

##############################
# Then 1D NMR Tools
##############################

class Show1D(HBox):
    """
    An interactive display, 1D NMR
        Show1D(spectrum)
    to be developped for peaks and integrals
    """
    def __init__(self, data, title=None, figsize=None, show=True, create_children=True):
        """
        data : to be displayed
        title : text for window
        figsize : size in inches (x,y)
        """
        super().__init__()
        self.data = data
        self.title = title
        # graphic set-up details
        self.yaxis_visible = False
        self.ax_spines = False
        self.canvas_header_visible = False
        if figsize is None:
            figsize = mpl.rcParams['figure.figsize'] #(10,5)
        if reverse_scroll:
            self.reverse_scroll = -1
        else:
            self.reverse_scroll = 1
        self.blay = Layout(width='80px')  # default layout for buttons
        self.blank = widgets.HTML("&nbsp;",layout=self.blay)
        self.done = Button(description="Done", button_style='success',layout=self.blay,
            tooltip="Stop interactivity and store figure")
        self.done.on_click(self.on_done)
        self.reset = Button(description="Reset", button_style='success',layout=self.blay,
            tooltip="Reset to full spectrum")
        self.savepdf = widgets.Button(description="Save figure", button_style='',layout=self.blay,
                tooltip='Save the current figure as pdf')
        def onsave(e):
            name = self.fig.get_axes()[0].get_title()
            name = name.replace('/','_')+'.pdf'
            if name.startswith('.'):
                name = 'Figure'+name
            self.fig.savefig(name)
            print('figure saved as: ',name)
        self.savepdf.on_click(onsave)
        self.reset.on_click(self.on_reset)
        self.scale = widgets.FloatSlider(description='scale:', value=1.0, min=0.1, max=200, step=0.1,
                            layout=Layout(width='80px', height=str(0.8*2.54*figsize[1])+'cm'), continuous_update=REACTIVE,
                            orientation='vertical')
        for widg in (self.scale,):
            widg.observe(self.obdisp)
        if create_children:
            plt.ioff()
            fi,ax = plt.subplots(figsize=figsize)
            fi.set_tight_layout(True)
            plt.ion()
            self.ax = ax
            self.fig = fi
#            self.xb = self.ax.get_xbound()
            self.children = [  VBox([self.reset, self.scale, self.savepdf, self.done]), self.fig.canvas ]
            if show:
                self.draw()
#                self.disp(redraw=True)
    # call backs
    def on_done(self, b):
        self.close()
        display(self.fig)   # shows spectrum
    def on_reset(self, b):
        self.scale.value = 1.0
        self.ax.set_xbound( (self.data.axis1.itoc(0),self.data.axis1.itoc(self.data.size1)) )
    def ob(self, event):
        "observe events and redraw"
        if event['name']=='value':
            self.draw()
    def obdisp(self, event):
        "observe events and display"
        if event['name']=='value':
            self.disp()
    def scale_up(self, step):
        self.scale.value *= 1.1892**(self.reverse_scroll*step) # 1.1892 is 4th root of 2.0
    def set_on_redraw(self):
        def on_scroll(event):
            self.scale_up(event.step)
        if Activate_Wheel:
            self.fig.canvas.capture_scroll = True
            cids = self.fig.canvas.mpl_connect('scroll_event', on_scroll)
    def draw(self):
        "builds and display the picture"
        self.ax.clear()
        self.data.display(new_fig=False, figure=self.ax, title=self.title)
        try:
            self.ax.set_xbound(self.xb)
            self.ax.set_ybound(self.yb0/self.scale.value)
        except AttributeError:               # the very first time
            self.xb = self.ax.get_xbound()
            self.yb0 = np.array(self.ax.get_ybound())
        self.fig.canvas.header_visible = False
        for s in ["left", "top", "right"]:
            self.ax.spines[s].set_visible(False)
        self.ax.yaxis.set_visible(False)
        visible_ticks = {   "top": False,   "right": False}
        self.ax.tick_params(axis="x", which="both", **visible_ticks)
        self.set_on_redraw()
    def disp(self):
        self.ax.set_ybound(self.yb0/self.scale.value)
class baseline1D(Show1D):
    def __init__(self, data, figsize=None, show=True):
        super().__init__( data, figsize=figsize, show=False)
        self.data.real()
        a, b = self.itoc3(0.0), self.itoc3(self.data.size1)   # spectral borders
        self.select = widgets.BoundedFloatText(description='select:',
                        min=min(a,b), max=max(a,b),
                        value=self.itoc3(0.5*self.data.size1), continuous_update=REACTIVE,
                        tooltip='position of the selected point in %s'%(data.unit[0]))
        self.smooth = widgets.IntSlider(description='smooth:', min=0, max=20,  layout=Layout(width='70%'),
                                        tooltip='apply a local smoothing for pivot points',
                                         value=1, readout=True, continuous_update=REACTIVE)
        self.toshow = widgets.Dropdown( options=['baseline', 'corrected', 'points'],  description='Display:')
        for w in [self.select, self.smooth, self.toshow]:
            w.observe(self.obdisp)

        self.bsl_points = []
        self.baseline = np.zeros_like(self.data.buffer)
        self.cancel = widgets.Button(description="Cancel", button_style='warning', layout=self.blay,
            tooltip='Exit without corrections')
        self.cancel.on_click(self.on_cancel)
        self.auto = widgets.Button(description="Auto", button_style='success', layout=self.blay,
            tooltip='Automatic set of points')
        self.auto.on_click(self.on_auto)
        self.set = widgets.Button(description="Add", button_style='success', layout=self.blay,
            tooltip='Add a baseline point at the selector position')
        self.set.on_click(self.on_set)
        self.unset = widgets.Button(description="Rem", button_style='warning', layout=self.blay,
            tooltip='Remove the baseline point closest to the selector')
        self.unset.on_click(self.on_unset)
        orig = self.children
        self.children = [VBox([
            HBox([widgets.HTML('Choose baseline points'), self.auto, self.select, self.set, self.unset]),
            HBox([self.cancel, self.toshow, self.smooth]),
            HBox(orig)])]
        def on_press(event):
            v = event.xdata
            self.select.value = round(v,3)
        cids = self.fig.canvas.mpl_connect('button_press_event', on_press)
        if show:
            self.draw()

    # def show(self):
    #     "create the widget and display the spectrum"
    #     display(self.Box)
    #     self.data.display(figure=self.ax)
    #     self.xb = self.ax.get_xbound()  # initialize zoom
    #     ppos = self.data.axis1.itop(self.select.value)
    #     self.ax.plot([ppos,ppos], self.ax.get_ybound())
    #     self.fig.canvas.header_visible = False
    # def close(self):
    #     "close all widget"
    #     for w in [ self.select, self.auto, self.set, self.unset, self.cancel, self.toshow, self.smooth, self.Box]:
    #         w.close()
    #     super().close()
    def itoc3(self,value):
        return round(self.data.axis1.itoc(value), 3)
    def on_done(self, e):
        if self.bsl_points == []:
            jsalert('Please define control points before applying baseline correction')
            return
        self.close()
        printlist = [round(i,3) for i in sorted(self.bsl_points)]
        print("Applied correction:\ndata.bcorr(xpoints={0}, nsmooth={1})".format(printlist, self.smooth.value))
        #print('Applied correction:\n', sorted(self.bsl_points))
        self.data.set_buffer( self.data.get_buffer() - self.correction() )
        self.selector.set_visible(False)
        self.selector.set_visible(False)
        super().disp()
        display(self.fig)   # shows spectrum
    def on_auto(self, e):
        "automatically set baseline points"
        self.bsl_points = [self.data.axis1.itoc(x) for x in bcorr.autopoints(self.data)]
        self.disp()
    def on_cancel(self, e):
        self.close()
        print('no baseline correction')
    def on_set(self, e):
        "add baseline points at selector"
        self.bsl_points.append( self.select.value )
        self.disp()
    def on_unset(self, e):
        "remove baseline points closest from selector"
        here = self.select.value
        distclose = np.inf
        pclose = np.NaN
        for i,p in enumerate(self.bsl_points):
            if abs(p-here)< distclose:
                pclose = p
                distclose = abs(p-here)
        self.bsl_points.remove(pclose)
        self.disp()
    def smoothed(self):
        "returns a smoothed version of the data"
        from scipy.signal import fftconvolve
        buf = self.data.get_buffer()
        mask = np.array([1,1,1,1,1])
        return fftconvolve(buf, mask, mode='same')
    def correction(self):
        "returns the correction to apply as a numpy array"
        ibsl_points = self.data.axis1.ptoi( np.array(self.bsl_points) ).astype(int)
        x = np.arange(self.data.size1)
        yy = self.data.get_buffer()
        if len(self.bsl_points) == 0 :
            return 0.0
        elif len(self.bsl_points) == 1 :
            value = self.data[ibsl_points[0]] * np.ones( self.data.size1 )
        elif len(self.bsl_points) < 4 :
            corr = bcorr._linear_interpolate(yy, ibsl_points, nsmooth=self.smooth.value)
            value = corr(x)
        else:
            corr = bcorr._spline_interpolate(yy, ibsl_points, nsmooth=self.smooth.value)
            value = corr(x)
        return value
    def corrected(self):
        value = self.data.copy()
        value.set_buffer( value.get_buffer() - self.correction() )
        return value
    def draw(self):
        "used to create the view with additional artists"
        super().draw()
        # baseline
        self.drbaseline = self.ax.plot(self.data.axis1.unit_axis(), self.baseline, color='r')[0]
        self.drbaseline.set_visible(False)
        # corrected spectrum
        self.drspectrum = self.ax.plot(self.data.axis1.unit_axis(), self.data.get_buffer()-self.correction(), color='r' )[0]
        self.drspectrum.set_visible(False)
        # selector
        ppos = self.select.value
        self.selector = self.ax.plot([ppos,ppos], self.ax.get_ybound())[0]   # visible !
        # pivot points
        y = bcorr.get_ypoints(  self.data.get_buffer(), 
                                self.data.axis1.ptoi( np.array(self.bsl_points)),
                                nsmooth=self.smooth.value )
        self.drpivot = self.ax.scatter(self.bsl_points, y,  c='r', marker='o')

    def disp(self):
        "used to refresh view"
        if len(self.bsl_points)>0:
            if self.toshow.value == 'baseline':
#                ( self.data.copy()-self.corrected() ).display(new_fig=False, figure=self.ax, color='r')
                self.drbaseline.set_visible(True)
                self.drspectrum.set_visible(False)
                self.drbaseline.set_ydata( self.correction() )
            elif self.toshow.value == 'corrected':
#                self.corrected().display(new_fig=False, figure=self.ax, color='r', scale=self.scale.value)
                self.drbaseline.set_visible(False)
                self.drspectrum.set_visible(True)
                self.drspectrum.set_ydata( self.data.get_buffer() - self.correction() )
        # selector
        ppos = self.select.value
        self.selector.set_xdata([ppos,ppos])
        self.selector.set_ydata(self.ax.get_ybound())
        # pivots
        y = bcorr.get_ypoints(  self.data.get_buffer(), 
                                self.data.axis1.ptoi( np.array(self.bsl_points)),
                                nsmooth=self.smooth.value )
        del(self.drpivot)
        self.drpivot = self.ax.scatter(self.bsl_points, y,  c='r', marker='o')
        # set zoom
        super().disp()

class SpforSuper():
    "a holder for one spectrum to SuperImpose"
    def __init__(self, i, filename='None', axref=False, ref=None):
        """
        filename is the dataset to show
        axref is the axes of the caller, where to display
        ref is a dataset in case filename is Self 
        """
        self.axref = axref
        self.ref = ref
        self.drartist = None
        self.filename = widgets.Text(value=filename, layout=Layout(width='30em'))
        def nactivate(e):
            if self.filename.value != 'None':
                self.direct.value = 'up'  # will call self.activate()
        self.filename.observe(nactivate)
        j = i%len(Colors)
        self.color = widgets.Dropdown(options=["steelblue"]+list(Colors),value=Colors[j],layout=Layout(width='90px'))
        self.direct = widgets.Dropdown(options=['up','down','off'], value='off', layout=Layout(width='60px'))
        self.direct.observe(self.activate)
        self.zmleft = widgets.FloatText(value=1000,layout=Layout(width='70px'))
        self.zmright = widgets.FloatText(value=-1000,layout=Layout(width='70px'))
        self.label = widgets.Checkbox(value=False, indent=False, layout=Layout(width='40px'))
        self.scale = widgets.FloatText(value=1.0,layout=Layout(width='70px')) 
        self.offset = widgets.FloatText(value=0.0, layout=Layout(width='60px'), tooltip="offset")
        self.splw = widgets.FloatText(value=1.0, step=0.1, layout=Layout(width='50px'))
        # for w in (self.color, self.direct, self.zmleft, self.zmright, self.label, self.scale, self.offset, self.splw):
        #     w.observe(self.ob)
        #control bar
        self.me = HBox([self.B(str(i)),
                        self.filename,
                        self.scale,
                        self.direct,
                        self.offset,
                        self.zmleft,
                        self.zmright,
                        self.splw,
                        self.color,
                        self.label,
                        ])
        self. header = HBox([self.B('Id','40px'),
                            self.B('file name','5em'),
                            self.I('or "Self" for spectrum excerpts', '26em'),
                            self.B('scale','60px'),
                            self.B('direction','60px'),
                            self.B('%offset','50px'),
                            self.B(' ','50px'), self.B('Zoom','90px'),
                            self.B('Linewidth','80px'),
                            self.B('color','70px'),
                            self.B('label')
                            ])
        self.activate(None)
    def ob(self, event):
        "observe events and display"
        if event['name']=='value':
            self.draw()
    def injectobserve(self, method):
        "used to to observe form outer object for draw"
        for w in (self.color,  self.zmleft, self.zmright, self.label, self.splw, self.direct,   self.scale, self.offset):
            w.observe(method)
    def activate(self, e):
        for w in (self.color, self.zmleft, self.zmright, self.label, self.scale, self.offset, self.splw):
            if self.direct.value == 'off':
                w.disabled = True
            else:
                w.disabled = False
        if self.filename.value == 'None':
            self.direct.disabled = True
        else:
            self.direct.disabled = False
    def B(self, text, size='auto'):
        return widgets.HTML(value="<b>%s</b>"%text, layout=widgets.Layout(width=size))
    def I(self, text, size='auto'):
        return widgets.HTML(value="<i>%s</i>"%text, layout=widgets.Layout(width=size))
    def draw(self, unit='ppm'):
        if self.filename.value == 'None' or self.direct.value == 'off':
            return  # nothing to do
        # else
        if self.drartist is not None:   # remove previous ones
            del(self.drartist)
        # label
        if self.label.value:
            lb = self.nmrname
        else:
            lb = None
        # load
        if self.filename.value == 'Self':
            d = self.ref.copy().set_unit(unit)
        else:
            try:
                d = NMRData(name=self.filename.value).set_unit(unit)
            except:
                jsalert('%s : File not found'%self.filename.value)
                return
        # draw
        zml = min(self.zmleft.value, d.axis1.itoc(0))
        zmr = max(self.zmright.value, d.axis1.itoc(d.size1))
        d.display(
            scale=1,
            figure=self.axref,
            zoom=(zml,zmr),
            color=self.color.value,
            linewidth=self.splw.value,
            label=lb)
        self.drartist = d.disp_artist[0]
        self.ydata = self.drartist.get_ydata()   # store spectral data
        # adapt to scale and direction
        self.drartist.set_ydata(  self.move_ydata()  )
        # text
        self.text = self.axref.text(zmr, self.ydata[-1], ' ', color=self.color.value) # will be changed by disp()
        self.disp()
    def move_ydata(self):
        "adapt self.ydata to direct, scale, and offset and return the value"
        # add offset
        off = self.axref.get_ybound()[1] * self.offset.value/100   # relative to screen in spectral unit
        return self.mult*self.ydata + off 
    def disp(self):
        " refresh "
        if self.filename.value == 'None' or self.direct.value == 'off':
            return  # nothing to do
        yd = self.move_ydata()
        self.drartist.set_ydata( yd )
        m = self.mult
        if m != 1:
            tt = "x%.1f"%(m,)
            self.text.set_text(tt)
            self.text.set_y(yd[-1])
        else:
            self.text.set_text(" ")

    @property
    def mult(self):
        if self.direct.value == 'up':
            mult = self.scale.value
        elif self.direct.value == 'down':
            mult = -self.scale.value
        else:
            mult = 1
        return mult
    @property
    def nmrname(self):
        fnm = self.filename.value
        if fnm not in ('None','Self'):
            return op.join(op.basename(op.dirname(fnm)), op.basename(fnm))
        elif fnm =='Self':
            return 'excerpt'
        else:
            return None
    def simpledisplay(self):
        "a simplified widget "
        lb = widgets.Label(self.nmrname)
        return HBox([lb, self.scale])
class Show1Dplus(Show1D):
    def __init__(self, data, base='/DATA', N=9, figsize=None, title=None, show=True):
        super().__init__( data, figsize=figsize, title=title, show=False)
        # spectrum control widgets
        self.sptitle = widgets.Text(description='Title',
                            value=self.title, layout=Layout(width='40%'))
        self.spcolor = widgets.Dropdown(description='Color',
                            options=["steelblue"]+list(Colors),value='steelblue',
                            layout=Layout(width='24%'))
        self.splw = widgets.FloatText(description='Linewidth',
                            value=1.0, step=0.1,
                            layout=Layout(width='24%'))
        self.showlogo = widgets.Checkbox(description="Logo", 
                            value = True, layout=widgets.Layout(left='200px'))
        self.axlogo = self.fig.add_axes([.92, .84, .08, .16], visible=True)
        self.axlogo.imshow(plt.imread(Logofile(),'rb'))
        self.axlogo.set_axis_off()
        def switchlogo(e):
            if self.showlogo.value:
                self.axlogo.set_visible(True)
            else:
                self.axlogo.set_visible(False)
        self.showlogo.observe(switchlogo)
        SP = VBox([ HBox([self.sptitle,self.showlogo]),
                    HBox([self.spcolor, self.splw])])
        # integral control widgets
        self.integ = widgets.Checkbox(description='Show',
                            value=False)
        self.scaleint = widgets.FloatSlider(description='Scale',
                            value=0.5, min=0.1, max=10, step=0.05,
                            layout=Layout(width='40%'),
                            continuous_update=REACTIVE)
        self.offset = widgets.FloatSlider(description='Offset',
                            value=0.4, min=0.0, max=1.0, step=0.01,
                            layout=Layout(width='45%'),
                            continuous_update=REACTIVE) 
        self.intlw = widgets.FloatText(description='Linewidth',
                            value=1.5, step=0.1,
                            layout=Layout(width='24%'))
        self.intcolor = widgets.Dropdown(description='Color',
                            options=Colors,value='red',
                            layout=Layout(width='45%'))
        self.intlabelcolor = widgets.Dropdown(description='Label',
                            options=Colors,value='crimson',
                            layout=Layout(width='40%'))
        self.labelx = widgets.FloatText(description='Label X pos',
                            value=1, step=0.1,
                            layout=Layout(width='24%'))
        self.labely = widgets.FloatText(description='Y pos',
                            value=-1, step=0.1,
                            layout=Layout(width='24%'))
        self.labelfont = widgets.BoundedFloatText(description='Size',
                            value=9, mini=1, maxi=128, step=1,
                            layout=Layout(width='24%'))
        self.labelrot = widgets.BoundedFloatText(description='Rot',
                            value=0, mini=-90, maxi=90, step=1,
                            layout=Layout(width='24%'))
        INT = VBox([HBox([widgets.HTML('<b>Integrals</b>'), self.integ]),
                    HBox([self.scaleint,self.offset, self.intlw]),
                    HBox([self.intcolor,self.intlabelcolor]),
                    HBox([self.labelx,self.labely, self.labelfont, self.labelrot])
                    ],layout=Layout(width='60%'))
        # peaks control widgets
        self.peaks = widgets.Checkbox(description='Show',
                value=False)
        markers = ['off','x','X','d','D','v','^', '<','>','|']
        self.marker = widgets.Dropdown(description='Marker',
                            options=markers,value="x",
                            layout=Layout(width='20%'))
        self.pkcolor = widgets.Dropdown(description='Color',
                            options=Colors,value='darkblue',
                            layout=Layout(width='40%'))
        self.pkvalues = widgets.Checkbox(description='Values',
                value=False, layout=Layout(width='30%'))
        self.pkrotation = widgets.BoundedFloatText(description='Rot',
                            value=45, mini=-90, maxi=90, step=1,
                            layout=Layout(width='20%'))
        self.pkfont = widgets.BoundedFloatText(description='Size',
                            value=5, mini=1, maxi=128, step=1,
                            layout=Layout(width='20%'))
        PK = VBox([ HBox([widgets.HTML('<b>Peaks</b>'), self.peaks]),
                    HBox([self.marker, self.pkcolor]),
                    HBox([self.pkvalues, self.pkfont, self.pkrotation])
                  ],layout=Layout(width='60%'))

        self.Chooser = FileChooser_code(path='/DATA/',filename='.gs1') #(base=base, filetype="*.gs1", mode='r', show=False)
        self.bsel = widgets.Button(description='Copy',layout=self.blay,
                 button_style='info', # 'success', 'info', 'warning', 'danger' or ''
                 tooltip='copy selected data-set to entry below')
        self.to = widgets.IntText(value=1,min=1,max=N,layout=Layout(width='10%'))
        self.bsel.on_click(self.copy)
        self.DataList = [SpforSuper(i+1, axref=self.ax, ref=self.data) for i in range(N)]
        for s in self.DataList:
            s.injectobserve(self.ob)
        self.DataList[0].color.value = 'red'

        for widg in (self.sptitle, self.spcolor, self.splw,
                    self.peaks, self.integ, self.scaleint, self.offset, self.intlw,
                    self.intcolor, self.intlabelcolor, self.labelx, self.labely, self.labelfont, self.labelrot,
                    self.pkvalues, self.marker, self.pkcolor, self.pkrotation, self.pkfont ):
            widg.observe(self.ob)

        orig = self.children
        controls = [SP]
        try:
            self.data.peaks
            controls.append(PK)
        except:
            pass
        try:
            self.data.integrals
            controls.append(INT)
        except:
            pass
        
        self.tabs = Tab()

        self.tabs.children = [
            VBox([  VBox(controls),
                    HBox(orig),
                    ]),
            VBox([  Label("Choose spectra to superimpose"),
                    HBox([self.Chooser, self.bsel, self.to]),
                    self.DataList[0].header] +
                    [sp.me for sp in self.DataList]
                )
            ]

        self.tabs.set_title(0, 'Spectrum')
        self.tabs.set_title(1, 'Superimpose')

        self.children = [self.tabs]
        if show:
            self.draw()


    def copy(self, event):
        if self.to.value <1 or self.to.value >len(self.DataList):
            jsalert('Destination is out of range !')
        else:
            entry = self.DataList[self.to.value-1] 
            entry.filename.value = self.Chooser.selected
            entry.direct.value = 'up'
            for w in entry.me.children:
                w.observe(self.ob)
            self.to.value = min(self.to.value, len(self.DataList)) +1
    def on_done(self, e):
        self.close()
        self.disp(zoom=True)
        display(self.fig)
    def draw(self, zoom=False):
        "refresh display - if zoom is True, display only in xbound"
        super().draw()
        self.xb = self.ax.get_xbound()
        if zoom:
            zoom = self.xb
        else:
            zoom = None
        self.data.display(new_fig=False, figure=self.ax, color=self.spcolor.value,
                            title=self.sptitle.value, linewidth=self.splw.value, zoom=zoom)
        if self.integ.value:
            try:
                if self.labely.value == -1:
                    labely = None
                else:
                    labely = self.labely.value
                self.data.display_integral(label=True, integscale=self.scaleint.value,
                    integoff=self.offset.value,
                    labelxposition=self.labelx.value,
                    labelyposition=labely,
                    curvedict = {'color':self.intcolor.value, 'linewidth':self.intlw.value},
                    labeldict = {'color':self.intlabelcolor.value, 'fontsize':self.labelfont.value, 'rotation':self.labelrot.value},
                    figure=self.ax, zoom=zoom)
            except:
                print('no or wrong integrals (have you clicked on "Done" in the Integration tool ?)')
                pass
        if self.peaks.value:
            try:
                self.data.display_peaks(peak_label=self.pkvalues.value, 
                                        color=self.pkcolor.value,
                                        markersize=self.pkfont.value,
                                        markerdict={'marker':self.marker.value},
                                        labeldict={'rotation':self.pkrotation.value,
                                                    'fontsize':self.pkfont.value},
                                        figure=self.ax, scale=self.scale.value, zoom=zoom)
            except:
                print('no or wrong peaklist (have you clicked on "Done" in the Peak-Picker tool ?)')
                pass
        # superimposition
        for s in self.DataList:
            s.draw()

        self.ax.set_xbound(self.xb)
    def disp(self, zoom=False):
        self.xb = self.ax.get_xbound()
        if zoom:
            zoom = self.xb
        else:
            zoom = None
        super().disp()
        for s in self.DataList:
            s.disp()

class Phaser1D(Show1D):
    """
    An interactive phaser in 1D NMR

        Phaser1D(spectrum, ...)

    requires %matplotlib widget

    """
    def __init__(self, data, figsize=None, title=None, maxfirstorder = 360, show=True, create_children=True):
        data.check1D()
        if data.itype == 0:
            jsalert('Data is Real - an Error will be generated \\n\\n Please redo Fourier Transform')
            data.phase(0,0)
        # we'll work on a copy of the data
        super().__init__( data, figsize=figsize, title=title, show=False, create_children=create_children)
        self.ydata = data.get_buffer()   # store (complex) buffer
    
        self.done.description = 'Apply'
        self.p0 = widgets.FloatSlider(description='P0:',min=-200, max=200, step=0.1,
                            layout=Layout(width='100%'), continuous_update=REACTIVE)
        self.p1 = widgets.FloatSlider(description='P1:',min=-maxfirstorder, max=maxfirstorder, step=1.0,
                            layout=Layout(width='100%'), continuous_update=REACTIVE)
        # self.pivot = widgets.FloatSlider(description='pivot:',
        #                 min=0.0, max=self.data.size1,
        #                 step=1, layout=Layout(width='80%'),
        #                 value=0.5*self.data.size1, readout=False, continuous_update=REACTIVE)
        pivl=self.data.axis1.itoc(self.data.size1)
        pivr=self.data.axis1.itoc(0)
        self.pivot = widgets.BoundedFloatText(description='Pivot',
                        value=round(self.data.axis1.itoc(0.5*self.data.size1),2), 
                        min=min(pivr, pivl),
                        max=max(pivr, pivl),
                        format='%.3f',
                        step=0.1, layout=Layout(width='20%'))
        self.cancel = widgets.Button(description="Cancel", button_style='warning')
        self.cancel.on_click(self.on_cancel)
        # remove done button and create an Apply one
        self.done.on_click(self.on_Apply)
        # draw HBox
        if create_children:
            orig = self.children
            self.children = [VBox([
                                HBox([self.cancel, self.pivot, widgets.HTML('<i>set with right-click on spectrum</i>')]),
                                self.p0,
                                self.p1,
                                HBox(orig)])]
        # add interaction
        for w in [self.p0, self.p1]:
            w.observe(self.ob)
        self.pivot.observe(self.on_movepivot)
        # add right-click event on spectral window
        def on_press(event):
            if event.button == 3:
                self.pivot.value = round(event.xdata,4)
        cids = self.fig.canvas.mpl_connect('button_press_event', on_press)
        self.lp0, self.lp1 = self.ppivot()
        if show:
            self.draw()

    def on_cancel(self, b):
        self.close()
        print("no applied phase correction")
    def on_Apply(self, b):
        self.close()
        lp0, lp1 = self.ppivot() # get centered values
        self.data.phase(lp0, lp1)
        print("Applied: data.phase(%.1f,  %.1f)"%(lp0, lp1))
    def ppivot(self):
        "converts from pivot values to centered ones"
        pp = 1.0-(self.data.axis1.ctoi(self.pivot.value)/self.data.size1)
        return (self.p0.value + (pp-0.5)*self.p1.value, self.p1.value)
    def ctopivot(self, p0, p1):
        "convert from centered to pivot values"
        pp = 1.0-(self.data.axis1.ctoi(self.pivot.value)/self.data.size1)
        p0p = p0 - (pp-0.5)*p1
        while p0p > 180:
            p0p -= 360
        while p0p < -180:
            p0p += 360
        return p0p, p1
    def on_movepivot(self, event):
        if event['name']=='value':
            if self.p1.value != 0:   # in that case values do not change 
                self.p0.value, self.p1.value = self.ctopivot(self.lp0, self.lp1)
            else:
                self.disp()
    @debounce(1.0)
    def rescale_p1(self):
        "shifts by +/- 180°  if P1 hits the border"
        if self.p1.value == self.p1.max:
            self.p1.max += 180
            self.p1.min += 180
        if self.p1.value == self.p1.min:
            self.p1.max -= 180
            self.p1.min -= 180
    def ob(self, event):
        "observe changes and start phasing"
        if event['name']=='value':
            self.rescale_p1()
            self.phase_n_disp()
    def draw(self):
        "copied from super() as it does not display the spectrum !"
        self.ax.clear()
        self.drspectrum = self.ax.plot(self.data.axis1.unit_axis()[::2] , self.ydata.real, lw=1 )[0]
        try:
            self.ax.set_xbound(self.xb)
        except AttributeError:               # the very first time
            xxb = self.ax.get_xbound()
            self.ax.set_xlim(xxb[1],xxb[0])
            y1,y2 = self.ax.get_ybound()
            my = max(abs(y1), abs(y2))
            self.yb0 = np.array([my,-my/2])
        self.ax.set_ybound(self.yb0/self.scale.value)
        self.fig.canvas.header_visible = False
        for s in ["left", "top", "right"]:
            self.ax.spines[s].set_visible(False)
        self.ax.yaxis.set_visible(False)
        self.set_on_redraw()
        self.phase()
        ppos = self.pivot.value
        self.drpivot = self.ax.plot([ppos,ppos], self.ax.get_ybound())[0]

    def phase(self):
        self.lp0, self.lp1 = self.ppivot()         # get centered values
        size = len(self.ydata)
        # compute correction in e
        if self.lp0==0:  # we can optimize a little
            e = np.exp(1J*m.radians(self.lp0)) * np.ones( size, dtype=complex)   # e = exp(j ph0)
        else:
            le = m.radians(self.lp0) + (m.radians(self.lp1))*np.linspace(-0.5, 0.5, size)
            e = np.cos(le) + 1J*np.sin(le)
        # then apply
        y = (e*self.ydata).real
        self.drspectrum.set_ydata( y )   # multiply and keep only real values
    def phase_n_disp(self):
        "apply phase and disp"
        self.phase()
        self.disp()
    def disp(self):
        "update display and display pivot"
        super().disp()
        ppos = self.pivot.value
        self.drpivot.set_xdata([ppos,ppos])
        self.drpivot.set_ydata(self.ax.get_ybound())



class AvProc1D:
    "Detailed 1D NMR Processing"
    def __init__(self, filename=""):
        print('WARNING this tool is not functional/tested yet')
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
            value=False, description='Baseline Correction', tooltip='Perform Baseline Correction')
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

from spike.plugins import Peaks
class NMRPeaker1D(Show1D):
    """
    a peak-picker for NMR experiments
    """
    # self.peaks : the defined peaklist, copyied in and out of data
    # self.temppk : the last computed pklist
    def __init__(self, data, figsize=None, show=True):
        super().__init__( data, figsize=figsize, show=show)
        self.data = data.real()
        try:
            self.peaks = self.data.peaks
        except AttributeError:
            self.peaks = Peaks.Peak1DList(source=self.data)
        self.temppk = Peaks.Peak1DList(source=self.data)
        self.thresh = widgets.FloatLogSlider(value=20.0,
            min=-1, max=2.0, base=10, step=0.01, layout=Layout(width='30%'),
            continuous_update=False, readout=True, readout_format='.2f')
        try: 
            self.thresh.value = 100*self.data.peaks.threshold/self.data.absmax  # if already peak picked
        except:
            self.thresh.value = 20.0
        self.thresh.observe(self.pickpeak)
        self.peak_mode = widgets.Dropdown(options=['marker', 'bar'],value='marker',description='show as')
        self.peak_mode.observe(self.ob)
        self.out = Output(layout={'border': '1px solid red'})
#        self.done = widgets.Button(description="Done", button_style='success')
#        self.done.on_click(self.on_done)
        self.badd = widgets.Button(description="Add", button_style='success', layout=self.blay)
        self.badd.on_click(self.on_add)
        self.brem = widgets.Button(description="Rem", button_style='warning', layout=self.blay)
        self.brem.on_click(self.on_rem)
        self.cancel = widgets.Button(description="Cancel", button_style='warning', layout=self.blay)
        self.cancel.on_click(self.on_cancel)
        self.selval = widgets.FloatText(
            value=0.0, description='selection', layout=Layout(width='20%'), step=0.001, disabled = True)
        self.newval = widgets.FloatText(
            value=0.0, description='calibration', layout=Layout(width='20%'), step=0.001, disabled = True)
        self.setcalib = widgets.Button(description="Set", layout=Layout(width='10%'),
                button_style='success', tooltip='Set spectrum calibration')
        self.setcalib.on_click(self.on_setcalib)
        def on_press(event):
            v = event.xdata
            iv = self.data.axis1.ptoi(v)  # store position in index (peaks are internally in index)
            distclose = np.inf     # search closest peak

            pclose = 0.0
            for p in self.data.peaks:
                if abs(p.pos-iv) < distclose:
                    pclose = p.pos
                    distclose = abs(p.pos-iv)
            self.selval.value = self.data.axis1.itop(pclose)  # back to ppm
            for w in (self.selval, self.newval):
                w.disabled = False
        cids = self.fig.canvas.mpl_connect('button_press_event', on_press)
        # redefine Box
        orig = self.children
        self.tabs = Tab()
        self.tabs.children = [
            VBox([
                HBox([self.badd, self.brem, Label('threshold - % largest signal'), self.thresh, self.peak_mode]),
                HBox([VBox([self.blank, self.reset, self.scale]), self.fig.canvas])
                ]),
            VBox([
                HBox([ Label('Select a peak with mouse and set calibrated values'), self.selval, self.newval, self.setcalib]),
                HBox([VBox([self.blank, self.reset, self.scale]), self.fig.canvas])
                ]),
            self.out]
        self.tabs.set_title(0, 'Peak Picker')
        self.tabs.set_title(1, 'calibration')
        self.tabs.set_title(2, 'Peak Table')

        self.children = [VBox([HBox([self.done, self.cancel]),self.tabs])]

        self.pickpeak({'name':'value'})
        self.disp()
    def on_add(self, b):
        self.peaks.pkadd(self.temppk)
        self.peaks = Peaks.peak_aggreg(self.peaks, distance=1.0)
        self.peaks.source = self.data
        self.temppk = Peaks.Peak1DList(source=self.data)
        self.disp()
    def on_rem(self, b):
        (up,down) = self.ax.get_xbound()
        iup = self.data.axis1.ptoi(up)
        idown = self.data.axis1.ptoi(down)
        iup,idown = (max(iup,idown), min(iup,idown))
        to_rem = []
        for pk in self.peaks:
            if pk.pos < iup and pk.pos>idown:
                to_rem.append(pk)
        for pk in to_rem:
            self.peaks.remove(pk)
        self.disp()
    def on_cancel(self, b):
        self.close()
        del self.data.peaks
        print("no Peak-Picking done")
    def on_done(self, b):
        self.temppk = Peaks.Peak1DList()  # clear temp peaks        
        self.close()
        # new figure
        self.disp()
        # and display        
        display(self.fig)
        display(self.out)
        self.data.peaks = self.peaks  # and copy

    def on_setcalib(self, e):
        off = self.selval.value-self.newval.value
        self.data.axis1.offset -= off*self.data.axis1.frequency   # off is in ppm, axis1.offset is in Hz
        self.selval.value = self.newval.value
#        self.pp()
        self.peaks.pos2label()
        self.temppk.pos2label()
        self.disp()
    def pkprint(self,event):
        self.out.clear_output(wait=True)
        with self.out:
            if len(self.temppk)>0:
                display(HTML("<p style=color:red> Transient peak list </p>"))
                self.data.peaks = self.temppk
                display(HTML(self.data.pk2pandas().to_html()))
            if len(self.peaks)>0:
                display(HTML("<p style=color:blue> Defined peak list </p>"))
                self.data.peaks = self.peaks
                display(HTML(self.data.pk2pandas().to_html()))
            else:
                display(HTML("<p style=color:blue> Defined peak list Empty</p>"))
    def _pkprint(self,event):
        self.out.clear_output(wait=True)
        with self.out:
            print(self.pklist())
    def pklist(self):
        "creates peaklist for printing or exporting"
        text = ["ppm\tInt.(%)\twidth (Hz)"]
        data = self.data
        intmax = max(data.peaks.intens)/100
        for pk in data.peaks:
            ppm = data.axis1.itop(pk.pos)
            width = 2*pk.width*data.axis1.specwidth/data.size1
            l = "%.3f\t%.1f\t%.2f"%(ppm, pk.intens/intmax, width)
            text.append(l)
        return "\n".join(text)
    def ob(self, event):
        if event['name']=='value':
            self.disp()        
    def disp(self):
        "interactive wrapper to peakpick"
        self.xb = self.ax.get_xbound()
        #self.yb = self.ax.get_ybound()
        self.ax.clear()
        #super().disp()
        self.data.display(scale=self.scale.value, new_fig=False, figure=self.ax, title=self.title)
        x = [self.data.axis1.itoc(z) for z in (0, self.data.size1) ]
        y = [self.data.absmax*self.thresh.value/100]*2
        self.ax.plot(x,y,':r')
        try:
            self.temppk.display(peak_label=False, peak_mode=self.peak_mode.value, f=self.data.axis1.itoc, figure=self.ax,color='red')
            self.peaks.display(peak_label=False, peak_mode=self.peak_mode.value, f=self.data.axis1.itoc, color='blue', figure=self.ax)
        except:
            rrr("problem")
        self.temppk.display(peak_label=True, peak_mode=self.peak_mode.value, color='red', figure=self.ax)
        self.ax.set_xbound(self.xb)
        self.ax.set_ylim(ymax=self.data.absmax/self.scale.value)
        self.pkprint({'name':'value'})  # send pseudo event to display peak table
    def pickpeak(self, event):
        "interactive wrapper to peakpick"
        if event['name']=='value':
            self.pp()
    def pp(self):
        "do the peak-picking calling pp().centroid()"
        #self.spec.clear_output(wait=True)
        th = self.data.absmax*self.thresh.value/100
        zm = self.ax.get_xbound()
        self.data.set_unit('ppm').peakpick(threshold=th, verbose=False, zoom=zm).centroid()
        self.temppk = self.data.peaks
        self.data.peaks = None
        self.disp()
        self.ax.annotate('%d peaks detected'%len(self.data.peaks) ,(0.05,0.95), xycoords='figure fraction')
        self.pkprint({'name':'value'})

from spike.plugins.NMR.Integrate import Integrals, Integralitem
class NMRIntegrate(Show1D):
    "an integrator for NMR experiments"
    def __init__(self, data, figsize=None, show=True):
        super().__init__( data, figsize=figsize, show=show)
        try:
            self.Integ = data.integrals
        except:
            self.Integ = Integrals(data, compute=False)   # initialize with empty list
        try:
            self.peaks_reserved = self.data.peaks
        except:
            print('no peaks')
            self.peaks_reserved = None
        self.thresh = widgets.FloatLogSlider(description='sensitivity',value=10.0,
            min=-1, max=2.0, base=10, step=0.01, layout=Layout(width='30%'),
            continuous_update=HEAVY, readout=True, readout_format='.1f',
            tooltip='sensitivity to weak signals')
        self.bias = widgets.FloatSlider(
            description='bias', layout=Layout(width='20%'),
            value=0.0,min=-10.0, max=10.0, step=0.1,
            continuous_update=HEAVY, readout=True, readout_format='.1f')
        self.sep = widgets.FloatSlider(
            description='separation', layout=Layout(width='30%'),
            value=3.0,min=0.0, max=20.0, step=0.1,
            continuous_update=HEAVY, readout=True, readout_format='.1f')
        self.wings = widgets.FloatSlider(
            description='extension', layout=Layout(width='30%'),
            value=5.0,min=0.5, max=20.0, step=0.1,
            continuous_update=HEAVY, readout=True, readout_format='.1f')
        for w in (self.bias, self.sep, self.wings):
            w.observe(self.integrate)
        self.thresh.observe(self.peak_and_integrate)
        self.cancel = widgets.Button(description="Cancel", button_style='warning', layout=self.blay)
        self.cancel.on_click(self.on_cancel)

        self.badd = widgets.Button(description="Add", layout=self.blay,
                button_style='success', tooltip='Add an entry from current zoom')
        self.badd.on_click(self.on_add)
        self.brem = widgets.Button(description="Rem", layout=self.blay,
                button_style='warning', tooltip='Remove all entries in current zoom')
        self.brem.on_click(self.on_rem)
        self.bauto = widgets.Button(description="Compute", layout=self.blay,
                button_style='success', tooltip='Automatic definition of integrals')
        self.bauto.on_click(self.peak_and_integrate)

        self.entry = widgets.IntText(value=0,description='Entry',min=0,layout=Layout(width='25%'),
                                    tooltip='Index of the calibrating integral')
        self.value = widgets.FloatText(value=100,description='Value',layout=Layout(width='25%'),
                                    tooltip='Value of the calibrating integral')
        self.set = widgets.Button(description="Set", button_style='success', layout=Layout(width='10%'))
        self.set.on_click(self.set_value)
        self.out = Output(layout={'border': '1px solid red'})
        # redefine children
        orig = self.children
        self.tabs = Tab()
        self.tabs.children = [
            VBox([
                HBox([self.badd, self.brem,
                    Label("Use buttons to add and remove integrals in the current zoom window")]), 
                HBox(orig),
                ]),
            VBox([HBox([self.bauto, Label("Define integral shapes using the sliders below")]),
                HBox([self.thresh,  self.sep, self.wings]),
                HBox(orig),
                ]),
            VBox([HBox([Label('Choose an integral for calibration'),
                                self.entry, self.value, self.set] ),
                self.out ])
            ]
        self.tabs.set_title(0, 'Manual integration')
        self.tabs.set_title(1, 'Automatic')
        self.tabs.set_title(2, 'Integral Table & Calibration')
        self.children = [VBox([HBox([self.done, self.cancel]),self.tabs])]
        # self.children = [VBox([HBox([self.done, self.cancel, self.bias, self.sep, self.wings]),
        #                     HBox([Label('Choose an integral for calibration'),
        #                         self.entry, self.value, self.set, self.blank, self.bprint] ),
        #                     HBox([VBox([self.blank, self.reset, self.scale]), self.fig.canvas]),
        #                     self.out ])]
    # def show(self):
    #     self.data.display(figure=self.ax)
    #     self.xb = self.ax.get_xbound()  # initialize zoom
    #     #self.int()
    #     self.data.display(figure=self.ax)
    #     self.fig.canvas.header_visible = False
    #     self.ax.set_xbound( (self.data.axis1.itop(0),self.data.axis1.itop(self.data.size1)) )
    #     self.disp()
    #     display(self)
        self.disp()
    def on_cancel(self, b):
        self.close()
        print("No integration")
    def on_done(self, b):
        self.close()
        self.disp(zoom=True)
        display(self.fig)
        display(self.out)
        self.data.integrals = self.Integ        # copy integrals
        if self.peaks_reserved:
            self.data.peaks = self.peaks_reserved
        else:
            del(self.data.peaks)
    def on_add(self, b):
        start, end = self.ax.get_xbound() 
        self.Integ.append( Integralitem(self.data.axis1.ptoi(start), self.data.axis1.ptoi(end), [], 0.0) )
        self.Integ.zonestocurves()
        self.disp()
        self.print(None)
    def on_rem(self, b):
        start, end = self.ax.get_xbound()
        start, end = self.data.axis1.ptoi(start), self.data.axis1.ptoi(end)
        start, end = (min(start, end)-1, max(start, end)+1)  #-1,+1 needed sometimes...
        to_rem = []
        for ii in self.Integ:
            if ii.start>start and ii.end<end:
                to_rem.append(ii)
        for ii in to_rem:
            self.Integ.remove(ii)
        self.disp()
        self.print(None)
    def set_value(self,b):
        self.Integ.recalibrate(self.entry.value, self.value.value)
        self.disp()
        self.print(None)
    def print(self,event):
        self.out.clear_output()
        self.Integ.sort(key=lambda x: x.start)   # sort integrals
        with self.out:
            display(HTML( self.Integ.to_pandas().to_html() ))
    def peak_and_integrate(self, event):
        self.data.pp(threshold=self.data.absmax*self.thresh.value/100,
            verbose=False, zoom=self.ax.get_xbound()).centroid()
        if len(self.data.peaks)>0:
            self.int()
    def integrate(self, event):
        #"integrate from event"
        #if event['name']=='value':
        self.int()
    def int(self):
        "do the automatic integration from peaks and parameters"
        self.on_rem(None) 
        try:
            calib = self.Integ.calibration
        except:
            calib = None
        self.data.set_unit('ppm')
        try:
            self.data.peaks
        except:
            self.data.pp(threshold=self.data.absmax*self.thresh.value/100,
                verbose=False, zoom=self.ax.get_xbound()).centroid()
        # appending to the list, self and Integrals are equivalent to list -
        self.Integ +=  Integrals(self.data, separation=self.sep.value, wings=self.wings.value,
                bias=self.data.absmax*self.bias.value/100)
        self.Integ.calibrate(calibration=calib)
        self.print(None)
        self.disp()
    def ob(self, event):
        if event['name']=='value':
            self.disp()
    def disp(self, zoom=False):
        "refresh display from event - if zoom is True, display only in xbound"
        self.xb = self.ax.get_xbound()
        # self.yb = self.ax.get_ybound()
        self.ax.clear()
        if zoom:
            zoom = self.xb
        else:
            zoom = None
        self.data.display(new_fig=False, figure=self.ax, scale=self.scale.value, zoom=zoom)
        try:
            self.Integ.display(label=True, figure=self.ax, labelyposition=None, regions=False, zoom=zoom)
        except:
            pass
        self.ax.set_xbound(self.xb)
        # self.ax.set_ybound(self.yb)


# if __name__ == '__main__':
#    unittest.main()    
