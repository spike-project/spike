#!/usr/bin/env python
# encoding: utf-8

"""
A set of utilities to use spike in NMR or FTMS within jupyter


First version MAD june 2017
Intermediate version MAD october 2019
Improvements MAD spring - summer - automn  2021 - winter 2022
"""

"""
The code contains several parts

first a bunch of general small gadgets used to improve a jupyter notebook (hidecode, jsalert, injectcss, etc...)

Then some NMR specific tools

The most important is the interactive tools (hence the name):

Show1D
    from which inherits
Show1Dplus
baseline1D
Phaser
NMRintegrate
NMRpeakpicker

all tools are not at the same level on completion, this is an evoluating code...

Naming Conventions are used
For interactive tools (phaser, show, baseline)
draw() creates static matplotlib "artists"
    should probably called once at creation 
    eventually is things change a lot
disp() 
    takes care of redisplaying after some modifications (x or y limits, eventually xdata / ydata ...)
    should note create new artists

"""


import asyncio
from spike.plugins import Peaks
from . import __path__
from spike.plugins.NMR.Integrate import Integrals, Integralitem
from .ipyfilechooser import FileChooser as FileChooser_code
from ..File.BrukerNMR import Import_1D
from ..NMR import NMRData
from .. import NMR
from .. import FTMS
from .. import version
from numpy.core.overrides import array_function_dispatch
import numpy as np
from IPython.display import display, HTML, Javascript, Markdown, Image
import ipywidgets as widgets
from ipywidgets import Layout, HBox, VBox, Label, Output, Button, Tab
from matplotlib.widgets import MultiCursor
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib as mpl
from datetime import datetime
import re
import warnings
import time
import math as m
import glob
import unittest
import sys
import os
import os.path as op
from pathlib import Path
if version.revision < "538":
    warnmsg = """
There is version missmatch between the core program and the interactive tools
You may experiment some difficulties or halting with this notebook
"""
    warnings.warn(warnmsg)


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
Colors = ('black', 'red', 'blue', 'green', 'orange',
          'blueviolet', 'crimson', 'turquoise', 'indigo',
          'magenta', 'gold', 'pink', 'purple', 'salmon', 'darkblue', 'sienna')
Cursor = dict(linestyle='--', color='darkblue',
              linewidth=0.8)   # used for running cursors
Contrast = dict(color="crimson")  # used for contrasting spectra
##############################
# First General Utilities
##############################


def initialize(verb=1):
    "initialize notebook interface for Spike"
    if verb > 0:
        print("\nInteractive module version,", __version__)
    if verb > 1:
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
        {   $('div.input').hide();
            $('div.prompt').hide();
            $('#but').val("show python code");
        } else {
            $('div.input').show();
            $('div.prompt').show();
            $('#but').val("hide python code");
        }
    code_show = !code_show
    } 
$(document).ready(code_toggle);
</script>
<form action="javascript:code_toggle()">
<input id="but" type="submit" style="border:1px solid black; background-color:#DDD">
%s
</form>''' % (init, message)))


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
</form>''' % (init, message)))


def UserLogofile():
    "return the address on disk of the User Logo file - located in $HOME/Spike/Logo.png"
    ufile = Path.home()/'.config'/'Spike'/'Logo.png'
    if ufile.exists():
        return ufile
    else:
        return Logofile()


def Logofile():
    "return the address on disk of the Logo file"
    return op.join(__path__[0], 'Logo.png')


def Logo(width=150):
    "display the Logo"
    display(Image(filename=Logofile(), width=width))


def jsalert(msg="Alert text"):
    "send a javascript alert"
    # use \\n if you want to add several lines
    display(Javascript("alert('%s')" % msg))


def _jsalert(msg="Alert text", title='Alert'):
    display(Javascript("""
    require(
        ["base/js/dialog"], 
        function(dialog) {
            dialog.modal({
                title: '%s',
                body: '%s',
                buttons: {
                    'Ok': {}
                }
            });
        }
    );
    """ % (title, msg)))


def space(width="80px"):
    "defines a layout of width, width should be a string in CSS syntax"
    return widgets.Layout(width=width)


def spacer(width="80px"):
    "defines a spacer to be used in widgets.Box - width should be a string in CSS syntax"
    return widgets.HTML("&nbsp;", layout=space(width))


class _FileChooser:
    """a simple file chooser for Jupyter - obsolete - not used any more"""

    def __init__(self, base=None, filetype=['fid', 'ser'], mode='r', show=True):
        if base is None:
            self.curdir = "/"
        else:
            self.curdir = base
        self.filetype = filetype
        self.mode = mode
        self.wfile = widgets.Text(layout=space(
            '70%'), description='File to load')
        self.ldir = widgets.Label(value="Chosen dir:  "+self.curdir)
        self.wdir = widgets.Select(
            options=self.dirlist(),
            description='Choose Dir',
            layout=space('70%'))
        if mode == 'r':
            self.wchooser = widgets.Select(
                options=self.filelist(),
                description='Choose File',
                layout=Layout(width='70%'))
            self.wchooser.observe(self.wob)
            self.wfile.disabled = True   # not mofiable in read mode
        elif mode == "w":
            self.wfile.description = 'File to create'
            self.wfile.disabled = False
            self.wchooser = widgets.Select(
                options=self.filelist(),
                description='Files',
                layout=space('70%'))
        else:
            raise Exception('Only "r" and  "w" modes supported')
        self.wsetdir = Button(description='❯', layout=space('20%'),
                              button_style='success',  # 'success', 'info', 'warning', 'danger' or ''
                              tooltip='descend in directory')
        self.wup = Button(description='❮', layout=space('20%'),
                          button_style='success',  # 'success', 'info', 'warning', 'danger' or ''
                          tooltip='up to previous directory')
        self.wsetdir.on_click(self.setdir)
        self.wup.on_click(self.updir)
        if show:
            self.show()

    def filelist(self):
        fl = []
        if self.mode == 'r':
            if type(self.filetype) is str:
                fl = glob.glob(op.join(self.curdir, self.filetype))
            elif type(self.filetype) in (tuple, list):
                for f in self.filetype:
                    fl += glob.glob(op.join(self.curdir, f))
            else:
                raise Exception(
                    'TypeError, filetype should be either a string or a list')
        else:   # 'w'
            fl = [f for f in glob.glob(
                op.join(self.curdir, '*')) if op.isfile(f)]
            self.wfile.value = op.join(self.curdir, self.filetype)
        if fl == []:
            fl = [" "]
        return fl

    def dirlist(self):
        base = self.curdir
        return [d for d in glob.glob(op.join(base, '*')) if (op.isdir(d) and not d.endswith('__pycache__'))]

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
            [self.ldir,
             HBox([self.wdir, VBox([self.wup, self.wsetdir])]),
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


# the following set of utilities is to display a nice summary of a SPIKE/Bruker dataset
# first parse parameters in dictionary derived from acqus or procs
# for the program to recognize arrays, the arrayname should be given below ( | separated, no space)
arraynames = 'D|P|PL|PCPD|SP'


def readplist(paramtoadd, paramdict):
    "parse lists from acqus files - only that ones defined in arraynames"
    m = (re.match('(%s)([0-9]+)' % (arraynames,), paramtoadd))
    if m:    # arrays are special !
        i = int(m.group(2))
        pp = paramdict['$'+m.group(1)]
        val = pp[i]
    else:
        val = paramdict['$%s' % paramtoadd]
    return val


ParamList = ['PULPROG', 'SFO1', 'NS', 'TE', 'TD', 'RG', 'SW', 'O1', 'D1', 'P1']


def param_line(acqus, paramlist=None, output='string'):
    """
    given a acqus dictionary, and a list of parameters to extract, 
    this functions return either a plain text or a HTML table depending on output value
    either "string" or "HTML" )
    """
    head = []
    lines = []
    if paramlist == None:
        paramlist = ParamList

    for k in paramlist:
        val = str(readplist(k, acqus))
        if val.startswith('<'):
            val = val[1:-1]
        lines.append(val)
        if k.startswith("$"):
            head.append(k[1:])
        else:
            head.append(k)
    if output == "string":
        res = "\t".join(head) + "\n" + "\t".join(lines)
    elif output == "HTML":
        Head = '<table class="dataframe" border="1"><thead><tr>'
        Head += '<th style="text-align: center;">' + \
            '</th><th style="text-align: center;">'.join(head) + '</th>'
        Head += '</tr></thead>'
        Lines = '<tbody><tr><td>' + \
            '</td><td>'.join(lines) + '</tr></tbody></table>'
        res = HTML(Head + Lines)
    return res


def param_table(data, output='string'):
    """ return the complete parameter table held into data
    output is either "string" or "HTML
    """
    lines = []
    acqus = data.params['acqu']
    if output == 'string':
        for k in acqus.keys():
            lines.append("%s : %s" % (k, acqus[k]))
        res = "\n".join(lines)
    if output == 'HTML':
        st = '<table class="dataframe" border="1"><thead><tr><th><b>Parameter</b></th><th style="text-align: left;">value</th></tr>'
        st += '</thead><tbody>\n'
        for k in acqus.keys():
            val = str(acqus[k]).replace("'", "") \
                .replace("&", '&amp;') \
                .replace("[", '') \
                .replace("]", '') \
                .replace("<", '&lsaquo;') \
                .replace(">", '&rsaquo;') \
                .replace('"', '&quot;')
            lines.append(
                '<tr><td>%s</td><td style="text-align: left;">%s</td></tr>' % (k, val))
        st += ('\n'.join(lines) + '\n</tbody></table>\n')
        res = st
    return res


def summary(data, param=True, output='string'):
    """ produces a summary of the data set
    if param is True, a short list of 1D parameters is given, using the globally defined Paramlist
    output is either "string" or "HTML"
    """
    if output == 'string':
        res = data.params['acqu']['title']
        if param:
            res += "\n" + param_line(data.params['acqu'], output='string')
    elif output == 'HTML':
        res = HTML('<b>' + data.params['acqu']['title'] + '</b>')
        if param:
            res2 = param_line(data.params['acqu'], output='HTML')
            res = HTML(res.data + res2.data)
    return res


def popup_param_table(data):
    css = """<style>
    table { border-collapse: collapse; border: 3px solid #eee; }
    table tr th:first-child { background-color: #eeeeee; color: #333; font-weight: bold }
    table thead th { background-color: #eee; color: #000; }
    tr, th, td { border: 1px solid #ccc; border-width: 1px 0 0 1px; border-collapse: collapse;
    padding: 3px; font-family: monospace; font-size: 10px }</style>
    """
    s = '<script type="text/Javascript">'
    s += 'var win = window.open("", "Parameters", "toolbar=no, location=no, directories=no, status=no, menubar=no, scrollbars=yes, resizable=yes, width=780, height=200, top="+(screen.height-400)+", left="+(screen.width-840));'
    s += 'win.document.body.innerHTML = \''
    s += (param_table(data, output='HTML') + css).replace("\n", '\\') + '\';'
    s += '</script>'
    return(HTML(s+css))
# ²


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


"""
There is a hierarchic structure of the tools
All these tools are Jupyter widgets which can be mixed with other regulat ipywidgets
- Show1D creates the basic environement
    it sets 2 functions called back by all actions
    - draw() which build the whole picture
    - disp() which refresh only the bounding box (in x and y), and eventually modify simple objects
    - Reset which redraw and reset to default display parameter
    - on-done() which freezes the current display as a static 
- Phaser1D, NMRPeaker1D, NMRIntegrate, baseline1D  are surcharging Show1d, 
    with new methods, and eventually surcharging draw() and disp()
- Show1Dplus() which presents all the options

"""


class Show1D(HBox):
    """
    An interactive display, 1D NMR
        Show1D(spectrum)
    """

    def __init__(self, data, title=None, figsize=None, show=True, create_children=True, yratio=21.0):
        """
        data : to be displayed
        title : text for window
        figsize : size in inches (x,y)
        yratio : inverse minimum neg display extension, 1=centered, 21:all positive
        """
        super().__init__()
        self.yratio0 = yratio            #  
        self.yratio = self.yratio0
        self.data = data
        self.drspectrum = None     # will store matplotlib artist
        try:
            (_, self.noise) = data.robust_stats()    # noise will be handy !
        except AttributeError:
            self.noise = 0
        self.title = title
        # graphic set-up details
        self.yaxis_visible = False
        self.ax_spines = False
        self.canvas_header_visible = False
        if figsize is None:
            figsize = mpl.rcParams['figure.figsize']  # (10,5)
        if reverse_scroll:
            self.reverse_scroll = -1
        else:
            self.reverse_scroll = 1
        self.blay = space('80px')       # default layout for buttons
        self.blank = widgets.HTML("&nbsp;", layout=self.blay)
        self.done = Button(description="Done", button_style='success', layout=self.blay,
                           tooltip="Stop interactivity and store figure")
        self.done.on_click(self.on_done)
        self.reset = Button(description="Reset", button_style='success', layout=self.blay,
                            tooltip="Reset to default display")
        self.savepdf = widgets.Button(description="Save figure", button_style='', layout=self.blay,
                                      tooltip='Save the current figure as pdf')

        self.savepdf.on_click(self.onsave)
        self.reset.on_click(self.on_reset)
        self.scale = widgets.FloatSlider(description='Scale:', value=1.0, min=0.1, max=200, step=0.1,
                                         layout=Layout(width='80px', height=str(0.8*2.54*figsize[1])+'cm'), continuous_update=REACTIVE,
                                         orientation='vertical')
        for widg in (self.scale,):
            widg.observe(self.obdisp)
        self.on_reset() # set to default
        if create_children:
            plt.ioff()
            fi, ax = plt.subplots(figsize=figsize)
            fi.set_tight_layout(True)
            plt.ion()
            self.ax = ax
            # def update_box(axes):
            #     self.xb = self.ax.get_xbound()
            #     self.ax.set_ybound(self.yb0/self.scale.value)
            # self.ax.add_callback(self.update_box)
            self.fig = fi
#            self.xb = self.ax.get_xbound()
            self.on_reset()
            self.children = [
                VBox([self.reset, self.scale, self.savepdf, self.done]), self.fig.canvas]
            if show:
                self.draw()

    # call backs
    def onsave(self, e=None, directory='.', verbose=True):
        "saves a pdf of the current display"
        name = self.fig.get_axes()[0].get_title()
        datedext = datetime.strftime(
            datetime.now(), " %y-%m-%d_%H:%M:%S.pdf")
        name = name.replace('/', '_') + datedext
        if name.startswith('.'):
            name = 'Figure'+name
        filename = op.join(directory,name)
        self.fig.savefig(filename)
        if verbose:
            print('figure saved as: ', filename)
        return filename
    def on_done(self, b):
        self.close()
        display(self.fig)   # shows spectrum

    def on_reset(self, e=None):
        self.scale.value = 1.0
        self.yratio = self.yratio0            # inverse minimum neg display extension 
        try:
            right,left = self.data.axis1.borders
        except AttributeError:   # if self.data not defined
            pass
            self.xb0 = (0.0, 10.0)
        else:
            a,b = (self.data.axis1.itoc(right), self.data.axis1.itoc(left))  # full box in x
            self.xb0 = (a,b) #(min(a,b), max(a,b))
            # /!\    xb and xb0 are always   (low, high)  so (Right - left)  in ppm
            self.ymax = self.data.absmax*1.1      # highest position
            self.yb0 = np.array( [-self.ymax/self.yratio, self.ymax] )       # full box in y
        self.xb = self.xb0                                                   # current x box
    def ob(self, event):
        "observe events and redraw"
        if event['name'] == 'value':
            self.draw()

    def obdisp(self, event):
        "observe events and display"
        if event['name'] == 'value':
            self.disp()

    @debounce(0.005)
    def scale_up(self, step):
        # 1.1892 is 4th root of 2.0
        self.scale.value *= 1.1892**(self.reverse_scroll*step)

    def set_on_redraw(self):
        "set-up mouse scroll control and click"
        @debounce(0.005)
        def on_draw(event):
            self.xb = self.ax.get_xbound()      # follow x and y box
            yb = self.ax.get_ybound()
            #  [-ymax/(yratio*scale), ymax/scale = [yb0, yb1]
            #  ymax/scale = yb1   =>  scale = ymax/yb1
            # -ymax/(yratio*scale) = yb0    =>   -ymax/(scale*yb0) = yratio
            ratio = -yb[1]/yb[0]
            if abs(ratio-self.yratio)/ratio > 0.1:      # more than 10% variation
                self.scale.value = self.ymax/yb[1]
                self.yratio = -self.ymax/(yb[0]*self.scale.value)
                self.yb0 = np.array( [-self.ymax/self.yratio, self.ymax] ) 
        self.cids_disp = self.fig.canvas.mpl_connect('draw_event', on_draw)
        def on_scroll(event):
            self.scale_up(event.step)
        if Activate_Wheel:
            self.fig.canvas.capture_scroll = True
            self.cids_scroll = self.fig.canvas.mpl_connect('scroll_event', on_scroll)

    def draw(self):
        "builds and display the picture"
        self.ax.clear()
        self.data.display(new_fig=False, figure=self.ax, title=self.title)
        self.drspectrum = self.data.disp_artist
        self.ax.set_xbound(*self.xb)
        self.ax.set_ybound(self.yb0/self.scale.value)   # self.yb0 is [low, up] np.array
        self.fig.canvas.header_visible = False
        for s in ["left", "top", "right"]:
            self.ax.spines[s].set_visible(False)
        self.ax.yaxis.set_visible(False)
        visible_ticks = {"top": False,   "right": False}
        self.ax.tick_params(axis="x", which="both", **visible_ticks)
        self.set_on_redraw()

    def disp(self):
        # scale = self.ax.get_ybound()
        # self.scale.value = self.yb0/scale
        # self.xb = self.ax.get_xbound()
        self.ax.set_ybound(self.yb0/self.scale.value)
        self.ax.set_xbound(*self.xb)


class baseline1D(Show1D):
    """
    An interactive baseline for 1D NMR

        baseline1D(spectrum, ...)

    """

    def __init__(self, data, figsize=None, show=True, create_children=True, yratio=10):
        try:
            self.BslPoints = data.BslPoints    # recover from previous correction
        except:
            self.BslPoints = []   # initialize with empty list

        super().__init__(data, figsize=figsize,
                        show=False, create_children=create_children,yratio=yratio)
        self.data.real()
        a, b = self.itoc3(0.0), self.itoc3(self.data.size1)   # spectral borders
        self.select = widgets.FloatSlider(description='Manually:',
                                          min=min(a, b), max=max(a, b),
                                          value=self.itoc3(0.5*self.data.size1), continuous_update=REACTIVE,
                                          layout=space('40%'),
                                          tooltip='position of the selected point in %s' % (data.unit[0]))
        self.smooth = widgets.IntSlider(description='Smoothing:', min=0, max=20,  layout=space('40%'),
                                        tooltip='apply a local smoothing for pivot points',
                                        value=1, readout=True, continuous_update=REACTIVE)
        self.toshow = widgets.Dropdown(
            options=['full', 'corrected', 'baseline', 'points'],  description='Display:')

        self.bsl_points = self.BslPoints
        self.baseline = np.zeros_like(self.data.buffer)
        self.cancel = widgets.Button(description="Exit", button_style='warning', layout=self.blay,
                                     tooltip='Exit without corrections')
        self.auto = widgets.Button(description="Auto", button_style='success', layout=self.blay,
                                   tooltip='Automatic set of points')
        self.set = widgets.Button(description="Add", button_style='success', layout=self.blay,
                                  tooltip='Add a baseline point at the selector position')
        self.unset = widgets.Button(description="Rem", button_style='warning', layout=self.blay,
                                    tooltip='Remove the baseline point closest to the selector')
        self.nbautopoints = widgets.IntText(value=8, layout=space('4em'))
        if create_children:
            orig = self.children
            self.children = [VBox([
                HBox([Label('Baseline points'),
                      self.nbautopoints,
                      Label("points"),
                      self.auto,
                      self.blank,
                      self.select,
                      self.set, self.unset]),
                HBox([self.cancel, self.toshow, self.smooth]),
                HBox(orig)])]
        # add interaction
        # done button is activated by base class
        self.cancel.on_click(self.on_cancel)
        self.auto.on_click(self.on_auto)
        self.set.on_click(self.on_set)
        self.unset.on_click(self.on_unset)
        for w in [self.select, self.smooth, self.toshow]:
            w.observe(self.ob)

        def on_press(event):
            if event.button == 3:       # right-click
                v = event.xdata
                self.select.value = round(v, 3)
                self.disp()
        self.cids_press = self.fig.canvas.mpl_connect('button_press_event', on_press)
        # finalize
        if show:
            self.draw()

    def disconnect_press(self):
        "should be called before close()"
        self.fig.canvas.mpl_disconnect(self.cids_press)

    def itoc3(self, value):
        return round(self.data.axis1.itoc(value), 3)

    def on_cancel(self, e):
        self.disconnect_press()
        self.close()
        print('no baseline correction')

    def on_done(self, e):
        if self.bsl_points == []:
            jsalert(
                'Please define control points before applying baseline correction')
            return
        printlist = [round(i, 4) for i in sorted(self.bsl_points)]
        if len(self.bsl_points) < 4:
            meth = 'linear'
        else:
            meth = 'spline'
        print(f"""Applied correction:
data.bcorr(method='{meth}', nsmooth={self.smooth.value}, xpunits='current',\\
           xpoints={printlist}))
""")
        sp = self.data     # short name for convenience
        sp.set_buffer(sp.get_buffer() - self.correction())
        self.drspectrum[0].set_xdata( sp.axis1.itoc( sp.axis1.points_axis() ) )    # I got once a shorten x display axis !!!  so restore here
        self.drspectrum[0].set_ydata(sp.get_buffer())
        self.clear()   # remove additional drawings
        self.drpivot.set_visible(True)   # but pivot points
        # self.ax.clear()
        # self.data.display(new_fig=False, figure=self.ax, title=self.title)
        self.drpivot.set_xdata(self.bsl_points)
        y = np.zeros(len(self.bsl_points))
        self.drpivot.set_ydata(y)
        # store for latter use
        self.data.BslPoints = self.bsl_points
        self.disconnect_press()
        self.close()
        display(self.fig)   # shows spectrum

    def on_auto(self, e):
        "automatically set baseline points"
        self.bsl_points = [self.data.axis1.itoc(x) for x in bcorr.autopoints(self.data, Npoints=self.nbautopoints.value)]
        self.disp()

    def on_set(self, e):
        "add baseline points at selector"
        self.bsl_points.append(self.select.value)
        self.disp()

    def on_unset(self, e):
        "remove baseline points closest from selector"
        if len(self.bsl_points) == 0:
            return
        here = self.select.value   # find
        distclose = np.inf
        pclose = np.NaN
        for i, p in enumerate(self.bsl_points):
            if abs(p-here) < distclose:
                pclose = p
                ipclose = i
                distclose = abs(p-here)
        self.bsl_points.remove(pclose)  # remove and clean display
        self.disp()

    def correction(self):
        "returns the correction to apply as a numpy array"
        ibsl_points = self.data.axis1.ptoi(
            np.array(self.bsl_points)).astype(int)
        x = np.arange(self.data.size1)
        yy = self.data.get_buffer()
        if len(self.bsl_points) == 0:
            return 0.0
        elif len(self.bsl_points) == 1:
            value = self.data[ibsl_points[0]] * np.ones(self.data.size1)
        elif len(self.bsl_points) < 4:
            corr = bcorr._linear_interpolate(
                yy, ibsl_points, nsmooth=self.smooth.value)
            value = corr(x)
        else:
            corr = bcorr._spline_interpolate(
                yy, ibsl_points, nsmooth=self.smooth.value)
            value = corr(x)
        return value

    def corrected(self):
        value = self.data.copy()
        value.set_buffer(value.get_buffer() - self.correction())
        return value

    def clear(self):
        "remove artists created by draw()"
#        self.drspectrum[0].set_visible(False)
        self.selector.set_visible(False)
        self.drhoriz.set_visible(False)
        self.drbaseline.set_visible(False)
        self.drpivot.set_visible(False)

    def draw(self):
        "used to create the view with additional artists"
        super().draw()
        # baseline
        self.drbaseline = self.ax.plot(
            self.data.axis1.unit_axis(), self.baseline, **Contrast)[0]
        # corrected spectrum
        corrected = self.data.get_buffer()-self.correction()
        # yoffset is used to make fig clearer
        yoffset = 0.2*self.ax.get_ybound()[1]
        self.drcorrected = self.ax.plot(
            self.data.axis1.unit_axis(), corrected+yoffset, alpha=0.5, **Contrast)
        # horizontal line at yoffset
        self.drhoriz = self.ax.axhline(yoffset, **Cursor)
        # selector
        ppos = self.select.value
        self.selector = self.ax.axvline(ppos, **Cursor)
        # pivots
        if len(self.bsl_points) > 0:
            y = bcorr.get_ypoints(self.data.get_buffer(),
                                  self.data.axis1.ptoi(
                                      np.array(self.bsl_points)),
                                  nsmooth=self.smooth.value)
            self.drpivot = self.ax.plot(
                self.bsl_points, y, marker='o', linestyle="None", **Contrast)[0]  # , markersize=12)
        else:
            self.drpivot = self.ax.plot(
                [], [], marker='o', linestyle="None", **Contrast)[0]  # , markersize=12)
        self.disp()

    def disp(self):
        "used to refresh view"
        # graphic objects: (drspectrum) drbaseline drpivot drcorrected drhoriz selector
        super().disp()
        if len(self.bsl_points) > 0:
            if self.toshow.value in ('baseline', 'full'):
                if self.toshow.value == 'baseline':
                    self.drcorrected[0].set_visible(False)
                    self.drhoriz.set_ydata(0.0)
                self.drhoriz.set_visible(True)
                self.drbaseline.set_ydata(self.correction())
                self.drbaseline.set_visible(True)
            if self.toshow.value in ('corrected', 'full'):
                yoffset = 0.2*self.ax.get_ybound()[1]
                self.drcorrected[0].set_ydata( self.data.get_buffer() - self.correction() + yoffset)
                self.drhoriz.set_ydata(yoffset)
                self.drhoriz.set_visible(True)
                self.drcorrected[0].set_visible(True)
                if self.toshow.value == 'corrected':
                    self.drbaseline.set_visible(False)
            if self.toshow.value == 'points':
                self.drbaseline.set_visible(False)
                self.drcorrected[0].set_visible(False)
                self.drhoriz.set_visible(False)
        else:
            self.drbaseline.set_visible(False)
            self.drcorrected[0].set_visible(False)
            self.drhoriz.set_visible(False)

        # selector
        ppos = self.select.value
        self.selector.set_xdata([ppos, ppos])
        # pivots
        # pivot points
        y = bcorr.get_ypoints(self.data.get_buffer(),
                              self.data.axis1.ptoi(np.array(self.bsl_points)),
                              nsmooth=self.smooth.value)
        #self.drpivot = self.ax.scatter(self.bsl_points, y,  c='r', marker='o')
        # self.drpivot = self.ax.plot(self.bsl_points, y, marker='o',linestyle="None", **Contrast)[0]#, markersize=12)
        self.drpivot.set_xdata(self.bsl_points)
        self.drpivot.set_ydata(y)
        # set zoom


def B(text, size='auto'):
    "bold HTML widget"
    return widgets.HTML(value="<b>%s</b>" % text, layout=widgets.Layout(width=size))

def I(text, size='auto'):
    "Italic HTML widget"
    return widgets.HTML(value="<i>%s</i>" % text, layout=widgets.Layout(width=size))
class SpforSuper():
    "a holder for one spectrum to SuperImpose"
    header = HBox([B('Id', '40px'),
                    B('file name', '5em'),
                    I('or "Self" for spectrum excerpts', '26em'),
                    B('scale', '60px'),
                    B('stretching', '60px'),
                    B('direction', '60px'),
                    B('xoffset', '50px'),
                    B('yoffset(%)', '50px'),
                    B(' ', '50px'), B('Zoom', '90px'),
                    B('Linewidth', '80px'),
                    B('color', '70px'),
                    B('label')
                    ])

    def __init__(self, i, filename='None', axref=False, ref=None):
        """
        filename is the dataset to show  (None or Self possible)
        axref is the axes of the caller, where to display
        ref is a dataset in case filename is Self 
        """
        self.Id = i
        self.axref = axref
        self.ref = ref
        self.drartist = None
        self.filename = widgets.Text(
            value=filename, layout=Layout(width='30em'))
        self.filename.observe(self.nactivate)
        j = i % len(Colors)
        self.color = widgets.Dropdown(
            options=["steelblue"]+list(Colors), value=Colors[j], layout=space('90px'))
        self.direct = widgets.Dropdown(
            options=['up', 'down', 'off'], value='off', layout=space('60px'))
        self.direct.observe(self.activate)
        self.zmleft = widgets.FloatText(
            value=1000, step=0.1, layout=space('70px'))
        self.zmright = widgets.FloatText(
            value=-1000, step=0.1, layout=space('70px'))
        self.label = widgets.Checkbox(
            value=False, indent=False, layout=space('40px'))
        self.scale = widgets.FloatText(value=1.0, step=0.1, layout=space(
            '70px'), tooltip="relative intensity scale")
        self.stretch = widgets.FloatText(
            value=1.0, step=0.1, layout=space('70px'), tooltip="stretching along x")
        self.xoffset = widgets.FloatText(value=0.0, step=0.1, layout=space(
            '60px'), tooltip="x offset in spec unit")
        self.yoffset = widgets.FloatText(
            value=0.0, step=1, layout=space('60px'), tooltip="y offset in %")
        self.splw = widgets.FloatText(
            value=1.0, step=0.1, layout=space('50px'))
        # for w in (self.color, self.direct, self.zmleft, self.zmright, self.label, self.scale, self.offset, self.splw):
        #     w.observe(self.ob)
        # control bar
        self.nactivate(None)
        self.xbox = np.array([None, None])   # will store the borders of the spectrum

    @property
    def me(self):
        return HBox([B(str(self.Id)),
                    self.filename,
                    self.scale,
                    self.stretch,
                    self.direct,
                    self.xoffset,
                    self.yoffset,
                    self.zmleft,
                    self.zmright,
                    self.splw,
                    self.color,
                    self.label,
                    ])
    def ob(self, event):
        "observe events and display"
        if event['name'] == 'value':
            self.draw()

    def injectobserve(self, method):
        "used to to observe form outer object for draw"
        for w in (self.color,  self.zmleft, self.zmright, self.label, self.splw, self.direct,
                  self.scale, self.stretch, self.xoffset, self.yoffset):
            w.observe(method)

    def nactivate(self, e):
        "called to check name validity, when name is entered"
        fileexists = (self.filename.value == 'Self' or
                      (self.filename.value not in ('None', '') and Path(self.filename.value).exists()))
        if fileexists:
            self.direct.value = 'up'  # will call self.activate()
            self.direct.disabled = False
        else:
            self.direct.value = 'off'  # will call self.activate()
            self.direct.disabled = True

    def load(self, name):
        """load a file as a NMRData from name - this one implements only native Gifa format
        if you want more, you need to overload it with your own version
        """
        return NMRData(name=name)

    def activate(self, e):
        "called when state changes"
        for w in (self.color, self.zmleft, self.zmright, self.label, self.scale, self.stretch,
                  self.xoffset, self.yoffset, self.splw):
            if self.direct.value == 'off':
                w.disabled = True
            else:
                w.disabled = False

    @property
    def todraw(self):
        "True if there is something to draw"
        if self.filename.value == 'None' or self.filename.value == '' or self.direct.value == 'off':
            return False
        else:
            return True

    def draw(self):
        if not self.todraw:    return  # nothing to do
        # else
        if self.drartist is not None:   # remove previous ones
            del(self.drartist)
        # label
        if self.label.value:
            m = self.mult
            if not (abs(m) == 1 or m == 0):
                lb = "%s (x%.1f)" % (self.nmrname,m)
            else:
                lb = self.nmrname
        else:
            lb = None
        # load
        unit = self.ref.axis1.currentunit
        if self.filename.value == 'Self':
            d = self.ref.copy().set_unit(unit)
        else:
            try:
                d = self.load(name=self.filename.value).set_unit(unit)
            except:
                self.direct.value = 'off'
                jsalert('%s : cannot open File' % self.filename.value)
                return
        # xbox
        self.xbox = np.array([d.axis1.itoc(0), d.axis1.itoc(d.size1)])  # full box in x

        # draw
        zml = min(self.zmleft.value, self.xbox[0])
        zmr = max(self.zmright.value, self.xbox[1])
        d.display(
            scale=1,
            figure=self.axref,
            zoom=(zml, zmr),
            color=self.color.value,
            linewidth=self.splw.value,
            label=lb)
        self.drartist = d.disp_artist[0]
        self.ydata = self.drartist.get_ydata()   # store spectral data
        self.xdata = self.drartist.get_xdata()
        # adapt to scale and direction
        # self.drartist.set_ydata(self.move_ydata())
        # self.drartist.set_xdata(self.move_xdata())
        # text
        # will be changed by disp()
        self.text = self.axref.text(
            zmr, self.ydata[-1], ' ', color=self.color.value)
        self.disp()

    def move_xdata(self):
        "adapt self.xdata to  offset and stretching, then  return the value"
        # get starting value
        start = self.xdata[0]
        # affine transform
        return self.xdata[0] + self.stretch.value*(self.xdata - start) + self.xoffset.value

    def move_ydata(self):
        "adapt self.ydata to direct, scale, and offset and return the value"
        # add offset
        # relative to screen in spectral unit
        off = self.axref.get_ybound()[1] * self.yoffset.value/100
        # and multiply
        return self.mult*self.ydata + off

    def disp(self):
        " refresh "
        if self.filename.value == 'None' or self.direct.value == 'off':
            return  # nothing to do
        yd = self.move_ydata()
        xd = self.move_xdata()
        self.drartist.set_ydata(yd)
        self.drartist.set_xdata(xd)
        self.text.set_text(" ")     # empty text
        if not self.label.value:    # unless... there is a scale != 1, and no label
            m = self.mult
            if not (abs(m) == 1 or m == 0):
                tt = "x%.1f" % (m,)
                self.text.set_text(tt)
                self.text.set_y(yd[-1])
                self.text.set_x(xd[-1])

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
        if fnm not in ('None', 'Self'):
            return op.join(op.basename(op.dirname(fnm)), op.basename(fnm))
        elif fnm == 'Self':
            return 'excerpt'
        else:
            return None

    def simpledisplay(self):
        "a simplified widget "
        lb = widgets.Label(self.nmrname)
        return HBox([lb, self.scale])


class Show1Dplus(Show1D):
    def __init__(self, data, base='/DATA', N=9, **kw):
        if isinstance(data, FTMS.FTMSData):  # pb is in xbox for Datalist, as MS goes reverse to NMR
            raise NotImplementedError('FTMS not implemented yet')
        show = kw.get('show', True)
        kw['show'] = False
        super().__init__(data, **kw)

        # initialise the DataList (used by spectral superposition)
        if N>0:
            self.DataList = [SpforSuper(i+1, axref=self.ax, ref=self.data) for i in range(N)]
            for s in self.DataList:
                s.injectobserve(self.ob)
            self.DataList[0].color.value = 'red'
        else:
            self.DataList = []
        
        self.NbMaxPeaks = Peaks.NbMaxDisplayPeaks
        # spectrum control widgets
        self.sptitle = widgets.Text(description='Title',
                                    value=self.title, layout=space('40%'))
        self.spcolor = widgets.Dropdown(description='Color',
                                        options=["steelblue"]+list(Colors), value='steelblue',
                                        layout=space('24%'))
        self.splw = widgets.FloatText(description='Linewidth',
                                      value=1.0, step=0.1,
                                      layout=space('24%'))
        self.showlogo = widgets.Checkbox(description="Logo",
                                         value=True, layout=widgets.Layout(left='200px'))
        self.axlogo = self.fig.add_axes([.92, .84, .08, .16], visible=True)
        self.axlogo.imshow(plt.imread(UserLogofile(), 'rb'))
        self.axlogo.set_axis_off()

        def switchlogo(e):
            if self.showlogo.value:
                self.axlogo.set_visible(True)
            else:
                self.axlogo.set_visible(False)
        self.showlogo.observe(switchlogo)
        SP = VBox([HBox([self.sptitle, self.showlogo]),
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
                                         options=Colors, value='red',
                                         layout=Layout(width='45%'))
        self.intlabelcolor = widgets.Dropdown(description='Label',
                                              options=Colors, value='crimson',
                                              layout=Layout(width='40%'))
        self.labelx = widgets.FloatText(description='Label X pos',
                                        value=1, step=0.1,
                                        layout=Layout(width='24%'))
        self.labely = widgets.FloatText(description='Y pos',
                                        value=-1, step=0.1,
                                        layout=Layout(width='24%'))
        self.labelfont = widgets.BoundedFloatText(description='Size',
                                                  value=8, mini=1, maxi=128, step=1,
                                                  layout=Layout(width='24%'))
        self.labelrot = widgets.BoundedFloatText(description='Angle',
                                                 value=0, mini=-90, maxi=90, step=1,
                                                 layout=Layout(width='24%'))
        INT = VBox([HBox([widgets.HTML('<b>Integrals</b>'), self.integ]),
                    HBox([self.scaleint, self.offset, self.intlw]),
                    HBox([self.intcolor, self.intlabelcolor]),
                    HBox([self.labelx, self.labely, self.labelfont, self.labelrot])
                    ], layout=Layout(width='60%'))
        # peaks control widgets
        self.peaks = widgets.Checkbox(description='Show',value=False, tooltip="show/hide peaks on the spectrum")
        markers = ['None', 'x', 'X', 'd', 'D', 'v', '^', '<', '>', '|']
        self.marker = widgets.Dropdown(description='Marker',
                                       options=markers, value="x",
                                       layout=Layout(width='20%'))
        self.pkcolor = widgets.Dropdown(description='Color',
                                        options=Colors, value='darkblue',
                                        layout=Layout(width='40%'))
        self.pkvalues = widgets.Checkbox(description='Values',
                                         value=False, layout=Layout(width='30%'))
        self.pkrotation = widgets.BoundedFloatText(description='Angle',
                                                   value=45, mini=-90, maxi=90, step=1,
                                                   layout=Layout(width='20%'))
        self.pkfont = widgets.BoundedFloatText(description='Size',
                                               value=8, mini=1, maxi=128, step=1,
                                               layout=Layout(width='20%'))
        PK = VBox([HBox([widgets.HTML('<b>Peaks</b>'), self.peaks]),
                   HBox([self.marker, self.pkcolor]),
                   HBox([self.pkvalues, self.pkfont, self.pkrotation])
                   ], layout=Layout(width='60%'))
        # Superposition control widgets
        # (base=base, filetype="*.gs1", mode='r', show=False)
        self.Chooser = FileChooser_code(path='/DATA/', filename='.gs1')
        self.bsel = widgets.Button(description='Copy', layout=self.blay,
                                   button_style='info',  # 'success', 'info', 'warning', 'danger' or ''
                                   tooltip='copy selected data-set to entry below')
        self.to = widgets.IntText(
            value=1, min=1, max=N, layout=Layout(width='10%'))
        self.bsel.on_click(self.copy)

        for widg in (self.sptitle, self.spcolor, self.splw,
                     self.peaks, self.integ, self.scaleint, self.offset, self.intlw,
                     self.intcolor, self.intlabelcolor, self.labelx, self.labely, self.labelfont, self.labelrot,
                     self.pkvalues, self.marker, self.pkcolor, self.pkrotation, self.pkfont):
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

        if len(self.DataList) > 0:
            self.superheader = self.DataList[0].header
        else:
            self.superheader = Label('Superimposition list empty')
        self.tabs.children = [
            VBox([VBox(controls),
                  HBox(orig),
                  ]),
            VBox([Label("Choose spectra to superimpose - first Select, then Copy to the numbered slot"),
                  HBox([self.Chooser, self.bsel, self.to]), self.superheader] + [sp.me for sp in self.DataList]
                 )
        ]

        self.tabs.set_title(0, 'Spectrum')
        self.tabs.set_title(1, 'Superimpose')

        self.children = [self.tabs]
        if show:
            self.draw()

    def on_reset(self, b=None):
        "overload reset() to compute the large xbox"
        super().on_reset(b)
        try:
            self.DataList
        except AttributeError:     # on_reset called  at initialisation  !
            return
        xb0 = self.xb0   # comput by super.reset
        for s in self.DataList:
            xb = s.xbox
            if all(xb != None):     # if defined  - all() because xb is an array
                xb0 = ( min(xb), max(xb) )  # always lower / upper
        self.xb0 = xb0
        self.xb = self.xb0                                                      # current x box
        self.disp()

    def copy(self, event):
        if self.to.value < 1 or self.to.value > len(self.DataList):
            jsalert('Destination is out of range !')
        else:
            entry = self.DataList[self.to.value-1]
            entry.filename.value = self.Chooser.selected
            entry.direct.value = 'up'
            for w in entry.me.children:
                w.observe(self.ob)
            self.to.value = min(self.to.value, len(self.DataList)) + 1

    def on_done(self, e):
        self.close()
        self.disp()
        display(self.fig)

    def draw(self, zoom=False):
        "refresh display - if zoom is True, display only in xbound"
        super().draw()
        self.xb = self.ax.get_xbound()
        if zoom:
            zoom = self.xb
        else:
            zoom = None
        # self.data.display(new_fig=False, figure=self.ax, color=self.spcolor.value,
        #                     title=self.sptitle.value, linewidth=self.splw.value, zoom=zoom)
        self.drspectrum[0].set_color(self.spcolor.value)
        self.drspectrum[0].set_linewidth(self.splw.value)
        self.ax.set_title(self.sptitle.value)
        self.title = self.sptitle.value
        if self.integ.value:
            if hasattr(self.data,'integrals') and self.data.integrals != []:
                if self.labely.value == -1:
                    labely = None
                else:
                    labely = self.labely.value
                self.data.display_integral(label=True, integscale=self.scaleint.value,
                                           integoff=self.offset.value,
                                           labelxposition=self.labelx.value,
                                           labelyposition=labely,
                                           curvedict={
                                               'color': self.intcolor.value, 'linewidth': self.intlw.value},
                                           labeldict={
                                               'color': self.intlabelcolor.value, 'fontsize': self.labelfont.value, 'rotation': self.labelrot.value},
                                           figure=self.ax, zoom=self.ax.get_xbound())
            else:
                # print('no or wrong integrals (have you clicked on "Done" in the Integration tool ?)')
                pass
        if self.peaks.value:
            if hasattr(self.data,'peaks') and self.data.peaks != []:
                self.data.display_peaks(peak_label=self.pkvalues.value,
                                        color=self.pkcolor.value,
                                        markersize=self.pkfont.value,
                                        NbMaxPeaks=self.NbMaxPeaks,
                                        markerdict={
                                            'marker': self.marker.value},
                                        labeldict={'rotation': self.pkrotation.value,
                                                   'fontsize': self.pkfont.value},
                                        figure=self.ax, scale=1.0, zoom=self.ax.get_xbound())
            else:
                #                print('no or wrong peaklist (have you clicked on "Done" in the Peak-Picker tool ?)')
                pass
        # superimposition
        for s in self.DataList:
            s.draw()
        if any([ s.label.value for s in self.DataList]):
            self.ax.legend(loc='upper left')

        self.ax.set_xbound(*self.xb)

    def disp(self, zoom=False):
    #     self.xb = self.ax.get_xbound()
    #     if zoom:
    #         zoom = self.xb
    #     else:
    #         zoom = None
        super().disp()
        for s in self.DataList:
#            print(type(s), s.filename)
            s.disp()


class Phaser1D(Show1D):
    """
    An interactive phaser in 1D NMR

        Phaser1D(spectrum, ...)

    """

    def __init__(self, data, figsize=None, title=None, maxfirstorder=360, show=True, warning=True, create_children=True, yratio=10):
        data.check1D()
        if data.itype == 0:
            if warning:
                jsalert('Data is Real - the Imaginary part is being reconstructed')
            data.real2cpx()
        # we'll work on a copy of the data
        super().__init__(data, figsize=figsize, title=title,
                         show=False, create_children=create_children, yratio=yratio)
        self.ydata = data.get_buffer()   # store (complex) buffer

        # change done button and create an Apply one
        self.done.description = 'Apply'
        self.p0 = widgets.FloatSlider(description='P0:', min=-200, max=200, step=0.1,
                                      layout=Layout(width='100%'), continuous_update=REACTIVE)
        self.p1 = widgets.FloatSlider(description='P1:', min=-maxfirstorder, max=maxfirstorder, step=1.0,
                                      layout=Layout(width='100%'), continuous_update=REACTIVE)
        # self.pivot = widgets.FloatSlider(description='pivot:',
        #                 min=0.0, max=self.data.size1,
        #                 step=1, layout=Layout(width='80%'),
        #                 value=0.5*self.data.size1, readout=False, continuous_update=REACTIVE)
        pivl = self.data.axis1.itoc(self.data.size1)
        pivr = self.data.axis1.itoc(0)
        self.pivot = widgets.BoundedFloatText(description='Pivot',
                                              value=round(self.data.axis1.itoc(
                                                  0.5*self.data.size1), 2),
                                              min=min(pivr, pivl),
                                              max=max(pivr, pivl),
                                              format='%.3f',
                                              step=0.1, layout=Layout(width='20%'))
        self.cancel = widgets.Button(
            description="Exit", button_style='warning', tooltip='Exit without corrections')
        # draw HBox
        if create_children:
            orig = self.children
            self.children = [VBox([
                HBox([self.cancel, self.pivot, widgets.HTML(
                    '<i>set with right-click on spectrum</i>')]),
                self.p0,
                self.p1,
                HBox(orig)])]
        # add interaction
        self.cancel.on_click(self.on_cancel)
        self.done.on_click(self.on_Apply)
        for w in [self.p0, self.p1]:
            w.observe(self.ob)
        self.pivot.observe(self.on_movepivot)
        # add click event on spectral window

        def on_press(event):
            if event.button == 3:       # right-click
                v = event.xdata
                self.pivot.value = round(v, 4)
                self.disp()
        self.cids_press = self.fig.canvas.mpl_connect('button_press_event', on_press)
        # finalize
        self.lp0, self.lp1 = self.ppivot()
        if show:
            self.draw()

    def disconnect_press(self):
        "should be called before close()"
        self.fig.canvas.mpl_disconnect(self.cids_press)

    def on_cancel(self, b):
        self.disconnect_press()
        self.close()
        print("no applied phase correction")

    def on_Apply(self, b):
        self.drpivot.set_visible(False)
        lp0, lp1 = self.ppivot()  # get centered values
        self.data.phase(lp0, lp1)
        self.disconnect_press()
        self.close()
        # display(self.fig)   # would show twice - no idea why ...
        print("Applied: data.phase(%.1f,  %.1f)" % (lp0, lp1))

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
        if event['name'] == 'value':
            if self.p1.value != 0:   # in that case values do not change
                self.p0.value, self.p1.value = self.ctopivot(
                    self.lp0, self.lp1)
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
        if event['name'] == 'value':
            self.rescale_p1()
            self.phase_n_disp()

    def draw(self):
        "copied from super() as it does not display the spectrum !"
#         self.ax.clear()
#         self.drspectrum = self.ax.plot(self.data.axis1.unit_axis()[::2] , self.ydata.real, lw=1 )[0]
        super().draw()
        self.phase()
        ppos = self.pivot.value
        self.drpivot = self.ax.axvline(ppos, **Cursor)

    def phase(self):
        self.lp0, self.lp1 = self.ppivot()         # get centered values
        size = len(self.ydata)
        # compute correction in e
        if self.lp0 == 0:  # we can optimize a little
            e = np.exp(1J*m.radians(self.lp0)) * \
                np.ones(size, dtype=complex)   # e = exp(j ph0)
        else:
            le = m.radians(self.lp0) + (m.radians(self.lp1)) * \
                np.linspace(-0.5, 0.5, size)
            e = np.cos(le) + 1J*np.sin(le)
        # then apply
        y = (e*self.ydata).real
        self.drspectrum[0].set_ydata(y)   # multiply and keep only real values

    def phase_n_disp(self):
        "apply phase and disp"
        self.phase()
        self.disp()

    def disp(self):
        "update display and display pivot"
        super().disp()
        ppos = self.pivot.value
        self.drpivot.set_xdata([ppos, ppos])
        self.drpivot.set_ydata(self.ax.get_ybound())


class AvProc1D:
    "Detailed 1D NMR Processing"

    def __init__(self, filename=""):
        print('WARNING this tool is not functional/tested yet')
        self.wfile = widgets.Text(
            description='File to process', layout=Layout(width='80%'), value=filename)
        self.wapod = widgets.Dropdown(
            options=['None', 'apod_sin (sine bell)', 'apod_em (Exponential)',
                     'apod_gm (Gaussian)', 'gaussenh (Gaussian Enhacement)', 'kaiser'],
            value='apod_sin (sine bell)',
            description='Apodisation')
        self.wpapod_Hz = widgets.FloatText(
            value=1.0,
            min=0,  # max exponent of base
            max=30,  # min exponent of base
            description='Width in Hz',
            layout=Layout(width='15%'),
            disabled=True)
        self.wpapod_enh = widgets.FloatText(
            value=2.0,
            min=0.0,  # max exponent of base
            max=5.0,  # min exponent of base
            description='strength',
            layout=Layout(width='15%'),
            step=1,
            disabled=True)
        self.wpapod_sin = widgets.FloatText(
            value=0.0,
            min=0,  # max exponent of base
            max=0.5,  # min exponent of base
            description='bell shape',
            layout=Layout(width='15%'),
            step=0.01,
            tooltip='value is the maximum of the bell, 0 is pure cosine, 0.5 is pure sine',
            disabled=False)
        self.wzf = widgets.Dropdown(
            options=[0, 1, 2, 4, 8],
            value=1,
            description='Zero-Filling')
        self.wphase0 = widgets.FloatText(
            value=0, description='Phase : P0', layout=Layout(width='20%'), disabled=True)
        self.wphase1 = widgets.FloatText(
            value=0, description='P1', layout=Layout(width='20%'), disabled=True)
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
        fi, ax = plt.subplots()
        self.ax = ax
        if os.path.exists(filename):
            self.load()
            # self.data.set_unit('sec')
            self.display()

    def apod_select(self, e):
        test = self.wapod.value.split()[0]
        self.wpapod_sin.disabled = True
        self.wpapod_Hz.disabled = True
        self.wpapod_enh.disabled = True
        if test == "apod_sin":
            self.wpapod_sin.disabled = False
        if test in ('apod_em', 'apod_gm', 'gaussenh'):
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
            todo = 'self.data.apod_sin(%f)' % (self.wpapod_sin.value,)
        elif func in ('apod_em', 'apod_gm'):
            todo = 'self.data.%s(%f)' % (func, self.wpapod_Hz.value)
        elif func == 'gaussenh':
            todo = 'self.data.gaussenh(%f,enhancement=%f)' % (
                self.wpapod_Hz.value, self.wpapod_enh.value)
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
            self.wphase0.value = round(self.data.axis1.P0, 1)
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
                  HBox([self.wapod, self.wpapod_sin, self.wpapod_Hz,
                        self.wpapod_enh, self.bapod]),
                  self.wzf,
                  HBox([self.wapmin, self.wphase0, self.wphase1]),
                  #                self.wbcorr,
                  self.bdoit]))


class NMRPeaker1D(Show1D):
    """
    a peak-picker for NMR experiments
    """
    # self.peaks : the defined peaklist, copyied in and out of data
    # self.temppk : the last computed pklist

    def __init__(self, data, figsize=None, show=True):
        self.data = data.real()
        try:
            if self.data.peaks is None:  # might happen
                self.peaks = Peaks.Peak1DList(
                    source=self.data)   # empty peak list
            else:
                self.peaks = self.data.peaks
        except AttributeError:
            self.peaks = Peaks.Peak1DList(source=self.data)   # empty peak list
        self.temppk = Peaks.Peak1DList(source=self.data)      # empty peak list
        self.thresh = widgets.FloatLogSlider(value=20.0,
                                             min=-1, max=2.0, base=10, step=0.01, layout=Layout(width='30%'),
                                             continuous_update=True, readout=True, readout_format='.2f')
        try:
            self.thresh.value = 100*self.data.peaks.threshold / \
                self.data.absmax  # if already peak picked
        except:
            self.thresh.value = 20.0
        self.thresh.observe(self.pickpeak)
        self.peak_mode = widgets.Dropdown(
            options=['marker', 'bar'], value='marker', description='show as')
        self.peak_mode.observe(self.ob)
        self.out = Output(layout={'border': '1px solid red'})
#        self.done = widgets.Button(description="Done", button_style='success')
#        self.done.on_click(self.on_done)
        self.badd = widgets.Button(
            description="Add", button_style='success', layout=space('80px'))
        self.badd.on_click(self.on_add)
        self.brem = widgets.Button(
            description="Rem", button_style='warning', layout=space('80px'))
        self.brem.on_click(self.on_rem)
        self.cancel = widgets.Button(
            description="Exit", button_style='warning', layout=space('80px'), tooltip='Exit without corrections')
        self.cancel.on_click(self.on_cancel)
        self.selval = widgets.FloatText(
            value=0.0, description='', layout=Layout(width='10%'), step=0.001, disabled=True)
        self.newval = widgets.FloatText(
            value=0.0, description='calibration', layout=Layout(width='20%'), step=0.001, disabled=True)
        self.setcalib = widgets.Button(description="Set", layout=Layout(width='10%'),
                                       button_style='success', tooltip='Set spectrum calibration')
        self.setcalib.on_click(self.on_setcalib)
        super().__init__(data, figsize=figsize, show=False)

        def on_press(event):
            if event.button == 3:       # right-click
                v = event.xdata
                # store position in index (peaks are internally in index)
                iv = self.data.axis1.ptoi(v)
                distclose = np.inf     # search closest peak

                pclose = 0.0
                for p in self.data.peaks:
                    if abs(p.pos-iv) < distclose:
                        pclose = p.pos
                        distclose = abs(p.pos-iv)
                self.selval.value = self.data.axis1.itop(pclose)  # back to ppm
                for w in (self.selval, self.newval):
                    w.disabled = False
        self.cids_press = self.fig.canvas.mpl_connect('button_press_event', on_press)
        # redefine Box
        orig = self.children
        self.tabs = Tab()
        self.tabs.children = [
            VBox([
                HBox([Label('threshold - % largest signal'),
                      self.thresh, self.badd, self.brem, self.peak_mode]),
                HBox(
                    [VBox([self.blank, self.reset, self.scale, self.done]), self.fig.canvas])
            ]),
            VBox([
                HBox([Label('Select a peak with mouse (or value) and set calibrated values',layout=Layout(width="40%")),
                      Label('selection (in ppm)'),self.selval, self.newval, self.setcalib]),
                HBox(
                    [VBox([self.blank, self.reset, self.scale, self.done]), self.fig.canvas])
            ]),
            self.out]
        self.tabs.set_title(0, 'Peak Picker')
        self.tabs.set_title(1, 'calibration')
        self.tabs.set_title(2, 'Peak Table')

        self.children = [VBox([HBox([self.cancel]), self.tabs])]

        self.pp()
        self.draw()

    def on_add(self, b):
        self.peaks.pkadd(self.temppk)
        self.peaks = Peaks.peak_aggreg(self.peaks, distance=1.0)
        self.peaks.source = self.data
        self.temppk = Peaks.Peak1DList(source=self.data)
        self.draw()

    def on_rem(self, b):
        (up, down) = self.ax.get_xbound()
        iup = self.data.axis1.ptoi(up)
        idown = self.data.axis1.ptoi(down)
        iup, idown = (max(iup, idown), min(iup, idown))
        to_rem = []
        for pk in self.peaks:
            if pk.pos < iup and pk.pos > idown:
                to_rem.append(pk)
        for pk in to_rem:
            self.peaks.remove(pk)
        self.draw()

    def on_cancel(self, b):
        self.close()
        del self.data.peaks
        print("no Peak-Picking done")

    def on_reset(self, b=None):
        self.thresh.value = 20.0
        super().on_reset()

    def on_done(self, b):
        self.temppk = Peaks.Peak1DList()  # clear temp peaks
        self.close()
        # new figure
        self.draw()
        # and display
        display(self.fig)
        display(self.out)
        self.data.peaks = self.peaks  # and copy

    def on_setcalib(self, e):
        off = self.selval.value-self.newval.value
        # off is in ppm, axis1.offset is in Hz
        self.data.axis1.offset -= off*self.data.axis1.frequency
        self.selval.value = self.newval.value
#        self.pp()
        self.peaks.pos2label()
        self.temppk.pos2label()
        self.draw()

    def pkprint(self, event):
        self.out.clear_output(wait=True)
        with self.out:
            if len(self.temppk) > 0:
                display(HTML("<p style=color:red> Transient peak list </p>"))
                self.data.peaks = self.temppk
                display(HTML(self.data.pk2pandas().to_html()))
            if len(self.peaks) > 0:
                display(HTML("<p style=color:blue> Defined peak list </p>"))
                self.data.peaks = self.peaks
                display(HTML(self.data.pk2pandas().to_html()))
            else:
                display(HTML("<p style=color:blue> Defined peak list Empty</p>"))

    def _pkprint(self, event):
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
            l = "%.3f\t%.1f\t%.2f" % (ppm, pk.intens/intmax, width)
            text.append(l)
        return "\n".join(text)

    def ob(self, event):
        if event['name'] == 'value':
            self.disp()

    def pickpeak(self, event):
        "interactive wrapper to pp"
        if event['name'] == 'value':
            self.disp()
            self.pp()

    @debounce(1.0)
    def pp(self):
        "do the peak-picking calling pp().centroid() within the current zoom"
        th = self.data.absmax*self.thresh.value/100
        zm = self.ax.get_xbound()
        self.data.set_unit('ppm').peakpick(
            threshold=th, verbose=False, zoom=zm).centroid()
        self.temppk = self.data.peaks
        self.data.peaks = None
        self.draw()
        self.ax.annotate('%d peaks detected' % len(
            self.data.peaks), (0.05, 0.95), xycoords='figure fraction')
        self.pkprint({'name': 'value'})

    def draw(self):
        "interactive wrapper to peakpick"
        super().draw()
        self.xb = self.ax.get_xbound()
        x = [self.data.axis1.itoc(z) for z in (0, self.data.size1)]
        y = [self.data.absmax*self.thresh.value/100]*2
        self.threshline = self.ax.plot(x, y, ':r')[0]
        if True:  # try:    # Pb should have been fixed
            self.temppk.display(peak_label=False, peak_mode=self.peak_mode.value,
                                f=self.data.axis1.itoc, figure=self.ax, color='red')
            self.peaks.display(peak_label=False, peak_mode=self.peak_mode.value,
                               f=self.data.axis1.itoc, color='blue', figure=self.ax)
        else:  # except:
            raise NotImplementedError("should never append")
        self.temppk.display(
            peak_label=True, peak_mode=self.peak_mode.value, color='red', figure=self.ax)
        self.ax.set_xbound(*self.xb)
        # self.ax.set_ylim(ymax=self.data.absmax/self.scale.value)
        # send pseudo event to display peak table
        self.pkprint({'name': 'value'})
        self.disp()

    def disp(self):
        y = [self.data.absmax*self.thresh.value/100]*2
        self.threshline.set_ydata(y)
        super().disp()


class NMRIntegrate(Show1D):
    "an integrator for NMR experiments"

    def __init__(self, data, figsize=None, show=True):
        try:
            self.Integ = data.integrals
        except:
            # initialize with empty list
            self.Integ = Integrals(data, compute=False)
        try:
            self.peaks_reserved = data.peaks
        except AttributeError:
            #print('no peaks')
            self.peaks_reserved = None
        self.idisp = 0
        self.idraw = 0
        super().__init__(data, figsize=figsize, show=show)
        self.thresh = widgets.FloatLogSlider(description='sensitivity', value=10.0,
                                             min=-1, max=2.0, base=10, step=0.01, layout=Layout(width='30%'),
                                             continuous_update=HEAVY, readout=True, readout_format='.1f',
                                             tooltip='sensitivity to weak signals')
        self.bias = widgets.FloatSlider(
            description='bias', layout=Layout(width='20%'),
            value=0.0, min=-10.0, max=10.0, step=0.1,
            continuous_update=HEAVY, readout=True, readout_format='.1f')
        self.sep = widgets.FloatSlider(
            description='separation', layout=Layout(width='30%'),
            value=3.0, min=0.0, max=20.0, step=0.1,
            continuous_update=HEAVY, readout=True, readout_format='.1f')
        self.wings = widgets.FloatSlider(
            description='extension', layout=Layout(width='30%'),
            value=5.0, min=0.5, max=20.0, step=0.1,
            continuous_update=HEAVY, readout=True, readout_format='.1f')
        for w in (self.bias, self.sep, self.wings):
            w.observe(self.integrate)
        self.thresh.observe(self.peak_and_integrate)
        self.cancel = widgets.Button(
            description="Exit", button_style='warning', layout=self.blay, tooltip='Exit without corrections')
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

        self.entry = widgets.IntText(value=0, description='Entry', min=0, layout=Layout(width='25%'),
                                     tooltip='Index of the calibrating integral')
        self.value = widgets.FloatText(value=100, description='Value', layout=Layout(width='25%'),
                                       tooltip='Value of the calibrating integral')
        self.set = widgets.Button(
            description="Set", button_style='success', layout=Layout(width='10%'))
        self.set.on_click(self.set_value)
        self.out = Output(layout={'border': '1px solid red'})
        # redefine children
        orig = self.children
        self.tabs = Tab()
        self.tabs.children = [
            VBox([
                HBox([widgets.HTML("Use buttons to add and remove integrals in the current zoom window&nbsp;&nbsp;"),
                      self.badd, self.brem]),
                HBox(orig),
            ]),
            VBox([HBox([widgets.HTML("Define integral shapes using the sliders below (erases the current integrals)&nbsp;&nbsp;"), self.bauto]),
                  HBox([self.thresh,  self.sep, self.wings]),
                  HBox(orig),
                  ]),
            VBox([HBox([Label('Choose an integral for calibration'),
                        self.entry, self.value, self.set]),
                  self.out])
        ]
        self.tabs.set_title(0, 'Manual integration')
        self.tabs.set_title(1, 'Automatic')
        self.tabs.set_title(2, 'Integral Table & Calibration')
        self.children = [VBox([self.cancel, self.tabs])]
        self.draw()
        self.print(None)

    def on_cancel(self, b):
        self.close()
        print("No integration")

    def on_done(self, b):
        self.draw()
        display(self.fig)
        display(self.out)
        self.data.integrals = self.Integ        # copy integrals
        if self.peaks_reserved:                # restore peaks
            self.data.peaks = self.peaks_reserved
        else:
            try:
                del(self.data.peaks)
            except AttributeError:    # if no peaks
                pass
        self.close()

    def on_add(self, b):
        start, end = self.ax.get_xbound()
        self.Integ.append(Integralitem(self.data.axis1.ptoi(
            start), self.data.axis1.ptoi(end), [], 0.0))
        self.Integ.zonestocurves()
        self.draw()
        self.print(None)

    def on_rem(self, b):
        start, end = self.ax.get_xbound()
        start, end = self.data.axis1.ptoi(start), self.data.axis1.ptoi(end)
        # -1,+1 needed sometimes...
        start, end = (min(start, end)-1, max(start, end)+1)
        to_rem = []
        for ii in self.Integ:
            if ii.start > start and ii.end < end:
                to_rem.append(ii)
        for ii in to_rem:
            self.Integ.remove(ii)
        self.draw()
        self.print(None)

    def set_value(self, b):
        self.Integ.recalibrate(self.entry.value, self.value.value)
        self.draw()
        self.print(None)

    def print(self, event):
        self.out.clear_output()
        self.Integ.sort(key=lambda x: x.start)   # sort integrals
        with self.out:
            display(HTML(self.Integ.to_pandas().to_html()))

    def peak_and_integrate(self, event):
        self.data.pp(threshold=self.data.absmax*self.thresh.value/100,
                     verbose=False, zoom=self.ax.get_xbound()).centroid()
        if len(self.data.peaks) > 0:
            self.int()

    def integrate(self, event):
        #"integrate from event"
        # if event['name']=='value':
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
        self.Integ += Integrals(self.data, separation=self.sep.value, wings=self.wings.value,
                                bias=self.data.absmax*self.bias.value/100)
        self.Integ.calibrate(calibration=calib)
        self.print(None)
        self.draw()

    def ob(self, event):
        if event['name'] == 'value':
            self.draw()

    def disp(self):
        self.ax.set_ybound(self.yb0/self.scale.value)
        self.idisp += 1

    def draw(self):
        "refresh display from event - if zoom is True, display only in xbound"
        # self.yb = self.ax.get_ybound()
        super().draw()
        self.idraw += 1
        self.Integ.display(label=True, figure=self.ax,
                           labelyposition=None, regions=False, zoom=self.ax.get_xbound())
#        self.ax.set_xbound(*self.xb)
        # self.ax.set_ybound(self.yb)


# if __name__ == '__main__':
#    unittest.main()
