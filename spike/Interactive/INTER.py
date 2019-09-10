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
import matplotlib.gridspec as gridspec
from ipywidgets import interact, interactive, fixed, interact_manual, Layout, HBox, VBox, Label, Output, Button
import ipywidgets as widgets
from IPython.display import display, HTML, Javascript
import numpy as np

from ..File.BrukerNMR import Import_1D
from .. import NPKData
try:
    import spike.plugins.bcorr as bcorr
except:
    print('Baseline correction plugins not installed !')

# REACTIVE modify callback behaviour
# True is good for inline mode / False is better for notebook mode
REACTIVE = True
HEAVY = False

def hidecode(initial='show', message=True):
    """
    this func adds a button to hide/show the code on a jupyter page
    initial is either 'show' or 'hide'
    see: https://stackoverflow.com/questions/27934885/how-to-hide-code-from-cells-in-ipython-notebook-visualized-with-nbviewer/28073228#28073228
    """
    from IPython.display import display, HTML, Markdown
    if initial == 'show':
        init = 'false'
    elif initial == 'hide':
        init = 'true'
    if message:
        msg = "<i>usefull to show/print a clean screen when processing is finished</i>"
    else:
        msg = ""
    display(HTML('''
<style>hr {height: 2px; border: 0;border-top: 1px solid #ccc;margin: 1em 0;padding: 0; }</style>
<script>
code_show=%s; 
function code_toggle() {
 if (code_show){
 $('div.input').hide();
 } else {
 $('div.input').show();
 }
 code_show = !code_show
} 
$( document ).ready(code_toggle);
</script>
<form action="javascript:code_toggle()">
<input type="submit" style="border:1px solid black; background-color:#DDD" value="hide/show the python code.">
%s
</form>'''%(init, msg)))

def jsalert(msg):
    "send a javascript alert"
    display(Javascript("alert('%s')"%msg))

class FileChooser:
    """a simple file chooser for Jupyter - obsolete - not used"""
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

class Show1D(object):
    """
    An interactive display, 1D NMR
        Show1D(spectrum)
    to be developped for peaks and integrals
    """
    def __init__(self, data, title=None, figsize=None):
        self.data = data
        self.title = title
        if figsize is None:
            figsize = (10,5)
        self.done = Button(description="Done", button_style='success',layout=Layout(width='10%'))
        self.scale = widgets.FloatSlider(description='scale:', value=1.0, min=1.0, max=100, step=0.1,
                            layout=Layout(width='50%'), continuous_update=REACTIVE)
        for widg in (self.scale,):
            widg.observe(self.ob)
        fi,ax = plt.subplots(figsize=figsize)
        self.ax = ax
        self.data.display(figure=self.ax, title=self.title)
        self.xb = self.ax.get_xbound()
        self.done.on_click(self.on_done)
        self.Box = HBox([self.scale, self.done])
        display( self.Box )
    def on_done(self, b):
        self.close()
    def close(self):
        for w in [ self.Box, self.done, self.scale]:
            w.close()
    def ob(self, event):
        "observe events and display"
        if event['name']=='value':
            self.disp()
    def disp(self):
        self.xb = self.ax.get_xbound()
        self.ax.clear()
        self.data.display(scale=self.scale.value, new_fig=False, figure=self.ax, title=self.title)
        self.ax.set_xbound(self.xb)

class baseline1D(Show1D):
    def __init__(self, data, figsize=None):
        super().__init__( data, figsize=figsize)
        self.data.real()
        ppos = self.data.axis1.itop(self.data.axis1.size/2)
        #self.ax.plot([ppos,ppos], self.ax.get_ybound())
        self.select = widgets.FloatSlider(description='select:', min=0.0, max=1.0, step=0.002, layout=Layout(width='100%'),
                                         value=0.5, readout=False, continuous_update=REACTIVE)
        self.selfine = widgets.FloatSlider(description='fine:', min=0.0, max=1.0, step=0.01, layout=Layout(width='70%'),
                                         value=0.5, readout=False, continuous_update=REACTIVE)
        self.bsl_points = []
        for w in [self.select, self.selfine]:
            w.observe(self.ob)
        bsize = '15%'
#        self.apply = widgets.Button(description="Apply", button_style='success', layout=Layout(width=2*bsize))
#        self.apply.on_click(self.on_apply)
        self.done.description = 'Apply'
        self.cancel = widgets.Button(description="Abort", button_style='warning', layout=Layout(width=2*bsize))
        self.cancel.on_click(self.on_cancel)
        self.Box.add_class("hidden")
        self.auto = widgets.Button(description="auto", button_style='success', layout=Layout(width=bsize))
        self.auto.on_click(self.on_auto)
        self.set = widgets.Button(description="add", button_style='success', layout=Layout(width=bsize))
        self.set.on_click(self.on_set)
        self.unset = widgets.Button(description="rem", button_style='warning', layout=Layout(width=bsize))
        self.unset.on_click(self.on_unset)
        self.toshow = widgets.Dropdown( options=['baseline', 'corrected', 'points'],  description='Display:')
        self.toshow.observe(self.ob)
        self.Box = HBox([self.scale, self.done,self.cancel])
        display(self.Box)
        display( HBox( [
                        VBox([self.select, self.selfine], layout=Layout(width='50%')),
                        VBox([  HBox([self.set, self.auto, self.toshow], layout=Layout(width='100%')),
                                self.unset
                              ], layout=Layout(width='50%')) 
                        ], layout=Layout(width='100%'))
                )
        self.disp()
    def close(self):
        for w in [ self.select, self.selfine, self.auto, self.set, self.unset, self.cancel, self.toshow]:
            w.close()
        super().close()
    def on_done(self, e):
        print('Applied correction:\n',self.bsl_points)
        self.data.set_buffer( self.data.get_buffer() - self.correction() )
        super().disp()
        self.close()
    def on_auto(self, e):
        "automatically set baseline points"
        self.bsl_points = [self.data.axis1.itop(x) for x in bcorr.autopoints(self.data)]
        self.disp()
    def on_cancel(self, e):
        self.close()
        self.ax.clear()
        self.data.display(new_fig=False, figure=self.ax)
    def on_set(self, e):
        "add baseline points at selector"
        self.bsl_points.append( self.selpos() )
        self.selfine.value = 0.5
        self.disp()
    def on_unset(self, e):
        "remove baseline points closest from selector"
        here = self.selpos()
        distclose = np.inf
        pclose = np.NaN
        for i,p in enumerate(self.bsl_points):
            if abs(p-here)< distclose:
                pclose = p
                distclose = abs(p-here)
        self.bsl_points.remove(pclose)
        self.disp()
    def selpos(self):
        "returns selector pos in ppm"
        pos = self.select.value + self.selfine.value/100
        return (1-pos)*(self.data.axis1.itop(0)-self.data.axis1.itop(self.data.size1))+self.data.axis1.itop(self.data.size1)
    def correction(self):
        "returns the correction to apply as a numpy array"
        ibsl_points = [int(self.data.axis1.ptoi(x)) for x in self.bsl_points]
        x = np.arange(self.data.size1)
        if len(self.bsl_points) == 1 :
            value = self.data[ibsl_points[0]] * np.ones( self.data.size1 )
        elif len(self.bsl_points) < 4 :
            corr = bcorr._linear_interpolate(self.data.get_buffer(), ibsl_points)
            value = corr(x)
        else:
            corr = bcorr._spline_interpolate(self.data.get_buffer(), ibsl_points)
            value = corr(x)
        return value
    def corrected(self):
        value = self.data.copy()
        value.set_buffer( value.get_buffer() - self.correction() )
        return value
    def disp(self):
        self.xb = self.ax.get_xbound()
        super().disp()
        if len(self.bsl_points)>0:
            if self.toshow.value == 'baseline':
                ( self.data.copy()-self.corrected() ).display(new_fig=False, figure=self.ax, color='r')
            elif self.toshow.value == 'corrected':
                self.ax.clear()
                self.corrected().display(new_fig=False, figure=self.ax, color='r', scale=self.scale.value)
        ppos = self.selpos()
        self.ax.plot([ppos,ppos], self.ax.get_ybound())# [0,self.data.absmax/self.scale.value])
        y = []
        for point in self.bsl_points:
            if self.toshow.value == 'corrected':
                y.append(0)
            else:
                loc = int(self.data.axis1.ptoi(point))
                y.append(self.data[loc])
        if y:
            self.ax.scatter(self.bsl_points, y, c='r', marker='o')
        self.ax.set_xbound(self.xb)

class Show1Dplus(object):
    """
    An interactive display, 1D NMR
        Show1D(spectrum)
    to be developped for peaks and integrals
    """
    def __init__(self, data, title=None):
        self.data = data
        self.title = title
        self.done = Button(description="Done", button_style='warning',layout=Layout(width='10%'))
        self.scale = widgets.FloatSlider(description='scale:', value=1.0, min=1.0, max=100, step=0.1,
                            layout=Layout(width='30%'), continuous_update=REACTIVE)
        self.scaleint = widgets.FloatSlider(value=0.5, min=0.1, max=10, step=0.05,
                            layout=Layout(width='20%'), continuous_update=REACTIVE)
        self.offset = widgets.FloatSlider(value=0.3, min=0.0, max=1.0, step=0.01,
                            layout=Layout(width='20%'), continuous_update=REACTIVE)
        self.zoom = widgets.FloatRangeSlider(value=[0, 100],
            min=0, max=100.0, step=0.1, layout=Layout(width='60%'), description='zoom (%):',
            continuous_update=REACTIVE, readout=True, readout_format='.1f')
        self.peaks = widgets.Checkbox(value=False, layout=Layout(width='15%'))
        self.integ = widgets.Checkbox(value=False, layout=Layout(width='15%'))

        for widg in (self.scale, self.zoom, self.scaleint, self.offset, self.peaks, self.integ):
            widg.observe(self.ob)
        fi,ax = plt.subplots()
        self.ax = ax
        self.data.display(figure=self.ax, zoom=self.zm, title=self.title)
        self.done.on_click(self.on_done)
        self.Box = VBox([self.zoom,
            HBox([self.scale,Label('Integral scale:'),self.scaleint,Label('offset:'),self.offset]),
            HBox([Label('Show Peaks'),self.peaks,Label('integrals'),self.integ,self.done])])
        display( self.Box )
    def on_done(self, b):
        self.Box.close()
    @property
    def zm(self):
        "returns the zoom window in ppm"
        z = self.zoom.value
        left=self.data.axis1.itop(z[0]*self.data.size1/100)
        right=self.data.axis1.itop(z[1]*self.data.size1/100)
        return (left,right)
    def ob(self, event):
        "observe events and display"
        if event['name']!='value':
            return
        self.ax.clear()
        self.data.display(scale=self.scale.value, new_fig=False, figure=self.ax, zoom=self.zm, title=self.title)
        if self.integ.value:
            try:
                self.data.display_integral(label=True, integscale=self.scaleint.value, integoff=self.offset.value, figure=self.ax, zoom=self.zm)
            except:
                print('no or wrong integrals')
                pass
        if self.peaks.value:
            try:
                self.data.display_peaks(peak_label=True, figure=self.ax, zoom=self.zm, scale=self.scale.value)
            except:
                print('no or wrong peaklist')
                pass
        self.ax.set_xlim(left=self.zm[0], right=self.zm[1])

class baseline2D_F2(baseline1D):
    def __init__(self, data, figsize=None):
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
        self.posview = widgets.Checkbox(value=True,description='Positive')
        self.negview = widgets.Checkbox(value=False,description='Negative')
        for w in (self.scale, self.posview, self.negview):
            w.observe(self.ob)
        grid = {'height_ratios':[1,4],'hspace':0,'wspace':0}
        if self.isDOSY:
            fsize = (10,5)
            grid['width_ratios']=[7,1]
        else:
            fsize = (8,8)
            grid['width_ratios']=[4,1]
#        fig, self.axarr = plt.subplots(2, 1, sharex=True, figsize=fsize, gridspec_kw=grid)
        fig = plt.figure(figsize=fsize, constrained_layout=False)
        spec2 = gridspec.GridSpec(ncols=2, nrows=2, figure=fig, **grid)
        axarr = np.empty((2,2), dtype=object)
        axarr[0,0] = fig.add_subplot(spec2[0, 0])
        axarr[1,0] = fig.add_subplot(spec2[1, 0],sharex=axarr[0, 0])
        axarr[1,1] = fig.add_subplot(spec2[1, 1],sharey=axarr[1, 0])
        self.top_ax = axarr[0,0]
        self.spec_ax = axarr[1,0]
        self.side_ax = axarr[1,1]
        self.Box = HBox( [self.scale, self.posview, self.negview])
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
        if self.posview.value:
            self.data.display(scale=self.scale.value, new_fig=False, figure=self.spec_ax)
        if self.negview.value:
#            m = -self.data.absmax/self.scale.value
#            level = sorted([m*0.05, m*0.1, m*0.25, m*0.5])  # correction for matplotlib 1.5.1
            self.data.display(scale=-self.scale.value, new_fig=False,
                figure=self.spec_ax, mpldic={'cmap':'Wistia'})
        self.spec_ax.set_xbound(xb)
        self.spec_ax.set_ybound(yb)

        # if self.isDOSY:   # does not seem needed
        #     self.side_ax.set_yscale('log')

class Show2Dplus(object):
    """
    A VERY preliminary minimimum interactive display for 2D NMR
        Show2D(spectrum)
    """
    def __init__(self, data, title=None):
        self.data = data
        self.isDOSY = isinstance(self.data.axis1, NPKData.LaplaceAxis)
        self.title = title
        self.scale = widgets.FloatLogSlider(description='scale:', value=1.0, min=-1, max=2,  base=10, step=0.01,
                            layout=Layout(width='80%'), continuous_update=HEAVY)
        self.phased = widgets.Checkbox(value=False, layout=Layout(width='15%'))
        self.peaks = widgets.Checkbox(value=False, layout=Layout(width='15%'))
        self.fullzoom()
        for widg in (self.scale, self.peaks, self.phased):
            widg.observe(self.ob)
        fi,ax = plt.subplots()
        self.ax = ax
        self.disp(new=True)
        self.Box = VBox([self.scale,
                        HBox([Label("Phase sensitive"), self.phased, Label("Showpeaks"),self.peaks]),
                        self.zoom_box()
                        ]) 
        display( self.Box )
    def fullzoom(self):
        ref = self.data
        if self.isDOSY:
            self._zoom =  (ref.axis1.itod(ref.size1), ref.axis1.itod(0), ref.axis2.itop(ref.size2),ref.axis2.itop(0))
        else:
            self._zoom =  (ref.axis1.itop(ref.size1), ref.axis1.itop(0), ref.axis2.itop(ref.size2),ref.axis2.itop(0))
    def zoom_box(self):
        wf = widgets.BoundedFloatText
        ref = self.data
        style = {'description_width': 'initial'}
        lay = Layout(width='80px', height='30px')
        if self.isDOSY:
            self.z1l = wf( value=self._zoom[0], max=ref.axis1.itod(ref.size1), min=ref.axis1.itod(0),
                description='Diff:', style=style, layout=lay)
            self.z1h = wf( value=self._zoom[1], max=ref.axis1.itod(ref.size1), min=ref.axis1.itod(0),
                description='..', style=style, layout=lay)
        else:
            self.z1l = wf( value=self._zoom[0], max=ref.axis1.itop(0), min=ref.axis1.itop(ref.size1),
                description='F1:', style=style, layout=lay)
            self.z1h = wf( value=self._zoom[1], max=ref.axis1.itop(0), min=ref.axis1.itop(ref.size1),
                description='..', style=style, layout=lay)
        self.z2l = wf( value=self._zoom[2], max=ref.axis2.itop(0), min=ref.axis2.itop(ref.size2),
            description='F2:', style=style, layout=lay)
        self.z2h = wf( value=self._zoom[3], max=ref.axis2.itop(0), min=ref.axis2.itop(ref.size2),
            description='..', style=style, layout=lay)
        def zupdate(b):
            self._zoom = (self.z1l.value, self.z1h.value, self.z2l.value, self.z2h.value)
        def set(b):
            self.disp()
        def reset(b):
            self.fullzoom()
            self.disp()
        for z in (self.z1l, self.z1h, self.z2l, self.z2h):
            z.observe(zupdate)
        self.bset = Button(description="Set",button_style='success')
        self.breset = Button(description="Reset",button_style='success')
        self.bset.on_click(set)
        self.breset.on_click(reset)
        #self.zmupd=('b_zupdate', 'Update', zupdate, layout=Layout(width='300px'), tooltip="Set zoom to values")
        return VBox([ Label('Zoom Window (in ppm)'),
                    HBox([self.z1l, self.z1h, self.z2l, self.z2h, self.bset, self.breset])])
    def on_done(self, b):
        self.Box.close()
    def ob(self, event):
        "observe events and display"
        if event['name']!='value':
            return
        self.disp()
    def disp(self,new=False):
        print(self._zoom)
        if not new:
            self.ax.clear()
        self.data.display(scale=self.scale.value, new_fig=new, figure=self.ax, zoom=self._zoom, title=self.title)
        if self.phased.value:
            self.data.display(scale=-self.scale.value, new_fig=new, figure=self.ax, zoom=self._zoom, color='red')
        self.ax.set_xlim(xmin=self._zoom[2], xmax=self._zoom[3])
        self.ax.set_ylim(ymin=self._zoom[1], ymax=self._zoom[0])

class SimpleZoom(HBox):
    def __init__(self, **kwargs):
        self.zoomw = widgets.FloatRangeSlider(value=[0, 100],
        min=0, max=100.0, step=0.1, layout=Layout(width='60%'), description='zoom (%):',
        continuous_update=REACTIVE, readout=False, readout_format='.1f',)
        self.zoomdisplay = widgets.Label('- / -')
        self.zoomfull = widgets.Button(description="Full", button_style='')
        self.zoomfull.on_click(self.zfull)
        self.Box = [HBox([self.zoomw,self.zoomdisplay,self.zoomfull])]
        self.value = [0,100]
        super().__init__( children=self.Box, layout=Layout(width='auto'), **kwargs)
        self.zoomw.observe(self.ob)
    def ob(self, event):
        "observe changes and set values"
        if event['name']=='value':
            self.value = self.zoomw.value
    @property
    def zm(self):
        "returns the zoom window in ppm"
        z = self.zoomw.value
        left = self.data.axis1.itop(z[0]*self.data.size1/100)
        right = self.data.axis1.itop(z[1]*self.data.size1/100)
        self.zoomdisplay.value = "%.2f / %.2f"%(left,right)
        return (left,right)
    def zfull(self, b):
        self.zoomw.value = [0,100]

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
        self.pivot = widgets.FloatSlider(description='pivot:', min=0.0, max=1.0, step=0.01, layout=Layout(width='30%'),
                                         value=0.5, readout=False, continuous_update=REACTIVE)
        for w in [self.p0, self.p1, self.scale, self.pivot]:
            w.observe(self.ob)
        self.button = widgets.Button(description="Ok", button_style='success')
        self.button.on_click(self.on_Apply)
        self.cancel = widgets.Button(description="Cancel", button_style='warning')
        self.cancel.on_click(self.on_cancel)
        display( VBox([ HBox([self.button, self.cancel]),
                        HBox([self.scale, self.pivot]),
                        self.p0, self.p1]) )
        fi,ax = plt.subplots()
        self.ax = ax
        self.data = data
        self.data.display(figure=self.ax)
        ppos = self.data.axis1.itop(self.data.axis1.size/2)
        self.ax.plot([ppos,ppos], self.ax.get_ybound())
        self.xb = self.ax.get_xbound()  # initialize zoom
    def close_all(self):
        for w in [self.p0, self.p1, self.scale, self.button, self.cancel, self.pivot]:
            w.close()
    def on_cancel(self, b):
        self.close_all()
        self.ax.clear()
        self.data.display(new_fig=False, figure=self.ax)
    def on_Apply(self, b):
        lp0, lp1 = self.ppivot() # get centered values
        if lp0 != 0 or lp1 != 0:
            print("Applied: phase(%.1f,  %.1f)"%(lp0, lp1))
            self.data.phase(lp0, lp1).display(new_fig=False, figure=self.ax)
        else:
            print("no applied phase")
        self.close_all()
#        self.ax.clear()
    def ppivot(self):
        "converts from pivot values to centered ones"
        pp = 1.0-self.pivot.value
        return (self.p0.value + (pp-0.5)*self.p1.value, self.p1.value)
    def ob(self, event):
        "observe changes and start phasing"
        if event['name']=='value':
            self.phase()
    def phase(self):
        "apply phase and display"
        self.xb = self.ax.get_xbound()   # get current zoom
        self.ax.clear()
        lp0, lp1 = self.ppivot() # get centered values
        self.data.copy().phase(lp0, lp1).display(scale=self.scale.value, new_fig=False, figure=self.ax)
        ppos = (1-self.pivot.value)*(self.data.axis1.itop(0)-self.data.axis1.itop(self.data.size1))+self.data.axis1.itop(self.data.size1)
        self.ax.plot([ppos,ppos], self.ax.get_ybound())
        self.ax.set_xbound( self.xb )
        #z = self.zm
        #self.ax.set_xlim(left=z[0], right=z[1] )

class Phaser2D(object):
    """
    An interactive phaser in 2D NMR

        Phaser2D(spec)

    best when in %matplotlib inline

    """
    def __init__(self, data):
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

class NMRPeaker(object):
    "a peak-picker for NMR experiments"
    def __init__(self, npkd):
        self.npkd = npkd.real()
        self.zoom = widgets.FloatRangeSlider(value=[0, 100],
            min=0, max=100.0, step=0.1, layout=Layout(width='60%'), description='zoom (%):',
            continuous_update=REACTIVE, readout=True, readout_format='.1f',)
        self.zoom.observe(self.display)
        self.thresh = widgets.FloatLogSlider(value=20.0,
            min=-1, max=2.0, base=10, step=0.01, layout=Layout(width='30%'),
            continuous_update=False, readout=True, readout_format='.1f')
        self.thresh.observe(self.pickpeak)
        self.peak_mode = widgets.Dropdown(options=['marker', 'bar'],value='marker',description='show as')
        self.peak_mode.observe(self.display)
        self.bprint = widgets.Button(description="Print", layout=Layout(width='10%'),
                button_style='success', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='Print to screen')
        self.bprint.on_click(self.pkprint)
        self.spec = Output(layout={'border': '1px solid black'})
        self.out = Output(layout={'border': '1px solid red'})
        self.Ok = widgets.Button(description="Ok", button_style='success')
        self.Ok.on_click(self.on_Apply)
        self.cancel = widgets.Button(description="Cancel", button_style='warning')
        self.cancel.on_click(self.on_cancel)
        self.Box = VBox([HBox([self.Ok, self.cancel]),
            self.zoom,
            HBox([Label('threshold - % largest signal'), self.thresh, self.peak_mode, self.bprint] )])
        display(self.Box )
        display(self.spec)
        display(self.out )
        with self.spec:
            self.fig, self.ax = plt.subplots()
            self.npkd.display(new_fig=False, figure=self.ax, zoom=self.zm)
    def on_cancel(self, b):
        self.Box.close()
        self.spec.close()
        self.out.close()
        del self.npkd.peaks
        self.npkd.display(new_fig=False, figure=self.ax)
    def on_Apply(self, b):
        pkm = self.peak_mode.value
        self.Box.close()
        self.spec.close()
        self.npkd.display()
        self.npkd.display_peaks(peak_mode=pkm,figure=self.npkd.mplfigure)
        self.npkd.mplfigure.annotate('%d peaks detected'%len(self.npkd.peaks) ,(0.05,0.95), xycoords='figure fraction')
    @property
    def zm(self):
        "returns the zoom window in ppm"
        z = self.zoom.value
        left=self.npkd.axis1.itop(z[0]*self.npkd.size1/100)
        right=self.npkd.axis1.itop(z[1]*self.npkd.size1/100)
        return (left,right)
    def pkprint(self,event):
        self.out.clear_output(wait=True)
        with self.out:
            display(HTML(self.npkd.pk2pandas().to_html()))
    def _pkprint(self,event):
        self.out.clear_output(wait=True)
        with self.out:
            print(self.pklist())
    def pklist(self):
        "creates peaklist for printing or exporting"
        text = ["ppm\tInt.(%)\twidth (Hz)"]
        data = self.npkd
        intmax = max(data.peaks.intens)/100
        for pk in data.peaks:
            ppm = data.axis1.itop(pk.pos)
            width = 2*pk.width*data.axis1.specwidth/data.size1
            #(data.axis1.itoh(pk.pos-pk.width) - data.axis1.itoh(pk.pos+pk.width))
            l = "%.3f\t%.1f\t%.2f"%(ppm, pk.intens/intmax, width)
            text.append(l)
        return "\n".join(text)
    def display(self, event):
        "interactive wrapper to peakpick"
        if event['name']=='value':
            self.ax.clear()
            self.npkd.display(new_fig=False, figure=self.ax, zoom=self.zm)
            try:
                self.npkd.display_peaks(peak_label=True, figure=self.ax, zoom=self.zm)
            except:
                pass
    def pickpeak(self, event):
        "interactive wrapper to peakpick"
        if event['name']=='value':
            self.pp()
    def pp(self):
        "do the peak-picking calling pp().centroid()"
        #self.spec.clear_output(wait=True)
        th = self.npkd.absmax*self.thresh.value/100
        self.npkd.set_unit('ppm').peakpick(threshold=th, verbose=False, zoom=self.zm).centroid()
        with self.spec:
            self.ax.clear()
            self.npkd.display(new_fig=False, figure=self.ax, zoom=self.zm)
            x = self.zm
            y = [self.npkd.peaks.threshold]*2
            self.ax.plot(x,y,':r')
            z = self.zoom.value
            left=self.npkd.axis1.itop(z[0]*self.npkd.size1/100)
            right=self.npkd.axis1.itop(z[1]*self.npkd.size1/100)
            self.npkd.display_peaks(peak_label=True, peak_mode=self.peak_mode.value,figure=self.ax, zoom=self.zm)
            self.ax.annotate('%d peaks detected'%len(self.npkd.peaks) ,(0.05,0.95), xycoords='figure fraction')
    def _pp(self):
        "do the peak-picking calling pp().centroid()"
        self.ax.clear()
        self.npkd.set_unit('m/z').peakpick(autothresh=self.thresh.value, verbose=False).centroid()
        self.npkd.display(new_fig=False, figure=self.ax)
        x = [self.npkd.axis1.lowmass, self.npkd.axis1.highmass]
        y = [self.npkd.peaks.threshold]*2
        self.ax.plot(x,y,':r')
        self.npkd.display_peaks(peak_label=True, figure=self.ax)
        self.ax.annotate('%d peaks detected'%len(self.npkd.peaks) ,(0.05,0.95), xycoords='figure fraction')

class NMRIntegrate(object):
    "an integrator for NMR experiments"
    def __init__(self, npkd):
        self.npkd = npkd.real().integrate()
        self.zoom = widgets.FloatRangeSlider(value=[0, 100],
            min=0, max=100.0, step=0.1, layout=Layout(width='60%'), description='zoom (%):',
            continuous_update=REACTIVE, readout=True, readout_format='.1f',)
        self.zoom.observe(self.display)
        self.bias = widgets.FloatSlider(
            description='bias', layout=Layout(width='30%'),
            value=0.0,min=-10.0, max=10.0, step=0.1,
            continuous_update=True, readout=True, readout_format='.1f')
        self.sep = widgets.FloatSlider(
            description='peak separation', layout=Layout(width='30%'),
            value=3.0,min=0.0, max=10.0, step=0.1,
            continuous_update=True, readout=False, readout_format='.1f')
        self.wings = widgets.FloatSlider(
            description='width', layout=Layout(width='30%'),
            value=5.0,min=0.5, max=10.0, step=0.1,
            continuous_update=True, readout=False, readout_format='.1f')
        self.bias.observe(self.integrate)
        self.sep.observe(self.integrate)
        self.wings.observe(self.integrate)
        self.Ok = widgets.Button(description="Ok", button_style='success')
        self.Ok.on_click(self.on_Apply)
        self.cancel = widgets.Button(description="Cancel", button_style='warning')
        self.cancel.on_click(self.on_cancel)
        self.bprint = widgets.Button(description="Print", layout=Layout(width='10%'),
                button_style='success', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='Print to screen')
        self.bprint.on_click(self.print)
        self.entry = widgets.IntText(value=0,description='Entry',min=0,layout=Layout(width='15%'))
        self.value = widgets.FloatText(value=100,description='Value',layout=Layout(width='15%'))
        self.set = widgets.Button(description="Set", button_style='success', layout=Layout(width='7%'))
        self.set.on_click(self.set_value)
        self.spec = Output(layout={'border': '1px solid black'})
        self.out = Output(layout={'border': '1px solid red'})
        self.Box = VBox([
            HBox([self.Ok, self.cancel]),
            self.zoom,
            HBox([self.bias, self.sep, self.wings,] ),
            HBox([self.entry, self.value, self.set, self.bprint] ),])
        display( self.Box)
        display(self.spec)
        display( self.out)
        with self.spec:
            self.fig, self.ax = plt.subplots()
        self.npkd.display(new_fig=False, figure=self.ax, zoom=self.zm)
        self.npkd.display_integral(label=True, zoom=self.zm, figure=self.ax, labelyposition=0.01, regions=False)
    def on_cancel(self, b):
        self.Box.close()
        self.spec.close()
        self.out.close()
        self.npkd.display()
    def on_Apply(self, b):
        self.Box.close()
        self.spec.close()
        self.npkd.display()
        self.npkd.display_integral(label=True, figure=self.npkd.mplfigure, labelyposition=0.01, regions=False)
    def set_value(self,b):
        self.npkd.integral_calibrate(self.entry.value, self.value.value)
        self.display({'name':'value'})
    @property
    def zm(self):
        "returns the zoom window in ppm"
        z = self.zoom.value
        left=self.npkd.axis1.itop(z[0]*self.npkd.size1/100)
        right=self.npkd.axis1.itop(z[1]*self.npkd.size1/100)
        return (left,right)
    def print(self,event):
        self.out.clear_output(wait=True)
        with self.out:
            display(HTML( self.npkd.integrals.to_pandas().to_html() ))
    def integrate(self, event):
        "integrate from event"
        if event['name']=='value':
            self.int()
    def int(self):
        "do the integration"
        self.npkd.set_unit('ppm').integrate(separation=self.sep.value, wings=self.wings.value,
            bias=self.npkd.absmax*self.bias.value/100)
        # self.value.value = 100.0
        # for i,val in enumerate(self.npkd.integvalues):
        #     if val == 100.0:
        #         self.entry.value = i
        #         break
        self.display({'name':'value'})
    def display(self, event):
        "refresh display from event"
        if event['name']=='value':
            self.ax.clear()
            self.npkd.display(new_fig=False, figure=self.ax, zoom=self.zm)
            try:
                self.npkd.display_integral(label=True, figure=self.ax, zoom=self.zm, labelyposition=0.01, regions=False)
            except:
                pass

#if __name__ == '__main__':
#    unittest.main()    