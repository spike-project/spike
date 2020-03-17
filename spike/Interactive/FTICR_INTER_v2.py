#!/usr/bin/env python 
# encoding: utf-8

"""
A tool to display FT-ICR data-sets

to be embedded in jupyter notebook

MAD Sept 2019

This version requires ipympl (see: https://github.com/matplotlib/jupyter-matplotlib )
and the notebook to be opened with %matplotlib widget
"""
import os.path as op
import tables
import matplotlib.pyplot as plt
from ipywidgets import interact, fixed, HBox, VBox, GridBox, Label, Layout, Output, Button
import ipywidgets as widgets
from IPython.display import display, Markdown, HTML, Image
import numpy as np

from .. import FTICR
from ..NPKData import flatten, parsezoom
from .ipyfilechooser import FileChooser
from ..FTMS import FTMSData
from ..FTICR import FTICRData

# REACTIVE modify callback behaviour
# True is good for inline mode / False is better for notebook mode
REACTIVE = True
HEAVY = False

SIZEMAX = 8*1024*1024       # largest zone to display
NbMaxDisplayPeaks = 1000      # maximum number of peaks to display at once

# TOOLS FOR 2D FTICR
class MR(object):
    "this class handles multiresolution datasets"
    def __init__(self, name, report=True, Debug=False):
        "name : filename of the msh5 multiresolution file"
        self.SIZEMAX = SIZEMAX
        self.name = name
        self.data = []                   # will contain all resolutions
        self.absmax = 0.0                # will accumulates global absmax
        self.load()
        if self.absmax == 0.0:
            self.compute_absmax()
        self.axis1 = self.data[0].axis1
        self.axis2 = self.data[0].axis2
#        self.col0 = self.data[0].col
#        self.row0 = self.data[0].row
        self.highmass = self.data[0].highmass     # a tuple (h1, h2)
        self.Debug = Debug
        if report: self.report()

    def load(self):
        "load from file"
        self.data = []
        for i in range(8):
            try:
                dl = FTICR.FTICRData(name=self.name, mode="onfile", group="resol%d"%(i+1))
            except tables.NoSuchNodeError:
                pass
            else:
                dl.unit='m/z'
                self.data.append(dl)
    def report(self):
        "report object content"
        print (self.name)
        print(self.data[0])
        print('=====================')
        print('multiresolution data:\n#: Size')
        for i,dl in enumerate(self.data):
            print ("%d: %d x %d : %.0f Mpix"%( i+1, dl.size1, dl.size2, dl.size1*dl.size2/1024/1024 ))
    def _report(self):
        "low level report"
        for i,dl in enumerate(self.data):
            print(i+1)
            print(dl)
    def colmz(self,i):
        "return a column with coordinate in m/z"
        return self.col(self.axis2.mztoi(i))
    def rowmz(self,i):
        "return a row with coordinate in m/z"
        return self.row(self.axis1.mztoi(i))
    def col(self,i):
        "return a column with coordinate in index"
        return self.col(int(round(i)))
    def row(self,i):
        "return a row with coordinate in index"
        return self.row(int(round(i)))
    def to_display(self,zoom=((0,FTICR.FTMS.HighestMass),(0,FTICR.FTMS.HighestMass)), verbose=False):
        """
        computes and return which dataset to display at a given zoom and scale level"
        in: zoom = ((F1low, F1up), (F2low,F2up))  - in m/z

        out: a tuple (data, zoomwindow), where data is a NPKData and zoomwindow an eventually recalibrated zoom window

        so, if DATA is an MR() object and Zoom a defined window in m/z ((F1low, F1up), (F2low,F2up))
        the sequence:
            datasel, zz = DATA.to_display(Zoom)
            datasel.display(zoom=zz, scale=...)
        will display the selected zone with the best possible resolution
        """
        z1,z2,z3,z4 = flatten(zoom)
#        print ("m/z F1: %.1f - %.1f  / F2 %.1f - %.1f"%(z1,z2,z3,z4))
        z1,z2 = min(z1,z2), max(z1,z2)
        z3,z4 = min(z3,z4), max(z3,z4)
        for reso,dl in enumerate(self.data):
            z11 = max(z1, dl.axis1.lowmass)
            z33 = max(z3, dl.axis2.lowmass)
            z22 = min(z2, dl.highmass[0])
            z44 = min(z4, dl.highmass[1])
            z1lo, z1up, z2lo, z2up = parsezoom(dl,(z11,z22,z33,z44))
            sz = (z1lo-z1up)* (z2lo-z2up)
            if sz < self.SIZEMAX:
                #print (dl.size1,dl.size2, z1lo, z1up, z2lo, z2up, sz)
                break
        zooml = (dl.axis1.itomz(z1lo), dl.axis1.itomz(z1up), dl.axis2.itomz(z2lo), dl.axis2.itomz(z2up))
        if verbose or self.Debug:
            print ("zoom: F1 %.1f-%.1f  /  F2 %.1f-%.1f"%zooml)
            print ("resolution level %d - %.1f Mpix zoom window"%(reso, sz/1024/1024))
        return dl, zooml
    def compute_absmax(self):
        "computes largest point from smaller resolution, and propagates"
        dsmall = self.data[-1]
        self.absmax = dsmall.absmax
        for dl in self.data:
            dl._absmax = self.absmax        

class MR_interact(MR):
    def __init__(self, name, figsize=None, report=True, show=True, Debug=False):
        """
        creates an interactive object.
        if display is True (default) the graphical tool will be displayed.
        """
        super(self.__class__, self).__init__(name, report=report, Debug=Debug)
        self.vlayout = Layout(width='60px')
        self.pltaxe = None
        self.pltaxe1D = None
        if figsize is None:
            self.figsize = (10,8)
        else:
            self.figsize = (figsize[0]/2.54, figsize[1]/2.54, )
#        self.reset_track()
        self.check_fig()
        if show:  self.show()

    def bb(self, name, desc, action, layout=None, tooltip=""):
        "build a button into self"
        if layout is None: layout = self.vlayout
        butt = widgets.Button(description=desc, layout=layout, tooltip=tooltip)
        butt.on_click(action)
        setattr(self, name, butt)
##################### 2D #######################
    def zoom_box(self):
        "defines the zoom box widget"
        # .  F1 .
        # F2 .  F2
        # .  F1 .
        self.scale = widgets.FloatLogSlider(description='scale:', value=1.0, min=-1, max=3,  base=10, step=0.01,
            layout=Layout(width='80%'), continuous_update=HEAVY)

        wf = widgets.BoundedFloatText
        ref = self.data[0]
        style = {'description_width': 'initial'}
        lay = Layout(width='100px', height='30px')
        self.z1l = wf(  min=ref.axis1.lowmass, max=ref.highmass[0],
            tooltip='vertical zoom', style=style, layout=lay)
        self.z1h = wf(  min=ref.axis1.lowmass, max=ref.highmass[0],
            style=style, layout=lay)
        self.z2l = wf(  min=ref.axis2.lowmass, max=ref.highmass[1],
            style=style, layout=lay)
        self.z2h = wf(  min=ref.axis2.lowmass, max=ref.highmass[1],
            style=style, layout=lay)
        self.fullzoom()
        self.bb('b_zupdate', 'Apply', lambda e : self.display(),
            layout=Layout(width='100px'), tooltip="Set zoom to values")
        self.bb('b_reset', 'Reset', self.reset,
            layout=Layout(width='100px'), tooltip="Set zoom to values")
        blank = widgets.HTML("&nbsp;",layout=lay)
        label =  widgets.HTML('Zoom Window (in <i>m/z</i>)')
        innerbox = VBox([ label,
                    HBox([blank,    self.z1h, blank]),  
                    HBox([self.z2l, self.b_zupdate,    self.z2h]),
                    HBox([blank,    self.z1l, blank]),  
                ]
            )
        action_box = VBox([blank ])
        box = HBox([ self.b_reset, innerbox])
        return innerbox
    def spec_box(self):
        "defines the spectral box widget"
        self.scale = widgets.FloatSlider(description='scale:', value=1.0, min=1.0, max=20, step=0.1,
                    tooltip='Set display scale', layout=Layout(width='100px', height='400px'), continuous_update=HEAVY,
                    orientation='vertical')
        self.scale.observe(self.ob)
        self.bb('b_redraw', 'Redraw', lambda e : self.display(),
            layout=Layout(width='100px'))
        box = VBox([self.b_reset, self.scale, self.b_redraw])
        return HBox([box, self.pltfig.canvas])
    def scale_up(self, step):
        self.scale.value *= 1.1892**step # 1.1892 is 4th root of 2.0
    def show(self):
        "actually show the graphical tool and the interactive spectrum"
        display(self.box)
        display(self.sbox)
        self.display()
    def ob(self, event):
        "observe events and display"
        if event['name'] != 'value':
            return
        self.display()
    def display(self, zoom=None, scale=None):
        "computes pictures (display in the SPIKE sense - not ipywidget one)"
        if zoom is not None:
            _zoom = zoom
        else:
            _zoom = (self.z1l.value, self.z1h.value, self.z2l.value, self.z2h.value)
        if scale is not None:
            self.scale.value = scale
        datasel, zz = self.to_display(_zoom)
        corner = (zz[1],zz[3])
        reso = [corner[i]/datasel.axes(i+1).deltamz(corner[i]) for i in range(2)]
#        print(corner, reso)
        self.check_fig()  # insure figure is there
#        with self.out:
        self.pltaxe.clear() # clear 
        datasel.display(zoom=zz, scale=self.scale.value, absmax=self.absmax, 
            xlabel='F2    m/z', ylabel='F1    m/z',
            show=False, figure=self.pltaxe)
        self.pltaxe.text(corner[1],corner[0],"D#%d R: %.0fx%.0f"%(self.data.index(datasel)+1, *reso))
    def check_fig(self):
        "create figure if missing"
        if self.pltaxe is None:
            plt.ioff()
#            plt.clf()
            fg,ax = plt.subplots(figsize=self.figsize)
            self.pltaxe = ax
            self.pltfig = fg
            self.set_on_redraw()
            self.box = self.zoom_box()
            self.sbox = self.spec_box()
            plt.ion()

    def update(self, e):
        "update internal zoom coordinates"
#        self.track.append((self._zoom, self.scale))
#        self.point = -1 # means last
        #print('update')
        xb = self.pltaxe.get_xbound()
        yb = self.pltaxe.get_ybound()
        #print( xb, yb)
        self.z1l.value = yb[0]
        self.z1h.value = yb[1]
        self.z2l.value = xb[0]
        self.z2h.value = xb[1]
    def set_on_redraw(self):
        def on_press(event):
            print('you pressed', event.button, event.xdata, event.ydata)
        def on_scroll(event):
            self.scale_up(event.step)
        cidd = self.pltfig.canvas.mpl_connect('draw_event', self.update)
        cids = self.pltfig.canvas.mpl_connect('scroll_event', on_scroll)
#        cidc = self.pltfig.canvas.mpl_connect('button_press_event', on_press)
    def reset(self, b):
        self.scale.value = 1.0
        self.fullzoom()
        self.reset_track()
        self.display()
    def fullzoom(self):
        self.z1l.value = self.axis2.lowmass+0.01
        self.z1h.value = self.highmass[0]
        self.z2l.value = self.axis2.lowmass+0.01
        self.z2h.value = self.highmass[1]

# ####### zoom track - experimental ########
    def reset_track(self):
        self.track = []
        self.point = -1
#     def back(self, *arg):
#         self.point -= 1
#         try:
#             self._zoom, self.scale = self.track[self.point]
#         except IndexError:
#             self.point = -len(self.track)
#             self._zoom, self.scale = self.track[0]
#         self.update()
#         self.display(redraw=True)
#     def forw(self, *arg):
#         self.point += 1
#         try:
#             self._zoom, self.scale = self.track[self.point]
#         except IndexError:
#             self._zoom, self.scale = self.track[-1]
#             self.point = -1
#         self.update()
#         self.display(redraw=True)
# ##################### 1D ##################
    def I1D(self):
        "show the 1D selector"
        self.r1D = None
        self.t1D = ''
        self.i1D = None
        self.check_fig1D()
        display(self.ext_box())

    def check_fig1D(self):
        if self.pltaxe1D is None:
            fg,ax = plt.subplots(figsize=(1.5*self.figsize[0], 0.75*self.figsize[0]))
            ax.text(0.1, 0.8, 'Empty - use "horiz" and "vert" buttons above')
            self.pltaxe1D = ax
            self.fig1D = fg

    def ext_box(self):
        "defines the interactive tools for 1D"
        wf = widgets.BoundedFloatText
        wi = widgets.BoundedIntText
        ref = self.data[0]
        style = {'description_width': 'initial'}
        lay = Layout(width='120px', height='30px')
        lay2 = Layout(width='50px', height='30px')
        self.z1 = wf( value=300.0, min=ref.axis1.lowmass, max=ref.highmass[0], description='F1', style=style, layout=lay)
        self.z2 = wf( value=300.0,  min=ref.axis2.lowmass, max=ref.highmass[1], description='F2', style=style, layout=lay)
        self.horiz = wi( value=1, min=1, max=20, style=style, layout=lay2)
        self.vert = wi( value=1, min=1, max=20, style=style, description="/", layout=lay2)
        def lrow(b, inc=0):
            rint = int(round(ref.axis1.mztoi(self.z1.value)))
            if inc !=0:
                rint += inc
            self.z1.value = ref.axis1.itomz(rint)
            if self.t1D == 'col':
                self.pltaxe1D.clear()
                self.r1D = None
            if self.b_accu.value == 'sum' and self.r1D is not None:
                self.r1D += ref.row(rint)
            else:
                self.r1D = ref.row(rint)
            self.t1D = 'row'
            self.i1D = rint
            self.display1D()
        def lrowp1(b):
            lrow(b,-1)
        def lrowm1(b):
            lrow(b,1)
        def lcol(b, inc=0):
            rint = int(round(ref.axis2.mztoi(self.z2.value)))
            if inc !=0:
                rint += inc
            self.z2.value = ref.axis2.itomz(rint)
            if self.t1D == 'row':
                self.pltaxe1D.clear()
                self.r1D = None
            if self.b_accu.value == 'sum' and self.r1D is not None:
                self.r1D += ref.col(rint)
            else:
                self.r1D = ref.col(rint)
            self.t1D = 'col'
            self.i1D = rint
            self.display1D()
        def lcolp1(b):
            lcol(b,-1)
        def lcolm1(b):
            lcol(b,1)
        def on_press(event):
            v = self.r1D.axis1.mztoi(event.xdata)
            if self.t1D == 'col':
                self.z1.value = v
            elif self.t1D == 'row':
                self.z2.value = v
        cids = self.fig1D.canvas.mpl_connect('button_press_event', on_press)

        self.bb('b_row', 'horiz', lrow, layout=Layout(width='60px'), tooltip='extract an horizontal row')
        self.bb('b_rowp1', '+1', lrowp1, layout=Layout(width='30px'), tooltip='next row up')
        self.bb('b_rowm1', '-1', lrowm1, layout=Layout(width='30px'), tooltip='next row down')
        self.bb('b_col', 'vert', lcol, layout=Layout(width='60px'), tooltip='extract a vertical col')
        self.bb('b_colp1', '+1', lcolp1, layout=Layout(width='30px'), tooltip='next col right')
        self.bb('b_colm1', '-1', lcolm1, layout=Layout(width='30px'), tooltip='next col left')
        self.b_accu = widgets.Dropdown(options=['off', 'graphic', 'sum'],
                value='off', description='Accumulate plots while scanning:', style=style)
        return VBox([ widgets.HTML('<h3>Extract 1D MS Spectrum going through given F1-F2 coordinates</h3>'),
                    HBox([widgets.HTML("<B>coord:</B>"),self.z1,
                          self.b_row, self.b_rowp1, self.b_rowm1, self.b_accu]),
                    HBox([widgets.HTML("<B>coord:</B>"),self.z2,
                          self.b_col, self.b_colp1, self.b_colm1])])
    def display1D(self):
        "display the selected 1D"
        if self.t1D == 'row':
            title = 'horizontal extract at F1=%f m/z (index %d)'%(self.z1.value, self.i1D)
            label = str(self.z1.value)
        elif self.t1D == 'col':
            title = 'vertical extract at F2=%f m/z (index %d)'%(self.z2.value, self.i1D)
            label = str(self.z2.value)
        if self.b_accu.value != 'graphic':
            self.pltaxe1D.clear()
            label = None
        self.r1D.display(xlabel='m/z', show=False, figure=self.pltaxe1D, new_fig=False, label=label, title=title)


class MSPeaker(object):
    "a peak-picker for MS experiments"
    def __init__(self, npkd, pkname):
        if not isinstance(npkd, FTMSData):
            raise Exception('This modules requires a FTMS Dataset')
        self.npkd = npkd
        self.pkname = pkname
        self.zoom = widgets.FloatRangeSlider(value=[npkd.axis1.lowmass, npkd.axis1.highmass],
            min=npkd.axis1.lowmass, max=npkd.axis1.highmass, step=0.1,
            layout=Layout(width='100%'), description='zoom',
            continuous_update=False, readout=True, readout_format='.1f',)
        self.zoom.observe(self.display)
        self.tlabel = Label('threshold (x noise level):')
        self.thresh = widgets.FloatLogSlider(value=20.0,
            min=np.log10(1), max=2.0, base=10, step=0.01, layout=Layout(width='30%'),
            continuous_update=False, readout=True, readout_format='.1f')
        self.thresh.observe(self.pickpeak)
        self.peak_mode = widgets.Dropdown(options=['marker', 'bar'],value='marker',description='show as')
        self.peak_mode.observe(self.display)
        self.bexport = widgets.Button(description="Export",layout=Layout(width='7%'),
                button_style='success', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='Export to csv file')
        self.bexport.on_click(self.pkexport)
        self.bprint = widgets.Button(description="Print", layout=Layout(width='7%'),
                button_style='success', tooltip='Print to screen')
        self.bprint.on_click(self.pkprint)
        self.bdone = widgets.Button(description="Done", layout=Layout(width='7%'),
                button_style='warning', tooltip='Fix results')
        self.bdone.on_click(self.done)
        #self.spec = Output(layout={'border': '1px solid black'})
        self.out = Output(layout={'border': '1px solid red'})
        display( VBox([self.zoom,
                      HBox([self.tlabel, self.thresh, self.peak_mode, self.bprint, self.bexport, self.bdone])
                      ]) )
        self.fig, self.ax = plt.subplots()
        self.npkd.set_unit('m/z').peakpick(autothresh=self.thresh.value, verbose=False, zoom=self.zoom.value).centroid()
        self.display()
        display(self.out)

    def pkprint(self,event):
        self.out.clear_output(wait=True)
        with self.out:
            display(HTML(self.npkd.pk2pandas().to_html()))
            #print(self.pklist())
    def pkexport(self,event):
        "exports the peaklist to file"
        with open(self.pkname,'w') as FPK:
            print(self.pklist(),file=FPK)
        print('Peak list stored in ',self.pkname)
    def pklist(self):
        "creates peaklist"
        text = ["m/z\t\tInt.(%)\tR\tarea(a.u.)"]
        data = self.npkd
        intmax = max(data.peaks.intens)/100
        for pk in data.peaks:
            mz = data.axis1.itomz(pk.pos)
            Dm = 0.5*(data.axis1.itomz(pk.pos-pk.width) - data.axis1.itomz(pk.pos+pk.width))
            area = pk.intens*Dm
            l = "%.6f\t%.1f\t%.0f\t%.0f"%(mz, pk.intens/intmax, round(mz/Dm,-3), area)
            text.append(l)
        return "\n".join(text)
    def display(self, event={'name':'value'}):
        "display spectrum and peaks"
        if event['name']=='value':  # event is passed by GUI - make it optionnal
            self.ax.clear()
            self.npkd.display(new_fig=False, figure=self.ax, zoom=self.zoom.value)
            try:
                self.npkd.display_peaks(peak_label=True, peak_mode=self.peak_mode.value, figure=self.ax, zoom=self.zoom.value, NbMaxPeaks=NbMaxDisplayPeaks)
                x = self.zoom.value
                y = [self.npkd.peaks.threshold]*2
                self.ax.plot(x,y,':r')
                self.ax.annotate('%d peaks detected'%len(self.npkd.peaks) ,(0.05,0.95), xycoords='figure fraction')
            except:
                pass
    def pickpeak(self, event):
        "interactive wrapper to peakpick"
        if event['name']=='value':
            self.pp()
    def pp(self):
        "do the peak-picking calling pp().centroid()"
        #self.spec.clear_output(wait=True)
        self.npkd.set_unit('m/z').peakpick(autothresh=self.thresh.value, verbose=False, zoom=self.zoom.value).centroid()
        self.display()
    def done(self, event):
        "exit GUI"
        for w in [self.zoom, self.thresh, self.peak_mode, self.bprint, self.bexport, self.bdone]:
            w.close()
        self.tlabel.value = "threshold %.2f noise level"%self.thresh.value
#        self.display()

class Calib(object):
    "a simple tool to show and modify calibration cste"
    def __init__(self, data):
        self.data = data
        self.res = [data.axis1.calibA, data.axis1.calibB, data.axis1.calibC]
        self.A = widgets.FloatText(value=data.axis1.calibA, description="A")
        self.B = widgets.FloatText(value=data.axis1.calibB, description="B")
        self.C = widgets.FloatText(value=data.axis1.calibC, description="C")
        self.bupdate = widgets.Button(description="Update",
                button_style='success', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='set current data-sets to displayed values')
        self.bupdate.on_click(self.update)
        self.bback = widgets.Button(description="Restore",
                button_style='success', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='restore dataset to initial values')
        self.bback.on_click(self.back)
        display(VBox([  HBox([self.A, widgets.Label('Hz/Th')]),
                        HBox([self.B, widgets.Label('Hz')]),
                        HBox([self.C, widgets.Label('Hz/Th^2')]),
                        HBox([self.bupdate, self.bback])]))
    def update(self,event):
        self.data.axis1.calibA = self.A.value
        self.data.axis1.calibB = self.B.value
        self.data.axis1.calibC = self.C.value
        print("Set - don't forget to rerun the Peak Picker")
    def back(self,event):
        self.data.axis1.calibA = self.res[0]
        self.data.axis1.calibB = self.res[1]
        self.data.axis1.calibC = self.res[2]
        self.A.value = self.res[0]
        self.B.value = self.res[1]
        self.C.value = self.res[2]

Colors = ('black','red','blue','green','orange',
'blueviolet','crimson','turquoise','indigo',
'magenta','gold','pink','purple','salmon','darkblue','sienna')

class SpforSuper(object):
    "a holder for SuperImpose"
    def __init__(self, i, name):
        j = i%len(Colors)
        self.name = widgets.Text(value=name, layout=Layout(width='70%'))
        self.color = widgets.Dropdown(options=Colors,value=Colors[j],layout=Layout(width='10%'))
        self.direct = widgets.Dropdown(options=['up','down','off'],value='off', layout=Layout(width='10%'))
        self.me = HBox([widgets.HTML(value="<b>%d</b>"%i),self.name, self.color,self.direct])
        self.fig = False
    def display(self):
        if self.name != 'None' and self.direct.value != 'off':
            scale = 1
            if self.direct.value == 'up':
                mult = 1
            elif self.direct.value == 'down':
                mult = -1
            else:
                return
            FTICRData(name=self.name.value).set_unit('m/z').mult(mult).display(
                new_fig=self.fig,
                scale=scale,
                color=self.color.value,
                label=op.basename(op.dirname(self.name.value)))

class SuperImpose(object):
    "a tool to superimpose spectra"
    def __init__(self, base=None, filetype='*.msh5', N=None):
        if N is None:
            N = int(input('how many spectra do you want to compare:  '))
        self.Chooser = FileChooser(base=base, filetype=filetype, mode='r', show=False)
        self.bsel = widgets.Button(description='Copy',layout=Layout(width='10%'),
                button_style='info', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='copy selected data-set to entry below')
        self.to = widgets.IntText(value=1,min=1,max=N,layout=Layout(width='10%'))
        self.bsel.on_click(self.copy)
        self.bdisplay = widgets.Button(description='Display',layout=Layout(width='10%'),
                button_style='info', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='display superimposition')
        self.bdisplay.on_click(self.display)
        self.spec = Output(layout={'border': '1px solid black'})
        self.DataList = [SpforSuper(i+1,'None') for i in range(N)]
        self.DataList[0].color.value = 'black'
        self.DataList[0].fig = True  # switches on the very first one
    def Show(self):
        display(widgets.Label('Select a file, and click on the Copy button to copy it to the chosen slot'))
        self.Chooser.show()
        display(HBox([self.bsel, widgets.Label('to'), self.to]))
        display(VBox([sp.me for sp in self.DataList]))
        display(self.bdisplay)
        display(self.spec)
    def copy(self, event):
        if self.to.value <1 or self.to.value >len(self.DataList):
            print('Destination is out of range !')
        else:
            self.DataList[self.to.value-1].name.value = self.Chooser.file
            self.DataList[self.to.value-1].direct.value = 'up'
        self.to.value = min(self.to.value, len(self.DataList)) +1
    def display(self, event):
        self.spec.clear_output(wait=True)
        for i,s in enumerate(self.DataList):
            with self.spec:
                s.display()

class MS2Dscene(HBox):
    "a widget to set all MS tools into one screen"
    def __init__(self, show=True, root='/'):
        super(self.__class__, self).__init__()
        # DATA
        self.MR2D = None # The MR2D object

        # GUI tools
        self.FC = FileChooser(root, filetype='*.msh5', mode='r')

        self.bdisp2D = Button(description='Load 2D',layout=Layout(width='20%'),
                button_style='success', tooltip='load and display interactive 2D map')
        self.bdisp2D.on_click(self.load2D)

        # GUI
        self.header = Output()
        with self.header:
            display(Markdown('---\n# Select a file, and choose a tool'))
            display(self.FC)
            display(self.bdisp2D)

        self.out2D = Output()  # the area where 2D is shown
        with self.out2D:
            display(HTML("<br><br><h3><i><center>No Data</center></i></h3>"))

        # self.outpp2D = Output()  # the area where peak picking is shown
        # with self.outpp2D:
        #     display(HTML("<br><br><h3><i><center>No Data</center></i></h3>"))

        self.out1D = Output()  # the area where 1D is shown
        with self.out1D:
            display(HTML("<br><br><h3><i><center>No Data</center></i></h3>"))

        self.outinfo = Output()  # the area where info is shown
        with self.outinfo:
            display(HTML("<br><br><h3><i><center>No Data</center></i></h3>"))

        # GUI set-up and show
        self.tabs = widgets.Tab()
        self.tabs.children = [ self.out2D, self.out1D, self.outinfo ]
        self.tabs.set_title(0, '2D Display')
        self.tabs.set_title(1, '1D extraction')
#        self.tabs.set_title(2, 'Peak Picking')
        self.tabs.set_title(2, 'info')

#        self.tabs = VBox([ self.out2D, self.outpp2D, self.out1D, self.outinfo ])
        self.box = VBox([   self.header,
                            self.tabs
                        ])
        # self.box = VBox([self.title,
        #                 self.FC,
        #                 HBox([self.bdisp2D, self.bpp2D, self.bdisp1D])
        #                 ])
        if show:
            display(self.box)

    def load2D(self, e):
        "create 2D object and display"
        self.out2D.clear_output(wait=True)
#        self.outpp2D.clear_output(wait=True)
        self.out1D.clear_output(wait=True)
        self.outinfo.clear_output(wait=True)
        try:
            self.MR2D = MR_interact(self.FC.selected, 
                report=False, show=False, figsize=(15,15), Debug=False)
        except FileNotFoundError:
            self.MR2D = None
            with self.out2D:
                display(HTML("<br><br><h3><i><center>No Data or wrong Data</center></i></h3>"))
        if self.MR2D is not None:     
            with self.out2D:
                self.MR2D.show()
#                display(self.MR2D.box)
#                display(self.MR2D.sbox)
            with self.out1D:
                self.MR2D.I1D()
            # with self.outpp2D:
            #     display(self.MR2D.box)  # copie of the main pane
            #     display(self.MR2D.sbox)
            with self.outinfo:
                self.MR2D.report()
                print()
