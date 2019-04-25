#!/usr/bin/env python 
# encoding: utf-8

"""
A tool to display 2D FT-ICR data-sets

to be embedded in jupyter notebook

First version MAD jan 2019
preliminary and not fully tested !

"""
import os.path as op
import tables
import matplotlib.pyplot as plt
from ipywidgets import interact, fixed, HBox, VBox, HTML, Label, Layout, Output
import ipywidgets as widgets
from IPython.display import display

from .. import FTICR
from ..NPKData import flatten, parsezoom
from .INTER import FileChooser
from ..FTMS import FTMSData
from ..FTICR import FTICRData


SIZEMAX = 8*1024*1024       # largest zone to display

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
        #self.check_fig()
        self._zoom = None
        self.fullzoom()   # in m/z
        self.scale = 1.0
        self.reset_track()
        self.bb('b_scup', "×", self.scaleup, tooltip="increase display scale")
        self.bb('b_scdown', "÷", self.scaledown, tooltip="decrease display scale")
        self.bb('b_screset', '≡', self.scale1, tooltip="reset display scale")

        self.bb('b_left', "◀︎", self.left, tooltip="move zoom window to the left")
        self.bb('b_right', "▶︎", self.right, tooltip="move zoom window to the right")
        self.bb('b_up', "▲", self.up, tooltip="move zoom window up")
        self.bb('b_down', "▼", self.down, tooltip="move zoom window down")
        self.bb('b_zin', "🔍 in", self.zoom_in, tooltip="zoom in")
        self.bb('b_zout', "🔍out", self.zoom_out, tooltip="zoom out")

        self.bb('b_reset', "⌘", self.reset, tooltip="reset display to default")
        self.bb('b_back', '⍇', self.back, tooltip="back in zoom list")
        self.bb('b_forw', '⍈', self.forw, tooltip="forward in zoom list")

        if show:  self.show()

    def bb(self, name, desc, action, layout=None, tooltip=""):
        "build button"
        if layout is None: layout = self.vlayout
        butt = widgets.Button(description=desc, layout=layout, tooltip=tooltip)
        butt.on_click(action)
        setattr(self, name, butt)
##################### 1D ##################
    def I1D(self):
        self.r1D = None
        self.t1D = ''
        self.i1D = None
        display(self.ext_box())
        self.check_fig1D()

    def check_fig1D(self):
        if self.pltaxe1D is None:
            fg,ax = plt.subplots(figsize=(self.figsize[0],self.figsize[0]/2))
            ax.text(0.1, 0.8, 'Empty - use "horiz" and "vert" buttons above')
            self.pltaxe1D = ax

    def ext_box(self):
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
        self.bb('b_row', 'horiz', lrow, layout=Layout(width='60px'), tooltip='extract an horizontal row')
        self.bb('b_rowp1', '+1', lrowp1, layout=Layout(width='30px'), tooltip='next row up')
        self.bb('b_rowm1', '-1', lrowm1, layout=Layout(width='30px'), tooltip='next row down')
        self.bb('b_col', 'vert', lcol, layout=Layout(width='60px'), tooltip='extract a vertical col')
        self.bb('b_colp1', '+1', lcolp1, layout=Layout(width='30px'), tooltip='next col right')
        self.bb('b_colm1', '-1', lcolm1, layout=Layout(width='30px'), tooltip='next col left')
        self.b_accu = widgets.Dropdown(options=['off', 'graphic', 'sum'],
                value='off', description='Accumulate plots while scanning:', style=style)
        return VBox([ HTML('Extract 1D MS Spectrum passing by F1-F2 coordinates'),
                    HBox([HTML("<B>coord:</B>"),self.z1, self.z2,
                          self.b_row, self.b_rowp1, self.b_rowm1,HTML("&nbsp;/&nbsp;"),
                          self.b_col, self.b_colp1, self.b_colm1, self.b_accu])])
    def display1D(self):
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
##################### 2D #######################
    def zoom_box(self):
        wf = widgets.BoundedFloatText
            # z11 = max(z1, dl.axis1.lowmass)
            # z33 = max(z3, dl.axis2.lowmass)
            # z22 = min(z2, dl.highmass[0])
            # z44 = min(z4, dl.highmass[1])
        ref = self.data[0]
        style = {'description_width': 'initial'}
        lay = Layout(width='80px', height='30px')
        self.z1l = wf( value=self._zoom[0], min=ref.axis1.lowmass, max=ref.highmass[0],
            description='F1:', style=style, layout=lay)
        self.z1h = wf( value=self._zoom[1], min=ref.axis1.lowmass, max=ref.highmass[0],
            description='..', style=style, layout=lay)
        self.z2l = wf( value=self._zoom[2], min=ref.axis2.lowmass, max=ref.highmass[1],
            description='F2:', style=style, layout=lay)
        self.z2h = wf( value=self._zoom[3], min=ref.axis2.lowmass, max=ref.highmass[1],
            description='..', style=style, layout=lay)
        def zupdate(b):
            self._zoom = (self.z1l.value, self.z1h.value, self.z2l.value, self.z2h.value)
            self.update()
            self.display()
        self.bb('b_zupdate', 'Update', zupdate, layout=Layout(width='300px'), tooltip="Set zoom to values")
        return VBox([ HTML('Zoom Window (in $m/z$)'),
                    HBox([self.z1l, self.z1h, self.z2l, self.z2h]),
                    self.b_zupdate],
                    layout=Layout(  ))
    def show(self):
        "actually show the graphical tool and the interactive spectrum"
        # .^..+
        # <=>.=
        # .v..-
        blank = widgets.Button(description=" ", layout=self.vlayout)
        zmbox1 = VBox([self.b_back,  self.b_left, self.b_zin])
        zmbox2 = VBox([self.b_up, self.b_reset, self.b_down])
        zmbox3 = VBox([self.b_forw,  self.b_right, self.b_zout])
        zmbox4 = VBox([self.b_scup, self.b_screset, self.b_scdown])
        spacer = HTML('&nbsp;&nbsp;&nbsp;&nbsp;')
        zbox = self.zoom_box()
        display(HBox([zmbox1,zmbox2, zmbox3, spacer, zmbox4, spacer, zbox],
            layout=Layout(height='100px')))
#        zmbox1 = HBox([blank,       self.b_up, blank, blank, self.b_scup])
#        zmbox2 = HBox([self.b_left, self.b_full, self.b_right, blank, self.b_sc1])
#        zmbox3 = HBox([blank,  self.b_down, blank,  blank, self.b_down])

#        display(VBox([zmbox1, zmbox2, zmbox3]))
        self.update()
        # fg,ax = plt.subplots()
        # self.pltaxe = ax
        self.display()
    def display(self, zoom=None, scale=None, redraw=None):
        if zoom is not None:
            self._zoom = zoom
        if scale is not None:
            self.scale = scale
        datasel, zz = self.to_display(self._zoom)
        corner = (zz[1],zz[3])
        reso = [corner[i]/datasel.axes(i+1).deltamz(corner[i]) for i in range(2)]
#        print(corner, reso)
        nf = not redraw
        self.check_fig()
        self.pltaxe.clear()
        # datasel.bokeh(zoom=zz, scale=self.scale, redraw=redraw, show=False)        datasel.display(zoom=zz, scale=self.scale, show=False, figure=self.pltaxe, new_fig=nf)
        datasel.display(zoom=zz, scale=self.scale, absmax=self.absmax, 
            xlabel='F2    m/z', ylabel='F1    m/z',
            show=False, figure=self.pltaxe, new_fig=nf)
        self.pltaxe.text(corner[1],corner[0],"D#%d R: %.0fx%.0f"%(self.data.index(datasel)+1, *reso))
    def reset_track(self):
        self.track = []
        self.point = -1
    def check_fig(self):
        if self.pltaxe is None:
            fg,ax = plt.subplots(figsize=self.figsize)
            self.pltaxe = ax
    def reset(self, b):
        self.scale = 1.0
        self.fullzoom()
        self.reset_track()
        self.update()
        self.display(redraw=True)
    def fullzoom(self):
        self._zoom=(self.axis1.lowmass+0.01, self.highmass[0],self.axis2.lowmass+0.01   , self.highmass[1])
    def update(self):
        self.track.append((self._zoom, self.scale))
        self.point = -1 # means last
        (self.z1l.value, self.z1h.value, self.z2l.value, self.z2h.value) = self._zoom
        if self.Debug: print(self.zoom, self.scale)
    @property
    def zoom(self):
        return "%.3f %.3f %.3f %.3f"%self._zoom
    def zoomwidth(self):
        return ( self._zoom[1]-self._zoom[0], self._zoom[3]-self._zoom[2])
    def zoom_in(self, b, factor=1.44):
        """
        z_in: waft = wbef/factor => waft*factor = wbef
                dw = wbef/factor-wbef = wbef*(1-factor)
                dw = waft*factor*(1-factor)
        """
        w1,w2 = self.zoomwidth()
        deltaw1 = w1*(factor-1)/2
        deltaw2 = w2*(factor-1)/2
        z1,z2,z3,z4 = self._zoom
        self._zoom = (z1+deltaw1, z2-deltaw1, z3+deltaw2, z4-deltaw2)
        self.update()
        self.display()
    def zoom_out(self, b, factor=1.44):
        """
        z_in: dw = waft*factor*(1-factor)
        z_out:dw = wbef*(factor*(1-factor))  
        """
        w1,w2 = self.zoomwidth()
        deltaw1 = w1*factor*(factor-1)/2
        deltaw2 = w2*factor*(factor-1)/2
        z1,z2,z3,z4 = self._zoom
        zoom = (z1-deltaw1, z2+deltaw1, z3-deltaw2, z4+deltaw2)
        self._zoom = (max(zoom[0],self.axis1.lowmass), min(zoom[1],self.highmass[0]),
                      max(zoom[2],self.axis2.lowmass), min(zoom[3],self.highmass[1]) )
        self.update()
        self.display()
    def up(self, b, factor=1.44):
        w1,w2 = self.zoomwidth()
        deltaw1 = w1*factor*(factor-1)/2
        z1,z2,z3,z4 = self._zoom
        zoom = [z1+deltaw1, z2+deltaw1]
        if zoom[1]>self.highmass[0]:
            zoom[1] = self.highmass[0]
            zoom[0] = zoom[1] - w1
        self._zoom = (zoom[0], zoom[1], z3, z4)
        self.update()
        self.display()
    def down(self, b, factor=1.44):
        w1,w2 = self.zoomwidth()
        deltaw1 = w1*factor*(factor-1)/2
        z1,z2,z3,z4 = self._zoom
        zoom = [z1-deltaw1, z2-deltaw1]
        if zoom[0]<self.axis1.lowmass:
            zoom[0] = self.axis1.lowmass
            zoom[1] = zoom[0] + w1
        self._zoom = (zoom[0], zoom[1], z3, z4)
        self.update()
        self.display()
    def left(self, b, factor=1.44):
        w1,w2 = self.zoomwidth()
        deltaw2 = w2*factor*(factor-1)/2
        z1,z2,z3,z4 = self._zoom
        zoom = [z3-deltaw2, z4-deltaw2]
        if zoom[0]<self.axis2.lowmass:
            zoom[0] = self.axis2.lowmass
            zoom[1] = zoom[0] + w2
        self._zoom = (z1, z2, zoom[0], zoom[1])
        self.update()
        self.display()
    def right(self, b, factor=1.44):
        w1,w2 = self.zoomwidth()
        deltaw2 = w2*factor*(factor-1)/2
        z1,z2,z3,z4 = self._zoom
        zoom = [z3+deltaw2, z4+deltaw2]
        if zoom[1]>self.highmass[1]:
            zoom[1] = self.highmass[1]
            zoom[0] = zoom[1] - w2
        self._zoom = (z1, z2, zoom[0], zoom[1])
        self.update()
        self.display()
    def scaleup(self, b):
        self.scale *= 1.4
        self.update()
        self.display(redraw=True)
    def scaledown(self, b):
        self.scale /= 1.4
        self.update()
        self.display(redraw=True)
    def scale1(self, b):
        self.scale = 1.0
        self.update()
        self.display(redraw=True)
    def back(self, *arg):
        self.point -= 1
        try:
            self._zoom, self.scale = self.track[self.point]
        except IndexError:
            self.point = -len(self.track)
            self._zoom, self.scale = self.track[0]
        self.update()
        self.display(redraw=True)
    def forw(self, *arg):
        self.point += 1
        try:
            self._zoom, self.scale = self.track[self.point]
        except IndexError:
            self._zoom, self.scale = self.track[-1]
            self.point = -1
        self.update()
        self.display(redraw=True)

class MR2(MR_interact):  # In development - do not use
    def back(self, *arg):
        self.point -= 1
        try:
            self._zoom, self.scale = self.track[self.point]
        except IndexError:
            self.point = -len(self.track)
            self._zoom, self.scale = self.track[0]
    def forw(self, *arg):
        self.point += 1
        try:
            self._zoom, self.scale = self.track[self.point]
        except IndexError:
            self._zoom, self.scale = self.track[-1]
            self.point = -1
    def _display(self, zoom=None, scale=None, redraw=None):
        if zoom is not None:
            self._zoom = zoom
        if scale is not None:
            self.scale = scale
        datasel, zz = self.to_display(self._zoom)
        corner = (zz[1],zz[3])
        reso = [corner[i]/datasel.axes(i+1).deltamz(corner[i]) for i in range(2)]
#        print(corner, reso)
        nf = not redraw
        self.pltaxe.clear()
        datasel.display(zoom=zz, scale=self.scale, show=False, figure=self.pltaxe, new_fig=nf)
        self.pltaxe.text(corner[1],corner[0],"R: %.0fx%.0f"%tuple(reso))

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
        self.thresh = widgets.FloatSlider(value=50.0,
            min=1, max=100.0, step=0.1, layout=Layout(width='30%'),
            continuous_update=False, readout=True, readout_format='.1f')
        self.thresh.observe(self.pickpeak)
        self.peak_mode = widgets.Dropdown(options=['marker', 'bar'],value='marker',description='show as')
        self.peak_mode.observe(self.display)
        self.bexport = widgets.Button(description="Export",layout=Layout(width='10%'),
                button_style='success', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='Export to csv file')
        self.bexport.on_click(self.pkexport)
        self.bprint = widgets.Button(description="Print", layout=Layout(width='10%'),
                button_style='success', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='Print to screen')
        self.bprint.on_click(self.pkprint)
        self.spec = Output(layout={'border': '1px solid black'})
        self.out = Output(layout={'border': '1px solid red'})
        display( VBox([self.zoom,
            HBox([Label('threshold (x noise level):'), self.thresh, self.peak_mode, self.bprint, self.bexport] ),
            self.spec, self.out ]) )
        with self.spec:
            self.fig, self.ax = plt.subplots()
            self.npkd.display(new_fig=False, figure=self.ax, zoom=self.zoom.value)

    def pkprint(self,event):
        self.out.clear_output(wait=True)
        with self.out:
            print(self.pklist())
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
    def display(self, event):
        "interactive wrapper to peakpick"
        if event['name']=='value':
            self.ax.clear()
            self.npkd.display(new_fig=False, figure=self.ax, zoom=self.zoom.value)
            try:
                self.npkd.display_peaks(peak_label=True, figure=self.ax, zoom=self.zoom.value)
            except:
                pass
    def pickpeak(self, event):
        "interactive wrapper to peakpick"
        if event['name']=='value':
            self.pp()
    def pp(self):
        "do the peak-picking calling pp().centroid()"
        #self.spec.clear_output(wait=True)
        with self.spec:
            self.npkd.set_unit('m/z').peakpick(autothresh=self.thresh.value, verbose=False, zoom=self.zoom.value).centroid()
            self.ax.clear()
            self.npkd.display(new_fig=False, figure=self.ax, zoom=self.zoom.value)
            x = self.zoom.value
            y = [self.npkd.peaks.threshold]*2
            self.ax.plot(x,y,':r')
            self.npkd.display_peaks(peak_label=True, peak_mode=self.peak_mode.value,figure=self.ax, zoom=self.zoom.value)
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

