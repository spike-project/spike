#!/usr/bin/env python 
# encoding: utf-8

"""
A tool to display 2D FT-ICR data-sets

to be embedded in jupyter notebook

First version MAD jan 2019
preliminary and not fully tested !

"""

import tables
import numpy as np

from .. import FTICR
from ..NPKData import flatten, parsezoom

import matplotlib.pyplot as plt
from ipywidgets import interact, interactive, fixed, interact_manual, HBox, VBox, HTML, Label, Layout
import ipywidgets as widgets
from IPython.display import display

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
        self.col = self.data[0].col
        self.row = self.data[0].row
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
        super(MR_interact, self).__init__(name, report=report, Debug=Debug)
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
