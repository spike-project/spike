#!/usr/bin/env python 
# encoding: utf-8

from __future__ import print_function
from .. util.debug_tools import* 
from .. Visu.Matplotlib_generictools import*
from .. Visu.Pyside_PyQt4 import*
import unittest

@dec_class_pr
@decclassdebugging
class DISPLAY():# class regrouping tools to change resolution with zoom
    '''
    Fill Canvas with Matplotlib figures.
    Handle connect disconnect.
    '''
    def __init__(self, QtMplCv, data, interface, paramz):
        self.QtMplCv = QtMplCv # Matplotlib Canvas
        self.data = data
        self.interface = interface
        self.currentd = self.data.d[len(self.data.d)-1]                       # current resolution is the lowest one.
        self.currsize = [self.currentd.size2, self.currentd.size1]          # list of size 1 and 2 of current resolution
        self.levfact = 4                                                    # level factor
        self.paramz = paramz                                                # 
        self.list_res()                                                     # list of the resolutions
        self.AREAMIN = 800000                                               # minimum number of pixels to be displayed
        self.curs_type = {'cross': QtCore.Qt.CrossCursor,
                    'hand': QtCore.Qt.OpenHandCursor,
                    'drag': QtCore.Qt.OpenHandCursor,
                    'zoom': QtCore.Qt.CrossCursor,
                    'cursor': QtCore.Qt.ArrowCursor,
                    'manual_profile': QtCore.Qt.IBeamCursor}
                    
    def disconnect(self, object_action, window = 'C'):
        '''
        Disconnect the canvas C and D
        '''
        if window == 'C': self.qmc.fig.canvas.mpl_disconnect(object_action) 
        if window == 'D': self.qmd.fig.canvas.mpl_disconnect(object_action)             
    
    def connect(self, event, action, window = 'C'):
        '''
        Connect the canvas C and D
        '''
        if window == 'C': self.qmc.fig.canvas.mpl_connect(event, action)
        if window == 'D': self.qmd.fig.canvas.mpl_connect(event, action)             
      
    def setcursor(self, name, window = 'C'):
        '''
        Set the type of cursor used.
        '''
        if window == 'C': self.qmc.setCursor(QtGui.QCursor(self.curs_type[name]))             
        if window == 'D': self.qmd.setCursor(QtGui.QCursor(self.curs_type[name]))             
        
    def set_canvasC(self):
        '''
        Instantiates the canvas C
        '''
        if debug(self): print("in display.set_canvasC")
        self.qmc = self.QtMplCv(self.interface.ui.centralwidget, self.paramz) 
        self.setcursor('cross') # Sets the cursor in the cross mode
        
    def set_canvasD(self):
        '''
        Instantiates the canvas D
        '''
        if debug(self): print("in display.set_canvasD")
        self.qmd = self.QtMplCv(self.interface.ui.centralwidget, self.paramz) # recreate a Qt for Mpl canvas for D
        self.setcursor('hand', 'D')    # Sets the cursor in openhand mode.
        
    def list_res(self):
        '''
        Objects for navigating in the different resolutions.
        self.LISTR, self.LISTRMAX, self.LISTD, self.RESMAX, self.LISTDNAME
        '''
        if debug(self): print("in display.list_res")
        self.LISTR = ['resol'+ str(i) for i in xrange(self.data.NBRES, 0, -1)]
        self.LISTRMAX = self.LISTR[0]                       # maximal resolution
        self.LISTD = self.data.d                            # list of the different resolutions
        self.RESMAX = self.LISTD[self.LISTR.index(self.LISTRMAX)]
        self.LISTDNAME = ['d' + str(i) for i in xrange(self.data.NBRES, 0, -1)] # d1 biggest resolution..
    
    def distrib(self, f, arg): # Use to make a more compact writing for llx,lly,urx,ury,
        '''
        Applying f to pairs of arg.
        '''
        if debug(self):  print("in display.distrib")
        return f(arg[0], arg[2]), f(arg[1], arg[3])
    
    def multiply_zoom_coord(self, alpha, beta):
        '''
        Changes the coordinates of the window with different
        factors (alpha, beta) according to the direction.
        '''
        zc = self.paramz.zoom_coord
        return [alpha*zc[0], beta*zc[1], alpha*zc[2], beta*zc[3]]
    
    def prepare_zoom(self, d, make_zoom):
        '''
        Calculates self.zoomtup before display with self.affichd()
        '''
        if len(self.paramz.zoom_coord) == 0 or not make_zoom : # If no coordinates or no zoom asked.. 
            coord = self.paramz.zoom_coord_old                      # Retrieves old coordinates.
            self.paramz.zoom_coord = coord
        else :
            coord = self.paramz.zoom_coord[:4]                     # Takes the current coordinates
        if d.axis1.currentunit == 'm/z':
            coord[1::2], coord[0::2] = coord[1::2][::-1], coord[0::2][::-1]                     # if m/z inverse the order
        self.zoomtup = (list(d.axis1.itoc(coord[1::2])), list(d.axis2.itoc(coord[0::2])))       # Transforming from point to current
        if d.axis1.currentunit == 'm/z':                    # Limit the m/z view with highmass values. 
            if self.zoomtup[0][1] > d.axis1.highmass: self.zoomtup[0][1] = d.axis1.highmass
            if self.zoomtup[1][1] > d.axis2.highmass: self.zoomtup[1][1] = d.axis2.highmass
    
    def affichd(self, canvas, d, make_zoom = True):             # print data in window C d is a fticrdata
        '''
        Making the display with NPKData display method. 
        '''
        self.prepare_zoom(d, make_zoom)                                 # Calculates self.zoomtup
        d.display(figure = canvas.axes, zoom = self.zoomtup, scale = self.paramz.scale, absmax = d.absmax) #  Display the zoom with current resolution

    def affi(self, canvas, d): # print data in window D
        '''
        Routine for printing the 2d in the D Window.
        It uses the coordinates of the zoom for a given resolution.
        '''
        if debug(self): print("in display.affi")
        d.display( figure = canvas.axes,  xlabel = "", ylabel = "")

    def change_resolution(self, layout1 = None, layout2 = None): # change resolution keeping the zoom.
        '''
        Changes the resolution and changing coordinates accordingly.
        self.currentd is selected according to self.data.resolu.
        '''
        if debug(self): print("in display.resol")
        for reso in self.LISTR:                                                                             # Loop on the resolutions. Changing resolution and zoom
            if self.data.resolu == reso :
                dd = self.LISTD[self.LISTR.index(reso)]                                                     # Retrieves the corresponding resolution
                rapp2, rapp1 = dd.size2/float(self.currentd.size2), dd.size1/float(self.currentd.size1)     # Calculates the ratio between old resolution and new resolution      
                self.paramz.zoom_coord = self.multiply_zoom_coord(rapp2, rapp1)                             # Changes the coordinates accordingly
                self.currentd = dd                                                                          # Changing the currrent resolution
                self.affd(self.currentd, self.data.resmin, layout1, layout2)                                # Visualizes with the current resolution
                llx , lly = self.distrib(min, self.paramz.zoom_coord)
                urx, ury  = self.distrib(max, self.paramz.zoom_coord)
                nbpix = (ury-lly)*(urx-llx)                                                                 # Number of pixels
    
    def right_order_coord(self, llx, lly, urx, ury):
        '''
        Restablishes right order for coordinates.
        '''
        llx, urx = min(llx, urx), max(llx, urx) # sort for making ptlx 
        lly, ury = min(lly, ury), max(lly, ury) # sort for making ptly
        return llx, lly, urx, ury
    
    def local_abs_max(self):
        '''
        Maximum peak height extracted from local view (self.zoomtup).
        '''
        from ..plugins.Peaks import peaks2d 
        x, y, z = peaks2d(self.res2dd(), threshold=0.1, zoom = self.zoomtup)
        self.interface.ui.label.setText("{0:.2e}".format(int(z.max()))) # shows absmax in the interface 

    def affd(self, d1, d2 , layout1, layout2 = None, make_zoom = True, message = None):# routine to visualize contours of 2D data.
        """
        Routine to show the Mpl in qt zoom in main window and global window
        d1 and d2 are the data to be plotted.
        """
        if debug(self): print("in display.affd")
        if debug(self): print("in display self.interface.clearlayout(self.interface.ui.layoutC)")
        self.interface.clearlayout(self.interface.ui.layoutC)       # Clear the data window C
        self.set_canvasC()                                          # makes the canvas C
        lay1 = self.interface.ui.layoutC.addWidget(self.qmc)
        if message:
            print("self.paramz.zoom_coord[2] ", self.paramz.zoom_coord[2])
            x, y = self.paramz.zoom_coord[3], self.paramz.zoom_coord[2]
            self.message(message, posx = x,  posy = y)
        if debug(self): print("in affd, using zoom ", zoom)
        self.affichd(self.qmc, d1, make_zoom = make_zoom)           # window C
        self.local_abs_max()                                    #  Indicates the maximal height of peaks in the local view. 
        if self.interface.ui.layoutD is not None :
            self.interface.clearlayout(self.interface.ui.layoutD)           # Clears the data window D
            self.set_canvasD()                                              # makes the canvas D
            self.interface.ui.layoutD.addWidget(self.qmd)
            self.affi(self.qmd, d2)                                         # window D

    def res2dd(self):
        '''
        Method for loading the resolution from the current self.data.resolu name
        '''
        self.currentd = self.LISTD[self.LISTR.index(self.data.resolu)] # Retrieves resolution corresponding to self.data.resolu
        return self.currentd
    
    def aff_listview_index(self):
        '''
        Displays the view at position self.paramz.listview_index.
        '''
        pos = self.paramz.listview_index                                        # Current index in self.paramz.listview
        self.paramz.zoom_coord = self.paramz.listview[pos]['zoom']              # Retrieves old coordinates
        self.paramz.scale = self.paramz.listview[pos]['scale']                  # Retrieves old scale
        self.data.resolu = self.paramz.listview[pos]['resolution']              # Retrieves old resolution
        self.paramz.sliderpos = self.paramz.listview[pos]['slider']             # Retrieves old slider position
        #self.interface.ui.horizontalSlider.setValue(self.sliderpos)
        self.affd(self.res2dd(), self.data.resmin, self.interface.ui.layoutC, self.interface.ui.layoutD) # Visualizes the current resolution in C and D
        self.aff_resolution()                                                   # Shows resolution in the Interface

    def zoom_area(self):
        '''
        Calculates area from zoom coordinates.
        '''
        zc = self.paramz.zoom_coord
        return (zc[0] - zc[2])*(zc[1] - zc[3])

    def register_coordinates(self):
        '''
        Stores in "self.paramz.listview" the zooms coordinates, the resolution, the scale and the slider position
        '''
        if debug(self): print("register coordinates ", self.paramz.zoom_coord)
        self.paramz.listview = self.paramz.listview[:int(self.paramz.listview_index) + 1]           # cuts self.paramz.listview current element
        self.paramz.listview += [{'zoom' : self.paramz.zoom_coord,
                                  'resolution': self.data.resolu,
                                  'scale' : self.paramz.scale,
                                  'slider': self.paramz.sliderpos
                                 }]                # keeps the zoom, the resolution
        self.paramz.listview_index += 1  # Increments the index of the list of views.

    def select_best_resolution(self):
        '''
        According to the size of the zoom chose the first resolution for which the visualized array is not too big.
        '''
        cnd = lambda : self.zoom_area() < self.AREAMIN and self.data.resolu != self.LISTRMAX    # If area too small and not resolution max..
        if  cnd(): 
            while cnd() :
                self.data.resolu = self.LISTR[self.LISTR.index(self.data.resolu)-1]         # Taking next resolution
                self.change_resolution(layout1 = self.interface.ui.layoutC, layout2 = self.interface.ui.layoutD) # Calculates new zoom with new resolution
        else :                                      # With the current resolution, the zoom size is OK. 
            self.change_resolution(layout1 = self.interface.ui.layoutC, layout2 = self.interface.ui.layoutD)
    
    def aff_resolution(self):
        '''
        Shows the resolution in the interface.
        '''
        if debug(self): print("in display.aff_resolution")
        curr_numres = str(self.data.NBRES - self.LISTR.index(self.data.resolu))
        self.interface.ui.label_2.setText("res :" + curr_numres +'/' + str(self.data.NBRES))         # setting the resolution in the ui
        self.interface.ui.lineEdit_3.setText(str(self.paramz.scale))    # setting the resolution in the ui
    
    def message(self, message, posx = None, posy = None, colorlab = (1.0, 0.7, 0.7)):
        '''
        function to label the peaks
        '''
        el = Ellipse((2, -1), 0.5, 0.5)
        self.qmc.axes.add_patch(el)
        color = colorlab
        if posx and posy :          
            print("posx, posy ", posx, posy)                                                                                
            ann = self.qmc.axes.annotate(message, xy = (posy, posx),  #xycoords='data',
                        xytext = (posy, posx), textcoords = 'offset points', size = 10,# va="center", #xytext=(posx, posy), 
                        bbox = dict(boxstyle = "round", fc = color, ec = "none"),)
                 
class SavingTests(unittest.TestCase):

    def test_display(self):
        "Testing display module"
        from .. Visu.canvas import Qt4MplCanvas as QtMplCv     # class for using matplotlib in Qt canvas.
        from .. Visu.Load import LOAD
        from .. Visu.interface import INTERFACE
        from .. Visu.paramzoom import PARAM_ZOOM
        data = LOAD(configfile = 'spike/Visu/visu2d_eg.mscf')
        interf = INTERFACE()
        paramz = PARAM_ZOOM(data)
        display = DISPLAY(QtMplCv, data, interf, paramz)  
        self.assertEqual(display.levfact, 4)
        self.assertEqual(display.AREAMIN, 800000)
        self.assertIsInstance(display.interface, INTERFACE)
        self.assertIsInstance(display.paramz, PARAM_ZOOM)

                                                  
if __name__ == '__main__':
    unittest.main()
         