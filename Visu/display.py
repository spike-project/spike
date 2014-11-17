from util.debug_tools import* 
from Matplotlib_generictools import*
from Pyside_PyQt4 import*

@dec_class_pr
@decclassdebugging
class DISPLAY():# class regrouping tools to change resolution with zoom
    '''
    Fill Canvas with Matplotlib figures.
    Handle connect disconnect.
    
    '''
    def __init__(self, QtMplCv, data, interface, paramz):
        self.QtMplCv = QtMplCv
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
        Disconnect
        '''
        if window == 'C':
            self.qmc.fig.canvas.mpl_disconnect(object_action) 
        if window == 'D':
            self.qmd.fig.canvas.mpl_disconnect(object_action)             # make the cursor in cross mode in window C
    
    def connect(self, event, action, window = 'C'):
        '''
        Connect
        '''
        if window == 'C':
            self.qmc.fig.canvas.mpl_connect(event, action)
        if window == 'D':
            self.qmd.fig.canvas.mpl_connect(event, action)             # make the cursor in cross mode in window C 
      
    def setcursor(self, name, window = 'C'):
        '''
        Set the type of cursor used.
        '''
      
        if window == 'C':
            self.qmc.setCursor(QtGui.QCursor(self.curs_type[name]))             # make the cursor in cross mode in window C
        if window == 'D':
            self.qmd.setCursor(QtGui.QCursor(self.curs_type[name]))             # make the cursor in cross mode in window C
        
    def set_canvasC(self):
        '''
        Makes the canvas C
        '''
        if debug(self):
            print "in display.set_canvasC"
        self.qmc = self.QtMplCv(self.interface.ui.centralwidget, self.paramz) # recreate a Qt for Mpl canvas for C
        self.setcursor('cross')
        
    def set_canvasD(self):
        '''
        Makes the canvas D
        '''
        if debug(self):
            print "in display.set_canvasD"
        self.qmd = self.QtMplCv(self.interface.ui.centralwidget, self.paramz) # recreate a Qt for Mpl canvas for D
        self.setcursor('hand', 'D')    # make the cursor in openhand mode.
        
    def list_res(self):
        '''
        Make list with the resolutions
        '''
        if debug(self):
            print "in display.list_res"
        self.LISTR = ['resol'+ str(i) for i in xrange(self.data.NBRES, 0, -1)]
        self.LISTRMAX = self.LISTR[0] # maximal resolution
        self.LISTD = self.data.d # list of the different resolutions
        self.RESMAX = self.LISTD[self.LISTR.index(self.LISTRMAX)]
        self.LISTDNAME = ['d' + str(i) for i in xrange(self.data.NBRES, 0, -1)] # d1 biggest resolution..
    
    def distrib(self, f, arg): # Use to make a more compact writing for llx,lly,urx,ury,
        '''
        Applying f to pairs of arg.
        '''
        if debug(self):
            print "in display.distrib"
        return f(arg[0], arg[2]), f(arg[1], arg[3])
    
    def multzoom_coord(self, alpha, beta):
        '''
        change the coordinates of the window with different
        factors (alpha, beta) according to the direction.
        '''
        if debug(self):
            print "in display.multzoom_coord"
        zc = self.paramz.zoom_coord
        return [alpha*zc[0], beta*zc[1], alpha*zc[2], beta*zc[3]]
    
    def affichd(self, canvas, d, zoom = True): # print data in window C d is a fticrdata
        '''
        Makes the display with NPKv2. 
        '''
        if debug(self):
            print "in display.affichd"
            print "in affichd, at beginning, self.paramz.zoom_coord ", self.paramz.zoom_coord
            print " zoom ", zoom 
        if len(self.paramz.zoom_coord) == 0 or zoom == False :
            print "using old coordinates"
            coord = self.paramz.zoom_coord_old
            self.paramz.zoom_coord = self.paramz.zoom_coord_old
        else :
            coord = self.paramz.zoom_coord
        if debug(self):
            print "in affichd, coord zoom are ", coord
        f1limits = (int(round(coord[1])), int(round(coord[3])))
        f2limits = (int(round(coord[0])), int(round(coord[2])))
        zoomtup = (f1limits, f2limits)
        d.display(figure = canvas.axes, zoom = zoomtup, scale = self.paramz.scale) # 

    def affi(self, canvas, d):# print data in window D
        '''
        Routine to print the 2d datas in Qt embedded environment
        It uses the coordinates of the zoom for a given resolution.
        '''
        if debug(self):
            print "in display.affi"
        d.display( figure = canvas.axes,  xlabel = "", ylabel = "") #, tick = False scale = res.scale,

    def change_resolution(self, layout1 = None, layout2 = None): # change resolution keeping the zoom.
        '''
        Changes the resolution.
        self.currentd is selected according to vis.resolu.
        '''
        if debug(self):
            print "in display.resol"
        for reso in self.LISTR:                                                         # loop on the resolutions. Change resolution and zoom
            if self.data.resolu == reso :
                dd = self.LISTD[self.LISTR.index(reso)]                                 # recuperation of the corresponding resolution
                # ration between old resolution and new resolution
                rapp2, rapp1 = (dd.size2/float(self.currentd.size2), dd.size1/float(self.currentd.size1))#        
                self.paramz.zoom_coord = self.multzoom_coord(rapp2, rapp1)              # changing the window with resolution
                zc = self.paramz.zoom_coord

                self.currentd = dd                                                      # current resolution is the resolution recuperated
                self.affd(self.currentd, self.data.resmin, layout1, layout2)             # visualize with the current resolution
                llx , lly = self.distrib(min, zc)
                urx, ury  = self.distrib(max, zc)
                nbpix = (ury-lly)*(urx-llx)                                             # number of pixels
    
    def right_order_coord(self, llx, lly, urx, ury):
        '''
        Restablishes right order for coordinates.
        '''
        llx, urx = min(llx, urx), max(llx, urx) # sort for making ptlx 
        lly, ury = min(lly, ury), max(lly, ury) # sort for making ptly
        return llx, lly, urx, ury
    
    def local_abs_max(self):
        '''
        maximum from local view
        '''
        d = self.res2dd()
        c = self.paramz.zoom_coord
        zc = self.paramz.zoom_coord
        x, y, z = d.peaks2d(zoom = [[zc[1], zc[3]], [zc[0], zc[2]]], value = True)
        if debug(self):
            print "maxi intensity is ", z.max()
        self.interface.ui.label.setText("{0:.2e}".format(int(z.max()))) # shows absmax in the interface 

    def affd(self, d1, d2 , layout1, layout2 = None, zoom = True, message = None):# routine to visualize contours of 2D data.
        """
        Routine to show the Mpl in qt zoom in main window and global window
        d1 and d2 are the data to be plotted.
        
        """
        if debug(self):
            print "in display.affd"
            print "in display self.interface.clearlayout(self.interface.ui.layoutC)"
        self.interface.clearlayout(self.interface.ui.layoutC) # Clear the data window C
        self.set_canvasC() # makes the canvas C
        lay1 = self.interface.ui.layoutC.addWidget(self.qmc)

        if message:
            print "self.paramz.zoom_coord[2] ", self.paramz.zoom_coord[2]
            x, y = self.paramz.zoom_coord[3], self.paramz.zoom_coord[2]
            #x, y = 100,100
            self.message(message, posx = x,  posy = y)
        if debug(self):
            print "in affd, using zoom ", zoom
        self.affichd(self.qmc, d1, zoom = zoom) # window C
        if debug(self):
            print "in display.affd show absmax"
        self.local_abs_max()
        if self.interface.ui.layoutD is not None :
            if debug(self):
                print "in display.affd"
                print "in display self.interface.clearlayout(self.interface.ui.layoutD)"
            self.interface.clearlayout(self.interface.ui.layoutD) # Clears the data window D
            self.set_canvasD() # makes the canvas D
            self.interface.ui.layoutD.addWidget(self.qmd)
            self.affi(self.qmd, d2) # window D

    def res2dd(self):
        '''
        Function to load the resolution from the current vis.resolu name
        '''
        if debug(self):
            print "in display.res2dd "
            print "vis.resolu ", vis.resolu
        dd = self.LISTD[self.LISTR.index(self.data.resolu)]#recuperation of the corresponding resolution
        self.currentd = dd # changing the resolution data
        return dd
    
    def afflistco(self):
        '''
        Prints the element of listco containing (zoom, resolution, scale)
        at position self.paramz.listview_index.
        '''
        if debug(self):
            print "in display.afflistco" 
        print "self.paramz.listview ", self.paramz.listview
        pos = self.paramz.listview_index
        self.paramz.zoom_coord = self.paramz.listview[pos]['zoom'] #retrieves old coordinates
        print "retrieved self.paramz.zoom_coord ",  self.paramz.zoom_coord
        self.paramz.scale = self.paramz.listview[pos]['scale'] #retrieves old scale
        print "retrieved self.paramz.scale ",  self.paramz.scale
        self.data.resolu = self.paramz.listview[pos]['resolution'] # retrieves old resolution
        self.paramz.sliderpos = self.paramz.listview[pos]['slider'] #retrieves old slider position
        #self.interface.ui.horizontalSlider.setValue(self.sliderpos)
        dd = self.res2dd() # loads current resolution
        self.affd(dd, self.data.resmin, self.interface.ui.layoutC,
                  self.interface.ui.layoutD) # visualizes the current resolution in C
        self.aff_resolution() # shows resolution in interface

    def zoom_area(self):
        '''
        Calculates area from zoom coordinates.
        '''
        zc = self.paramz.zoom_coord
        return (zc[0] - zc[2])*(zc[1] - zc[3])

    def register_coordinates(self):
        '''
        Keeps in a list "self.paramz.listview" the zooms coordinates and the associated resolutions.
        '''
        if debug(self):
            print "register coordinates ", self.paramz.zoom_coord
        self.paramz.listview = self.paramz.listview[:int(self.paramz.listview_index) + 1]           # cuts self.paramz.listview current element
        self.paramz.listview += [{'zoom' : self.paramz.zoom_coord,\
            'resolution': self.data.resolu, 'scale' : self.paramz.scale, 'slider': self.paramz.sliderpos }]                # keeps the zoom, the resolution
        self.paramz.listview_index += 1

    def select_best_resolution(self):
        '''
        In function of the size of the zoom chose the best resolution.
        '''
        if debug(self):
            print "in display.select_best_resolution"
        if self.zoom_area()  < self.AREAMIN and self.data.resolu != self.LISTRMAX :
            if debug(self):
                print "self.zoom_area()  < self.AREAMIN and self.data.resolu != self.LISTRMAX ",\
                 self.zoom_area()  < self.AREAMIN and self.data.resolu != self.LISTRMAX
            while self.zoom_area()  < self.AREAMIN and self.data.resolu != self.LISTRMAX :
                self.data.resolu = self.LISTR[self.LISTR.index(self.data.resolu)-1] # next resolution
                self.change_resolution(layout1 = self.interface.ui.layoutC,
                           layout2 = self.interface.ui.layoutD) # Calculates new zoom with new resolution
        else :
            print "in changeres_listzoom refresh"
            self.change_resolution(layout1 = self.interface.ui.layoutC, layout2 = self.interface.ui.layoutD)
    
    def aff_resolution(self):
        '''
        Shows the resolution in the interface.
        '''
        if debug(self):
            print "in display.aff_resolution"
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
            print "posx, posy ", posx, posy                                                                                
            ann = self.qmc.axes.annotate(message, xy = (posy, posx),  #xycoords='data',
                        xytext = (posy, posx), textcoords = 'offset points', size = 10,# va="center", #xytext=(posx, posy), 
                        bbox = dict(boxstyle = "round", fc = color, ec = "none"),)
                 
if __name__ == '__main__':
    print "hello"