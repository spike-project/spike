'''
Created by Lionel Chiron  02/10/2013 
Copyright (c) 2013 __NMRTEC__. All rights reserved.
'''
import sys, os
from interfUi import Ui_MainWindow as Ui
from PySide import QtCore, QtGui
from PySide.QtCore  import SIGNAL 
#import the NavigationToolbar Qt4Agg widget
from matplotlib.backends.backend_qt4 import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.patches import PathPatch, Ellipse, Rectangle
import matplotlib.lines as mlines
from matplotlib.pyplot import figure, show
from matplotlib.figure import Figure #
from matplotlib import pyplot as plt

class Qt4MplCanvas(FigureCanvas): 
    """
    Class for integrating Matplotlib in PySide
    Drop events capabilities.
    """ 
    def __init__(self, parent):
        self.fig = Figure(figsize=(16, 4))
        self.axes = self.fig.add_subplot(111)
        FigureCanvas.__init__(self, self.fig) # set the parent widget 
        self.setParent(parent)
        FigureCanvas.updateGeometry(self)
        try:
            hasattr(self, droplist)
            print "self.droplist ",self.droplist
        except Exception:
            self.droplist = []

    def dropEvent(self, event):
        '''
        If drop event is an url, append the address to self.droplist
        '''
        if event.mimeData().hasUrls():
            addr_file = str(event.mimeData().urls()[0].toLocalFile())
        self.droplist.append(addr_file)
        self.emit(QtCore.SIGNAL("dropped"))
        
    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.accept()
             
class interf(): #interface
    '''
    Creation of the graphic window
    Methods :
        init_interf
        makelayout
        makecanvas
        clearlayout
        makecontextmenu
        maketoolbar
        copyclipboard
        run
    '''
    def __init__(self):
        '''
        Two ways for launching Qt application, normal way and Canopy way.
        '''
        try:
            whichapp =  'app is QtGui.QApplication(sys.argv)'
            self.app = QtGui.QApplication(sys.argv)
        except :
            whichapp = 'app is QtCore.QCoreApplication.instance()'
            self.app = QtCore.QCoreApplication.instance() # Correction for Canopy 
        self.window = QtGui.QMainWindow()
        self.ui = Ui()# 
        self.ui.setupUi(self.window)#set up of the main window
        self.show_toolbar = True
        self.init_interf()
        
    def init_interf(self):
        '''
        Initialization of the interface
        '''
        self.window.setCursor(QtCore.Qt.ArrowCursor)
        try:
            self.clearlayout()
        except :
            self.makelayout()
        self.makecanvas()
        self.makecontextmenu()
        if self.show_toolbar:
            self.maketoolbar()
        self.makemenubar()
        
    def makemenubar(self):
        '''
        Add elements to the menubar
        '''
        self.menubar = self.window.menuBar()
        
        #self.menubar = QtGui.QMenuBar()
        #menu = QtGui.QMenu(self.menubar)
        #menu.setTitle("Testing")
        
        #print dir(self.menubar)
        if not hasattr(self,'menubar_file'):
            self.menubar_file = self.menubar.addMenu('&File')
        if not hasattr(self,'menubar_edit'):
            self.menubar_edit = self.menubar.addMenu('&Edit')
        if not hasattr(self,'menubar_tools'):
            self.menubar_tools = self.menubar.addMenu('&Tools')
        if not hasattr(self,'menubar_navigation'):
            self.menubar_navigation = self.menubar.addMenu('&Navigation')
        if not hasattr(self,'menubar_help'):
            self.menubar_help = self.menubar.addMenu('&Help')
        
        #self.window.setMenuBar(self.menubar)
        
    def makelayout(self):
        '''
        make the layout in the centralwidget
        '''
        print "makes and fills the layout"
        self.layout = QtGui.QVBoxLayout(self.ui.centralwidget)#creation of the layout widget.
        self.layout.setContentsMargins(280, 10, 5, 70)# Margin for the picture left,right, up, dow
        
    def makecanvas(self):
        '''
        Put the canvas in the layout
        '''
        self.canvas = Qt4MplCanvas(self.ui.centralwidget)#
        self.canvas.setAcceptDrops(True)
        self.layout.addWidget(self.canvas)
        self.canvas.setCursor(QtCore.Qt.OpenHandCursor)# make the cursor in openhand mode.
        
    def clearlayout(self):#   Delete item widgets in the layout to allow the refresh.
        '''
        Clear the layout of the centralwidget
        '''
        print "clear layout"
        while self.layout.count() > 0:
            item = self.layout.takeAt(0)
            if not item: continue
            w = item.widget()
            if w: w.deleteLater()

    def context_menu(self,event):
        '''
        context menu
        '''
        self.menu.popup(self.canvas.mapToGlobal(event))
        
    def makecontextmenu(self):
        '''
        Creates the context menu
        '''
        self.menu = QtGui.QMenu(self.canvas)
        self.canvas.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.canvas.customContextMenuRequested.connect(self.context_menu)
        
    def maketoolbar(self):
        '''
        using the native toolbar
        '''
        self.navi_toolbar = NavigationToolbar(self.canvas, self.window)
        self.layout.addWidget(self.navi_toolbar)
    
    def copy_clipboard(self):
        '''
        Copy the canvas picture to clipboard
        '''
        print "copy to clipboard"
        name_pic = "clip.jpg"
        clipboard = self.app.clipboard()
        self.canvas.fig.savefig(name_pic)
        loadedImage = QtGui.QImage()
        loadedImage.load(name_pic)
        clipboard.setImage(loadedImage)
        os.remove(name_pic)
        print "copy done"
    
    def run(self):
        self.window.show()
        sys.exit(self.app.exec_()) # Correction for Canopy
        