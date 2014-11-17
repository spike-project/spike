'''
Created by Lionel Chiron  02/10/2013 
Copyright (c) 2013 __NMRTEC__. All rights reserved.
'''
import sys, os
from util.debug_tools import*  
from Visu.init.fticrvisuUi import Ui_MainWindow as Ui
from Visu.init.fticrvisupopUi import Ui_Form as Pop
from Pyside_PyQt4 import*

@dec_class_pr
@decclassdebugging
class INTERFACE(): #interface
    '''
    Creation of the graphic window
    Methods :
        init_interf
        makelayout
        makecanvas
        clearlayout
        run
    '''
    def __init__(self):
        '''
        Two ways for launching Qt application, normal way and Canopy way.
        '''
        print "in INTERFACE"
        try:
            whichapp =  'app is QtGui.QApplication(sys.argv)'
            self.app = QtGui.QApplication(sys.argv)
        except :
            whichapp = 'app is QtCore.QCoreApplication.instance()'
            self.app = QtCore.QCoreApplication.instance() # Correction for Canopy
        #print whichapp 
        self.window = QtGui.QMainWindow()
        self.ui = Ui()# 
        self.ui.setupUi(self.window)#set up of the main window
        #self.show_toolbar = False
        self.init_interf()
        
    def init_interf(self):
        '''
        Initialization of the interface
        '''
        self.window.setCursor(QtCore.Qt.ArrowCursor)
        self.makelayout()
           
    def makelayout(self):
        '''
        make the layout in the centralwidget
        '''
        print "makes and fills the layout"
        self.ui.layoutC = QtGui.QVBoxLayout(self.ui.centralwidget)#creation of the layout widget.
        self.ui.layoutC.setContentsMargins(280, 50, 10, 70) # Margin for the picture left, up, right, bottom
        
    def clearlayout(self, layout):#   Delete item widgets in the layout to allow the refresh.
        '''
        Clear the layout of the centralwidget
        '''
        if debug(self):
            print "in interface.clearlayout"
        while layout.count() > 0:
            item = layout.takeAt(0)
            if not item: continue
            w = item.widget()
            if w: w.deleteLater()

    def run(self):
        self.window.show()
        sys.exit(self.app.exec_()) # Correction for Canopy
        
