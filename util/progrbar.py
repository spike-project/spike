#!/usr/bin/env python 
# encoding: utf-8

from __future__ import print_function
from PyQt4 import QtCore, QtGui
from PyQt4.QtCore  import SIGNAL 
from PyQt4.QtCore import QObject as Qobj
from progrbarUi import Ui_MainWindow as Ui
import sys


def interfGraph(): #Initialisation of the user interface
    '''
    Creation of the graphic window
    ''' 
    app = QtGui.QApplication(sys.argv)
    progressbar = QtGui.QMainWindow()
    ui = Ui()# 
    ui.setupUi(progressbar)#set up of the main window
    
    # Defining the Signals
    #Qobj.connect(ui.pushButton, SIGNAL("clicked()"), zo.backzoo)# going back
    progressbar.show()
    app.exec_()

if __name__ == '__main__':
    interfGraph()#executing Main loop. 