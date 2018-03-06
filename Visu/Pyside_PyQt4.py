#!/usr/bin/env python 
# encoding: utf-8

from __future__ import print_function

import os
from matplotlib import pyplot as plt

try: 
    #import PySide
    import PySide.QtGui 
    qtbib = 'pyside'
    print("using PySide for visu2D")
except ImportError:
    try:
        import PyQt4.QtGui 
        qtbib = 'pyqt4'
        print("using PyQt4 for visu2D")
    except ImportError:
        try:
            import PyQt5.QtGui 
            qtbib = 'pyqt5'
            print("using PyQt5 for visu2D")
        except ImportError:
            print("neither PySide nor PyQt4 nor PyQt5 work.")
    
try:
    if qtbib == 'pyside':
        from PySide.QtGui import QMainWindow, QApplication, QWidget, QToolButton, QAction, QVBoxLayout
        os.environ["QT_API"] = "pyside"
        from PySide import QtCore, QtGui
        from PySide.QtCore  import SIGNAL 
        from PySide.QtCore import QObject as Qobj
    elif qtbib == 'pyqt4' :
        from PyQt4.QtGui import QMainWindow, QApplication, QWidget, QToolButton, QAction, QVBoxLayout
        from PyQt4 import QtCore, QtGui
        from PyQt4.QtCore  import SIGNAL 
        from PyQt4.QtCore import QObject as Qobj
    elif qtbib == 'pyqt5' :
        from PyQt5.QtWidgets import QMainWindow, QApplication, QToolButton, QAction, QVBoxLayout
        # from PyQt5.QtGui import QMainWindow, QApplication, QWidget, QToolButton, QAction, QVBoxLayout
        from PyQt5 import QtCore, QtGui
        from PyQt5 import QtWidgets as QWidget 
        from PyQt5.QtCore  import SIGNAL 
        from PyQt5.QtCore import QObject as Qobj
    try:
        _fromUtf8 = QtCore.QString.fromUtf8
    except AttributeError:
        _fromUtf8 = lambda s: s
except:
    print("neither Pyside not PyQt4 nor PyQt5 are imported.")

