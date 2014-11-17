import os

from matplotlib import pyplot as plt

try: 
    #import PySide
    import PySide.QtGui 
    pyside = True
    print "using PySide for visu2D"
except:
    try:
        import PyQt4.QtGui 
        pyside = False
        print "using PyQt4 for visu2D"
    except:
        print "neither PySide nor PyQt4 works."
    
try:
    if pyside:
        from PySide.QtGui import QMainWindow, QApplication, QWidget, QToolButton, QAction, QVBoxLayout
        os.environ["QT_API"] = "pyside"
        from PySide import QtCore, QtGui
        from PySide.QtCore  import SIGNAL 
        from PySide.QtCore import QObject as Qobj
    else:
        from PyQt4.QtGui import QMainWindow, QApplication, QWidget, QToolButton, QAction, QVBoxLayout
        from PyQt4 import QtCore, QtGui
        from PyQt4.QtCore  import SIGNAL 
        from PyQt4.QtCore import QObject as Qobj
    try:
        _fromUtf8 = QtCore.QString.fromUtf8
    except AttributeError:
        _fromUtf8 = lambda s: s
except:
    print "neither Pyside not PyQt4 are imported."

