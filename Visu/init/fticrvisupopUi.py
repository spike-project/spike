# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '/Users/chiron/QtSDK/fticrpopup.ui'
#
# Created: Thu Oct 27 11:34:34 2011
#      by: PyQt4 UI code generator 4.8.4
#
# WARNING! All changes made in this file will be lost!


import sys
from ..Pyside_PyQt4 import*


try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s
    

class Ui_Form(object):
    def __init__(self):
        pass
        #loadbibqt()

    
    def setupUi(self, Form):
        #self.loadbibqt()
        Form.setObjectName(_fromUtf8("Form"))
        Form.resize(394, 388)
        self.treeView = QtGui.QTreeView(Form)
        self.treeView.setGeometry(QtCore.QRect(20, 30, 351, 311))
        self.treeView.setObjectName(_fromUtf8("treeView"))
        self.label = QtGui.QLabel(Form)
        self.label.setGeometry(QtCore.QRect(180, 10, 131, 16))
        self.label.setObjectName(_fromUtf8("label"))
        self.pushButton = QtGui.QPushButton(Form)
        self.pushButton.setGeometry(QtCore.QRect(150, 350, 101, 32))
        self.pushButton.setObjectName(_fromUtf8("pushButton"))
        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(QtGui.QApplication.translate("Form", "Form", None, QtGui.QApplication.UnicodeUTF8))
        self.treeView.setToolTip(QtGui.QApplication.translate("Form", "files", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("Form", "FTICR files", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton.setText(QtGui.QApplication.translate("Form", "Select", None, QtGui.QApplication.UnicodeUTF8))

