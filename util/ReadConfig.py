
# encoding: utf-8
"""
ReadConfig.py
Small programm which reads, presents and adjust configuration file

Created by mac on 2012-03-05.
Copyright (c) 2012 __NMRTEC__. All rights reserved.
"""

import sys
import os
from configobj import ConfigObj     # standard dans Enthought !!!
from PySide import QtGui
from PySide import QtCore
import SetConfig

    
class readConfig():
    """
    Read configuration from filename
    """
    def __init__(self, filename):
        self.filename = filename
        self.cp = self.get_parser()
        
    def get_parser(self):
        config = ConfigObj(self.filename)  # load it as a dictionnary
        print config.keys()
        return config

class InfoWindow(QtGui.QWidget):
    """
    Window presenting the configuration read : 
    Two columns per section, one for the label and one for the line edit.
    Each LineEdit is connected to adjust_config which builds the dictionnary of modification
    """
    def __init__(self,filename):
        super(InfoWindow, self).__init__()
        self.filename = filename
        self.config = readConfig(self.filename)
        self.modif = {}
        self.create_widgets()
    def create_widgets(self):
        """
        Presenting informations
        """
        ok_button = QtGui.QPushButton("Use")
        cancel_button = QtGui.QPushButton("Cancel")
        # 
        gridall = QtGui.QGridLayout()
        gridall.setSpacing(10)

        sections = self.config.cp.keys()
        self.tabWidget = QtGui.QTabWidget(self)
        self.tabWidget.setGeometry(0,20,400,400)
        for j in range(len(sections)):
            self.pn = QtGui.QWidget(self.tabWidget)
            items = self.config.cp[sections[j]]
            grid = QtGui.QGridLayout(self.pn)
            for i in range(len(items)):
                item_label = QtGui.QLabel(items.keys()[i])
                self.item_edit = QtGui.QLineEdit(text = items.values()[i])
                self.item_edit.textChanged.connect(self.adjust_config)                  # when text is changed, adjust_config is called 
                self.item_edit.setObjectName(items.keys()[i])
                grid.addWidget(item_label,i,0)
                grid.addWidget(self.item_edit,i,1)
            self.pn.setLayout(grid)
            self.tabWidget.addTab(self.pn,sections[j])
        self.setWindowTitle('Read Configuration')
        gridall.addWidget(self.tabWidget,0,1,0,2)
        gridall.addWidget(ok_button,1,1)
        QtCore.QObject.connect(ok_button, QtCore.SIGNAL("clicked()"), self.on_ok_clicked)
        gridall.addWidget(cancel_button,1,2)
        QtCore.QObject.connect(cancel_button, QtCore.SIGNAL("clicked()"), self.on_cancel_clicked)
        self.setLayout(gridall)
        self.show() 
    def on_ok_clicked(self):
        SetConfig.set_config(self.modif,self.filename)
        self.close()
    def on_cancel_clicked(self):
        self.close()
    # def closeEvent(self,event): 
    #     """
    #     On close
    #     """
    #     reply = QtGui.QMessageBox.question(self
    #         , 'Close Confirmation'
    #         , 'Do you want to use this configuration ? :'
    #         , QtGui.QMessageBox.Yes | QtGui.QMessageBox.No
    #         , QtGui.QMessageBox.No)
    #     if reply == QtGui.QMessageBox.Yes:
    #         SetConfig.set_config(self.modif,self.filename)      # update the configuration file
            
    def adjust_config(self,text):
        """
        build the dictionnary storing modification
        """
        sender = self.sender()
        key = sender.objectName()
        self.modif[key] = text
        

   
if __name__ == '__main__':
    #config = getConfig("/Volumes/XeonData/Developement/Novalix/sscreen/ConfigFile.txt")
    app = QtGui.QApplication(sys.argv)
    main_window = InfoWindow("process.mscf")
    #main_window = InfoWindow("/Volumes/XeonData/Developement/NPK_v2/draft/visu2D.mscf")
    sys.exit(app.exec_())
