from __future__ import print_function
from os.path import *

import os
import shutil
import Tkinter
import tkFileDialog as Selector
import sys
from os.path import join

"""
Passer en Enthought avant de lancer le code.
Routine to recuperate the .py code of the interface created with QtDesigner as an .ui file.
It places the .py code next to the .py file we want to interface.
The addresses of the .py file to be interfaced and the .ui are asked in this order with a Tkinter interface.
"""
gene_interf = False
igbmc = True
home = False
 
class py2ui(object):
   def __init__(self, name_py = None, name_ui = None, extf = None ):
      self.addressUi = os.getcwd()
      self.extf = extf
      self.path_py = name_py
      self.path_ui = name_ui
   
   def make_py_from_ui(self):
      '''
      Create the py module to be import for the User Intrface
      '''
      if igbmc:
         addpyuic = '/Library/Frameworks/EPD64.framework/Versions/7.1/lib/python2.7/site-packages/PyQt4/uic/pyuic.py'
      if home or gene_interf:
         addpyuic = '/usr/lib/python2.7/dist-packages/PyQt4/uic/pyuic.py'     # PyQt conversion routine from ui to py
      dir_py = os.path.dirname(self.path_py)
      name_py = os.path.basename(self.path_py)
      #self.add_py = self.path_py[:-3] + self.extf # Python address for the python resulting interface module, name of python routine to be interface+"Ui"
      self.add_py = join(dir_py,"Visu", "init", name_py[:-3] + self.extf)
      print("self.add_py ",self.add_py)
      os.system('python ' + addpyuic + ' ' + self.path_ui+ '> ' + self.add_py)# create .py file for the interface next to the routine to be interfaced.
   
   def name_path_py(self):
      currentpath = sys.path[0] # address for the program to interface
      self.path_py = Selector.askopenfilename(title = 'enter name of the main .py file', initialdir = currentpath)
      print("pathfile ", self.path_py)
   
   def name_path_ui(self):
      self.path_ui = Selector.askopenfilename(title = 'enter name of the .ui file', initialdir = self.addressUi)
   
   def findpath(self):
      ''' 
      path of .py file and .ui file with Tkinter window asking
      '''
      if not self.path_py and not self.path_py:
         self.name_path_py()
         self.name_path_ui()
   
   def open_files(self):
      self.f = open(self.add_py, 'r')
      self.name_corr = self.path_py[:-3] + '_corr' + self.extf
      self.g = open(self.name_corr , 'w')
   
   def close_files_and_rename(self):
      self.f.close()
      self.g.close()
      os.rename(self.name_corr, self.add_py)
   
   def test_change_line(self, line):
      '''
      Correcting lines
      '''
      list_err_line = ['from PyQt4 import QtCore, QtGui']
      for elem in list_err_line:
         if elem in line:
            print("elem found is ", elem)
            #self.g.write('from PySide import QtCore, QtGui')
            self.g.write('from ..Pyside_PyQt4 import*')
            return True
      return False
   
   def test_remove_line(self, line):
      list_rem = ['self.verticalLayout.setMargin(0)',
                  'self.layoutD.setMargin(0)',
                  'self.horizontalLayout.setMargin(0)',
                  'self.gridLayout.setMargin(0)',
                  'self.gridLayout_2.setMargin(0)',
                  'self.horizontalLayout_2.setMargin(0)',
                  'self.verticalLayoutWidget.setGeometry(QtCore.QRect(250, 90, 621, 541))'
                  'self.layoutDWidget.setGeometry(QtCore.QRect(250, 90, 621, 541))']
                  #'self.pushButton_9 = QtGui.QPushButton(self.centralwidget)',
                  #'self.pushButton_9.setGeometry(QtCore.QRect(180, 420, 71, 61))',
                  #'self.pushButton_9.setObjectName(_fromUtf8("pushButton_9"))',
                  #'self.pushButton_9.setText(QtGui.QApplication.translate("MainWindow", "PushButton", None, QtGui.QApplication.UnicodeUTF8))']
      for elem in list_rem:
         if elem in line:
            print("elem found is ", elem)
            self.g.write('')
            return True
      return False
         
   def correc_py(self):
      self.open_files()
      for line in self.f.readlines():
         test = False
         if self.test_change_line(line) or self.test_remove_line(line):
            test = True
         #print "test ", test
         if not test:
            self.g.write(line)
      self.close_files_and_rename()
         
if __name__ == '__main__':

   if gene_interf:
      name_py = '/home/lio/Bureau/Encours/General_Interf/interf.py'
      name_ui = '/home/lio/Bureau/Encours/General_Interf/interf.ui'
   elif igbmc:
      name_py = '/Users/chiron/bitbuck/draft/fticrvisu.py'
      name_ui = '/Users/chiron/bitbuck/draft/Visu/util/fticrResol.ui'
   elif home:
      name_py = '/home/lio/Documents/bitbucket/draft/fticrvisu.py'
      name_ui = '/home/lio/Documents/bitbucket/draft/Visu/util/fticrResol.ui'
   pyui = py2ui(name_py = name_py, name_ui = name_ui, extf = 'Ui.py')
   pyui.findpath()
   pyui.make_py_from_ui()#
   pyui.correc_py()

