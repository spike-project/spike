from os.path import *
import os
import shutil
import Tkinter
import tkFileDialog as Selector
import sys

"""
Routine to recuperate the .py code of the interface created with QtDesigner as an .ui file.
It places the .py code next to the .py file we want to interface.
The addresses of the .py file to be interfaced and the .ui are asked in this order with a Tkinter interface.

"""
addressUi = os.getcwd() #address where is the file.ui, here address of the generator.

popup = False

def recupUipyth(pthpy,pthui,extf):
   '''
   Create the py module to be import for the User Intrface
   '''
   
   addpyuic='/usr/lib/python2.7/dist-packages/PyQt4/uic/pyuic.py'     # PyQt conversion routine from ui to py
   addpy=pthpy[:-3] + extf # Python address for the python resulting interface module, name of python routine to be interface+"Ui"
   os.system('python ' + addpyuic + ' ' + pthui + '> ' + addpy)# create .py file for the interface next to the routine to be interfaced.


def findpath(both):
    ''' 
    path with Tkinter window asking
    '''
    
    if both==True :
        currentpath=sys.path[0]# address for the program to interface
        # recuperate the address of program to interface
        pathtrgpy=Selector.askopenfilename(title='enter the name of the program',initialdir=currentpath)
        print "pathfile ", pathtrgpy
    currentpath = addressUi# address of the Ui files
    print currentpath
    # recuperate the address of the interface module
    pathui = Selector.askopenfilename(title='enter name of the Ui files',initialdir=currentpath)
    print "pathui ", pathui
    print "current path ", os.getcwd()
    print 'base ',os.getcwd()+'/'+os.path.basename(pathui)
    local_name_ui = os.getcwd()+'/'+os.path.basename(pathui)
    if pathui != local_name_ui :
    	shutil.copyfile(pathui, local_name_ui)
    if both==True :
        return pathtrgpy, pathui
    else :
        return pathui

pathtrgpy,pathui = findpath(both = True)
recupUipyth(pathtrgpy,pathui,'Ui.py')# create the user interface at the address of the file for which is imported the interface module.

if popup :
    print "enter the pop up"
    pathui = findpath(both=False)
    recupUipyth(pathtrgpy,pathui,'popUi.py')# create the user interface at the address of the file for which is imported the interface module.
