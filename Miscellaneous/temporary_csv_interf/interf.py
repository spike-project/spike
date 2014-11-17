'''
Created by Lionel Chiron  02/10/2013 
Copyright (c) 2013 __NMRTEC__. All rights reserved.
'''
import os,sys
os.environ["QT_API"] = "pyside"
from orbitrap_fticr import orbitrap_fticr
from  interface import interf
from interface_actions import interact

'''
Program for viewing Orbitrap and FTICR datasets..
Interactions are done in interface_actions.py
The interface graphics is produced by interfOrbitUi.py.
The main functions of the interface are done in interface.py.
'''
def run(namef = None):
    
    interface = interf() # instantiates the interface
    orb_fti = orbitrap_fticr(interface) # instantiates the library for Orbitrap datasets.
    inter = interact(interface, orb_fti) # interaction of the interface with orbitrap library
    if namef :
        interface.canvas.droplist.append(namef)
        inter.drop_refresh()
    interface.run() # launches the interface. 

if __name__ == '__main__':
    
    #namef = '/home/lio/Bureau/Encours/Orbitrap/results/20130523_ubi_5_scan_res_7500.pkl'
    #run(namef)
    run()
