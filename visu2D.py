# encoding: utf-8
"""
Created by Marc Andre Delsuc & Lionel Chiron on 2011-05-19.
Copyright (c) 2011 IGBMC. All rights reserved.
###
Program for visualizing FTICR2D data. 
"""
import sys, os
from Visu.Load import LOAD #
from Visu.Saving import SAVE
from Visu.graphic_tools import GRAPHTOOLS       # class contiaining the graphic tools for doing profiles etc.. 
from Visu.display import DISPLAY            # class to handle diplay of the dataset.
from Visu.interface_actions import INTERACT     # class for to handle interaction with the interface.
from Visu.canvas import Qt4MplCanvas as QtMplCv     # class for using matplotlib in Qt canvas.
from Visu.paramzoom import PARAM_ZOOM               # class object for the zoom parameters
from Visu.zooming import ZOOMING            # class for taking care of the zooms
from Visu.move_window import MOVE_WINDOW    # class to handle the zoom windows positions
from Visu.zooming import ZOOM3D             # class for viewing peaks in 3D
from Visu.interface import INTERFACE    # class for completing the interface. 
from Visu.convert import CONVERT            # Class for conversion operations between mz and point
from Visu.single.select_tools import SELECT_TOOLS # class to handle the selected tool used in the interface.
#from util.log_all import Logger # Logger for saving stdoud and sterr in a log file. 

def debugs_activate(*args):
    '''
    Debugging class methods.
    Classes debugged are
    -interface
    -display
    -convert
    -gtools
    -zooming
    -move window
    -interact
    '''
    interf, display, convert, gtools, zoom, mwind, interact = args
    ### interface
    interf.clearlayout_debug = False
    interf.affichd_debug = False
    interf.affi_debug = False
    ### display
    display.select_best_resolution_debug = False
    display.affd_debug = False
    display.affi_debug = False
    display.affichd_debug = False
    display.change_resolution_debug = False
    display.set_displayasC_debug = False
    display.distrib_debug = False
    display.set_canvasD_debug = False
    ### convert
    convert.pass_to_pt_debug = False
    ### gtools
    gtools.pass_to_curr_mode_debug = False
    gtools.itomzcoorr_debug = False
    gtools.make_ptl_mz_debug = False
    gtools.param_zoom_right_order_debug = True
    ### zooming
    zoom.change_view_debug = False
    zoom.change_view_from_list_debug = False
    zoom.on_press_C_event_debug = False
    zoom.on_press_debug = False
    zoom.on_press_zoom_debug = False
    zoom.press_zoomready_and_out_debug = False
    zoom.press_zoom_debug = False
    zoom.changeres_listzoom_debug = False
    ### move_window
    mwind.moverectC_debug = False
    mwind.moverectD_debug = False
    mwind.moverect_debug = False
    ### interface_actions
    interact.coord_profile_y_debug = False
    interact.coord_profile_x_debug = False
    interact.coord_profile_debug = False
    interact.selectL_debug = False
    interact.afffile_debug = False
    interact.take_lineEdit_xy_format_debug = False

def main(argv = None):
    """ creates and runs """
    #####################
    ######### read arguments
    # logger_activated = False
    # if logger_activated:
    #     sys.stderr = Logger(erase = True)
    #     sys.stdout = Logger()
    configfile = None
    msh5file = None
    if not argv:
        argv = sys.argv
    try:                        # First try to read config file from arg list
        root, ext = os.path.splitext(argv[1])
        if ext == '.mscf':
            configfile = argv[1]
        elif ext == '.msh5': 
            msh5file = argv[1]
    except IndexError:          # then assume standard name
        configfile = "Visu/visu2D.mscf"
    if msh5file:
        data = LOAD(msh5file = msh5file)  
    else:
        if configfile:
            print "using %s as configuration file" %configfile
        data = LOAD(configfile = configfile)                                                                  # loads msh5 file
    save = SAVE(data)                                                                               # saves 2D, 3D, profiles.
    paramz = PARAM_ZOOM(data)                                                                       # takes the parameters for zoom. 
    interf = INTERFACE()                                                                            # instantiate the interface 
    display = DISPLAY(QtMplCv, data, interf, paramz)                                                # control the display, zoom/resolution. etc.. 
    convert = CONVERT(display, data, paramz)                                                                # conversion mz/point
    gtools = GRAPHTOOLS(paramz, display, data, save, convert)                                                # graphic tools
    stools = SELECT_TOOLS(display, paramz, interf)                                                          # orthogonally select tools 
    mwind = MOVE_WINDOW(display, interf, data, paramz, gtools, convert, stools)   # moving zoom window, drag etc.. 
    zoom = ZOOMING(display, interf, data, paramz, gtools, convert, stools, mwind)                          # zooming tools.
    
    zoom3d = ZOOM3D()                                                                       # zoom 3D on localalized area.
    interact = INTERACT(zoom, mwind, data, display, interf, paramz,\
        gtools, zoom3d, stools, convert, save)                                                      # interactions with the interface.
    debugs_activate(interf, display, convert, gtools, zoom, mwind, interact)                               # activates all the debugs.
    interact.interfGraph()                                                                          #executing Main loop. 
    interf.run()

if __name__ == '__main__':
    print "python version ", sys.version
    main()
