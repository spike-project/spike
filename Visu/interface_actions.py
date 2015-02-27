from util.debug_tools import*  
import sys, os, re, time
from Visu.init.fticrvisuUi import Ui_MainWindow as Ui
from Visu.init.fticrvisupopUi import Ui_Form as Pop 
from Pyside_PyQt4 import*
from canvas_event import CANVAS_EVENT
from functools import partial

@dec_class_pr
@decclassdebugging
class INTERACT(object):
    def __init__(self, zoom, mwind, data, display, interf, paramz, gtools, zoom3d, stools, convert, save):
        print "in INTERACT"
        self.zoom = zoom
        self.data = data            # data loaded
        self.display = display
        self.interface = interf     # interface.. 
        self.paramz = paramz        # parameters for zoom
        self.gtools = gtools        # graphic tools
        self.select_tools = stools  # tools selection 
        self.zoom3d = zoom3d        # 3D zoom
        self.convert = convert      # conversions m/z point
        self.save = save 
        ########
        self.canv_event = CANVAS_EVENT(display, interf, data, paramz, gtools, convert, stools, mwind, zoom)
        self.list_res()             # list of resolutions, made from self.data
        self.lowmass_shift = 1e-7   # very small shift for profiles.
        self.pt_shift = 1e-7        # very small shift for doing the profiles
        #self.list_res()  
        self.select_tools.change_to('zoom') # position the tool to zoom tool.
        self.interface.ui.pushButton_4.setText("go to m/z")

    def list_res(self):
        '''
        resolutions list made from self.data
        '''
        self.LISTR = ['resol'+ str(i) for i in xrange(self.data.NBRES, 0, -1)]
        self.LISTRMAX = self.LISTR[0]                                                               # maximal resolution
        self.LISTD = self.data.d                                                                    # list of the different resolutions
        self.LISTDNAME = ['d' + str(i) for i in xrange(self.data.NBRES, 0, -1)]                     # d1 biggest resolution..
    
    def _name_number(self, name, nb):
        '''
        Used in self.button() and self.lineEdit() for making the names used in those methods.
        '''
        if nb == 0:
            name = name
        else:
            name = name + '_' + str(nb)
        return name
    
    def button(self, nb, action, name_icon = None, icon_size = None):
        '''
        General definition for the buttons
        nb : number for the button
        action : method associated to the action.
        name_icon : name of the icon in the directory Visu/iconsUi/, must be in png format.
        icon_size : size of the icon for fitting to the button size.
        '''
        name = self._name_number('pushButton', nb)
        button = getattr(self.interface.ui, name)
        Qobj.connect(button, SIGNAL("clicked()"), action)                   # attach an acton to the button.
        if name_icon:           
            icon = QtGui.QIcon("Visu/iconsUi/" + name_icon + ".png")      # use an icon picture for the button.                                
            button.setIcon(icon) 
        if icon_size:                                                       
            button.setIconSize(QtCore.QSize(icon_size[0], icon_size[1]))   # set the size of the icon
    
    def lineEdit(self, nb, action):
        '''
        General lineEdit
        nb : number of the lineEdit
        action : associated action.
        '''
        name = self._name_number('lineEdit', nb)
        button = getattr(self.interface.ui, name)
        button.returnPressed.connect(action)
    
    def lineEdit_and_button(self, name_lineEdit, name_button, action):
        '''
        Defines a lineEdit and an associated button.
        Uses self.lineEdit() and self.button()
        '''
        self.lineEdit(name_lineEdit, action)
        self.button(name_button, action)            
    
    def interfGraph(self):                                                                          
        '''
        Defines the buttons, lineEdits and actions associated.
        '''
        ### Navigation
        self.button(0, self.backzoo, 'arrowl')                  # go backward
        self.button(2, self.forwzoo, 'arrowr')                      # go forward
        self.button(3, self.backhome, 'homeee', icon_size = (80, 100))          # home button 
        self.button(15, partial(self.select_tools.change_to,'zoom'), 'zoom_icon', icon_size = (30,30))   # zoom button         
        ### 2D view savings.
        self.button(5, self.savefigurepng)              # save png
        self.button(6, self.savefigurepdf)                  # save pdf
        ### mz/point and scale
        self.button(4, self.swap_pt_mz)                     # point or m/z        
        self.lineEdit_and_button(3, 8, self.manual_scale)   # lineEdit 3 and button 8, changes scale
        ###
        self.lineEdit_and_button(2, 13, self.manual_zoom)   # lineEdit 2 and button 13, makes manual zoom
        self.lineEdit_and_button(4, 9, self.manual_profile)  # lineEdit 4 and button 9, makes manual profile
        self.button(7, self.zoom3D)   # zoom3D
        #self.button(10, self.select_curs, 'cursor_arrow')   # cursor button
        ### Tools
        self.button(10, partial(self.select_tools.change_to,'cursor'), 'cursor_arrow')   # cursor button
        self.button(11, self.select_drag, 'cursor_hand')   # drag button
        self.button(12, self.gtools.createpeaks)            # make peakpicking                             
        if self.data.param.multresfile :
            self.afffile()                                                                          # visualization
            
    def permute_states (self, test, etat0, etat1):                                                  # pour les boutons dont a fonction est un va et vient. 
        '''
        function for swapping between two states..
        have to define first value of state before.
        '''
        if debug(self): print "etat is : ", test
        if test == etat0:  test = etat1
        else : test = etat0
        return test
    
    def drag_connect(self):
        '''
        Drag the main image.
        '''
        self.paramz.mouse_motion = self.display.connect(\
                'motion_notify_event', self.canv_event.on_motion)
        #self.display.disconnect(self.paramz.mouse_motion)                               # sets the mouse free from C window

    def select_drag(self):
        '''
        Select the drag function (hand)
        '''
        if debug(self):   
            print "in interface_actions.select_drag" 
            print self.select_tools.drag
        self.select_tools.change_to('drag')                                                         # passses to drag mode
        self.drag_connect()
        if debug(self):   
            self.select_tools.report()

    def select_curs(self):
        '''
        Select the arrow cursor
        '''
        if debug(self):   
            print "in interface_actions.select_curs"
        self.select_tools.change_to('curs')
        if self.select_tools.curs :
            self.paramz.mouse_motion = self.display.connect(\
                    'motion_notify_event', self.canv_event.on_motion) 
            self.display.setcursor('cross')   # make the cursor in cross mode
            # make the cursor a hand if out of corner
        else:
            self.display.disconnect(self.paramz.mouse_motion)                            # frees the mouse from C window 
            self.display.setcursor('cross')   # make the cursor in cross mode

    def select_manual_profile(self):
        '''
        Select the function manual_profile
        '''
        if debug(self):   
            print "in interface_actions.select_manual_profile"
        self.select_tools.change_to('manual_profile') 
        if self.select_tools.manual_profile :
            self.paramz.mouse_motion = self.display.connect(\
                    'motion_notify_event', self.canv_event.on_motion)  
            self.display.qmc.setCursor(QtGui.QCursor(QtGui.QPixmap(\
                    'Visu/iconsUi/pencil-iconsm.png'))) 
        else:
            self.display.disconnect(self.paramz.mouse_motion)                            # frees the mouse from C window
            self.display.setcursor('cross')   # make the cursor in cross mode 
        if debug(self):   
            self.select_tools.report()

    def whichres(self):                                                                             # function to know the current resolution
        for delem in self.LISTD :
            if self.display.currentd == delem :
                print self.LISTDNAME[self.LISTD.index(delem)]

    def scale_control(self):
        '''
        Show scale in the interface.
        '''
        if debug(self):   
            print "in interface_actions.scale_control"
        self.interface.ui.label.setText(str(int(self.scale)))                                       # prints value of the scale
        if len(self.paramz.zoom_coord) == 0 :
             self.paramz.zoom_coord = self.paramz.zoom_coord_old
        if debug(self):
            print "len(self.paramz.zoom_coord)", len(self.paramz.zoom_coord)
        self.display.affd(self.display.currentd, self.data.resmin, self.interface.ui.layoutC)                       # changes scale in real time only in C window
        
    def afffile(self):
        '''
        Show dataset in C window and addresses of used files
        vis.resmin : resolution for window D
        '''
        if debug(self):   
            print "##                                                                               # in interface_actions.afffile"
        sizefticr = os.path.getsize(self.data.param.multresfile)                                    # Takes the size of the multiresolution file
        self.interface.ui.label_8.setText(str(round(sizefticr/1e6)) + " MB")                        # show the size of the data in the interface 
        self.interface.ui.label_6.setText(self.data.param.multresfile)                              # Path to the multiresolution file
        if debug(self):
            print "in interface_actions.afffile makes self.display.affd "
        self.display.affd(self.data.d[len(self.data.d)-1], self.data.resmin,
                      self.interface.ui.layoutC, self.interface.ui.layoutD)                                                 #initialisation with resolution 1
        if debug(self):
            print "in interface_actions.afffile makes self.canv_event.interact_with_canvasC "
        self.canv_event.interact_with_canvasC()
        absmx = self.data.d[len(self.data.d)-1].absmax
        for dd in self.data.d:
            dd.absmax = absmx
        self.data.resolu = 'resol1'                                                                 # initial resolution
        self.display.aff_resolution()                                                                     # update resolution viewing in the interface
        
    def swap_from_proint(self):
        '''
        Passes from point to m/z
        '''
        if debug(self):   
            print "in interface_actions.swap_from_proint"
        self.data.mode_point = False
        self.interface.ui.pushButton_4.setText("go to pt")
        for dd in self.data.d:
            if debug(self):   
                print "dd.units = m/z"
            dd.units = "m/z"
    
    def swap_from_mz(self):
        '''
        Passes from m/z to point
        '''
        if debug(self): print "swap_from_mz "
        self.data.mode_point = True
        self.interface.ui.pushButton_4.setText("go to m/z")
        for dd in self.data.d:  # change the mode in all the resolutions
            if debug(self): print "dd.units = pt"
            dd.units = "points"

    def swap_pt_mz(self):                                                                         
        '''
        Function for passing from "point mode" to "m/z mode" and inversely.
        Change only in C window
        '''
        if debug(self):
            print " in swap_pt_mz, self.paramz.zoom_coord ", self.paramz.zoom_coord
            print "self.data.mode_point ", self.data.mode_point
        if self.data.mode_point :                                                                   # if mode point is true
            self.swap_from_proint()
        elif not self.data.mode_point :
            self.swap_from_mz()
        if debug(self):
            print "after swap_from_mz,  self.paramz.zoom_coord ", self.paramz.zoom_coord
        self.display.affd(self.display.currentd, self.data.resmin,\
            self.zoom.layoutC, self.zoom.layoutD, zoom = not(self.paramz.zoomready) ) #not(self.paramz.zoomready)
        self.select_tools.change_to('zoom')
        self.canv_event.interact_with_canvasC()
        self.canv_event.aff_param()
      
    def coord_profile_y(self, yval):
        '''
        Makes coordinates profile for y profile 
        Called by coord_profile
        Returns the extreme coordinates in m/z or point format
        '''
        dd = self.display.currentd
        if debug(self):
            print "in interface_actions.coord_profile_y"
            print "in interface_actions.coord_profile_y yval is ", yval
        if self.data.mode_point :
            '''
            point mode
            '''
            llx, lly, urx, ury = 1., yval + self.pt_shift, dd.size2, yval
        else:
            '''
            m/z mode
            '''
            if debug(self):
                print "in interface_actions.coord_profile_y  self.data.mode_point = True"
                print "dd.axes(2).highmass ", dd.axes(2).highmass
            llx, lly, urx, ury = dd.axes(2).highmass, yval, dd.axes(2).lowmass + self.lowmass_shift, yval
        if debug(self):
            print 'in interface_actions.coord_profile_y coordinates for profile are ', llx, lly, urx, ury
        return llx, lly, urx, ury
    
    def coord_profile_x(self, xval):
        '''
        Makes coordinates profile for x profile.
        Called by coord_profile.
        Returns the extreme coordinates in m/z or point format.
        '''
        dd = self.display.currentd
        if debug(self):
            print "in interface_actions.coord_profile_x"
            print "in interface_actions.coord_profile_x xval is ", xval
        if self.data.mode_point :
            '''
            point mode
            '''
            llx, lly, urx, ury = xval + self.pt_shift, 1, xval, dd.size1
        else:
            '''
            m/z mode
            '''
            if debug(self):
                print "in interface_actions.coord_profile_x  self.data.mode_point = True"
                print "dd.axes(1).highmass ", dd.axes(1).highmass
            llx, lly, urx, ury = xval , dd.axes(1).highmass , xval, dd.axes(1).lowmass + self.lowmass_shift
        if debug(self):
            print "in interface_actions.coord_profile_x llx, lly, urx, ury is ", llx, lly, urx, ury
        return llx, lly, urx, ury
    
    def prepare_coord_profile(self, values, ct = None):
        '''
        Makes coordinates profile
        Called by take_lineEdit_xy_format
        '''
        if debug(self):
            print "###########                                                                      # in interface_actions.coord_profile "
        val = float(values.split(ct)[1])
        if ct == 'x':
            self.zoom.profile_type = 'x'
            if debug(self):
                print "in interface_actions.coord_profile case x "
            llx, lly, urx, ury = self.coord_profile_x(val)                                          # Completes x coordinates
        elif ct == 'y':
            self.zoom.profile_type = 'y'
            if debug(self):
                print "in interface_actions.coord_profile case y "
            llx, lly, urx, ury = self.coord_profile_y(val)                                          # Completes y coordinates
        if debug(self):
            print 'in interface_actions.coord_profile coordinates for profile are ', llx, lly, urx, ury
        return llx, lly, urx, ury
           
    def take_lineEdit_xy_format(self, lineEdit_value):
        '''
        Takes the profile with format "y200" for horizontal line y=200
        ct is the coordinates type
        '''
        if debug(self):
            print "in interface_actions.take_lineEdit_xy_format values is ", lineEdit_value
            print "will check kind of coordinates "
            print "type(lineEdit_value) ", type(lineEdit_value)
            print "lineEdit_value ", lineEdit_value
        lcoord_type = ['x','y']
        for ct in lcoord_type:                                                                      # search if the lines begins with x or y
            if debug(self):
                print "lineEdit_value ", lineEdit_value
                print "ct ", ct
            if lineEdit_value.find(ct) == 0:                                                        # if x or y is found
                if debug(self):
                    print "coord type is ", ct
                self.name_profile = lineEdit_value.replace('.','-')
                llx, lly, urx, ury =  self.prepare_coord_profile(lineEdit_value, ct)                # returns profile coordinates from lineEdit value.
                if debug(self):
                    print "name profile is ", self.name_profile
                    print "extracted values for profile are ", llx, lly, urx, ury
        return llx, lly, urx, ury, ct

    def take_lineEdit(self, ledit):
        '''
        Takes values from lineEdit for zoom and profile. 
        Passes the coordinates in "point mode".
        Called by manual_profile.
        '''
        ct = None
        try:
            if debug(self):
                print "using PyQt4 lineEdit"
            values = str(ledit.text())              # coordinates from QlineEdit PyQt4
        except:
            if debug(self):
                print "using PySide lineEdit"
            values = ledit.text()                                   # coordinates from QlineEdit PySide
        if debug(self):
            print "values is ", values
        print "type(values) ", type(values)
        try:
            if self.data.mode_point:
                llx, lly, urx, ury = map(float, values.split(','))                                  # from string to float
            else:
                urx, ury, llx, lly = map(float, values.split(','))                                  # from string to float
            self.zoom.profile_type = 'diag'
            if debug(self):
                print "in take_lineEdit profile coordinates are ", llx, lly, urx, ury
            #self.name_profile = 'diagonal'
            self.name_profile = 'diag_'+ values.replace(',','-')
        except:
          try:
            print "values = ", values
            print "x/y profile"
            llx, lly, urx, ury, ct = self.take_lineEdit_xy_format(values)                           # takes the profiles in format y345 etc.. 
          except:
            print "format incorrect"
            print "values ", values                                                                  # passing to point
        if llx ==  urx:                                                                             # correction for avoiding division by 0 in profile
            if debug(self): print "adding 1e-7 to urx"                                                                                   
            urx += 1e-7
        elif lly == ury:                                                                            # correction for avoiding division by 0 in profile
            if debug(self): print "adding 1e-7 to ury"
            ury += 1e-7
        return llx, lly, urx, ury

    def manual_scale(self):
        '''
        Manual scaling.
        '''
        self.paramz.scale = float(self.interface.ui.lineEdit_3.text())
        self.select_tools.change_to('scale')
        self.zoom.change_view(change_layoutD = True)
        self.canv_event.release_refrechC()                                                          # Makes the zoom
        self.zoom.change_zoom()                                                                     # makes zoom and saves the zoom, resolution and scale
        if debug(self):
            print "#### changing to zoom mode   "
        self.select_tools.change_to('zoom') 
        if debug(self):
            print "self.paramz.zoomready ", self.paramz.zoomready
        
    def manual_profile(self):
        '''
        Manual profile
        '''
        if debug(self):
            print "#### in manual_profile "
        self.select_tools.change_to('profile')                                                      # set mode to profile.
        self.select_tools.report()
        llx, lly, urx, ury = self.take_lineEdit(self.interface.ui.lineEdit_4)                       # retrieves from interface values for profile.
        if debug(self):
            print "pass to current mode"
            print "drawline"
            print "llx, lly, urx, ury ", llx, lly, urx, ury
        self.gtools.drawline(llx, lly, urx, ury, layout1 = self.interface.ui.layoutC)               # draws the line in the 2D on which is made the profile.
        llx, lly, urx, ury = self.convert.maxres(llx, lly, urx, ury)
        self.gtools.do_profile_coord(llx, lly, urx, ury)                                            # Makes the coordinates for plotting the profile in a new window
        self.canv_event.release_refrechC(self.name_profile)                                         # makes the profile.

    def manual_zoom(self):
        '''
        Manual zoom
        Coordinates are entered manually with format llx, lly, urx, ury.
        '''
        llx, lly, urx, ury = self.take_lineEdit(self.interface.ui.lineEdit_2)
        if self.data.mode_point:
            self.convert.set(llx, lly, urx, ury)                                                       # saves coordinates in self.paramz.zoom_coord
        else:
            llx, lly, urx, ury = self.convert.pass_to_pt(llx, lly, urx, ury)
            self.convert.set(llx, lly, urx, ury)
        if debug(self): print "in manual_zoom, self.paramz.zoom_coord ", self.paramz.zoom_coord
        self.zoom.change_view(change_layoutD = True)                                                # Refreshes resolution
        self.canv_event.release_refrechC()                                                          # prepares the zoom
        self.zoom.change_zoom()                                                                     # makes zoom and put the zoom and resolution in a list
    
    def savefigure(self, kind):
        '''
        Function to save figure from window C.
        '''
        path_save  = self.save.prep_path_save(self.data.namefile + str(self.interface.ui.lineEdit.text()) + '.' + kind)
        self.display.qmc.fig.savefig(path_save)                                                     # saves with matplotlib
        self.interface.ui.label_9.setText("picture saved")

    def savefigurepng(self):
        '''
        Save the main view in png
        '''
        self.savefigure('png')
        
    def savefigurepdf(self):
        '''
        Save the main view in pdf
        '''
        self.savefigure('pdf')

    def backhome(self):                                                                             # retrieving last zoom and resolution
        '''
        going to the original view
        '''
        self.paramz.listview_index = 0
        self.zoom.change_view_from_list()                                                           # changing the zoom with history
        self.select_tools.change_to('zoom')
        
    def backzoo(self):                                                                              # retrieving last zoom and resolution
        '''
        going back in the zooms
        '''
        if len(self.paramz.listview) > 1 and self.paramz.listview_index > 0:                        # if there is more than one element in list of coordinates.
            self.paramz.listview_index += -1
        self.zoom.change_view_from_list()                                                           # changing the zoom with history  

    def forwzoo(self):                                                                              # retrieving last zoom and resolution
        '''
        going forward in the zooms
        '''
        if self.paramz.listview_index < len(self.paramz.listview) - 1 :                                # if there is more than one element in list of coordinates.
            self.paramz.listview_index += 1
        self.zoom.change_view_from_list()                                                              # changing the zoom with history
    
    def zoom3D(self):
        '''
        From zoom coordinates, calculates the zoom area. 
        If the area is small enough, makes the 3D reprensentation in m/z coordinates.
        '''                                                                                   # in interface_actions, 3D zoom"
        f2min, f1min, f2max, f1max = map(int, self.paramz.zoom_coord)
        if debug(self): print "f2min, f2max, f1min, f1max ",f2min, f2max, f1min, f1max
        area = abs(f1min - f1max)*abs(f2min - f2max)                                                # calculates the zoom area.
        if debug(self): print area
        if area < 1e9: #self.paramz.area3Dmax                                                                   # test if area is small enough
            if self.data.mode_point:
                print "point"
            else:
                print "in m/z mode"
            d = self.display.data.d[len(self.display.data.d)-1]                                       # lowest resolution
            #d = self.display.data.currentd 
            self.zoom3d.trunc = 1
            self.zoom3d.plotregion3d(d, f1min, f1max, f2min, f2max, visible = True)                 # plot the 3D view
            plt.show()
        else :
            message = " the zoom area is too large"
            print message
            self.display.affd(self.display.currentd, self.data.resmin,\
             self.interface.ui.layoutC, zoom = False, message = message)
        self.select_tools.change_to('zoom')
