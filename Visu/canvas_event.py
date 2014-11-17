from util.debug_tools import*  
from Matplotlib_generictools import*
import numpy as np
import re
from zoom_plot import ZOOM_PLOT 
from Pyside_PyQt4 import*

@dec_class_pr
@decclassdebugging
class CANVAS_EVENT(object):
    '''
    Gestion of the events with the interface Canvas.
    '''
    def __init__(self, display, interf, data, paramz, gtools, convert, stools, mwind, zoom):
        self.display = display 
        self.interface = interf
        self.data = data
        self.paramz = paramz # zoom parameters
        self.gtools = gtools  # graphic tools
        self.convert = convert # conversions mz/point
        self.select_tools = stools # gestion of the selected tools.
        self.move_wind = mwind # zoom window movements after drag or zoom or anything else.
        self.zoom = zoom
        self.zplot = ZOOM_PLOT(display, interf, data, paramz, gtools)
        self.rectd = None                                                                           # graphical object for the zoom rectangle of the small window
        self.rectc = None                                                                           # graphical object for the zoom rectangle of the main window
        self.layoutC = self.interface.ui.layoutC
        self.layoutD = self.interface.ui.layoutD

    def recupxy(self, event):
        '''
        Take the event coordinates and transform from "m/z" to "point" if necessary.
        '''
        if debug(self):
            print "in zooming.recupxy"
            self.paramz.report()
        dx = event.xdata; dy = event.ydata
        if not self.data.mode_point and event.canvas != self.display.qmd.fig.canvas :               # if mode point false and not taken in D
            dd = self.display.res2dd()                                                              # retrieves fticrdata
            if dx is not None and dy is not None:
                dx = self.convert.mztoi(dd, dx, 2)                                                  # if in mode m/z translates x from point format to m/z format
                dy = self.convert.mztoi(dd, dy, 1)                                                  # if in mode m/z translates y from point format to m/z format
        return dx, dy, event.xdata, event.ydata
    
    def on_motion(self, event):
        '''
        Activate on mouse motion.
        '''
        zc = self.paramz.zoom_coord
        dx, dy, dxcurr, dycurr = self.recupxy(event)                      # retrieves dx, dy (point) and dxcurr dycurr (current mode)
        if dx is not None and  dy is not None:
            self.lastdx, self.lastdy = dx, dy
        if debug(self):
            print "self.select_tools.zoom ", self.select_tools.zoom
        if dx is not None and dy is not None and len(zc) == 2 :
            if self.zplot.rectc is None and self.select_tools.zoom :                         # if no rectangle in C window and no line tool
                self.zplot.drawrect(zc[0], zc[1], dx, dy, self.layoutC, self.layoutD)
            if event.canvas == self.display.qmd.fig.canvas :                                        # if in canvas D
                ####                                                                                # Stretch zoom in window D
                if self.zoom.stretch :                                                              # if stretch activated because on corner and press
                    self.zoom.stretchrect(dx,dy)                                                    # stretchs rectangle in canvas D
                else :                                                                              # motion canvas D 
                    self.zplot.rectd.remove()                                                       # erase the zoom rectangle in canvas d
                    self.move_wind.moverect([dx, dy], "D")                                               # move the zoom rectagle to position dx, dy                                                                                    # Draw zoom in C and D windows on motion
            #############################                                                           # window C
            if event.canvas == self.display.qmc.fig.canvas and event.button == 1 and self.select_tools.zoom:                   # if in canvas C
                self.move_wind.move_refreshC(dx, dy)
            self.interface.ui.label_16.setText(str(dxcurr)) # print in the interface x coordinate in current mode
            self.interface.ui.label_18.setText(str(dycurr))  # print in the interface y coordinate in current mode
    
    def on_press_D_event(self, event, xpress, ypress): 
        '''
        Permits to interact with the zoom in window D. 
        '''
        if debug(self):   
            print "in zooming.on_press_D_event"
        if len(self.paramz.zoom_coord) == 0 :
            self.paramz.zoom_coord = self.paramz.zoom_coord_old
        if event.inaxes:                                                                            # when inside window D
            self.paramz.mouse_motion = self.display.connect('motion_notify_event', self.on_motion, 'D')
            self.zoom.scro = self.display.connect('scroll_event', self.zoom.on_scroll, 'D')
        #######                                                                                     # detect position for stretching window
        zc = self.paramz.zoom_coord
        if len(zc) == 4 :                                                                           # recuperating corner to enlarge the zoom from global zoom window.
            self.zoom.listrch = [[zc[0], zc[1]], [zc[2], zc[3]]]                                    # list for stretch.
            for elemch in self.zoom.listrch :
                if np.sqrt((xpress-elemch[0])**2 + (ypress-elemch[1])**2) < 60 :                    # if mouse near corner at less than 60 points.
                     self.zoom.plstrch = self.zoom.listrch.index(elemch)
                     self.zoom.corner = self.display.qmd.axes.plot(\
                            elemch[0], elemch[1], self.zoom.colpt)                                  # plot a plot on the corner
                     self.display.qmd.fig.canvas.draw()                                             # draw the stretched zoom
                     #self.display.qmd.setCursor(QtGui.QCursor(QtCore.Qt.CrossCursor))               
                     self.display.setcursor('cross','D')   # makes the cursor in cross mode
                     self.zoom.stretch = True

    def on_press_C_event(self, event, xpress, ypress):
        '''
        Makes the drawing of the zoom window in the layout C.
        If a click is produced again in the window, it makes the zoom.
        '''
        if debug(self):
            print "in zooming.on_press_C_event"
        if event.inaxes  :                                                    # if left button:# when in window C
            self.paramz.zoom_coord_old = self.paramz.zoom_coord                                     # keep in memory the window
            if debug(self):
                print "in on_press_C_event keeping old coord self.paramz.zoom_coord_old ", self.paramz.zoom_coord_old
            if len(self.paramz.zoom_coord) >  2:                                                    # if there is yet one corner the new press eliminates it.
                if debug(self):
                    print "in zooming.on_press_C_event, len(self.paramz.zoom_coord) >  2 "
                self.zoom.press_zoom(xpress, ypress)
                if debug(self):
                    print "in zooming.on_press_C_event ", xpress, ypress
            if not self.paramz.zoomready :                                                          # if zoom window is not here.. add a point and connect to motion
                if debug(self):
                    print "zoom not ready "
                self.paramz.zoom_coord += [xpress, ypress]              # add a new point to the zoom window
                if debug(self):                                         
                    print "in on_press_C_event, self.paramz.zoom_coord ", self.paramz.zoom_coord
                self.paramz.mouse_motion = self.display.connect('motion_notify_event', self.on_motion)                                      # 
            self.paramz.zoomready = False
            if debug(self):
                print "self.paramz.zoomready = False "

    def on_press(self, event):                                                                      # 
        '''
        When pressed, triggers the rectangle drawing
        Waits for the mouse button's release to make the zoom. 
        Calls self.on_press_D_event and self.on_press_C_event
        '''
        if debug(self):   
            print "in zooming.on_press"
        dx, dy, dxcurr, dycurr = self.recupxy(event)                                                        # gets the x and y coords in mode point.
        if event.canvas == self.display.qmd.fig.canvas and event.button == 1 :                      # if canvas D selected and click left
           self.on_press_D_event(event, dx, dy)  
        if event.canvas == self.display.qmc.fig.canvas and event.button == 1 :                      # if canvas C selected and click left
           if debug(self):
               print "pressed in canvas C ", dx, dy
           self.on_press_C_event(event, dx, dy)
           if debug(self):
               print "in on_press "
               self.paramz.report()
    ##########                                                                                      # on release
    
    def release_refrechC(self, name_profile = None):
        '''
        When mouse is released, refreshes layout C.
        '''
        self.display.disconnect(self.paramz.mouse_motion)                             # sets free the mouse from C window
        if debug(self):   
           print "in zooming.release_refrechC"
        zc = self.paramz.zoom_coord
        if debug(self):
            print "self.paramz.zoom_coord ", self.paramz.zoom_coord
        if not self.select_tools.drag:
            self.convert.set(zc[0], zc[1], zc[2], zc[3])
        self.zoom.zoom_check_size()
        if self.select_tools.profile:                                                               # if line tool selected plot profile 
           self.gtools.plotprofile(self.zoom.profile_type, name_profile = name_profile)             # makes the profile in an other window
           self.zoom.profile_type = None                                                            # reinitializes profile type
           self.select_tools.change_to('zoom')                                                      # changes from profile to zoom.
        if self.select_tools.drag  and not self.paramz.zoomready :                                  # if zoom not ready and drag mode
           print "in release_refrechC self.paramz.zoom_coord ", self.paramz.zoom_coord
           self.move_wind.moverect([1,1], "C")
           self.zoom.change_view()     
           print "self.paramz.zoom_coord ", self.paramz.zoom_coord                                                                  # 
        if not self.select_tools.drag:
            self.paramz.zoomready = True # ready for performing a zoom by clicking inside.
        if debug(self):
            print "####### end of release_refrechC"

    def on_release_D_event(self): 
        '''
        Release D event
        '''
        if debug(self):   
           print "in zooming.on_release_D_event"
        self.stretch = False
        self.display.disconnect(self.paramz.mouse_motion, 'D')                                # sets free the mouse from D window
        if self.select_tools.zoom:
            self.register_coordinates()
            self.display.change_resolution(layout1 = self.layoutC, layout2 = self.layoutD)         
            self.plot_zooms()                                                                           # draw rectangle of new zoom

    def aff_param(self):
        '''
        Show resolution and print coordinates in lineEdit.
        '''
        self.display.aff_resolution()
        self.make_coord_manual()
    
    def make_coord_manual(self):
        '''
        Mouse coordinates are automatically written in the interface for manual interaction.
        '''
        if debug(self):   
           print "in zooming.make_coord_manual"
           print "in make_coord_manual self.paramz.zoom_coord", self.paramz.zoom_coord
        llx, lly, urx, ury = self.gtools.convert.pass_to_curr_mode(*self.paramz.zoom_coord)                 # pass from zoom.coord to current format pt or m/z                                                                                         #     llx, lly, urx, ury = conv.mztoi_all(self, d, llx, lly, urx, ury)
        if debug(self): 
            print "llx, lly, urx, ury ", llx, lly, urx, ury
        if self.data.mode_point:
            strzoom_coord = str(map(int, [llx, lly, urx, ury]))                                         # transform to intergers 
        else:
            strzoom_coord = str(map(int, [urx, ury, llx, lly]))                  
        zmcoord = re.findall('(\d*[,].*[^]]\d*)', strzoom_coord)[0]                                 # filters the string 
        if debug(self):
            print "in canvas_event zmcoord ", zmcoord
        self.interface.ui.lineEdit_2.setText(zmcoord)                                               # Edits coordinates in LineEdit.
        
    def on_release_C_event(self, dx, dy):  
        '''
        Release C event
        
        '''
        if debug(self):   
           print "in zooming.on_release_C_event"                                                    # if len(self.paramz.zoom_coord) >=  4 : #if nb of points too high                                                                                   #     self.paramz.zoom_coord = []  # reinitialize the zoom window zoom_coord
        if self.zoom.newrectcready :                                         # ready to draw new rectangle and not dragging.
           print "add new point to self.paramz.zoom_coord "
           self.paramz.zoom_coord += [dx, dy]                                                       # add a new point to zoom_coord
        if len(self.paramz.zoom_coord) == 2 :
           print "self.zoom.newrectcready ", self.zoom.newrectcready
           self.paramz.zoom_coord = []                                                              # reinitialize zoom_coord if rectangle has a too small size.
           self.display.disconnect(self.paramz.mouse_motion)
        if len(self.paramz.zoom_coord) == 4 :   
           print "four coordinates so refreshes "                                                    
           self.release_refrechC()                                                                  # Refreshes layout C
           self.aff_param()                                                                         # writes zoom coordinates and resolution                                                                 

    def on_release(self, event):
        '''
        On release keep the rectangle and makes the zoom
        on_released is used to make lines .
        self.paramz.listview, list of all the zooms. 
        '''
        if debug(self):   
           print "in zooming.on_release"
        dx, dy, dxcurr, dycurr = self.recupxy(event)                                                 # retrieves x and y when event is detected.
        if event.canvas == self.display.qmd.fig.canvas :
            if event.inaxes:                                                                           # when in window D
                self.on_release_D_event()
        if event.canvas == self.display.qmc.fig.canvas :                                            # if event in window C
            if event.inaxes:                                                                           # when in window C
                self.on_release_C_event(dx, dy)
            else:
                self.on_release_C_event(self.lastdx, self.lastdy)   # release outside axes, so use last dx and dy

    def detect_corner(self,event):
        '''
        detect position for stretching window and put the good mouse's shape
        '''
        if len(self.paramz.zoom_coord) == 4 :                                                       # recuperating corner to enlarge the zoom from global zoom window.
           dx, dy, dxcurr, dycurr = self.recupxy(event)                                                         # get the x and y pixel coords
           zc = self.paramz.zoom_coord
           self.listrch = [[zc[0], zc[1]], [zc[2], zc[3]]]                                          # list for stretch.
           if dx :
               for elemch in self.listrch :
                   if np.sqrt((dx-elemch[0])**2 + (dy-elemch[1])**2) < 60 :                     # if mouse near corner at less than 60 points.
                       self.plstrch = self.listrch.index(elemch)                                    # index of coordinates found in listrech
                       self.display.setcursor('cross','D')   # make the cursor in cross mode if on corner
                       break                                                                        # to avoid to check second element in case of low left corner
                   else :
                       self.display.qmd.setCursor(QtGui.QCursor(QtCore.Qt.OpenHandCursor))          
                       self.display.setcursor('hand','D')   # make the cursor a hand if out of corner

    def interact_with_canvasC(self):
        '''
        Connects event to canvas C. Press, release, corner detection.
        '''
        if debug(self):
           print "in zooming.interact_with_canvasC"
        self.display.connect('button_press_event', self.on_press)                # connecting the press event 
        self.display.connect('button_release_event', self.on_release)            # connecting the release event
        self.display.connect('motion_notify_event', self.detect_corner) 
        #self.fig.canvas.mpl_connect('scroll_event', zo.on_scroll)