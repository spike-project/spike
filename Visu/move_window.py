from util.debug_tools import*  
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Ellipse, Rectangle
import matplotlib.lines as mlines
import numpy as np
from mpl_toolkits.mplot3d import axes3d
from Pyside_PyQt4 import*
from zoom_plot import ZOOM_PLOT 

@dec_class_pr
@decclassdebugging
class MOVE_WINDOW():                                                                                    # class regrouping methods about drag function
    '''
    Drag the zoomed area or just makes the zooming area.
    The coordinates are in "point" format at the root of the treatment
    so as to simplify all the procedures. 
    '''
    def __init__(self, display, interf, data, paramz, gtools, convert, stools):
        print "in MOVE_WINDOW"
        self.display = display
        self.interface = interf
        self.layoutC = self.interface.ui.layoutC
        self.layoutD = self.interface.ui.layoutD
        self.data = data   # data loaded
        self.paramz = paramz
        self.gtools = gtools # graphic tools
        self.select_tools = stools # selected tools
        self.convert = convert
        self.zplot = ZOOM_PLOT(display, interf, data, paramz, gtools)
        self.stretch = False                                                                        # switch to allow stretching of windows by taking the corners
        self.selline = False 
        self.vecshift = []
        self.drag = False
        self.listrch = None
        self.scroll = None
     
    def move_dragC(self, inzoomdrag, dx, dy):
        '''
        Moving the rectangle in C
        Dragging the window
        '''
        try:
            self.zplot.rectc.remove()                                                               # erase the zoom rectangle in canvas c
            self.zplot.rectd.remove()                                                               # erase the zoom rectangle in canvas d
            #self.gzoo.remove()                                                                     # erase grey blue area
        except :
            pass
        if inzoomdrag :                                                                             # move the zoom rectangle
            self.moverect([dx, dy],"Czoom")                                                         # dragging the zoom window
        elif  self.paramz.zoomready == False:
            self.moverect([dx, dy], "C")                                                            # dragging all the picture

    def move_refreshC_line(self, dx, dy):
        '''
        Moves  the profile line with mouse.
        '''
        if debug(self):
            print 'in move_refreshC_line '
        zc = self.paramz.zoom_coord
        self.display.qmc.setCursor(QtGui.QCursor(QtGui.QPixmap(\
                'Visu/iconsUi/pencil-iconsm.jpg')))
        if len(zc) == 2 :                                                                           # draw intermediate line
            try :
                self.linec.remove()                                                                 # erase the line in canvas C
            except :
                pass 
            self.gtools.drawline(zc[0], zc[1], dx, dy, self.layoutC, self.layoutD)                  # draw line in windows C

    def move_refreshC_rect(self, dx, dy):
        '''
        Moves the rectangle in C window.
        '''
        if self.zplot.rectc is not None :
            try:
                self.zplot.rectc.remove()                                                               # erases the zoom rectangle in canvas c
            except:
                print "no self.zplot.rectc"
        if self.zplot.rectd is not None :
            try:
                self.zplot.rectd.remove()                                                               # erases the zoom rectangle in canvas d
            except:
                print "no self.zplot.rectd"
        zc = self.paramz.zoom_coord
        self.zplot.drawrect(zc[0], zc[1], dx, dy, self.layoutC, self.layoutD)                       # draw zoom rectangles in windows C and D

    def move_refreshC(self, dx, dy):
        '''
        Moves the C zoom window or draw a line for profile. 
        if drag mode, drags the whole window else changes its size.
        '''
        if debug(self):
             print "in move_refreshC "
        zc = self.paramz.zoom_coord
        try :
            llx , lly = self.display.distrib(min, zc)
            urx, ury = self.display.distrib(max, zc)
            inzoomdrag = dx > llx and dx < urx and dy > lly and dy < ury and self.paramz.zoomready
        except :
            pass
        if  self.drag :                                                                             # right mouse button  "hand cursor"
            self.move_dragC(inzoomdrag, dx, dy)
        ########                                                                                    # Cursor mode
        if not self.drag :                                                                          # left mouse button   "arrow cursor"
            if self.selline :                                                                       # if line selection # and not and not vis.mode_point
                self.move_refreshC_line(dx, dy)
            else : 
                if debug(self):
                     print "will make self.move_refreshC_rect "                                                                                 # draw zoom window in C and D
                self.move_refreshC_rect(dx, dy)
                    
    def move_zoom_win(self, poswin, llx , lly, urx, ury):
        '''
        Move zoom window C with mouse
        '''
        pwx = float(poswin[0]); pwy = float(poswin[1])
        winx = (urx-llx)/2; winy = (ury-lly)/2
        #print "winx, winy",winx, winy
        newllx, newlly, newurx, newury = pwx - winx, pwy - winy, pwx + winx, pwy + winy
        self.paramz.zoom_coord = [newllx, newlly, newurx, newury ]#
        self.zplot.drawrect(newllx, newlly, newurx,newury , layout1 = self.layoutC)                 # draw zoom rectangle in C
        self.zplot.drawrect(newllx, newlly, newurx, newury , layout2 = self.layoutD)                # draw zoom rectangle in D

    # def moverectC(self, poswin, llx , lly, urx, ury):
    #     '''
    #     Drags the zoom area with mouse
    #     '''
    #     if debug(self):
    #         print "in moverectC "
    #     pwx = float(poswin[0]); pwy = float(poswin[1])
    #     self.vecshift += [pwx, pwy]
    #     winx = (urx - llx)/2; winy = (ury - lly)/2
    #     centx = (urx + llx)/2
    #     centy = (ury + lly)/2
    #     vs = self.vecshift
    #     if len(vs) == 4 :
    #         depx = -10*(vs[2] - vs[0])
    #         depy = -10*(vs[3] - vs[1])
    #         newllx, newurx = centx - winx + depx, centx + winx + depx
    #         newlly, newury = centy - winy + depy, centy + winy + depy
    #         self.paramz.zoom_coord = [newllx, newlly, newurx, newury ] #
    #         self.zplot.drawrect(newllx, newlly, newurx, newury, layout1 = self.layoutC)             # draws zoom rectangle in C
    #         self.zplot.drawrect(newllx, newlly, newurx, newury, layout2 = self.layoutD)             # draws zoom rectangle in D
    #         self.vecshift = []
    #         #self.change_view()
    
    def moverectC(self):
        '''
        Drags the zoom area with mouse
        '''
        if debug(self):
            print "in moverectC "
        llx, lly, urx, ury = self.paramz.zoom_coord_old
        # self.gtools.param_zoom_right_order()
        
        print "before shift llx, lly, urx, ury ", llx, lly, urx, ury
        self.vecshift = list(self.paramz.zoom_diag_vector())
        print "self.vecshift = ", self.vecshift
        newllx, newlly, newurx, newury  = llx + self.vecshift[0], lly + self.vecshift[1], urx + self.vecshift[0], ury + self.vecshift[1]
        print "after shift newllx, newlly, newurx, newury ", newllx, newlly, newurx, newury
        self.paramz.zoom_coord = [newllx, newlly, newurx, newury] #
        self.zplot.drawrect(newllx, newlly, newurx, newury, layout1 = self.layoutC)             # draws zoom rectangle in C
        self.zplot.drawrect(newllx, newlly, newurx, newury, layout2 = self.layoutD)             # draws zoom rectangle in D
        self.vecshift = []

    def moverectD(self, poswin, llx , lly, urx, ury):
        '''
        Moves zoom window D
        Takes care of the ration due to the changing resolution.
        '''
        if debug(self):
            print "in moverectD "
        rapp2, rapp1 = (self.data.resmin.size2/float(self.display.currentd.size2),
                        self.data.resmin.size1/float(self.display.currentd.size1))                  # 
        winx = (urx - llx)/2; winy = (ury-lly)/2
        pwx = float(poswin[0])/rapp2; pwy = float(poswin[1])/rapp1
        newllx, newlly, newurx, newury = pwx - winx, pwy - winy, pwx + winx, pwy + winy
        self.paramz.zoom_coord = [newllx, newlly, newurx, newury ]#
        self.zplot.drawrect(newllx, newlly, newurx, newury, layout2 = self.layoutD)                 # draws zoom rectangle in D
                                                                         # Calculates new position in window C
    def moverect(self, poswin, winselect):
        '''
        Function to move zoom area in window "C" and window "D"
        winselect is the window in which we want to move the zoom rectangle it is "C " or "D"
        '''
        if debug(self):
            print "in moverect "
        llx, lly, urx, ury = self.gtools.param_zoom_right_order()
        if winselect == "D" :
             self.moverectD(poswin, llx , lly, urx, ury)
        if winselect == "C" :
            self.moverectC()
        if winselect == "Czoom" :
            self.move_zoom_win(poswin, llx , lly, urx, ury) # 