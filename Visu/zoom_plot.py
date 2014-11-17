from util.debug_tools import*  
from Matplotlib_generictools import*
import numpy as np
from Pyside_PyQt4 import*

@dec_class_pr
@decclassdebugging
class ZOOM_PLOT(object):                                                                            # class regrouping methods about zoom
    '''
    Makes the zoom rectangles for C and D window and the greyarea.
    '''
    def __init__(self, display, interf, data, paramz, gtools):
        self.display = display                              # handle the displays in the canvas.
        self.interface = interf                             # object for interactinf with the interface.
        self.data = data                                    # data loaded
        self.paramz = paramz
        self.gtools = gtools
        self.rectd = None                                                                           # graphical object for the zoom rectangle of the small window
        self.rectc = None                                                                           # graphical object for the zoom rectangle of the main window
        self.drawgreyzoom = True

    def greyzoom(self, llx, lly, urx, ury):
        '''
        Grey area around the zoom
        '''
        if debug(self):
            print "in zooming.greyzoom"
        patches = []
        hei = self.display.currentd.buffer.shape[0]
        wid = self.display.currentd.buffer.shape[1]
        x_edge = 0
        y_edge = 0
        if not self.data.mode_point:
            dax = self.display.currentd.axes(2)
            day = self.display.currentd.axes(1)
            wid = dax.lowmass                                                                       # if in mode m/z translates x from point format to m/z format
            hei = day.lowmass                                                                       # if in mode m/z translates y from point format to m/z format
            x_edge = dax.highmass
            y_edge = dax.highmass
        verts = np.array([
              (llx, lly),
            (llx, ury),
            (urx, ury),
            (urx, lly),
            (wid, lly),
            (wid, hei),
            (x_edge, hei),
            (x_edge, lly),
            (x_edge, y_edge),
            (wid, y_edge),
            (wid, lly),
            (llx, lly),
            ])
        codes = [Path.MOVETO,
                 Path.LINETO, Path.LINETO, Path.LINETO, Path.MOVETO,
                 Path.LINETO, Path.LINETO, Path.LINETO,
                 Path.LINETO, Path.LINETO, Path.LINETO, Path.MOVETO]
        pathgrey = Path(verts, codes)
        pathgreyzoo = PathPatch(pathgrey, alpha = 0.1)
        self.greyzoo = self.display.qmc.axes.add_patch(pathgreyzoo)
        self.display.qmc.fig.canvas.draw()                                                          # draw rectangle in C

    def _draw_grey_zone(self, llx, lly, urx, ury):
        '''
        Draws a grey area around the zoom. 
        '''
        if debug(self):
            print "in zooming._draw_grey_zone"
        try :
            self.greyzoo.remove()                       # remove the greyzoom first
        except :
            if debug(self):
                print "no grey zoom yet"
        llx, lly, urx, ury = self.gtools.right_order_coord(llx, lly, urx, ury)
        self.greyzoom(llx, lly, urx, ury)                                                           # draw a grey zone around the zoom
    
    def _makes_poly(self, llx, lly, urx, ury, rapp = None):
        '''
        Makes codes and the four points pt0, pt1, pt2, pt3, pt4
        for drawing the zoom rectangle.
        Called in drawrectC and drawrectD.
        '''
        codes = []
        margx = 0; margy = 0  
        if rapp:
            rapp1, rapp2 = rapp             # ratios for drawing rectangle in window D. 
        else:
            rapp1, rapp2 = 1, 1             # no change.                                                                 
        codes +=  [Path.MOVETO] + [Path.LINETO]*3 + [Path.CLOSEPOLY]
        pt0 = (int(llx*rapp2) + margx, int(lly*rapp1) + margy)
        pt1 = (int(llx*rapp2) + margx, int(ury*rapp1) - margy)
        pt2 = (int(urx*rapp2) - margx, int(ury*rapp1) - margy)
        pt3 = (int(urx*rapp2) - margx, int(lly*rapp1) + margy)
        pt4 = (0, 0)
        return codes, pt0, pt1, pt2, pt3, pt4
    
    def drawrectC(self, llx, lly, urx, ury, layout1 = None):
        '''
        Draws the zoom rectangle in the C windows
        '''
        if debug(self):
            print "in zooming.drawrectC"
        verticesC = []
        llx, lly, urx, ury = self.gtools.convert.pass_to_curr_mode(llx, lly, urx, ury)                      # if mode m/z restablish mz coordinates
        codes, pt0, pt1, pt2, pt3, pt4 = self._makes_poly(llx, lly, urx, ury)                                                                                        # make vertices
        verticesC +=  [pt0, pt1, pt2, pt3, pt4]                                                     # for window C
        verticesC = np.array(verticesC, float)
        pathC = Path(verticesC, codes)
        self.gtools.corner_color()                                                                  # change color for corner point red for m/z , green for point format. 
        pathrectc = PathPatch(pathC, facecolor = 'None', edgecolor = self.gtools.colzoo)
        if layout1:                                                                                 # if layout 1 is asked to be changed
            self.rectc = self.display.qmc.axes.add_patch(pathrectc)
            self.display.qmc.fig.canvas.draw()                                                      # draws rectangle in C
        if debug(self):
            print "coordinates for greyzoom are ", llx, lly, urx, ury
        if abs(llx - urx) > 5 and abs(lly - ury) > 5:                                                       # makes the grey zone only if rectangle large enough
            if self.drawgreyzoom:
                self._draw_grey_zone(llx, lly, urx, ury)                                                     # draw grey zone around the zoom.

    def drawrectD(self, llx, lly, urx, ury, layout2 = None):
        '''
        Draws the zoom rectangle in the D windows
        '''
        if debug(self):
            print "in zooming.drawrectD"
        rapp2, rapp1 = (self.data.resmin.size2/float(self.display.currentd.size2),
                        self.data.resmin.size1/float(self.display.currentd.size1))                  # 
        verticesD = []
        codes, pt0, pt1, pt2, pt3, pt4 = self._makes_poly(llx, lly, urx, ury, rapp = [rapp1, rapp2])
        verticesD +=  [pt0, pt1, pt2, pt3, pt4]                                                     # for window D
        verticesD = np.array(verticesD, float)
        pathD = Path(verticesD, codes)
        pathrectd = PathPatch(pathD, facecolor = 'None', edgecolor = self.gtools.colzoo)
        if layout2:                                                                                 # if layout 2 is asked to be changed
            if debug(self):
                print "shows zoom in D "
            self.rectd = self.display.qmd.axes.add_patch(pathrectd)
            self.display.qmd.fig.canvas.draw()                                                      # draw rectangle in D

    def drawrect(self, llx, lly, urx, ury, layout1 = None , layout2 = None):

        '''
        Draws rectangles with absolute coordinates llx, lly, urx, ury in C and D windows 
        '''
        if debug(self):
            print "in zooming.drawrect"
        self.drawrectC(llx, lly, urx, ury, layout1 = layout1)   
        self.drawrectD(llx, lly, urx, ury, layout2 = layout2)
    
    