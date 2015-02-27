from util.debug_tools import*
from Matplotlib_generictools import*
import numpy as np
import scipy.interpolate as spinterp
from mpl_toolkits.mplot3d import axes3d
from Pyside_PyQt4 import*
from zoom_plot import ZOOM_PLOT 
from canvas_event import CANVAS_EVENT
import json
import io

@dec_class_pr
@decclassdebugging
class ZOOMING():                                                                                    # class regrouping methods about zoom
    '''
    Zoom and draw a rectangle around the zooming area
    The coordinates are in "point" format at the root of the treatment
    so as to simplify all the procedure. 
    '''
    def __init__(self, display, interf, data, paramz, gtools, convert, stools, mwind):
        print "in ZOOMING"
        self.display = display
        self.interface = interf
        self.layoutC = self.interface.ui.layoutC
        self.layoutD = self.interface.ui.layoutD
        self.data = data
        self.paramz = paramz
        self.gtools = gtools
        self.convert = convert
        self.select_tools = stools
        self.move_wind = mwind
        self.zplot = ZOOM_PLOT(display, interf, data, paramz, gtools)
        self.canv_event = CANVAS_EVENT(display, interf, data, paramz,\
            gtools, convert, stools, mwind, self)
        self.stretch = False                                                                        # switch to allow stretching of windows by taking the corners
        self.selline = False
        self.profile_type = None
        self.newrectcready = True                                                                   # switch to indicate that a new rectangle can be drawn or not. 
        self.vecshift = []
        self.listrch = None
        self.corner = None
        self.scro = None
        self.plstrch = None
        self.nbpix = None
        self.margin = 20 # margin for zoom window rectangle
    
    def debug_trig(self):
        '''
        Debug for canvas event
        '''
        self.canv_event.recupxy_debug = True
    
    def change_view(self, change_layoutD = False):
        '''
        Changes views.
        '''
        if debug(self):
            print "in zooming.change_view"
        if not change_layoutD:
            self.display.change_resolution(layout1 = self.layoutC)                                              # 
        else: 
            self.display.change_resolution(layout1 = self.layoutC, layout2 = self.layoutD)                      # 
        self.canv_event.interact_with_canvasC()                         # restablishes interaction with canvas.
    
    def change_view_from_list(self):
        '''
        Refreshes views from history.
        '''
        if debug(self):
            print "in zooming.change_view_from_list"
        self.display.afflistco()                                                                    # 
        self.canv_event.interact_with_canvasC()                                                     # restablishes interaction with canvas.
     
    def clearprofile(self):
        plt.clf()

    def on_scroll(self, event):
        '''
        Use of scrolling function to control the level
        '''
        scrollstep = event.step*10
        self.display.scale = 1 +\
            self.display.levfact*float(self.interface.ui.horizontalSlider.value())
        self.interface.ui.horizontalSlider.setValue(self.display.scale + scrollstep) #
        #self.scroll = res.qmc.fig.canvas.mpl_connect('scroll_event', self.paramz.on_motion)
        
    def stretchrect(self, dx, dy):  
        '''
        function to stretch the zoom by taking the corners
        '''
        zc = self.paramz.zoom_coord
        self.zplot.rectd.remove()                                                                   # erase the zoom rectangle in canvas d
        if self.plstrch == 1:                                                                       # if plot stretch is 1 

            self.set(zc[0], zc[1], dx, dy)
        elif self.plstrch == 0 :                                                                    # if plot stretch is 0
            zc = self.paramz.zoom_coord
            self.paramz.zoom_coord = [dx, dy, zc[2], zc[3]]
        self.change_view()                                                                          # calculate new zoom with new resolution only C window changed
        self.corner.pop(0).remove()
        self.corner = self.display.qmd.axes.plot(dx, dy, self.colpt)
        zc = self.paramz.zoom_coord
        if self.plstrch == 1:                                                                       # draw zoom rectangles in windows C and D, stretching from down left corner
            self.zplot.drawrect(zc[0], zc[1], dx, dy,\
                layout1 = self.layoutC, layout2 = self.layoutD)
        elif self.plstrch == 0 :                                                                    # draw zoom rectangles in windows C and D, stretching from up right corner
            self.zplot.drawrect(dx, dy, zc[2], zc[3],\
                layout1 = self.layoutC, layout2 = self.layoutD)
    
    def find_numbpix(self):
        '''
        Numper of pixels in the zoom area.
        '''
        llx, lly, urx, ury = self.gtools.param_zoom_right_order()
        self.nbpix = (ury-lly)*(urx-llx)                                                            # number of pixels

    def plot_zooms(self):
        '''
        Calculates the coordinates of the zoom and prints in both windows
        '''
        if debug(self):
            print "in zooming.cal_coor"
        llx, lly, urx, ury = self.gtools.param_zoom_right_order()
        self.paramz.zoom_coord = [llx, lly, urx, ury]                                            # updates coordinates
        self.zplot.drawrect(llx - self.margin, lly - self.margin, urx + self.margin, ury + self.margin,
                layout1 = self.layoutC, layout2 = self.layoutD)                                     # draws zoom rectangle in both windows
        self.display.disconnect(self.paramz.mouse_motion)                                     # sets free the mouse from C window
        self.display.aff_resolution()                                                   # resolution and zoom coordinates
    
    def zoom_check_size(self):
        '''
        If zoom too small reinitializes ie self.paramz.zoom_coord = [], disconnects the mouse
        and sets flag self.newrectcready to True
        '''
        if debug(self):
            print "in zooming.zoom_check_size"
        zc = self.paramz.zoom_coord
        cnd1 = (abs(zc[2] - zc[0]) < 10)                                                            # condition along x axis                 
        cnd2 = (abs(zc[3] - zc[1]) < 10)                                                            # condition along y axis   
        if self.data.resolu == 'resol1' and (cnd1 or cnd2) :
            self.newrectcready = True

    def press_zoomready_and_in(self):
        '''
        Activated when zoom is ready and pressed is in the zoom window.
        '''
        if debug(self):   
            print "in zooming.press_zoomready_and_in"
        self.paramz.zoom_coord_old = self.paramz.zoom_coord                                         # keeps last window coordinates
        if debug(self): print "in and zoom ready so makes zoom"
        self.change_zoom()                                                                          # makes zoom and put the zoom and resolution in a list
        self.newrectcready = True   # ready for new zoom                                                                 # ready to draw an other zoom window
        self.canv_event.aff_param()

    def press_zoomready_and_out(self):
        '''
        Activated when zoom is ready an pressed is outside the zoom window
        '''
        if debug(self):   
            print "in zooming.press_zoomready_and_out"                                                         
            print "remove grey area in draw C cursor"
        self.newrectcready = True
        try:
            if debug(self):
                print "removing rectangles C and D"
            self.zplot.rectc.remove()                                                                # clicked not in the area .. new rectangle is possible
            self.zplot.rectd.remove()                                                                   # removes zoom window in window D
            #self.zplot.greyzoo.remove() # remove the grey area.. 
        except:
            pass
        self.paramz.zoomready = False
        if debug(self):
            print "self.paramz.zoomready set to False"
        self.zplot.rectc = None
        if debug(self):
            print "Redraw"
        self.display.qmc.fig.canvas.draw()                                                          # refreshes layout C
        self.display.qmd.fig.canvas.draw()                                                          # refreshes layout D
        if debug(self):
            print "keep last zoom"
        last_zoom = self.paramz.listview[int(self.paramz.listview_index)]['zoom']
        self.paramz.zoom_coord = last_zoom                                                          # reinitializes coordinates.
        self.paramz.zoom_coord_old = last_zoom
        if debug(self):
            print "self.paramz.zoom_coord_old ", self.paramz.zoom_coord_old
            print "in zooming.press_zoomready_and_out"
            self.paramz.report()
        #self.paramz.zoom_coord = []

    def press_zoom(self, xpress, ypress):
        '''
        Zoom activated or not after pressing the left mouse button
        '''
        if debug(self):
            print " in zooming.press_zoom "
        llx, lly, urx, ury = self.gtools.param_zoom_right_order()
        inzoom = xpress > llx and xpress < urx and ypress > lly and ypress < ury
        if self.paramz.zoomready: # zoom is ready
            if inzoom :
                self.press_zoomready_and_in()  # pressed in
            else:
                self.press_zoomready_and_out()  # pressed out
                self.change_view()
                #self.paramz.zoom_coord = []
        else :
            
            self.paramz.zoom_coord = []
        if debug(self): print "self.paramz.zoom_coord ", self.paramz.zoom_coord
           
    def change_zoom(self):
        '''
        Makes zoom and resolution if necessary, condition = self.areazoom() >= AREAMIN
        Keeps both zoom and resolution in a list.
        vis.resolu is the name of the current resolution. 
        '''  
        if debug(self):
            print "in zooming.change_zoom "  
        if self.select_tools.zoom:                            
            self.find_numbpix()
            self.display.select_best_resolution()                                                       # changes to the best resolution.
            self.plot_zooms()                                                                           # Calculates new zoom and draw rectangle
            self.display.register_coordinates()                                                         # keeps the coordinates.
            self.canv_event.interact_with_canvasC()

@dec_class_pr
@decclassdebugging
class ZOOM3D(object):
    '''
    Zoom 3D 
    '''
    def __init__(self):
        self.plot_3D_ok = True
        self.trunc = 4
    
    def makemeshfromirreg(self, x, y, z, sizef1, sizef2):
        '''
        Makes the meshgrid from irregular mesh (x, y) in m/z coordinates.
        Called after make_mz_xyz()
        '''
        if debug(self):
            print "in makemeshfromirreg "
            print "x.size, y.size, z.size ", x.size, y.size, z.size
        spline = spinterp.Rbf(x, y, z, function = 'cubic')                                          # Spline interpolation function.
        xi = np.linspace(x.min(), x.max(), sizef1)                                                  # m/z1
        yi = np.linspace(y.min(), y.max(), sizef2)                                                  # m/z2
        xi, yi = np.meshgrid(xi, yi)                                                                # Makes the regular meshgrid in m/z
        zi = spline(xi, yi)                                                                         # Spline interpolation for z component.
        
        return xi, yi, zi
    
    def make_mz_xyz(self, d, pt1min, pt1max, pt2min, pt2max):
        '''
        From the resolution d and frequency coordinates,
        returns the m/z coordinates x, y, z for making the irregular meshgrid.
        '''
        if debug(self):
            print "in make_mz_xyz "
        x, y, z = [], [], []                                                                        # initializes the m/z coordinates x, y, z
        for i in range(pt1min, pt1max) : 
            for j in range(pt2min, pt2max) :                                                        # scans the area in mode point and makes x, y, z in m/z coordinates.
                x.append(d.axes(1).itomz(i))                                                        # pass to m/z
                y.append(d.axes(2).itomz(j))                                                        # pass to m/z
                z.append(d.col(j)[i])
        x, y, z = np.array(x), np.array(y), np.array(z)
        if debug(self):
            print "x.min(), x.max(), y.min(), y.max() ", x.min(), x.max(), y.min(), y.max()
            print "x, y", x, y
        return x, y, z
                
    def makemesh(self, d, pt1min, pt1max, pt2min, pt2max):
        '''
        from frequencies limits pt1min, pt1max, pt2min, pt2max, makes the 3D mesh X, Y, Z.
        Calls make_mz_xyz() then makemeshfromirreg().
        '''
        if debug(self):
            print "makemesh"
            print "pt1min, pt1max, pt2min, pt2max ", pt1min, pt1max, pt2min, pt2max
        x, y, z = self.make_mz_xyz(d, pt1min, pt1max, pt2min, pt2max)                               # Makes the regular m/z mesh from limits in point coordinates.
        sizef1 = pt1max - pt1min                                                                    # number of points in f1 direction
        sizef2 = pt2max - pt2min                                                                    # number of points in f2 direction
        X, Y, zs = self.makemeshfromirreg(x, y, z, sizef1, sizef2)                                  # Makes interpolation 
        np.clip(zs, zs.min(), zs.max()/self.trunc, out = zs)                                        # truncates 
        if debug(self):
            print "X.shape", X.shape
        Z = zs.reshape(X.shape)
        return X, Y, Z
    
    def drawrect(self, llx, lly, urx, ury):
        '''
        Draws rectangle for showing the area where the 3D is performed. 
        '''
        verticesC = []
        codes = []
        margx = 0; margy = 0                                                                        # margin for zoom
        codes +=  [Path.MOVETO] + [Path.LINETO]*3 + [Path.CLOSEPOLY]
        pt0 = (llx + margx, lly + margy); pt1 = (llx + margx, ury - margy)
        pt2 = (urx - margx, ury - margy); pt3 = (urx - margx, lly + margy); pt4 = (0, 0)
        ###                                                                                         # make vertices
        verticesC +=  [pt0, pt1, pt2, pt3, pt4]                                                     # for window C
        verticesC = np.array(verticesC, float)
        pathC = Path(verticesC, codes)
        pathrectc = PathPatch(pathC, facecolor = 'None', edgecolor = 'r')
        
    def plotregion3d(self, d, pt1min, pt1max, pt2min, pt2max, visible = False):
        '''
        Makes the 3D plot from the meshgrid.
        Makes the json that is read by FTICR2D_3d.html
        '''
        if debug(self):
            print "pt1min, pt1max, pt2min, pt2max",pt1min, pt1max, pt2min, pt2max
        X, Y, Z = self.makemesh(d, pt1min, pt1max, pt2min, pt2max)                                   # makes the meshgrid
        data = []
        dimx, dimy = X.shape[0], X.shape[1]
        xflat = np.array([np.arange(dimx) for i in range(dimy)]).flatten()
        yflat = np.array([np.ones(dimx)*i for i in range(dimy)]).flatten()
        zflat = Z.flatten()
        sizemat = xflat.size
        with io.open('Visu/3d.json', 'w', encoding ='utf-8') as f:
              for i in range(sizemat):
                  print xflat[i]
                  data.append({'x':int(xflat[i]),'y':int(yflat[i]),'z': zflat[i]*1e-5 })
              f.write(unicode('data = '))
              f.write(unicode(json.dumps(data, ensure_ascii = False)))
            