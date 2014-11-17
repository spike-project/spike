'''
Created by Lionel Chiron  02/10/2013 
Copyright (c) 2013 __NMRTEC__. All rights reserved.
'''
import os, sys
import numpy as np
from Matplotlib_generictools import*
from util.debug_tools import*
from FTICR import FTICRData
from profile_popup import PROFILE 
from label_2D import PEAKPICK

'''
Graphic tools
Peakpicking, Profiles, Corners.
'''

@dec_class_pr
@decclassdebugging
class GRAPHTOOLS():
    '''
    Class containing tools like peak peaking, profiles etc.. 
    Additionnal tools to the basic zoom one.
    '''
    def __init__(self, paramz, display, data, save, convert):
        print "in GRAPHTOOLS"
        self.paramz = paramz    # zoom parameters
        self.display = display
        self.data = data
        self.save = save
        self.convert = convert  # conversions point/mz
        self.selpeak = False
        self.xpeak = None                                                                   # Numpy list containing the x coordinates of the peaks
        self.ypeak = None                                                                   # Numpy list containing the y coordinates of the peaks
        self.deltfact = None
        self.colzoo = None                                                                  # color zoom
        self.colpt = None #
        self.profile_list = []
        self.maxnblabels = 10   # maximum number of lables for peakpicking.
        self.roundmz = 4    # round value for labels of mz peaks.
        self.namefile = data.namefile
    
    def corner_color(self):
        '''
        Change corner point's color
        '''
        if not(self.data.mode_point) :                                                      # if mode point false
            self.colzoo = 'red'
            self.colpt = 'ro'
        else :
            self.colzoo = 'green'
            self.colpt = 'go'
    
    def createpeaks(self):
        '''
        Peak peaking routine working for "point" and "m/z" formats
        '''
        pp = PEAKPICK(self.data, self.display, self.convert, self.paramz) # , pat.data2D, range_axis = range_axis, ax = ax)
        if len(self.paramz.zoom_coord) == 0 :
            self.paramz.zoom_coord = self.paramz.zoom_coordold
        absmax = self.display.currentd.absmax
        pp.find_peaks(thresh = absmax/self.paramz.scale*0.05, zoom = self.convert.to_npk(), maxpeaks = self.maxnblabels) 
        pp.make_labels(method = 'barycenter') 

    def data2screencoor(self, x, y, llx, lly, urx, ury):
        '''
        passes from data coordinates to screen coordinates. 
        normalizes both on x and y. 
        '''
        x = np.float64(x - llx)/abs(llx - urx)
        y = np.float64(y - lly)/abs(lly - ury)
        return x, y
        
    def make_ptl_mz(self, llx, lly, urx, ury, a, b):
        '''
        Makes the lists ptlx and ptly for m/z mode for profiles.
        Passes by the point domain.
        '''
        if debug(self):
            print "in m/z mode "
        ###############
        mztoi2 = self.display.RESMAX.axes(2).mztoi
        mztoi1 = self.display.RESMAX.axes(1).mztoi
        ##
        hlimx = np.int64(mztoi2(llx))                                                 # low bound in "point" format
        llimx = np.int64(mztoi2(urx))                                                 # high bound in "point" format
        ###
        hlimy = np.int64(mztoi1(lly))                                                 # low bound in "point" format
        llimy = np.int64(mztoi1(ury))                                                 # high bound in "point" format
        
        if debug(self):
            print "hlimx, llimx, hlimy, llimy", hlimx, llimx, hlimy, llimy
        if (hlimx - llimx) > (hlimy - llimy):
            if debug(self):
                print "case (hlimx - llimx) > (hlimy - llimy)"
                print "wider than higher"
            for i in range(llimx, hlimx):                                                           # moving in "point" space with self.canvolution taken on x
                self.ptlx.append(i)  # makes the x coordinates
                yval = int(self.display.RESMAX.axes(1).mztoi(a*self.display.RESMAX.axes(2).itomz(i) + b))
                self.ptly.append(yval)
        else :
            if debug(self):
                print "case  (hlimx - llimx) < (hlimy - llimy)" 
                print "higher than wider"
            for i in range(llimy, hlimy):                                               # moving in "point" space with resolution taken on y
                self.ptly.append(i) # makes the y coordinates
                xval = int(self.display.RESMAX.axes(2).mztoi((self.display.RESMAX.axes(1).itomz(i) - b)/a))
                self.ptlx.append(xval)
        if debug(self):
            print "self.ptlx, self.ptly ", self.ptlx, self.ptly

    def make_ptl_pt(self, llx, lly, urx, ury, a, b):
        '''
        Makes the ptlx and ptly for point mode
        '''
        print "in point mode "
        if (urx-llx) > (ury - lly): 
            for i in range(int(llx), int(urx)):                                                                 # moving in "point" space with resolution taken on x
                self.ptlx.append(i)
                self.ptly.append(int(a*i + b))
        else :      
            for i in range(int(lly), int(ury)):                                                                 # moving in "point" space with resolution taken on y
                self.ptly.append(i)
                self.ptlx.append(int((i - b)/a))

    def right_order_coord(self, llx, lly, urx, ury):
        '''
        Put llx, lly, urx, ury in the right order.
        '''
        if debug(self):
            print "in right_order_coord llx, lly, urx, ury  = ", llx, lly, urx, ury
        llx, lly, urx, ury = min(llx, urx), min(lly, ury), max(llx, urx), max(lly, ury)                         # sort for making ptlx and ptly
        return llx, lly, urx, ury
        
    def param_zoom_right_order(self):
        '''
        retrieve llx, lly, urx, ury from self.paramz.zoom_coord in the right order.
        '''
        if debug(self):
            "##### self.paramz.zoom_coord ", self.paramz.zoom_coord
        try:
            llx , lly, urx, ury =  self.paramz.zoom_coord
        except:
            llx , lly, urx, ury = 0, 0, 1024, 1024
            
        return self.right_order_coord(llx , lly, urx, ury)
            
    def point_mz_right_order(self):
        '''
        for a given mode(pt or m/z) and resolution and self.paramz.zoom_coord returns llx, lly, urx, ury
        '''
        zc = self.paramz.zoom_coord
        if not self.data.mode_point :# if in m/z view
            urx, ury, llx, lly = self.convert.itomz_all(self.display.currentd, zc[0], zc[1], zc[2], zc[3])
        else :
            llx, lly, urx, ury = zc[0], zc[1], zc[2], zc[3]
        return self.right_order_coord(llx, lly, urx, ury)

    def do_profile_coord(self, llx, lly, urx, ury):
        """
        from llx, lly, urx, ury in m/z
        builds the self.ptlx/y coordinates in point of the profile
        If in m/z takes resolution on x axis or y axis according to the one giving the best resolution.
        At the final we have self.ptlx and self.ptly, coordinates of the profile.
        """
        llx, lly, urx, ury = self.right_order_coord(llx, lly, urx, ury)
        a = float(lly - ury)/(llx - urx)  # 
        b = (llx*ury - lly*urx)/(llx - urx) #
        self.ptlx, self.ptly = [], []                                                                           # coordinates line in "point"
        if self.data.mode_point: # 
            self.make_ptl_pt(llx, lly, urx, ury, a, b)                                                          # Makes the ptlx and ptly for point mode
        else:
            self.make_ptl_mz(llx, lly, urx, ury, a, b)                                                          # Makes the ptlx and ptly for m/z mode
    
    def makeFTICRData(self, pts, typeprof):
        '''
        Make profiles with FTICRData format.
        returns data_profile
        '''
        data_profile = FTICRData(dim = 1, buffer = pts)
        data_profile.axes = self.display.currentd.axes
        
        if self.data.mode_point:
            data_profile.type = 'points'
            data_profile.axis1.units = 'points'
        else:
            data_profile.type = 'm/z'
            data_profile.axis1.units = 'm/z'
        data_profile.along = typeprof
        return data_profile

    def plotproftype(self, name_profile, typeprof = None):
        '''
        Makes the profile in m/z.
        typreprof can be "diag", "x" or "y"
        '''
        pts = self.display.RESMAX.buffer[self.ptly, self.ptlx]
        data_profile = self.makeFTICRData(pts, typeprof)
        ax = None
        if typeprof == 'diag':
            ax = np.array(self.ptlx)
        self.profile_list.append(PROFILE(data_profile, self.save, name_profile, self.namefile, ax))
        self.profile_list[-1].show()
        
    def plotprofile(self, profile_type, name_profile = None):
        '''
        Plot profile in m/z mode
        Called in canv_event.release_refrechC
        '''
        if debug(self):
            print "############## in plotprofile "
        if profile_type == 'diag':
            self.plotproftype(name_profile, typeprof = 'diag')
        if profile_type == 'x':
            self.plotproftype(name_profile, typeprof = 'x')
        if profile_type == 'y':
            self.plotproftype(name_profile, typeprof = 'y')
        self.coorrect = []                                                                          # reinitializes the line coordinates

    def clearprofile(self):
        '''
        Remove line of the profile from canvas.
        '''
        plt.clf()

    def drawline(self, llx, lly, urx, ury, layout1 = None , layout2 = None):
        '''
        Draw line with absolute coordinates llx, lly, urx, ury
        Called in interface_action.manual_profile
        '''
        print "in drawline llx, lly, urx, ury ", llx, lly, urx, ury
        vertices = []
        codes = []
        margx = 0; margy = 0
        codes +=  [Path.MOVETO] + [Path.LINETO] + [Path.CLOSEPOLY]
        pt0 = (llx + margx, lly + margy); pt1 = (urx + margx, ury - margy)
        pt2 = (0, 0)
        vertices +=  [pt0, pt1, pt2]
        vertices = np.array(vertices, float)
        path = Path(vertices, codes)
        pathlineC = PathPatch(path, facecolor = 'None', edgecolor = 'black')
        if layout1:                                                                                 # if layout 1 is asked to be changed
            print "add patch"
            print 'pathlineC ', pathlineC
            self.linec = self.display.qmc.axes.add_patch(pathlineC)                                 # adds line with patch.. 
            self.display.qmc.fig.canvas.draw()
