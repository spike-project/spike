#!/usr/bin/env python
"""
Illustrate simple contour plotting, contours on an image with
a colorbar for the contours, and labelled contours.

See also contour_image.py.
"""
from util.debug_tools import*  
import numpy as np
import time
from scipy.linalg import norm
from math import sqrt
from itertools import cycle
from Matplotlib_generictools import*
import matplotlib
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
colors = cycle(['bx', 'mx', 'kx', 'rx', 'gx']) # cycling iterator

@dec_class_pr
@decclassdebugging
class PEAKPICK(object):
    '''
    Peakpicking
    '''
    def __init__(self, data, display, convert, paramz): #, data2D,ax):
        self.data = data
        self.display = display
        self.convert = convert
        self.paramz = paramz
        self.nbiter_getfar = 6                  # nb of iteration for getting the labels far from each other and far from peaks.
        self.divdir = 5
        self.nb_wrong = {'label_label': 0, 'peak_label': 0}
        self.list_wrong = {'label_label': [], 'peak_label': []}
        self.colorlab = (1.0, 0.7, 0.7)         # label color
        self.nbchsign = 4       # number of significant decimals for the m/z values

    def test_correct(self, pt0, pt1, correct = False, kind = None):
        '''
        correct permits to chose between correction (True) and simple test (False).. 
        kind be set to 'label_label' or 'peak_label'
        '''
        if debug(self):
            print "in test_correct "
            print "kind is ", kind
        dist_min = min(self.lx, self.ly)
        if debug(self):
            print "pt1.shape ", pt1.shape
        for i in range(pt1.shape[1]):
            if kind == 'label_label':
                deb = i+1
            elif kind == 'peak_label':
                deb = 0
            for j in range(deb, pt1.shape[1]):
                dist  = norm(pt1[:, i]- pt0[:, j])
                if dist < dist_min  :
                    self.nb_wrong[kind] += 1
                    self.list_wrong[kind].append([i, j])
                    if correct:  # if correct is True, makes correction
                        v = pt1[:, i] - pt0[:, j]
                        vdir = v/norm(v)
                        pt0[:, i] += (dist_min - dist)*vdir
        return pt0  

    def correc_label_label(self, correct):
        '''
        Correction between labels
        '''
        if debug(self):
            print "dist_min = ", dist_min
        self.ptfar = self.test_correct(self.ptfar, self.ptfar, correct = correct, kind = 'label_label')
        
    def correc_peak_label(self, correct):
        '''
        Correction between peaks and labels
        '''
        if debug(self):
            print "dist_min = ", dist_min
        self.ptfar = self.test_correct(self.ptfar, self.peaks, correct = correct, kind = 'peak_label')
    
    def decrossing(self):
        '''
        Avoidance of arrow crossing.
        '''
        for i in range(self.peaks.shape[1]):
            for j in range(i+1, self.peaks.shape[1]):
                p0 = self.ptfar[:, i]
                p1 = self.ptfar[:, j]
                v0 = self.peaks[:, i] - p0
                v1 = self.peaks[:, j] - p1
                r = p1 - p0
                mat = np.empty((2,2))
                mat[:,0] = v0; mat[:, 1] = - v1
                mata = mat.copy(); mata[:, 0] = r
                matb = mat.copy(); matb[:, 0] = r
                detglob = np.linalg.det(mat)
                if detglob != 0:
                    a = np.linalg.det(mata)/detglob
                    b = np.linalg.det(matb)/detglob
                    if 0 < a < 1 and  0 < b < 1:
                    
                        print "crossing", i, j
                        print "inversing the labels "
                        for p in [self.ptfar]: #
                            p[:, i], p[:, j] = p[:, j], p[:, i] # swapping the position of the labels
                    
    def baryc_far(self):
        '''
        Makes the points far from the self.barycenter.
        self.lx, self.ly defined from a small division fo the window dimensions.
        '''
        for i in range(self.peaks.shape[1]):
            v = self.peaks[:,i] - self.barycenter
            vdir = v/norm(v) # vector director
            if debug(self):  
                print i
            self.ptfar[0, i] = self.peaks[0, i] + self.lx*vdir[0]
            self.ptfar[1, i] = self.peaks[1, i] + self.ly*vdir[1]
        plt.plot(self.ptfar[0], self.ptfar[1], 'rx')
    
    def barycenter_method(self):
        '''
        Takes the labels far one from another and all the labels far from the peaks
        1) Begins by  making a self.barycenter and getting the labels at the opposite of the vector 'peak to self.barycenter' 
        at the distance dist_min.
        2) Makes a correction on the labels ovelapping the peaks
        3) Makes a correction on the labels ovelapping each other.
        '''
        print "barycenter method "
        self.baryc_far() # takes the points far from self.barycenter.
        
        for i in range(self.nbiter_getfar):
            #self.decrossing()
            if debug(self):
                print "#### doing peak_label"
            self.correc_peak_label(correct = True) # make correction peak-label
            if debug(self):
                print "peak_label finished "
                print "#### doing label_label"
            self.correc_label_label(correct = True) # make correction label-label
            if debug(self):
                print "self.nb_wrong_peak_label", self.nb_wrong['peak_label']
                print "self.nb_wrong_label_label ", self.nb_wrong['label_label']
                print "self.list_wrong_peak_label", self.list_wrong['peak_label']
                print "self.list_wrong_label_label ", self.list_wrong['label_label']
            
    def get_far(self, method = 'barycenter'):
        '''
        Take the labels far one from another and from peaks
        Two methods: barycenter 
        '''
        self.ptfar = np.empty(self.peaks.shape) # Initialize positoin of the labels on the peaks.
        if method == 'barycenter':
            self.barycenter_method()
      
    def peaklabel(self, i):
        '''
        function to label the peaks.
        Plots the labels if the labels are not calculated outside the limited range.
        '''
        absc = str(round(self.peaks[0,i], self.nbchsign ))
        ordo = str(round(self.peaks[1,i], self.nbchsign ))
        labeltxt = str(i) + ': ' + absc + ', ' + ordo
        cnd1 =   self.xmin < self.ptfar[0, i] < self.xmax
        cnd2 =  self.ymin < self.ptfar[1, i] < self.ymax
        if cnd1 and cnd2:
            ann = self.display.qmc.axes.annotate(labeltxt, xy = (self.peaks[0, i], self.peaks[1, i]), xycoords = 'data', 
                        xytext = (self.ptfar[0, i], self.ptfar[1, i]), size = 10, textcoords = 'data',
                        bbox = dict(boxstyle = "round", fc = self.colorlab, ec = "none"),
                        arrowprops = dict(arrowstyle = "->"))

    def find_peaks(self, thresh, zoom, maxpeaks):
        '''
        Find the peaks in the 2D dataset.
        '''
        if debug(self):
            print "find peaks "
            print "zoom is ", zoom
            print "self.paramz.zoom_coord ", self.paramz.zoom_coord 
        y, x, z = self.display.currentd.peaks2d(threshold = thresh, zoom = zoom, value = True)       # 
        if not self.data.mode_point :                                               # if in m/z view pass to m/z coordinates
            x = self.display.currentd.axes(2).itomz(x)
            y = self.display.currentd.axes(1).itomz(y)
        if debug(self):
            print "x, y ", x, y
        ind = np.lexsort((x,z))[::-1]
        x, y = x[ind][:maxpeaks], y[ind][:maxpeaks]     # take only a maximum number of labels self.maxnblabels
        self.peaks = np.array([x, y])           #self.convert_array(self.peaks_raw)              
        for i in range(self.peaks.shape[1]):
            print "self.peaks[0, i], self.peaks[1, i] ", self.peaks[0, i], self.peaks[1, i]
   
    def make_labels(self, method = 'barycenter'):
        '''
        Finds positions of the labels and plots.
        self.lx, self.ly defined from dimension of the window.
        '''
        if debug(self):
            print "##### make labels "
        ########
        self.xmax, self.ymax, self.xmin, self.ymin  =  self.convert.itomz_all(self.display.currentd, *self.paramz.zoom_coord)
        self.lx, self.ly  = (self.xmax - self.xmin )/self.divdir,   (self.ymax - self.ymin )/self.divdir
        ######
        t0 = time.time()
        self.barycenter = np.sum(self.peaks, axis = 1)/self.peaks.shape[1]      # self.barycenter of the peaks in mpl range format.
        if debug(self):
            print "self.barycenter ", self.barycenter
        plt.plot(self.barycenter[0], self.barycenter[1], 'yo')                  # plots self.barycenter as a yellow point
        self.get_far(method = method)                              # go far from barycenter and makes correction to avoid ovelaps (label-label and peak-label)
        for i in range(self.peaks.shape[1]):
            self.display.qmc.axes.plot(self.peaks[0, i], self.peaks[1, i],'r+') 
            self.peaklabel(i)                                   # plots the labels
        t1 = time.time()
        if debug(self):
            print 'time elapsed ', t1-t0
        self.display.qmc.axes.set_xlim(self.xmin, self.xmax)                                                                               # corrects on x
        self.display.qmc.axes.set_ylim(self.ymin, self.ymax)                                                                               # corrects on y
        self.display.qmc.fig.canvas.draw()
        