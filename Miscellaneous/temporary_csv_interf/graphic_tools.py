'''
Created by Lionel Chiron  02/10/2013 
Copyright (c) 2013 __NMRTEC__. All rights reserved.
'''

import os,sys
from matplotlib import pyplot as plt
import numpy as np
import pickle as pkl
import glob
import scipy as sp
import scipy.interpolate as spinterp
import scipy.fftpack as fft
import time
from time import sleep
from scipy.linalg import norm
from matplotlib.patches import PathPatch, Ellipse, Rectangle
from debug_tools import*

'''
Tools for making peakpicking, finding centroids etc..
    Label
    Centroids
    Peakpicking.

'''
@dec_class_pr
@decclassdebugging
class Label(object):

    '''
    
    '''
    def  __init__(self, ax, ax_name, thresh):
        self.ax = ax
        self.name = ax_name
        self.color_label = (1.0, 0.7, 0.7)
        self.pt = polygons()
        self.thresh = thresh
        self.nb_print_labels = 10
        self.ratio_dist_min = 3.5
        
    def plot_lim_poly(self, pm, p_inside = False, p_whole = False, p_divise = False):
        '''
        Plot the limits of the polygons.
        '''
        if p_inside:
            x, y = pm.vec_to_2lists(pm.poly_inside) 
            self.ax['pkl'].plot(x, y, '--') # plot the poly in which labels must be contained.
        if p_whole:
            xw, yw = pm.vec_to_2lists(pm.whole_poly)
            self.ax['pkl'].plot(xw, yw, '--') # plot the whole poly
        if p_divise:
            for poly in pm.list_poly: # plot all the polys which are divisiion of the whole polygon.
                x, y = pm.vec_to_2lists(poly)
                self.ax['pkl'].plot(x, y, '--')
                
    def label_above_arrow(self):
        '''
        Plot label from bottom to up to avoid overlap of the arrows on the labels.
        '''
        mix_list = np.hstack([self.peaks_vec, self.label_vec]) #b = a[a[:,1].argsort()]
        self.mix_list_sorted = mix_list[mix_list[:,3].argsort()][::-1] # sort and reverse for plotting labels on arrows.. 
    
    def remove_doublon_begin_end(self):
        for poly_dr in self.list_poly_drawn:
            for axis in ['x', 'y']:
                poly_dr[axis].pop() # removing last point equal to first for test.
    
    def in_polys(self, pt):
        '''
        Check if point is in the polys
        '''
        print "###### pt observed ",pt
        inside = False
        for poly_dr in self.list_poly_drawn:
            poly = np.array(zip(poly_dr['x'], poly_dr['y']))
            if self.pt.inside(pt, poly): # test if point is inside
                print "inside poly !!!! ", poly
                inside = True
                break
        if not inside :
            print "outside the polys "
        return inside
        
    def reduce_nb_labels(self, l):
        '''
        Avoids having too much labels in the plot
        '''
        if debug(self):
            print "type(l) ", type(l)
        if len(l) > self.nb_print_labels:
            self.pr('l')
            #sys.exit()
            indices = l[:,1].argsort()
            print indices
            l = l[indices][-self.nb_print_labels:]
        else:
            pass
        return l
    
    def prepare_peaks(self, xmz, curv, use_poly):
        if use_poly:
            self.list_poly_drawn = list_poly_drawn
            self.remove_doublon_begin_end()
            self.pr('self.list_poly_drawn')
        index_peaks = peaks1d(curv, self.thresh) # extraction of the peaks
        spec = np.array(zip(xmz, curv)) # making pairs of coordinates
        self.size_max = [spec[:,0].max(), spec[:,1].max()]
        print "########### making list of interesting peaks "
        if use_poly:
            l = np.array([ list(spec[i]) for i in index_peaks if curv[i] > self.thresh and self.in_polys(spec[i])])[::-1] #
        else :
            l = np.array([ list(spec[i]) for i in index_peaks if curv[i] > self.thresh ])[::-1] #
        if debug(self):
            self.pr('l')
        l = self.reduce_nb_labels(l) # reduces the nb of label to self.nb_print_labels
        if debug(self):
            print "type(l) ",type(l)
        self.peaks_vec = l # peaks to be treated with phsyical motor.
        if debug(self):
            self.pr('self.peaks_vec')
    
    def plot(self, xmz, curv, list_poly_drawn = [], use_poly = False, plot_poly = False):
        '''
        Plot the labels
        Seach best position for labels with physical motor.
        '''
        self.prepare_peaks(xmz, curv, use_poly)
        self.pm = phys_motor(self.peaks_vec, fact_repulse = 0.1,
            min_dist = int(xmz.max()/self.ratio_dist_min), size_max = self.size_max)
        self.pm.search_debug = True
        self.pm.search() # search best position for labels
        if plot_poly :
            self.plot_lim_poly(self.pm, p_whole = True) # plot the limit of each polygon and of the polygon around the peaks.
        self.label_vec = self.pm.label_vec # saving the solution from the physical motor.. 
        self.label_above_arrow() # avoid having labels behind arrows
        self.peaklabel()
    
    def make_ellipse(self, ax, xcenter, ycenter):
        '''
        Ellipse around label, showing the limit for collision.
        '''
        a, b = (self.pm.min_dist, self.pm.min_dist/self.pm.anis_norm[1])
        ell = Ellipse((xcenter, ycenter), a, b, angle = 0, linewidth = 1, fill = False, zorder = 2)
        ax.add_patch(ell)
    
    def test_in(self, x, y, xylim):
        print "in test_in "
        print "xylim ",xylim
        print "x > xylim[0][0] and x < xylim[0][1] ", x > xylim[0][0] and x < xylim[0][1]
        print "y > xylim[1][0] and x < xylim[1][1] ", y > xylim[1][0] and y < xylim[1][1]
        return x > xylim[0][0] and x < xylim[0][1] and y > xylim[1][0] and y > xylim[1][1]
    
    def change_xtxt(self, x, xylim):
        if x < xylim[0][0]:
            x = xylim[0][0]
        else:
            x = xylim[0][1]
        return x
    
    def annotate(self, name_label,x, y, xtxt, ytxt):
        self.ax[self.name].annotate(name_label, xy = (x, y), xytext = (xtxt, ytxt),
            bbox = dict(boxstyle = "round", facecolor = self.color_label, 
            edgecolor = "0.5", alpha = 0.9, clip_on = True),
            arrowprops = dict(arrowstyle = "-", clip_on = False), )
    
    def make_label(self, x, y, xtxt, ytxt, xylim):
        '''
        called by peaklabel
        '''
        if debug(self):
            self.make_ellipse(self.ax[self.name], xtxt, ytxt)
        print "make_label"
        name_label = str(round(x, 4))
        print "in make label xylim ", xylim
        print "point is ",x, y
        print "label as coordinates ", xtxt, ytxt
        if xylim is not None : # limit zoom exists
            print "xylim is not empty"
            if self.test_in(x, y, xylim):
                print "x, y ", x,y, " is inside "
                if not self.test_in(xtxt, ytxt, xylim):
                    ytxt = y
                    xtxt = self.change_xtxt(xtxt, xylim)
                    print "annotate with new xtxt, ytxt ", xtxt,ytxt
                self.annotate(name_label, x, y, xtxt, ytxt)
        else:
            self.annotate(name_label,x, y, xtxt, ytxt)

    def peaklabel(self, xylim = None):
        '''
        plot the sorted and reverse list for avoiding to see arrows on labels..
        '''
        print "peaklabel"
        print "xylim is ",xylim
        for p in self.mix_list_sorted: # read the list containing peaks and labels
            x, y, xtxt, ytxt = p # 
            self.make_label(x, y, xtxt, ytxt, xylim = xylim)

def peaks1d(fid, threshold = 0.1):
    '''
    Extract peaks from 1d from FID
    '''
    print "in peaks1d"
    listpk = np.where(( (fid > threshold*np.ones(fid.shape))&# thresholding
                    (fid > np.roll(fid,  1, 0)) &     
                    (fid > np.roll(fid, -1, 0)) )) # roll 1 and -1 on axis 0
    print "nb peaks ", len(listpk[0])
    print type(listpk[0])
    return listpk[0]

class Subplotnum():
    '''
    Iterator for returning automatically
    the number of the subplot iteratively.
    '''
    def __init__(self, nbsubpltvertic, nbsubplthoriz, fig , xlabel = None, ylabel = None):
        self.numplot = 0
        self.nbsubpltvertic = nbsubpltvertic
        self.nbsubplthoriz = nbsubplthoriz
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.fig = fig
            
    def __iter__(self):
        return self
    
    def make_xlabel(self, ax, xlabel):
        if self.xlabel is None and xlabel is None : 
            ax.set_xticklabels([])
        elif xlabel is not None:
            ax.set_xlabel(xlabel)   
        else:
            ax.set_xlabel(self.xlabel)
        ax.locator_params(axis = 'x', tight = False, nbins = 6)
    
    def make_ylabel(self, ax, ylabel):
        if self.ylabel is None and ylabel is None: 
            ax.set_yticklabels([])
        elif ylabel is not None :
            ax.set_ylabel(ylabel)
        else :
            ax.set_ylabel(self.ylabel)
        ax.locator_params(axis = 'y', tight = False, nbins = 5)
    
    def next(self, xlabel = None, ylabel = None):
        if self.numplot < self.nbsubpltvertic:
            self.numplot += 1
            print "self.numplot ",self.numplot
            ax = self.fig.add_subplot(self.nbsubpltvertic*100+self.nbsubplthoriz*10+self.numplot)
            ## x axis
            self.make_xlabel(ax, xlabel)
            ## y axis
            self.make_ylabel(ax, ylabel)
            return ax

class Centroid(object):
    '''
    Plot centroids
    Methods
        plot
    '''
    def __init__(self, ax, axname, thresh):
        self.ax = ax
        self.name = axname
        self.thresh = thresh
        
    def plot(self, xmz, spec):
        self.xmz = xmz
        self.spec = np.real(spec)
        self.findcentroid()
        self.plotparab()
       
    def __calccoeffs(self, i):
        '''
        '''
        xfit0 = [self.xmz[self.index_peaks[i]-1], self.xmz[self.index_peaks[i]],
                 self.xmz[self.index_peaks[i]+1]] # indices of curv to be taken
        xindex = [self.index_peaks[i]-1, self.index_peaks[i], self.index_peaks[i]+1]
        yfit = self.spec[xindex]
        coeffs = sp.polyfit(xfit0, yfit, 2)
        pospeak = -coeffs[1]/(2*coeffs[0])#/factexp# position of maximum of parabola in m/z
        def parab(x):
            return coeffs[0]*x**2+coeffs[1]*x+coeffs[2]
        return self.index_peaks[i], pospeak, coeffs, parab

    def __calcresolution(self, coeffs):
        '''
        Calculus of resolution
        coeffs : coeeficient of fit
        factexp dilatation coeeficient for more accuracy 
        '''
        c = coeffs[2]/2 + coeffs[1]**2/(8*coeffs[0]) #half height 
        b = coeffs[1]
        a = coeffs[0]
        delt = b**2-4*a*c
        width =  np.sqrt(delt)/(a) # width at half height
        resol = -b/(2*np.sqrt(delt))
        return resol

    def __init_find(self):
        '''
        Initializes lists for findcentroid
        '''
        self.lcentroid = [] # list of peaks corrected
        self.lparab = []
        self.lindexc = []
        self.lheight = []# list of peaks height
        self.lresol = [] # list of resolutions
        
    def __find_append(self, i, l_to_append):
        '''
        Append indices, coordinates and functions to different lists.
        '''
        indexcentr, centroid, coeffs, parab = l_to_append
        self.lindexc.append(indexcentr) # index of centroids in curv
        self.lcentroid.append(centroid) # mz centroids
        self.lheight.append(self.spec[self.index_peaks[i]]) # peaks height
        self.lresol.append(self.__calcresolution(coeffs) ) # list resolutions
        self.lparab.append(parab)
    
    def findcentroid(self):
        """ 
        centroid correction
        Find the maximum from five points fitting parabola by leastsquares. 
        """
        self.__init_find()
        self.index_peaks = peaks1d(self.spec, threshold = self.thresh)
        for i in range(self.index_peaks.shape[0]):# read all the peaks 
            l_to_append = self.__calccoeffs(i) # find coeff of parabolic fit m/z
            self.__find_append(i, l_to_append)

    def plotparab(self):
        '''
        plot the parabolae on the spectrum.
        '''
        for ind, ptxf in enumerate(self.lcentroid) :
            if ptxf > 0 and ptxf < self.xmz.max():
                ptsparab = np.linspace(self.xmz[self.lindexc[ind]-1], self.xmz[self.lindexc[ind]+1], 20)
                self.ax[self.name].plot(ptsparab, self.lparab[ind](ptsparab),'r')
                
class polygons(object):
    '''
    Tools for playing with polygons.
    '''
    def __init__(self, lpeaks = []):
        self.lpeaks = lpeaks
        self.list_poly = []
        
    def test_trigo(self, v1, v2):
        '''
        test if the vectorial product is positive or negative.
        '''
        return (v1[0]*v2[1] - v1[1]*v2[0]) < 0
        
    def inside(self, p, poly):
        '''
        Test if point p is inside poly
        '''
        inside = None      
        for j in range(len(poly)):
            p1, p2 = poly[j-1], poly[j]
            test = self.test_trigo(p2-p, p1-p)
            if j > 0 and  test != test_old: # outside
                inside = False
                break
            elif j == len(poly)-1: # inside
                inside = True
            test_old = test
        return inside
    
    def vec_to_2lists(self, vec):
        '''
        Passing from vector to tuple
        '''
        dez = zip(*vec)
        return np.array(dez[0]), np.array(dez[1])

    def make_poly(self, lpts):
        '''
        Makes a polygon from the list of points lpts.
        close with vertical edges at the beginning and end of the list of points.
        '''
        k = np.array([1,0])
        interm = np.vstack([lpts[0]*k, lpts]) # add the first point
        poly = np.vstack([interm, lpts[-1]*k]) # add the last point
        return poly
        
    def close_and_save_poly(self):
        '''
        '''
        newp = self.make_poly(np.array(self.new_poly)) # close the poly
        self.list_poly.append(newp) # add the resulting poly to the list
        self.nb_poly += 1
        
    def add_pt_poly(self, pt):
        print "convex save add pt ", pt, "to ", self.new_poly
        self.new_poly = np.vstack([self.new_poly, pt]) # continue the poly
    
    def begin_new_poly(self):
        self.new_poly = np.vstack([self.lpeaks[self.pt_index+1], self.lpeaks[self.pt_index+2]])
        
    def make_diff(self):
        ind = self.pt_index
        p1, p2, self.lastpt = self.lpeaks[ind], self.lpeaks[ind+1], self.lpeaks[ind+2]
        self.diff1 = (p2-p1)[1]/float((p2-p1)[0]) # y axis
        self.diff2 = (self.lastpt-p2)[1]/float((self.lastpt-p2)[0]) # y axis

    def if_end_close(self):
        '''
        
        '''
        if self.pt_index == len(self.lpeaks)-2 :
            self.close_and_save_poly()
            self.nb_poly += 1
            
    def case_convex(self):
        '''
        '''
        print "convex"
        self.add_pt_poly(self.lastpt)
        self.pt_index += 1
        self.if_end_close()
    
    def case_concave(self):
        '''
        '''
        print "concave"
        self.close_and_save_poly()
        self.begin_new_poly()
        self.pt_index += 1
        self.if_end_close()
    
    def divide_in_convex(self):
        '''
        Divide concave polygon in convex polygons
        '''
        self.new_poly = [self.lpeaks[0], self.lpeaks[1]]
        #print "in divide_in_convex self.new_poly ",self.new_poly
        self.pt_index = 0
        self.nb_poly = 0
        while self.pt_index < len(self.lpeaks)-2: # read the list of peaks
            self.make_diff()
            if self.diff2 < self.diff1 : # convex
                self.case_convex()
            elif self.diff2 > self.diff1 : #concave.
                self.case_concave()              

@decclassdebugging
class phys_motor(object):
    def __init__(self, lpeaks, fact_repulse = 0.5, min_dist = 5, size_max = None): #poly_inside,
        '''
        '''
        print "############### in phys motor"
        self.label_vec = lpeaks + 10 # begins a little far from peak
        self.fact_repulse = fact_repulse # repulsion factor.. speed
        self.win_size = size_max[1] # size of the observation scene.. 
        self.rep_anis = np.array([1,1]) # anisotropic repulsion
        self.lpeaks = lpeaks # list of peaks in the spectrum.
        self.div_arrow = 10.0
        self.cross_tolerance = 3 # number minimal of poly to be crossed between making mirror effect.
        self.xmax, self.ymax  = size_max
        self.ratio_label = 5 # five time longer than larger.
        self.ratioxy = np.array([1, self.xmax/self.ymax])# ratio width/height
        self.orient = np.array([1,0.85])
        self.min_dist = min_dist
        self.nb_dist_lab = 1
        
    def including_rectangle(self):
        '''
        Makes a rectangle containing the peaks
        '''
        self.whole_poly = self.pt.make_poly(self.lpeaks)
        xmin, ymin = self.whole_poly.min(axis = 0)
        xmax, ymax = self.whole_poly.max(axis = 0)
        gap = 3*self.min_dist
        xmi, xma, yma = xmin - gap, xmax + gap, ymax + gap/self.anis_norm[1]
        include_poly = np.array([[xmi,0], [xmi, yma], [xma, yma], [xma, 0]])
        return include_poly
        
    def near_out_in_side(self, p, poly):
        '''
        return the vector point-middle of the nearest edge.
        Works for stay_inside and repulsion_poly
        '''
        dist = 1e15
        for j in range(len(poly)):
            p1, p2 = poly[j-1], poly[j]
            new_dist = norm(np.cross(p-p2,p1-p2))/norm(p1-p2) # distance from point to edge..
            vec_out = p1 + p2 - 2*p # vector twice p to the middle of [p1,p2]
            if debug(self):
                print "p1, p2 ", p1, p2
                print "new_dist ", new_dist
            if j > 0 and new_dist < dist : # keep the shortest vector
                dist = new_dist
                min_vec_out = vec_out
        if debug(self):
            print " min_vec_out ",min_vec_out
       
        return min_vec_out
    
    def calc_bary(self):
        '''
        barycentrum of the cloud of points
        '''
        self.bary = self.label_vec.mean(axis=0)
    
    def calc_std(self):
        '''
        Standard deviation of the cloud of points
        '''
        self.std = self.label_vec.std(axis=0) 
        
    def vec_to_2lists(self, vec):
        '''
        Passing from vector to tlist
        '''
        return np.array(zip(*vec)[0]), np.array(zip(*vec)[1])
    
    def repulsion_pairs(self, i):
        '''
        Separates each pairs (i,j)
        '''
        for j in range(i+1, len(self.label_vec)):
            if not self.dist_ok(i, j): 
                self.separate_pairs(i, j)
    
    def separate_pairs(self, i, j):
        '''
        self.fact_repulse : repulsion coefficient
        self.anisotropy : anisotropy coefficient
        '''
        self.label_vec[i] += self.dij*self.fact_repulse*self.rep_anis # moves particle i
        self.label_vec[j] += -self.dij*self.fact_repulse*self.rep_anis # moves particle j
        
    def repulsion_poly(self, i):
        '''
        check if particule i is inside on of the  polygons
        '''
        if debug(self):
            print 'in repulsion poly '
        l = self.label_vec[i]
        for poly in self.list_poly:
            if self.pt.inside(l, poly):
                self.nboutpeaks += 1
                self.move_out_in_poly(l, self.whole_poly) #displaces the point out of the whole_poly
    
    def not_crossing_poly(self, i):
        '''
        check if particule i is inside on of the  polygons
        '''
        if debug(self):
            print 'in not_crossing_poly'
        p = self.lpeaks[i]
        l = self.label_vec[i]
        if debug(self):
            if i == 7 : print 'label ',p[0], ' has position ',l
        arrow = [p + n_arr*(l-p)/self.div_arrow  for n_arr in range(int(self.div_arrow))]#arrow vector
        inside = False
        nb_pt_inside = 0
        for pt_arrow in arrow:
            for num_poly, poly in enumerate(self.list_poly):
                if self.pt.inside(pt_arrow, poly):
                    nb_pt_inside += 1
                    if debug(self):
                        print "point label ", p[0], "has a crossing arrow in poly number: ",num_poly
                    break
            if nb_pt_inside > self.cross_tolerance:
                inside = True
                if debug(self):
                    print "nb pts inside is ", nb_pt_inside
                    print "position  before ", self.label_vec[i]
                self.label_vec[i][0] += 2*(p-l)[0]
                if debug(self):
                    print "position  after", self.label_vec[i]
                    print "adding ", 2*(p-l)[0]
                break
            
    def stay_above_peak(self, i):
        '''
        check if particule i is inside on of the  polygons
        '''
        if debug(self):
            print "in repulsion above peak"
            print "len(self.label_vec) ",len(self.label_vec)
            print "len(self.lpeaks[i]) ",len(self.lpeaks)
        l = self.label_vec[i]
        p = self.lpeaks[i]
        if l[1] < p[1]:
            self.label_vec[i][1] += (p[1]-l[1])*1.1
            if debug(self):          
                print "self.label_vec ",self.label_vec

    def keep_dist_to_peak(self, i):
        '''
        Keeps the labels at a standard horizontal distance from peaks
        '''
        if debug(self):
            print "############## correction on length label-peak"
            print "i ",i
        l = self.label_vec[i]
        p = self.lpeaks[i]
        length = self.nb_dist_lab*self.min_dist # length of arrows
    
        line = (l-p)*self.ratioxy
        if debug(self):
            print "before position is ", self.label_vec[i]
        for ax in [0,1]:
            if line[ax] < 0.2*length or line[ax] > length: #line[ax] < 0.5*length or
                if debug(self):
                    print "###########################norm is wrong !!!  "
                corr = self.orient[ax]/self.ratioxy[ax]*np.sign(line[ax])
                self.nblengthwrong += 1
                l[ax] = p[ax] + length*0.6*corr  # taking back the label
        if debug(self):
            print "after position is ", self.label_vec[i]

        if debug(self):
            l = self.label_vec[i]
            p = self.lpeaks[i]
            print "final  norm(l-p)",  norm(l-p)
            
    def move_out_in_poly(self, p, poly):
        '''
        Displaces the particle i outside the polygon.
        '''
        shift = self.near_out_in_side(p, poly)*self.poly_out_gap
        p += shift 
        
    def stay_inside(self, i):
        '''
        keep particle inside the polygon self.poly_inside
        '''
        if debug(self):
            print "in repulsion inside.. "
        p = self.label_vec[i]
        if not self.pt.inside(p, self.poly_inside):
            self.nboutside += 1
            self.move_out_in_poly(p, self.poly_inside)  # displace the point inside.
    
    def stay_upper_plane(self, i):
        '''
        Forbids to particle i having negative values.
        Mirrored position.
        '''
        if debug(self):
            print "in repulsion base"
        if self.label_vec[i][1] < 0:
            self.label_vec[i][1] *= -1

    def calc_dij(self, i, j):
        '''
        distance between point i and j
        calculates also self.dijx and self.dijy
        '''
        self.dij = self.label_vec[i] - self.label_vec[j]
        anis_dij = self.dij*self.anis_norm
        return norm(anis_dij)

    def dist_ok(self, i, j):
        '''
        tests if distance between pairs is ok
        If pair is superimposed, returns False.
        '''
        distij = self.calc_dij(i, j)
        if distij < self.min_dist : # if superimposed
            self.nbsuperp += 1 # a pair more counted.
            return False
        else:
            return True
    
    def init_polygons(self):
        '''
        Prepares polygons for search
        '''
        self.poly_out_gap = 1.1 # factor applied to the vector point-mid_edge for extracting points inside polygons.
        self.pt = polygons(self.lpeaks)
        self.pt.divide_in_convex() # transform polygon made with peaks to a set of convex polygons. 
        self.list_poly = self.pt.list_poly # list of the convex polygons.
        self.poly_inside = self.including_rectangle() # external polygon
    
    def keep_list_nb_notgood(self):
        '''
        Saves number of points not respecting the rules
        '''
        if debug(self):
            print "self.nbsuperp, self.nboutside, self.nboutpeaks ",self.nbsuperp, self.nboutside, self.nboutpeaks
        self.list_nbsuperp.append(self.nbsuperp)
        self.list_nboutside.append(self.nboutside)
        self.list_nboutpeaks.append(self.nboutpeaks)
        self.list_nblengthwrong.append(self.nblengthwrong)
            
    def init_notgood_lists(self):
        '''
        Initialize the list of the nb of points respecting the rules.
        '''
        self.list_nbsuperp = [] # list of the number of superposition at each iteration.
        self.list_nboutside = [] # list of the number outside the ultimate limit at each iteration.
        self.list_nboutpeaks =[] # list of the number outside the polygone defined with the peaks limit at each iteration.
        self.list_nblengthwrong =[] # list of the number of wrong length at each iteration.
        
    def init_notgood(self):
        '''
        Initialize the nb of points respecting the rules
        '''
        self.nbsuperp = 0 # nb  of superpositions
        self.nboutside = 0 # nb  of outside ultimate limit
        self.nboutpeaks = 0 # nb outside the limit defined by the peaks.
        self.nblengthwrong = 0
    
    def search(self, do_show_iter = False, do_show_result = False):
        '''
        Search a solution with no superpositions.
        Makes a repulsion on pairs, a repulsion out of the polygon,
        and repulsion for keeping positive y values.
        '''
        self.anis_norm = np.array([1, self.ratio_label])*self.ratioxy # anisotropy factor for the repulsion 
        print "self.anis_norm  ",self.anis_norm
        print "in search()"
        self.init_polygons()
        self.init_notgood_lists()
        self.nbiter = 0
        while 1 :
            self.init_notgood()
            for i in range(len(self.label_vec)):
                #print "treating label number ", i
                self.repulsion_poly(i) # repulsiton outside polygons
                self.stay_upper_plane(i) # repulsion outside negative y-plane
                self.stay_above_peak(i)
                self.not_crossing_poly(i)
                self.repulsion_pairs(i) # repulsion between pairs
                if self.nbiter > 1 :
                    self.keep_dist_to_peak(i)
                self.stay_inside(i) # repulsion toward inside polygon.
            if debug(self):
                print "self.label_vec ",self.label_vec
                print "self.nbiter  ",self.nbiter 
            self.keep_list_nb_notgood()
            self.nbiter += 1
            #print '#### new iteration number ', self.nbiter
            if do_show_iter:
                self.show_iter()
            if (self.nbsuperp == 0 and self.nboutside == 0
                and  self.nboutpeaks == 0  and self.nblengthwrong == 0) or (self.nbiter == 100): #or (self.nbiter == 50)
                #print "self.nbiter ",self.nbiter
                self.calc_bary()
                self.calc_std()
                break
        if do_show_result:
            self.show_result()

    def show_result(self):
        '''
        Show final result
        '''
        plt.ioff()
        self.show_iter()
        plt.figure()
        plt.plot(self.list_nbsuperp)
        plt.show()
        
    def wind_lim(self):
        '''
        setting limit of the scene
        '''
        plt.ylim(-self.win_size, self.win_size)
        plt.xlim(-self.win_size, self.win_size)
        
    def show_poly(self, poly, line = '-'):
        '''
        extracts x,y and plot
        '''
        x, y = self.vec_to_2lists(poly) # 
        plt.plot(x, y, line) 
        
    def plot_iter(self):
        '''
        Plot points and polygones
        '''
        print "###plot iter"
        self.label_vecx, self.label_vecy = self.vec_to_2lists(self.label_vec) # 
        plt.plot(self.label_vecx, self.label_vecy,'ro')
        self.show_poly(self.whole_poly, line = '--') # the whole polygon defined by the peaks.
        for poly in self.list_poly: # all the small polygons
           self.show_poly(poly)
        self.show_poly(self.poly_inside, line = '--') # including rectangle
    
    def show_iter(self):
        '''
        Show_iterations
        '''
        print "in show iter"
        plt.clf() #clear plot
        self.wind_lim() # set the size of the window.
        self.plot_iter() # plot iterations
        plt.draw()
        
@decclassdebugging
class Mouse_line:
    '''
    Draw polygons with mouse
    '''
    def __init__(self, ax, plt, canvas):
        print "init Mouse_line"
        self.line, = ax.plot([0], [0])
        self.canvas = canvas
        print "dir(self.canvas) ",dir(self.canvas)
        self.ax = ax
        self.plt = plt
        self.list_polygon = [{'x' :[], 'y' :[]}]
        self.cid = self.canvas.mpl_connect('button_press_event', self.click)
        self.rid = self.line.figure.canvas.mpl_connect('button_release_event', self.click)
        self.omid = self.canvas.mpl_connect('motion_notify_event', self.on_motion)
        
        self.pressed = False
    
    def test_trigo(self, v1, v2):
        '''
        test if the vectorial product is positive or negative.
        '''
        return (v1[0]*v2[1] - v1[1]*v2[0]) < 0
    
    def test_trigo_ij(self, poly, i, j):
        '''
        Returns test trigo for pair i,j
        '''
        p1 = np.array([poly['x'][i+1]-poly['x'][i], poly['y'][i+1]-poly['y'][i]])
        p2 = np.array([poly['x'][j+1]-poly['x'][j], poly['y'][j+1]-poly['y'][j]])
        return self.test_trigo( p1, p2)
    
    def test_sense_ij(self, poly, i):
        print "poly['x'][i+1]-poly['x'][i] ",poly['x'][i+1]-poly['x'][i]
        if poly['x'][i+1]-poly['x'][i] < 0 :
            return False
        else:
            return True
    
    def reduce_concave(self, test_deb, poly):
        i = 0
        nberror = 0
        while 1:
            i += 1
        
            if self.test_trigo_ij(poly, i, i+1) != test_deb: # if segment as same test as beginning.. 
                nberror += 1
                for axis in ['x', 'y']:
                    poly[axis][i+1] = (poly[axis][i] + poly[axis][i+2])/2
            if i == len(poly['x'])-3:
                print "nberror",nberror
                break
        return nberror

    def show_last_first(self, poly):
        print 'last is ', poly['x'][-1], poly['y'][-1]
        print 'first is ', poly['x'][0], poly['y'][0]
    
    def shift(self, l):
        for i in ['x', 'y']:
            l[i].pop()
            l[i].insert(0,l[i].pop())
            l[i].append(l[i][0])
            
    def circ(self, poly, step):
        for i in range(step):
            self.shift(poly)
      
    def correc_poly(self, poly):
        if debug(self): print "in correc_poly"
        i = 0
        test_deb = self. test_sense_ij(poly, 0)
        print "test_deb is ", test_deb
        nberror = 1
        for i in range(50):
            nberror = self.reduce_concave(test_deb, poly)
        self.circ(poly, 5)
        for i in range(50):
            nberror = self.reduce_concave(test_deb, poly)
            
    def line_draw(self, event):
        '''
        save and draw polygons
        '''
        if event.inaxes != self.line.axes: return
        poly = self.list_polygon[-1]
        poly['x'].append(event.xdata)
        poly['y'].append(event.ydata)
        self.draw_polygons(poly)
    
    def init_line(self, event):
        '''
        Begin new polygon
        '''
        for i in ['x', 'y']:
            self.list_polygon[-1][i] = []
    
    def add_last_point(self, poly):
        '''
        Adding last point when button is released
        '''
        for i in ['x', 'y']:
            poly[i].pop()
            poly[i].append(poly[i][0])
            
    def plotpoint(self, poly, i, col):
        self.ax.plot(poly['x'][i], poly['y'][i], col)
    
    def plot_last_points(self, poly):
        '''
        Debug tool
        '''
        col = iter(['ko','ko','go','ro'])
        pts = iter([-1,-2,-3,-4])
        for i in pts:
            self.plotpoint(poly, i, col.next())
    
    def plot_first_points(self, poly):
        '''
        Debug tool
        '''
        col = iter(['kx','rx'])
        pts = iter([0,1])
        for i in pts:
            self.plotpoint(poly, i, col.next())
        
    def draw_last_polygons(self, poly):
        '''
        final drawing
        '''
        self.correc_poly(poly)
        self.ax.plot(poly['x'], poly['y'])
        print "poly", poly
        if debug(self):
            self.plot_last_points(poly)
            self.plot_first_points(poly)
        self.plt.draw()
      
    def draw_polygons(self, poly):
        '''
        Real time draw
        '''
        self.line.set_data(poly['x'], poly['y'])
        self.line.figure.canvas.draw()
    
    def finish_line(self, event):
        '''
        Close the polygon
        '''
        poly = self.list_polygon[-1]
        self.add_last_point(poly)
        self.draw_last_polygons(poly)
        self.list_polygon.append({'x' :[], 'y' :[]})
        print "self.list_polygon ",self.list_polygon
        self.draw_polygons(poly)

    def click(self, event):
        '''
        First click begin the polygon, second click finish the polygon
        '''
        if debug(self): print "in click "
        if not self.pressed: # if not in a current drawing
            self.init_line(event)
        self.line_draw(event)
        self.pressed = not self.pressed
        if not self.pressed and len(self.list_polygon[-1]) > 1: # if end of drawing
            self.finish_line(event)
    
    def on_motion(self, event):
        '''
        when moving mouse save data in polygon list
        '''
        if debug(self):
            print "in on_motion Mouse_line"
        if self.pressed:
            self.line_draw(event)

