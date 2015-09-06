#!/usr/bin/env python 
# encoding: utf-8

"""
set of function for Peak detections and display

Very First functionnal - Not finished !

Sept 2015 M-A Delsuc
"""

from __future__ import print_function
import numpy as np
import unittest

from spike import NPKError
from spike.NPKData import NPKData_plugin, NPKData
from spike.util.counter import timeit
#from spike.util.signal_tools import findnoiselevel, findnoiselevel_2D

def _identity(x):
    "the null function - used as default function argument"
    return x
class Peak(object):
    """ a generic class to store peaks"""
    def __init__(self, Id, label, intens):
        self.Id = Id
        self.label = label
        self.intens = intens

class Peak1D(Peak):
    """a class to store a single 1D peak"""
    def __init__(self, Id, label, intens, pos):
        super(Peak1D, self).__init__(Id, label, intens)
        self.pos = pos
class Peak2D(Peak):
    """a class to store a single 2D peak"""
    def __init__(self, Id, label, intens, posF1, posF2 ):
        super(Peak2D, self).__init__(Id, label, intens)
        self.posF1 = posF1
        self.posF2 = posF2
        
class Peak1DList(object):
    """
    store a list of peaks
    contains the array version of the Peak1D object :
    self.pos is the numpy array of the position of all the peaks
    and self[k] is the kth Peak1D object of the list
    """
    def __init__(self):
        self.id = None
        self.pos = None
        self.intens = None
    def len(self):
        return len(self.id)
    def __len__(self):
        return len(self.id)
    def __getitem__(self,i):
        return Peak1D(Id = self.id[i],
                label = self.label[i],
                intens = self.intens[i], 
                pos = self.pos[i])
    def export(self, f=_identity):
        """
        print the peak list
        f is a function used to transform respectively the coordinates
        indentity function is default,
        for instance you can use something like
        d.peaks.export(f=s.axis1.itop)    to get ppm values on a NMR dataset
        """
        print ("# Peak list")
        print( "id, label, position, intensity")
        for i in range(len(self)):
            print("%d, %s, %.5f, %.5f" % \
                (self.id[i], self.label[i], f(self.pos[i]), self.intens[i]) )
class Peak2DList(object):
    """
    store a list of peaks
    contains the array version of the Peak2D object :
    self.posF1 is the numpy array of the position of all the peaks
    and self[k] is the kth Peak2D object of the list
    """
    def __init__(self):
        self.id = None
        self.posF1 = None
        self.posF2 = None
        self.intens = None
    def len(self):
        return len(self.id)
    def __len__(self):
        return len(self.id)
    def __getitem__(self,i):
        return Peak2D(Id = self.id[i],
                label = self.label[i],
                intens = self.intens[i], 
                posF1 = self.posF1[i],
                posF2 = self.posF2[i] )
    def export(self, f1=_identity, f2=_identity):
        """
        print the peak list
        f1, f2 are two functions used to transform respectively the coordinates in F1 and F2
        indentity function is default,
        for instance you can use something like
        d.peaks.export(f1=s.axis1.itop, f2=s.axis2.itop)    to get ppm values on a NMR dataset
        """
        print ("# Peak list")
        print( "id, label, posF1, posF2, intensity")
        for i in range(len(self)):
            print("%d, %s, %.5f, %.5f, %.5f" % \
                (self.id[i], self.label[i], f1(self.posF1[i]), f2(self.posF2[i]), self.intens[i]) )
#     npkd.peaks_ordered = npkd.peaks[np.argsort(a_buf[npkd.peaks])[::-1]]        # list from maximum to             
def _peaks2d(npkd, threshold = 0.1, zoom = None, value = False, zones=0):
    '''
    Extract peaks from 2d Array dataset
    if value is True, return the magnitude at position (x,y)
    '''
    if not zones == 0:
        _peaks2d(npkd, threshold = threshold, zoom = zoom, value = value)
    else:
        print("** assuming no zoom predefined **")
        for z in range(zones):
            print (z)

#----------------------------------------------------------
def peakpick(npkd, threshold = None, zoom = None):
    """
    performs a peak picking of the current experiment
    threshold is the level above which peaks are picked
        None (default) means that 3*(noise level of dataset) will be used - using d.std() as proxy for noise-level
    zoom defines the region on which detection is made (same syntax as in display)
        None (default) means whole data
    """
    if npkd.dim == 1:
        if threshold is None:
            threshold = 3*npkd.std()
        listpkF1, listint = peaks1d(npkd, threshold=threshold, zoom=zoom)
        pp = Peak1DList()
        pp.id = range(len(listpkF1))
        pp.pos = listpkF1
        pp.intens = listint
        pp.label = ["%d"%i for i in pp.id]
        npkd.peaks = pp
    elif npkd.dim == 2:
        if threshold is None:
            threshold = 3*npkd.std()
        listpkF1, listpkF2, listint = peaks2d(npkd, threshold=threshold, zoom=zoom)
        pp = Peak2DList()
        pp.id = range(len(listpkF1))
        pp.posF1 = listpkF1
        pp.posF2 = listpkF2
        pp.intens = listint
        pp.label = ["%d"%i for i in pp.id]
        npkd.peaks = pp
    else:
        raise NPKError("Not implemented of %sD experiment"%npkd.dim)
    return npkd
        
#----------------------------------------------------------
#@timeit
def peaks2d(npkd, threshold, zoom):
    '''
    code for NPKData 2D peak picker
    '''
    npkd.check2D()
    #print threshold
    if zoom:        # should ((F1_limits),(F2_limits))
        z1lo=zoom[0][0]
        z1up=zoom[0][1]
        z2lo=zoom[1][0]
        z2up=zoom[1][1]
        if z1up < 0 : z1up = npkd.size1 + z1up
        if z2up < 0 : z2up = npkd.size2 + z2up
    else:
        z1lo=0
        z1up=npkd.size1-1
        z2lo=0
        z2up=npkd.size2-1
    buff = npkd.get_buffer()[z1lo:z1up, z2lo:z2up]            # take the zoom window
    if npkd.itype != 0:
        buff = abs(buff)
    
    # listpk=np.where(((buff > threshold*np.ones(buff.shape))&            # thresholding
    #                 (buff > np.roll(buff,  1, 0)) &         # laplacian - kind of
    #                 (buff > np.roll(buff, -1, 0)) &
    #                 (buff > np.roll(buff,  1, 1)) &
    #                 (buff > np.roll(buff, -1, 1))))
    # listpkF1 = int(z1lo) + listpk[0]
    # listpkF2 = int(z2lo) + listpk[1]        # absolute coordinates

    # this second version is about 2x faster on my machine
    tbuff = buff[1:-1,1:-1]  # truncated version for shifting
    listpk=np.where(((tbuff > threshold*np.ones(tbuff.shape))&            # thresholding
                    (tbuff > buff[:-2 , 1:-1]) &         # laplacian - kind of
                    (tbuff > buff[2: , 1:-1]) &
                    (tbuff > buff[1:-1 , :-2]) &
                    (tbuff > buff[1:-1 , 2:])))
    listpkF1 = int(z1lo) + listpk[0] +1
    listpkF2 = int(z2lo) + listpk[1] +1        # absolute coordinates

    listint = npkd.buffer[listpkF1, listpkF2]
    return listpkF1, listpkF2, listint
def peaks1d(npkd, threshold, zoom):
    """
    code for NPKData 2D peak picker
    """
    npkd.check1D()
    if zoom:        # should (left,right)
        z1lo=zoom[0]
        z1up=zoom[1]
        if z1up < 0 : z1up = npkd.size1 + z1up
    else:
        z1lo=0
        z1up=npkd.size1-1
    buff = npkd.get_buffer()[z1lo:z1up]            # take the zoom window
    if npkd.itype == 1:
        buff = abs(buff)
    tbuff = buff[1:-1]
    listpk = np.where(((tbuff > threshold*np.ones(tbuff.shape))&# thresholding
                    (tbuff > buff[:-2]) &     
                    (tbuff > buff[2:]) )) # roll 1 and -1 on axis 0
#    return listpk[0]
    listpkF1 = int(z1lo) + listpk[0] +1
    listint = npkd.buffer[listpkF1]
    return listpkF1, listint
#-------------------------------------------------------
def centroid1d(npkd, npoints=3):
    """
    from peak lists determined with peak()
    realize a centroid fit of the peak summit and width,
    computes Full width at half maximum
    creates lists npkd.centered_peaks and npkd.width_peaks
    
    Temporary  so far,
        only based on regular sampling, not unit axis.
        ah-hoc structure, waiting for a real PEAK object
    """
    from scipy import polyfit
    noff = (int(npoints)-1)/2
    npkd.centered_peaks = []
    npkd.width_peaks = []
    if (2*noff+1 != npoints) or (npoints<3):
        raise NPKError("npoints must odd and >2 ",data=npkd)
    for i, pk in enumerate(npkd.peaks):
        xdata = np.arange(pk-noff, pk+noff+1)
        ydata = npkd.get_buffer()[xdata]
        coeffs = polyfit(xdata, ydata, 2)
        pospeak = -coeffs[1]/(2*coeffs[0])#/factexp# position of maximum of parabola in m/z
        c = coeffs[2]/2 + coeffs[1]**2/(8*coeffs[0]) #half height 
        b = coeffs[1]
        a = coeffs[0]
        delt = b**2-4*a*c
        width =  np.sqrt(delt)/abs(a) # width at half height
        npkd.centered_peaks.append(pospeak)
        npkd.width_peaks.append(width)
    npkd.centered_peaks = np.array(npkd.centered_peaks)
    npkd.width_peaks = np.array(npkd.width_peaks)
#-------------------------------------------------------
def display_peaks(npkd, axis = None, peak_label = False, zoom = None, show = False):
    if npkd.dim == 1:
        return display_peaks1D(npkd, axis=axis, peak_label=peak_label, zoom=zoom, show=show)
    elif npkd.dim == 2:
        return display_peaks2D(npkd, axis=axis, peak_label=peak_label, zoom=zoom, show=show)
    else:
        raise Exception("to be done")

def display_peaks1D(npkd, axis = None, peak_label = False, zoom = None, show = False):
    """displays 1D peaks"""
    import spike.Display.testplot as testplot
    plot = testplot.plot()
    pl = npkd.peaks
    if zoom:
        z0=zoom[0]
        z1=zoom[1]
        if z1 < 0 : z1 = npkd.size1 + z1
        pk = [i for i,p in enumerate(pl) if p.pos>=z0 and p.pos<=z1]
    else:
        z0 = 0
        z1 = npkd.size1
        pk = xrange(pl.len())
    if axis is None:
        plot.plot(pl.pos[pk], pl.intens[pk], "x")
        if peak_label:
            for p in pk:
                plot.text(pl.pos[p], 1.05*pl.intens[p], pl.label[p])
    else:
        raise Exception("to be done")
    if show: plot.show()
    return npkd
def display_peaks2D(npkd, axis = None, peak_label = False, zoom = None, show = False):
    """displays 2D peaks"""
    import spike.Display.testplot as testplot
    plot = testplot.plot()
    pl = npkd.peaks
    if zoom:        # should be ((F1_limits),(F2_limits))
        z1lo=zoom[0][0]
        z1up=zoom[0][1]
        z2lo=zoom[1][0]
        z2up=zoom[1][1]
        pk = []
        for p in xrange(pl.len()):
            plp = pl[p]
            if plp.posF1>=z1lo and plp.posF1<=z1up and plp.posF2>=z2lo and plp.posF2<=z2up:
                pk.append(p)
    else:
        z1lo=0
        z1up=npkd.size1-1
        z2lo=0
        z2up=npkd.size2-1
        pk = xrange(pl.len())
    if axis is None:
        plot.plot(pl.posF2[pk], pl.posF1[pk], "x")
        if peak_label:
            for p in pk:
                plot.text(1.01*pl.posF2[p], 1.01*pl.posF1[p], pl.label[p])
    else:
        raise Exception("to be done")
    if show: plot.show()
    return npkd

class PeakTests(unittest.TestCase):
    def test_peaks2d(self):
        "test 2D peak picker"
        print(self.test_peaks2d.__doc__)
        M=np.zeros((30, 30))
        M[5,7] = 20
        M[10,12] = 20
        d = NPKData(buffer = M)
        d.pp() #threshold = 10)  3*d.std is just right
        self.assertEqual(list(d.peaks.posF1) , [ 5, 10])
        self.assertEqual(list(d.peaks.posF2) , [ 7, 12])
        self.assertEqual(list(d.peaks.intens) , [ 20.0, 20.0])
        
    def test_peaks1d(self):
        "test 1D peak picker"
        print(self.test_peaks1d.__doc__)
        M = np.zeros((30))
        M[5] = 20
        M[7] = 8
        M[15] = 11
        M[10] = 20
        d = NPKData(buffer = M)
        d.pp(threshold=3)
        self.assertEqual(list(d.peaks.pos) , [ 5, 7, 10, 15])
        self.assertEqual(list(d.peaks.intens) , [ 20.0, 8.0, 20.0, 11.0])

NPKData_plugin("display_peaks", display_peaks)
NPKData_plugin("pp", peakpick)

