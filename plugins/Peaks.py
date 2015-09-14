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
from scipy.optimize import curve_fit
#from spike.util.signal_tools import findnoiselevel, findnoiselevel_2D

def _identity(x):
    "the null function - used as default function argument"
    return x
class Peak(object):
    """ a generic class to store peaks
    defines :
    Id          a unique integer
    intens      The intensity (the height of the largest point)
    area        The area/volume of the peak
    label       a string 
    intens_err  The uncertainty of the previous values
    area_err    ..
    
    """
    def __init__(self, Id, label, intens):
        self.Id = Id
        self.label = label
        self.intens = float(intens)    # intensity
        self.area = 0.0         # area under the curve
        self.intens_err = 0.0      # uncertainty on intensity
        self.area_err = 0.0     # uncertainty on area
        

class Peak1D(Peak):
    """a class to store a single 1D peak
    defines in addition to Peak
    pos         position of the peak in index
    width       width of the peak in index
    pos_err     uncertainty of the previous values
    width_err   ...
    """
    def __init__(self, Id, label, intens, pos):
        super(Peak1D, self).__init__(Id, label, intens)
        self.pos = float(pos)
        self.pos_err = 0.0
        self.width = 0.0
        self.width_err = 0.0
    def report(self, f=_identity):
        """
        return a string with the peak description
        f is a function used to transform the coordinate
        indentity function is default,
        for instance you can use something like
        peaks.report(f=d.axis1.itop)    to get ppm values on a NMR dataset
        order is "id, label, position, intensity"
        """
        return  "%d, %s, %.5f, %.5f, %.2f" % \
            (self.Id, self.label, f(self.pos), self.intens, self.width) 
class Peak2D(Peak):
    """a class to store a single 2D peak
    defines in addition to Peak
    posF1 posF2         positions in F1 and F2 of the peak in index
    widthF1 ...F2       widthes of the peak in index
    posF1_err ...       uncertainty of the previous values
    widthF1_err   ...
    """
    def __init__(self, Id, label, intens, posF1, posF2 ):
        super(Peak2D, self).__init__(Id, label, intens)
        self.posF1 = posF1
        self.posF2 = posF2
        self.widthF1 = 0.0
        self.widthF2 = 0.0
        self.posF1_err = 0.0
        self.posF2_err = 0.0
    def report(self, f1=_identity, f2=_identity):
        """
        print the peak list
        f1, f2 are two functions used to transform respectively the coordinates in F1 and F2
        indentity function is default,
        for instance you can use something like
        peaks.report(f1=s.axis1.itop, f2=s.axis2.itop)    to get ppm values on a NMR dataset
        order is "id, label, posF1, posF2, intensity"
        """
        return "%d, %s, %.5f, %.5f, %.5f" % \
                (self.Id, self.label, f1(self.posF1), f2(self.posF2), self.intens)
        
class Peak1DList(list):
    """
    store a list of peaks
    contains the array version of the Peak1D object :
    self.pos is the numpy array of the position of all the peaks
    and self[k] is the kth Peak1D object of the list
    """
    @property
    def pos(self):
        return np.array( [pk.pos for pk in self] )
    @property
    def intens(self):
        return np.array( [pk.intens for pk in self] )
    @property
    def label(self):
        return [pk.label for pk in self]
    def report(self, f=_identity, file=None):
        """
        print the peak list
        f is a function used to transform respectively the coordinates
        indentity function is default,
        for instance you can use something like
        d.peaks.export(f=s.axis1.itop)    to get ppm values on a NMR dataset
        """
        print ("# Peak list")
        print( "id, label, position, intensity")
        for pk in self:
            print(pk.report(f=f), file=file)
class Peak2DList(list):
    """
    store a list of peaks
    contains the array version of the Peak2D object :
    self.posF1 is the numpy array of the position of all the peaks
    and self[k] is the kth Peak2D object of the list
    """
    @property
    def posF1(self):
        return np.array( [pk.posF1 for pk in self] )
    @property
    def posF2(self):
        return np.array( [pk.posF2 for pk in self] )
    @property
    def intens(self):
        return np.array( [pk.intens for pk in self] )
    @property
    def label(self):
        return [pk.label for pk in self]
    # def __getitem__(self,i):
    #     return Peak2D(Id = self.id[i],
    #             label = self.label[i],
    #             intens = self.intens[i], 
    #             posF1 = self.posF1[i],
    #             posF2 = self.posF2[i] )
    def report(self, f1=_identity, f2=_identity, file=None):
        """
        print the peak list
        f1, f2 are two functions used to transform respectively the coordinates in F1 and F2
        indentity function is default,
        for instance you can use something like
        d.peaks.export(f1=s.axis1.itop, f2=s.axis2.itop)    to get ppm values on a NMR dataset
        the file keyword allows to redirect the output to a file object 
        """
        print ("# Peak list", file=file)
        print( "id, label, posF1, posF2, intensity", file=file)
        for pk in self:
            print(pk.report(f1=f1, f2=f2), file=file)
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
                            #     Id, label, intens, pos        
        pkl = Peak1DList( [ Peak1D(i, str(i), intens, pos) \
            for i, pos, intens in zip( range(len(listint)), list(listpkF1), list(listint) ) ] )
        npkd.peaks = pkl
    elif npkd.dim == 2:
        if threshold is None:
            threshold = 3*npkd.std()
        listpkF1, listpkF2, listint = peaks2d(npkd, threshold=threshold, zoom=zoom)
                            #     Id, label, intens, posF1, posF2 
        pkl = Peak2DList( [ Peak2D(i, str(i), intens, posF1, posF2) \
            for i, posF1, posF2, intens in zip( range(len(listpkF1)), listpkF1, listpkF2, listint ) ] )
        npkd.peaks = pkl
    else:
        raise NPKError("Not implemented of %sD experiment"%npkd.dim, data=npkd)
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
# centroid measure
def center(x, xo, intens, width):
    """
    the centroid definition, used to fit the spectrum
    x can be a nparray
    FWMH is   sqrt(2) x width.
    """
    return intens*(1 - ((x-xo)/width)**2)
def center2d(yx, yo, xo, intens, widthy, widthx):
    """
    the 2D centroid, used to fit 2D spectra - si center()
    xy is [x0, y_0, x_1, y_1, ..., x_n-1, y_n-1] - is 2*n long for n points,
    returns [z_0, z_1, ... z_n-1]
    """
    y = yx[::2]
    x = yx[1::2]
    return intens*(1 - ((x-xo)/widthx)**2)*(1 - ((y-yo)/widthy)**2)
## previous version - based on polyfit - different API
def _centroid1d(npkd, npoints=3):
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
    npkd.check1D()
    noff = (int(npoints)-1)/2
    if (2*noff+1 != npoints) or (npoints<3):
        raise NPKError("npoints must odd and >2 ",data=npkd)
    for i, pk in enumerate(npkd.peaks):
        xdata = np.arange(pk.pos-noff, pkpk.pos+noff+1)
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
def centroid1d(npkd, npoints=3):
    """
    from peak lists determined with peak()
    realize a centroid fit of the peak summit and width,
    computes Full width at half maximum
    updates in data peak list
    
    TODO : update uncertainties
    """
    from scipy import polyfit
    npkd.check1D()
    noff = (int(npoints)-1)/2
    if (2*noff+1 != npoints) or (npoints<3):
        raise NPKError("npoints must odd and >2 ",data=npkd)
    for pk in npkd.peaks:
        xdata = range(int(round(pk.pos-noff)), int(round(pk.pos+noff+1)))
        ydata = npkd.get_buffer()[xdata]
        try:
            popt, pcov = curve_fit(center, xdata, ydata, p0=[pk.pos, pk.intens, 1.0] ) # fit
        except RuntimeError:
            print (pk.label, "centroid could not be fitted")
        pk.pos = popt[0]
        pk.intens = popt[1]
        pk.width = popt[2]
def centroid2d(npkd, npoints_F1=3, npoints_F2=3):
    """
    from peak lists determined with peak()
    realize a centroid fit of the peak summit and width,
    computes Full width at half maximum
    updates in data peak list
    
    TODO : update uncertainties
    """
    from scipy import polyfit
    npkd.check2D()
    nF1 = npoints_F1
    nF2 = npoints_F2
    noff1 = (int(nF1)-1)/2
    noff2 = (int(nF2)-1)/2
    if (2*noff1+1 != nF1) or (nF1<3) or (2*noff2+1 != nF2) or (nF2<3):
        raise NPKError("npoints must odd and >2 ",data=npkd)
    for pk in npkd.peaks:
        st1 = int(round(pk.posF1-noff1))
        end1 = int(round(pk.posF1+noff1+1))
        st2 = int(round(pk.posF2-noff2))
        end2 = int(round(pk.posF2+noff2+1))
        yxdata = np.array( [(y,x) for y in range(st1, end1) for x in range(st2, end2)] ).ravel()
        # yx = np.array([[y,x] for y in range(1,4) for x in range(10,12)]).ravel()
        # => array([ 1, 10,  1, 11,  2, 10,  2, 11,  3, 10,  3, 11])  will be decoded by center2d
        zdata = npkd.get_buffer()[st1:end1, st2:end2].ravel()
        try:
                                            # yx,           yo, xo, intens, widthy, widthx
            popt, pcov = curve_fit(center2d, yxdata, zdata, p0=[pk.posF1, pk.posF2, pk.intens, 1.0, 1.0] ) # fit
        except RuntimeError:
            print (pk.label, "centroid could not be fitted")
        pk.posF1 = popt[0]
        pk.posF2 = popt[1]
        pk.intens = popt[2]
        pk.widthF1 = popt[3]
        pk.widthF2 = popt[4]
#-------------------------------------------------------
def centroid(npkd, *arg, **kwarg):
    if npkd.dim == 1:
        centroid1d(npkd, *arg, **kwarg)
    elif npkd.dim == 2:
        centroid2d(npkd, *arg, **kwarg)
    else:
        raise Exception("to be done")
    return npkd
#-------------------------------------------------------
def display_peaks(npkd, axis = None, peak_label = False, zoom = None, show = False, f=_identity, f1=_identity, f2=_identity):
    if npkd.dim == 1:
        return display_peaks1D(npkd, axis=axis, peak_label=peak_label, zoom=zoom, show=show, f=f)
    elif npkd.dim == 2:
        return display_peaks2D(npkd, axis=axis, peak_label=peak_label, zoom=zoom, show=show, f1=f1, f2=f2)
    else:
        raise Exception("to be done")
def display_peaks1D(npkd, axis = None, peak_label = False, zoom = None, show = False, f=_identity):
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
        pk = range(len(pl))
    if axis is None:
        plot.plot(f(pl.pos[pk]), pl.intens[pk], "x")
        if peak_label:
            for p in pk:
                plot.text(f(pl.pos[p]), 1.05*pl.intens[p], pl.label[p])
    else:
        raise Exception("to be done")
    if show: plot.show()
    return npkd
def display_peaks2D(npkd, axis = None, peak_label = False, zoom = None, show = False, f1=_identity, f2=_identity):
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
        pk = range(len(pl))
    if npkd.debug>0:  print("plotting %d peaks"%len(pk))
    if axis is None:
        plF1 = pl.posF1 # these are ndarray !
        plF2 = pl.posF2
        plot.plot(f2(plF2[pk]), f1(plF1[pk]), "x")
        if peak_label:
            for p in pk:
                plp = pl[p]
                plot.text(1.01*plp.posF2, 1.01*plp.posF1, plp.label)
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
        self.assertEqual(list(d.peaks.pos) , [ 5.0, 7.0, 10.0, 15.0])
        self.assertEqual(list(d.peaks.intens) , [ 20.0, 8.0, 20.0, 11.0])
    def test_center1d(self):
        x = np.arange(5)
        y = center(x, 2.2, 10.0, 1.2)
        # y == [ -2.36111111e+01  -4.44089210e-15   9.72222222e+00   5.55555556e+00 -1.25000000e+01]
        self.assertAlmostEqual(y[2], 9.72222222)
        d = NPKData(buffer = np.maximum(y,0.0))
        d.peaks = Peak1DList()
        d.peaks.append(Peak1D(0, "0", 9.7, 2 ))
        d.peaks[-1].width = 1.0
        d.centroid(npoints=3)
        self.assertAlmostEqual(d.peaks[0].pos,2.2)
        self.assertAlmostEqual(d.peaks[0].intens,10.0)
        self.assertAlmostEqual(d.peaks[0].width,1.2)
    def test_center2d(self):
        M=np.zeros((20, 20))
        # add one peak at (F1,F2) 5.3, 7.9 with widthes (5.0,1.3) 
        for y in range(1,10):
            for x in range(6,11):
                #                 yo, x0 intens, widthy, widthx
                M[y,x] = center2d(np.array([y,x]), 5.3, 7.9, 20.0, 5.0, 1.3)
        #print (M[1:10,6:11])
        self.assertAlmostEqual(M[2,7],5.87777515)
        d = NPKData(buffer = np.maximum(M,0.0))
        d.peaks = Peak2DList()
        # self, Id, label, intens, posF1, posF2 
        d.peaks.append(Peak2D(0, "0", 18.0, 5, 8 ))
        d.centroid(npoints_F1=5)
        self.assertAlmostEqual(d.peaks[0].posF1,5.3)
        self.assertAlmostEqual(d.peaks[0].posF2,7.9)
        self.assertAlmostEqual(d.peaks[0].intens,20.0)
        self.assertAlmostEqual(d.peaks[0].widthF1,5.0)
        self.assertAlmostEqual(d.peaks[0].widthF2,1.3)

NPKData_plugin("display_peaks", display_peaks)
NPKData_plugin("pp", peakpick)
NPKData_plugin("centroid", centroid)

