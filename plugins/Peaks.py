"""
set of function for Peak detections and display

Very Sloppy - Not finsihed !
"""

from __future__ import print_function
import numpy as np
import unittest

from spike import NPKError
from spike.NPKData import NPKData_plugin, NPKData


#----------------------------------------------------------
def peaks2d(npkd, threshold = 0.1, zoom = None, value = False):
    '''
    Extract peaks from 2d Array dataset
    if value is True, return the magnitude at position (x,y)
    '''
    npkd.check2D()
    #print threshold
    if zoom:        # should ((F1_limits),(F2_limits))
        z1lo=zoom[0][0]
        z1up=zoom[0][1]
        z2lo=zoom[1][0]
        z2up=zoom[1][1]
    else:
        z1lo=0
        z1up=npkd.size1-1
        z2lo=0
        z2up=npkd.size2-1
    buff = npkd.buffer[z1lo:z1up, z2lo:z2up]            # take the zoom window
    listpk=np.where(((buff > threshold*np.ones(buff.shape))&            # thresholding
                    (buff > np.roll(buff,  1, 0)) &         # laplacian - kind of
                    (buff > np.roll(buff, -1, 0)) &
                    (buff > np.roll(buff,  1, 1)) &
                    (buff > np.roll(buff, -1, 1))))
    listpk = [int(z1lo) + listpk[0], int(z2lo) + listpk[1]]         # absolute coordinates
    if value: 
        return listpk[0], listpk[1], npkd.buffer[listpk[0], listpk[1]]
    else:
        return listpk[0], listpk[1] # list f1, f2

def peak(npkd, pos_neg = 1, threshold = 0.1, offset = None):
    """
    first trial for peakpicker
    1D only
    pos_neg = 1 / -1 / 0   :  type of peaks positive / negative / 
    threshold = minimum level, as absolute value
    npkd.peaks : index of the peaks
    npkd.peaks_ordered : index of the ordered peaks from maximum to minimum.
    """
    npkd.check1D()
    a_buf = npkd.get_buffer()
    if npkd.itype == 1:
        a_buf = abs(a_buf)
    a_rest = a_buf[1:-1]
    #valmin = threshold*np.max(np.abs(a_rest))            # minimumvalue to be a peak
    valmin = threshold
    diff = (a_rest-a_buf[:-2])*(a_rest-a_buf[2:])       # positive at extremum
    peaks = diff > 0                              # true at extremum
    if pos_neg == 1:
        peaks = np.logical_and(peaks , (a_rest-a_buf[:-2])>0)       # true at positive extremum
        peaks = np.logical_and(peaks, (a_rest>=valmin))
    elif pos_neg == -1:
        peaks = np.logical_and(peaks , (a_rest-a_buf[:-2])<0)       # true at negative extremum
        peaks = np.logical_and(peaks, (a_rest<=valmin))
    elif pos_neg == 0:
        peaks = np.logical_and(peaks, np.abs(a_rest)>=valmin)
    (ipeaks,) = np.where(peaks)     # returns index list
    npkd.peaks = ipeaks+1   # +1 because to compensate initial shifts
    npkd.peaks_ordered = npkd.peaks[np.argsort(a_buf[npkd.peaks])[::-1]]        # list from maximum to minimum
    return npkd.peaks
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
    """displays peaks generated with peak()"""
    import spike.Display.testplot as testplot
    plot = testplot.plot()
    
    if zoom:
        z0=zoom[0]
        z1=zoom[1]
        if z1 == -1 : z1=npkd.size1
        pk = npkd.peaks[np.where(np.logical_and(npkd.peaks>=z0,npkd.peaks<=z1))]
    else:
        z0 = 0
        z1 = npkd.size1
        pk = npkd.peaks
    if axis is None:
        plot.plot(pk-z0,npkd.buffer[pk], "o")
        if peak_label:
            for p in pk:
                plot.text(1.05*npkd.buffer[p], "%.2f"%axis[p])
    else:
        plot.plot(axis[pk], npkd.buffer[pk], "o")
        if peak_label:
            for p in pk:
                plot.text(axis[p], 1.05*npkd.buffer[p], "%.2f"%axis[p])
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
        thresh = 10
        x,y = d.peaks2d(threshold = thresh)
        print("hou 2D",list(x),list(y))
        self.assertEqual(list(x) , [ 5, 10])
        self.assertEqual(list(y) , [ 7, 12])
        
    def test_peaks1d(self):
        "test 1D peak picker"
        print(self.test_peaks1d.__doc__)
        M = np.zeros((30))
        M[5] = 20
        M[7] = 8
        M[15] = 11
        M[10] = 20
        d = NPKData(buffer = M)
        thresh = 10
        x = d.peak(threshold = thresh)
        self.assertEqual(list(x) , [ 5, 10, 15])
        #self.assertEqual(list(y) , [ 20, 20, 11])

NPKData_plugin("display_peaks", display_peaks)
NPKData_plugin("peak", peak)
NPKData_plugin("centroid1d", centroid1d)
NPKData_plugin("peaks2d", peaks2d)

