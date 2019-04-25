#!/usr/bin/env python 
# encoding: utf-8

"""set of function for Peak detections and display - 1D and 2D

Very First functionnal - Not finished !


Peak1D and Peak2D are simple objects
    with attributes like Id, label, intens(ity), pos(ition), or width
    the only added method is report() (returns a string)

Peak1DList and Peak2DList are python list, with a few added methods
    - report (to stdio or to a file)
    - largest sort in decreasing order of intensity
        other sorts can simply done by peaklist.sort(key = lambda p: p.XXX)
            where XXX is any peak attribute  (see largest code)

Example of usage:

# assuming d is a 2D NPKData / 1D will be just as simple
d.pp()          # computes a peak picking over the whole spectrum using 3 x standard_deviation(d)
                # This is just a detection of all local maxima

# We can be more specific:
d.pp(threshold=5E5, zoom=((700,1050),(300,350)) )     # zoom is always in the currently active unit, defined with d.unit

# this attached the peak list to the dataset as d.peaks,
# it is a list of Peaks2D objects, with some added properties
print( "number of detected peaks: %d" % len(d.peaks))

p0 = d.peaks[0]     # peaks have label, intensity and positions attributes

print( p0.report() )        # and a report method
                            # report has an additional format parameter which enables control on the output

# we can call centroid to improve the accuracy and move the position to center of a fitted (2D) parabola
d.centroid()     


# The peak list can be displayed on screen as simple crosses
d.display_peaks()

# The label can be modififed for specific purposes:
for p in d.peaks:
    if 150 < p.posF2 < 1500 :
        p.label = "%.2f x %.f"%(p.posF1,p.posF2)    # for instance changing to the coordinates for a certain zone
    else:
        p.label = ""                                # and removing elsewhere

d.display_peaks(peak_label=True)

# peak lists can also be reported
d.report_peak()

# but also as a formatted stream, and redirected to a file:
output = open("my_peak_list.csv","w")      # open the file
output.write("# LABEL, INTENSITY, F1, Width, F2, width")
d.report_peak(file=output, format="{1}, {4:.2f}, {2:.7f}, {5:.2f}, {3:.7f}, {6:.2f}")
        # arguments order order is   id, label, posF1, posF2, intensity, widthF1, widthF2
output.close()


Sept 2015 M-A Delsuc
"""

from __future__ import print_function
import numpy as np
import unittest
###
from spike import NPKError
from spike.NPKData import NPKData_plugin, NPKData, flatten, parsezoom
from spike.util.counter import timeit
from scipy.optimize import curve_fit
import warnings
#from spike.util.signal_tools import findnoiselevel, findnoiselevel_2D

debug = 0
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
    pos         position of the peak in index relative to the typed (real/complex) buffer
    width       width of the peak in index
    pos_err     uncertainty of the previous values
    width_err   ...
    """
    report_format = "{}, {}, {:.2f}, {:.2f}"
    full_format = "{}, "*8
    def __init__(self, Id, label, intens, pos):
        super(Peak1D, self).__init__(Id, label, intens)
        self.pos = float(pos)
        self.pos_err = 0.0
        self.width = 0.0
        self.width_err = 0.0
    def report(self, f=_identity, format=None):
        """
        print the peak list
        f is a function used to transform the coordinate
        indentity function is default,
        for instance you can use something like
        peaks.report(f=s.axis1.itop)    to get ppm values on a NMR dataset
        order is "id, label, position, intensity"

        parameters are : 
        Id label positions intens width intens_err pos_err width_err 
        in that order.

        By default returns only the 4 first fields with 2 digits, but the format keyword can change that.
        format values:
        - None or "report", the standard value is used:  "{}, {}, {:.2f}, {:.2f}" 
                  (so only the four first parameters are shown)
        - "full" is all parrameters at full resolution  ( "{}; "*8 )
        - any othe string following the format syta will do.
            you can use any formating syntax. So for instance the following format
            "{1} :   {3:.2f}  F1: {2:.7f} +/- {4:.2f}"
            will remove the Id, show position with 7 digits after the comma, and will show width

        you can change the report and full default values by setting
        pk.__class__.report_format   and   pk.__class__.full_format  which class attributes
        """
        if format is None or format == "report":
            format = self.report_format
        elif format == "full":
            format = self.full_format
        return format.format( self.Id, self.label, f(self.pos), self.intens, self.width,  self.pos_err, self.intens_err, self.width_err)
    def _report(self, f=_identity):
        """
        full report for 1D Peaks
        list is : Id, label, pos intens, widthF1, width, pos_err, intens_err, width_err
        """
        return self.report(f=f, format=self.full_format)
        
class Peak2D(Peak):
    """a class to store a single 2D peak
    defines in addition to Peak
    posF1 posF2         positions in F1 and F2 of the peak in index relative to the typed (real/complex) axis
    widthF1 ...F2       widthes of the peak in index
    posF1_err ...       uncertainty of the previous values
    widthF1_err   ...
    """
    report_format = "{}, {}, {:.2f}, {:.2f}, {:.2f}"
    full_format = "{}, "*12
    def __init__(self, Id, label, intens, posF1, posF2 ):
        super(Peak2D, self).__init__(Id, label, intens)
        self.posF1 = posF1
        self.posF2 = posF2
        self.widthF1 = 0.0
        self.widthF2 = 0.0
        self.posF1_err = 0.0
        self.posF2_err = 0.0
        self.widthF1_err = 0.0
        self.widthF2_err = 0.0
    def report(self, f1=_identity, f2=_identity, format=None):
        """
        print the peak list
        f1, f2 are two functions used to transform respectively the coordinates in F1 and F2
        indentity function is default,
        for instance you can use something like
        peaks.report(f1=s.axis1.itop, f2=s.axis2.itop)    to get ppm values on a NMR dataset
        order is "id, label, posF1, posF2, intensity, widthF1, widthF2"

        printed parameters are : 
        Id label posF1 posF2 intens widthF1 widthF2 posF1_err posF2_err intens_err widthF1_err widthF2_err 

        By default returns only the 5 first fields with 2 digits, but the format keyword can change that.
        format values:
        - None or "report", the standard value is used:  "{}, {}, {:.2f}, {:.2f}, {:.2f}" 
                  (so only the four first parameters are shown)
        - "full" is all parrameters at full resolution  ( "{}; "*12 )
        - any othe string following the format syntxa will do.
            you can use any formating syntax. So for instance the following format
            "{1} :   {4:.2f}  F1: {2:.7f} +/- {5:.2f}  X  F2: {3:.7f} +/- {6:.2f}"
            will remove the Id, show position with 7 digits after the comma, and will show widthes    
        """
        if format is None or format == "report":
            format = self.report_format
        elif format == "full":
            format = self.full_format
        return format.format( self.Id, self.label, f1(self.posF1), f2(self.posF2), self.intens, self.widthF1, self.widthF2, self.posF1_err, self.posF2_err, self.intens_err, self.widthF1_err, self.widthF2_err )
    def _report(self, f1=_identity, f2=_identity):
        """
        full report for 2D Peaks
        list is : Id, label, posF1, posF2, intens, widthF1, widthF2, posF1_err, posF2_err, intens_err, widthF1_err, widthF2_err
        """
        return self.report(f1=f1, f2=f2, format=self.full_format)

class PeakList(list):
    """
    the class generic to all peak lists
    """
    def __init__(self, *arg, **kwds): #threshold=None, source=None):
        super(PeakList, self).__init__(*arg)
        # I can't figure out how to explictly specify a keyword arg with *args:
        #   def __init__(self, *arg, threshold=None, source=None): ...
        # so I use **kwds and sqauwk if something unexpected is passed in.
        # code taken from   lib/python2.7/pstats.py
        #
        # additional kw are source: the originating dataset, threshold: actual value
        if "threshold" in kwds:
            self.threshold = kwds["threshold"]
            del kwds["threshold"]
        else:
            self.threshold = None
        if "source" in kwds:
            self.source = kwds["source"]
            del kwds["source"]
        else:
            self.source = None
        if kwds:
            keys = kwds.keys()
            keys.sort()
            extras = ", ".join(["%s=%s" % (k, kwds[k]) for k in keys])
            raise(ValueError, "unrecognized keyword args: %s" % extras)
    @property
    def intens(self):
        "returns a numpy array of the intensities"
        return np.array( [pk.intens for pk in self] )
    @property
    def label(self):
        "returns an array of the labels"
        return [pk.label for pk in self]
    def len(self):
        return len(self)
    def largest(self):
        "sort the peaklist in decresing order of intensities"
        self.sort(reverse = True, key = lambda p:p.intens)
    def getitem(self,i,j):
        print('getitem')
        return type(self)(list.getitem(self,i,j))
class Peak1DList(PeakList):
    """
    store a list of 1D peaks
    contains the array version of the Peak1D object :
    self.pos is the numpy array of the position of all the peaks
    and self[k] is the kth Peak1D object of the list
    """
    @property
    def pos(self):
        "returns a numpy array of the positions in index"
        return np.array( [pk.pos for pk in self] )
    def report(self, f=_identity, file=None, format=None):
        """
        print the peak list
        f is a function used to transform respectively the coordinates
        indentity function is default,
        for instance you can use something like
        d.peaks.report(f=d.axis1.itop)    to get ppm values on a NMR dataset

        check documentation for Peak1D.report() for details on output format
        """
        print ("# %d in Peak list"%len(self), file=file)
        for pk in self:
            print(pk.report(f=f, format=format), file=file)
    def _report(self, f=_identity, file=None):
        """return full report for 1D peak list
        list is : Id, label, pos, intens, width,  pos_err, intens_err, width_err
        """
        lst = [pk._report(f=f) for pk in self ]
        return "\n".join(lst)
    def display(self, peak_label=False, peak_mode="marker", zoom=None, show=False, f=_identity, color = 'red', markersize=None, figure=None):
        """
        displays 1D peaks
        zoom is in index
        peak_mode is either "marker" or "bar"
        """
        from spike.Display import testplot
        plot = testplot.plot()
        if figure is None:
            fig = plot.subplot(111)
        else:
            fig = figure
        if zoom:
            z0=zoom[0]
            z1=zoom[1]
            pk = [i for i,p in enumerate(self) if p.pos>=z0 and p.pos<=z1]
        else:
            pk = range(len(self))
        if peak_mode == "marker":
            fig.plot(f(self.pos[pk]), self.intens[pk], "x", color=color)
        elif peak_mode == "bar":
            for p in self:
                fig.plot( [f(p.pos),f(p.pos)], [0,p.intens], '-', color=color)
        else:
            raise Exception("wrong peak_mode")
        if peak_label:
            for p in self:
                fig.annotate(p.label,(f(p.pos), p.intens),
                    color='red', xycoords='data', xytext=(0, 10), textcoords='offset points',rotation=40,
                    arrowprops=dict(arrowstyle='-'), horizontalalignment='left', verticalalignment='bottom')
        if show: plot.show()
    def pos2label(self):
        "for FTMS: use pos in current unit, using converion f and set it as label for each peak"
        try:
            f = self.source.axis1.itomz
        except:
            return
        for pk in self:
            pk.label = "%.4f"%(f(pk.pos),)

class Peak2DList(PeakList):
    """
    store a list of 2D peaks
    contains the array version of the Peak2D object :
    self.posF1 is the numpy array of the position of all the peaks
    and self[k] is the kth Peak2D object of the list
    """
    @property
    def posF1(self):
        "returns a numpy array of the F1 positions in index"
        return np.array( [pk.posF1 for pk in self] )
    @property
    def posF2(self):
        "returns a numpy array of the F2 positions in index"
        return np.array( [pk.posF2 for pk in self] )
    def report(self, f1=_identity, f2=_identity, file=None, format=None):
        """
        print the peak list
        f1, f2 are two functions used to transform respectively the coordinates in F1 and F2
        indentity function is default,
        for instance you can use something like
        d.peaks.export(f1=s.axis1.itop, f2=s.axis2.itop)    to get ppm values on a NMR dataset
        the file keyword allows to redirect the output to a file object 
        
        check documentation for Peak2D.report() for details on output format
        """
        print ("# %d in Peak list"%len(self), file=file)
        for pk in self:
            print(pk.report(f1=f1, f2=f2, format=format), file=file)
    def _report(self, f1=_identity, f2=_identity, file=None):
        """return full report for 2D peak list
        list is : Id, label, posF1, posF2, intens, widthF1, widthF2, posF1_err, posF2_err, intens_err, widthF1_err, widthF2_err
        """
        lst = [pk._report(f1=f1, f2=f2) for pk in self ]
        return "\n".join(lst)
    def display(self, axis = None, peak_label=False, zoom=None, show=False, f1=_identity, f2=_identity, color=None, markersize=6):
        """
        displays 2D peak list
        zoom is in index
        """
        import spike.Display.testplot as testplot
        plot = testplot.plot()
        if zoom is not None:
            (z1lo, z1up, z2lo, z2up) = flatten(zoom)
            pk = []
            for p in range(len(self)):
                plp = self[p]
                if plp.posF1>=z1lo and plp.posF1<=z1up and plp.posF2>=z2lo and plp.posF2<=z2up:
                    pk.append(p)
        else:
            pk = range(len(self))
        if debug>0:  print("plotting %d peaks"%len(pk))
        if axis is None:
            plF1 = self.posF1 # these are ndarray !
            plF2 = self.posF2
            plot.plot(f2(plF2[pk]), f1(plF1[pk]), "x", color=color, markersize=markersize)
            if peak_label:
                for p in pk:
                    plp = self[p]
                    plot.text(1.01*plp.posF2, 1.01*plp.posF1, plp.label, color=color, markersize=markersize)
        else:
            raise Exception("to be done")
        if show: plot.show()

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
def peakpick(npkd, threshold = None, zoom = None, autothresh=3.0, verbose=True):
    """
    performs a peak picking of the current experiment
    threshold is the level above which peaks are picked
        None (default) means that autothresh*(noise level of dataset) will be used - using d.std() as proxy for noise-level
    zoom defines the region on which detection is made
        zoom is in currentunit (same syntax as in display)
        None means the whole data
    """
    if threshold is None:
        threshold = autothresh*np.std( npkd.get_buffer().real )
    if npkd.dim == 1:
        listpkF1, listint = peaks1d(npkd, threshold=threshold, zoom=zoom)
                            #     Id, label, intens, pos
        pkl = Peak1DList( ( Peak1D(i, str(i), intens, pos) \
            for i, pos, intens in zip( range(len(listint)), list(listpkF1), list(listint) ) ), \
                    threshold=threshold, source=npkd )
        # use pos as labels - used only for FTMS
        pkl.pos2label()
        npkd.peaks = pkl
    elif npkd.dim == 2:
        listpkF1, listpkF2, listint = peaks2d(npkd, threshold=threshold, zoom=zoom)
                            #     Id, label, intens, posF1, posF2 
        pkl = Peak2DList( ( Peak2D(i, str(i), intens, posF1, posF2) \
            for i, posF1, posF2, intens in zip( range(len(listpkF1)), listpkF1, listpkF2, listint ) ), \
                                    threshold=threshold, source=npkd )
        npkd.peaks = pkl
    else:
        raise NPKError("Not implemented of %sD experiment"%npkd.dim, data=npkd)
    if threshold is None:
        if verbose: print ('PP Threshold:',threshold)
    if verbose: print('PP: %d detected'%(len(npkd.peaks),))
    return npkd

def peaks2d(npkd, threshold, zoom):
    '''
    math code for NPKData 2D peak picker
    '''
    npkd.check2D()
    #print threshold
    print("########## in peaks2d ")
    print("zoom ", zoom)
    z1lo, z1up, z2lo, z2up = parsezoom(npkd, zoom)
    
    print("z1lo, z1up, z2lo, z2up ", z1lo, z1up, z2lo, z2up)
    buff = npkd.get_buffer()[z1lo:z1up, z2lo:z2up]            # take the zoom window
    if npkd.itype != 0:
        buff = buff.real
    
    # listpk=np.where(((buff > gold*np.ones(buff.shape))&            # thresholding
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
def peaks1d(npkd, threshold, zoom=None):
    """
    math code for NPKData 1D peak picker
    """
    npkd.check1D()
    z1, z2 = parsezoom(npkd, zoom)
    buff = npkd.get_buffer()[z1:z2]            # take the zoom window
    if npkd.itype == 1:
        buff = buff.real
    tbuff = buff[1:-1]
    listpk = np.where(((tbuff > threshold*np.ones(tbuff.shape))&# thresholding
                    (tbuff > buff[:-2]) &     
                    (tbuff > buff[2:]) )) # roll 1 and -1 on axis 0
#    return listpk[0]
    listpkF1 = int(z1) + listpk[0] +1
    listint = npkd.get_buffer()[listpkF1]
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
def centroid1d(npkd, npoints=3, reset_label=True):
    """
    from peak lists determined with peak()
    realize a centroid fit of the peak summit and width,
    will use npoints values around center  (npoints has to be odd)
    computes Full width at half maximum
    updates in data peak list
    reset_label when True (default) reset the labels of FTMS datasets
    TODO : update uncertainties
    """
    from scipy import polyfit
    npkd.check1D()
    noff = (int(npoints)-1)/2
    if (2*noff+1 != npoints) or (npoints<3):
        raise NPKError("npoints must odd and >2 ",data=npkd)
    buff = npkd.get_buffer().real
    for pk in npkd.peaks:
        xdata = range(int(round(pk.pos-noff)), int(round(pk.pos+noff+1)))
        ydata = buff[xdata]
        try:
            popt, pcov = curve_fit(center, xdata, ydata, p0=[pk.pos, pk.intens, 1.0] ) # fit
        except RuntimeError:
            print ( "peak %d (id %s) centroid could not be fitted"%(pk.Id, pk.label) )
        pk.pos = popt[0]
        pk.intens = popt[1]
        pk.width = popt[2]
    if reset_label:
        npkd.peaks.pos2label()

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
            print ( "peak %d (label %s) centroid could not be fitted"%(pk.Id, pk.label) )
        pk.posF1 = popt[0]
        pk.posF2 = popt[1]
        pk.intens = popt[2]
        pk.widthF1 = popt[3]
        pk.widthF2 = popt[4]
#-------------------------------------------------------
def centroid(npkd, *arg, **kwarg):
    if npkd.dim == 1:
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', message='.*ovariance.*')
            centroid1d(npkd, *arg, **kwarg)
    elif npkd.dim == 2:
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', message='.*ovariance.*')
            centroid2d(npkd, *arg, **kwarg)
    else:
        raise Exception("Centroid yet to be done")
    return npkd
#-------------------------------------------------------
def display_peaks(npkd, peak_label=False, peak_mode="marker", zoom=None, show=False, color=None, markersize=6, figure=None):
    """
    display the content of the peak list, 
    peak_mode is either "marker" (default) or "bar" (1D only)
    zoom is in current unit.
    """
    if npkd.dim == 1:
        z1, z2 = parsezoom(npkd, zoom)
        if npkd.axis1.itype == 0:  # if real
            ff1 = npkd.axis1.itoc
        else:
            ff1 = lambda x : npkd.axis1.itoc(2*x)
        return npkd.peaks.display( peak_label=peak_label, peak_mode=peak_mode, zoom=(z1,z2), show=show, f=ff1, color=color, markersize=markersize, figure=figure)
    elif npkd.dim == 2:
        z1lo, z1up, z2lo, z2up = parsezoom(npkd, zoom)
        if npkd.axis1.itype == 0:  # if real
            ff1 = npkd.axis1.itoc
        else:
            ff1 = lambda x : npkd.axis1.itoc(2*x)
        if npkd.axis2.itype == 0:  # if real
            ff2 = npkd.axis2.itoc
        else:
            ff2 = lambda x : npkd.axis2.itoc(2*x)
        return npkd.peaks.display( peak_label=peak_label, zoom=((z1lo,z1up),(z2lo,z2up)), show=show, f1=ff1, f2=ff2, color=color, markersize=markersize, figure=figure)
    else:
        raise Exception("to be done")
#-------------------------------------------------------
def report_peaks(npkd, file=None, format=None):
    """
    print the content of the peak list, using the current unit
    
    if file should be an already opened writable file stream. 
    if None, output will go to stdout
    
    for documentation, check Peak1D.report() and Peak2D.report()
    """
    if npkd.dim == 1:
        if npkd.axis1.itype == 0:  # if real
            ff1 = npkd.axis1.itoc
        else:
            ff1 = lambda x : npkd.axis1.itoc(2*x)
        return npkd.peaks.report(f=ff1, file=file, format=format)
    elif npkd.dim == 2:
        if npkd.axis1.itype == 0:  # if real
            ff1 = npkd.axis1.itoc
        else:
            ff1 = lambda x : npkd.axis1.itoc(2*x)
        if npkd.axis2.itype == 0:  # if real
            ff2 = npkd.axis2.itoc
        else:
            ff2 = lambda x : npkd.axis2.itoc(2*x)
        return npkd.peaks.report( f1=ff1, f2=ff2, file=file, format=format)
    else:
        raise Exception("to be done")

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

NPKData_plugin("pp", peakpick)
NPKData_plugin("peakpick", peakpick)
NPKData_plugin("centroid", centroid)
NPKData_plugin("report_peaks", report_peaks)
NPKData_plugin("display_peaks", display_peaks)
NPKData_plugin("peaks2d", peaks2d)
