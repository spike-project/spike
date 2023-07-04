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
import warnings
import unittest 
from collections import UserList
import numpy as np
from scipy.optimize import curve_fit
###
from spike import NPKError
from spike.NPKData import NPKData_plugin, _NPKData, flatten, parsezoom
from spike.util.counter import timeit
#from spike.util.signal_tools import findnoiselevel, findnoiselevel_2D

debug = 0
NbMaxDisplayPeaks = 1000      # maximum number of peaks to display at once
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
    def __init__(self, Id, label, intens, pos, pos_err=0.0, width=0.0, width_err=0.0 ):
        super(Peak1D, self).__init__(Id, label, intens)
        self.pos = float(pos)
        self.pos_err = float(pos_err)
        self.width = float(width)
        self.width_err = float(width_err)
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

class PeakList(UserList):
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
            raise ValueError("unrecognized keyword args: %s" % extras)
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
    # def __getitem__(self,i):
    #     #print('getitem')
    #     return type(self)(list.__getitem__(self,i))
def peak_aggreg(pklist, distance, maxdist=None, method='max'):
    """
    aggregates 1D peaks in peaklist if peaks are closer than a given distance in pixel
    distance : if two consecutive peaks are less than distance (in points),  they are aggregated
    if maxdist is not None, this is the maximal distance to the largest peak 
    method is either 'max' or 'mean'

    """
    if len(pklist) == 0:
        return pklist
    pkl = sorted(pklist, key = lambda p: p.pos)   # first, sort the list in position
    newlist = Peak1DList()
    prev = pkl[0]
    ipk = 1                    # runs on peak index
    maxgrp = (prev.intens, 0)  # maintains maxval and index of max peak
    newgrp = [0,]              # current group
    if maxdist is None:
        dmax = 10*distance
    else:
        dmax = maxdist
    while ipk < len(pkl):
        current = pkl[ipk]
        # print (current.pos, int(current.pos-prev.pos), distance)
        if (abs(current.pos-prev.pos) <= distance) and   \
            (abs(pkl[maxgrp[1]].pos-current.pos) <= dmax):  # aggregates
            # print(ipk, current.pos)
            newgrp.append(ipk)                     # add to list
            if current.intens>maxgrp[0]:           # check summit
                maxgrp = (current.intens, ipk)
        else:       # stop aggregating
            if len(newgrp) == 1:
                inewpk = newgrp[0]   # isolated peak
            else:
                inewpk =  maxgrp[1]  # larger peak of the aggregate
                if method == 'max':
                    pass   # keep position of max peak
                elif method == 'mean':   # put position in middle of the group
                    pkl[inewpk].pos = 0.5*(pkl[newgrp[0]].pos + pkl[newgrp[-1]].pos)
                else:
                    raise Exception('wrong method, use either "mean" or "max"')
            newlist.append( pkl[inewpk] )
            newgrp = [ipk]
            maxgrp = (current.intens, ipk)
        ipk += 1
        prev = current

    # finished - left-overs
    if len(newgrp) == 1:
        inewpk = newgrp[0]   # isolated peak
    else:
        inewpk =  maxgrp[1]  # larger peak of the aggregate
        if method == 'mean':   # put position in middle of the group
            pkl[inewpk].pos = 0.5*(pkl[newgrp[0]].pos + pkl[newgrp[-1]].pos)
    newlist.append( pkl[inewpk] )
    return newlist 
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
    def report(self, f=_identity, file=None, format=None, NbMaxPeaks=NbMaxDisplayPeaks):
        """
        print the peak list
        f is a function used to transform respectively the coordinates
        indentity function is default,
        for instance you can use something like
        d.peaks.report(f=d.axis1.itop)    to get ppm values on a NMR dataset

        check documentation for Peak1D.report() for details on output format
        """
        if NbMaxPeaks<len(self):
            print("# %d in Peak list - reporting the %d largest"%(len(self), NbMaxPeaks), file=file)
            indices = list(np.argpartition(self.intens,len(self)-NbMaxPeaks)[-NbMaxPeaks:])
            for i,pk in enumerate(self):
                if i in indices:
                    print(pk.report(f=f, format=format), file=file)
        else:
            print ("# %d in Peak list"%len(self), file=file)
            for pk in self:
                print(pk.report(f=f, format=format), file=file)
    def _report(self, f=_identity, file=None):
        """return full report for 1D peak list
        list is : Id, label, pos, intens, width,  pos_err, intens_err, width_err
        """
        lst = [pk._report(f=f) for pk in self ]
        return "\n".join(lst)
    def pkadd(self, pkl2):
        """adds the 1D peaklist pkl2 to self  - (no check on pkl2)"""
        if len(self)>0:
            maxId = max((pk.Id for pk in self))
            inc = maxId+1
        else:
            inc = 0
        for pk2 in pkl2:
            self.append( Peak1D(pk2.Id+inc, pk2.label, pk2.intens, pk2.pos, pk2.pos_err, pk2.width, pk2.width_err) )

    def peakaggreg(self, distance, maxdist=None, method='max'):
        """
        aggregates 1D peaks in peaklist if peaks are closer than a given distance in pixel
        check peak_aggreg() for detailed doc
        """
        self.data = peak_aggreg(self.data, distance=distance, maxdist=maxdist, method=method)

    def display(self, peak_label=False, peak_mode="marker", zoom=None, show=False, f=_identity, color='red',
            marker='x', markersize=6, figure=None, scale=1.0, NbMaxPeaks=NbMaxDisplayPeaks,
            markerdict=None, labeldict=None):
        """
        displays 1D peaks
        zoom is in index
        peak_mode is either "marker"  "bar" or "None" (to be used to just add labels)
        NbMaxPeaks is the maximum number of peaks to display in the zoom window (show only the largest)
        f() should be a function which converts from points to current display scale - typically npk.axis1.itoc
        """
        if len(self) == 0:
            return  # nothing to display
        from spike.Display import testplot
        plot = testplot.plot()
        if figure is None:
            fig = plot.subplot(111)
        else:
            fig = figure
        # create and filter list
        if zoom:
            z0=zoom[0]
            z1=zoom[1]
            pkl = [i for i,p in enumerate(self) if p.pos>=z0 and p.pos<=z1]  # index of peaks on zoom window
        else:
            pkl = range(len(self))
        if peak_mode == "marker" and len(self)>0:   # in marker mode, removes peaks too high
            mmax = max(self.intens)/scale
            pkl = list(filter(lambda i:self.intens[i]<=mmax, pkl))
        if len(pkl)>NbMaxPeaks:     # too many to display
            pkl.sort(reverse=True, key=lambda i: self[i].intens)
            pkl = pkl[:NbMaxPeaks]
        # create arg for display
        # default
        mark = {'linestyle':'',
                'marker':'x',
                'markersize':markersize,
                'color':color}
        label = {'color':color,
                 'fontsize':7,
                 'rotation':40}
        # then update with args
        if markerdict is None: markerdict = {}
        mark.update(markerdict)
        if labeldict is None: labeldict = {}
        label.update(labeldict)
        # now display
        if peak_mode == "marker":
            fig.plot(f(self.pos[pkl]), self.intens[pkl], **mark)
        elif peak_mode == "bar":
            for i in pkl:
                p = self[i]
                fig.plot( [f(p.pos),f(p.pos)], [0,p.intens], '-', color=color)
        elif peak_mode != "None":                 # no plot - 
            raise Exception("wrong peak_mode")
        if peak_label:
            for i in pkl:
                p = self[i]
                fig.annotate(p.label,(f(p.pos), p.intens),
                    xycoords='data', xytext=(0, 10), textcoords='offset points',
                    arrowprops=dict(arrowstyle='-'), horizontalalignment='left', verticalalignment='bottom', **label)
        if show: plot.show()
    # def merge(self, pklist, tolerance=1.0):
    #     """
    #     merge self and additional pklist by removing redundant peaks (closer than tolerance - in points)
    #     """
    #     self = peak_aggreg(self+pklist, distance=tolerance)  # that simple !
    def pos2label(self):
        "use pos in current unit, using converion f and set it as label for each peak"
        f = self.source.axis1.itoc
        for pk in self:
            pk.label = "%.4f"%(f(pk.pos),)
    def from_csv(self, filename):
        """import peak list from csv file,
        as would be created by pk2pandas.to_csv() coded in ppm, full=False
        """
        import pandas as pd
        df = pd.read_csv( filename, index_col='Id', dtype={"Label":"str"})
        self.clear()
        for i,v in enumerate(df.to_numpy()):
            self.append(Peak1D( Id=df.index[i],
                                label=v[0],
                                intens=v[2],
                                pos=self.source.axis1.itoix(self.source.axis1.ptoi(v[1])),
                                width=v[3])
                        )
        return self

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
    def report(self, f1=_identity, f2=_identity, file=None, format=None, NbMaxPeaks=NbMaxDisplayPeaks):
        """
        print the peak list
        f1, f2 are two functions used to transform respectively the coordinates in F1 and F2
        indentity function is default,
        for instance you can use something like
        d.peaks.export(f1=s.axis1.itop, f2=s.axis2.itop)    to get ppm values on a NMR dataset
        the file keyword allows to redirect the output to a file object 
        
        check documentation for Peak2D.report() for details on output format
        """
        if NbMaxPeaks<len(self):
            print("# %d in Peak list - reporting the %d largest"%(len(self), NbMaxPeaks), file=file)
            indices = list(np.argpartition(self.intens,len(self)-NbMaxPeaks)[-NbMaxPeaks:])
            for i,pk in enumerate(self):
                if i in indices:
                    print(pk.report(f1=f1, f2=f2, format=format), file=file)
        else:
            print ("# %d in Peak list"%len(self), file=file)
            for pk in self:
                print(pk.report(f1=f1, f2=f2, format=format), file=file)
    def _report(self, f1=_identity, f2=_identity, file=None):
        """return full report for 2D peak list
        list is : Id, label, posF1, posF2, intens, widthF1, widthF2, posF1_err, posF2_err, intens_err, widthF1_err, widthF2_err
        """
        lst = [pk._report(f1=f1, f2=f2) for pk in self ]
        return "\n".join(lst)
    def display(self, axis = None, peak_label=False, zoom=None, show=False, f1=_identity, f2=_identity, color=None,
            markersize=6, figure=None, NbMaxPeaks=NbMaxDisplayPeaks,
            markerdict=None, labeldict=None):
        """
        displays 2D peak list
        zoom is in index
        f1 and f2 should be functions which convert from points to current display scale - typically npk.axis1.itoc npk.axis2.itoc
        """
        import spike.Display.testplot as testplot
        plot = testplot.plot()
        if figure is None:
            fig = plot.subplot(111)
        else:
            fig = figure
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
        # create arg for display
        # default
        mark = {'markersize':markersize,
                'color':color}
        label = {'color':color,
                 'fontsize':7,
                 'rotation':40}
        # then update with args
        if markerdict is None: markerdict = {}
        mark.update(markerdict)
        if labeldict is None: labeldict = {}
        label.update(labeldict)
        if axis is None:
            plF1 = self.posF1 # these are ndarray !
            plF2 = self.posF2
            fig.plot(f2(plF2[pk]), f1(plF1[pk]), "x", **mark)
            if peak_label:
                for p in pk:
                    plp = self[p]
                    fig.text(1.01*plp.posF2, 1.01*plp.posF1, plp.label, **label)
        else:
            raise Exception("to be done")
        if show: fig.show()

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
def peakpick(npkd, threshold = None, zoom = None, autothresh=3.0, verbose=False):
    """
    performs a peak picking of the current experiment
    threshold is the level above which peaks are picked
        None (default) means that autothresh*(noise level of dataset) will be used - using d.robust_stats() as proxy for noise-level
    zoom defines the region on which detection is made
        zoom is in currentunit (same syntax as in display)
        None means the whole data
    """
    if threshold is None:
        mu, sigma = npkd.robust_stats()
        threshold = autothresh*sigma
        if mu >0:
            threshold += mu
    if npkd.dim == 1:
        listpkF1, listint = peaks1d(npkd, threshold=threshold, zoom=zoom)
                            #     Id, label, intens, pos
        pkl = Peak1DList( ( Peak1D(i, str(i), intens, pos) \
            for i, pos, intens in zip( range(len(listint)), list(listpkF1), list(listint) ) ), \
                    threshold=threshold, source=npkd )
        # use pos as labels - 
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
    # print("########## in peaks2d ")
    # print("zoom ", zoom)
    z1lo, z1up, z2lo, z2up = parsezoom(npkd, zoom)
    
    # print("z1lo, z1up, z2lo, z2up ", z1lo, z1up, z2lo, z2up)
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
    listint = npkd.get_buffer()[listpkF1].real
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
def centroid1d(npkd, npoints=3, reset_label=True, cure_outliers=True):
    """
    from peak lists determined with peak()
    realize a centroid fit of the peak summit and width,
    will use npoints values around center  (npoints has to be odd)
    computes Full width at half maximum
    updates in data peak list
    reset_label when True (default) reset the labels of FTMS datasets
    cure_outliers : restore peaks with pathological parameters
    TODO : update uncertainties
    """
    from scipy import polyfit
    npkd.check1D()
    noff = (int(npoints)-1)/2
    if (2*noff+1 != npoints) or (npoints<3):
        raise NPKError("npoints must odd and >2 ",data=npkd)
    buff = npkd.get_buffer().real
    for pk in npkd.peaks:
        xdata = np.arange(int(pk.pos-noff), int(pk.pos+noff+1))
        ydata = buff[xdata]
        try:
            popt, pcov = curve_fit(center, xdata, ydata, p0=[pk.pos, pk.intens, 1.0] ) # fit
        except RuntimeError:
            print ( "peak %d (id %s) centroid could not be fitted"%(pk.Id, pk.label) )
        if cure_outliers and (popt[0]>npkd.size1 or popt[0]<0):  # peak is outside
            continue
        else:
            pk.pos = popt[0]
            pk.intens = popt[1]
            pk.width = np.sqrt(2.0)*popt[2]
            errors = np.sqrt(np.diag(pcov))
            pk.pos_err = errors[0]
            pk.intens_err = errors[1]
            pk.width_err = np.sqrt(2.0)*errors[2]

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
        raise NPKError("npoints must be odd and >2 ",data=npkd)
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
        pk.widthF1 = np.sqrt(2.0)*popt[3]
        pk.widthF2 = np.sqrt(2.0)*popt[4]
        errors = np.sqrt(np.diag(pcov))
#        print(errors)
        pk.posF1_err = errors[0]
        pk.posF2_err = errors[1]
        pk.intens_err = errors[2]
        pk.widthF1_err = np.sqrt(2.0)*errors[3]
        pk.widthF2_err = np.sqrt(2.0)*errors[4]
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
def display_peaks(npkd, peak_label=False, peak_mode="marker", zoom=None, show=False, color=None,
markersize=6, figure=None, scale=1.0, NbMaxPeaks=NbMaxDisplayPeaks, markerdict=None, labeldict=None):
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
        return npkd.peaks.display( peak_label=peak_label, peak_mode=peak_mode, zoom=(z1,z2), show=show, f=ff1, color=color,
                markersize=markersize, figure=figure, scale=scale, NbMaxPeaks=NbMaxPeaks, markerdict=markerdict, labeldict=labeldict)
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
        return npkd.peaks.display( peak_label=peak_label, zoom=((z1lo,z1up),(z2lo,z2up)), show=show, f1=ff1, f2=ff2, color=color,
               markersize=markersize, figure=figure, NbMaxPeaks=NbMaxPeaks, markerdict=markerdict, labeldict=labeldict)
    else:
        raise Exception("to be done")
#-------------------------------------------------------
def report_peaks(npkd, file=None, format=None, NbMaxPeaks=NbMaxDisplayPeaks):
    """
    print the content of the peak list, using the current unit
    
    if file should be an already opened writable file stream. 
    if None, output will go to stdout
    
    for documentation, check Peak1D.report() and Peak2D.report()
       (or d.peaks[0].report() documentation)
    """
    if npkd.dim == 1:
        if npkd.axis1.itype == 0:  # if real
            ff1 = npkd.axis1.itoc
        else:
            ff1 = lambda x : npkd.axis1.itoc(2*x)
        return npkd.peaks.report(f=ff1, file=file, format=format, NbMaxPeaks=NbMaxPeaks)
    elif npkd.dim == 2:
        if npkd.axis1.itype == 0:  # if real
            ff1 = npkd.axis1.itoc
        else:
            ff1 = lambda x : npkd.axis1.itoc(2*x)
        if npkd.axis2.itype == 0:  # if real
            ff2 = npkd.axis2.itoc
        else:
            ff2 = lambda x : npkd.axis2.itoc(2*x)
        return npkd.peaks.report( f1=ff1, f2=ff2, file=file, format=format, NbMaxPeaks=NbMaxPeaks)
    else:
        raise Exception("to be done")
#-------------------------------------------------------
def pk2pandas(npkd, **kw):
    """export extract of current peak list to pandas Dataframe - in current unit
    if full is False (default), the uncertainties are not listed
    uses nmr or ms version depending on data_type
    """
    import spike.FTMS
    if isinstance(npkd, spike.FTMS.FTMSData):
        return pk2pandas_ms(npkd, **kw)
    else:
        return pk2pandas_nmr(npkd, **kw)
def pk2pandas_ms(npkd, full=False):
    "export extract of current peak list to pandas Dataframe for MS datasets"
    import pandas as pd
    if npkd.dim == 1:
        width = np.array([pk.width for pk in npkd.peaks])   # width have to be computed in m/z
        width_array = 0.5*abs(npkd.axis1.itoc(npkd.peaks.pos+width) - npkd.axis1.itoc(npkd.peaks.pos-width))
        width_array = np.where(width_array==0, np.NaN, width_array)
        if full:
            err = np.array([pk.pos_err for pk in npkd.peaks])
            pos_err_array = 0.5*abs(npkd.axis1.itoc(npkd.peaks.pos+err) - npkd.axis1.itoc(npkd.peaks.pos-err))
            P1 = pd.DataFrame({
                'Id':[ pk.Id for pk in npkd.peaks],
                'Label':npkd.peaks.label,
                'm/z':npkd.axis1.itoc(npkd.peaks.pos),
                u'Δm/z':width_array,
                'R':  np.around(npkd.axis1.itoc(npkd.peaks.pos)/width_array, 3),
                'Intensity':npkd.peaks.intens,
                u'ε_m/z': pos_err_array,
                u'ε_Intensity': [pk.intens_err for pk in npkd.peaks]
            }).set_index('Id')
        else:
            P1 = pd.DataFrame({
                'Id':[ pk.Id for pk in npkd.peaks],
                'Label':npkd.peaks.label,
                'm/z':npkd.axis1.itoc(npkd.peaks.pos),
                u'Δm/z':width_array,
                'R':  np.around(npkd.axis1.itoc(npkd.peaks.pos)/width_array, 3),
                'Intensity':npkd.peaks.intens
            }).set_index('Id')

    elif npkd.dim == 2:
        # width in m/z
        width1 = np.array([pk.widthF1 for pk in npkd.peaks])   # width have to be computed in m/z
        width1_array = 0.5*abs(npkd.axis1.itoc(npkd.peaks.posF1+width1)
                             - npkd.axis1.itoc(npkd.peaks.posF1-width1))
        width2 = np.array([pk.widthF2 for pk in npkd.peaks])   # width have to be computed in m/z
        width2_array = 0.5*abs(npkd.axis2.itoc(npkd.peaks.posF2+width2)
                             - npkd.axis2.itoc(npkd.peaks.posF2-width2))

        # then R
        R1 =  npkd.axis1.itoc(npkd.peaks.posF1)/width1_array
        R2 =  npkd.axis2.itoc(npkd.peaks.posF2)/width2_array
        if full:
            # then for position_error to current unit
            err1 = np.array([pk.posF1_err for pk in npkd.peaks])
            err2 = np.array([pk.posF2_err for pk in npkd.peaks])
            pos1_err_array = 0.5*abs(npkd.axis1.itoc(npkd.peaks.posF1+err1) - npkd.axis1.itoc(npkd.peaks.posF1-err1))
            pos2_err_array = 0.5*abs(npkd.axis2.itoc(npkd.peaks.posF2+err2) - npkd.axis2.itoc(npkd.peaks.posF2-err2))
            P1 = pd.DataFrame({
                'Id':[ pk.Id for pk in npkd.peaks],
                'Label':npkd.peaks.label,
                'm/z F1':npkd.axis1.itoc(npkd.peaks.posF1),
                'm/z F2':npkd.axis2.itoc(npkd.peaks.posF2),
                'Intensity':npkd.peaks.intens,
                u'Δm/z F1':width1_array,
                u'Δm/z F2':width2_array,
                'R1':  np.around(R1, 3),
                'R2':  np.around(R2, 3),
                'R *1E6':  np.around(R1*R2/1E6, 3),
                u'ε_m/z F1': pos1_err_array,
                u'ε_m/z F2': pos2_err_array,
                u'ε_Intensity': [pk.intens_err for pk in npkd.peaks]
            }).set_index('Id')
        else:
            P1 = pd.DataFrame({
                'Id':[ pk.Id for pk in npkd.peaks],
                'Label':npkd.peaks.label,
                'm/z F1':npkd.axis1.itoc(npkd.peaks.posF1),
                'm/z F2':npkd.axis2.itoc(npkd.peaks.posF2),
                'Intensity':npkd.peaks.intens,
                u'Δm/z F1':width1_array,
                u'Δm/z F2':width2_array,
                'R1':  np.around(R1, 3),
                'R2':  np.around(R2, 3),
                'R *1E6':  np.around(R1*R2/1E6, 3)
            }).set_index('Id')
    else:
        raise Exception('Not implemented in 3D yet')
    return P1
def pk2pandas_nmr(npkd, full=False, unit='current'):   
    """export extract of current peak list to pandas Dataframe for NMR datasets, 
    unit is either "current", "points", "ppm" or "Hz" .
    """
    import pandas as pd          # LAZY IMPORT
    if npkd.dim == 1:
        # ŧhis part could be buggy in certain strange cases
        if unit == "current":
            conv = npkd.axis1.ixtoc
        if unit == "points":
            conv = lambda x: x
        if unit == "ppm":
            conv = lambda x: npkd.axis1.ixtoi(npkd.axis1.itop(x))
        if unit == "Hz":
            conv = lambda x: npkd.axis1.ixtoi(npkd.axis1.itoh(x))
        w_itohz = 2*npkd.axis1.specwidth/npkd.cpxsize1   # for width: points to Hz
        # then for position_error to given unit
        err = np.array([pk.pos_err for pk in npkd.peaks])
        pos_err_array = conv(npkd.peaks.pos+err) - conv(npkd.peaks.pos-err)
        pos_err_array = 0.5*np.abs(pos_err_array)
        if full:
            P1 = pd.DataFrame({
                'Id':[ pk.Id for pk in npkd.peaks],
                'Label':npkd.peaks.label,
                'Position':conv(npkd.peaks.pos),
                'Intensity':npkd.peaks.intens,
                'Width':[ pk.width*w_itohz for pk in npkd.peaks],   # width are in Hz
                'Position_err': pos_err_array,
                'Intensity_err': [pk.intens_err for pk in npkd.peaks],
                'Width_err': [pk.width_err*w_itohz for pk in npkd.peaks]   # width are in Hz
            }).set_index('Id')
        else:
            P1 = pd.DataFrame({
                'Id':[ pk.Id for pk in npkd.peaks],
                'Label':npkd.peaks.label,
                'Position':conv(npkd.peaks.pos),
                'Intensity':npkd.peaks.intens,
                'Width':[ pk.width*w_itohz for pk in npkd.peaks]   # width are in Hz
            }).set_index('Id')
    elif npkd.dim == 2:
        w_itohz1 = 2*npkd.axis1.specwidth/npkd.cpxsize1   # for width: points to Hz
        w_itohz2 = 2*npkd.axis2.specwidth/npkd.cpxsize2
        # then for position_error to current unit
        err1 = np.array([pk.posF1_err for pk in npkd.peaks])
        err2 = np.array([pk.posF2_err for pk in npkd.peaks])
        pos1_err_array = npkd.axis1.itoc(npkd.peaks.posF1+err1) - npkd.axis1.itoc(npkd.peaks.posF1-err1)
        pos2_err_array = npkd.axis2.itoc(npkd.peaks.posF2+err2) - npkd.axis2.itoc(npkd.peaks.posF2-err2)
        pos1_err_array = 0.5*np.abs(pos1_err_array)
        pos2_err_array = 0.5*np.abs(pos2_err_array)
        if full:
            P1 = pd.DataFrame({
                'Id':[ pk.Id for pk in npkd.peaks],
                'Label':npkd.peaks.label,
                'Position F1':npkd.axis1.itoc(npkd.peaks.posF1),
                'Position F2':npkd.axis2.itoc(npkd.peaks.posF2),
                'Intensity':npkd.peaks.intens,
                'Width F1':[ pk.widthF1*w_itohz1 for pk in npkd.peaks],   # width are in Hz
                'Width F2':[ pk.widthF2*w_itohz2 for pk in npkd.peaks],
                'Position_err F1': pos1_err_array,
                'Position_err F2': pos2_err_array,
                'Intensity_err': [pk.intens_err for pk in npkd.peaks],
                'Width_err F1': [pk.widthF1_err*w_itohz1 for pk in npkd.peaks],   # width are in Hz
                'Width_err F2': [pk.widthF2_err*w_itohz1 for pk in npkd.peaks]
            }).set_index('Id')
        else:
            P1 = pd.DataFrame({
                'Id':[ pk.Id for pk in npkd.peaks],
                'Label':npkd.peaks.label,
                'Position F1':npkd.axis1.itoc(npkd.peaks.posF1),
                'Position F2':npkd.axis2.itoc(npkd.peaks.posF2),
                'Intensity':npkd.peaks.intens,
                'Width F1':[ pk.widthF1*w_itohz1 for pk in npkd.peaks],   # width are in Hz
                'Width F2':[ pk.widthF2*w_itohz2 for pk in npkd.peaks]
            }).set_index('Id')
    else:
        raise Exception('Not implemented in 3D yet')
    return P1

class PeakTests(unittest.TestCase):
    def test_peaks2d(self):
        "test 2D peak picker"
        print(self.test_peaks2d.__doc__)
        M=np.zeros((30, 30))
        M[5,7] = 20
        M[10,12] = 20
        d = _NPKData(buffer = M)
        d.pp(threshold = 10) # 3*d.std was just right - but does not work anymore with robust_stats() !
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
        d = _NPKData(buffer = M)
        d.pp(threshold=3)
        self.assertEqual(list(d.peaks.pos) , [ 5.0, 7.0, 10.0, 15.0])
        self.assertEqual(list(d.peaks.intens) , [ 20.0, 8.0, 20.0, 11.0])
        d.peaks.report(NbMaxPeaks=2)
    def test_peaksaggreg1d(self):
        "test 1D peak aggregator"
        M = np.zeros((300))
        j = 10
        for i in range(2,20):   # makes 18 peaks with increasing distance
            M[j] = 10
            j += i
        d = _NPKData(buffer = M)
        d.pp(threshold=3)
        self.assertEqual(len(d.peaks), 18)
        self.assertEqual(d.peaks[1].pos, 12)
        # aggreg 2 first peaks
        d.peaks.peakaggreg(distance=2)   # method = max
        self.assertEqual(len(d.peaks), 17, msg="length")
        self.assertEqual(d.peaks[1].pos, 15)
        self.assertEqual(d.peaks[2].pos, 19, msg="pk1")
        # aggreg next using mean
        d.peaks.peakaggreg(distance=4, method="mean")
        self.assertEqual(len(d.peaks), 16, msg="length")
        self.assertEqual(d.peaks[1].pos, 17, msg="pk1")
        # test maxdist
        d.peaks.peakaggreg(distance=10, method="mean", maxdist=20)
        self.assertEqual(len(d.peaks), 11, "length")
        self.assertEqual(d.peaks[1].pos, 45.5, "pk1")
    def test_center1d(self):
        x = np.arange(5)
        y = center(x, 2.2, 10.0, 1.2)
        # y == [ -2.36111111e+01  -4.44089210e-15   9.72222222e+00   5.55555556e+00 -1.25000000e+01]
        self.assertAlmostEqual(y[2], 9.72222222)
        d = _NPKData(buffer = np.maximum(y,0.0))
        d.peaks = Peak1DList(source=d)
        d.peaks.append(Peak1D(0, "0", 9.7, 2 ))
        d.peaks[-1].width = 1.0
        d.centroid(npoints=3)
        self.assertAlmostEqual(d.peaks[0].pos,2.2)
        self.assertAlmostEqual(d.peaks[0].intens,10.0)
        self.assertAlmostEqual(d.peaks[0].width,1.2*np.sqrt(2))
        d.peaks.report(NbMaxPeaks=10)
    def test_center2d(self):
        M=np.zeros((20, 20))
        # add one peak at (F1,F2) 5.3, 7.9 with widthes (5.0,1.3) 
        for y in range(1,10):
            for x in range(6,11):
                #                 yo, x0 intens, widthy, widthx
                M[y,x] = center2d(np.array([y,x]), 5.3, 7.9, 20.0, 5.0, 1.3)
        #print (M[1:10,6:11])
        self.assertAlmostEqual(M[2,7],5.87777515)
        d = _NPKData(buffer = np.maximum(M,0.0))
        d.peaks = Peak2DList(source=d)
        # self, Id, label, intens, posF1, posF2 
        d.peaks.append(Peak2D(0, "0", 18.0, 5, 8 ))
        d.centroid(npoints_F1=5)
        self.assertAlmostEqual(d.peaks[0].posF1,5.3)
        self.assertAlmostEqual(d.peaks[0].posF2,7.9)
        self.assertAlmostEqual(d.peaks[0].intens,20.0)
        self.assertAlmostEqual(d.peaks[0].widthF1,5.0*np.sqrt(2))
        self.assertAlmostEqual(d.peaks[0].widthF2,1.3*np.sqrt(2))

NPKData_plugin("pp", peakpick)
NPKData_plugin("peakpick", peakpick)
NPKData_plugin("centroid", centroid)
NPKData_plugin("report_peaks", report_peaks)
NPKData_plugin("display_peaks", display_peaks)
NPKData_plugin("peaks2d", peaks2d)
NPKData_plugin("pk2pandas", pk2pandas)
