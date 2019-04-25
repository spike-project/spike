#!/usr/bin/env python
# encoding: utf-8
"""Autmomatic phase correction for 1D NMR spectra

based on an earlier version from NPK

works by minimizing the negative part of the spectrum
WILL NOT WORK on positive/negative spectra (JMOD, W-LOGSY, etc.)

Created by DELSUC Marc-André on 2016-05-23.
Copyright (c) 2016 IGBMC. All rights reserved.
"""
import numpy as np
#import matplotlib.pyplot as plt

from spike.NPKData import NPKData_plugin
#from spike.Algo.BC import correctbaseline as cbl

def neg_wing(d, bcorr=False, inwater=False, apt=False):
    """ measure negative wing power of NPKData d 
    
    if bcorr == True, a baseline correction is applied
    if inwater == True, the 10% central zone is just zeroed
    if apt == False, computes the std() of the negative points  (distance to mean)
           == True, computes the sum(abs()) of all the points (l_1 norm)
    """
    dd = d.copy().real()
    if inwater:
        dd[int(0.45*d.size1):int(0.55*d.size1)] = 0.0
    if bcorr:   # complete baseline corr
        dd.bcorr(method="spline",  xpoints=4)

    data = dd.get_buffer()
    lendata = len(data)

    if not bcorr:       # simple linear corr
        wind = [int(0.05*lendata), int(0.95*lendata)]
        hwidth = int(0.01*lendata/2)
        y0 = data[wind[0]-hwidth:wind[0]+hwidth].mean()
        y1 = data[wind[1]-hwidth:wind[1]+hwidth].mean()
        bl1 = np.poly1d(np.polyfit(wind, [y0,y1], 1))
        data -= bl1(np.arange(lendata))

    if not apt:
        data[data>0.0] = 0.0    # set positive to 0.0
        val = data.std()        # and compute std() as l_2 norm
    else:
        val = np.sum( np.abs(data) )  # simply l_1 norm
    return val

def phase_pivot(d, p0, p1, pivot=0.5):
    """ three parameter phasing routine
        pivot = 0 is on left side
        pivot = 1 is on right side
        all intermediate values are possible
        returns actual (P0, P1)
    """
    lp0, lp1 = p0+(0.5-pivot)*p1, p1
    d.phase(lp0, lp1)
    return (lp0, lp1)

def apmin(d, first_order=True, inwater=False, baselinecorr=True, apt=False, debug=False):
    """automatic 1D phase correction
    phase by minimizing the negative wing of the 1D spectrum
    
    first_order = False   inhibit optimizing 1st order phase
    inwater = True   does not look to the central zone of the spectrum
    baselinecorr = True, an advanced baseline correction is applied on the final steps
    apt = True   (Attached proton test) performs the phasing on up-down spectra, such as APT / DEPT 13C spectra.
    
    performs a grid/simplex search on P0 first then on (P0 P1)    
    the dataset is returned phased and the values are stored in d.axis1.P0 and d.axis1.P1

    P1 is kept to 0 if first_order=False
    
    note that if baselinecorr is True
        - the algo becomes quite slow !
        - a simple linear baseline correction is always applied anyway anyhow

    adapted from NPK v1 
    MAD, may 2016
    """
    d.check1D()
    if d.axis1.itype != 1:
        raise Exception("On complex data only")

# find largest point and use as Pivot
    im = float(abs(d.get_buffer()).argmax())
    pivot = (im/d.cpxsize1)
    if debug: print ("Pivot:",pivot)
# initialize
    valmin = neg_wing(d, bcorr=False, inwater=inwater, apt=apt)
    P0min, P1min = 0,0
    P0minnext, P1minnext = 0,0
    neval=0
    bcorr = baselinecorr #False
# first coarse
    P0step=10.0
    if first_order:
        P1step=40.0
    else:
        P1step=0.0
    while( abs(P0step)>=1.0 ):      # stops when increment is 1 degree phase
        moved=1
        while(moved):
            moved=0
            if first_order:
                PPlist = ((P0step,0),(-P0step,0),(0,P1step),(0,-P1step))
                # if abs(P1min)>360:
                #     print("1st order correction too large - reseting algo")
                #     P1min = 0.0
                #     valmin = np.inf
            else:
                PPlist = ((P0step,0),(-P0step,0))
            for (dP0,dP1) in PPlist:
                neval = neval+1
                dd = d.copy()
                phase_pivot(dd, P0min+dP0, P1min+dP1,pivot)
                pw = neg_wing(dd, bcorr=bcorr, inwater=inwater, apt=apt)
                if debug: print (" %.2f %.2f %g"%(P0min+dP0, P1min+dP1, pw))
                if (pw<valmin):
                    moved=1
                    valmin=pw
                    P0minnext = P0min+dP0
                    P1minnext = P1min+dP1
                    if (P0step*dP0 <0):     # try to remember variation direction
                        P0step = -P0step
                    if (P1step*dP1 <0):
                        P1step = -P1step
                    break
            P0min = P0minnext
            P1min = P1minnext
            if debug:
                dd = d.copy()
                P0,P1 = phase_pivot(dd, P0min, P1min, pivot)
                print ("*** P0 P1 :", P0, P1)
                color_sequence = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c',
                                  '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5',
                                  '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f',
                                  '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5']
                dd.display(new_fig=False, label="%.0f %.0f"%(P0,P1), color=color_sequence[neval%len(color_sequence)])
        P0step=P0step/2.0
        P1step=P1step/2.0
        if P0step < 5.0:
            bcorr = baselinecorr     # bcorr is expensive, so we make it only at the end if needed
            if debug: print ('bcorr = True')
    (P0,P1) = phase_pivot(d, P0min, P1min, pivot)
    if debug:
        print ("**FINAL** %.2f %.2f   in %d evaluations"%(P0, P1, neval))
    d.axis1.P0 = P0
    d.axis1.P1 = P1
    return d

NPKData_plugin("apmin", apmin)
