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

def findnoiselevel(fid, nbseg=10):
    """
    Routine for determining the noise level
    cut the data in segments then make the standard deviation 
    for each segment compare the std dev with the average one and eliminates 
    segments ginving a std dev above the mean
    finally makes the average on the sorted segments. 
    nbseg=10   nb of segment to cut the spectrum
    """

    less = len(fid)%nbseg     # rest of division of length of data by nb of segment
    restpeaks = fid[:-less]   # remove the points that avoid to divide correctly the data in segment of same size.
    print restpeaks.size
    print "less",less
    newlist = np.array(np.hsplit(restpeaks,nbseg))    #Cutting in segments
    noisestd = newlist.std(axis=1).min()
    noisemean = restpeaks.mean()
    print "noisestd ",noisestd
    print "noisemean ",noisemean
    noiselev = noisemean + noisestd
    return noiselev

def zerfapod(self, data, zerofill):
    '''
    Zerofilling and apodisation
    Zerofilling on self.refresol
    '''
    print "self.refresol ",self.refresol
    if zerofill :
        spec = data.apod_sin(maxi = 0.5).chsize(self.refresol).rfft().modulus().buffer
    else :
        spec = data.apod_sin(maxi = 0.5).rfft().modulus().buffer
    print "spec.size",spec.size
    return spec