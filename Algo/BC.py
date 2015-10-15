#!/usr/bin/env python 
# encoding: utf-8

"""
This module contains several utilities for baseline correction of spectra

Created by Marc-Andre on 2015-03-26.

Modification by Lionel 2015-07-10

"""

from __future__ import print_function
from scipy.optimize import minimize
import numpy as np
import unittest
import multiprocessing as mp
from scipy import interpolate

def poly(x,coeff):
    "computes the polynomial over x with coeff"
    lcoeff = list(coeff)
    y = np.zeros_like(x)
    y += lcoeff.pop()    # a
    while lcoeff:
        y *= x          # ax + b
        y += lcoeff.pop()
    return y

def fitpolyL1(x, y, degree=2, power=1, method="Powell"):
    "fit with L1 norm a polynome to a function y over x, returns coefficients"
    coeff0 = [0]*(degree+1)
    #pmin = lambda c, x, y: np.sum(np.abs(y-poly(x,c)) )
    pmin = lambda c, x, y: np.sum(np.power(np.abs(y-poly(x,c)), power) )
    res = minimize(pmin, coeff0, args=(x,y), method=method)
    return res.x

def bcL1_init(y, degree=2, power=1, method="Powell"):
    "compute a baseline on y using fitpolyL1"
    x = np.arange(1.0*y.size)
    coeff = fitpolyL1(x, y, degree=degree, power=power, method=method)
    return poly(x,coeff)

def bcL1(args): #y, degree=2, power=1, method="Powell"
    # print(args)
    y, degree, power, method  = args
    "compute a baseline on y using fitpolyL1"
    x = np.arange(1.0*y.size)
    coeff = fitpolyL1(x, y, degree=degree, power=power, method=method)
    return poly(x,coeff)

def interp_bl1(bl, lsize, nbseg = 17, kind = "linear"):
    '''
    Takes the whole baseline, divides it in smaller chunks than the initial chunks. Then makes an interpolation.
    '''
    bl_reduced = np.concatenate((bl[::lsize/nbseg], bl[-1:]))
    x_reduced = np.concatenate((np.arange(bl.size)[::lsize/nbseg], np.array([bl.size-1]))) 
    bl_interp = interpolate.interp1d(x_reduced, bl_reduced, kind=kind) # interpolate.splrep(x_reduced, bl_reduced, s=0)
    bl_final = bl_interp(np.arange(bl.size))
    return bl_final


def baseline0(y, degree=2, power=1, method="Powell", chunksize=2000, nbcores=None):
    """
    compute a piece-wise baseline on y using fitpolyL1
    degree is the degree of the underlying polynome
    chunksize defines the size of the pieces
    a cosine roll-off is used to smooth out chunks junctions

    y - baseline(y) produces a baseline corrected spectrum
    """
    nchunk = y.size/chunksize
    if nchunk <2:
        bl = bcL1_init(y, degree=degree, power=power, method=method)
    else:
        lsize = y.size/nchunk
        recov = lsize/10  # recovering parts
        corr = np.linspace(0.0,1.0,2*recov)
        corr = np.sin( np.linspace(0,np.pi/2,2*recov) )**2  # cosine roll-off
        corrm1 = 1.0-corr
        bl = np.zeros_like(y)
        bl[0:lsize+recov] = bcL1_init(y[0:lsize+recov], degree=degree, power=power)
        i = 0 # if nchunk == 2 !
        ###
        if nbcores: # Parallelization
            p = mp.Pool(nbcores) # Multiprocessing
            args = iter([[y[i*lsize-recov:(i+1)*lsize+recov], degree, power, method] for i in range(1,nchunk-1)])
            res = p.imap(bcL1, args)
            for k, estimate in enumerate(res):
                s0 = slice((k+1)*lsize-recov,(k+2)*lsize+recov)
                srecov0 = slice((k+1)*lsize-recov,(k+1)*lsize+recov)
                srecov1 = slice((k+1)*lsize+recov,(k+2)*lsize+recov)
                tbl = estimate
                bl[srecov0] = bl[srecov0]*corrm1 + tbl[:2*recov]*corr
                bl[srecov1] = tbl[2*recov:]
                
            p.close() 
        else:
            for i in range(1,nchunk-1):
                tbl = bcL1_init(y[i*lsize-recov:(i+1)*lsize+recov], degree=degree, power=power)
            
                bl[i*lsize-recov:i*lsize+recov] = bl[i*lsize-recov:i*lsize+recov]*corrm1 + tbl[:2*recov]*corr
                bl[i*lsize+recov:(i+1)*lsize+recov] = tbl[2*recov:]
        i = i+1
        tbl = bcL1_init(y[i*lsize-recov:-1], degree=degree)
        bl[i*lsize-recov:i*lsize+recov] = bl[i*lsize-recov:i*lsize+recov]*corrm1 + tbl[:2*recov]*corr
        bl[i*lsize+recov:] = tbl[2*recov-1:]
    return bl

def baseline1(y, degree=2, power=1, method="Powell", chunksize=2000, nbcores=None):
    """
    compute a piece-wise baseline on y using fitpolyL1
    degree is the degree of the underlying polynome
    chunksize defines the size of the pieces

    y - baseline(y) produces a baseline corrected spectrum
    """
    nchunk = y.size/chunksize
    if nchunk <2:
        bl = bcL1(y, degree=degree, power=power, method=method)
    else: # Parallelize
        lsize = y.size/nchunk
        bl = np.zeros_like(y)
        p = mp.Pool(nbcores) # Multiprocessing
        args = iter([[y[i*lsize:(i+1)*lsize], degree, power, method] for i in range(nchunk)])
        res = p.imap(bcL1, args)
        for i, estimate in enumerate(res):
            bl[i*lsize:(i+1)*lsize] = estimate
        p.close()
    bl_final = interp_bl1(bl, lsize, kind = "linear") # Smoothing to avoid segmentation of the baseline.
    
    return bl_final

def correctbaseline(y, iterations=1, chunksize=100, firstpower=0.3,
                        secondpower=7, degree=2,  chunkratio=1.0,
                        interv_ignore = None, method="Powell",
                        nbcores= 10,
                        debug = False, choiceBL = 0):
    '''
    Find baseline by using low norm value and then high norm value to attract the baseline on the small values.
    iterations : number of iterations for convergence toward the small values. 
    chunksize : size of each chunk on which is done the minimization. Typically, it must be larger than the peaks. 
    firstdeg : degree used for the first minimization 
    degree : degree of the polynome used for approaching each signal chunk. 
    chunkratio : ratio for changing the chunksize inside main loop
    interv_ignore : ignore a given intervall in the spectrum (eg : avoids issues with water pick)
    method : Algorithm used for minimization on each chunk
    nbcores : number of cores used for minimizing in parallel on many chunks. 
    '''
    if choiceBL == 0:
        baseline = baseline0
    elif choiceBL == 1:
        baseline = baseline1
    if interv_ignore:
        ii = interv_ignore
        delta = ii[1]-ii[0]
        y[ii[0]:ii[1]] = y[ii[0]] + np.arange(delta)/float(delta)*(y[ii[1]]-y[ii[0]]) # linear interpolation on the intervall.
    
    bl = baseline(y, degree=degree, power=firstpower, chunksize = chunksize, nbcores=nbcores, method="Powell")
    bls = {'bl':[], 'blmin':[]}
    for i in range(iterations):
        blmin = np.minimum.reduce([bl, y])
        bl = baseline(blmin, degree=degree, power=secondpower, chunksize = int(chunksize*chunkratio), nbcores=nbcores, method=method)
        bls['bl'].append(bl)
        bls['blmin'].append(blmin)
    if debug:
        return bl, bls
    else:
        return bl

class BC_Tests(unittest.TestCase):
    def test_poly(self):
        "tests the poly function"
        p = poly(np.arange(10.0),(.1,.2,.3,.4))
        self.assertEqual(p[6], 98.5)
        self.assertAlmostEqual(sum(p), 905.5)
    def test_baseline(self):
        N = 100000
        x = np.linspace(0,10,N)
        y = np.sin(x/2) + 0.2*np.random.randn(N)
        b = baseline(y,chunksize=N/20)
        corr = y-b
        self.assertTrue(np.std(corr) < 0.21)
    def test_correctbaseline(self):
        N = 100000
        x = np.linspace(0,10,N)
        y = np.sin(x/2) + 0.2*np.random.randn(N)
        b = correctbaseline(y, iterations=10, chunksize=N/20)
        corr = y-b
        self.assertTrue(np.std(corr) < 0.25)
        
if __name__ == '__main__':
    unittest.main()

