#!/usr/bin/env python 
# encoding: utf-8

"""
poulding.py

Apply the Poulding denoising algorithm

Poulding et al.
Removal of t(1) noise from metabolomic 2D (1)H-(13)C HSQC NMR spectra by Correlated Trace Denoising. J Magn Reson (2007) vol. 189 (2) pp. 190-9

Created by Marc-André on 2010-10-16.
Copyright (c) 2010 IGBMC. All rights reserved.
"""

from __future__ import print_function
import numpy as np
import math
import matplotlib.pyplot as plt
import NPKData as npkd

#----------------------------------
def dofilter(spec,threshold=2.921):
    """
    separate signal from noise, see 3.2 in the paper
    """
    def smooth(buf,t):
        "used by thresholding"
        return buf - (t**2)/buf
    # compute r2
    if spec.axis1.itype == 0:
        print("reel !")
        r2 = spec.buffer
    else:
        d = npkd.as_cpx(spec.buffer)
        r2 = np.real(np.sqrt(d*d.conj()))
    # compute threshold value
    med = np.median(r2)
    thresh = threshold*med
    # apply threshold and compute rpeak and rnoise
    rpeak = np.where(r2>thresh, smooth(r2,thresh), 0.0)    # set to zero or smooth
    rnoise = r2-rpeak
    # compute noise trace
    if spec.axis1.itype == 0:
        pnoise = rnoise.copy()
    else:
        pnoise = rnoise*(d/r2)  # d/r2 is phase only
    return pnoise
#----------------------------------
def corr(ph1, ph2):
    """
    computes complex correlation between complex vectors ph1 ad ph2  - eq(4)
    returns (x,y) as a complex pair BEWARE !
    """
    ph1c = ph1 - ph1.mean()
    ph2c = ph2 - ph2.mean()
    num = (ph1c*ph2c.conj()).sum()
    denom = math.sqrt( np.sum( ph1c*ph1c.conj() ) ) * math.sqrt( np.sum( ph2c*ph2c.conj() ) )
    c = num/denom
    return c
#    return complex( float(c*c.conj()), math.atan2(c.imag,c.real) )
#----------------------------------
def filter_matrix(data2D,clist):
    """
    given data2D, a 2D NPKData 
    builds a matrix (n x size1) containing the n "filtered" columns which list is found in clist
    """
    n = len(clist)
    mat = np.zeros((n,data2D.size1/2),dtype=complex)
    for i in xrange(len(clist)):
        #print i,clist[i]
        dd = data2D.col(clist[i])
        df = dofilter(dd)
        mat[i,:] = df
    return mat
#----------------------------------
def corr_matrix(data):
    """computes correlation matrix for data filtered matrix"""
    n = data.shape[0]
    mat = np.eye(n, dtype=complex)
    for i in xrange(n):
        for j in xrange(i+1,n):
            mat[i,j] = corr(data[i,:], data[j,:])
            mat[j,i] = mat[i,j].conj()
#            print i, j, mat[i,j], mat[j,i]
    return mat
#----------------------------------
def mask_matrix(corrmat, filtered, thresh=0.2):
    """
    compute noise mask from the complex correlation matrix, and the filtered data
    eqs 5 and 6
    """
    from operator import itemgetter
    def medabs(v):  # compute the median of the absolute value of the vector v
         return np.median( np.abs( v ) )
    n, p = filtered.shape
    # compute medians once
    med = []
    for i in xrange(n):
        med.append(  medabs( filtered[i,:] ) )   # computes before hand weights of eq 5
    cc = corrmat-np.eye(n)
    rcorr = np.sqrt( np.real(cc*cc.conj()) )
    anglecorr = np.arctan2( cc.real,  cc.imag)
    mat = np.zeros((n,p),dtype=complex)
    for i in xrange(n):     # for all entries
        dic = {}
        for j in xrange(n):     # copy columns
            dic[j] = rcorr[i,j]
        pkm = []
        # keep all entries such that rcorr>thres, and at least 3, accumulate into pkm
        for (pk,val) in sorted(dic.items(), key=itemgetter(1), reverse=True):     # sort on content
            if val<thresh:
                if len(pkm)>3:
                    break
            pkm.append(pk)
        print("**", i, pkm)
        for p in pkm:
            mat[i,:] +=  ( rcorr[i,p]*np.exp(1j*anglecorr[i,p])/med[p] ) * filtered[p,:]    # eq 5
        mat[i,:] *= ( med[i]/medabs(mat[i,:]) )                                  # eq 6
    return mat
#----------------------------------
def testcorr():
    c1 = npkd.as_cpx(np.random.rand(2000))
    c2 = npkd.as_cpx(np.random.rand(2000))
    c = corr(c1, c2)
    r =  math.sqrt(c.real**2 + c.imag**2)
    teta = math.atan2(c.real, c.imag)
    print(r, teta)
    c = corr( c1, c2*np.exp(-1j*teta) ) # ( math.cos(-teta)+1j*math.sin(-teta) ) )
    r =  math.sqrt(c.real**2 + c.imag**2)
    teta = math.atan2(c.real, c.imag)
    print(r, teta)
#----------------------------------
def test():
    #load
    # try:
    #     c1 = npkd.NPKData(name="col8200.gs1")
    #     c2 = npkd.NPKData(name="col8248.gs1")
    #     c3 = npkd.NPKData(name="col12104.gs1")
    d1 = npkd.NPKData(name="/DATA/Marie van/data_raw.gf2")  # load 2D
    d1.apod_sin(axis=2,maxi=0.5).rfft(axis=2)   # process
    d1.apod_sin(axis=1,maxi=0.5).rfft(axis=1)
    d1.extract(0, d1.size1, 2000, d1.size2)
    d2 = d1.copy()
    d2.modulus()
    proj = d2.proj(axis='F2', projtype="m")     # get proj
    proj.display( label="somme de tout")
    med = np.median(proj.buffer)    # and peak-pick
    proj.addbase(-med)
    proj.peak(threshold=0.05)
    proj.addbase(med)
    print("%d peaks"%len(proj.peaks))
    proj.display_peaks()
    # compute filtered spectra
    filtered = filter_matrix(d1,2*proj.peaks)
    plt.figure()
    plt.plot(filtered.T)
    # corr matrix
    corrmat = corr_matrix(filtered)
    rcorr = np.sqrt( np.real(corrmat*corrmat.conj()) )
    anglecorr = np.arctan2( corrmat.real,  corrmat.imag)
    plt.figure()
    plt.contour(rcorr-np.eye(len(proj.peaks)))
    # apply eq 7
    mask = mask_matrix(corrmat, filtered)
    for i in xrange(len(proj.peaks)):
        d1.col(2*proj.peaks[i]).display(label="avant")
        plt.plot(-mask[i,:],label="mask")
        d1.buffer[:,2*proj.peaks[i]] -= npkd.as_float( mask[i,:])
        d1.col(2*proj.peaks[i]).display(new_fig=False, label="apres")
    d1.modulus()
    d2.display(scale=2)
    d1.display(scale=2)
    
if __name__ == '__main__':
    test()
    plt.show()
