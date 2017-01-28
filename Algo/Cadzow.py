#!/usr/bin/env python 
# encoding: utf-8

"""
cadzow.py

Created by Marc-André on 2010-03-18.
Copyright (c) 2010 IGBMC. All rights reserved.
"""

from __future__ import print_function, division

__version__ = "0.1.3"

"0.1.3 : added full_matrices=False in truncated mode"


import numpy as np
import scipy.linalg as lin
import math
import sys
debug = False
import unittest
import time
if sys.version_info[0] < 3:
    pass
else:
    xrange = range

truncated = True    # When True, a flag is used which allow a slightly faster algorithm for SVD to be used.
# from NPKdata
def as_cpx(arr):
    return arr.view(dtype="complex")
def as_float(arr):
    return arr.view(dtype="float")

#-----------------------------------------------
def cadzow(data_arg, n_of_line=5, n_of_iter=5, orda=10):
    """ apply the cadzow procedure to remove noise from the current FID.
    should be followed with an FT or an LP-SVD analysis
    
    will return a dataset of the same type (float or cpx) as the input
    
    will do nothing if the data is null
    """
#    ord = max(4*n_of_line,(data.size)/3)
#    ord = min(4*n_of_line, (data.size+1)/2, 1000)
    data = data_arg
    if np.allclose(data,0.0):
        return data
    orda = min(orda, (data.size+1)/2, 4000)
    if debug:
        print("orda :", orda)
        print("processing", data.dtype)
    fid1 = data
    for i in xrange(n_of_iter):
        if debug: print("fid1", fid1.dtype, fid1.shape)
        try:
            (U, S, Vh) = dt2svd(fid1, orda=orda)
        except np.linalg.LinAlgError: # this may happend sometime ???
            print("SVD error, exiting iteration loop")
            break   # usually not at first iteration, so just break
        if debug: print("U S V", U.dtype, S.dtype, Vh.dtype)
        S1 = svdclean(S, keep=n_of_line, remove=1)
        if debug: print("S1", S1.dtype)
        fid1 = svd2dt(U, S1, Vh)
        if debug: print("fid1", fid1.dtype, fid1.shape)
    
    if data.dtype == "float":
        fid1 = np.real(fid1)
#    print "ret", fid1.dtype, fid1.shape
    return fid1
    

#-----------------------------------------------
def cadzow1d(d1D, n_of_line=5, n_of_iter=5, orda=100):
    """
    applies the cadzow procedure to a 1D dataset
    """
    if d1D.axis1.itype == 0:    # if real
        d1D.buffer = cadzow(d1D.buffer, n_of_line, n_of_iter, orda)
    else:    # if complex
        d1D.buffer = as_float(cadzow(as_cpx(d1D.buffer), n_of_line, n_of_iter, orda))
#-----------------------------------------------
def cadfun(iterelem):
    "utility for cadzow2d - has to be at top level"
    fid, n_of_line, n_of_iter, orda = iterelem
    if fid.axis1.itype == 0:    # if real
        corr = cadzow(fid.buffer, n_of_line, n_of_iter, orda)
    else:    # if complex
        corr = as_float( cadzow(as_cpx(fid.buffer), n_of_line, n_of_iter, orda) )
    return corr
#-----------------------------------------------
def cadzow2d(d2D, n_of_line=5, n_of_iter=5, orda=100, mp=True, N_proc=None, verbose=1):
    """
    applies the cadzow procedure to all columns of the 2D

    if mp,  does it in a multiprocessing fashion using multiprocessing.Pool()
    if N_proc is None, finds the optimum numer itself.
    verbose = 1 is minimum output, verbose 0 is no output
    """
    # ver = sys.version_info
    # if ver[0] <= 2 and ver[1] < 4:
    #     raise Exception("MP requires python 2.4 or higher")
    # if ver[0] == 2 and (ver[1] == 4 or ver[1] == 5):
    #     import processing as mp
    # else:
    #     import multiprocessing as mp
    if mp:
        import multiprocessing as mproc
    import itertools as it
    d2D.check2D()
    d1D = d2D.col(0)            # sets all parameters
    iterlist = it.izip(d2D.xcol(), it.repeat(n_of_line), it.repeat(n_of_iter), it.repeat(orda))    
    if mp:
        pool = mproc.Pool(processes=N_proc)  #         start N_proc independent processes  
        result = pool.imap( cadfun, iterlist )     # and apply the computation on them
    else:
        result = it.imap( cadfun, iterlist)
    for i in xrange(d2D.size2):     # do and copy the results
        if verbose>0: print("processing column %d / %d"%(i+1, d2D.size2))
        d1D.buffer = as_float(result.next())
        d1D.check()
        d2D.set_col(i, d1D)
    if mp:
        pool.close()    # closes the pool completely (useless ?)
        pool.join()

cadzow2dmp = cadzow2d   # for compatibility with previous version...

#----------------------

def dt2svd(data, orda=5):

    """computes SVD at order 'orda' from current 1D data (FID)
    data should be a numpy complex 1D array
    
    Calculates the rectangular matrix X(size-order,order) from the data and perform its singular value decomposition.
    This is the first step of the LP-SVD spectral analysis method.
    The same singular value decomposition can be used to perform forward and backward analysis.

    data is untouched

    returns (U, S, Vh)

    see also : svd2dt svdlist svdclean1 svd2ar
    """
    # check values
    (n,) =  data.shape
    n1 = n-orda+1
    if n1 < orda:
        raise Exception("orda is too large for this data, orda max : %d"%n/2)
    # build matrix
    X = np.empty((n1, orda),'complex_')   # matrix is complex
    if debug: print("Hankel matrix (%d,%d)"%X.shape)
    for i in xrange(n1):
        X[i,:] = data[i:(i+orda)].copy() # build hankel matrix
    # solve it
    if truncated :     # careful with full_matrices=False - modifies use of lin.diagsvd
        (U, S, Vh) = lin.svd(X, full_matrices = False)
    else:
        (U, S, Vh) = lin.svd(X)
    return (U, S, Vh)
#-----------------------------------------------
def svd2dt(U, S, V):
    """
    Recalculate the data which would produce the SVD decomposition U S V,
    Does this by approximating the Hankel data matrix x from the svd and u,v

    U S V are untouched,
    
    return data

    see also : dt2svd svdclean svd2ar
    """
    # determines sizes
    if debug: print("########in svd2dt")
    M = U.shape[0]  # 14
    N = V.shape[0]  # 6
    size = M+N-1
    data = np.empty((size,), dtype="complex_")
    if debug: print("M, N, size",M, N, size)
    # recompute matrix
    if truncated:       # take min, because matrices are truncated - due to the full_matrices=False in svd step
        MN = min(M,N)
        Sig = np.mat(lin.diagsvd(S, MN, MN))
    else:
        Sig = np.mat(lin.diagsvd(S, M, N))
    if debug: print("U S V : (%d x %d)  (%d x %d) (%d x %d)"%(U.shape+Sig.shape+V.shape))
    X = U*Sig*V
    Xt = X[::-1,:]
    for k in xrange(size):  # 0..19     # rebuild data from Hankel matrix,  a bit painful...
        data[k] = np.diag(Xt,k-M+1).mean()
#         temp = complex(0.0, 0.0)
#         l = max(0, k-M+1)   
#         s = min(N, k+1)
# #        print k,l,s, len(range(l,s))
#         for i in xrange(l, s):
#             temp = temp + X[k-i, i]
#         data[k] = temp/(s-l)
    if debug : print("len(data) at the end of svd2dt ",len(data))
    return data
#-----------------------------------------------
def svdclean(svd, keep=0, threshold=0, remove=1):
    """
    removes noise related features from a svd spectrum
    two methods available (uses whichever is available (eventually both)):
    keep != 0 : the number of svd to keep
    threshold != 0 : the minimum value to keep
    when remove == 1, the mean power of the removed SVD is removed from the remaining ones
    
    svd is untouched
    """
    svdn = svd.copy()
    if keep != 0:
        svdn[keep:] = 0.0
    if threshold != 0:
        svdn = np.where(svdn<threshold, 0.0, svdn)
    if remove == 1:
        power =  np.sum( (svd-svdn)**2)               # total power removed
        n = len(svd) - len(np.nonzero(svdn)[0])     # number of entry removed
        if debug: print("reduced power : %f"%math.sqrt(power/n))
        if n > 0 :
            svdn = np.sqrt( (svdn**2 - power/n).clip(0.0, svd[0]**2) ) # do not forget 
    return svdn
#-----------------------------------------------
class CadzowTest(unittest.TestCase):
    "  - Testing Cadzow mathematics - "
    def assertAlmostList(self, a, b, places=7):
        "apply asserAlmostEqual on two list of numbers"
        for ia,ib in zip(a,b):
            self.assertAlmostEqual(ia, ib, places=places)
    
    def mfft(self,v):
        "utility that returns the modulus of the fft of v"
        import scipy.fftpack as fft
        s0 = fft.fft(v)
        s0 = np.real(np.sqrt(s0*s0.conj()))   # ref spectrum
        return s0
        
    def test1D(self):
        """
        ============= a  demo / test of cadzow noise cleaning technique ==============
        applied on an additive noise
        """
        from ..Display import testplot
        plt = testplot.plot()
        print(self.test1D.__doc__)
        # simulate fid
        LB = 3      # linewidth
        N = 2000    # length of FID
        noise = 30.0 # noise level
        x = np.arange(N*1.0)/N          # 1000 Hz on N points
        fid0 = 1j*np.zeros_like(x)      # complex fid
        for i in range(1, 6):
            fid0 +=  i*20*np.exp(2*i*432.1j*x)*np.exp(-LB*x)   # that's 5 lines altogether
#        fid0 = np.real(fid0)
        s0 = self.mfft(fid0)
        plt.subplot(3,1,1)
        plt.plot( s0 , label="initial spectrum")
        plt.legend()
        fid = fid0 + noise*(np.random.randn(x.size)) # + 1j*np.random.randn(x.size))       # add some noise
        plt.subplot(3,1,2)
        plt.plot( self.mfft(np.exp(-LB*x)*fid), label="noised (filtered) spectrum")                              # adapted filter - the best one can classically do.
        plt.legend()
        # now try to clean it
        fidn = np.zeros_like(fid)
        NN = 1
#        t0 = time.time()
        for i in range(1, NN+1):  # makes NN computation with varying order
            o = N//4  # N*i/10
            fid1 = cadzow(fid, n_of_line=10, n_of_iter=2, orda=o)
            fidn = fidn+fid1
#        print "CALCUL : ",time.time()-t0
        spn = self.mfft(fidn)*(1.0/(NN))
        plt.subplot(3,1,3)
        plt.plot(spn, label="Cadzow cleaned")
        plt.legend()
        plt.show()
        diff = np.abs(s0-spn)   # diff between noise-free and cadzow-filtered
        self.assertTrue( np.sum(diff)/N < 2*noise)  # diff should be small enough
        self.assertTrue( max(diff) < 10000 )
    def _test2D(self):
        """
        ==============test for multiprocessing in cazow2d()===============
        This test might fail because svd is multithreaded on MKL, so mp version may actually be slower !!!!
        """
        from .. import NPKData
        import multiprocessing as mproc
        print(self.test2D.__doc__)
        d1 = NPKData.NPKData(buffer=np.random.rand(500,200))  # create fake data
        d1.axis1.itype=1
        d2 = d1.copy()
        print("one processor")
        t0 = time.time()
        cadzow2d(d2, n_of_line=5, n_of_iter=3, orda=20, mp=False, verbose=0)
        tmono = time.time()-t0
        print("Time : ",tmono)
        d2 = d1.copy()
        N = mproc.cpu_count()
        if N>1:
            print(N,"processors")
            if N>8:
                print("Limiting to 8 proc")
                N = 8
            t0 = time.time()
            cadzow2d(d2, n_of_line=5, n_of_iter=3, orda=20, mp=True, verbose=0, N_proc=N)
            tduo = time.time()-t0
            print("Time : ",tduo)
            self.assertTrue(tduo<tmono)
        else:
            print("test not valid as you have only one processor !")
            

   
            
    
#################################################
if __name__ == '__main__':
    unittest.main()

