#!/usr/bin/env python 
# encoding: utf-8

"""

SL0 from http://ee.sharif.ir/~SLzero/

Authors: Massoud Babaie-Zadeh and Hossein Mohimani
Version: 1.3
Last modified: 4 August 2008.

adapted to numpy by Marc-André on 2011-11-18.
Copyright (c) 2011 IGBMC. All rights reserved.
"""

from __future__ import print_function
import numpy as np
import scipy.fftpack as fft
import unittest
#    import numpy.fft as fft

def SL0(A, x, sigma_min, sigma_decrease_factor=0.5, mu_0=2, L=3, A_pinv=None, true_s=None):
    """
    SL0(A, x, sigma_min, sigma_decrease_factor, mu_0, L, A_pinv, true_s)
   
      Returns the sparsest vector s which satisfies underdetermined system of
      linear equations  A*s=x, using  Smoothed L0  (SL0) algorithm. Note that 
      the matrix  A should  be a 'wide' matrix  (more columns than rows). The 
      number of the rows of  matrix A  should  be  equal to the length of the 
      column vector x.
   
        The first 3 arguments should necessarily be provided by the user. The 
      other parameters have defult values  calculated within the function, or
      may be provided by the user.
   
      Sequence of Sigma (sigma_min and sigma_decrease_factor):
        This is a decreasing geometric sequence of positive numbers:
          - The  first  element   of  the  sequence  of  sigma is  calculated 
        automatically. The last  element  is  given  by  'sigma_min', and the 
        change factor for decreasing sigma is given by 'sigma_decrease_factor'. 
          - The default value of 'sigma_decrease_factor' is 0.5. Larger value 
        gives better results  for less sparse sources, but it uses more steps 
        on   sigma   to  reach  sigma_min,  and  hence  it  requires   higher 
        computational cost.
          - There is no default  value for  'sigma_min',  and  it  should  be 
        provided  by  the  user (depending  on his/her estimated source noise 
        level,  or  his/her  desired  accuracy).  By `noise' we mean here the
        noise in the sources, that is, the energy of the inactive elements of
        's'.   For example,  by  the  noiseless  case,  we  mean the inactive
        elements of 's' are exactly equal to zero. As a rule of tumb, for the
        noisy case,  sigma_min should be about 2 to 4  times  of the standard
        deviation of this noise.  For the noiseless case, smaller 'sigma_min'
        results in  better estimation of the sparsest solution, and hence its
        value is determined by the desired accuracy.
    
      mu_0: 
           The  value  of  mu_0  scales  the sequence of mu. For each vlue of 
        sigma, the value of  mu is chosen via mu=mu_0*sigma^2. Note that this 
        value effects Convergence.
           The default value is mu_0=2 (see the paper).
   
      L: 
           number  of  iterations of the internal (steepest ascent) loop. The
        default value is L=3.
   
      A_pinv: 
           is the  pseudo-inverse of matrix A defined by A_pinv=A'*inv(A*A'). 
        If it is not provided, it will be calculated within the function.  If
        you use this function for solving x(t)=A s(t) for different values of
        't', it would be a good idea to calculate A_pinv outside the function
        to prevent its re-calculation for each 't'.
   
      true_s: 
           is the  true value of the  sparse  solution.  This argument is for
        simulation purposes. If it is provided by the user, then the function
        will  calculate the SNR of the estimation for each value of sigma and
        it provides a progress report.
   
    Authors: Massoud Babaie-Zadeh and Hossein Mohimani
    Version: 1.4
    Last modified: 4 April 2010.
    Web-page:
    ------------------
       http://ee.sharif.ir/~SLzero

    Code History:
   --------------
   MODIF MAD : FT if true, spacity is evaluated in Fourier space
   Version 2.0: 4 April 2010
         Doing a few small modifications that enable the code to work also
         for complex numbers (not only for real numbers).
  
   Version 1.3: 13 Sep 2008
         Just a few modifications in the comments
  
   Version 1.2: Adding some more comments in the help section
   
    Version 1.1: 4 August 2008
       - Using MATLAB's pseudo inverse function to generalize for the case
         the matrix A is not full-rank.
   
    Version 1.0 (first official version): 4 July 2008.
   
    First non-official version and algorithm development: Summer 2006
    """
    if A_pinv is None:
        A_pinv = np.linalg.pinv(A)

    s = np.dot(A_pinv, x)
    sigma = 2*max(np.abs(s))

    # Main Loop
    while sigma>sigma_min:
        for i in range(L):
            delta = OurDelta(s, sigma)
            s -=  mu_0*delta
            s -=  np.dot(A_pinv, np.dot(A,s)-x)   # Projection
    
        if true_s is not None:
            print('     sigma=%f, SNR=%f' % (sigma,estimate_SNR(s,true_s)))
    
        sigma = sigma * sigma_decrease_factor
    return s

def SL0FT(A, x, sigma_min, sigma_decrease_factor=0.5, mu_0=2, L=3, A_pinv=None, true_s=None):
    """
    Simplistic implementation of FT in SL0
    same argument as SL0,

    A is a sampling/scrambling matrix
    if A is identity, trunctation or random sampling, then  A_pinv should be provided as A.T
    """
    def trans(s):
        "transform s to x"
        return np.dot(A, fft.ifft(s) )
    def ttrans(x):
        "ttransform x to s"
        return fft.fft(np.dot(A_pinv, x))
    if A_pinv is None:
        A_pinv = np.linalg.pinv(A)

    s = ttrans(x)
    sigma = 2*max(np.abs(s))

    # Main Loop
    while sigma>sigma_min:
        for i in range(L):
            delta = OurDelta(s, sigma)
            s -=  mu_0*delta
            s -=  ttrans(trans(s)-x)   # Projection
    
        if true_s is not None:
            print('     sigma=%f, SNR=%f' % (sigma,estimate_SNR(s,true_s)))
    
        sigma = sigma * sigma_decrease_factor
    return s
class transformations(object):
    """
    this class contains methods which are tools to generate transform and ttransform functions needed by
    convergence algorithms
    """
    def __init__(self, size_image, size_mesure, debug=0):
        """
        size_image and size_mesure are sizes of x and s space
        all other fields are meant to be overloaded after creation
        
        direct transform refers to S => X  // image => Data transform
        """
        self.size_image = size_image    
        self.size_mesure = size_mesure
        self.pre_ft = self.Id       # applied before FT (in the direct transform)
        self.tpre_ft = self.Id     # its transpose (so applied after FT in ttransform)
        self.post_ft = self.Id      # applied after FT (in the direct transform)
        self.tpost_ft = self.Id     # its transpose (so applied before FT in ttransform)
        self.ft = fft.ifft          # FT in the direct transform
        self.tft = fft.fft          # FT in the inverse transform
        self.sampling = None        # sampling vector if not none
        self.debug = debug
    def Id(self, x):
        return x
    def check(self):
        if self.debug:
            print("""
size_image: %d - size_mesure: %d
sampling
%s"""%(self.size_image, self.size_mesure, str(self.sampling)))
        if self.sampling is not None:
            assert(len(self.sampling) == self.size_mesure)
            assert(max(self.sampling) <= self.size_image)
        # assert( (self.pre_ft == self.Id and self.tpre_ft == self.Id) or \
        #         (self.pre_ft != self.Id and self.tpre_ft != self.Id) )                  # if one, then both !
        # assert( (self.post_ft == self.Id and self.tpost_ft == self.Id) or \
        #         (self.post_ft != self.Id and self.tpost_ft != self.Id) )                # if one, then both !
            
    def sample(self, x):
        """
        apply a sampling function - using self.sampling
        """
        return x[self.sampling]
    def tsample(self, x):
        """
        transpose of the sampling function
        """
        xx = np.zeros(self.size_image,'complex128')
        #print xx.dtype, x.dtype, self.sampling.dtype
        xx[self.sampling] = x
        return xx
    def transform(self, s):
        """
        transform s (image) to x (data)
        pre_ft() : s->s
        ft()     : s->x - fft.ifft by default - should not change size
        post_ft() : x->x - typically : broadening, sampling, truncating, etc...
        """
        if self.debug: print('entering trans', s.shape)
        if self.pre_ft != self.Id:
            if self.debug: print('trans pre_ft')
            s = self.pre_ft(s)
        x = self.ft(s)
        if self.post_ft != self.Id:
            if self.debug: print('trans post_ft')
            x = self.post_ft(x)
        if self.sampling is not None:
            if self.debug: print('trans sample')
            x = self.sample(x)
        if self.size_mesure != len(x):      # eventually truncate
            if self.debug: print('trans trunc')
            x = x[0:self.size_mesure]
        if self.debug: print('exiting trans', x.shape)
        return x
    def ttransform(self, x):
        """
        the transpose of transform
        transform x to s
        """
        if self.debug: print('entering ttrans', x.shape, x.dtype)
        if self.sampling is not None:
            if self.debug: print('ttrans sample')
            x = self.tsample(x)
        elif self.size_image != len(x):      # eventually zerofill
            if self.debug: print('ttrans zerofill',len(x),self.size_image)
            xx = np.zeros(self.size_image, dtype='complex')
            xx[:len(x)] = x[:]
            x = xx
        if self.tpost_ft != self.Id:
            if self.debug: print('ttrans tpost_ft')
            x = self.tpost_ft(x)
        s = self.tft(x)
        if self.tpre_ft != self.Id:
            if self.debug: print('ttrans tpre_ft')
            s = self.tpre_ft(s)
        if self.debug: print('exiting ttrans', s.shape)
        return s
def SL0gene(transf, ttransf, x, sigma_min, sigma_decrease_factor=0.5, mu_0=2, L=3, true_s=None):
    """
    general implementation of SL0
    most arguments as SL0
    
    transform and ttransform are generalized linear operation
    simplistic implementation could be :
    def trans(s):
        "transform s to x"
        return np.dot(A, fft.ifft(s) )
    def ttrans(x):
        "ttransform x to s"
        return fft.fft(np.dot(A_pinv, x))
    
    just write your own (eventually algorithmic) or used provided ones
    """
    s = ttransf(x)
    sigma = 2*max(np.abs(s))
    print("sigma initial:", sigma)
    # Main Loop
    while sigma>sigma_min:
        for i in range(L):
            delta = OurDelta(s, sigma)
            s -=  np.dot(mu_0, delta)
            s -=  ttransf(transf(s)-x)   # Projection
    
        if true_s is not None:
            print('     sigma=%f, SNR=%f' % (sigma,estimate_SNR(s,true_s)))
        else:
            print('     sigma=%f' % (sigma))
        sigma = sigma * sigma_decrease_factor
    return s
def OurDelta(s,sigma):
#    return s*np.exp(-(s*np.conj(s))/sigma**2)
    return s*np.exp(-abs(s)**2/sigma**2)

def estimate_SNR(estim_s, true_s):
    err = true_s - estim_s
#    return 10*np.log10(np.sum(true_s*np.conj(true_s))/sum(err*np.conj(err)))
    return 10*np.log10(sum(abs(true_s)**2)/sum(abs(err)**2))

def norm(v):
    return np.real(np.sqrt(np.sum(v * v.conj())))


class SL0_Tests(unittest.TestCase):
    def test_igor1(self):
        """
        from http://nuit-blanche.blogspot.com/2011/11/how-to-wow-your-friends-in-high-places.html
        code from Igor Carron
        """
        from ..Display import testplot
        plt = testplot.plot()
        M = 20
        N = 500
        noise = 0.01
        A = np.random.randn(M,N)
        x = np.zeros(N)
        x[10] = 1
        x[140] = 1
        y1 = np.dot(A, x) + noise*np.random.randn(M)
        plt.subplot(2,2,1)
        plt.plot(x,'o')
        plt.title(' original')

        x2 = np.linalg.lstsq(A, y1)[0]
        plt.subplot(2,2,2)
        plt.plot(x2,'*')
        plt.title(' using lstsq')

        A3 = np.dot(A.T.conj(), np.linalg.inv(np.dot(A, A.T.conj())))
    #    A3 = linalg.pinv(A)
        x3 = np.dot(A3, y1)
        plt.subplot(2,2,3)
        plt.plot(x3,'*')
        plt.title(' using pinv')

        plt.subplot(2,2,4)
        # following parameters for all my computations: sigma_min = 0.00004, sdf= 0.95 )
        sigma_min = 0.0001
        sdf= 0.95
        mu = 2
        L = 3
        x4 = SL0(A, y1, sigma_min, sdf, mu, L, true_s=x)
        plt.plot(x4,'*')
        plt.title(' using SL0 solver')
        plt.show()

if __name__ == '__main__':
    unittest.main()