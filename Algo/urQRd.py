"""
urQRd.py
#########
Algorithm for denoising time series, named urQRd (standing for "uncoiled random QR denoising")

main function is 
urQRd(data, rank)
data : the series to be denoised
rank : the rank of the analysis

Copyright (c) 2013 IGBMC. All rights reserved.
Marc-Andr\'e Delsuc <madelsuc@unistra.fr>
Lionel Chiron <lionel.chiron@gmail.com>

This software is a computer program whose purpose is to compute urQRd denoising.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

Created by Lionel Chiron and Marc-Andr\'e on 2013-10-13.

version 2.0 
28/oct/2013
"""
# import sys
# sys.path.append('util')
import numpy as np
import numpy.linalg as linalg
from numpy.fft import fft, ifft
import unittest
import time
from Algo.urQRd_optk import OPTK

debug = 1 # put to 1 for debuging message

def urQRd(data, k, orda = None, iterations = 1, optk = False):
    
    """ 
    urQRd algorithm. Name stands for uncoiled random QR denoising.
    From a data series return a denoised series denoised
    data : the series to be denoised - a (normally complex) numpy buffer
    k : the rank of the analysis
    orda : is the order of the analysis
        internally, a Hankel matrix (M,N) is constructed, with M = orda and N = len(data)-orda+1
        if None (default) orda = (len(data)+1)/2
    iterations : the number of time the operation should be repeated

    values are such that
    orda <= (len(data)+1)/2
    k < orda
    N = len(data)-orda+1
    Omega is (N x k)
    ##########
    BECAREFUL datasize must be different from a product of primes !!!!!!..
    a processing with a datasize of 120022 for example will be 50 times longer than
    a procesing of a datasize of 120000.
    ##########
    """
    if optk:                                                            #  optimal rank, ensures that all the signal part is retrieved.
        optrank = OPTK(data, orda)
        k = optrank.find_best_rank()  
    if np.allclose(data,0.0):                                           # dont do anything if data is empty
        return data
    if not orda:
        orda = (data.size)/2                                            # defining orda if not given.
    if (2*orda > data.size):                                            # checks if orda not too large.
        raise(Exception('order is too large'))
    if (k >= orda):                                                     # checks if rank not too large
        raise(Exception('rank is too large'))
    N = len(data)-orda+1
    dd = data
    for i in range(iterations):
        Omega = np.random.normal(size=(N,k))                            # Omega random real gaussian matrix Nxk
        Q, QstarH = urQRdCore(dd, Omega)                                # H = QQ*H
    dd = Fast_Hankel2dt(Q,QstarH)
    denoised = dd
    if data.dtype == "float":                                           # this is a kludge, as a complex data-set is to be passed - use the analytic signal if your data are real
        denoised = np.real(denoised)
    return denoised

def urQRdCore(data, Omega):
    '''
    Core of urQRd algorithm
    '''
    Y =  FastHankel_prod_mat_mat(data, Omega)
    Q,r = linalg.qr(Y)                                                  # QR decompsition of Y
    del(r)                                                              # we don't need it any more
    QstarH = FastHankel_prod_mat_mat(data.conj(), Q).conj().T# 
    return Q, QstarH                                                    # H approximation given by QQ*H    

def vec_mean(M,L):
    '''
    Vector for calculating the mean from the sum on the antidiagonal.
    data = vec_sum*vec_mean
    '''
    vec_prod_diag = [1/float((i+1)) for i in range(M)]
    vec_prod_middle = [1/float(M) for i in range(L-2*M)]
    vec_mean_prod_tot = vec_prod_diag + vec_prod_middle + vec_prod_diag[::-1]
    return np.array(vec_mean_prod_tot)

def FastHankel_prod_mat_mat(gene_vect, matrix):
    '''
    Fast Hankel structured matrix matrix product based on FastHankel_prod_mat_vec
    '''
    N,K = matrix.shape 
    L = len(gene_vect)
    M = L-N+1
    data = np.zeros(shape = (M, K), dtype = complex)
    for k in range(K):
        prod_vect = matrix[:,k] 
        data[:,k] = FastHankel_prod_mat_vec(gene_vect, prod_vect) 
    return data

def FastHankel_prod_mat_vec(gene_vect, prod_vect):
    """
    Compute product of Hankel matrix (gene_vect)  by vector prod_vect.
    H is not computed
    M is the length of the result
    """
    L = len(gene_vect)                                                  # length of generator vector
    N = len(prod_vect)                                                  # length of the vector that is multiplied by the matrix.
    M = L-N+1
    prod_vect_zero = np.concatenate((np.zeros(M-1), prod_vect[::-1]))   # prod_vect is completed with zero to length L
    fft0, fft1 = fft(gene_vect), fft(prod_vect_zero)                    # FFT transforms of generator vector and 
    prod = fft0*fft1                                                    # FFT product performing the convolution product. 
    c = ifft(prod)                                                      # IFFT for going back 
    return np.roll(c,+1)[:M]

def Fast_Hankel2dt(Q,QH):
    '''
    returning to data from Q and QstarH
    Based on FastHankel_prod_mat_vec.
    '''
    M,K = Q.shape 
    K,N = QH.shape 
    L = M+N-1
    vec_sum = np.zeros((L,), dtype = complex)
    for k in range(K):
        prod_vect = QH[k,:]
        gene_vect = np.concatenate((np.zeros(N-1), Q[:, k], np.zeros(N-1))) # generator vector for Toeplitz matrix
        vec_k = FastHankel_prod_mat_vec(gene_vect, prod_vect[::-1])         # used as fast Toeplitz
        vec_sum += vec_k 
    datadenoised = vec_sum*vec_mean(M,L)                                    # from the sum on the antidiagonal to the mean
    return datadenoised


class urQRd_Tests(unittest.TestCase):
    def test_urQRd(  self,
                    lendata = 10000,
                    rank = 100,
                    orda = 4000,
                    nbpeaks = 20,
                    noise = 200.0,
                    noisetype = "additive"):
        """
        ============== example of use of urQRd on a synthetic data-set ===============
        """
        import Display.testplot as testplot
        plt = testplot.plot()
        from util.dynsubplot import subpl
        from util.signal_tools import fid_signoise_type, SNR_dB, mfft
        ###########
        print "=== Running rQR algo ===",
        print "lendata:", lendata,
        print " orda:", orda,
        print ' rank:', rank
        data = fid_signoise_type(nbpeaks,lendata, noise, noisetype)                 # noisy signal
        fdatanoise = mfft(data)
        noise = 0
        data0 = fid_signoise_type(nbpeaks,lendata, noise, noisetype)                # clean signal
        fdata = mfft(data0)
        
        iSNR = SNR_dB(data,data0)
        print "Initial Noisy Data SNR: %.2f dB - noise type : %s"%(iSNR, noisetype)
        t0 = time.time()
        dataurqrd = urQRd(data, k = rank, orda = orda, optk =True)                  # denoise signal with urQRd
        turQRd = time.time()-t0
        fdataurqrd  = mfft(dataurqrd )# FFT of urQRd denoised signal
        print "=== Result ==="
        fSNR = SNR_dB(dataurqrd, data0)
        print "Denoised SNR: %.2f dB  - processing gain : %.2f dB"%( fSNR, fSNR-iSNR )
        print "processing time for urQRd : %f sec"%turQRd
        print fSNR-iSNR
        #####
        sub = subpl(2, plt)
        sub.next()
        #######################
        sub.plot(data0,'b',label = "clean signal")                                  # original signal
        sub.title('data series')
        sub.next()
        sub.plot(fdata,'b',label = "clean spectrum")                                # original spectrum
        sub.title('FFT spectrum')
        ######################
        sub.next()
        sub.plot(data,'k', label = "noisy signal")                                  # plot the noisy fid
        sub.next()
        sub.plot(fdatanoise,'k', label = "noisy spectrum")                          # plot the noisy spectrum
        #######################
        sub.next()
        sub.plot(dataurqrd ,'r', label = 'urQRd filtered signal')                   # plot the fid denoised with urQRd
        sub.next()
        sub.plot(fdataurqrd ,'r', label = 'urQRd filtered spectrum')                # plot the spectrum denoised with urQRd
        sub.title("Noise type : " + noisetype)
        ############################
        sub.show()
    
if __name__ == '__main__':
    unittest.main()
