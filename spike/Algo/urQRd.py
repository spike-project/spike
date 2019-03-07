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
import numpy as np
import numpy.linalg as linalg
from scipy.linalg import norm
from numpy.fft import fft, ifft


def urQRd(data, k, orda=None, iterations=1):
    
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
    """
    if np.allclose(data,0.0):   # dont do anything if data is empty
        return data
    if not orda:
        orda = (data.size+1)/2
    if (2*orda > data.size):
        raise(Exception('order is too large'))
    if (k >= orda):
        raise(Exception('rank is too large'))
    N = len(data)-orda+1
    dd = data
    for i in range(iterations):
        Omega = np.random.normal(size=(N,k))
        Q, QstarH = urQRdCore(dd, orda, Omega) # H = QQ*H
        dd = Fast_Hankel2dt(Q,QstarH)
    denoised = dd
    if data.dtype == "float":           # this is a kludge, as a complex data-set is to be passed - use the analytic signal if your data are real
        denoised = np.real(denoised)
    return denoised

def urQRdCore(data, orda, Omega):
    '''
    Core of urQRd algorithm
    '''
  
    Y =  FastHankel_prod_mat_mat(data, Omega)
    Q,r = linalg.qr(Y) # QR decomopsition of Y
    del(r)          # we don't need it any more

    QstarH = FastHankel_prod_mat_mat(data.conj(), Q).conj().T# 
    return Q, QstarH # H approximation given by QQ*H    

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
    L = len(gene_vect)
    N = len(prod_vect)
    M = L-N+1
    prod_vect_zero = np.concatenate((np.zeros(M-1), prod_vect[::-1]))   # prod_vect is completed with zero to length L
    fft0, fft1 = fft(gene_vect), fft(prod_vect_zero)      # FFT transforms of generator vector and 
    prod = fft0*fft1                          # FFT product performing the convolution product. 
    c = ifft(prod)                        # IFFT for going back 
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
        vec_k = FastHankel_prod_mat_vec(gene_vect, prod_vect[::-1]) # used as fast Toeplitz
        vec_sum += vec_k 
    datadenoised = vec_sum*vec_mean(M,L) # from the sum on the antidiagonal to the mean
    return datadenoised


def test_urQRd(  
    lendata = 10000,
    rank = 100,
    orda = 4000,
    noise = 200.0,
    iterations=1,
    noisetype = "additive"):
    """
    ============== example of use of urQRd on a synthetic data-set ===============
    """
    import time
    from numpy import pi
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator

    def plot_param(fig,fignumb):
        ax = fig.add_subplot(fignumb)
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.yaxis.set_major_locator(MaxNLocator(4))
    def mfft(v):
        "utility that returns the modulus of the fft of v"
        import scipy.fftpack as fft
        s0 = fft.fft(v)    ###
        s0 = np.real(np.sqrt(s0*s0.conj()))   # ref spectrum
        return s0
    def SNR(noisy,target):
        "computes and return SNR value, in dB"
        return 10*np.log10(sum(abs(target)**2)/sum(abs(noisy - target)**2))
    ################################################
    # Data built for tests
    ################################################ Create the data
    nbpeaks = 8     # number of simulated signals
    LB = 1.11       # linewidth
    Freq = [(i+1+np.sqrt(10))*pi*500.0j for i in range(nbpeaks)]  # frequencies
    Amp = [(i+1)*20 for i in range(nbpeaks)]    # amplitudes

    data0 = np.zeros(lendata,dtype=complex)

    if noisetype == "additive":
        x = np.arange(lendata*1.0)/lendata          #  time series
        for i in range(nbpeaks):
            data0 +=  Amp[i] * np.exp(Freq[i]*x) * np.exp(-LB*x)

        dataadd = data0 + noise*(np.random.randn(x.size)+1j*np.random.randn(x.size))  # additive complex noise
        data=dataadd

    elif noisetype == "multiplicative":
        x = np.arange(lendata*1.0)/lendata          #  time series
        for i in range(nbpeaks):
            data0 +=  Amp[i] * np.exp(Freq[i]*x) * np.exp(-LB*x)

        data = np.zeros(lendata,dtype=complex)
        Anoise = noise/2
        Fnoise = noise/200
        for i in range(nbpeaks):
            nAmp = Amp[i] + Anoise*np.random.randn(x.size)
            nFreq = Freq[i] + Fnoise*np.random.randn(x.size)
            data +=  nAmp * np.exp(nFreq*x) * np.exp(-LB*x)
        
    elif noisetype == "sampling":
        x = np.arange(lendata*1.0)/lendata          #  time series
        xn = x + 0.5*np.random.randn(x.size)/lendata          #  time series with noisy jitter
        for i in range(nbpeaks):
            data0 +=  Amp[i] * np.exp(Freq[i]*x) * np.exp(-LB*x)
        data = np.zeros(lendata,dtype=complex)
        for i in range(nbpeaks):
            data +=  Amp[i] * np.exp(Freq[i]*xn) * np.exp(-LB*xn)
        
    elif noisetype == "missing points":
        x = np.arange(lendata*1.0)/lendata          #  time series
        for i in range(nbpeaks):
            data0 +=  Amp[i] * np.exp(Freq[i]*x) * np.exp(-LB*x)
        miss = np.random.randint(2, size=len(x))
        
        dataadd = data0*miss
        data=dataadd
    else:
        raise Exception("unknown noise type")

    iSNR = SNR(data,data0)
    print("Initial Noisy Data SNR: %.2f dB - noise type : %s"%(iSNR,noisetype))

    ###########----
    fdata = mfft(data0) # FFT of noiseless signal
    fdatanoise = mfft(data)# FFT of noisy signal 

    ###########
    print('''
    === Running urQR algo ===",
    lendata : {0}
    orda : {1}
    rank : {2}
    '''.format(lendata, orda, rank))
    
    t0 = time.time()
    datarqrd = urQRd(data, k=rank, orda=orda, iterations=iterations) # denoise signal with urQRd
    trQRd = time.time()-t0
    
    fdatarqrd  = mfft(datarqrd )# FFT of urQRd denoised signal
    # normrQR = norm(fdatarqrd -fdata)/norm(fdata)
    # print "= normratio ",normrQR
    print("=== Result ===")
    fSNR = SNR(datarqrd, data0)
    print("Denoised SNR: %.2f dB  - processing gain : %.2f dB"%( fSNR, fSNR-iSNR ))
    print("processing time for urQRd : %.2f sec"%trQRd)
    ################################################################# Plotting
    fig = plt.figure()
    plot_param(fig,321)
    plt.plot(data0.real,'b',label="clean signal")# plot the clean data
    plt.legend()
    plt.title('data series')
    plot_param(fig,323)
    plt.plot(data.real,'k', label="noisy signal")# plot the noisy data
    plt.legend()
    plot_param(fig,325)
    plt.plot(datarqrd.real ,'r', label='urQRd filtered signal') # plot the signal denoised with urQRd
    plt.legend()
    plot_param(fig,322)
    plt.plot(fdata,'b',label="clean spectrum")
    plt.legend()
    plt.title('FFT spectrum')
    plot_param(fig,324)
    plt.plot(fdatanoise,'k', label="noisy spectrum")# plot the noisy data
    plt.legend()
    plot_param(fig,326)
    plt.plot(fdatarqrd ,'r', label= 'urQRd filtered spectrum') # plot the signal denoised with urQRd
    plt.suptitle("Noise type : "+noisetype)
    plt.legend()
    plt.show()


if __name__ == '__main__':
    test_urQRd()

