"""
rQRd.py
#########
Algorithm for denoising time series, named rQRd (standing for random QR denoising)

main function is 
rQRd(data, rank)
data : the series to be denoised
rank : the rank of the analysis


Copyright (c) 2013 IGBMC and CNRS. All rights reserved.

Marc-Andr\'e Delsuc <madelsuc@unistra.fr>
Lionel Chiron <lionel.chiron@gmail.com>

This software is a computer program whose purpose is to compute rQRd denoising.

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


Created by Lionel Chiron and Marc-Andr\'e on 2012-04-04.
version 1.0

"""

import numpy as np
import numpy.linalg as linalg
import unittest

def rQRd(data, k, orda=None, iterations=1):
    
    """ 
    rQRd algorithm. Name stands for random QR denoising.
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
        H = dt2Hankel(dd, orda) # Data to Hankel
        H = rQRdCore(H, k, Omega)# H=QQ*H
        dd = Hankel2dt(H) # Hankel to data
    denoised = dd
    if data.dtype == "float":           # this is a kludge, as a complex data-set is to be passed - use the analytic signal if your data are real
        denoised = np.real(denoised)
    return denoised

def rQRdCore(H, k, Omega):
    '''
    Core of rQRd algorithm
    '''
    Y = np.dot(H, Omega) # prodcut to obtain random basis
    Q,r = linalg.qr(Y) # QR decomopsition of Y
    del(r)          # we don't need it any more
#    QstarH = np.dot(Q.conj().T, H) # Q*H
    QQH = np.dot(Q, np.dot(Q.conj().T, H))# QQ*H product
    return QQH# Approximation of H given by QQ*H    

#-----------------------------------------------
def dt2Hankel(data, orda):
    '''
    constructs the Hankel H matrix from the data
    Build the matrix by sticking shifted column vectors. 
    '''
    (lengthdata,) =  data.shape
    N = lengthdata-orda+1
    H = np.empty((orda,N),'complex_')   # matrix is complex
    for i in xrange(orda):
        H[i,:] = data[i:(i+N)].copy() # build hankel matrix
    return H

#-----------------------------------------------
def Hankel2dt(QQH):
    '''
    Goes from QQH to the datadenoised 
    '''
    M,N = QQH.shape
    size = M+N-1
    datadenoised = np.empty((size,), dtype=complex)
    Qt = QQH[::-1,:]
    for k in xrange(size):
        datadenoised[k] = np.diag(Qt,k-M+1).mean()
    return datadenoised
#-----------------------------------------------


def test_rQRd(  lendata = 10000,
                rank=100,
                orda=4000,
                noise = 200.0,
                iterations=1,
                noisetype = "additive"):
    """
    ============== example of use of rQRd on a synthetic data-set ===============
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
        s0 = fft.fft(v)
        s0 = np.real(np.sqrt(s0*s0.conj()))   # ref spectrum
        return s0
    def SNR(noisy,target):
        "computes and return SNR value, in dB"
        return 10*np.log10(sum(abs(target)**2)/sum(abs(noisy - target)**2))
    ########--------
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
    print ("Initial Noisy Data SNR: %.2f dB - noise type : %s"%(iSNR,noisetype))

    ###########----
    fdata = mfft(data0) # FFT of noiseless signal
    fdatanoise = mfft(data)# FFT of noisy signal 

    ###########
    print "=== Running rQR algo ===",
    print "lendata:",lendata,
    print " orda:",orda,
    print ' rank:',rank

    t0 = time.time()
    datarqrd = rQRd(data, k=rank, orda=orda, iterations=iterations) # denoise signal with rQRd
    trQRd = time.time()-t0

    fdatarqrd  = mfft(datarqrd )# FFT of rQRd denoised signal
    # normrQR = norm(fdatarqrd -fdata)/norm(fdata)
    # print "= normratio ",normrQR
    print "=== Result ==="
    fSNR = SNR(datarqrd, data0)
    print "Denoised SNR: %.2f dB  - processing gain : %.2f dB"%( fSNR, fSNR-iSNR )
    print "processing time for rQRd : %.2f sec"%trQRd
    ################################################################# Plotting
    fig = plt.figure()
    plot_param(fig,321)
    plt.plot(data0.real,'b',label="clean signal")# rQR normal rQR
    plt.legend()
    plt.title('data series')
    plot_param(fig,323)
    plt.plot(data.real,'k', label="noisy signal")# plot the noisy data
    plt.legend()
    plot_param(fig,325)
    plt.plot(datarqrd.real ,'r', label='rQRd filtered signal') # plot the signal denoised with rQRd
    plt.legend()
    plot_param(fig,322)
    plt.plot(fdata,'b',label="clean spectrum")# rQR normal rQR
    plt.legend()
    plt.title('FFT spectrum')
    plot_param(fig,324)
    plt.plot(fdatanoise,'k', label="noisy spectrum")# plot the noisy data
    plt.legend()
    plot_param(fig,326)
    plt.plot(fdatarqrd ,'r', label= 'rQRd filtered spectrum') # plot the signal denoised with rQRd
    plt.legend()
    plt.suptitle("Noise type : "+noisetype)
    plt.show()

if __name__ == '__main__':
    test_rQRd()

    
