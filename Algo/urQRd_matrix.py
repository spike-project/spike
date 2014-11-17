"""
rQRd.py

Created by Lionel Chiron and Marc-Andr\'e on 2012-04-04.
Copyright (c) 2012 IGBMC. All rights reserved.

version 1.0 

#########
Algorithm rQRd method for  denoising time series, named rQRd (standing for random QR denoising)

main function is 
rQRd(data, rank)
data : the series to be denoised
rank : the rank of the analysis

urQRd of the 27/09/2012

"""

import numpy as np
import numpy.linalg as linalg
import unittest


def rQRd(data, k, orda = None, iterations = None):
    
    """ 
    rQRd algorithm. Name stands for random QR denoising.
    From a data series return a denoised series denoised
    data : the series to be denoised
    k : the rank of the analysis
    orda : is the order of the analysis
        internally, a Hankel matrix (M,N) is constructed, with M = orda and N = len(data)-orda+1
        if None (default) orda = (len(data)+1)/2

    values are such that
    orda <= (len(data)+1)/2
    k < orda
    Omega is (orda x k)
    """
    if np.allclose(data,0.0):   # dont do anything if data is empty
        return data
    if not orda:
        orda = (data.size+1)/2
    if (2*orda > data.size):
        raise(Exception('order is too large'))
    if (k >= orda):
        raise(Exception('rank is too large'))
   
    Omega = np.random.normal(size=(orda,k))

    H = dt2Hankel(data, orda) # Data to Hankel
    H = rQRdCore(H, k, Omega)# H=QQ*H
#    denoised = Hankel2dtmp(H) # Hankel to dat
    denoised = Hankel2dt(H) # Hankel to data

    if data.dtype == "float":
        denoised = np.real(denoised)
    return denoised

def rQRdCore(H, k, Omega):
    '''
    Core of rQRd algorithm
    '''
    Y = np.dot(H, Omega) # prodcut to obtain random basis
    Q,r = linalg.qr(Y) # QR decomopsition of Y
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
    M = lengthdata-orda+1
    H = np.empty((M, orda),'complex_')   # matrix is complex
    for i in xrange(orda):
        H[:,i] = data[i:(i+M)].copy() # build hankel matrix
    return H

#-----------------------------------------------
def Hankel2dt(QQH):
    '''
    Goes from QQH to the datadenoised 
    '''
    M,N = QQH.shape
    size = M+N-1
    datadenoised = np.empty((size,), dtype="complex_")
    for k in xrange(size):
        tmean = complex(0.0, 0.0)
        l = max(0, k-M+1)
        s = min(N, k+1)
        if k>= M  :
            antidiag = np.diag(QQH[::-1,:],k-M+1)[l-(k-M+1):s-(k-M+1)]
        else : 
            antidiag = np.diag(QQH[::-1,:],k-M+1)[l:s]
        tmean+=antidiag.mean()
        datadenoised[k]=tmean

    return datadenoised
#-----------------------------------------------
def inject(arg):
    """
    this functions injects in the global dictionnary a new object
    arg is (name, value)
    
    used to set-up global before Pool
    generates many errors for pylint which does not understand the techniques
    """
    (nom, val) = arg
    globals()[nom]=val
    return None
def Hankel2dtmp_f(k):
    global QQH
    M,N = QQH.shape
    size = M+N-1
    tmean = complex(0.0, 0.0)
    l = max(0, k-M+1)
    s = min(N, k+1)
    if k>= M  :
        antidiag = np.diag(QQH[::-1,:],k-M+1)[l-(k-M+1):s-(k-M+1)]
    else : 
        antidiag = np.diag(QQH[::-1,:],k-M+1)[l:s]
    tmean+=antidiag.mean()
    return tmean
def Hankel2dtmp(QQH):
    '''
    Goes from QQH to the datadenoised 
    '''
    import multiprocessing
    inject(("QQH",QQH))
    M,N = QQH.shape
    size = M+N-1
    datadenoised = np.empty((size,), dtype="complex_")
    pool = multiprocessing.Pool()
    res = pool.imap(Hankel2dtmp_f, xrange(size))
    for i,v in enumerate(res):
        datadenoised[i] = v
    pool.close()
    pool.join()
    return datadenoised
    
#-----------------------------------------------

class rQRd_Tests(unittest.TestCase):
    def test_rQRd(  self,
                    lendata = 10000,
                    rank=100,
                    orda=4000,
                    noise = 200.0,
                    do_plot=True,
                    noisetype = "additive",
                    debug=0):
        """
        ============== example of use of rQRd on a synthetic data-set ===============
        """
        import time

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
        Freq = [(i+1+np.sqrt(10))*np.pi*500.0j for i in range(nbpeaks)]  # frequencies
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
        print "Initial Noisy Data SNR: %.2f dB - noise type : %s"%(iSNR,noisetype)

        ###########----
        fdata = mfft(data0) # FFT of noiseless signal
        fdatanoise = mfft(data)# FFT of noisy signal 

        ###########
        print "=== Running rQR algo ===",
        print "lendata:",lendata,
        print " orda:",orda,
        print ' rank:',rank

        t0 = time.time()
        datarqrd = rQRd(data, k=rank, orda=orda) # denoise signal with rQRd
        trQRd = time.time()-t0

        fdatarqrd  = mfft(datarqrd )# FFT of rQRd denoised signal
        # normrQR = norm(fdatarqrd -fdata)/norm(fdata)
        # print "= normratio ",normrQR
        print "=== Result ==="
        fSNR = SNR(datarqrd, data0)
        print "Denoised SNR: %.2f dB  - processing gain : %.2f dB"%( fSNR, fSNR-iSNR )
        print "processing time for rQRd : %f sec"%trQRd
        if noisetype == "additive":
            self.assertTrue( fSNR-iSNR>17.0 )   # should give 18 dB on 100000 points
        ################################################################# Plotting
        if do_plot:
            import matplotlib.pyplot as plt
            plt.subplot(321)
            plt.plot(data0,'b',label="clean signal")# rQR normal rQR
            plt.title('data series')
            plt.legend()
            plt.subplot(323)
            plt.plot(data,'k', label="noisy signal")# plot the noisy data
            plt.legend()
            plt.subplot(325)
            plt.plot(datarqrd ,'r', label='rQRd filtered signal') # plot the signal denoised with rQRd
            plt.legend()
            plt.subplot(322)
            plt.plot(fdata,'b',label="clean spectrum")# rQR normal rQR
            plt.title('FFT spectrum')
            plt.legend()
            plt.subplot(324)
            plt.plot(fdatanoise,'k', label="noisy spectrum")# plot the noisy data
            plt.legend()
            plt.subplot(326)
            plt.plot(fdatarqrd ,'r', label= 'rQRd filtered spectrum') # plot the signal denoised with rQRd
            plt.suptitle("Noise type : "+noisetype)
            plt.legend()
    #        plt.savefig("rQrd.png")
            plt.show()
    
if __name__ == '__main__':
    unittest.main()

