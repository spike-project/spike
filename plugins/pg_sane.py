#!/usr/bin/env python 
# encoding: utf-8

"""plugin for PG_Sane algorithm,  used for NUS processing

It takes a NUS acquired transient, and fill it by estimating the missing values.

associated publications

- Lionel Chiron, Afef Cherni, Christian Rolando, Emilie Chouzenoux, Marc-André Delsuc
  Fast Analysis of Non Uniform Sampled DataSets in 2D-FT-ICR-MS. - in progress

- Bray, F., Bouclon, J., Chiron, L., Witt, M., Delsuc, M.-A., & Rolando, C. (2017).
  Nonuniform Sampling Acquisition of Two-Dimensional Fourier Transform Ion Cyclotron Resonance Mass Spectrometry for Increased Mass Resolution of Tandem Mass Spectrometry Precursor Ions.
  Anal. Chem., acs.analchem.7b01850. http://doi.org/10.1021/acs.analchem.7b01850

- Chiron, L., van Agthoven, M. A., Kieffer, B., Rolando, C., & Delsuc, M.-A. (2014).
  Efficient denoising algorithms for large experimental datasets and their applications in Fourier transform ion cyclotron resonance mass spectrometry.
  PNAS , 111(4), 1385–1390. http://doi.org/10.1073/pnas.1306700111

"""

from __future__ import print_function

import unittest
import numpy as np
from numpy.fft import fft, ifft, rfft, irfft

from spike.NPKData import NPKData_plugin,  as_cpx, as_float, _base_fft,\
            _base_ifft, _base_rfft, _base_irfft
from spike.Algo.sane import sane
from spike.util.signal_tools import filtering

import sys #
if sys.version_info[0] < 3:
    pass
else:
    xrange = range


def noise(data, iterations=10, tresh=3.0):
    "Simple noise evaluation "
    b = data.copy()
    for i in range(iterations):
        b = b[ (b-b.mean()) < tresh*b.std() ]
    return b.std()

def HT(x, thresh):
    """
    returns the Hard Thresholding of x,
    i.e. all points such as x_i <= thresh are set to 0.0
    """
    ax = abs(x)
    tx = np.where(ax>thresh, x, 0.0)
    return tx

def HTproj(x, k):
    """
    returns the Hard Thresholding of x, on the ball of radius \ell_o = k
    i.e. the k largest values are kept, all other are set to 0
    """
    ax = abs(x)
    N = len(ax)-k
    tx = np.argpartition(ax, N)
    hpx = np.zeros_like(x)
    hpx[tx[N:]] = x[tx[N:]]
    return hpx

def pg_sane(npkd, HTmode='projection', axis=1, iterations=10, HTratio=None, rank=20, Lthresh=2.0, sampling=None, size=None, final='SANE'):
    """
    Papoulis Gershberg algorithm - stabilized with SANE
    
    This function takes FID with a partial sampling, and fills the holes with estimated values.
    The FID can then be processed normally, as it was completely acquired.
    
    mode (threshold or projection) determine PG algorithm
    iterations : the number of iterations used for the program
        pg_sane converges quite rapidely, and usually 5 to 20 iterations are suffisant.
        when the SNR is high, or for large data-sets (>8k) more iterations might be needed
    HTratio:    number of lines left by HT/projection - usually a few % of the whole spectrum
        default is 0.01 ( 1% )
    rank : a rough estimate of the number of lines present in the FT spectrum of the FID - not stringent -
            if in doubt, use a slightly larger value.
    Lthresh : the multiplier for HT/threshold, usually 2.0
        lower recovers more weak signals but requires more iterations
        higher requires less iteration but look only to larger peaks.
    sampling : if provided, the sampling index list and npkd is supposed to be zerofilled (missing values set to zero)
               if None, npkd.axis1.sampled should be True, and sampling it will be fetched via npkd.axis1.get_sampling()
    final: the final step after iterations,
           default is 'SANE': better (noise-less) data-sets, the best reconstruction quality in simulations
           'PG' to reapply PG step - produces the cleanest/more compressible data-sets
           'Reinject': the closest to acquired data - to use if there is very little noise.
    size: if different from None, pg_sane will also extrapolate to size
    """
    def PG(fid):
        "local version of the PG algorithm"
        b = directFT(fid)            # directFT() defined locally below
        if HTmode == 'threshold':
            bf = HT(b, Lthresh*noise(b))
        else:
            bf = HTproj(b, Ndata)
        return idirectFT(bf)

    if npkd.dim == 1:
        if sampling is None:
            if not npkd.axis1.sampled :
                raise Exception("this function works only on NUS datasets")
            fiditer = npkd.copy().zf().get_buffer()         # initialize
            lsampling = npkd.axis1.get_sampling()
        else:
            fiditer = npkd.get_buffer()                     # initialize
            lsampling = sampling

        GoodOnes = fiditer[lsampling]
        RATIO = float(len(GoodOnes))/len(fiditer)       # evaluate mean power
        power = np.linalg.norm(fiditer)/np.sqrt(RATIO)
        if size is not None and size>len(fiditer):       # extend is size is required
            fidzf = np.zeros(size, dtype=fiditer.dtype)
            fidzf[:len(fiditer)] = fiditer[:]
            fiditer = fidzf
        if HTratio is None:
            Ndata = len(fiditer)//100                    # 1% for HT projection 
        else:
            Ndata = int(HTratio * len(fiditer))

        cpx = (fiditer.dtype == 'complex128')           # assign FT
        if cpx:
            directFT = fft
            idirectFT = ifft
        else:
            directFT = rfft
            idirectFT = lambda sp: irfft(sp, len(fiditer))  # check irfft() doc for the len() trick
        
        
        for i in range(iterations): # then loop
            fiditer[lsampling] = GoodOnes                            #Reinject
            fiditer = sane(fiditer, rank)                           # Apply Sane
            if i == 0:
                fiditer *= power/(np.linalg.norm(fiditer))   # normalize first SANE step
            fiditer[lsampling] = GoodOnes                            #Reinject
            fiditer = PG(fiditer)                                   # Apply PG
            if i == 0:
                fidnorm = np.linalg.norm(fiditer)
                if fidnorm>0.0:                  # happends in some pathological cases
                    fiditer *= power/(fidnorm)   # normalize first PG step

        if final == 'SANE':
            fiditer = sane(fiditer, rank)
        elif final == 'PG':
            pass  # nothing special to do
        elif final == 'Reinject':
            fiditer[lsampling] = GoodOnes
        else:
            raise Exception('wrong mode for "final", choose "SANE", "PG", or "Reinject"')
        # we're done
        npkd.set_buffer(fiditer)
    elif npkd.dim == 2:
        todo = npkd.test_axis(axis)
        if todo == 2:
             for i in xrange(npkd.size1):
                 r = npkd.row(i).pg_sane(rank=rank, iterations=iterations, Lthresh=Lthresh, sampling=sampling, final=final, size=size)
                 npkd.set_row(i,r)
        elif todo == 1:
             for i in xrange(npkd.size2):
                 r = npkd.col(i).pg_sane(rank=rank, iterations=iterations, Lthresh=Lthresh, sampling=sampling, final=final, size=size)
                 npkd.set_col(i,r)
    elif npkd.dim == 3:
        raise Exception("not implemented yet")
    npkd.axes(axis).sampling = None
    return npkd
    

class sane_pgTests(unittest.TestCase):
    def test_NUS_sampling(self):
        '''
        NUS example
        removing the sampling noise 
        '''
        from ..Tests import filename, directory
        import spike.Display.testplot as testplot
        plt = testplot.plot()
        import spike.util.signal_tools as u
        from numpy.fft import fft
        from spike.Tests import filename
        from spike.NPKData import NPKData
        samplingfile = filename("Sampling_file.list")
        e = NPKData(dim = 1)
        e.axis1.load_sampling(samplingfile)
        size = 20000
        signal = u.SIGNAL_NOISE(lenfid=size, nbpeaks=10, amplitude=100, noise=50, shift = 7000)
        signal.fid
        echant = signal.fid[e.axis1.sampling]
#        print "echant.size ",echant.size
        f = NPKData(buffer = echant)
        f.axis1.load_sampling(samplingfile)
        h = f.copy()
        h.pg_sane()
        pgdb = u.SNR_dB(signal.fid0,h.get_buffer())
        print("PG_SANE reconstruction is %.1fdB should be greater than 19dB"%(pgdb))
        f.zf().fft().modulus().display(label = "NUS : FFT with sampling noise")
        h.fft().modulus().display(label = "NUS : processed with PG_SANE", new_fig = False, show = True)
        self.assertTrue(pgdb>19.0)
    def test_NUS_sampling2(self):
        '''
        NUS larger example
        removing the sampling noise 
        '''
        from spike.NPKData import NPKData
        import spike.Display.testplot as testplot
        plt = testplot.plot()
        import time
        # First tools
        def build(tp):
            """
            build a synthetic transient signal
            containing 10 frequencies with intensities from 1 to 10
            return the time axis and the signal
            """
            # set random signal lines
            np.random.seed(123)
            freq = 3E5*np.random.rand(10)   # 10 random lines
            amp = range(1,11)
            # build transient
            fid = np.zeros(N, dtype=complex)
            for a,nu in zip(amp, freq):
                fid += a*np.exp(-tp/tau)*np.exp(2j*np.pi*nu*tp)
            return fid
        # compute spectrum
        def FT(v):
            " Fourier transform, with a simple cosine-bell apodisation"
            vapod = v*np.cos(np.linspace(0,np.pi/2, len(v)))
            return  np.fft.fftshift(np.fft.fft(vapod))
        # generate sampling list
        def gene_sampling(ratio):
            np.random.seed(1234)
            perm = np.random.permutation(N-1)  # generate a permutation 
            sampling = sorted(perm[:int(round(N*ratio))-2])
            sampling.insert(0,0)  # insure the first point is set
            sampling.append(N-1)  # insure the last point is set
            return sampling
        def noise(data):
            " estimate noise in the spectrum by iteratively removing all signals above 3 sigma "
            b = data.copy()
            for i in range(10):
                b = b[ b-b.mean()<3*b.std() ]
            return b.std()

        # Then generate data
        N = 64000           # Size of the complete sampling
        SR = 1000000.0      # Spectral Width
        tau = .1           # ion life-time - make it short enough to reduce FT artifacts
        dt = 1/SR           # Nyquist sampling, the example are presented in complex
        tp = np.arange(0,N)*dt     # time axis
        fq = np.linspace(0,SR,N)   # frequency axis
        fid0 = build(tp)  # noise-free data
        NOISE = 5   # noise level, same unit as the amp in the first cells.
                    # here, noise level is half of the largest amplitude.
        np.random.seed(12345)
        # noisy data
        nfid = fid0 + NOISE*np.random.randn(len(fid0)) + 1j*NOISE*np.random.randn(len(fid0))
        # generate sampling
        RATIO = 1./8
        sampling = gene_sampling(RATIO)
        # prepare
        f = NPKData(buffer = nfid[sampling])
        f.axis1.sampling = sampling
        # do it
        t0 = time.time()
        g = f.copy().pg_sane( iterations=20, rank=15)
        elaps = time.time() - t0
        SNR = -20*np.log10(np.linalg.norm(g.get_buffer()-fid0) / np.linalg.norm(fid0) )
        print ("test_NUS_sampling2: elaps %.2f sec  SNR: %.1f dB should be larger than 30dB"%(elaps, SNR))
        self.assertTrue(SNR>30.0)
        ax1 = plt.subplot(211)
        f.copy().apod_sin().zf(2).fft().display(title='spectrum original data with sampling noise', figure=ax1)
        ax2 = plt.subplot(212)
        g.copy().apod_sin().zf(2).fft().display(title='spectrum after pg_sane cleaning', figure=ax2)



NPKData_plugin("pg_sane", pg_sane)
