#!/usr/bin/env python 
# encoding: utf-8

'''
Created by Lionel Chiron  18/10/2013 
Copyright (c) 2013 __NMRTEC__. All rights reserved.
Various tools for performing signal processing (SNR calculation,
    synthetic Fids with noise, noise level etc... )
'''

from __future__ import print_function, division
import numpy as np
from numpy.fft import rfft as nprfft, irfft as npirfft
from numpy import pi
import math
import scipy.fftpack as fft
from scipy.signal import firwin2, firwin, convolve, lfilter
from scipy.optimize import minimize
#from matplotlib import pyplot as plt
from ..Display import testplot
plt = testplot.plot()

import unittest
import time
#from .. import NPKData

def filtering(sig, window = None):
    '''
    Filter for real signal
    window : list [f1, f2], eg: [0.1, 0.3], 0 < f < 1
    '''
    def FIR_filter(sig, nbpoles = 80, window = None ):
        filt = firwin(nbpoles, window, pass_zero = False) 
        filtered_sig = lfilter(filt, 1.0, sig)
        return filtered_sig
    filtered_sig = FIR_filter(sig, nbpoles = 100, window = window)
    return filtered_sig

def SNR_dB(noisy,target):
    "computes and return SNR value, in dB"
    return 10*np.log10(sum(abs(target)**2)/sum(abs(noisy - target)**2))

def mfft(v):
     "utility that returns the modulus of the fft of v"
     
     s0 = fft.fft(v)
     s0 = np.sqrt(s0*s0.conj())   # ref spectrum
     return s0

def mrfft(v):
      "utility that returns the modulus of the rfft of v"
      s0 = nprfft(v)
      s0 = np.sqrt(s0*s0.conj())   # ref spectrum
      return s0

class SIGNAL_NOISE(object):
    '''
    Tool for generating SIGNAL with/without NOISE. 
    
    input:
        lenfid : length of the Fid
        nbpeaks : nb of frequencies
        amplitude : reference amplitude for signal
        LB : exponential apodisation factor
        noise : noise amplitude
        shape : shape of the spectrum. Available : triangular, list. 
        noisetype : kind of noise, ie : additive, multiplicative etc...
        shift : shift of the frequencies, the shift is expressed in points.
        trunc : truncaton ratio for the Fid. 
        seed : for using a precise seed for random processes.
        baseline : boolean for specifying if a baseline is added to the spectrum
    output:
        Noisy signal : self.fid
        Signal without noise :  self.fid0
        Noisy spectrum : self.spec
        Spectrum without noise :  self.spec0
    
    eg1: # simplest usage
        sig = SIGNAL_NOISE( lenfid= 1000, nbpeaks = 20, amplitude = 27, noise = 5)
    eg2:  # with a list of amplitudes
        lamplit = list(np.random.randn(20)*10)
        sig = st.SIGNAL_NOISE( lenfid= 10000, nbpeaks = 20,  amplitude = lamplit, noise = 4, shape = "list")
    
    '''
    
    def __init__(self, lenfid, nbpeaks, amplitude, noise, LB = 1.11,\
        shape = "triangular", noisetype = 'additive', shift = 0,\
                realfid = False, trunc = None, seed = None, baseline=None):
        self.LB = LB # linewidth
        self.nbpeaks = nbpeaks   # number of lines in the signal
        self.lenfid = lenfid     # length of the fid
        self.noise = noise            # noise amplitude
        self.noisetype = noisetype      # kind of noise
        self.shift = shift  # shift the spectrum in frequency domain.
        self.realfid = realfid
        self.trunc = trunc # factor of truncation of the Fid. 
        if not seed:
            np.random.seed(seed)                    # frequency shift
        self.baseline = baseline
        self.x = np.arange(self.lenfid*1.0)/self.lenfid          #  time series
        self.fid = np.zeros(self.lenfid, dtype = complex) # fid without noise initialization
        self.fid0 = np.zeros(self.lenfid, dtype = complex) # noisy fid initialization
        self.shape = shape # shape of spectrum
        self.slope = 20 # slope for the triangular shape
        self.mult_Anoise = noise/2      # amplitude for multiplicative noise
        self.mult_Fnoise = noise/200  # drift in frequency for multiplicative noise
        self.amplitude = amplitude # amplitude of the signal
        #######
        self.prepare_amplitudes()          # make the list of amplitudes
        self.prepare_frequencies()         # make the list of frequencies
        self.make_pure_signal()         # make the signal without noise
        self.make_noisy_signal()        # make the noisy signal
        if self.realfid:
            self.mfft = mrfft
            self.fid = np.real(self.fid)
            self.fid0 = np.real(self.fid0)
        else: 
            self.mfft = mfft
        self.fid_orig = self.fid.copy() # keep initial Fid when performing truncation etc.. 
        self.fid0_orig = self.fid0.copy() # keep initial Fid when performing truncation etc.. 
        
    def cut(self):
        '''
        Uses Fid truncated
        '''
        self.fid[self.lenfid//self.trunc:] = 0 # truncate
        self.fid0[self.lenfid//self.trunc:] = 0 # truncate
    
    def full():
        '''
        Uses full Fid.
        '''
        self.fid = self.fid_orig.copy() # retrieve original
        self.fid0 = self.fid_orig.copy() # retrieve original
        
    def prepare_amplitudes(self):
        '''
        Makes the list of amplitudes according to the case. 
        triangular: makes a triangular shape beginning at amplitude "self.amplitude" with slope "self.slope"
        list: retrieves amplitudes from given list of amplitudes.
        '''
        if self.shape == 'triangular':
            self.Amp = [self.amplitude + (i+1)*self.slope for i in range(self.nbpeaks)]    # amplitudes
        elif self.shape == 'list':
            self.Amp = self.amplitude    # list of amplitudes
    
    def prepare_frequencies(self):
        '''
        Frequencies for the signal
        '''
        self.Freq = [(i+1+np.sqrt(10))*pi*500.0j for i in range(self.nbpeaks)]  # frequencies
        
    def make_pure_signal(self):
        '''
        Signal without noise
        '''
        for i in range(self.nbpeaks):
            self.fid0 +=  self.Amp[i] * np.exp(self.Freq[i]*self.x)
            self.fid0 *= np.exp(-self.LB*self.x)
            self.fid0 *= np.exp(1j*self.shift*self.x)
    
    def make_noisy_signal(self):
        '''
        '''
        if self.noisetype == "additive":
            self._make_additive_noise()
        elif self.noisetype == "multiplicative":
            self._make_multiplicative_noise()
        elif self.noisetype == "sampling":
            self._make_sampling_noise()
        elif self.noisetype == "missing points":
            self._make_missing_point_noise()
    
    def _make_additive_noise(self):
        '''
        Makes complex additive noise.
        '''
        real = np.random.randn(self.x.size)
        imag = np.random.randn(self.x.size)
        self.fid = self.fid0 + self.noise*(real + 1j*imag)  # additive complex noise
    
    def _make_multiplicative_noise(self, ampl = True, freq = False):
        '''
        Adding multiplicative noise on amplitude and frequency. 
        '''
        for i in range(self.nbpeaks):
            nAmp = self.Amp[i]
            nFreq = self.Freq[i]
            if ampl:
#                nAmp *= self.mult_Anoise*np.random.randn(self.x.size)
                nAmp += self.mult_Anoise*np.random.randn(self.x.size)
            if freq:
                nFreq += self.mult_Fnoise*np.random.randn(self.x.size) 
            self.fid +=  nAmp * np.exp(nFreq*self.x)
            self.fid *= np.exp(-self.LB*self.x)
            self.fid *= np.exp(1j*self.shift*self.x)
    
    def _make_sampling_noise(self):
        '''
        Takes point with a random step.
        '''
        xn = self.x + 0.5*np.random.randn(self.x.size)/self.lenfid          #  time series with noisy jitter
        for i in range(self.nbpeaks):
            self.fid +=  self.Amp[i] * np.exp(self.Freq[i]*xn) * np.exp(-self.LB*xn) * np.exp(1j*self.shift*self.x)
    
    def _make_missing_point_noise(self):
        '''
        Removes randomly some points after simulation of a regular set of points.
        '''
        miss = np.random.randint(2, size = len(self.x))
        self.fid = self.fid0*miss
    
    def _spec_baseline(self):
        maxi = np.abs(self.mfft(self.fid).max())
        x = np.linspace(0,10,self.lenfid)
        return maxi*(0.1*np.sin(1.5*x)+0.2*x/x.max())
    
    def _spec_from_fid(self, fid):
        '''
        Passing from fid to spectrum
        '''
        spec = self.mfft(fid)
        if self.baseline : spec += self._spec_baseline()
        return spec
        
    @property
    def spec(self):
        '''
        Noisy spectrum
        '''
        return self._spec_from_fid(self.fid)
    
    @property
    def spec0(self):
        '''
        Spectrum wihtout noise
        '''
        return self._spec_from_fid(self.fid0)
    
    def zerof(self, fid, fid_ref):
        '''
        Zerofilling
        '''
        zeropadded = np.concatenate((fid, np.zeros(len(fid_ref)-len(fid))))
        if self.realfid:
            fid_zerf = npirfft(nprfft(zeropadded))
        else:
            fid_zerf = npifft(npfft(zeropadded))
        return fid_zerf
        
    @property
    def iSNR(self):
        '''
        Reference SNR
        '''
        return SNR_dB(self.fid, self.fid0)
    
    def SNR_dB(self, processed_fid):
        '''
        SNR of the denoised signal
        '''
        if processed_fid.size == self.fid0.size:
            fid_ref = self.fid0
        else:
            fid_ref = self.zerof(self.fid0, processed_fid)
            
        return SNR_dB(processed_fid, fid_ref)

def fid_signoise_type(nbpeaks,lendata, noise, noisetype):
    "Obsolete"
    LB = 1.11       # linewidth
    Freq = [(i+1+np.sqrt(10))*pi*500.0j for i in range(nbpeaks)]  # frequencies
    Amp = [(i+1)*20 for i in range(nbpeaks)]    # amplitudes

    data0 = np.zeros(lendata,dtype=complex)
    x = np.arange(lendata*1.0)/lendata          #  time series
    for i in range(nbpeaks):
        data0 +=  Amp[i] * np.exp(Freq[i]*x) * np.exp(-LB*x)

    if noisetype == "additive":
        dataadd = data0 + noise*(np.random.randn(x.size) + 1j*np.random.randn(x.size))  # additive complex noise
        data = dataadd

    elif noisetype == "multiplicative":
        Anoise = noise/2
        Fnoise = noise/200
        for i in range(nbpeaks):
            nAmp = Amp[i] + Anoise*np.random.randn(x.size)
            nFreq = Freq[i] + Fnoise*np.random.randn(x.size)
            data +=  nAmp * np.exp(nFreq*x) * np.exp(-LB*x)
        
    elif noisetype == "sampling":
        xn = x + 0.5*np.random.randn(x.size)/lendata          #  time series with noisy jitter
        data = np.zeros(lendata,dtype = complex)
        for i in range(nbpeaks):
            data +=  Amp[i] * np.exp(Freq[i]*xn) * np.exp(-LB*xn)
        
    elif noisetype == "missing points":
        miss = np.random.randint(2, size=len(x))
        dataadd = data0*miss
        data = dataadd
    else:
        raise Exception("unknown noise type")
    return data

def fid_signoise(nbpeaks, ampl, lengthfid, noise, shift = 0, shape = "triangular", seed = True):
    '''
    Obsolete

    Build fid with triangular spectrum from number of peaks : nbpeaks, amplitude : ampl,
    the length of the fid : lengthfid, and the noise level : noise.
    if shape is "triangular", uses ampl as minimum amplitude.
    if shape is "list", uses the given list to make the amplitudes.
    The seed of random generator can be activated or deactivated with boolean 'seed'
    '''
    if seed:
        np.random.seed(11232)
    LB = 1  # linewidth
    x = np.arange(lengthfid*1.0)/lengthfid        
    fid0 = 1j*np.zeros_like(x)      # complex fid
    omeg = 430.1*1j
    for i in range(1, nbpeaks + 1):
        if shape == "triangular":
            fid0 +=  i*ampl*np.exp(omeg*(i)*x)*np.exp(-LB*x)*np.exp(1j*shift*x)   #
        elif shape == 'list':
            fid0 +=  ampl[i-1]*np.exp(omeg*(i)*x)*np.exp(-LB*x)*np.exp(1j*shift*x)
    fid = fid0 + noise*(np.random.randn(x.size) + 1j*np.random.randn(x.size)) # additive noise 
    return fid

def findnoiselevel(fid, nbseg = 20):
    """
    Routine for determining the noise level from a numpy buffer "fid"
    cut the data in segments then make the standard deviation
    return makes the average on the 25% lowest segments
    
    nbseg=10   nb of segment to cut the spectrum
    """
    less = len(fid)%nbseg     # rest of division of length of data by nb of segment
    restpeaks = fid[less:]   # remove the points that avoid to divide correctly the data in segment of same size.
    newlist = np.array(np.split(restpeaks, nbseg))    #Cutting in segments
    levels = newlist.std(axis = 1)
    levels.sort()
    if nbseg < 4:
        noiselev = levels[0]
    else:
        noiselev = np.mean(levels[0:nbseg//4])
    return noiselev

def findnoiselevel_2D(data_array, nbseg=20, nbrow=200):
    """
    measure noise level on a 2D spectrum
    """
    si1 = data_array.shape[0]
    if si1 > nbrow:
        plist = np.random.permutation(si1)[:nbrow]
    else:
        plist = range(si1)
    levels = [findnoiselevel(data_array[p,:], nbseg=nbseg) for p in plist]
    return np.mean(levels)

def _findnoiselevel_2D(data_array, nbseg = 20):
    """
    obsolete
    """
    dim_array = data_array.size
    x = np.random.randint(dim_array, size = dim_array[0])
    y = np.random.randint(dim_array, size = dim_array[1])
    random2D = zip(x, y)
    noiselev = data_array[random2D].std()

    return noiselev

def findnoiselevel_2D_slice(data_array, ratio = 1):
    """
    Routine for finding noise level by dividing the 2D spectrum in vertical slices. 
    """
    lstd = []
    for i in range(int(data_array.shape[1]*ratio)):
        randcol = np.random.randint(data_array.shape[1])
        cut = data_array[:,randcol]
        noisecut = findnoiselevel(cut, nbseg=30)
        lstd.append(noisecut)
    noiselev = np.array(lstd).mean()
    return noiselev

def findnoiselevel2(fid, nbseg=20):
    """
    (findnoiselevel seems more reliable)
    Routine for determining the noise level from a numpy buffer "fid"
    cut the data in segments then make the standard deviation 
    and returns the smallest deviation

    nbseg=10   nb of segment to cut the spectrum
    """
    less = len(fid)%nbseg     # rest of division of length of data by nb of segment
    restpeaks = fid[less:]   # remove the points that avoid to divide correctly the data in segment of same size.
    newlist = np.array(np.hsplit(restpeaks,nbseg))    #Cutting in segments
    noiselev = newlist.std(axis = 1).min()
    return noiselev

def findnoiselevel_offset(fid, nbseg = 10):
    """
    Routine for determining the noise level and mean value (considering a flat signal with gaussian noise and some outliers in the signal)
    cut the data in segments then make the standard deviation 
    for each segment compares the std dev with the average one and eliminates 
    segments giving a std dev above the mean
    finally makes the average on the sorted segments. 
    nbseg=10   nb of segment to cut the spectrum
    """
    def fun(x, t, data):
        return np.sum(np.abs(x[0] + x[1]* t - data)) # L1 norm for the mean.
    less = len(fid)%nbseg     # rest of division of length of data by nb of segment
    if less !=0:
        restpeaks = fid[:-less]   # remove the points that avoid to divide correctly the data in segment of same size.
    else: 
        restpeaks = fid
    newlist = np.array(np.hsplit(restpeaks,nbseg))    #Cutting in segments
    noiselev = newlist.std(axis=1).min()
    noisemean = newlist.mean(axis=1)
    t = np.arange(len(noisemean))
    res = minimize(fun, x0 = [0,0], args=(t,noisemean), method="Powell") # Find the mean with L1 norm
    return noiselev, res.x[0]

def signal_to_noise(fid, nbseg = 20, peak = "max"):
    """
    determines the signal to noise of a numpy buffer, by dividing le largest peak witht the noiselevel
    see findnoiselevel
    """
    if peak == "max":
        sn = max(fid)/findnoiselevel(fid, nbseg)
    else:
        sn = fid[peak]/findnoiselevel(fid, nbseg)
    return sn

def SNR(dataset, nbseg = 20, peak = "max", unit = None):
    """
    determines the signal to noise of NPKData, by dividing le largest peak witht the noiselevel
    handles complex spectra, taking the noiselevel and the peak on the real part
    noise is evaluated by segmenting the data-set in nbseg peices, check in  findnoiselevel, 
    peak is either max or an index on which to evaluate the signal intensity
    if unit is None, value is in ration, if uint is dB, value is in dB
    """
    noise = findnoiselevel(dataset.buffer[::dataset.axis1.itype+1], nbseg = nbseg)
    if dataset.axis1.itype == 1:
        buff = dataset.buffer[::2]
        try:
            peak = peak/2
        except:
            pass # peak is "max"
    else:
        buff = dataset.buffer
    if peak == "max":
        pk = max(buff)
    else:
        pk = buff[peak]
    snr = pk/noise
    if unit == 'dB':  snr = dB(snr)
    return snr

def dB(val):
    """
    converts a real value to dB
    20*math.log(val)/math.log(10)
    """
    return 20*math.log(val)/math.log(10)

def em_wiener(x,lb):
    "multiply x(t) by exp(-lb*t)"
    t = 1.0*np.arange(len(x))
    return x*np.exp(-t*lb)
    
def apod(lb = 1):
    '''
    apodisation for improving sparsity
    '''
    return em_wiener(x,lb)


class Test_Units(unittest.TestCase):
    '''
    Unittests 
    test_signoise
    test_signoise_truncated
    test_noiselevel
    test_multitests_urQRdtrick
    test_multitests_urQRd_npk
    '''      
    def ttest_signoise(self):
        '''
        Testing the class for generating noisy signal.
        '''
        from ..Display import testplot
        plt = testplot.plot()
        lendata = 10000
        nbpeaks  = 10
        amplitude = 10
        noise = 20
        signoise = SIGNAL_NOISE(lendata, nbpeaks, amplitude, noise)
        plt.plot(signoise.spec)
        plt.show()
    
    def ttest_signoise_truncated(self):
        '''
        Testing the class for generating noisy signal.
        '''
        from ..Display import testplot
        plt = testplot.plot()
        lendata = 10000
        nbpeaks  = 10
        amplitude = 10
        noise = 20
        signoise = SIGNAL_NOISE(lendata, nbpeaks, amplitude, noise, trunc = 10)
        plt.plot(signoise.spec) # spectrum of the full signal
        signoise.cut() # takes the truncated signal with zerofilling
        plt.plot(signoise.spec) # spectrum of the truncated signal
        plt.show()
    
    def ttest_noiselevel(self):
        '''
        Testing noiselevel on experimental data.
        '''
        from ..Tests import filename
        from .. import FTICR 
        import math
        d = FTICR.FTICRData(name = filename("ubiquitine_2D_000002_mr.msh5"))
        e = d.col(11033)
        e.display(show = True)
        for n in (2,5,10,20,50,100,200):
            print(findnoiselevel(e.buffer, n))
        snr = SNR(e)
        print("total SNR: %.0f  %.2f dB"%( snr, 20*math.log(snr)/math.log(10)))
        snr = SNR(e, peak = 453)
        print("pk 453 SNR: %.0f  %.2f dB"%( snr, 20*math.log(snr)/math.log(10)))
    
    def ttest_multitests_urQRdtrick(self):
        '''
        Multitests on urQRd_trick
        '''
        from ..Algo.urQRd_trick import urQRd
        mt = MULTITESTS(algo = urQRd, ampl_sig = 50, ampl_noise = 120.0,\
                    nbpeaks = 6, nb_tries = 1, plot_fft = True)
        algo_param = '2*self.nbpeaks'
        mt.run(algo_param)
    
    def ttest_multitests_urQRd_npk(self):
        '''
        Multitests on urQRd_trick
        '''
        mt = MULTITESTS(algo = 'urqrd', ampl_sig = 50, ampl_noise = 30.0,\
                nbpeaks = 6, nb_tries = 3, npk = True)
        algo_param = 'k = 2*self.nbpeaks, orda = self.lenfid/2'
        mt.run(algo_param)

    def ttest_findnoiselevel_2D_slice(self):  
        arr = np.random.randn(1000,2000)
        arr[25:75,50:150]+=1e9 # Square with high values
        nl = findnoiselevel_2D_slice(arr, ratio=0.1) # noise level
        self.assertAlmostEqual(nl, 0.825, 2)
        
if __name__ == '__main__':
    unittest.main()
    
