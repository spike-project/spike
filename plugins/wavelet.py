#!/usr/bin/env python 
# encoding: utf-8

"""A plugin which install wavelet denoising

This plugin is based on the PyWavelet library, which should be installed independently before trying to use this plugin
It can be found at: http://www.pybytes.com/pywavelets/

M-A Delsuc april 2016 from an idea by L Chiron
"""

from __future__ import print_function
import numpy as np
import unittest

try:
    import pywt
    ok = True
except:
    raise Exception('This plugin requires the installation of the PyWavelet library ( www.pybytes.com/pywavelets )')
    print('*** The wavelet plugin requires the installation of the PyWavelet library ( www.pybytes.com/pywavelets )')
    ok = False

from spike import NPKError
from spike.NPKData import NPKData_plugin
from spike.util.signal_tools import findnoiselevel

########################################################################
def ilog2(x):
    "integer log2 definition"
    return int(np.floor(np.log(x)/np.log(2)))
def denoise1D(data, noiseSigma, wavelet='db3'):
    """performed the 1D denoising
    data : a 1D numpy array
    wavelet : the wavelet basis used, 
    
    """
    levels = ilog2(data.shape[0])
    WC = pywt.wavedec(data,wavelet,level=levels)
    threshold=noiseSigma*np.sqrt(2*ilog2(data.size))
    NWC = map(lambda x: pywt.thresholding.soft(x,threshold), WC)
    return pywt.waverec( NWC, wavelet)

def denoise2D(data, noiseSigma, wavelet='db3'):
    """performed the 2D denoising
    data : a 2D numpy array
    wavelet : the wavelet basis used
    
    """
    levels = ilog2(data.shape[0])
    WC = pywt.wavedec2(data,wavelet,level=levels)
    threshold=noiseSigma*np.sqrt(2*ilog2(data.size))
    NWC = map(lambda x: pywt.thresholding.soft(x,threshold), WC)
    return pywt.waverec2( NWC, wavelet)

def wavelet(npkd, nsigma=1.0, wavelet='db3'):
    """
    Performs the wavelet denoising of a 1D or 2D spectrum.
    
    nsigma  the threshold is nsigma times the estimate noise level,
        default 1.0 - corresponds to a relatively strong denoising
    wavelet the wavelet basis used, default 'db3' (Daubechies 3)
        check pywt.wavelist() for the list of possible wavelet
    
    eg:
    d.wavelet(nsigma=0.5)  # d is cleaned after execution
    
    ref: Donoho DL (1995) De-noising by soft-thresholding. IEEE Trans Inf Theory 41:613–621.
    
    Based on the PyWavelet library
    """
    noise = findnoiselevel(npkd.get_buffer())
    if npkd.dim == 1:
        z = denoise1D(npkd.get_buffer(), nsigma*npkd.get_buffer().std(), wavelet=wavelet)
    elif npkd.dim == 2:
        z = denoise2D(npkd.get_buffer(), nsigma*npkd.get_buffer().std(), wavelet=wavelet)
    else:
        raise NPKError("not implemented")
    npkd.set_buffer(z)
    return npkd

if ok:
    NPKData_plugin("wavelet", wavelet)

class WaveLetTest(unittest.TestCase):
    """ - Testing Wavelet plugin- """
    def setUp(self):
        self.verbose = 1    # verbose > 0 switches messages on
    def announce(self):
        if self.verbose >0:
            print("\n========",self.shortDescription(),'===============')
    def test_wave(self):
        """ - testing wavelet - """
        from spike.util.signal_tools import findnoiselevel
        from spike.NPKData import NPKData
        if not ok:
            print("wavelet plugin not installed")
            return
        self.announce()        
        si = 8000
        b = np.random.randn(si)
        d = NPKData(buffer=b)
        d[1234]=100.0
        d[4321]=30.0
        n = findnoiselevel(d.get_buffer())
        self.assertTrue(abs(n-0.95)<0.1)
        self.assertAlmostEqual(d[1234],100.0)
        self.assertAlmostEqual(d[4321],30.0)
        d.wavelet(nsigma=0.5)        # relatively strong denoising
        n = findnoiselevel(d.get_buffer())
        print (n, d[1234], d[4321])
        self.assertTrue(n<0.01)
        self.assertTrue(d[1234]>80.0)
        self.assertTrue(d[4321]>10.0)

