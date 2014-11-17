from util.signal_tools import findnoiselevel, mfft, mrfft
import numpy as np
from numpy.fft import fft
from scipy.linalg import norm
from math import sqrt
import unittest

'''
Module for finding the best rank for urQRd.
Passed parameters are the Fid "fid", the estimated number of lines "estim_nbpeaks" and the order "orda"
'''
class OPTK(object):

    def __init__(self, fid, orda):
        self.norm_show = False                  # Show the evolution of the norm
        self.list_xis = []                      # Subspace signal
        self.list_xin = []                      # Subspace noise
        self.psig = 0                           # signal power
        self.pnoise = 0                         # noise power
        self.fid = fid                          # initial Fid. 
        self.xis = 0                            # initial signal dimension retrieved
        self.xin = 0                            # initial noise dimension retrieved
        self.above_noise = 4                    # correction for noise level.. threshold = noiselevel*above_noise
        self.step_k = 20                        # rank step for calculating the norm
        self.list_urqrd_norm = []               # list of norm values with urQRd.
        self.estim_nbpeaks = None               # estimation of the number of lines.
        self.orda =  orda
        self.prec = 0.9
        
    def subspace_filling_with_rank(self):
        '''
        Estimates subpsaces filling of signal subspace and noise subspace.
        '''
        self.separate_signal_from_noise()									# 
        M = self.orda - self.estim_nbpeaks
        self.psig = norm(self.spec_trunc)**2						                    # Calculation of the signal power
        self.pnoise = norm(self.spec - self.spec_trunc)**2				                # Calculation of the noise power
        for i in range(self.orda):
            empty = (self.psig*(1-self.xis)**2 + self.pnoise*(1 - self.xin)**2)		    # power not yet retrieved.
            self.xis += self.psig*(1-self.xis)**2/empty/self.estim_nbpeaks		        # part of signal dimension retrieved
            self.xin += self.pnoise*(1-self.xin)**2/empty/M				                # part of noise dimension retrieved
            self.list_xis.append(self.xis)
            self.list_xin.append(self.xin)

    def separate_signal_from_noise(self):
        '''
        SIgnal is separated from noise.
        signal is in kept in self.spec
        '''
        if self.fid.dtype == 'complex':
            self.spec = mfft(self.fid)
            print "complex fid"
        elif self.fid.dtype == 'float':
            self.spec = mrfft(self.fid)
            print "real fid"
        self.spec_trunc = self.spec.copy()
        self.noiselev = findnoiselevel(self.spec, nbseg = 10)				            # finds noise level
        self.spec_trunc[self.spec < self.above_noise*self.noiselev] = 0			        # truncates spectrum
        peaks = self.peaks1d(self.spec, threshold = self.above_noise*self.noiselev)
        self.estim_nbpeaks = len(peaks)
        print "self.estim_nbpeaks ", self.estim_nbpeaks
        print "self.above_noise*self.noiselev ", self.above_noise*self.noiselev
        
    def find_best_rank(self):
        '''
        Finds the optimal rank
        '''
        self.subspace_filling_with_rank()                                          	   # Subspaces filling
        diff = abs(np.array(self.list_xis) - self.prec)                            
        #print "diff ", diff
        #print "self.list_xis ", self.list_xis
        minval = (diff).min()
        optk = list(diff).index(minval)
        print "optimal rank is ", optk
        return optk
    
    def peaks1d(self, fid, threshold = 0.1):
        '''
        Extracts peaks from 1d from FID
        '''
        listpk = np.where(((fid > threshold*np.ones(fid.shape))&            # thresholding
                        (fid > np.roll(fid,  1, 0)) &     
                        (fid > np.roll(fid, -1, 0)) ))                      # roll 1 and -1 on axis 0
        return listpk[0]

class optim_urQRd_Tests(unittest.TestCase):  
    def test_optim(self):
        from matplotlib import pyplot as plt
        from util.signal_tools import fid_signoise
        nbpeaks = 15                                                                       # number of peaks
        sigmas = 2                                                                         # amplitude for the peaks
        lengthfid = 2000                                				   # length of the Fid.
        noise = 20                                      				   # white noise amplitude
        fid = fid_signoise(nbpeaks, ampl , lengthfid = lengthfid, noise = noise) 		           # builds the signal
        ########
        orda = lengthfid/4
        optrk = OPTK(fid, orda = orda)                        # optimal rank estimation.
        optk = optrk.find_best_rank()                                             
        print "optk", optk

if __name__ == '__main__':
    unittest.main()
    
