"""
Adaptation of code from :

file              CollombBurg.py
author/translator Ernesto P. Adorio
                  UPDEPP (UP Clark)
                  ernesto.adorio@gmail.com

Version           0.0.1 jun 11, 2010 # first release.
References        Burg's Method, Algorithm and Recursion, pp. 9-11

Created by Lionel on 2011-09-18.
Removed the "for loops" so as to speed up using numpy capabilities.
Copyright (c) 2010 IGBMC. All rights reserved.
"""
import numpy as np
import random
import time
import unittest

def predict(data, ar, length):
    """
    returns a vector with additional points, predicted at the end of "data" up to total size "length", using "ar" polynomial
    """
    
    m = len(ar)
    n = len(data)
    if length<=n:
        raise Exception('values do not allow data extention')
    pred = np.zeros(length)
    pred[0:n] = data[0:n]
    for i in xrange(n,length):
        x = np.sum(- ar * pred[i-1:i-m-1:-1])
        pred[i] = x
    return pred

def denoise(data, ar):
    """
    returned a denoised version of "data", using "ar" polynomial
    first len(ar) points are untouched.
    """
    m = len(ar)
    filtered = data.copy()
    for i  in xrange(m+1, len(data)):
        filtered[i] = predict(data[:i], ar, i+1)[i]
    return filtered

def  burg(m,  x):
    """
    Based on Collomb's C++ code, pp. 10-11
    Burgs Method, algorithm and recursion
      m - number of lags in autoregressive model.
      x  - data vector to approximate.
    """
    N = len(x)-1
    f  = np.array(x[:])
    b = np.array(x[:])
    Ak = np.zeros(m+1); Ak[0] = 1.0 
    Dk= (2.0 * f[:N+1]**2).sum()-(f[ 0 ] **2+ b[ N ] **2)

    for k in range(m):
        mu= (np.roll(f,-(k+1))*b)[:N-k].sum() 
        mu *= -2.0 / Dk
        # update Ak
        t1 = Ak+mu*np.roll(Ak[::-1],k+2)
        t2 = mu*Ak+np.roll(Ak[::-1],k+2)
        limk = (k + 1)/2
        limkr = int(round((k + 1)/2.0))
        Ak[:limk+1] = t1[:limk+1]
        Ak[limk+1:k+2] = np.roll(t2[::-1],limkr)[:len(Ak[limk:k+1])]
        #update f and b
        t1 = (np.roll(f,-(k+1))+mu*b)[:N-k]
        t2 = (b+mu*np.roll(f,-(k+1)))[:N-k]
        f[k+1:N+1] = t1
        b[:N-k] = t2
        Dk = ( 1.0 - mu **2 ) * Dk - f[ k + 1 ] **2 - b[ N - k - 1 ] **2
    return Ak[1:]

class LinpredTests(unittest.TestCase):
    """ - Testing linear prediction , Burg algorithm- """
    def setUp(self):
        self.verbose = 1    # verbose > 0 switches messages on
    def announce(self):
        if self.verbose >0:
            print "\n========",self.shortDescription(),'==============='
    def curvetest(self):
        longtest=128
        x = np.arange(longtest)
        original = np.cos( x * 0.01 ) + 0.75*np.cos( x * 0.2 ) + 0.5*np.cos( x * 0.05 ) + 0.25*np.cos( x * 0.11 ) + 0.1*np.random.randn(longtest)
        return original

    def test_burg(self):
        """ - testing burg algo - """
        import Display.testplot as testplot
        plt = testplot.plot()
        self.announce()        
        mlength=12      # ideal is 2*number of signal - here should 8, 12 gives "a little room"" to noise, and seems better
        original=self.curvetest()
        plt.plot(original,label='noisy original')
        coeffs = burg(mlength,  original)
        m = len(coeffs)
        predicted = predict(original, coeffs, 180)
        self.assertEqual(180, len(predicted))
        plt.plot(predicted, label='extended')
        denoised = denoise(original, coeffs)
        self.assertTrue( np.sum(abs(original-denoised)) / np.sum(abs(original)) <0.2 )
        plt.plot(denoised, label='denoised')
        plt.legend()
        plt.show()

if __name__ == "__main__":
    unittest.main()