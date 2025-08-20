#!/usr/bin/env python
# encoding: utf-8
""""
Fista.py

FISTA: A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems
From Beck and Treboulle 2009. SIAM J. IMAGING SCIENCES. Vol. 2, No. 1, pp. 183–202
Adapted for Python by Lionel and Marc-André on 2013-3-5.
Copyright (c) 2011 IGBMC. All rights reserved.

Solves the LASSO problem : min x {F(x) = ||Ax-b||**2 + lambda||x||_1}
|| ||_1 is the norm l1.
For a given estimated lambda in function of the noise the algorithm tries to converge toward a solution.
If the Chi2 test is ok it stops otherwise the parameter lambda is decreased until to find a correct solution in term of Chi2 test. 
The algorithm is protected for convergence by lambda values always decreasing so as to avoid noise contamination.
"""
from __future__ import print_function, division
import sys
import numpy as np
from math import sqrt
from scipy.linalg import norm
from ..Display import testplot
plt = testplot.plot()
from spike.util.counter import counting, timeit # decorators for timing, counting etc.. 
from spike.util.debug_tools import dec_class_pr, decclassdebugging, debug, pr # debugging
from spike.util.signal_tools import findnoiselevel_offset as findnoiselevel # find the noise level
#libraries for test. 
import unittest

#from util.log_stdout import Logger

__version__ = 0.2
__date__ = "august 2025"

#sys.stdout = Logger()

def hypcpxsqsum(buffer):
    """
    compute the sum of the square of 4 quadrants for hypercomplex 2D data-set, provided as a 2*si1, 2*si2 float matrix 
    return value as si1, si2 floats, 
    """
    a = buffer[::2,::2]**2 + \
                    buffer[1::2,::2]**2 +\
                    buffer[::2,1::2]**2 +\
                    buffer[1::2,1::2]**2
    return a

def hypcpxabs(buffer):
    """
    compute abs() for hypercomplex 2D data-set, provided as a 2*si1, 2*si2 float matrix 
    return value as si1, si2 floats,
    """
    return np.sqrt( hypcpxsqsum(buffer) )

def hypcpxnorm(buffer):
    """
    compute the cartesian norm of hypercomplex 2D data-set, provided as a 2*si1, 2*si2 float matrix 
    return value a float
    """
    return np.sqrt( hypcpxsqsum(buffer).sum() )

def hypcpxabsmask(buffer):
    """
    compute abs() for hypercomplex 2D data-set, provided as a 2*si1, 2*si2 float matrix 
    return abs() as 2*si1, 2*si2 float, where the abs value are in the 4 fields of the hypercpx values
    usefull to compute masking on hypercomplex thresholding.
    """
    aPP = hypcpxabs(buffer)
    a = np.zeros_like(buffer)
    a[::2,::2] = aPP
    a[1::2,::2] = aPP
    a[::2,1::2] = aPP
    a[1::2,1::2] = aPP
    return a

@dec_class_pr
@decclassdebugging
class Fista(object):
    '''
    b = Ax
    b is the result of measurement - real, complex or hypercomplex
    A is the linear transformation for passing from x to b
    x is the solution - real
    ####
    transf : from x to result b, from image to observed
    ttransf : from result b to x, from observed to image. 
    evalnoise : estimation of the noise by the user.

    if hypcpxdatashape is defined, then it should the numpy shape (si1,si2) of the hypcpx data (here as 1D float)
    Then thresholding is performed on the hypcpxabs basis.

    self.lambdaparam : lambda parameter.. regularization parameter.
    iterat : number of iteration for a given lambda. 
    lengthx : length of the vector to be recovered
    prec : precision for convergence for the Chi2 test.
    General structure is : solve() calling fistaloop() calling corealgo()

    '''
    def __init__(self, transf, ttransf, b, scale_noise = 1.0,\
            hypcpxdatashape = None, 
            iterat = 4, miniterat = 300, prec = 0.1, noise_level = None, pow = 1, mfista = False):
        self.transf = transf                                # transformation from the searched vector x to the observed data b
        self.ttransf = ttransf                              # transformation from the data b to searched vector x
        self.b = b                                          # observed data
        self.b_dtype = b.dtype
        self.hypcpxdatashape = hypcpxdatashape
        if hypcpxdatashape is not None:   # if hypercomplex data, modify the norm !
            si1, si2 = hypcpxdatashape
            self.norm = hypcpxnorm
        else:
            self.norm = norm

        temp = ttransf(b)
        print(len(b), len(temp))
        self.x_dtype = temp.dtype
        self.lengthx = len(temp)                             # length of the vector x to be recovered
        del(temp)
#        divide_noiselevel = 1./scale_noise        # divides the noiselevel.                      
        self.prec = prec                                    # precision for the rest of chi2
        self.pow = pow                                      # power for the slight generalization introduced by Rick Chantrand (paper with Voronin 2013)
        ######
        '''
        Initializes Fista parameters
        '''
        self.iternb = 0                                     # counter of iterations in the main loop over labda
        self.iterat = iterat                                # number of iterations for the main loop with lambda fitting. 
        self.miniterat = miniterat                          # number of iterations for the inner loop at lambda constant; if None, loops until
        self.miniteratmax = 100*np.sqrt(b.size)             # maximal number of iterations for Fista for a given lambda 
        #### chi2 Test
        self.chi2 = 1.0                                     # chi2 initialization
        self.limsupstddev = 0.0001                          # upper limit for standard deviation
        #### Fista initialization
        self.xtargm1 = None                                 #
        if self.x_dtype == float:
            self.xtarg = np.ones(self.lengthx)              # initialization of the vector to be recovered.
        elif self.x_dtype == complex:
            self.xtarg = (1+1j)*np.ones(self.lengthx)
        self.yk = self.xtarg                                # initialization of the solution.
        self.tk = 1.0                                       # initialization of the step for the Gradient step.
        self.tkp1 = 1.0
        self.ykp1 = None # 
        #### Lambda
        if noise_level is None:
            self.noise_level = findnoiselevel(abs(self.ttransf(self.b)), nbseg = 11)[0]*scale_noise # The noise level used for finding the first lambda.
        else: 
            self.noise_level = noise_level
        self.recomputlamb = False                           # Recomputes lambda if the chi2 test is negative
        self.nrandip = 5                                    # parameter for calculating the lipschitz constant randomly.
        self.lipsch = self.lipschitz_coeff()                # lipschitz factor from transformation matrix, from Gram matrix.
        self.noiselevelb = self.noise_level/self.lipsch     # noiselevel for the observed data b is equal to noise level of the image x divided by the lipschitz coefficient.
        self.lambdaparam = self.noiselevelb                 # self.lambdaparam initialization to noiselevel of observed data.
        self.t = 0.1/self.lipsch                            # inverse of lipschitz coefficient
        ###### Debugs
        '''
        Initializes the monitoring of chi2 and results.
        '''
        self.residual = None
        self.monitor_chi2 = []                              # list monitoring chi2 on one lambda value, takes the chi2 values for the last lambda value.
        self.monitor_all_final_chi2 = []                    # list for monitoring the final chi2 for each lambda value.
        self.monitor_lambda = []                            # list for monitoring lambda values
        self.all_monitor_chi2 = []                          # list for monitoring the evolution of chi2 for each lambda value.
        self.monitor_result = []                            # list for monitoring the result for each lambda value.
        self.plot_controls = False                          # plot chi2 test and Fista result.
        self.name_plot = ''
        self.best_precision = None
        self.best_yk = self.xtarg                   # Best solution
        self.mfista = mfista                # activate  MFISTA for performing a monotone Fista.

        
    def report(self):
        "dumps content"
        for i in dir(self):
            if not i.startswith('_') and not callable(getattr(self,i)):
                print (i, getattr(self,i))
    
    def F(self, y):
        '''
        Value of the LASSO expression.
        '''
        if self.hypcpxdatashape:
            val =     hypcpxnorm(self.transf(y) - self.b)**2 - self.lambdaparam*(hypcpxabs(y).sum())
        else:
            val =     norm(self.transf(y) - self.b)**2 - self.lambdaparam*abs(y).sum()
        return val 
    
    def lipschitz_coeff(self):
        '''
        Lipfschitz coeff calculated from maximal eigenvalues
        self.lengthx is the length of the observation
        nrandip is the number of trials
        '''
        rr = np.zeros((self.lengthx, self.nrandip), dtype = self.x_dtype)
        for i in range(self.nrandip):                                    # sample N times the trf()
            if self.x_dtype == float:
                r = np.random.randn(self.lengthx)
            elif self.x_dtype == complex:
                r = np.random.randn(self.lengthx) + 1j*np.random.randn(self.lengthx)
            rr[:, i] = self.ttransf(self.transf(r))                      # go back and forth
        rrc = np.dot(np.conj(rr.T), rr)                                  # QQ*
        w, v = np.linalg.eig(rrc)                                        # computes the eigenvalues of the squared matrix
        w = sqrt(max(abs(w)))                                            # then computes largest eigenvalue of (trf o ttrf) 
        return 2*w
    
    def Thresh_shrink(self, x):
        '''
        Thresholding is made on the absolute value abs(x) instead of x as in Fista original algorithm.
        It is an adaptation for complex values.
        for hypercomplex values, see Thresh_shrink_hypcpx
        '''
        a = abs(x)  # this works on real and on complex values - for hypercomplex, see Thresh_shrink_hypcpx
        return np.where(a >= self.lambdaparam, x*(1 - self.lambdaparam/a), 0.0)

    def Thresh_shrink_hypcpx(self, x):
        '''
        Thresholding is made on the absolute value abs(x) instead of x as in Fista original algorithm.
        It is an adaptation for complex values.
        for hypercomplex values

        has to first compute abs() for hypercomplex 2D data-set, provided as a 2*si1, 2*si2 float matrix 
        return abs() as 2*si1, 2*si2 float, where the abs value is copied into the 4 fields of the hypercpx values
        '''
        buffer = x.reshape(self.hypcpxdatashape)
        a = hypcpxabsmask(buffer).ravel()
        return np.where(a >= self.lambdaparam, x*(1 - self.lambdaparam/a), 0.0)
        
    def Thresh_shrink_pow(self, x):
        '''
        Thresholding is made on the absolute value abs(x) instead of x as in Fista original algorithm.
        It is an adaptation for complex values.
        power for the slight generalization introduced by Rick Chantrand (paper with Voronin 2013)
        '''
        a = abs(x)
        return np.where(a >= self.lambdaparam*pow(a, self.pow-1), x*(1 - self.lambdaparam*pow(a, self.pow-1)/a), 0.0)
    
    def chi2_test(self):
        '''
        Performs the test of chi2 on the residual.
        Compares the relative error with the precision defined by self.prec
        Returns the test value.
        '''
        current_precision = abs(self.chi2-float(self.b.size))/float(self.b.size)                 # Current precision for the resisual
        test = (current_precision <= self.prec)                                                  # chi2 test True is current precision is under self.prec.
        if debug(self):
             print ("#### self.chi2, test", self.chi2, test)
        if self.best_precision is None or current_precision < self.best_precision : # test of precision 
             self.best_precision = current_precision
             self.best_yk = self.yk
        return test
    
    def monitor_parameters(self):
        '''
        Follows parameters for each lambda.
        Registers Fista result, chi2 normalized and lambda parameter.
        '''
        self.monitor_all_final_chi2.append(self.chi2/self.b.size)                   # List of the chi2 for each lambda parameter
        self.monitor_result.append(self.yk)                                         # List of Fista results for each lambda parameter
        self.monitor_lambda.append(self.lambdaparam)                                # list of all the lambda parameters used.

    def fista_next_step(self):
        '''
        FISTA updating step
        '''
        self.tkp1 = (1 + sqrt(1 + 4*self.tk**2))/2
        if self.mfista:                                 # Applying MFISTA
            if debug(self):
                print ("self.ztarg.size ",self.ztarg.size )
            if self.F(self.ztarg) < self.F(self.xtargm1):
                self.xtarg  = self.xtargm1
            else:
                self.xtarg  = self.ztarg
            self.ykp1 = self.xtarg + (self.tk)/(self.tkp1)*(self.ztarg - self.xtarg) + (self.tk - 1)/(self.tkp1)*(self.xtarg - self.xtargm1) 
        else:
            self.ykp1 = self.xtarg + (self.tk - 1)/(self.tkp1)*(self.xtarg - self.xtargm1) 
        self.tk = self.tkp1
        self.yk = self.ykp1
        
    def adaptlamb(self):
        '''
        Changes lambda if the chi2 test is wrong ( parameter self.recomputlamb )
        '''
        if debug(self): 
            print ("beginning of adaptlamb ############# : self.chi2/self.b.size ", self.chi2/self.b.size)
            print ("in self.adaptlamb self.recomputlamb ", self.recomputlamb)
        randfact = 1                                      
        if abs(self.chi2/self.b.size-1) > self.prec :     # we are far from 1 
            if debug(self):
                print ("case abs(self.chi2/self.b.size-1) > 0.1")
            if self.recomputlamb : 
                lambfact = self.chi2/self.b.size/2                              # Overshoot we increase strongly lambda again.
                if debug(self):
                    print ("## recompute lambda")
                    print ("## lambfact ", lambfact)
            else :
                if debug(self):
                    print ("##square root for lambfact")
                lambfact = np.sqrt(self.chi2/self.b.size)*randfact              # No overshoot we diminish lambda slowly using square root
        if debug(self):
            print ("ratio self.chi2/self.b.size ", self.chi2/self.b.size)
            print ("decrement lamb ", lambfact)
            print ("self.lambdaparam before ", self.lambdaparam)
        self.lambdaparam /= lambfact                                            # dividing or multiplying factor.
        if debug(self): 
            print ("self.lambdaparam after ", self.lambdaparam)
            
    @timeit
    def solve(self):
        '''
        Repeats fistaloop until abs(self.chi2-float(self.b.size))/float(self.b.size) < self.prec. 
        Other stopping criterium is divergence break
        Adapt the lambda parameter for reaching the noise level.
        '''
        if debug(self):        print ("#################in Fista solve ##########")
        while 1:
            self.iternb += 1
            if debug(self): print ('Iteration ', self.iternb)
            self.fistaloop()                                    # Fista algorithm for a given lambda, if lambda not good go out and make adaptlamb
            self.all_monitor_chi2.append(self.monitor_chi2)    # saves all the monitoring of chi2 with iterations for given lambda.
            if not self.miniterat :                             # if miniterat is None, breaks the main loop, only one lambda used.
                break   
            self.monitor_parameters()                                   #  monitors parameters, list of chi2 for lambda, lambda and list of yk            
            if self.chi2_test() :                               
                if debug(self):
                    print ("self.test_noise break ")
                break                                                   # if chi2 test is ok, stops the algoritm.
            else:
                self.adaptlamb()                                                        # lambda was not good chi2 too high or chi2 too low find another one better.. 
                if self.lambdaparam > 10*self.noiselevelb*self.lipsch:                  # diverging solution case
                    self.yk = np.zeros(self.xtarg.size, dtype=self.x_dtype)             # return a null solution, yk is set equal to 0
                    if debug(self):
                        print ("lambda diverging, not sparse, yk = 0")
                        if debug(self):
                             print ("divergence break ")
                        break                                                                   # divergence, breaks the loop.
            if self.iternb >= self.iterat:
                if debug(self):
                    print ("Exit on maximum iteration number")
                break
        self.residual = self.b-self.transf(self.yk)
        if debug(self):
            print ("self.b.dtype  ",self.b.dtype)
            print ("Done in %d iterations"%self.iternb)
            print ("Final normalized Chi2 :", self.chi2/self.b.size)
            (self.b-self.transf(self.yk)).tofile("residual.npy")
        return self.best_yk                                                                          # return the solution 
    
    def fistaloop(self):
        '''
        Performs Fista for a given lambda.
        The loops breaks if the number of iterations is above max,
        or if the chi2 has a low variation or if its derivative is positive.
        '''
        self.miniternb = 0                                                                      # iteration index inside Fista for a given lambda.
        self.monitor_chi2 = []                                                                  # list for keep track of the chi2                   
        if debug(self):
            print ("now in fistaloop ")
            print ("self.noiselevelb ", self.noiselevelb)
            print ("self.lambdaparam  at beginning of fistaloop", self.lambdaparam)
            print ("chi2 test just before fistaloop  test and result ", self.chi2_test())
        while 1 : 
            if self.miniterat:                                                                  # if the number of iterations is given
                if self.miniternb > self.miniterat:                                             # if miniternb > number of iterations then breaks. 
                    break                                       
            self.recomputlamb = self.corealgo()                                                 # if convergence gives overshoot recomputes lambda. 
            if not self.recomputlamb :
                if debug(self):
                    print ("does not recompute lambda ")
                self.miniternb += 1                                                             #  Counting new iteration of Fista
                # calculates the chi2.
                ir = (self.b-self.transf(self.yk))
                if self.hypcpxdatashape:
                    ir.reshape(self.hypcpxdatashape)
                    self.chi2 = (hypcpxnorm(ir)**2)/self.noiselevelb**2
                else:
                    self.chi2 = (norm(ir)**2)/self.noiselevelb**2          # calculates the chi2.
                if debug(self):
                    print ("self.miniternb, chi2 normalized ", self.miniternb, self.chi2/self.b.size)
                self.monitor_chi2.append(self.chi2/self.b.size)                                 # adds the chi2 ratio of last solution.. 
                ######## Conditions for getting out of Fista loop.
                cnd1 = self.miniternb > self.miniteratmax                                       # self.miniteratmax  is reached
                try:
                    cnd2 = self.chi2_stdv() < self.limsupstddev                                 # variation of the chi2 is too weak
                    cnd3 = self.chi2deriv_max() > 0                                             # the chi2 has a positive derivative.
                except: 
                    cnd2 = False
                    cnd3 = False
                if cnd1 or cnd2 or cnd3: 
                    if debug(self):
                        print ("#### above max, low variation, derivative positive, ", cnd1, cnd2, cnd3)
                    break                                                                       # if nb of iterations > nbconverg or chi2 variation to weak breaks. 
            else :                                                                              # chi2 test is wrong, changing the value of lambda
                if debug(self): print ("go out of fistaloop ")
                break                                                                           # computes solution for other lambda
            
    def corealgo(self):
        '''
        Main part of FISTA algorithm.
        Returns True if lambda is bad and False to contiue if lambda is good.
        '''
        if debug(self):
            print ("now in self.corealgo")
            print ("self.yk.size ", self.yk.size)
            print ("self.b.size ", self.b.size)
        self.xtargm1 = self.xtarg
        interm_residual = self.transf(self.yk) - self.b                    # intermediate residual
        print("interm_residual", interm_residual.dtype )
        if debug(self):
            print ("self.b.size ", self.b.size)
        #     print "interm_residual.size ", interm_residual.size
        # MAD on doit pouvoir localiser self.grad
        grad = self.ttransf(interm_residual)                               # Gradient part, returning to image domain
        if debug(self):
            print ("grad.size ", grad.size)
            print ("self.t ", self.t)
            print ("self.yk.size ", self.yk.size)
        '''
        Intermediate chi2 test
        '''
        new = self.yk - 2.0*self.t*grad                                    # new step before Thresholding
        if self.hypcpxdatashape:
            irhyp = interm_residual.reshape(self.hypcpxdatashape)
            chi2interm = (hypcpxnorm(irhyp)**2)/self.noiselevelb**2        # calculate an intermediate chi2.
        else:
            chi2interm = (norm(interm_residual)**2)/self.noiselevelb**2    # calculate an intermediate chi2.
        if debug(self):
            print ("intermediate chi2 in corealgo", chi2interm/self.b.size)
        if chi2interm/self.b.size > 1 :                                    # if not overshooting continue the calculus 
            if self.mfista:
                # calculation of self.ztarg for MFISTA
                if self.hypcpxdatashape:    # hypercomplex
                    self.ztarg = self.Thresh_shrink_hypcpxmask(new)        # Thresholding since the precedent step gives a good residual with Chi2 test. 
                else:        # real or complex
                    self.ztarg = self.Thresh_shrink(new)                   # Thresholding since the precedent step gives a good residual with Chi2 test. 
            else:
                if self.hypcpxdatashape:    # hypercomplex
                    self.xtarg = self.Thresh_shrink_hypcpxmask(new)         # Thresholding since the precedent step gives a good residual with Chi2 test. 
                else:        # real or complex
                    self.xtarg = self.Thresh_shrink(new)                    # Thresholding since the precedent step gives a good residual with Chi2 test. 
            self.fista_next_step()
            return False                                                    # do not recompute lambda
        else:
            return True                                                     # overshoot, needs to recompute lambda

    def chi2_stdv(self, nbpts = 5):                                         # criterium for convergence
        '''
        Returns the standard deviation on chi2 over few values.
        This value is then compared outside to a up limit.
        '''
        if self.miniternb > nbpts + 1 :
            end_chi2 = np.array(self.monitor_chi2)[-nbpts:]                     # takes the last nblastpoints in monitor_chi2
        return end_chi2.std()/end_chi2.mean() 
    
    def chi2deriv_max(self):                                                    # maximum of the derivative of chi2 values for a given lambda
        '''
        Returns the derivative of evolution of the normalized chi2 value.
        '''
        return max(np.diff(np.array(self.monitor_chi2))) 
    ############# PLots
    def plot_lambda(self):
        '''
        Plots the evolution of lambda
        '''
        print ("self.monitor_lambda", self.monitor_lambda)
        plt.title("evolution of lambda " + self.name_plot)
        plt.semilogy(self.monitor_lambda)                                           # plots values of lambda during Fista execution.

    def plot_chi2(self):
        '''
        Plots the evolution of the normalized chi2 for each lambda
        '''
        for i, fol in enumerate(self.all_monitor_chi2):
            plt.title("evolution of the chi2 for the lambda values " + self.name_plot)

            plt.xlabel("iterations")
            plt.ylabel("$\chi$2")
            plt.semilogy(fol, label = "$\lambda=$" + str(round(self.monitor_lambda[i], 2)))  # plots evolution of chi2 for one lambda value.
        plt.legend()
        
    def plot_residual(self):
        '''
        Plots the residual at the end of Fista.
        '''
        plt.title('residual at the end of Fista ' + self.name_plot)           
        plt.plot(self.residual)                                                     # plots values of the residual during Fista execution.
   
    def plot_result(self, fft_apart = True):
        '''
        Plots the result of Fista in modulus.
        '''
        plt.title('FFT ' + self.name_plot)
        plt.plot(abs(self.ttransf(self.b)), label = 'FFT')                          # plots FFT modulus
        if fft_apart: 
            plt.figure()
        plt.title('Fista ' + self.name_plot)
        plt.plot(abs(self.yk), label = 'FISTA')                                     # plots Fista result modulus
        plt.legend()
    
    def plot_monitorings(self):
        '''
        Controls of Fista parameters and results.
        plots lambda evolution, chi2 evolution, residual and Fista result.
        '''
        self.plot_lambda()                              #
        plt.figure()
        self.plot_chi2()                                                            # Shows the evolution of chi2 for each lambda.
        plt.figure()
        self.plot_residual()                                                        # Shows the residual.
                                                            
class FistaTests(unittest.TestCase):
    '''
    Unittests for Fista with synthetic data.
    '''
    
    def fista_debug(self, fista):
        '''
        Triggering debugs in Fista.
        '''
        fista.adaptlamb_debug = True                                        # debugs the adaptation of lambda parameter.
        fista.solve_debug = True
        fista.corealgo_debug = False                                        # debugs The core algorithm of Fista.
        fista.fistaloop_debug = False                                       # debugs Fista loop
    
    def _test_superresol_synthetic_data_Fista(self):  # ok
        '''
        Fista on synthetic data for testing superresolution.
        '''
        print ("### in test_superresol_synthetic_data_Fista ### ")
        import spike.Algo.CS_transformations as cstr
        from spike.util.signal_tools import fid_signoise
        np.random.seed(1234)
        nbpeaks, lengthx = 20, 20000                                                         # number of peaks and image length.
        lengthb = lengthx//4                                                                 # length of the observed data.
        b = fid_signoise(nbpeaks, ampl = 3, lengthfid = lengthb, noise = 10)               # synthetic noisy Fid of reduced size
        tr = cstr.transformations(lengthx, lengthb, debug = 0)                               # linear transformations for Fista.
        trans = tr.transform                                                                 # Transformation from image to observed.
        ttrans = tr.ttransform                                                               # Transformation from observed to image.
        tr.check()
        ####
        fista = Fista(trans, ttrans, b)                                                      # Instantiation of Fista
        ####
        fista.report()
        fista.solve()
        print ("#### plot results ")
        fista.name_plot = sys._getframe().f_code.co_name
        fista.plot_monitorings()                                                             # Monitoring of Fista parameters
        plt.figure()
        fista.plot_result(fft_apart = True)                                                  # Shows the solution found by Fista
        plt.show()
    
    def _test_superresol_DATA_test_Fista(self):
        '''
        Fista on synthetic data for testing superresolution.
        Thresholding modified according to the paper of Rick Chantrand.
        '''
        print ("### in test_superresol_DATA_test_Fista ### ")
        import spike.Algo.CS_transformations as cstr
        from spike.File.Apex import Import_2D
        from spike.Tests import filename
        from time import time
        ###
        superresol = 4
        f = Import_2D(filename('cytoC_2D_000001.d'))
        f.report()
        f.currentunit = 'm/z'
        sig = f.col(9500).rfft().irfft().buffer
        fsize = sig.size
        tr = cstr.transformations(fsize*superresol, fsize, debug = 0)                               # linear transformations for Fista.
        trans = tr.transform                                                                 # Transformation from image to observed.
        ttrans = tr.ttransform                                                               # Transformation from observed to image.
        tr.check()
        ####
        fista = Fista(trans, ttrans, sig)                                                      # Instantiation of Fista
        ####
        fista.report()
        t0 = time()
        fista.solve()
        t1 = time()
        print ("#### plot results ")
        print ("time elapsed ", t1-t0)
        fista.name_plot = sys._getframe().f_code.co_name
        fista.plot_monitorings()                                                             # Monitoring of Fista parameters
        plt.figure()
        fista.plot_result(fft_apart = False)                                                  # Shows the solution found by Fista
        plt.show()
    
    def _test_superresolution_modif_Thresholding_Fista(self):  #ok
        '''
        Fista on synthetic data for testing superresolution.
        '''
        print ("### in test_superresolution_modif_Thresholding_Fista ### ")
        import spike.Algo.CS_transformations as cstr
        from spike.util.signal_tools import fid_signoise
        np.random.seed(1234)
        nbpeaks, lengthx = 20, 20000                                                         # number of peaks and image length.
        lengthb = lengthx//4                                                                  # length of the observed data.
        b = fid_signoise(nbpeaks, ampl = 3, lengthfid = lengthb, noise = 10)               # synthetic noisy Fid of reduced size
        tr = cstr.transformations(lengthx, lengthb, debug = 0)                               # linear transformations for Fista.
        trans = tr.transform                                                                 # Transformation from image to observed.
        ttrans = tr.ttransform                                                               # Transformation from observed to image.
        tr.check()
        ####
        fista = Fista(trans, ttrans, b)                                                      # Instantiation of Fista
        ####
        fista.pow = 1.05
        fista.report()
        fista.solve()
        print ("#### plot results ")
        fista.name_plot = sys._getframe().f_code.co_name
        fista.plot_monitorings()                                                             # Monitoring of Fista parameters
        plt.figure()
        fista.plot_result(fft_apart = True)                                                  # Shows the solution found by Fista
        plt.show()
        
    def test_sampling_Fista(self):
        '''
        Fista on synthetic data with NUS (Non Uniform Sampling).
        Creates a synthetic Fid, samples it with sampling file at ./DATA_test/Sampling_file.list . 
        Then rebuilds the spectrum with Fista.
        '''
        import spike.Algo.CS_transformations as cstr
        from spike.NPKData import _NPKData 
        print ("### in test_sampling_Fista ### ")
        from ..util.signal_tools import fid_signoise                                                  # import method for producing synthetic data.
        np.random.seed(1234)
        nbpeaks, lengthx = 20, 20000                                                                # number of peaks and length of synthetic data
        b = fid_signoise(nbpeaks, ampl = 3, lengthfid = lengthx, noise = 3)                       # sytnthetic data : noisy Fid
        plt.subplot(311)
        plt.plot(b)                                                                                 # plots the synthetic data. 
        sampling, param = cstr.sampling_load('../DATA_test/Sampling_file.list')                     # retrieves sampling parameters from sampling file.
        print (param)
        plt.subplot(312)
        s = np.zeros((lengthx))
        s[sampling] = 1.0                                                                           # Sampling vector obtained from sampling informations.
        plt.plot(s)                                                                                 
        ####### Linear Transformations for Fista
        tr = cstr.transformations(lengthx, len(sampling), debug = 0)                                # linear transformations for Fista.
        trans = tr.transform                                                                        # Transformation from image to data
        ttrans = tr.ttransform                                                                      # Transformation from data to image
        tr.sampling = sampling                                                                      # insertion of sampling information in the transformations.
        tr.report()
        plt.subplot(313)
        plt.plot(b*s)                                                                               # Plots the synthetic observed data sampled. 
        ######
        c = _NPKData(buffer = b[sampling].copy())                                                    # produces an NPKdata sampled with sampling informations.
        c.axis1.load_sampling('../DATA_test/Sampling_file.list')                                    # associated sampling scheme to data
        c.zf()                                                                                      # zerofill the sampled data
        noise_level = findnoiselevel(abs(ttrans(c.get_buffer())), nbseg = 11)                             # calculates the noiselevel on the zerofilled sampled data.
        #tr.check()
        ####
        fista = Fista(trans, ttrans, b[sampling], noise_level = noise_level)     # Fista NUS
        ####
        fista.report()
        fista.solve()                                                                               # Solves Fista.
        print ("####### plot results ")
        fista.name_plot = sys._getframe().f_code.co_name                                            # Takes the name of the current method.
        fista.plot_monitorings()                                                                    # plots the monitoring of Fista parameters.
        plt.figure()
        fista.plot_result(fft_apart = True)                                                        # Shows the solution found by Fista
        plt.show()
        
if __name__ == '__main__':
    unittest.main()
    
