#!/usr/bin/env python 
# encoding: utf-8

"""
maxent.py

Created by Marc-André on 2012-03-06.
Copyright (c) 2012 IGBMC. All rights reserved.
"""
from __future__ import print_function
import sys
import os
import math
import numpy as np
import scipy
import scipy.optimize
import unittest
from ..util.counter import counting, timeit
from ..Display import testplot
plt = testplot.plot()

__version__ = 0.3
__date__ = "march 2013"

########################################################################
# UTILITIES for testing
def br():
    plt.show()
    os._exit(0)
def initial_scene_delta(size = 100,sc=1.0):
    """draw test scene"""
#    Ideal = 0.1*sc*np.ones((size,))
    Ideal = np.zeros((size,))
    Ideal[25]=sc
    return Ideal
#
def initial_scene_2delta(size = 100,sc=1.0):
    """draw test scene"""
#    Ideal = 0.1*sc*np.ones((size,))
    Ideal = np.zeros((size,))
    Ideal[25]=sc
    Ideal[40]=sc
    return Ideal
def initial_scene(size = 100,sc=1.0):
    """draw test scene"""
#    Ideal = 0.1*sc*np.ones((size,))
    Ideal = np.zeros((size,))
    Ideal[15]=5*sc
    for i in range(10):
        Ideal[i+30]=sc
    for i in range (20):
        Ideal[i+65]=i*sc/20.0
    return Ideal
figdef = 0
def plot(buf, text=None, fig=0, logy=False):
    "plot with a text and return the figure number"
    global figdef
    if fig == 0:
        figdef += 1
        lfig = figdef
    else:
        lfig = fig
    plt.figure(lfig)
    if logy:
        plt.semilogy(buf,label=text)
    else:
        plt.plot(buf,label=text)
    if text and fig == 0:
        plt.title(text)
    plt.legend()
    return lfig
#
def estimate_SNR(estim_s, true_s):
    err = true_s - estim_s
    return 10*np.log10(sum(abs(true_s)**2)/sum(abs(err)**2))

################################################################################
class ExpData(object):
    """
    Implement an experimental data to be analysed
    Combined with the definition of the transform ( TransferFunction ), defines the scene for the MaxEnt solver
    """
    def __init__(self, data):
        """
        defines
        data : an array that contains the experimental data
        noise : the value of the noise, in the same unit as data
        window : allows to define error on data for each data point - initialized to 1.0
                error in data[i] is noise/window[i]
        """
        self.data = data
        self.window = np.ones(data.shape)
        self.noise = 0

class TransferFunction(object):
    """ defines the transfert functions for a given problem to be solved by MaxEnt
    data    <-- transform --  image         The experimental transfer function
    data     -- t_transform -->  image      The transpose of `transform'

    must implement
        transform()         // defined above
        t_transform()
        init_transform()
        norm                // the norm of the transform operator, defined as ||tr(x)|| / ||x||
    """
    #-------------------------------
    def __init__(self):
        """
        sets-up the default function, which should be overloaded by a sub classer
        """
        self.init_transform()
    #-------------------------------
    #@counting
    def t_transform(self,datain):
        """
        default method for transform
        here, simple convolution by a 5 pixel rectangle
        """
        return np.convolve(datain, self.convol_kernel, mode="full")     # increases size
    #-------------------------------
    #@counting
    def transform(self,datain):
        """
        default method for transform
        here, simple convolution by a 5 pixel rectangle
        """
        return np.convolve(datain, self.convol_kernel, mode="valid")    # reduces size
    #-------------------------------
    def init_transform(self,size=5):
        """each transform/ttransform pair should have its initialisation code"""
#        self.convol_kernel = (1.0/size)*np.ones(size)
#        self.convol_kernel = np.arange(size)+1.0
#        k = np.asarray([1.0, 2.0, 3.0, 4.0, 3.0, 2.0, 1.0])
#        k = np.asarray([1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0])
#        k = np.asarray([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0])
        k = np.concatenate( (np.arange(float(size)), np.arange(float(size))[::-1]))
        self.convol_kernel = k/np.sum(k)
        self.norm = 1.0

#########################################################
#-------------------------------
def entropy(F, sc):
    """
    compute the entropy of the ndarray a, given scaling value sc
    S =  sum of -Z(i)*LOG(Z(i)), where  Z(i) = a(i)/sc )
    """
    l = np.abs(F) / sc
    S = - np.sum(l*np.log(l))
    return S
#-------------------------------
def entropy_prior(F, A):
    """
    compute the entropy of the ndarray F, given a prior knowledge in ndarray A
    S =  - sum  [F(i) - A(i) + F(i)*LOG(F(i)/A(i)) ]
    """
    l = np.abs(F/A)
    S = - np.sum(F - A + F*np.log(l))
    return S
#-------------------------------
def d_entropy_Gifa(F, S, sc):
    """
    compute the derivative of entropy of the ndarray a, given entropy S and scaling value sc
    Gifa  : dS/dFi = -1/sc (S+log(Fi/A))
    """
    l = np.abs(F) / sc
    return  (-S - np.log(l))/sc
#-------------------------------
def d_entropy_Skilling(F, sc):   # to change - MAD
    """
    compute the derivative of entropy of the ndarray a, given entropy S and scaling value sc
    Skilling  : dS/dFi = -1/sc (1+log(Fi/A))
    """
    l = np.abs(F)/sc
    return (-1.0/sc) * (1.0 + np.log(l))
#-------------------------------
def d_entropy_prior(F, A):   # to change - MAD
    """
    compute the derivative of entropy of the ndarray F, given a prior knowledge in ndarray A
    """
    return (-np.log(F/A))
########################################################
class MaxEnt(object):
    """
    implements the Maximum Entropy solver
    
    given M 
    finds the minimum of Q = S(M) - lamb*Chi2(M)
    with lamb such that Chi2(M) is normalized
    """
    #-------------------------------
    def __init__(self, transferfunction, expdata, iterations=10, debug=0 ):
        """
        initialize a MaxEnt minimiser
        
        sets expdata and transferfunction attributes from arguments
        
        expdata should be a ExpData instance
        transferfunction should be a TransferFunction instance
        
        debug == 0 : production code
        debug == 1 : display evolution
        debug == 2 : full debug        
        
        """
        self.debug = debug
        self.expdata = expdata
        if self.expdata.noise > 0:
            self.sigma = self.expdata.noise
        else:
            self.sigma = 0.01*max(self.expdata.buffer)     # assume a 1% noise level
        self.tf = transferfunction
        self.iterations = iterations
        # then initialisation
        self.setup()
        self.iter = 0               # current iteration counter
        self.lchi2 = []               # will be used to store avancement
        self.lentropy = []
        self.llamb = []
        # initialize image
        tt = self.tf.t_transform(self.expdata.data)
        mm = max(tt)/tt.size
        self.image = mm*np.ones(tt.shape)
        self.sum = np.sum(self.image)
        if self.debug>2:
            plot(self.image,"image initiale")
        self.chi2 = self.do_chisquare(self.image, mode='none')
        self.S = self.do_entropy(self.image, mode='none')
        self.lamb = self.S / self.chi2
        self.true_s = None
    def report(self):
        """report internal state values"""
        for i in dir(self):
            att = getattr(self,i)
            if not callable(att) and not i.startswith('_'):
                print("%s\t:\t%s"% (i, att))
    def setup(self, mode="default"):
        """ 
        Set-up internal MaxEnt parameters so that a kind of optimal state is given
        
        mode can be (case insensitive)
            "Descent"
            "Newton"
            "Gifa"
            "default"   : currently sets to Descent
        modifies (M the MaxEnt instance):
            M.algo M.exptang M.lambdacont M.lambdamax M.lambdasp M.miniter M.stepmax

        """
        if mode.lower() in ('descent', 'default'):
            # convergence control
            self.algo = "steepest"  # convergence algorithm :  steepest  cg  ncg  bfgs  Gifa
            self.stepmax = 2.0      # maximum step along search axis
            self.miniter = 10       # number of iterations for 1D search 
            self.exptang = 1.0          # point at which exp() function switches to linear - used by Gifa algo
            self.chi2target = 1.0       # stopping criterium - stops when normalized chi2 reaches 1.0
            # lambda control
            self.lambdacont = "increment"   # algo for driving lambda :  none  increment  angle  cosine 
            self.lambdasp = 1.1             # how hard we push on lambda
            self.lambdamax = 1.3            # how fast we modify it
        elif mode.lower() == 'gifa' :
            self.algo = "Gifa"
            self.stepmax = 2.0
            self.miniter = 10
            self.exptang = 0.0
            self.chi2target = 1.0
            self.lambdacont = "cosine"
            self.lambdasp = 2.0
            self.lambdamax = 1.5
        else:
            raise (Exception("Wrong mode in MaxEnt.setup()"))
        return
    #-------------------------------
    def do_chisquare(self,image,mode='derivative'):
        """
        computes chi2 (scalar) and its derivative (array),
        computed between     the current image
                    and     self.expdata which is the experimental data
        returned as (chi2, dchi2)
        if mode != 'derivative'  then only chi2 is computed, and returned
        """
        #chi2 = Sum (TF(image)-data)**2)*window / sigma
        residu = (self.tf.transform(image)-self.expdata.data)/self.sigma
        chi2 = np.dot(residu,residu) / residu.size
        if mode == 'derivative':              # then 
            #     dCHI = 2*tTF(TF(image)-data))/sigma
            dchi2 = (2.0/residu.size) * self.tf.t_transform( residu )
        if mode == 'derivative':
            return (chi2, dchi2)
        else:
            return chi2
    def do_entropy(self,image,mode='derivative'):
        """
        computes the entropy S (scalar) and its derivative (array),
        computed on image
        returned as (S, dS)
        if mode != 'derivative'  then only S is computed, and returned
        """
        S = entropy(image, self.sum)
        if mode == 'derivative':            #     compute  derivative of  entropy
            dS = d_entropy_Gifa(self.image, S, self.sum)
            return (S, dS)
        else:
            return S
    #-------------------------------
    #@counting
    def Q(self, im):
        """
        the Q function returns the value of -Q for a given image
        -Q so that it can be minimized
        """
        imclp = im.clip(1.0E-8)
        S = self.do_entropy(imclp, mode="entropy")
        Chi2 = self.do_chisquare(imclp, mode="Chi2")
        Q = S - self.lamb*Chi2
        if self.debug>2: print("in Qfunc     S=%f, Chi2=%f,  Q=%f"%( S, Chi2, Q))
        return -Q
    #-------------------------------
    def sexp(self,x):
        """
        a local and linearized version of the exponential
        exptang is the tangent point
        """
        return np.where(x<self.exptang, np.exp(x), math.exp(self.exptang)*(x+1.0-self.exptang))
    #-------------------------------
    def drive_conv(self,dchi,dS):
        """
        returns (conv, lambda_optim)
        given the derivative dchi and dS, and the current settings
        """
#     convergence test cosine(dS,dC) = dS.dC/(|dS| |dC|)
        if self.iter == 0:
            conv=0.0
            lambda_optim = self.lamb
        else:
            dSsq = np.dot(dS,dS)
            dCsq = np.dot(dchi,dchi)
            dSdC = np.dot(dchi,dS)
            if dSsq*dCsq != 0.0 :   # We'd be in trouble !
                conv = dSdC / math.sqrt(dSsq*dCsq)   # conv is the cos(angle) between dC and dS
                SCratio = math.sqrt(dSsq/dCsq)
            else:
                if self.debug>0: print("*** WARNING in drive_conv, dSsq*dCsq == 0.0")
                conv = 0.0
                SCratio = 1.0
# lambda control
            if self.lambdacont == "none":
                lambda_optim = self.lamb
            elif self.lambdacont == "increment":
                lambda_optim = self.lamb*self.lambdamax
    #if lambcont == "angle"
    #     the best lambda could be such as dQ.dS = 0.0 
    #     i.e. descending Q does not modify the entropy
    #         => lambda = dS dS / dS dC
    #     but requires that conv.gt.0.0  (dC dS angle larger than pi/2)
            elif self.lambdacont == "angle":
                if (conv > 0.0):
                    lambda_optim = dSsq / dSdC
                else:
                    lambda_optim = 0.0
                    print("Warning : Inverse geometry in Angle")
    #     or such as cos(dQ,-dC) = cos(dQ.dS)           if lambcont == "cosine"
            elif self.lambdacont == "cosine":
                lambda_optim = (SCratio*dSdC + dSsq)/(SCratio*dCsq + dSdC)
                sqdSsq = np.sqrt(dSsq)
                sqdCsq = np.sqrt(dCsq)
#                lambda_optim = (sqdCsq*dSsq + sqdSsq*dSdC) / (sqdCsq*dSdC - sqdSsq*dCsq)
            else:
                raise Exception("error in lambdacont value")
            if lambda_optim < 0.0:
                if self.debug>1:
                    print("lambda_optim neg",lambda_optim)
                lambda_optim = 0.0
        return (conv, lambda_optim)
    #-------------------------------
    def update_lambda(self, newlambda, inc=0.5):
        """
        moves self.lamb to new value, using new value newlambda as internal control
        typically : l = l (1-inc) * l inc
        restricted to +/- lambdamax
        """
        if self.iter>0:
            new_lamb = (1.0-inc)*self.lamb + inc*self.lambdasp*newlambda
            self.lamb = min(new_lamb,self.lambdamax*self.lamb)
        else:
            self.lamb = newlambda
        if self.debug>0:
            self.llamb.append(self.lamb)
    #@counting
    def concentring(self):
        "used when no valid step is found"
        # set as a function to allow counting
        self.image *= 0.99
        if self.debug>0: print("concentring")
    #-------------------------------
    @timeit
    def solve(self):
        """
        main entry point for MaxEnt
        after having been initialised, and all attributes set,
        calling solve will realize the MaxEnt convergence.
        
        after having called solve(), the following fields are updated : (see self.report_convergence())
        self.iter   The number of iteration performed so far
        self.image  The final image obtained on the last iteration
        self.chi2   The chi2 obtained on the last iteration
        self.S      The entropy obtained on the last iteration
        self.sum    The integral of self.image
        self.lamb   The lambda obtained on the last iteration
        self.dchi2      The derivative of chi2 obtained on the last iteration
        self.dS         The derivative of S obtained on the last iteration
        the following will be populated only if self.debug>0 :
        self.lchi2      The list of chi2 values since the beginning
        self.lentropy   The list of S values since the beginning
        self.llamb      The list of lambda values since the beginning
        """
        def clp(im):
            "force positivity of im"
            im = im.clip(1e-8)
            return im
        def foptim(x):
            """
            return -Q( image + x.dQ )
            this function get minimized during 1D line-search
            """
            im = self.image + x*dQ
            im = clp(im)
            return self.Q(im)
        def Qprime(im):
            """returns the derivative of -Q"""
            lim = im
            (chi2,dchi) =  self.do_chisquare(lim)
            (S,dS) =  self.do_entropy(lim)
            return self.lamb * dchi - dS
        def Qhess_p(im, p):
            """
            the Qhess_p function returns the result of the matrix product of the hessian of -Q at given image with vector p
            """
            im = np.clip(im,0.001,np.Inf)
            H = (1.0/im) + 2*self.lamb * self.tf.norm**2     # estimate of the Hessian - Cornwell et al algorithm
            return  p/H
        #-------------------------------
        if self.debug >1:
            fig = plot(self.image,'iter 0')
        while (self.iter <= self.iterations and self.chi2 > self.chi2target):
#            print "self.chi2 > self.chi2target",self.chi2, self.chi2target, self.chi2 > self.chi2target
            # prepare
            self.sum = np.sum(self.image)
            print("----------- ITER:",self.iter)
            if self.debug>1:
                print("Sum : %f"%self.sum)
            # compute values and derivatives
            (chi2,dchi) =  self.do_chisquare(self.image)
            self.chi2 = chi2
            (S,dS) =  self.do_entropy(self.image)
            self.S = S # / math.log(self.image.size)
            self.dchi2 = dchi
            self.dS = dS
            if self.debug>0:
                print("S : %f   Chi2 : %f"%(self.S, self.chi2))
                self.lchi2.append(chi2)
                self.lentropy.append(S)
            # compute convergence and lambda
            conv, lambda_optim = self.drive_conv(dchi,dS)
            self.update_lambda(lambda_optim)
            print("convergence : %F   lambda : %f"%(conv,self.lamb))
            step = 0.0
            if self.debug>0:
                print("lambda_optim : %f,    new_lambda : %f"%(lambda_optim, self.lamb))
                Q0 = self.Q(self.image)    # store for later
            # make one step
            if self.algo == "steepest":            #here steepest descent !
                dQ =  self.lamb * dchi - dS
                step =  scipy.optimize.brent(foptim, brack=(0.0, self.stepmax), maxiter=self.miniter, full_output=0)
                if step != 0:
                    self.image = clp(self.image + step*dQ)
                else:
                    self.concentring()
            elif self.algo == "Gifa":               #Gifa is a kind of fixed point
                dQ = self.sum * self.sexp(-self.sum*dchi - self.S) - self.image      # sexp() is an ad-hoc version of exp()
                step =  scipy.optimize.brent(foptim, brack=(0.0, self.stepmax), maxiter=self.miniter, full_output=0)
#                step =  scipy.optimize.fminbound(foptim, x1=0.0, x2=self.stepmax, maxfun=self.miniter, full_output=False)
                if step != 0:
                    self.image = clp(self.image + step*dQ)
                else:
                    self.concentring()
            # These are increasing complexity minimiser. ncg seems to perform best
            elif self.algo == "cg": # conjugated gradient
                self.image = scipy.optimize.fmin_cg(self.Q, self.image, fprime=Qprime, gtol=1e-4, maxiter=self.miniter, disp=0, callback=clp)
            elif self.algo == "bfgs":   # Broyden, Fletcher, Goldfarb, and Shanno
                self.image = scipy.optimize.fmin_bfgs(self.Q, self.image, Qprime, gtol=1e-4, maxiter=self.miniter, disp=0, callback=clp)
            elif self.algo == "ncg":    # newton conjugated gradient
                self.image = scipy.optimize.fmin_ncg( self.Q, self.image, Qprime, fhess_p=Qhess_p, fhess=None, maxiter=self.miniter, disp=1, callback=clp)
            elif self.algo == "slsqp":  # not debuged yet !  Sequential Least-SQuare Programing
                self.image = scipy.optimize.fmin_bfgs(self.Q, self.image, Qprime, gtol=1e-4, maxiter=self.miniter, disp=0, callback=clp)
            else:
                raise Exception('Error with algo')
            # and update
            if self.debug>0:
                print("Q avant %f    Q apres %f    inc %f"%(Q0, self.Q(self.image), step))
            self.iter += 1
            if self.debug >1:
                plot(self.image,"iter %d"%self.iter,fig=fig)
            if self.true_s is not None:
                print('SNR ',estimate_SNR(self.image, self.true_s))
        if self.debug >0:
            plot(self.lentropy,"Entropy")
            plot(self.lchi2, "$\chi^2$",logy=True)
            plot(self.llamb, "$\lambda$ evolution")
        return
    def report_convergence(self):
        """
        returns a set of array, which describe the way convergence was handled
        returns
        self.iter   The number of iteration performed so far
        self.image  The final image obtained on the last iteration
        self.chi2   The chi2 obtained on the last iteration
        self.S      The entropy obtained on the last iteration
        self.sum    The integral of self.image
        self.lamb   The lambda obtained on the last iteration
        self.dchi2      The derivative of chi2 obtained on the last iteration
        self.dS         The derivative of S obtained on the last iteration
        (following are [] if self.debug==0 )
        self.lchi2      The list of chi2 values since the beginning
        self.lentropy   The list of S values since the beginning
        self.llamb      The list of lambda values since the beginning
        """
        if self.iter>0:
            r = [self.iter, self.image, self.chi2, self.S, self.sum, self.lamb, self.dchi2, self.dS, self.lchi2, self.lentropy, self.llamb]
        else:
            r = [self.iter, self.image]
        return r
#-------------------------------
# TESTS

class maxent_Tests(unittest.TestCase):
    def setup(self, size = 100, noise = 0.1, sc = 1):
        "set-up scene"
        # initialize TransferFunction
        T = TransferFunction()
        # create starting image
        Ideal = initial_scene(size,sc)
        Exp = T.transform(Ideal)
        Exp += noise*np.random.randn(Exp.size)
        # initialize ExpData and load data
        D = ExpData(Exp)
        D.noise = noise
        return Ideal, D, T
    def test_sexp(self):
        "draw sexp()"
        Ideal, D, T = self.setup()
        M = MaxEnt(T, D)
        x = np.linspace(-5,5)
        for i in range(5):
            M.exptang = float(i)
            plt.plot(x, M.sexp(x), label = "exptang=" + str(i))
        plt.legend()
    def setup_test(self, M):
        "set-up M fr tests"
        M.algo = "Gifa"
    #    M.chi2target = 0.9
        M.exptang = 0.0
        M.lambdacont = "cosine" #"increment"# "cosine"
    #    M.lamb = 80
        M.lambdamax = 1.5
        M.lambdasp = 2.0
        M.miniter = 100
    def test_MaxEnt(self):
        "run MaxEnt test"
        Ideal, D, T = self.setup(noise = 0.03)
        M = MaxEnt(T, D, debug = 1, iterations = 500)
        M.true_s = Ideal
        self.setup_test(M)
        M.report()
        M.solve()
        f = plot(Ideal,'ideal')   # dummy, just to create a new image
        plt.subplot(411)
        plot(Ideal,'ideal',fig = f)
        plt.subplot(412)
        plot(D.data,'experimental (convoluted, noised)',fig = f)
        plt.subplot(413)
        plot(M.image,'finale',fig = f)
        residu = T.transform(M.image) - D.data
        plt.subplot(414)
        plot(residu, 'residu',fig = f)
        print("Noise in residu :", np.sum(residu**2), M.chi2)
        return  M.report_convergence()
if __name__ == '__main__':
    unittest.main()
