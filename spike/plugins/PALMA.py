#!/usr/bin/env python 
# -*- coding: utf-8 -*-
"""complete DOSY processing, using the PALMA algorithm

This program uses the PALMA algorithm, presented in the manuscript

Cherni, A., Chouzenoux, E., & Delsuc, M.-A. (2017).
PALMA, an improved algorithm for DOSY signal processing.
Analyst, 142(5), 772–779. http://doi.org/10.1039/c6an01902a

see manuscript for details.

Authors:
Afef Cherni, Emilie Chouzenoux, and Marc-André Delsuc

Licence: CC-BY-NC-SA   https://creativecommons.org/licenses/by-nc-sa/4.0/
"""
from __future__ import print_function, division

import sys
import os.path as op
import unittest
import re

import numpy as np
import scipy

from spike import NPKError
from spike.NPKData import NPKData_plugin
from spike.NPKData import LaplaceAxis
from spike.File.BrukerNMR import Import_2D, Import_2D_proc
import spike.Algo.savitzky_golay as sgm
import spike.util.signal_tools

debug = True

#version = 1.0
# february 2018 - added an randomisation of colonne processing for a cleaner progress bar
version = 1.1

################# PPXA+ Algo ######################
def residus(x, K, y):
    "computes distance between y and Kx"
    return np.linalg.norm(np.dot(K,x)-y,2)
def L1(x):
    "computes L1 norm of x"
    absx = np.abs(x)
    return np.sum(absx)
def ent(x,a):
    "computes Entropy of x"
    xpos = x[x>0]/a           #if x = 0: rent sera infini (log (0)) ==> erreur!
    return -np.sum(xpos*np.log(xpos)) 

def criterion(x, K, y, lamda, a):
    """
    Compute regularization function, (not used during iteration without full_output)
    """
    f = 0.5*residus(x,K,y)**2
    if (1-lamda) > 0:
        rl1 = (1-lamda)*L1(x)
    else:
        rl1 = 0
    if lamda > 0:
        rent = -lamda*ent(x,a)
    else:
        rent = 0
    return f + rl1 + rent

def lambert_w(x):
    """
    W Lambert function
    """
    d = scipy.special.lambertw(x, k=0, tol=1e-8)
    return np.real(d)

def approx_lambert(x):
    """
    approximation of W( exp(x) )
    no error below 50, and  less than 0.2% error for x>50 and converges toward 0 for x->inf
    does not overflow !     does not NaN
    """
    limit = 50
    nz = np.nonzero( (x < limit) )[0]
    A = 0.00303583748046
    s = x*(1 - np.log(x)/(1+x)) + A/(1+x-limit)**0.75
    s[nz] = lambert_w(np.exp(x[nz]))
    return np.nan_to_num(s)

# def approx_lambert(x):
#     """ 
#     Lambert approximation of exp(x)
#     """
#     limit = 100
#     nz = np.nonzero(x < limit)[0]
#     s = x - np.log(x)
#     s[nz] = lambert_w(np.exp(x[nz]))
#     return s 

def prox_l1(x, w) :
    """
    Compute proximity operator of L1 norm"
    """
    p = np.zeros_like(x)
    pos = np.nonzero(x>w)[0]
    p[pos] = x[pos] - w
    neg = np.nonzero(x<-w)[0]
    p[neg] = x[neg] + w
    return p

def prox_l2(x, dx, eta) :
    """
    Compute projection of x onto l2 ball ||z-dx||<=eta
    x and dx are image vectors  
    """
    t = x-dx
    s = t*np.minimum(eta/np.linalg.norm(t),1)
    return x + s - t

def prox_l1_Sent(x, lamda, a):
    """
    Compute the proximity operator of L1 + Shannon Entropy
    """    
    if lamda == 0:
        p = prox_l1(x, 1)
    else:
        loga = np.log(a)
        loglamda  = np.log(lamda)
        c = (a*x - a*(1-lamda))/lamda - 1 - loglamda + 2*loga
        p = (lamda/a)*approx_lambert(c)
    return p

debug=False
def PPXAplus(K, Binv, y, eta, nbiter=1000, lamda=0.1, prec=1E-12, full_output=False):
    r"""
    performs the PPXA+ algorithm
    K : a MxN matrix which transform from data space to image space
    Binv : inverse of (Id + K.t K)
    y : a M vector containing the data
    a : an estimate of $\sum{x}$ where x is final image - used as a bayesian prior of x
    eta : an estimate of the standard deviation of the noise in y
    nbiter : maximum number of iteration to do
    lamda: is in [0..1], the weigth of l1 vs MaxEnt regularisation
        lamda = 0 is full L1
        lamda = 1 is full MaxEnt
    prec: precision of the result, algo will stop if steps are below this evel
    full_output: if True, will compute additional terms during convergence (slower):
        parameters =  (lcrit, lent, lL1, lresidus)
            with lcrit: the optimized criterion
            len: evolution of -entropy
            lL1: evolution of L1(x)
            lresidus: evolution of the distance ||Kx-y||
        if False, returns the number of performed iterations
    returns
    (x, parameters), where x is the computed optimal image
    
    """

    # 1 - scaling step
    scale = y[0]
    y = y / scale
    eta = eta / scale
    a = y[0]   
    # 2 - PPXA+ 
    # 2.1 preparation
    gamma = 1.99
    M,N = K.shape
    x0 = np.ones((N,1))
    x0 *= np.sum(y) / (M*N) #x0 = y[0] /N
    lcrit = []
    lent = []
    lL1 = []
    lresidus = []
    Kt = K.T
    x_n_old = x0.copy()
    tmp1 = x0.copy()
    tmp2 = np.dot(K,x0)
    x_n = np.dot(Binv,tmp1 + np.dot(Kt,tmp2))
    # 2.2 loop
    n = 0
    for n in range(0,nbiter):
        xx1 = prox_l1_Sent(tmp1, lamda, a)   # L1 + Shannon
        xx2 = prox_l2(tmp2, y, eta)
        c = np.dot(Binv, xx1 + np.dot(Kt,xx2))
        cmxn = c - x_n
        c2mxn = c + cmxn
        tmp1 += gamma*(c2mxn - xx1)
        tmp2 += gamma*(np.dot(K, c2mxn) - xx2)
        x_n += gamma*cmxn
        
        n_x_n = np.linalg.norm(x_n-x_n_old,2) / np.linalg.norm(x_n)
        if np.isnan( x_n.sum() ):  # Nan appear sometimes in pathological cases
            break
        if n_x_n < prec:
            break
        x_n_old[:,:] = x_n[:,:]
        if full_output is True:
            lcrit.append(criterion(x_n, K, y, lamda, a))
            lent.append(ent(x_n,a))
            lL1.append(L1(x_n))
            lresidus.append(residus(x_n,K,y))
    #  3 - eliminate scaling step
    x_n = x_n * scale
    if full_output:
        comp = [lcrit, lent, lL1, lresidus]
    else:
        comp = n        # number of iterations
    return x_n, comp

def eval_dosy_noise(x, window_size=9, order=3):
    """
    we estimate the noise in x by computing difference from polynomial fitting
    input: x - a real vector
    return: the noise level
    """
    m = sgm.sgolay_coef(window_size, order, deriv=0)
    noise = (sgm.sgolay_comp(x, m, window_size)- x).std()

    return noise

############# DOSY set-up ###########################
def auto_damp_width(d):
    """
    uses the tab buffer to determine the optimum dmin and dmax for ILT processing
    """
    import math
    mn = min(d.axis1.qvalues)
    mx = max(d.axis1.qvalues)

    idmax = (d.axis1.dfactor*2/mn)
    idmin = (d.axis1.dfactor*0.5/mx)
    logidmin = (int(math.log(idmin)/math.log(10)))
    edmin = (math.exp(math.log(10)*logidmin))
    dmin = (int((int(idmin/edmin))*edmin))

# calculate dmax from idmax by taking it's upper 10th decimal
    logidmax = (int(math.log(idmax)/math.log(10)))
    edmax = (math.exp(math.log(10)*logidmax))

    dmax = (int(((int(idmax/edmax))+1)*edmax))

    return dmin, dmax

#----------------------------------------------------
def calibdosy(litdelta, bigdelta, recovery=0.0, seq_type='ste', nucleus='1H', maxgrad=50.0, maxtab=50.0, gradshape=1.0, unbalancing=0.2, os_tau=None, os_version=1):
    """
    returns the DOSY calibrating factor from the parameters
    
      bigdelta float
       "Big Delta"  : diffusion delay in msec

      litdelta float
       "little Delta"  : gradient duration in msec

      seq_type enum "pgse","ste","bpp_ste","ste_2echoes","bpp_ste_2echoes","oneshot" / default ste
       the type of DOSY sequence used
        pgse : the standard hahn echoe sequence
        ste : the standard stimulated echoe sequence
        bpp_ste : ste with bipolar gradient pulses
        ste_2echoes : ste compensated for convection
        bpp_ste_2echoes : bpp_ste compensated for convection
        oneshot : the oneshot sequence from Pelta, Morris, Stchedroff, Hammond, 2002, Magn.Reson.Chem. 40, p147
            uses unbalancing=0.2, os_tau=None, os_version=1
            unbalancing is called alpha in the publication
            os_tau is called tau in the publication
            os_version=1 corresponds to equation(1) / os_version=2 to (2)
      nucleus enum "1H","2H","13C","15N","17O","19F","31P" / default 1H
       the observed nucleus

     recovery float
       Gradient recovery delay

     maxgrad float
       Maximum Amplificator Gradient Intensity, in G/cm    / default 50.0

     maxtab float
       Maximum Tabulated Gradient Value in the tabulated file. / default 50.0
       Bruker users with gradient list in G/cm (difflist) use maxgrad here
       Bruker users with gradient list in % use 100 here
       Varian users use 32768 here

     gradshape float
       integral factor depending on the gradient shape used / default 1.0
       typical values are :
           1.0 for rectangular gradients
           0.6366 = 2/pi for sine bell gradients
           0.4839496 for 4% truncated gaussian (Bruker gauss.100 file)
       Bruker users using difflist use 1.0 here, as it is already included in difflist

    """
# MAD : modified august-sept 2007 - corrected ste; added oneshot; added PGSE
# MAD : modified june 2018 - corrected oneshot


    g = (maxgrad / maxtab)*1E-4 # now in Tesla/cm
    aire = g*gradshape*litdelta
    if nucleus == "1H":
        gama = 2.675E8                           # rad.s-1.T-1
    elif nucleus == '2H':
        gama = 0.411E8                           # rad.s-1.T-1
    elif nucleus =='13C':
        gama = 0.673E8                           # rad.s-1.T-1
    elif nucleus == '15N':
        gama = -0.271E8                          # rad.s-1.T-1
    elif nucleus == '17O':
        gama = -0.363E8                          # rad.s-1.T-1
    elif nucleus =='19F':
        gama = 2.517E8                           # rad.s-1.T-1
    elif nucleus =='31P':
        gama = 1.083E8                           # rad.s-1.T-1
    else:
        raise Exception('Unknown nucleus')

    K = ((gama * aire)**2)            # Calcul de q^2

# equation references are in    Jerschow,A.;Muller,N.;JMR;125;1997;Suppresion of convection artifacts
    if seq_type == 'ste':
        K = (K * (bigdelta + ((2 * litdelta)/3) + recovery))                 # cm 2 sec-1 pour Q
    elif seq_type == 'bpp_ste':
        K = (K * (bigdelta + ((2 * litdelta)/3) + ((3 * recovery)/4))) # cm 2 sec-1 pour Q
    elif seq_type == 'ste_2echoes':
        K = (K * (bigdelta + ((4 * litdelta)/3) + (2 * recovery)))     # cm 2 sec-1 pour Q
    elif seq_type == 'bpp_ste_2echoes':
        K = (K * (bigdelta + ((4 * litdelta)/3) + ((3 * recovery)/2))) # cm 2 sec-1 pour Q
    elif seq_type == 'oneshot':
        if os_version == 1:
            K = (K * (bigdelta + litdelta * (unbalancing**2 - 2) / 6 + os_tau * (unbalancing**2 - 1) / 2))
        elif os_version == 2:
            K = (K * (bigdelta + litdelta * (unbalancing**2 + 3*unbalancing - 2) / 6 + os_tau * (unbalancing**2 + 2*unbalancing - 1) / 2))
        else:
            raise Exception('Unknown sequence: '+str(seq_type)+str(os_version))
    elif seq_type == 'pgse':
        K = (K * bigdelta + (2 * litdelta)/3)
    else:
        raise Exception('Unknown sequence: '+str(seq_type))

    K = (K * 1e-8)      # from cm^2 to um^2
    return(1/K)

def determine_seqtype(pulprog):
    """
    given the PULPROG name, determines which seq_type is to be used
    PULPROG should be follow the standard Bruker naming scheme
    """
# Bruker avance sequences
    if (re.search('dstebp',pulprog)):
        sequence = 'bpp_ste_2echoes'
    elif re.search('dstegp',pulprog):
        sequence = 'ste_2echoes'
    elif re.search('stebp|stegpbp|ledbp',pulprog):
        sequence = 'bpp_ste'
    elif re.search('stegp|led',pulprog):
        sequence = 'ste'
    elif re.search('oneshot',pulprog):
        sequence = 'oneshot'
    else:
        print("<%s> : Unsupported pulse program."%pulprog)
        sequence = "None"
    print (sequence)
    return sequence

def dcalibdosy(npk, nucleus='1H'):
    """use stored parameters to determine correct DOSY calbiration"""
    d20 = float(npk.params['acqu']['$D'][20])
    d16 = float(npk.params['acqu']['$D'][16])
    d17 = float(npk.params['acqu']['$D'][17])
    p1 =  float(npk.params['acqu']["$P"][1])*1e-6
    p19 = float(npk.params['acqu']["$P"][19])*1e-6
    p30 = float(npk.params['acqu']['$P'][30])*1e-6

    nuc1 = npk.params['acqu']["$NUC1"]
    if nucleus is None:
        if (nuc1 == '1H' or nuc1 == '15N' or nuc1 == '13C' or nuc1 == '31P' or nuc1 == '19F' or nuc1 == '17O'):
            nucleus = nuc1
        else:
            nucleus = '1H'
    print ("DOSY performed on %s"%(nucleus,))
    pulprog = npk.params['acqu']['$PULPROG']
    seq_type = determine_seqtype(pulprog[1:-1])
    
    # STEBP_2echos Bruker avance sequences
    if seq_type == 'bpp_ste_2echoes':
        litdelta = (2*p30)
        bigdelta = (d20-(10*p1)-(8*p30)-(8*d16)-(8*d17)-(2*p19))
        recovery = d16
    # STE_2echos Bruker avance sequences
    elif seq_type == 'ste_2echoes':
        litdelta = p30
        bigdelta = (2*(d20-(2*p1)-(p30)-(2*d16)-(p19)))
        recovery = d16
    # BPP_LED NMRtec and Bruker Avance sequences
    elif seq_type == 'bpp_ste':
        litdelta = 2*p30
        bigdelta = d20-(4*p1)-(2*p30)-(3*d16)-(p19)
        recovery = 2*d16
    # LEDgp/STEgp Bruker Avance sequence
    elif seq_type == 'ste':
        litdelta = p30
        bigdelta = d20-(2*p1)-(p30)-(2*d16)-(p19)
        recovery = d16
    #Doneshot from Morris and Nillson
    elif seq_type == 'oneshot':
        litdelta = 2*p30
        bigdelta = d20-(4*p1)-(2*p30)-(3*d16)-(p19)
        recovery = 2*d16
    else:
        litdelta = p30
        bigdelta = d20
        recovery = d16
    
    print (litdelta, bigdelta, recovery, seq_type, nucleus)
    npk.axis1.dfactor = calibdosy(litdelta, bigdelta, recovery, seq_type=seq_type, nucleus=nucleus)
    return npk

def Import_DOSY(fname, nucleus=None, verbose=False):
    """
    Import and calibrate DOSY data-set from a Bruker ser file
    """
    d = Import_2D(fname)
    d.axis1 = LaplaceAxis(size=d.size1)
    dire=op.dirname(fname)
    d.axis1.load_qvalues(op.join(dire,"difflist"))
    if d.axis1.size != len(d.axis1.qvalues):
        l = min(d.axis1.size, len(d.axis1.qvalues))
        print ("WARNING in Import_DOSY(), size missmatch data is %d while difflist is %d"%(d.axis1.size, len(d.axis1.qvalues)))
        print("truncating to %d"%(l,))
        d.chsize(sz1=l)
        d.axis1.qvalues = d.axis1.qvalues[:l]
    d.calibdosy()
    if verbose:
        print("imported 2D DOSY, size = %d x %d\n%s"%(d.axis1.size, d.axis2.size, d.params['acqu']['title']))
    return d

def Import_DOSY_proc(fname, nucleus='1H', verbose=False):
    """
    Import and calibrate DOSY data-set from a Bruker 2rr file
    """
    d = Import_2D_proc(fname)
    d.axis1 = LaplaceAxis(size=d.size1)
    dire=op.dirname(op.dirname(op.dirname(fname)))   # up to expno
    d.axis1.load_qvalues(op.join(dire,"difflist"))
    # if d.axis1.size != len(d.axis1.qvalues):
    #     l = min(d.axis1.size, len(d.axis1.qvalues))
    #     print ("WARNING in Import_DOSY_proc(), size missmatch data is %d while difflist is %d"%(d.axis1.size, len(d.axis1.qvalues)))
    #     print("truncating to %d"%(l,))
    #     d.chsize(sz1=l)
    #     d.axis1.qvalues = d.axis1.qvalues[:l]
    d.calibdosy()
    # In Topspin, the diffusion axis parameters are faked as ppm - lets decode them from proc2
    Dmax = 10**(float(d.params['proc2']['$OFFSET'])+12)  # 12 to change from m2/s to µm2/s
    width = float(d.params['proc2']['$SW_p'])/float(d.params['proc2']['$SF'])
    Dmin = Dmax*(10**(-width))
    d.reverse(axis=1)
    d.axis1.dmin = Dmin
    d.axis1.dmax = Dmax
    if verbose:
        print("imported 2D DOSY spectrum, size = %d x %d\n%s"%(d.axis1.size, d.axis2.size, d.params['acqu']['title']))
    return d

#################### PALMA setup ###########################
def process(param):
    " do the elemental processing, used by loops"
    icol, c, N, valmini, nbiter, lamda, precision, uncertainty = param
    if (c[0] > valmini):
        y = c.get_buffer()
        c = c.palma(N, nbiter=nbiter, lamda=lamda, precision=precision, uncertainty=uncertainty)
        lchi2 = np.linalg.norm(y-np.dot(c.axis1.K,c.get_buffer()))
    else:
        c = c.set_buffer(np.zeros(N))
        lchi2 = 0
    return (icol, c, lchi2)

def do_palma(npkd, miniSNR=32, mppool=None, nbiter=1000, lamda=0.1, uncertainty=1.2, precision=1E-8):
    """
    realize PALMA computation on each column of the 2D datasets
    dataset should have been prepared with prepare_palma()
    
    the noise in the initial spectrum is analysed on the first DOSY increment
    then each column is processed with palma() if its intensity is sufficient

    miniSNR: determines the minimum Signal to Noise Ratio of the signal for allowing the processing
    mppool: if passed as a multiprocessing.Pool, it will be used for parallel processing
    
    the other parameters are transparently passed to palma()

    """
    import multiprocessing as mp
    from spike.util import progressbar as pg
    from spike.util import widgets
    import sys #
    if sys.version_info[0] < 3:
        import itertools
        imap = itertools.imap
    else:
        imap = map

    # local functions
    def palmaiter(npkd):
        "iterator for // processing around palma() using mp.pool.imap()"
        #for c in npkd.xcol(): #
        for icol in np.random.permutation(npkd.size2):  # create a randomized range
            c = npkd.col(icol)
            yield (icol, c, N, valmini, nbiter, lamda, precision, uncertainty)
        
    # prepare
    if mppool is not None:
        if isinstance(mppool, mp.pool.Pool):
            paral = True
        else:
            raise Exception("parameter mpool should be either None or of multiprocessing.Pool type")
    else:
        paral = False
    npkd.check2D()
    K = npkd.axis1.K
    M,N = K.shape
    output = npkd.copy()
    output.chsize(sz1=N)
    chi2 = np.zeros(npkd.size2)   # this vector contains the final chi2 for each column
    noise = spike.util.signal_tools.findnoiselevel(npkd.row(0).get_buffer())
    valmini = noise*miniSNR
    # loop
    xarg = palmaiter(npkd)
    wdg = ['PALMA: ', widgets.Percentage(), ' ', widgets.Bar(marker='-',left='[',right=']'), widgets.ETA()]
    pbar= pg.ProgressBar(widgets=wdg, maxval=npkd.size2).start() #, fd=sys.stdout)
    if paral:
        result = mppool.imap(process, xarg)
    else:
        result = imap(process, xarg)
    # collect
    for ii, res in enumerate(result):
        # if icol%50 == 0 :
        #     print ("DOSY # %d / %d"%(icol,npkd.size2))
        pbar.update(ii+1)
        sys.stdout.flush()
        icol, c, lchi2 = res
        chi2[icol] = lchi2
        output.set_col(icol, c)                                                                                                         
    
    # for icol in range(npkd.size2):
    #     #if icol%10 ==0: print (icol, "iteration")
    #     c = npkd.col(icol)
    #     if (c[0] > miniSNR*noise):
    #         y = c.get_buffer()
    #         c = c.palma(N, nbiter=nbiter, lamda!!, precision=precision, uncertainty=uncertainty)
    #         chi2[icol] = np.linalg.norm(y-np.dot(K,c.get_buffer()))
    #     else:
    #         c = c.set_buffer(np.zeros(N))
    #     result.set_col(icol, c)                                                                                                         
    pbar.finish()
    output.axis1.chi2 = chi2

    return output


def prepare_palma(npkd, finalsize, Dmin, Dmax):
    """
    this method prepares a DOSY dataset for processing
    - computes experimental values from imported parameter file
    - prepare DOSY transformation matrix
    """
    npkd.check2D()
    M = npkd.size1
    N = finalsize
    # computes t / direct space sampling
    t = npkd.axis1.qvalues**2
    t /= npkd.axis1.dfactor
    t = t.reshape((M,1))
    # compute T / Laplace space sampling
    npkd.axis1.dmin = Dmin
    npkd.axis1.dmax = Dmax
    targetaxis = LaplaceAxis(size=N)
    targetaxis.dmin = Dmin
    targetaxis.dmax = Dmax
    T = targetaxis.itod( np.arange(N) )
    T = T.reshape((1,N))
    K = np.exp(-np.kron(t, T))
    npkd.axis1.K = K

    #Stepsize parameter 
    Kt = np.transpose(K)
    KtK = np.dot(Kt,K)
    
    #Inversion of linear operator
    B = np.identity(N)
    B = B + KtK
    Binv = np.linalg.inv(B)
    npkd.axis1.Binv = Binv

    return npkd


def palma(npkd, N, nbiter=1000, uncertainty=1.0, lamda=0.1, precision=1E-8, full_output=False):
    """
    realize PALMA computation on a 1D dataset containing a decay
    dataset should have been prepared with prepare_palma
    on each column, the noise is estimated, then PPXA+ algorithm is applied
    nbiter: maximum iteration number of the PALMA algo
    
    uncertainty: noise is estimated on the dataset, and then multiplied by this value
        so uncertainty=1 is full confidence in the noise evaluation algo
        uncertainty>1 allows more room for algo when data quality is poor
    lamda is the balance between entropy and L1
        lamda = 0 is full L1
        lamda = 1 is full Ent
    
    precision: is the required precision for the convergence
    full_output is used for debugging purposes, do not use in production
        check PPXAplus() doc for details
    """
    NaN_found = 0
    npkd.check1D()
    K = npkd.axis1.K            # the measure matrix
    Binv = npkd.axis1.Binv
    y = npkd.get_buffer()       # the mesasured values
    M, Nk = K.shape
    if debug:
        print (npkd.size1, M, Nk ,N)
    if (npkd.size1 != M) or (Nk != N):
        raise Exception("Size missmatch in palma : %d x %d  while data is %d x %d" % (M, Nk, npkd.size1, N))

    evald_noise = eval_dosy_noise(y)
    y = y.reshape((M,1))
    Ok = False
    while not Ok:  # this is to force positivity or not NaN; 
        eta = uncertainty*np.sqrt(M)*evald_noise
        if debug:
            print(" noise: %f  uncertainty: %f  eta: %f"%(evald_noise,uncertainty,eta))
        x, c = PPXAplus(K, Binv, y, eta, nbiter=nbiter, lamda=lamda, prec=precision, full_output=full_output)
        Ok = not np.isnan( x.sum() )  #  the current algo sometimes produces NaN values
        if not Ok:
            NaN_found += 1
            uncertainty *= 1.4
    npkd.set_buffer(x[:,0])
    npkd.noise = eta
    if full_output:
        npkd.full_output = c
    if NaN_found >0:
        print ("%d NaN conditions encountered during PALMA processing"%NaN_found)
    return npkd

def test(npkd):
    print('Not implemented')
NPKData_plugin("palma", palma)
NPKData_plugin("do_palma", do_palma)
NPKData_plugin("prepare_palma", prepare_palma)
NPKData_plugin("calibdosy", dcalibdosy)

