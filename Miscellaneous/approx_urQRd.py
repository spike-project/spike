'''
#########
Algorithm urQRd method for long fid denoising, name "urQRd" stands for unfolded random QR denoising
urQRd(data, k)
data are the data to be denoised
k is the rank of the analysis 

##############
Parameters :
trunc_Omega defines the truncature used to approximate HxOmeg product.. 
trunc_QQstar defines the truncature used to approximate HxOmeg product.. 
trunc_antidiag defines the truncature used to approximate the antidiagonals averaging. 
'''
import time
import numpy as np
import numpy.linalg as linalg
from scipy.linalg import norm
import functools
import multiprocessing as mproc
import unittest


global trunc_Omega
global trunc_QQstar
global trunc_antidiag
trunc_Omega = 200
trunc_QQstar = 1000
trunc_antidiag = 1000

debug=0

# flags for Hankel2dtAlg - both should be true for optimum speed
# on the devel machine (Core i7, 2.2GHz, 4 core hyperthreaded to 8) I get : (lendata = 10000, rank=100, orda=4000)
#//     flat    seconds
#T      T       19 (24 with Nbproc=4)
#T      F       28
#F      T       50
#F      F       85
parall = True    # chose between parallel processing or linear processing
flatrand = True  # chose if doing the averaging on the antidiagonals by using the flatten method. 

def urQRd(data, k, orda = None, nbproc=None, trunc_1=None, trunc_2=None, trunc_3=None):
    '''
    rQRd algorithm. Name stands for random QR denoising.
    From a data series return a denoised series denoised
    data : the series to be denoised - a (normally complex) numpy buffer
    k : the rank of the analysis
    orda : is the order of the analysis
        internally, a Hankel matrix (M,N) is constructed, with M = orda and N = len(data)-orda+1
        if None (default) orda = (len(data)+1)/2

    runs urQRCore, 
    if no orda given it takes orda = (data.size+1)/2
    the matrix Omega is done locally.
    nbproc is the number of cores to use for // processing. None means maximum possible, 1 is monoprocessor.
    trunc_1 trunc_2 trunc_3 allow to modulate truncation level (see documentation in source)
    '''

    global trunc_Omega  # these are given default values at import time
    global trunc_QQstar
    global trunc_antidiag
    if trunc_1 : trunc_Omega = trunc_1
    if trunc_2 : trunc_QQstar = trunc_2
    if trunc_3 : trunc_antidiag = trunc_3
    
    if np.allclose(data,0.0):
        return data
    if not orda:
        orda = (data.size-1)/2
    if (2*orda > data.size):
        raise(Exception('orda is too large'))
    if (k >= orda):
        raise(Exception('rank is too large'))
        
    N = len(data)-orda+1
    
    Omega = np.random.normal(size = (N,k))
        
    rdata = urQRdCore(data, orda = orda, Omega=Omega, nbproc=nbproc) # makes  H=QQ*H algorithmicly
    if debug>0 : print "fid.shape after rebuild fid ", data.shape
    if rdata.dtype == "float":
        rdata = np.real(rdata)
    return rdata

#---------------------------------------------
def urQRdCore(fid, orda, Omega, nbproc=None):
    '''
    Makes the Y product, then QR decomposition 
    and the averaging on the antidiagonals 
    Calls Yprod, QstarHprod and Hank2dt
    It returns the data denoised. 
    '''
    k = Omega.shape[1]
    if debug>0: print "=== in urQRdCore ==="
    lengthfid = len(fid)
    Y = Yprod(fid, orda, Omega, nbproc=nbproc)###### Y= H Omeg
    Q,r = linalg.qr(Y)###### QR decomposition
    del(r)
    QstarH = QstarHprod(fid, orda, Q.conj().T, nbproc=nbproc)######  Q*H
    data = Hank2dt(fid = fid, orda = orda, Q = Q, QstarH = QstarH, nbproc=nbproc)  ####### QQ*H + Hankel to fid
    return data

def fY_r(i, lengthfid, orda):
    """
    Calculate each column Ycol for the product Y = H Omeg
    For constructing Y and then doing QR decomposition
    It is called by Yprod(fid,orda,Omeg,k) 
    trunc_Omega define the truncature used to approximate HxOmeg product.. 
    """
    global Omeg, fid
    global trunc_Omega
#    lenfiminorda = lengthfid-orda+1
    N = lengthfid-orda+1
    M = orda
    Ycol = np.zeros(M, dtype = complex)
    l = min(trunc_Omega, N)
    if l<N:     # l smaller than M - do sampling
        perm = np.random.permutation(N)
        p = perm[:l]
        p.sort()
    else:
        p = np.arange(N)   
    for j in range(M):# spans lines of H
        ff = fid[j:j+N]
        OmegOmeg = Omeg[:N,i]# ith column
        Ycol[j] = np.dot(ff[p],OmegOmeg[p])
    Ycol *= float(N)/float(l)
    return Ycol
    


def Yprod(fid, orda, Omeg, nbproc=None):   # product Y = M*Omeg, first product for producing the basis
    '''
    Makes the random basis Y from the hankel matrix H and the random matrix Omega
    product Y = H x Omeg 
    Calculus parallelized on the columns
    '''
    k = Omeg.shape[1]
    if debug>0: print "== in Yprod =="
    lengthfid= len(fid)
    inject(('fid',fid))
    inject(('Omeg',Omeg))
    Y = np.empty((orda, k), dtype=complex)
    todo = range(k) 
    pool = mproc.Pool(nbproc)
    res = pool.imap(functools.partial(fY_r,lengthfid = lengthfid,orda = orda), todo)
    for i in range(k):
        col = res.next()
        Y[:,i] = col
    pool.close()
    pool.join()
    return Y

#---------------------------------------------
def fqstarh_r(i,lengthfid,orda,k):
    """
    used by QstarHprod to compute randomly the product Q*H, without storing H
    performs product random sampling to speed-up computation
    trunc_QQstar defines the truncature used to approximate product. 
    """
    global Qstar, fid
    global trunc_QQstar
    colqstarh = np.empty(k, dtype=complex)
#    lenfidminorda = len(fid)-orda+1
    N = len(fid)-orda+1
    M = orda
    l = min(trunc_QQstar, M)
    if l<M:
        perm = np.random.permutation(M)
        p = perm[:l]
        p.sort()
    else:
        p = np.arange(M)
    for j in range(k):
        Qs = Qstar[j,:]
        f = fid[i:i+M]
        colqstarh[j] = np.dot(Qs[p], f[p])
    colqstarh *= float(M)/float(l)
    return colqstarh



def QstarHprod(fid, orda, Qstar, nbproc=None):   # Q*H
    '''
    Makes the intermediate product Q*H before the complete product QQ*H
    Calculus parallelized on the columns. 
    '''
    if debug>0: print "== in QstarHprod =="
    lengthfid = len(fid)
    k = Qstar.shape[0]
    N = len(fid)-orda+1
    M = orda
    QstarH = np.empty((k,N),dtype=complex)##
    todo = range(N)
    inject(('Qstar',Qstar))
    pool = mproc.Pool(nbproc)
    res = pool.imap(functools.partial(fqstarh_r, lengthfid = lengthfid, orda = orda,k = k),todo)
    for i in range(N):
        col = res.next() 
        QstarH[:,i] = col  
       
    pool.close()
    pool.join()
    return QstarH

####################
def antidiag1flatrnd(p,M,N):
    '''
    Makes the averaging of each antidiagonal in the left part of the matrix QQ*H
    It uses random sampling plus flatten function f numpy
    called by Hank2dt(fid,orda,rank,Q,QstarH) through paralldiag(inf, sup, func, fid, M, N, pool)
    '''
    global fid, Q, QstarH
    global trunc_antidiag
    l = min(p+1,trunc_antidiag)  # truncating the mean
    z=np.random.permutation(l)
    c=Q[::-1,:][(M-p-1):,:][z,:].ravel()
    d=QstarH[:,:(p+1)][:,z].T.ravel()
    mean=np.dot(c,d)

    return mean/l

def antidiag2flatrnd(p,M,N):
    '''
    Makes the averaging of each antidiagonal in the middle part of the matrix QQ*H
    It uses random sampling plus flatten function f numpy
    called by Hank2dt(fid,orda,rank,Q,QstarH) through paralldiag(inf, sup, func, fid, M, N, pool)
    '''
    global fid, Q, QstarH
    global trunc_antidiag
    l = min(M,trunc_antidiag)
    z=np.random.permutation(l)
    c=Q[::-1,:][z,:].ravel()
    d=QstarH[:,(p-M+1):(p+1)][:,z].T.ravel()
    mean=np.dot(c,d)
    return mean/l

def antidiag3flatrnd(p,M,N): 
    '''
    Makes the averaging of each antidiagonal in the right part of the matrix QQ*H
    It uses random sampling plus flatten function f numpy
    called by Hank2dt(fid,orda,rank,Q,QstarH) through paralldiag(inf, sup, func, fid, M, N, pool)
    '''
    global fid, Q, QstarH   
    global trunc_antidiag
    l = min(N+M-p-1,trunc_antidiag)
    z=np.random.permutation(l)
    c=Q[::-1,:][:(M+N-1-p),:][z,:].ravel()
    d=QstarH[:,(p-M+1):N][:,z].T.ravel()
    mean=np.dot(c,d) 
    return mean/l

def antidiag1(p,M,N):
    '''
    used by paralldiag to compute averaging antidiagonals on the 1st part of QQ*H
    called by Hank2dt(fid,orda,rank,Q,QstarH) through paralldiag(inf, sup, func, fid, M, N, pool)
    '''
    global fid, Q, QstarH
    global trunc_antidiag
    mean=0
    l = min(p+1,trunc_antidiag)  # truncating the mean
    for j in np.random.permutation(p+1)[:l]:     
        QstarHj = QstarH[:,j] # matrix k x N
        Qi = Q[p-j,:] # matrix M x k
        mean += np.dot(Qi,QstarHj) # direct calculus of Hij
    return mean/l

def antidiag2(p,M,N):
    '''
    used by paralldiag to compute averaging antidiagonals on the 2nd part of QQ*H
    called by Hank2dt(fid,orda,rank,Q,QstarH) through paralldiag(inf, sup, func, fid, M, N, pool)
    '''
    global fid, Q, QstarH
    global trunc_antidiag
    mean=0
    l = min(M,trunc_antidiag)
    for j in np.random.permutation(range(p-M+1,p+1))[:l]:
        QstarHj=QstarH[:,j] # matrix k x N
        Qi=Q[p-j,:] # matrix M x k
        mean += np.dot(Qi,QstarHj) # direct calculus of Hij
    return mean/l

def antidiag3(p,M,N): 
    '''
    used by paralldiag to compute averaging antidiagonals on the 3rd part of QQ*H
    called by Hank2dt(fid,orda,rank,Q,QstarH) through paralldiag(inf, sup, func, fid, M, N, pool)
    '''
    global fid, Q, QstarH   
    global trunc_antidiag
    mean=0
    l = min(N+M-p+1,trunc_antidiag)
    for j in np.random.permutation(range(p-M+1,N))[:l]:
        QstarHj=QstarH[:,j] # matrix k x N
        Qi=Q[p-j,:] # matrix M x k
        mean += np.dot(Qi,QstarHj) #  direct calculus of Hij
    return mean/l
    
def inject(arg):
    """
    this functions injects in the global dictionnary a new object
    arg is (name, value)
    
    used to set-up global before Pool
    generates many errors for pylint which does not understand the techniques
    """
    (nom, val) = arg
    globals()[nom]=val
#    print "inject ",nom, val.dtype
    return None

def paralldiag(inf, sup, func, fidres, M, N, pool):
    """ 
    call antidiagx utilities through a pool f workers
    requires python 2.7 as partial is broken in 2.6
    """
    todo = np.arange(inf,sup)
    res = pool.imap(functools.partial(func, M=M, N=N),todo)
    for p in range(inf,sup):
        fidres[p] = res.next()
#    return fidres

def Hank2dt(fid, orda, Q, QstarH, nbproc=None):
    """
    computes fid from Hankel matrix, taking the mean of antidiagonals
    """
 
    
    if debug>0: print "== in Hank2dt =="
    M = orda #number of lines of Hankel martix
    N = len(fid)-orda+1 #number of columns of Hankel martix
    fidresult = np.zeros(len(fid), dtype=complex)# initialize the fid
    listnoflat = [antidiag1, antidiag2, antidiag3]
    listflat = [antidiag1flatrnd, antidiag2flatrnd, antidiag3flatrnd]
    if flatrand :
        antid=listflat# list of the method to call for doing the flatten method randomly. 
    else :
        antid=listnoflat## list of the method to call for doing multiple dot product randomly
    if parall :
        if debug>0: print "use parallel"
        inject(('Q',Q))
        inject(('QstarH',QstarH))
        pool = mproc.Pool(nbproc)
        paralldiag(0, M, antid[0], fidresult, M, N, pool)
        paralldiag(M, N, antid[1], fidresult, M, N, pool)
        paralldiag(N, M+N-1, antid[2], fidresult, M, N, pool)
        pool.close()
        pool.join()
    else :
        inject(('fid',fid))
        inject(('Q',Q))
        inject(('QstarH',QstarH))
        for p in range(M):
            fidresult[p] = antid[0](p,M,N)
        for p in range(M,N):
            fidresult[p] = antid[1](p,M,N)
        for p in range(N,M+N-1):
            fidresult[p] = antid[2](p,M,N)
    return np.array(fidresult)


class urQRd_Tests(unittest.TestCase):
    def test_rQRd(  self,
        lendata = 10000,
        rank = 100,
        orda = 4000,
        noise = 200.0,
        noisetype = "additive"):
        """
        ============== example of use of rQRd on a synthetic data-set ===============
        """
        import spike.Display.testplot as testplot
        plt = testplot.plot()

        def mfft(v):
            "utility that returns the modulus of the fft of v"
            import scipy.fftpack as fft
            s0 = fft.fft(v)    ###
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
        print "=== Running urQR algo ===",
        print "lendata:",lendata,
        print " orda:",orda,
        print ' rank:',rank

        t0 = time.time()
        datarqrd = urQRd(data, k=rank, orda=orda) # denoise signal with rQRd
        trQRd = time.time()-t0
        
        fdatarqrd  = mfft(datarqrd )# FFT of rQRd denoised signal
        # normrQR = norm(fdatarqrd -fdata)/norm(fdata)
        # print "= normratio ",normrQR
        print "=== Result ==="
        fSNR = SNR(datarqrd, data0)
        print "Denoised SNR: %.2f dB  - processing gain : %.2f dB"%( fSNR, fSNR-iSNR )
        print "processing time for rQRd : %f sec"%trQRd
        if noisetype == "additive":
            self.assertTrue( fSNR-iSNR>16.0 )   # should give 18 dB on 100000 points
        ################################################################# Plotting
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
        plt.plot(fdatarqrd ,'r', label= 'urQRd filtered spectrum') # plot the signal denoised with rQRd
        plt.legend()
#        plt.savefig("rQrd.png")
        plt.show()

     
if __name__ == '__main__':
    debug = 1
    unittest.main()