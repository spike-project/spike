#!/usr/bin/env python 
# encoding: utf-8

"""
Authors: Marc-André
First version: 2011/11/18.
adapted by Lionel on 2013-3-6.
new version by M-A 04/2025

"""

import numpy as np
#import scipy.fftpack as fft
import numpy.fft as fft # numpy scipy seem equivalent for ifft and fft.. 


class transformations():
    """
    this class contains methods which are tools to generate transform and ttransform for 1D signal
    ttrans form data to image
    trans from image to data.
    handles real or complex data
    """
    def __init__(self, size_image, size_mesure, sampling = None, debug=0):
        """
        size_image and size_mesure are sizes of x and s space
        sampling is a nuslist [...i,..] of the index of the experimental points
        all other fields are meant to be overloaded after creation
        
        direct transform refers to S => X  // image => Data transform
        """
        self.size_image = size_image    
        self.size_mesure = size_mesure
        #print "data size is ",self.size_mesure
        #print "image size is ",self.size_image 
        self.pre_ft = self.Id       # applied before FT (in the direct transform)
        self.tpre_ft = self.Id     # its transpose (so applied after FT in ttransform)
        self.post_ft = self.Id      # applied after FT (in the direct transform)
        self.tpost_ft = self.Id     # its transpose (so applied before FT in ttransform)
        self.ft = fft.ifft          # trans: from image to data
        self.tft = fft.fft          # ttrans: from data to image.
        self.sampling = sampling       # sampling vector if not none
        self.debug = debug
    def report(self):
        "dumps content"
        for i in dir(self):
            if not i.startswith('_') :
                if not callable(i):
                    print(i, getattr(self,i))
    def Id(self, x):
        return x
    def Error(self,x):
        raise RuntimeError('The step is not initialized yet')
    def check(self):
        if True: #self.debug:
            print(f"""
size_image: {self.size_image} - size_mesure: {self.size_mesure}
sampling {len(self.sampling)} points
""")
        if self.sampling is not None:
            assert(len(self.sampling) == self.size_mesure)
            assert(max(self.sampling) <= self.size_image)
        assert( (self.pre_ft == self.Id and self.tpre_ft == self.Id) or \
                (self.pre_ft != self.Id and self.tpre_ft != self.Id) )                  # if one, then both !
        assert( (self.post_ft == self.Id and self.tpost_ft == self.Id) or \
                (self.post_ft != self.Id and self.tpost_ft != self.Id) )                # if one, then both !
    
    def zerofilling(self,x):
        # eventually zerofill
#        xx = np.zeros(self.size_image, dtype = 'complex')
        xx = np.zeros(self.size_image, dtype=x.dtype)
        xx[:len(x)] = x[:]
        x = xx
        return x
    
    def sample(self, x):
        """
        apply a sampling function - using self.sampling
        """
        #print self.sampling
        return x[self.sampling]
    def tsample(self, x):
        """
        transpose of the sampling function
        """
#        xx = np.zeros(self.size_image,'complex128')
        xx = np.zeros(self.size_image, dtype=x.dtype)
        #print xx.dtype, x.dtype, self.sampling.dtype
        xx[self.sampling] = x
        return xx
    def transform(self, s):
        """
        transform to data.
        Passing from s (image) to x (data)
        pre_ft() : s->s
        ft()     : s->x - fft.ifft by default - should not change size
        post_ft() : x->x - typically : broadening, sampling, truncating, etc...
        """
        if self.debug: print('entering trans', s.shape, s.dtype)
        if self.pre_ft != self.Id:
            s = self.pre_ft(s)
            if self.debug: print('trans pre_ft', s.shape, s.dtype)
        x = self.ft(s)
        if self.post_ft != self.Id:
            x = self.post_ft(x)
            if self.debug: print('trans post_ft', x.shape, x.dtype)
        if self.sampling is not None:
            x = self.sample(x)
            if self.debug: print('trans sample', x.shape, x.dtype)
        if self.size_mesure != len(x):      # eventually truncate
            x = x[0:self.size_mesure]
            if self.debug: print('trans trunc', x.shape, x.dtype)
        if self.debug: print('exiting trans', x.shape, x.dtype)
        return x
        
    def ttransform(self, x):
        """
        the transpose of transform
        Passing from x to s (data to image)
        """
        if self.debug: print('entering ttrans', x.shape, x.dtype)
        if self.sampling is not None:
            if self.debug: print('ttrans sample')
            x = self.tsample(x)
        elif self.size_image != len(x):      # eventually zerofill
            if self.debug: print('ttrans zerofill',len(x),self.size_image)
            x = self.zerofilling(x)
 
        if self.tpost_ft != self.Id:
            if self.debug: print('ttrans tpost_ft')
            x = self.tpost_ft(x)
        s = self.tft(x)
        if self.tpre_ft != self.Id:
            if self.debug: print('ttrans tpre_ft')
            s = self.tpre_ft(s)
        if self.debug: print('exiting ttrans', s.shape, s.dtype)
        return s

# algebra for hypercomplex 2D version
# internally, hypcpx are handled as float 2D arrays, with imag parts interleaved, on both axes, as in Spike#

from spike.NPKData import conj_ip, _conj_ip_from_float, _conj_ip_to_float, as_cpx, as_float
# def ftf2(buffer):
#     "out of place FFT-F2 HyperCpx - buffer and value are real"
#     # following works because axis=1 is contiguous
#     return _conj_ip_to_float( np.fft.fft( _conj_ip_from_float(buffer), axis=1 ) )
# def ftf1(buffer):
#     "out of place FFT-F1 HyperCpx - buffer and value are real"
#     # here axis=0 is not contiguous, we have to build the cpx, and then get each parts
#     si1,si2 = buffer.shape
#     outbuf = np.zeros((si1,si2))
#     for i in range(si2):
#         v = buffer[::2,i]  -  1j*buffer[1::2,i]
#         V =  np.fft.fft( v )
#         outbuf[::2,i] = V.real
#         outbuf[1::2,i] = V.imag
#     return outbuf
def ftf2(buffer):
    "out of place FFT-F2 HyperCpx - buffer should be complex"
    # following works because axis=1 is contiguous
    return np.fft.fft(np.conjugate(buffer))
def ftf1(buffer):
    "out of place FFT-F1 HyperCpx - buffer should be complex"
    # here axis=0 is not contiguous, we have to build the cpx, and then get each parts
    si1,si2 = buffer.shape
    outbuf = np.zeros((si1, 2*si2))
    for i in range(si2):
        V =  np.fft.fft( buffer[::2,i].real  -  1j*buffer[1::2,i].real )
        outbuf[::2, 2*i] = V.real                  # RR
        outbuf[1::2, 2*i] = V.imag                 # RI
        V =  np.fft.fft( buffer[::2,i].imag  -  1j*buffer[1::2,i].imag )
        outbuf[::2, 2*i+1] = V.real                  # IR
        outbuf[1::2, 2*i+1] = V.imag                 # II
    return as_cpx(outbuf)
def iftf2(buffer):
    "out of place iFFT-F2 HyperCpx - buffer and value are real"
    # following works because axis=1 is contiguous
    return _conj_ip_to_float( np.fft.ifft( _conj_ip_from_float(buffer), axis=1 ) )
def iftf1(buffer):
    "out of place iFFT-F1 HyperCpx - buffer and value are real"
    # here axis=0 is not contiguous, we have to build the cpx, and then get each parts
    si1,si2 = buffer.shape
    outbuf = np.zeros((si1,si2))
    for i in range(si2):
        v = buffer[::2,i]  -  1j*buffer[1::2,i]
        V =  np.fft.ifft( v )
        outbuf[::2,i] = V.real
        outbuf[1::2,i] = V.imag
    return outbuf

def fthyp2(x):
    "direct 2D-hyp-FFT"
    return ftf1(ftf2(x))
def ifthyp2(x):
    "inverse 2D-hyp-FFT"
    return iftf2(iftf1(x))
def hypabs(buffer):
    """
    compute abs() for hypercomplex 2D data-set, provided as a 2*si1, 2*si2 float matrix 
    return abs() as 2*si1, 2*si2 float, where the abs value are in the 4 fields of the hypercpx values
    """
    abuffer = np.sqrt(buffer[::2,::2]**2 + \
                      buffer[1::2,::2]**2 +\
                      buffer[::2,1::2]**2 +\
                      buffer[1::2,1::2]**2)
    a = np.zeros_like(buffer)
    a[::2,::2] = aPP
    a[1::2,::2] = aPP
    a[::2,1::2] = aPP
    a[1::2,1::2] = aPP
    return a

############ following code is utilities for compressing and decompressing sampled 2D points according to a nuslist
def samplecomp(data, sampling, datatype='complex'):
    """
    apply a sampling function - using sampling
    here we do compress the sampled 2D to a 1D compacted vector
    all buffers real or complex, and datatype tells how they should be considered

    datatype: 'real'/'complex'/'hypercomplex'  
         "complex" corresponds to F1 real and F2 complex - the other setup is not provided , and is handled as phase modulated dataset
    data : numpy real 2D buffer to sample
        if complex/hypercomplex,
            data should be complex
            if hypercomplex 
                nuslist entry handle all F1 real and imaginary parts, Bruker style
                F1 imaginary part follows the F1 real part 
    sampling is a nuslist [...(i1,i2)..] of the index of the experimental 2D points
            each entry corresponds to 1, or 2 data points depending on datatype being hypercomplex
    returns a 1D np buffer
    """
    # create
    if datatype == 'hypercomplex':
        dataout = 1j*np.zeros(2*len(sampling))  # has to be complex !  - double to accept hyperimag
    else:
        dataout = np.zeros(len(sampling), dtype=data.dtype)   # depends on data

    # then store
    for i,(i1,i2) in enumerate(sampling):
        if datatype == 'hypercomplex':
            dataout[2*i]   = data[2*i1,   i2]
            dataout[2*i+1] = data[2*i1+1, i2]
        else:
            dataout[i] = data[i1, i2]
    return dataout
def tsamplecomp(data, outshape, sampling, datatype='complex'):
    """
    transpose of the samplecomp function

    apply a sampling function 

    here we decompress the 1D compacted vector to a sampled 2D
    all buffers real, and datatype tells how they should be considered
        if complex/hypercomplex, imaginary parts are interleaved with real parts 
        in which case, one nuslist entry handle both real and imaginary parts
    data : numpy 1D real buffer to unsample
    outshape : numpy shape of the returned buffer
    datatype: 'real'/'complex'/'hypercomplex'  in which case imaginary parts is interleaved on the second axis or on both axis
    sampling is a nuslist [...(i1,i2)..] of the index of the experimental 2D points
            each entry corresponds to 1, 2 or 4 depending on datatype
    returns a 2D np buffer 
    """
    si1,si2 = outshape

    if datatype == 'real':
        dataout = np.zeros((si1,si2))
        for j,(i1,i2) in enumerate(sampling):
            dataout[i1,  i2  ] = data[j]
    if datatype == 'complex':
        dataout = np.zeros((si1, 2*si2))
        for j,(i1,i2) in enumerate(sampling):
            dataout[i1,  2*i2  ] = data[2*j]
            dataout[i1,  2*i2+1] = data[2*j+1]
    if datatype == 'hypercomplex':
        dataout = np.zeros((2*si1, 2*si2))
        for j,(i1,i2) in enumerate(sampling):
            dataout[2*i1,  2*i2  ] = data[4*j]
            dataout[2*i1,  2*i2+1] = data[4*j+1]
            dataout[2*i1+1,2*i2  ] = data[4*j+2]
            dataout[2*i1+1,2*i2+1] = data[4*j+3]
    return dataout

class transformations2D(transformations):
    """
    this class contains methods which are tools to generate transform and ttransform for 2D signal
    ttrans form data to image
    trans from image to data.

    all computationq in 2D, but everything is linearized outside (using ravel / reshape)
        agnostic Linear algo can use it

    handles real, complex or hypercomplex data
    for hypercomplex, the F2 rows are python complex and 0,2,...2n rows contain the F1 real and 1,2,..2n+1 contains the F1 imag

    """
    def __init__(self, size_image, size_mesure, datatype='complex', sampling=None, compressed=False, debug=0):
        """
        size_image and size_mesure are numpy shapes of x and s 2D space
        datatype is either 'real', 'complex', 'hypercomplex' 
        sampling is a nuslist [...(i1,i2)..] of the index of the experimental 2D points

        compressed defines how the 1D representation of the raw data is handled if there is sampling
            False: with 0 on the unmeasured points positions
            True: compressed to a vector of the measured values
        
        all other fields are meant to be overloaded after creation
        
        direct transform refers to S => X  // image => Data transform
        """
        if len(size_image) != 2: raise TypeError( "size_image should be the shape of a 2D")
        if len(size_mesure) != 2: raise TypeError( "size_mesure should be the shape of a 2D")
        self.datatype = datatype
        self.compressed = compressed

        super(transformations2D, self).__init__(size_image=size_image, size_mesure=size_mesure, sampling=sampling, debug=debug)
        if datatype == 'real':
            self.ft = fft.irfft2          # trans: from image to data
            self.tft = fft.rfft2          # ttrans: from data to image.
        elif datatype == 'complex':
            self.ft = fft.ifft2          # trans: from image to data
            self.tft = fft.fft2          # ttrans: from data to image.
        elif datatype == 'hypercomplex':
            self.ft = ifthyp2          # trans: from image to data
            self.tft = fthyp2          # ttrans: from data to image.
        else:
            raise ValueError("datatype is either 'real', 'complex', or 'hypercomplex'")

        if sampling is not None:    # build mask
            if datatype == 'hypercomplex':
                self.mask = np.zeros( (self.size_mesure[0], self.size_mesure[1]) ) # mask is real valued with room for all hyperimaginaries
                for i,j in sampling:
                    self.mask[2*i,j] = 1.0
                    self.mask[2*i+1,j] = 1.0
            else:  # real of complex
                self.mask = np.zeros(size_mesure)  # mask is real
                for i,j in sampling:
                    self.mask[i,j] = 1.0
            #compute limits
            self.min1 = min(self.sampling[:,0])
            self.min2 = min(self.sampling[:,1])
            self.max1 = max(self.sampling[:,0])
            self.max2 = max(self.sampling[:,1])

        else:
            self.mask = None
    
    def report(self):
        print(f"size_image: {self.size_image} - size_mesure: {self.size_mesure}")
        if self.sampling is not None:
            n1,n2,m1,m2 = (self.min1, self.min2, self.max1, self.max2)
            print(f"sampling {len(self.sampling)} pairs - datatype : {self.datatype} - missing values : {'compressed' if self.compressed else 'zero padded'}")
            
            NUS_ratio = len(self.sampling)/(m1+1)/(m2+1)
            print(f"""
Sampling : {n1} ... {m1} x {n2} ... {m2}  in {len(self.sampling)} pairs
NUS ratio, {100*NUS_ratio}""")

    def check(self):
        self.report()
        if self.sampling is not None:
            assert(self.max1 <= self.size_mesure[0])
            assert(max(self.sampling[0]) <= self.size_image[0])
            assert(self.max2 <= self.size_mesure[1])
            assert(max(self.sampling[1]) <= self.size_image[1])
        assert( (self.pre_ft == self.Id and self.tpre_ft == self.Id) or \
                (self.pre_ft != self.Id and self.tpre_ft != self.Id) )                  # if one, then both !
        assert( (self.post_ft == self.Id and self.tpost_ft == self.Id) or \
                (self.post_ft != self.Id and self.tpost_ft != self.Id) )                # if one, then both !

    def zerofilling(self,x):
        # eventually zerofill
        xx = np.zeros(self.size_image, dtype=x.dtype)
        si1,si2 = x.shape
        xx[:si1,:si2] = x
        return xx
    
    def truncate(self,x):
        # transpose of zerofill
        xx = np.zeros(self.size_mesure, dtype=x.dtype)
        si1,si2 = xx.shape
        xx = x[:si1,:si2]
        return xx

    ############ following code is for not compressed sampling, based on masking
    def sample(self, x):
        """
        apply a sampling function - using self.sampling
        in 2D we do not compress
        """
        #print self.sampling
        return x*self.mask
    def tsample(self, x):
        """
        transpose of the sample function
        equiv to sample in this set-up
        """
        return x*self.mask
    def transform(self, sflat):
        """
        transform to data.
        Passing from s (image) to x (data)
        pre_ft() : s->s   
        ft()     : s->x - fft.ifft by default - should not change size
        post_ft() : x->x - typically : broadening, sampling, truncating, etc...
        """
        s = sflat.reshape(self.size_image)
        if self.debug: print('entering trans', s.shape, s.dtype)
        if self.pre_ft != self.Id:
            s = self.pre_ft(s)
            if self.debug: print('  trans pre_ft', s.shape, s.dtype)
        x = self.ft(s)
        if self.debug: print('  trans ft', x.shape, x.dtype)
        if self.post_ft != self.Id:
            x = self.post_ft(x)
            if self.debug: print('  trans post_ft', x.shape, x.dtype)
        if self.size_mesure != self.size_image:      # eventually truncate
            x = self.truncate(x)
            if self.debug: print('  trans trunc', x.shape, x.dtype)
        if self.sampling is not None:
            x = self.sample(x)
            if self.debug: print('  trans sample', x.shape, x.dtype)
        if self.debug: print('exiting trans', x.shape, x.dtype)
        return x.ravel()

    def ttransform(self, xflat):
        """
        the transpose of transform
        Passing from x to s (data to image)
        """
        x = xflat.reshape(self.size_mesure)
        if self.debug: print('entering ttrans', x.shape, x.dtype)
        if self.sampling is not None:
            x = self.tsample(x)
            if self.debug: print('  ttrans sample', x.shape, x.dtype)
        if self.size_image != self.size_mesure:      # eventually zerofill
            x = self.zerofilling(x)
            if self.debug: print('  ttrans zerofill',x.shape, x.dtype)
        if self.tpost_ft != self.Id:
            x = self.tpost_ft(x)
            if self.debug: print('  ttrans tpost_ft', x.shape, x.dtype)
        s = self.tft(x)
        if self.debug: print('  ttrans ft', s.shape, s.dtype)
        if self.tpre_ft != self.Id:
            s = self.tpre_ft(s)
            if self.debug: print('  ttrans tpre_ft', s.shape, s.dtype)
        if self.debug: print('exiting ttrans', s.shape, s.dtype)
        return s.ravel()

def sampling_load(addr_sampling_file):
    '''
    Loads a sampling protocole from a list of indices stored in a file named addr_sampling_file
    returns an nparray with the sampling scheme.
    
    i.e. if b is a full dataset, data[sampling] is the sampled one
    '''
    with open(addr_sampling_file, 'r') as F:
#        print "### reads the sampling file ", addr_sampling_file
        param = read_param(F)
        F.seek(0)
        sampling = read_data(F)
#        print "sampling[0], sampling[len(sampling)//2], sampling[-1]",sampling[0], sampling[len(sampling)//2], sampling[-1]
    return sampling, param

def read_data(F):
    '''
    Reads data from the sampling file, used by sampling_load()
    '''
    data = []
    for l in F:
        if not l.startswith("#"):
            if l.strip() == "":
                continue
            data.append(int(l))
    return np.array(data)

def read_param(F):
    '''
    Reads the sampling parameters. used by sampling_load()
    '''
    """
        given F, an opend file , retrieve all parameters found in file header
        read_param returns  values in a plain dictionnary
    """
    dic = {}
    for l in F:
#        print l.rstrip()
        if not l.startswith("#"): 
            break
        v = l.rstrip().split(':')       # remove trailing chars and split around :
        if len(v)<2:    # comment lines
            pass #print l
        else:
            entry = v[0][1:].strip()
            dic[entry] = v[1].lstrip()    # metadata lines
    return dic

def sizend(shapenD):
    " computes the total size of a buffer from its shape" 
    m = 0
    for i in shapenD:
        m *= i
    return m

if __name__ == "__main__":
    tr = transformations(2000, 1000)
    tr.report()
