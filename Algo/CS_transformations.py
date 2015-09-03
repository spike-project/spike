#!/usr/bin/env python 
# encoding: utf-8

"""
Authors: Marc-André
Last modified: 2011/11/18.

adapted by Lionel on 2013-3-6.
Copyright (c) 2011 IGBMC. All rights reserved.
"""

from __future__ import print_function
import numpy as np
#import scipy.fftpack as fft
import numpy.fft as fft # numpy scipy seem equivalent for ifft and fft.. 


class transformations(object):
    """
    this class contains methods which are tools to generate transform and ttransform.
    ttrans form data to image
    trans from image to data.
    """
    def __init__(self, size_image, size_mesure, sampling = None, debug=0):
        """
        size_image and size_mesure are sizes of x and s space
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
                print(i, getattr(self,i))
    def Id(self, x):
        return x
    def check(self):
        if self.debug:
            print("""
size_image: %d - size_mesure: %d
sampling %s
"""%(self.size_image, self.size_mesure, str(self.sampling)))
        if self.sampling is not None:
            assert(len(self.sampling) == self.size_mesure)
            assert(max(self.sampling) <= self.size_image)
        # assert( (self.pre_ft == self.Id and self.tpre_ft == self.Id) or \
        #         (self.pre_ft != self.Id and self.tpre_ft != self.Id) )                  # if one, then both !
        # assert( (self.post_ft == self.Id and self.tpost_ft == self.Id) or \
        #         (self.post_ft != self.Id and self.tpost_ft != self.Id) )                # if one, then both !
    
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
    
def sampling_load(addr_sampling_file):
    '''
    Loads a sampling protocole from a list of indices stored in a file named addr_sampling_file
    returns an nparray with the sampling scheme.
    
    i.e. if b is a full dataset, b[sampling] is the sampled one
    '''
    with open(addr_sampling_file, 'r') as F:
#        print "### reads the sampling file ", addr_sampling_file
        param = read_param(F)
        F.seek(0)
        sampling = read_data(F)
#        print "sampling[0], sampling[len(sampling)/2], sampling[-1]",sampling[0], sampling[len(sampling)/2], sampling[-1]
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

if __name__ == "__main__":
    tr = transformations(2000, 1000)
    tr.report()
