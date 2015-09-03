#!/usr/bin/env python 
# encoding: utf-8

from __future__ import print_function
import numpy as np
import unittest
import json

'''
Utilities for reading and writing
'''

def write_json(data, file_name):
    '''
    Save numpy data
    '''
    with open(file_name, mode='w') as f:  
        json.dump(data.tolist(), f)
        
def read_json(file_name):
    with open(file_name) as f:
        d = np.array(json.load(f))
    return d
    
def read(f):
    """
    reads Orbitrap .dat FID files
    returns a numpy array with the FID
    """
    with open(f,'rb') as F:
        for l in F:
            if l.startswith("Data Points"): # find the number of data points
                print(re.findall(r'\d+', l)[0]) # extract this number
            if l.startswith("Data:"): # find where begin the data storage.
                break
        buf = np.fromfile(F, dtype = 'f4') # read the data from the preceding stop point.
    return np.array(buf)

class readwrite_test(unittest.TestCase):
    def test_json_numpy1D(self):
        print("######## numpy 1D ")
        a = np.arange(53)/2.23
        print(a.dtype)
        write_json(a, 'test_save_a.js')
        b = read_json('test_save_a.js')
        print("(b-a)[0] ",(b-a)[0])
        print(b.dtype)
        
    def test_json_numpy2D(self):
        print("######## numpy 2D ")
        a = np.random.randn(10,5)
        print("a.dtype ",a.dtype)
        print("a.shape ", a.shape)
        write_json(a, 'test_save_a.js')
        b = read_json('test_save_a.js')
        print("(b-a)[0] ",(b-a)[0])
        print("b.dtype ",b.dtype)
        print("b.shape ", b.shape)


if __name__ == '__main__':
    unittest.main()
