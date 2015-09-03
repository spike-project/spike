#!/usr/bin/env python 
# encoding: utf-8
from __future__ import print_function
#Toggle line numbers
def peaks2d(data, threshold=0.1):
    '''
    Extract peaks from 1d from FID
    '''
    listpk=np.where(((data > threshold*np.ones(data.shape))&# thresholding
                        (data > np.roll(data,  1, 0)) &     # laplacian - kind of
                        (data > np.roll(data, -1, 0)) &
                        (data > np.roll(data,  1, 1)) &
                        (data > np.roll(data, -1, 1))))
    listpk=[listpk[0],listpk[1]]
    return listpk[0],listpk[1]# list f1, f2