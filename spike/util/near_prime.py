#!/usr/bin/env python 
# encoding: utf-8

from __future__ import print_function

def nearfactor(num, thresh):
    '''
    find the product of prime factors of num the nearest from thresh 
    '''
    def primefactors(x):
        '''
        decomposition in prime factors
        '''
        factorlist = []
        loop = 2
        while loop <= x:
            if x%loop == 0:
                x /= loop
                factorlist.append(loop)
            else:
                loop += 1
        return factorlist
    lprime =  primefactors(num)
    # from decompositon in primes, find the nearer product from thresh
    div = 1
    for i in lprime[::-1]:
        div*=i
        if num/div != 0 and div <= thresh:
            pass
        else :
            div /= i
    return div