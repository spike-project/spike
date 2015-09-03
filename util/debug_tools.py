#!/usr/bin/env python 
# encoding: utf-8

'''
Created by Lionel Chiron  07/10/2013 
Copyright (c) 2013 __NMRTEC__. All rights reserved.
decorate the class with decclassdebugging and in each method to be debugged make :
if debug(self) : instructions.. 
'''

from __future__ import print_function
import sys, os
import inspect
import unittest

def pr(self, var, message = ''):
    '''
    print
    '''
    frame = inspect.currentframe()
    loc = frame.f_back.f_locals
    def show(value):
        if message == '':
            print(var + ' ', value) 
        else:
            print(message + ' ' + var + ' ', value) 
    try:
        show(eval(var))
    except Exception:
        show(loc[var])

def dec_class_pr(cl):
    '''
    Decorates the class for printing the local variables
    '''
    setattr(cl, 'pr', pr)
    return cl

def name_up():
    frame = sys._getframe(2)
    return frame.f_code.co_name

def debug(cl):
    '''
    returns the name of the class followed by "_debug"
    '''
    return getattr(cl, name_up() + '_debug')

def decclassdebugging(cl):
    '''
    decorator for debugging classes.
    It allows to trigger specifically which method in a class will be observed.
    if debug(self) : " where debug is needed
    After instantiation of the class, add class.method_name_debug = True to activate the debugging.
    '''
    for m in dir(cl):
        if len(m.split('__')) == 1: 
            setattr(cl, m + '_debug', False)
    return cl

@dec_class_pr
@decclassdebugging
class Test_debug_tools(unittest.TestCase):
        
    def hip(self):
        print("passing by hip")
        ici = 'here'
        if debug(self) :
            print("hip is being debugged..!")
            self.pr('ici')
        
    def hop(self):
        print("you are in hop")
        variable = 0
        if debug(self) :
            self.pr('variable', 'the unknown is ')
            
    def test_debug(self):
        '''
        On this example as just hip is activated, 
        debug(self) is True only in hip.
        '''
        self.hip_debug = True
        self.hip()
        self.hop()

if __name__ == '__main__':
    unittest.main()
    
