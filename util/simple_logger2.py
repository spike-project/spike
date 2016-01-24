#!/usr/bin/env python 
# encoding: utf-8

'''
Created by Lionel Chiron  18/10/2013 
Copyright (c) 2013 __NMRTEC__. All rights reserved.

Utility for copying Standard Error (stderr) and Standard Output (stdout) to a file

'''

from __future__ import print_function
import sys, os
from datetime import datetime
import os.path as op
import unittest
import threading

def strdate():
    n = datetime.now()
    #return n.strftime('%Y-%h-%d-%Hh%M') #CR error at runtime on Windows
    return n.strftime('%a-%d-%b-%Y-%H-%M-%S')

class _Writer(object):
    '''
    Writes two fluxes. One to the standard console(flux) and the other one in a file(log)
    input: 
        flux: stderr or stdout
        log: file descriptor of the log file
    '''
    def __init__(self, flux, log, prefix=""):
        self.log = log
        self.flux = flux
        self.prefix = prefix
        self.newline = True

    def write(self, message):
        '''
        Writes the flux in console and in file
        '''
        self.flux.write(message)
        if self.newline:
            self.log.write(self.prefix+message)
        else:
            self.log.write(message)
        self.newline = message.endswith("\n")

    def close(self):
        '''
        do nothing
        '''
        # self.flux.close()
        # self.log.close()
        pass
    def flush(self):
        '''
        flush both flux
        '''
        self.flux.flush()
        self.log.flush()

class TeeLogger(object):
    '''
    Simple logger.

    TeeLogger(log_name = "log_file.log", erase = False, date_in = True, err_prefix = "___")
    or 
    TeeLogger()

    copies standard Error (stderr) and Standard Output (stdout) to a file
    
    log_name : name of the log file, if no name given, it is called "log_file.log"
    erase : if True erase the existing log file with the same name.
    date_in: in True (defualt) date is indicated in log file 
    err_prefix stderr is prefixed with this string
    '''
    def __init__(self, log_name = "log_file.log", erase = False, date_in = True, err_prefix = "___"):
        self.log_name = log_name
        if erase :
            choice_write = 'w'
        else:
            choice_write = 'a'
        self.log = open(self.log_name, choice_write) # open log file with name self.log_name
        if date_in:
            self.log.write("""
========================================
TeeLogging information from stdout and stderr
%s

"""%(strdate()))
        self.stdout = _Writer(sys.stdout, self.log) # implements the double flux redirection for stdout
        self.stderr = _Writer(sys.stderr, self.log, prefix=err_prefix) # implements the double flux redirection for stderr
        sys.stderr  = self.stderr # reassigning sys.stderr 
        sys.stdout = self.stdout # reassigning sys.stdout

if __name__ == '__main__':
    unittest.main()

    
