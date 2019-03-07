#!/usr/bin/env python 
# encoding: utf-8

'''
Created by Lionel Chiron  18/10/2013 
Copyright (c) 2013 __NMRTEC__. All rights reserved.

Utility for reporting error from standard error output via gmail. 
When sys.stderr performs a "write", a mail is sent with a report. 
Typical syntax is: 
    import sys
    sys.stderr = Logger()
    f = open("fff.jpg")
If the picture "fff.jpg" doesn't exist, an arror message is sent by mail. 
'''
from __future__ import print_function
import sys, os
from time import sleep, localtime, strftime
import os.path as op


class Logger(object):
    '''
    Standard error logger.
    '''
    def __init__(self, erase = False):
        try:
            applic = sys.modules['__main__'].__file__ # retrieve the name of the module.
        except Exception:
            print("no __file__")
            applic = ''
        self.stderr = sys.stderr
        date = self.datetime() # takes the date 
        app = op.splitext(op.basename(applic))[0]+'_' # shorten the application name
        #print "shorten name is ", app
        self.log_name = "log_" + app + date + ".dat"
        if os.path.exists(self.log_name) and  erase :
            os.remove(self.log_name)
        if os.path.exists(self.log_name):
            kindopen = 'a'
        else:
            kindopen = 'w'
        self.log = open(self.log_name, kindopen) # open log file
    
    def datetime(self):    
        return strftime('%Y-%h-%d-%Hh%M', localtime())
    
    def write(self, message):
        '''
        If error, sys.stderr will call this method
        '''
        #print "writes message"
        self.stderr.write(message)
        self.log.write(message)


if __name__ == '__main__':
    sys.stderr = Logger(erase = True)
    sys.stdout = Logger()
    print("hello toto")
    f = open("fff.jpg")

    
