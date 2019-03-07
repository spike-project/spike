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
import threading

from time import sleep, localtime, strftime
from uuid import getnode as get_mac
import os.path as op


class Logger(object):
    '''
    Standard error logger.
    '''
    def __init__(self):
        try:
            applic = sys.modules['__main__'].__file__ # retrieve the name of the module.
        except Exception:
            print("no __file__")
            applic = ''
        self.terminal = sys.stderr
        date = self.datetime() # takes the date 
 
        app = op.splitext(op.basename(applic))[0]+'_' # shorten the application name
        #print "shorten name is ", app
        self.log_name = "log_" + app + date + ".dat"
        self.log = open(self.log_name, "w") # open log file
    
    
    def mac(self):
        mac_addr = str(get_mac())
        #print "mac address is ", mac_addr
        mac_dic = {'149885691548389': 'kartikeya'} #149885691548389
        if mac_addr in mac_dic :
            return mac_dic[mac_addr]
        else:
            return mac_addr
    
    def datetime(self):    
        return strftime('%Y-%m-%d-%H-%M', localtime()) 
    
    def write(self, message):
        '''
        If error, sys.stderr will call this method
        '''
        #print "writes message"
        self.terminal.write(message)
        self.log.write(message)


if __name__ == '__main__':
    sys.stderr = Logger()
    f = open("fff.jpg")

    
