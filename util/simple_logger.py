#!/usr/bin/env python 
# encoding: utf-8

'''
Created by Lionel Chiron  18/10/2013 
Copyright (c) 2013 __NMRTEC__. All rights reserved.

Utility for reporting Standard Error (stderr) and Standard Output (stdout)
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
import unittest
from mail import GMAIL
import threading

class Writer(object):
    '''
    Writes two fluxes. One to the standard console(flux) and the other one in a file(log)
    If mail is given, send a mail with attached report. 
    input: 
        flux: stderr or stdout
        log: address of the log file
        mail_address: mail address for mail alerts
        log_name: log file for mail alerts
        date_in: date indicated in log file
    '''
    def __init__(self, flux, log, mail_address = None,
                    log_name = None, date_in = False):
        self.log = log
        self.flux = flux
        self.mail_address = mail_address
        self.log_name = log_name
    
    def send_mail(self):
        #print "sending mail"
        a = threading.Thread(None, self.mail_if_error, None, )
        a.start()
    
    def mail_if_error(self):
        print("sending mail")
        sleep(1)
        gm = GMAIL()
        gm.send(to = self.mail_address, subject = 'logger report',
                text = "here is lof file attached", attach = self.log_name)
    
    def write(self, message):
        '''
        Writes the flux in console and in file
        '''
        self.flux.write(message)
        self.log.write(message)

    def flush(self):
        '''
        flushes all 
        '''
        # self.flux.flush()
        # self.log.flush()
        self.flux.close()
        self.log.close()
        if self.mail_address:
            self.send_mail()

class Logger(object):
    '''
    Simple logger.
    input:
        erase : by default False. if True erase the existing log file with the same name.
        log_name : name of the log file, if no name given, it calls it anonymously "log_file"
        date : if True adds the date after the name
        date_in: date indicated in log file 
        mail: mail address for mail alerts  
    output:
        log file with name log_name
    '''
    def __init__(self, erase = False, log_name = None, date = False, date_in = False,  mail = None):
        if not log_name:
            self.log_name = 'log_file'
        else:
            self.log_name = log_name
        if date:
            self.log_name += self.datetime()
        self.log_name +=".log"
        if erase :
            choice_write = 'w'
        else:
            choice_write = 'a'
        self.log = open(self.log_name, choice_write) # open log file with name self.log_name
        self.stderr = Writer(sys.stderr, self.log,
                mail_address = mail, log_name = self.log_name) # implements the double flux redirection for stderr
        self.stdout = Writer(sys.stdout, self.log,
                mail_address = mail, log_name = self.log_name) # implements the double flux redirection for stdout
        sys.stderr  = self.stderr # reassigning sys.stderr 
        sys.stdout = self.stdout # reassigning sys.stdout
        
    def datetime(self):    
        return strftime('%Y-%h-%d-%Hh%M', localtime())

class Test(unittest.TestCase):
    Logger(erase = True)#log_name = 'blabla', mail = 'gmalert67@gmail.com'
    print("hello toto")
    f = open("fff.jpg") # opening an inexisting file
        
if __name__ == '__main__':
    unittest.main()

    
