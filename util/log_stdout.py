#!/usr/bin/env python 
# encoding: utf-8

from __future__ import print_function
import sys

class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open("log_stdout.dat", "w")
    
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
     