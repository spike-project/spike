#!/usr/bin/env python
# encoding: utf-8
"""
A utility that wraps ConfigParser for NPK

Typical use is :

cp = NPKConfigParser()
cp.readfp(open(configfilename))     # where configfilename is the name of the config file
var1 =    cp.get( section, varname1)
var2 =    cp.get( section, varname2, default_value)
...
you can use
get()  getint()  getfloat()  getboolean()
see details in methods.

Created by Marc-André on 2011-10-14.
Copyright (c) 2011 IGBMC. All rights reserved.

modif on 21 - may 2012 - added getword and removing trailing comments
"""

import sys
import os
import unittest
from ConfigParser import SafeConfigParser
import re

    
class NPKConfigParser(SafeConfigParser):
    """
    this is a subclass of ConfigParser.SafeConfigParser, providing default value for values
    will never raise an error on missing values
    """
    def get(self, section, option, default = None, raw = 0, vars = None):
        """read a value from the configuration, with a default value"""
        if self.has_option(section, option):
            vv = SafeConfigParser.get(self, section, option, raw = raw, vars = vars) # returns the string after "option = "
            vl = re.split('\s*#',vv)        # this removes trailing comments
            return vl[0]
        else:
            print "Using default value for {} : {}".format(option,default) # send message if option not in configfile
            return default
    def getword(self, section, option, default = None, raw = 0, vars = None):
        "read a value from the configuration, with a default value - takes the first word of the string"
        vv = self.get(section, option, default = str(default), raw = raw, vars = vars)
        v = vv.split()[0]
        return v
    def getint(self, section, option, default = 0, raw = 0, vars = None):
        """
        read a int value from the configuration, with a default value
        understands 123 16k 32M 4G 2T ...
        """
        v = self.getword(section, option, default = str(default), raw = raw, vars = vars).lower()
        if v.endswith('k'):
            val = 1024*int(v[:-1])
        elif v.endswith('m'):
            val = 1024*1024*int(v[:-1])
        elif v.endswith('g'):
            val = 1024*1024*1024*int(v[:-1])
        elif v.endswith('t'):
            val = 1024*1024*1024*1024*int(v[:-1])
        else:
            val = int(v)
        return val
    def getfloat(self, section, option, default=0.0, raw=0, vars=None):
        """read a float value from the configuration, with a default value"""
        return float(self.getword(section, option, default=str(default), raw=raw, vars=vars))
    def getboolean(self, section, option, default="OFF", raw=0, vars=None):
        """read a boolean value from the configuration, with a default value"""
        v = self.getword(section, option, default=str(default), raw=raw, vars=vars)
        if v.lower() not in self._boolean_states:
            raise ValueError, 'Not a boolean: %s' % v
#        print self._boolean_states
        return self._boolean_states[v.lower()]

class Tests(unittest.TestCase):
    def setUp(self):
        self.verbose = 1    # verbose >0 switches on messages
    def announce(self):
        if self.verbose >0:
            print self.shortDescription()
    def test_def(self):
        "testing configparser default values"
        self.announce()
        cp = NPKConfigParser()
        cp.read("None")
        self.assertEqual("def_apex", cp.get("section", "apex", "def_apex"))              # string
        self.assertEqual(123,        cp.getint("section", "int_value", 123))   # int
        self.assertEqual(2048,       cp.getint("section", "int_value", "2k"))   # int
        self.assertEqual(3*1024*1024,  cp.getint("section", "int_value", "3M"))   # int
        self.assertEqual(123.45,     cp.getfloat("section", "f_value", 123.45))   # float
        self.assertEqual(0.0,        cp.getfloat("section", "f_value"))   # float
        self.assertFalse(   cp.getboolean("section", "b_value"))   # boolean
        self.assertTrue(   cp.getboolean("section", "b_value", "true"))   # boolean
    def test_read(self):
        "testing configparser - getting values from file"
        self.announce()
        cp = NPKConfigParser()
        cp.read("test.mscf")
        self.assertEqual("../DATA_test/ubiquitine_2D_000002.d", cp.get("import", "apex", "def_apex"))              # string
        self.assertEqual(2500.0, cp.getfloat("import", "highmass", 123.45))   # float
        self.assertEqual(True, cp.getboolean("processing", "do_modulus", False))   # bool
        self.assertEqual(16*1024*1024*1024, cp.getint("processing", "largest_file", 123))   # int

if __name__ == '__main__':
    unittest.main()