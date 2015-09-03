#!/usr/bin/env python 
# encoding: utf-8

"""
Tests for v1 compatibility

Created by Marc-André on 2012-10-09.
Copyright (c) 2012 IGBMC. All rights reserved.
"""

from __future__ import print_function
import unittest
from v1.Kore import *

class v1Tests(unittest.TestCase):
    def setUp(self):
        self.verbose = 1    # verbose >0 switches on messages
    def announce(self):
        if self.verbose >0:
            print(self.shortDescription())
    def test_BrukerImport(self):
        "tests Bruker 2D Import"
        from ..v1 import Bruker
        from ..Tests import filename
        name2D = filename("Lasalocid-Tocsy/dataset/ser")
        self.announce()
        rep = Bruker.Import_2D(name2D)
        self.assertTrue(rep == (328, 2048))
        d = get_Kore_2D()
        # make a few checks
        self.assertAlmostEqual(d.axis1.specwidth, 5201.560468, 6)
        self.assertAlmostEqual(d[100,100],80)
        self.assertAlmostEqual(d[4,400], -3186.)
    
if __name__ == '__main__':
    unittest.main()