#!/usr/bin/env python 
# encoding: utf-8

"""
Few examples of how to use NPK v2. 
1 simple FID
    default display parameters:
        show = False, label = None, new_fig = True, axis = None, xlabel="_def_", ylabel="_def_"
2 FFT with zerofilling
3 RECTIAL
4 urQRd

"""

from __future__ import print_function
from File.Thermo import Import_1D

import numpy as np                  # pour faire des calculs
import matplotlib.pyplot as plt     # pour afficher les données brutes

eg = 1

if eg == 1:
    # Show FID
    o = Import_1D("../DATA_test/ubiquitin_5_scan_res_30000_1.dat")
    o.report()
    oo = o.copy()
    oo.chsize(len(o.buffer)/2)
    oo.buffer = oo.buffer[:len(o.buffer)/2]/2
    o.display(label = "FID")
    oo.display(label = "FID cut", new_fig = False)

if eg == 2:
    # FFT with zerofilling
    o = Import_1D("../DATA_test/ubiquitin_5_scan_res_30000_1.dat")
    o.report()
    o.currentunit = 'm/z'
    o.apod_sin(maxi = 0.5).chsize(o.buffer.size*2).rfft().modulus().display(label = "zerofill x2", show = True)

if eg == 2.1:
    # FFT with zerofilling, processing cutting the pipes.
    o = Import_1D("../DATA_test/angio_ms_000005.d")
    o.currentunit = 'm/z'
    o.apod_sin(maxi = 0.5)
    o.chsize(o.buffer.size*2).rfft()
    o.modulus().display(label = "zerofill x2", show = True)

if eg == 3:
    # urQRd
    o = Import_1D("../DATA_test/ubiquitin_5_scan_res_30000_1.dat")
    o.report()
    o.currentunit = 'm/z'
    o.urqrd(k = 300).rfft().modulus().display(label = "urQRd, rank = 300", show = True)

plt.show()
