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
from File.Apex import Import_1D

import numpy as np                  # pour faire des calculs
import matplotlib.pyplot as plt     # pour afficher les données brutes
from FTICR import FTICRData

eg = 2.1

if eg == 1:
    # Show FID
    f = Import_1D("../DATA_test/angio_ms_000005.d")
    f.report()
    f.display(label = "FID")
    ff = f.copy()
    ff.chsize(len(o.buffer)/2)
    ff.buffer = ff.buffer[:len(f.buffer)/2]/2
    f.display(label = "FID")
    ff.display(label = "FID cut", new_fig = False)

if eg == 2:
    # FFT with zerofilling
    f = Import_1D("../DATA_test/angio_ms_000005.d")
    f.report()
    f.currentunit = 'm/z'
    f.apod_sin(maxi = 0.5).chsize(f.buffer.size*2).rfft().modulus().display(label = "zerofill x2", show = True)

if eg == 2.1:
    # FFT with zerofilling, processing cutting the pipes.
    f = Import_1D("../DATA_test/angio_ms_000005.d")
    f.currentunit = 'm/z'
    f.apod_sin(maxi = 0.5)
    f.chsize(f.buffer.size*2).rfft()
    f.modulus().display(label = "zerofill x2", show = True)

# 
# if eg == 2.1:
#     # FFT with zerofilling
#     f = FTICRData(name = "../DATA_test/1D_test.msh5")
#     f.display(label = "FID of 1D_test.msh5", show = True)
#     # f.report()
#     # f.apod_sin(maxi = 0.5).chsize(f.buffer.size*2).rfft().modulus().display(label = "zerofill x2", show = True)
    
if eg == 3:
    # urQRd
    f = Import_1D("../DATA_test/angio_ms_000005.d")
    f.units = 'm/z'
    f.urqrd(k = 300).rfft().modulus().display(label = "urQRd, rank = 300", show = True)

plt.show()
