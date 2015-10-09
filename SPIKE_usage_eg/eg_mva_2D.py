#!/usr/bin/env python 
# encoding: utf-8

"""
quelques examples de traitements, en utilisant la base du nouveau NPK v2

"""

from __future__ import print_function
import NPKData as npkd
from FTICR import FTICRData
import numpy as np                  # pour faire des calculs supplémentaires  - pas utilisé ici
import matplotlib.pyplot as plt     # pour afficher les données brutes - pas utilisé ici
import Apex as A
######## chargement


# d = FTICRData(name="D:/donnees Marie/11_05_26/insuline_2D_000002.gf2")

#ou

d = A.Import_2D("../../MS-FTICR/subP_2D_19juillet.d")

######### calibration

d.axis1.specwidth = 1000*500.0      # spectral width in Hz
d.axis2.specwidth = 1000*1000.0
d.ref_mass = 344.0974          # reference mass
d.ref_freq =  1000*419.62      # freq for the reference mass
d.axis1.highmass = 1500.0            # highest mass of interest in this axis, zoom will be automatic!
d.axis2.highmass = 1500.0

########### maintenant le calcul complet

# en F2 d'abord

d.apod_sin(axis=2,maxi=0.5).chsize(d.size1,2*d.size2).rfft(axis=2)
d.flip()
d.phase(0.0, 90*d.axis1.mztoi(d.axis1.highmass), axis = 1) 
d.flop()

# en F1 ensuite

d.apod_sin(axis=1,maxi=0.5).chsize(2*d.size1,d.size2).rfft(axis=1)

d.modulus()

######### recalibration

d.axis1.left_point = d.axis1.mztoi(d.axis1.highmass)

d.axis1.specwidth = d.axis1.specwidth + d.axis1.specwidth*d.axis1.mztoi(d.axis1.highmass)/d.size1

d.currentunit = "m/z"

######### affichage

d.display(scale=3, show=True, xlabel="m/z", ylabel="m/z")
