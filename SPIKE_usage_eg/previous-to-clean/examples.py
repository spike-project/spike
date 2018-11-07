#!/usr/bin/env python 
# encoding: utf-8

"""
quelques examples de traitements, en utilisant la base du nouveau NPK v2

il manque encore le module Kore qui fait le lien avec NPK v1
"""

from __future__ import print_function

import File.GifaFile as gf
import NPKData as npkd
import numpy as np                  # pour faire des calculs
import matplotlib.pyplot as plt     # pour afficher les données brutes

eg = 3

if eg == 1:
    # crée le et charge les données
    R = npkd.NPKData(name = "mad2-zoom.gf2")
    print(R.report())
    R.axis1.itype = 1

    # affiche les premières lignes en module
    R.row(0).modulus().display(label="premiere ligne")
    R.row(1).modulus().display(label="deuxieme ligne")

    # additionne les 10 premiers
    r = R.row(0).modulus()
    for i in range(1,10):
        r.buffer += R.row(i).modulus().buffer
    r.display(label="somme 10 premiers")

    # additionne tous
    r = R.row(0).modulus()
    for i in range(1,R.size1):
         r.buffer += R.row(i).modulus().buffer
    r.display(label="somme de tous")
    r.save("sptout4k.gf1")

    # on peut aussi
    R.col(4249).apod_sin(maxi = 0.5).chsize(4096).fft().modulus().display(label = "4249 depart",show=True)
    # l'option show force l'affichage
    # inutile sous ipython, mais utile sinon

elif eg == 2:
    # charge un spectre de masse
    r = npkd.NPKData(name="sum_tot.gs1")
    # specw = 1667 kHz, le pic à 419.62kHz est à m/z 344.0974
    masses = 344.0974*419.62/(1667.0* (1+np.arange(r.size1))/r.size1 )
    
    z = 2000    # pour zoomer
    # affiche en masse
    r.display(axis = masses, zoom = (z,-1))
    # calcul les pics et affiche les, (bien mettre le meme zoon et le meme axe)
    r.peak().display_peaks(axis = masses, zoom = (z,-1))  
    # nomme les plus grands
    r.peak(threshold = 0.2).display_peaks(peak_label = True,axis=masses,zoom=(z,-1))  

elif eg == 3:
    # simule 
    r = npkd.NPKData(shape = (100,100))   # pas besoin de 2^n !
    r.fill(1.0).apod_sin(axis = 1,maxi=0).apod_sin(axis = 2,maxi = 0)
    # et fft
    r.axis1.itype = 0
    r.axis2.itype = 0
    r.chsize(2*r.size1,2*r.size2).revf(2).rfft(2).revf(1).rfft(1).modulus().display()
    print(r.report())

# l'option show affiche tout les display() créés
plt.show()
