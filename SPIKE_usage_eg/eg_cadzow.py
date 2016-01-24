#!/usr/bin/env python 
# encoding: utf-8


"""
quelques examples de traitements, en utilisant la base du nouveau NPK v2

"""

from __future__ import print_function
import NPKData as npkd
import Cadzow
import numpy as np                  # pour faire des calculs supplémentaie  - pas utilisé ici
import matplotlib.pyplot as plt     # pour afficher les données brutes - pas utilisé ici


d1 = npkd.NPKData(name="/DATA/Marie van/data_raw.gf2")  # je charge les données

print(d1.report())       # état des lieux

d1.apod_sin(axis=2,maxi=0.5).rfft(axis=2)               # en F2 d'abord

n_of_line=10
n_of_iter=5
order=500
Cadzow.cadzow2dmp(d1,n_of_line,n_of_iter, order)        # cadzow2dmp fait tout le travail en multiprocesseur
# pas testé ! trop long !!

d1.apod_sin(axis=1,maxi=0.5).rfft(axis=1)               # puis FT en F1

d1.modulus()                    # divise la taille par 4, je t'expliquerai
d1.display()
# en 2D, display n'a pas encore beaucoup d'option.
d1.display(scale=3,show=True)   # 3 fois plus bas. show=True affiche tout et bloque.
d1.save("data_cadzow_2D.gs2")   
# pourra être relue dans NMRnotebook


