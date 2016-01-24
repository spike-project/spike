#!/usr/bin/env python 
# encoding: utf-8


"""
quelques examples de traitements, en utilisant la base du nouveau NPK v2

"""
from __future__ import print_function

import NPKData as npkd
import numpy as np                  # pour faire des calculs supplémentaires  - pas utilisé ici
import matplotlib.pyplot as plt     # pour afficher les données brutes - pas utilisé ici


######## chargement
d1 = npkd.NPKData(name="/DATA/Marie van/data_raw.gf2")  # je charge les données

print(d1.report())       # état des lieux

###### calculs 1D
d1.row(0).rfft().modulus().display(label="FT de la premiere ligne")    # je montre la 1ere ligne

d1.row(0).apod_sq_sin(maxi=0.5).rfft().modulus().display(new_fig=False,label="avec apodisation")   # new_fig=False dessine sur la même image
# l'apodisation permet d'améliorer l'aspect des raies après modulus()

# on peut calculer un spectre 1D, somme de toutes les lignes :
sum = npkd.NPKData(shape=(2*d1.axis2.size,))     # on crée une données vide, *4 pour chsize /2 pour modulus
for i in range(d1.size1):
        ri = d1.row(i).apod_sin(maxi=0.5).chsize(4*d1.size2).rfft().modulus()      # on fait les calculs
        sum.add(ri)                                             # et on aditionne
sum.save("sum_1D.gs2")  # pourra être relue dans NMRnotebook
sum.display(label="somme de tout")
sum.peak()
sum.display_peaks()
# ici j'ai fait du "zerofilling" 2 fois : chsize(4*d1.size2)  pour augmenter la résolution.

# je peut même y mettre des axes :
# specw = 1667 kHz, le pic à 419.62kHz est à m/z 344.0974
masses = 344.0974*419.62/(1667.0 * (1+np.arange(sum.size1))/sum.size1 )  # axe est un tableau numpy contenant les coordonnées
sum.display(axis=masses,zoom=(8000,-1)) # axis : défini l'axe à utliser; zoom=(debut,fin) (-1 c'est le bout) enlève le bruit à fréquence faible
sum.display_peaks(axis=masses,zoom=(8000,-1))


########### maintenant le calcul complet
#d1.apod_sin(axis=2,maxi=0.5).chsize(d1.size2,2*d1.size2).rfft(axis=2)               # en F2 d'abord
#d1.apod_sin(axis=1,maxi=0.5).chsize(2*d1.size2,d1.size2).rfft(axis=1)               # puis en F1
# ces chsize me dise "dimensions too large" sur mon systeme 32bits - dommage - pourtant ca devrait passer
d1.apod_sin(axis=2,maxi=0.5).rfft(axis=2)               # en F2 d'abord
d1.apod_sin(axis=1,maxi=0.5).rfft(axis=1)               # puis en F1
d1.modulus()                    # divise la taille par 4, je t'expliquerai
d1.display()
# en 2D, display n'a pas encore beaucoup d'option.
d1.display(scale=3,show=True)   # 3 fois plus bas. show=True affiche tout et bloque.
d1.save("data_2D.gs2")  # pourra être relue dans NMRnotebook


