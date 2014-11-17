#!/usr/bin/env python
# encoding: utf-8
"""
rem_ridge.py

Created by Marc-André on 2011-08-15.
Copyright (c) 2011 IGBMC. All rights reserved.
"""

import NPKData

def rem_ridge(data):
    """
    cette fonction soustrait un ridge en F1 en l'évaluant sur les derniers 10% d'une 2D
    prévue pour être injectée dans NPKData.NPKData au run-time
    """
    data.check2D()
    deb = int(0.9*data.size1)   # debut et fin de l'évaluation
    fin = data.size1
    r = data.row(deb)
    for i in xrange(deb+1, fin):    # je calcule la moyenne
        r.add(data.row(i))
    r.mult(-1.0/(fin-deb))
    for i in xrange(data.size1):
        data.set_row(i, data.row(i).add(r) )
    return data     # et garde la syntaxe standard NPKData
    
if __name__ == '__main__':
    # fait l'injection
    NPKData.NPKData.rem_ridge = rem_ridge
    print u"""
injection de rem_ridge()
maintenant il suffit de faire 

data.rem_ridge()

pour réaliser la correction de ldb
    """

