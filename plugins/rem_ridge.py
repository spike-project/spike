#!/usr/bin/env python 
# encoding: utf-8

"""removes ridges in 2D

Created by Marc-André on 2011-08-15.
Copyright (c) 2011 IGBMC. All rights reserved.
"""

from __future__ import print_function
from spike import NPKError
from spike.NPKData import NPKData_plugin
import sys #
if sys.version_info[0] < 3:
    pass
else:
    xrange = range

def rem_ridge(data):
    """
    This function removes a F1 ridge by evaluating a mean avalue over the last 10% data of each column of a 2D
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
    
NPKData_plugin("rem_ridge", rem_ridge)

"""
rem_ridge() injection 
now on (in this running version)

data.rem_ridge()

will realize a baseline ridge correction
    """

