#!/usr/bin/env python 
# encoding: utf-8

"""
testplot

allow to import either matplotlib.pyplor or fakeplot depending on the PLOT flag

usage:


import Display.testplot as testplot
testplot.PLOT = False    # eventually
plt = testplot.plot()


Created by Marc-André on 2012-10-03.
Copyright (c) 2012 IGBMC. All rights reserved.
"""
from __future__ import print_function
import sys
PLOT = True

def plot():
    """
    import the current plotpackage
    usage:

    import spike.Display.testplot as testplot
    plt = testplot.plot()
    
    then use plt as matplotlib
    """
#    print "importing ",plotname()
    __import__(plotname())
    return sys.modules[plotname()]

def plotname():
    "returns current plot package name"
    if PLOT:
        return "matplotlib.pyplot"
    else:
        return "spike.Display.fakeplot"
