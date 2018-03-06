#!/usr/bin/env python 
# encoding: utf-8


"""
fakeplot.py

This library allows to fake the matplotlib.plot module where it has not been installed.

Created by Marc-André on 2010-03-22.
Copyright (c) 2010 IGBMC. All rights reserved.
"""

from __future__ import print_function

__version__ = "0.2.1"


import sys
import os

current = 0
FAKE = True

class fake(object):
    """fake figure object"""
    def __init__(self, *args, **key):
        global current
        current += 1
        print("******* figure %d *****"%current)
    def figure(self, *args, **key):
        return fake(*args, **key)
    def plot(self, *args, **key):
        x = args[0]
        label = key.get("label","")
        plotname = key.get("plot_name","plot")
        print("-- %s : %d points  %s"%(plotname, len(x), label))
    def hist(self, *args, **key):
        key['plot_name'] = 'histogram'
        plot(*args, **key)
        try:
            return (args[1],args[0],[3,4])
        except:
            return (len(args[0]),args[0],[3,4])
    def errorbar(self, *args, **key):
        key['plot_name'] = 'errorbar plot'
        plot(*args, **key)
    def semilogx(self, *args, **key):
        plot(*args, **key)
        xscale('log')
    def semilogy(self, *args, **key):
        plot(*args, **key)
        yscale('log')
    def loglog(self, *args, **key):
        plot(*args, **key)
        xscale('log')
        yscale('log')
    def contour(self, *args, **key):
        if len(args) in (1,2):
            (n,m) = args[0].shape
        elif len(args) in (3,4):
            (n,m) = ( len(args[0]), len(args[1]))
            if (n,m) != args[2].shape:
                raise "size missmatch in contour"
        print("-- contour-plot : %d x %d"%(n,m))
    contourf = contour
    def yscale(self, *args, **key):
        print("-- yscale",args)
    def xscale(self, *args, **key):
        print("-- xscale",args)
    def set_xscale(self, *args, **key):
        self.xscale(*args, **key)
    def set_yscale(self, *args, **key):
        self.yscale(*args, **key)
    def text(self, *args, **key):
        print("-- text :",args)
    def legend(self, *args, **key):
        pass
    def show(self, *args, **key):
        print("******** show")
    def xlabel(self, *args, **key):
        print("-- xlabel :",args[0])
    set_xlabel = xlabel
    def ylabel(self, *args, **key):
        print("-- ylabel :",args[0])
    set_ylabel = ylabel
    def title(self, *args, **key):
        print("-- title :",args[0])
    def set_title(self, *args, **key):
        self.title(*args, **key)    
    def suptitle(self, *args, **key):
        print("-- suptitle :",args[0])
    def grid(self, *args, **key):
        if args[0]: 
            print('grid on')
    def add_subplot(self, *args, **key):
        print("*** subplot :",args[0])
        return self
    def set_xticklabels(self, *args, **key):
        print("-- xtick :",args[0])
    def set_yticklabels(self, *args, **key):
        print("-- ytick :",args[0])
    def xlim(self, *args, **key):
        print("-- xlim :",args[0])
    def ylim(self, *args, **key):
        print("-- ylim :",args[0])
def figure(*args, **key):
    return gca().figure(*args, **key)
def subplot(*args, **key):
    return gca().add_subplot(*args, **key)
def plot(*args, **key):
    return gca().plot(*args, **key)
def hist(*args, **key):
    return gca().hist(*args, **key)
def errorbar(*args, **key):
    return gca().errorbar(*args, **key)
def semilogx(*args, **key):
    return gca().semilogx(*args, **key)
def semilogy(*args, **key):
    return gca().semilogy(*args, **key)
def loglog(*args, **key):
    return gca().loglog(*args, **key)
def contour(*args, **key):
    return gca().contour(*args, **key)
contourf = contour
def yscale(*args, **key):
    return gca().yscale(*args, **key)
def xscale(*args, **key):
    return gca().xscale(*args, **key)
def text(*args, **key):
    return gca().text(*args, **key)
def legend(*args, **key):
    return gca().legend(*args, **key)
def show(*args, **key):
    return gca().show(*args, **key)
def xlabel(*args, **key):
    return gca().xlabel(*args, **key)
def ylabel(*args, **key):
    return gca().ylabel(*args, **key)
def xlim(*args, **key):
    return gca().xlim(*args, **key)
def ylim(*args, **key):
    return gca().ylim(*args, **key)
def title(*args, **key):
    return gca().title(*args, **key)
def suptitle(*args, **key):
    return gca().suptitle(*args, **key)
def grid(*args, **key):
    return gca().grid(*args, **key)

print("-- matplotlib not available, using fake plot instead")
print("-- no graphic will be available")

def gca():
    global gca0
    if gca0 is None:
        gca0 = fake()
    return gca0
gca0 = None


if __name__ == '__main__':
    print("This library creates fake entry methods allowing to fake minimum use of matplotlib.pyplot\neg:")
    print("""
    x=range(10)
    plot(x)
    f = figure()
    f.semilogx(x+x,label="with label")
    f.title('title')
    f.text(10,20,'infigure text')
    show()

    produces :
    """)
    x=range(10)
    plot(x)
    f = figure()
    f.semilogx(x+x,label="with label")
    f.title('title')
    f.text(10,20,'infigure text')
    show()
    