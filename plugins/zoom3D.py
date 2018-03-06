"""Code to display 2D data-sets in interactive 3D, using Mayavi

Original code from Alice Lynch - mai 2016

adapted my M-A Delsuc
"""
from __future__ import print_function

import numpy as np

try:
    from mayavi import mlab
    ok = True
except:
    raise Exception('This plugin requires the installation of Mayavi (http://docs.enthought.com/mayavi/mayavi/installation.html) )')
    print('*** The zoom3D plugin requires the installation of Mayavi (http://docs.enthought.com/mayavi/mayavi/installation.html)')
    ok = False

from spike import NPKError
from spike.NPKData import NPKData_plugin
import spike.NPKData as npk

def zoom3D(npkd, zoom, fontsize=0.7, font='times', colormap='blue-red', showaxes=True):
    """
    use the zoom region and display it in interactive 3D
    
    zoom : f1lo, f1up, f2lo, f2up - expressed in current unit

    remark, use a small zoom region !

    x = axis 2
    y = axis 1
    z = intensity
    """
    # z1lo = int(npkd.axis1.ctoi(zoom[0]))
    # z1up = int(npkd.axis1.ctoi(zoom[1]))
    # z2lo = int(npkd.axis2.ctoi(zoom[2]))
    # z2up = int(npkd.axis2.ctoi(zoom[3]))
    z1lo, z1up, z2lo, z2up = npk.parsezoom(npkd, zoom)
    #print (z1lo, z1up, z2lo, z2up)
    # d2 = npkd.extract(*zoom)  # extract modifies the data - 
    # d3 = d2.get_buffer()
    # d4 = np.array(d3)
    # d5 = d4.transpose()
    # get_buffer() loads everything in memory - so we'll do it line by line
    d5 = np.zeros((z2up-z2lo+1, z1up-z1lo+1))  # tranposed matrix compared to npkd
    for i in range(z2lo, z2up+1):
        cc = npkd.col(i)              # going column wise is probably faster ...
        d5[i-z2lo,:] = cc[z1lo:z1up+1]  # taking a slice out of a npkdata returns a np.array
    zmax = np.amax(d5)
    zmin = np.amin(d5) # 0  - some data-set are negative
    xmin = zoom[2]
    xmax = zoom[3]
    ymin = zoom[0]
    ymax = zoom[1]
    mlab.figure(bgcolor=(1., 1., 1.), fgcolor = (0., 0., 0.))
    mlab.surf(d5, extent = [0, 1000, 0, 1000, 0, 1000], warp_scale='auto', colormap=colormap)
    ax = mlab.axes(x_axis_visibility=showaxes, y_axis_visibility=showaxes, z_axis_visibility=showaxes, xlabel="F2 "+npkd.axis2.currentunit, ylabel="F1 "+npkd.axis1.currentunit, zlabel='Intensity', ranges=[xmin, xmax, ymin, ymax, zmin, zmax], nb_labels=5)
    ax.label_text_property.font_family = font
    ax.title_text_property.font_family = font
    ax.axes.font_factor = fontsize

if ok:
    NPKData_plugin("zoomwindow", zoom3D)