#!/usr/bin/env python 
# encoding: utf-8

"""
displaying using bokeh

"""

from __future__ import print_function
from spike import NPKError
from spike.NPKData import NPKData_plugin, parsezoom
import numpy as np

import bokeh
import bokeh.plotting as bk
from bokeh.models import Range1d
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

BokehMax = 60000000

########################################################################
def get_contour_data(ax):
    """
    Get informations about contours created by matplotlib.
    ax is the input matplotlob contour ax (cf. fig,ax produced by matplotlib)
    xs and ys are the different contour lines got out of the matplotlib. col is the color corresponding to the lines.
    """
    xs = []
    ys = []
    col = []
    isolevelid = 0
    for isolevel in ax.collections:
        isocol = isolevel.get_color()[0]
        thecol = 3 * [None]
        theiso = str(ax.collections[isolevelid].get_array())
        isolevelid += 1
        for i in range(3):
            thecol[i] = int(255 * isocol[i])
        thecol = '#%02x%02x%02x' % (thecol[0], thecol[1], thecol[2])
        for path in isolevel.get_paths():
            v = path.vertices
            x = v[:, 0]
            y = v[:, 1]
            xs.append(x.tolist())
            ys.append(y.tolist())
            col.append(thecol)
    return xs, ys, col

def bokeh_display(npkd, scale=1.0, autoscalethresh=3.0, absmax=None, show=False, title=None, label=None, xlabel="_def_", 
        ylabel="_def_", axis=None, image=False, mode3D=False, zoom=None, mpldic={}, 
        dbkdic={}, dfigdic={}, linewidth=1, color=None, plot_width=600, plot_height=400, 
        sizing_mode=None, redraw=False, tools="pan, box_zoom, box_select, reset, save"):
    """
        Display using bokeh instead of matplotlib
        
        scale       allows to increase the vertical scale of display
        absmax      overwrite the value for the largest point, which will not be computed
                display is scaled so that the largest point is first computed (and stored in absmax),
                and then the value at absmax/scale is set full screen 
        show        will call bk.show() at the end, allowing every declared display to be shown on-screen
                useless in ipython/jupyter notebook
        title       add a title to the bokeh plot
        label       add a label text to plot
        xlabel, ylabel    axes label (default is self.currentunit - use None to remove)
        axis        used as axis if present, axis length should match experiment length
                in 2D, should be a pair (xaxis,yaxis)
        image       if True, the function will generate the 2D NMR FID of data,
                if False (Default), the function present contour plots.
        mode3D      not implemented
        zoom        is a tuple defining the zoom window (left,right) or   ((F1_limits),(F2_limits))
                defined in the current axis unit (points, ppm, m/z etc ....)
        mpldic      a dictionnary passed as is to the matplotlib plot command 
        dbkdic      a dictionnary passed as is to populated the parameters of the bokeh graph 
        dfigdic     a dictionnary passed as is to populated the content of the bokeh figure
        linewidth   linewidth for the plots (useful for example when using seaborn)
        color       if 1D is the color of the curve, 
                if 2D FID is the palette name to be used,
                if 2D contour is the color set to be used by matplotlib.
        plot_width, plot_height  width and height of the plot
        sizing_mode if provided, resize plot according to the window chosen sizes.
                e.g. "scale_width", "scale_height", "scale_both"
        tools       a string containing the tools to be available for bokeh interactivity.
                e.g. "pan, box_zoom, box_select, reset, save" (see bokeh doc for more info)
    """
    #Settings global graphical parameters
    dbk = {}
    dbk['tools'] = tools
    dbk['title'] = title
    if sizing_mode is not None:
        dbk['sizing_mode'] = sizing_mode
    else:
        dbk['plot_width'] = plot_width
        dbk['plot_height'] = plot_height

    if npkd.dim == 1:
        if not absmax:  # absmax is the largest point on spectrum, either given from call, or handled internally
            if not npkd.absmax:     # compute it if absent
                #print "computing absmax...",
                npkd.absmax = np.nanmax( np.abs(npkd.buffer) )
        else:
            npkd.absmax = absmax
        mmin = -npkd.absmax/scale
        mmax = npkd.absmax/scale
        step = npkd.axis1.itype+1
        z1, z2 = parsezoom(npkd,zoom)
        while (z2-z1)//step > BokehMax:
            step += 2
        if axis is None:
            ax = npkd.axis1.unit_axis()
        else:
            ax = axis
#        fig.set_xscale(npkd.axis1.units[npkd.axis1.currentunit].scale)  # set unit scale (log / linear)
        # if npkd.axis1.units[npkd.axis1.currentunit].reverse:           # set reverse mode
        #     ax.invert_xaxis()
        if xlabel == "_def_":
            xlabel = npkd.axis1.currentunit
        if ylabel == "_def_":
            ylabel = "a.u."
        if xlabel is not None:
            dbk['x_axis_label'] = xlabel
        if ylabel is not None:
            dbk['y_axis_label'] = ylabel
        dbk.update(dbkdic)
        dbk['x_range'] = (npkd.axis1.itoc(z1), npkd.axis1.itoc(z2))

        p = bk.figure(**dbk)
        dfig = {}
        dfig['x'] = ax[z1:z2:step]
        dfig['y'] = npkd.buffer[z1:z2:step].clip(mmin,mmax)
        dfig['legend'] = label
        dfig['line_width'] = linewidth
        dfig['color'] = color
        dfig.update(dfigdic)
        p.line(**dfig)
        npkd.bokeh_fig = dfig
        npkd.bokeh_plot = p

    elif npkd.dim == 2 and image:
    # If image is set to True (Default False) - used to display 2D NMR FID as images
        bout1 = npkd.axis1.itoc(npkd.size1-1)
        bout2 = npkd.axis2.itoc(npkd.size2-1)
        deb1 =  npkd.axis1.itoc(0)
        deb2 =  npkd.axis2.itoc(0)
        dbk['x_range'] = (deb2, bout2)
        dbk['y_range'] = (deb1, bout1)
        if xlabel == "_def_":
            xlabel = npkd.axis2.currentunit
        if ylabel == "_def_":
            ylabel = npkd.axis1.currentunit
        if xlabel is not None:
            dbk['x_axis_label'] = xlabel
        if ylabel is not None:
            dbk['y_axis_label'] = ylabel
        dbk.update(dbkdic)
        p = bk.figure(**dbk)
        dfig = {}
        dfig['image'] = [npkd.get_buffer().real]
        if color is None:
            dfig['palette'] = bokeh.palettes.brewer['OrRd'][9][:-1] + bokeh.palettes.brewer['Blues'][9][::-1]
        else:
            dfig['palette'] = color #ex. "Magma256"
        dfig['x'] = deb2
        dfig['y'] = deb1
        dfig['dw'] = bout2
        dfig['dh'] = bout1
        dfig['legend'] = label
        dfig.update(dfigdic)
        p.image(**dfig)
        npkd.bokeh_fig = dfig
        npkd.bokeh_plot = p

    elif npkd.dim == 2 and not image:
    # When image is False, contour plots are generated for display.
        fig, ax = plt.subplots(figsize=(9,9)) #The figure is created as if it was a classical matplotlib figure
        z1lo, z1up, z2lo, z2up  = parsezoom(npkd,zoom)
        print(z1lo, z1up, z2lo, z2up)
        dbk['y_axis_type'] = npkd.axis1.units[npkd.axis1.currentunit].scale # set unit scale (log / linear) for bokeh
        dbk['y_range'] = (npkd.axis1.itoc(z1lo), npkd.axis1.itoc(z1up))
        dbk['x_axis_type'] = npkd.axis2.units[npkd.axis2.currentunit].scale # set unit scale (log / linear) for bokeh
        dbk['x_range'] = (npkd.axis2.itoc(z2lo), npkd.axis2.itoc(z2up))
        npkd.display(scale=scale, autoscalethresh=autoscalethresh, absmax=absmax, show=False, new_fig = False, axis=axis,
            zoom=zoom, title=None, figure=ax, color=color, mpldic=mpldic)
        if npkd.axis1.units[npkd.axis1.currentunit].reverse:
                ax.invert_yaxis()
        if npkd.axis2.units[npkd.axis2.currentunit].reverse:
                ax.invert_xaxis()
        if xlabel == "_def_":
            xlabel = npkd.axis2.currentunit
        if ylabel == "_def_":
            ylabel = npkd.axis1.currentunit
        if xlabel is not None:
            dbk['x_axis_label'] = xlabel
        if ylabel is not None:
            dbk['y_axis_label'] = ylabel
        dbk.update(dbkdic)
        xs, ys, col = get_contour_data(ax)
        if redraw:
            npkd.bokeh_fig['xs'] = xs
            npkd.bokeh_fig['ys'] = ys
        else:
            p = bk.figure(**dbk)
            xs, ys, col = get_contour_data(ax)
            dfig = {}
            dfig['xs'] = xs
            dfig['ys'] = ys
            dfig['color'] = col
            dfig.update(dfigdic)
            p.multi_line(**dfig)
            npkd.bokeh_fig = dfig
            npkd.bokeh_plot = p
        del fig, ax
    if show:
        bk.show(npkd.bokeh_plot)
        return npkd

NPKData_plugin("bokeh", bokeh_display)