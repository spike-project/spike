
import numpy as np

try:
    from mayavi import mlab
    ok = True
except:
    print('*** This plugin requires the installation of Mayavi (http://docs.enthought.com/mayavi/mayavi/installation.html)')
    ok = False

from spike import NPKError
from spike.NPKData import NPKData_plugin
import spike.NPKData as npk
import spike.FTICR
import numpy as np
def zoom3D(npkd, zoom, fontsize = 0.7, font = 'times', colormap = 'blue-red', showaxes = True):
    #assumes zoom is in format displayed by visu2D
    #x = axis2
    #y = axis 1
    #z = intensity
    z1lo = int(npkd.axis1.mztoi(zoom[1]))
    z1up = int(npkd.axis1.mztoi(zoom[3]))
    z2lo = int(npkd.axis2.mztoi(zoom[0]))
    z2up = int(npkd.axis2.mztoi(zoom[2]))
    d2 = npkd.extract(z1lo, z1up, z2lo, z2up)
    d3 = d2.get_buffer()
    d4 = np.array(d3)
    d5 = d4.transpose()
    zmax = np.amax(d5)
    xmin = zoom[2]
    xmax = zoom[0]
    ymin = zoom[3]
    ymax = zoom[1]
    zmin = 0
    mlab.figure(bgcolor=(1., 1., 1.), fgcolor = (0., 0., 0.))
    mlab.surf(d5, extent = [0, 1000, 0, 1000, 0, 1000], warp_scale = 'auto', colormap = colormap)
    ax = mlab.axes(x_axis_visibility = showaxes, y_axis_visibility = showaxes, z_axis_visibility = showaxes, xlabel = 'm/z (h)', ylabel = 'm/z (v)', zlabel = 'Intensity', ranges = [xmin, xmax, ymin, ymax, zmin, zmax], nb_labels = 5)
    ax.label_text_property.font_family = font
    ax.title_text_property.font_family = font
    ax.axes.font_factor = fontsize
if ok:
    NPKData_plugin("zoomwindow", zoom3D)