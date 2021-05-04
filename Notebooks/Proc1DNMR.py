# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown] class="test_MAD"
# # 1D NMR Processing and Display
#
# a simplified environment for processing 1D Bruker NMR datasets with `SPIKE`.
#
# Run each python cell in sequence by using the ⇥Run button above (or typing *shift* Enter).
#
# Cells are meant to be used in order, taking you to the complete analysis, but you can go back at any time.
#
# The SPIKE code used for processing is visible in the cells, and can be used as a minimal tutorial.
# You can hide it when done to present a clean NoteBook.
#
# `reverse_scroll` inverse the direction of the mouse wheel, whether it is `True` or `False`
#
# ***Remark*** *to use this program, you should have installed the following packages:*
#
# - *a complete scientific python environment ( tested with python 3.7 / [anaconda](https://www.anaconda.com/)  with no support for python 2.7 )*
# - [`spike`](https://www.bitbucket.org/delsuc/spike) ( *version 0.99.21 minimum* )
# - [`ipywidgets`](https://ipywidgets.readthedocs.io/en/latest/)  ( *tested with version 7.6* )
# - [`ipympl`](https://github.com/matplotlib/jupyter-matplotlib)  ( *adds interactivity in the notebook* )
#
# ## Initialization
# the following cell should be run only once, at the beginning of the processing

# %%
# load all python and interactive tools - should be run only once
from IPython.display import display, HTML, Markdown, Image
display(Markdown('## STARTING Environment...'))
import matplotlib as mpl
# %matplotlib widget
import spike
from spike.File.BrukerNMR import Import_1D
from spike.Interactive import INTER as I
from spike.Interactive.ipyfilechooser import FileChooser
print("\nInteractive module version,",I.__version__)
from datetime import datetime
print('Run date:', datetime.now().isoformat() )
I.initialize()
display(Markdown('## ...program is Ready'))
from importlib import reload  # this line is debugging help

# configurable items
mpl.rcParams['figure.figsize'] = (8,4)   # (X,Y) default figure size
I.Activate_Wheel = True                  # True/False wheel control in the graphic cells 

# %% [markdown]
# ### Choose the file
# The `FileChooser()` tool creates a dialog box which allows to choose a file on your disk
#
# - use the `Select` button, click on the directory to move around `..` is for going back in the directory tree
# - click on the file to select it
# - modify the ( *optional* ) `path` argument, to start the exploration on a given location
# - changing to `path='.'` will start the browsing from the current directory 
# - After the selection, the selected filename is found in `FC.selected`

# %%
# FileChooser
FC = FileChooser(path='/DATA/',filename='fid')
display(FC)

# %% [markdown]
# ### Import dataset and display FID
#
# This is simply done with the `Import_1D()` tool, which returns a `SPIKE` object.
#
# We store the dataset into a variable, here called d1. 

# %%
print('Reading file ',FC.selected)
d1 = Import_1D(FC.selected)                    # Import_1D creates a SPIKE NMRData object, from which everything is available
d1.set_unit('sec')                             # it can be acted upon
d1.filename = FC.selected                      # and be extended at will
print(d1)                                      # print() of the dataset shows a summary of the parameters
display(HTML('<b>title: </b>'+ d1.params['acqu']['title']))    # d1.params is a dictionary which contains the whole 'acqu' and 'proc' Bruker parameters
I.Show1D(d1, title=FC.selected)
# alternatively you can use the low-level tool below:
# d1.display()  # with many options, plus access to matplotlib details

# %% [markdown]
# Spectra can be interactively explored with the jupyter tools displayed  on the side of the dataset:
#
# - zoom with <button class="jupyter-matplotlib-button jupyter-widgets jupyter-button" href="#" title="Zoom to rectangle" style="outline: currentcolor none medium;"><i class="center fa fa-square-o"></i></button>
# - shift and resize
# <button class="jupyter-matplotlib-button jupyter-widgets jupyter-button" href="#" title="Pan axes with left mouse, zoom with right" style="outline: currentcolor none medium;"><i class="center fa fa-arrows"></i></button>
#  (with left and right click)
# - <button class="jupyter-matplotlib-button jupyter-widgets jupyter-button" href="#" title="Back to previous view" style="outline: currentcolor none medium;"><i class="center fa fa-arrow-left"></i></button>
# and
# <button class="jupyter-matplotlib-button jupyter-widgets jupyter-button" href="#" title="Forward to next view" style="outline: currentcolor none medium;"><i class="center fa fa-arrow-right"></i></button>
# allow to navigate in the zoom history
#
# - <button class="jupyter-matplotlib-button jupyter-widgets jupyter-button" href="#" title="Download plot" style="outline: currentcolor none medium;"><i class="center fa fa-fw fa-floppy-o"></i></button> is used to store a `png` graphic file of the current display.
#
# The drawing zone can be resized using the little grey triangle on the lower-right corner
#
# In addition, the `I.Show1D()` is a high level tool which adds a scale slider, a
# <button class="p-Widget jupyter-widgets jupyter-button widget-button" style="width: 80px;" >Save figure</button>
#  button, to store a high quality pdf version,
#  and a
# <button class="p-Widget jupyter-widgets jupyter-button widget-button mod-success" style="width: 80px;" >Done</button>
# button to freeze the picture.
#
# ## Basic Processing
# The following cell applies a basic processing, check the documentation for more advanced processing

# %%
# Basic Processing
LB = 0.1                        # you can adapt LB to your means, in Hz
D1 = d1.copy()                  # copy the imported data-set to another object for processing
D1.apod_em(LB).zf(4).ft_sim().bk_corr().apmin()  # chaining  apodisation - zerofill - FT - Bruker correction - autophase
D1.set_unit('ppm')              # set to ppm unit ('Hz' and 'point' also available)
                                # all Spike command can be pipelined at will - these 3 lines could be piped as one.
I.Show1D(D1, title=FC.nmrname)  #  and display

# %% [markdown]
# <hr/>
#
# **Following steps are for special operations**
#
# ### rephasing
# If the spectrum requires rephasing, use the interactive phaser down here.
#
# Use `scale` (or mouse wheel) and the zoom box to tune the display; then use `P0, P1, pivot` to optimize the phase.
#
# `pivot` is where the $1^{st}$ order correction is not acting - and can be moved with the slider or by right-clicking on the spectrum
#
# Once finished, click on `Done`

# %%
# rephasing
reload(I)
I.Phaser1D(D1, title=FC.nmrname)

# %% [markdown]
# ### Baseline correction
# A simple interactive baseline correction tool
#
# Choose positions on the baseline of the spectrum with the `select` slider or  by clicking on the baseline.
# `Add` a control point and see its effect either on the spectrum, or the computed baseline.
#
# You can also `Rem`ove the control point closest to the selector.
# You can also try the `Auto` button for a set of selector points, a set you can optimize by adding and removing points.

# %%
# Baseline Correction
I.baseline1D(D1)

# %% [markdown]
# ## Peak-Picker
# - first detect peaks (<span style=color:red>in red</span>) in the zoom window by moving the `threshold`
# - `Add` the detected peaks to the permanent peak list (<span style=color:blue>in blue</span>)
# - recalibrate the spectrum in the `calibration` tab by clicking on a peak and setting the correct ppm value
# - get the peak-list in the `Peak Table` tab

# %%
# Peak Picker
reload(I)
I.NMRPeaker1D(D1)

# %% [markdown]
# ## Integrate
# Integration zones are computed from the peaks detected with the Peak-Picker above **required**
#

# %%
# Integration
reload(I)
D1.real()
I.NMRIntegrate(D1)

# %% [markdown]
# ## Interactive composite display
# Convenient to set-up your own figure
# (spectral superposition is not operational)
#
# label y pos -1

# %%
# Composite display
reload(I)
I.Show1Dplus(D1, title=FC.nmrname)

# %% [markdown]
# ---
# optional steps
#
# ## Save the data-set
# either as stand alone native SPIKE files, (there are other formats)

# %%
D1.save('example1.gs1')

# %% [markdown]
# or as a `csv` text file, - in which case, it is probably better to remove the imaginary part, not useful there.
#
# The file contains some basic informations in addition to the spectral data

# %%
D1.copy().real().save_csv('example.csv')

# %% [markdown]
# ### Save the peak list to a csv file

# %%
D1.pk2pandas().to_csv('peaklist.csv')

# %% [markdown]
# ### Save the integrals to a csv file

# %%
D1.integrals.to_pandas().to_csv('integrals.csv')

# %% [markdown]
# ---
#
# This part adds the bucket list tool 
#
# ## Export a buckelist

# %%
# adapt the parameters below
Zoom = (0.5,8)                    # zone to bucket       - (start, end) in ppm
BucketSize = 0.04                 # width of the buckets - in ppm
Output = 'screen'                 # 'screen'  or  'file'  determines output
BucketFileName = 'bucket.csv'     #  the filename if Output (above) is 'file'  - don't forget the .csv extension.

# %%
# the following cell executes the bucketing
if Output == 'file':
    with open(BucketFileName,'w') as F:
        D1.bucket1d(zoom=Zoom, bsize=BucketSize, pp=True, file=F)
    print('buckets written to %s\n'%op.realpath(BucketFileName))
else:
    D1.bucket1d(zoom=Zoom, bsize=BucketSize, pp=True);

# %% [markdown]
# *Tools in this page is under intensive development - things are going to change rapidly.*
#
# to come/finish:
#
# - spectral superposition
# - annotations

# %%
