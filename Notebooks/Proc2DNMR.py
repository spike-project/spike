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

# %% [markdown]
#
# # 2D NMR Processing and Display
#
# a simplified environment for processing 2D Bruker NMR datasets with `SPIKE`.
#
# Run each python cell in sequence by using the ⇥Run button above (or typing *shift* Enter).
#
# Cells are meant to be used in order, taking you to the complete analysis, but you can go back at any time.
#
# The SPIKE code used for processing is visible in the cells, and can be used as a minimal tutorial.
#
# ***Remark*** *to use this program, you should have installed the following packages:*
#
# - *a complete scientific python environment ( tested with python 3.7 - [anaconda](https://www.anaconda.com/))*
# - [`spike`](https://www.bitbucket.org/delsuc/spike) ( *version 0.99.21 minimum* )
# - [`ipywidgets`](https://ipywidgets.readthedocs.io/en/latest/)  ( *tested with version 7.6* )
# - [`ipyml`](https://github.com/matplotlib/jupyter-matplotlib)  ( *adds interactivity in the notebook* )
#
# ## Initialization
# the following cell should be run only once, at the beginning of the processing

# %%
# load all python and interactive tools
from __future__ import print_function, division
from IPython.display import display, HTML, Markdown, Image
display(Markdown('## STARTING Environment...'))
import matplotlib.pyplot as plt
# %matplotlib widget
import os.path as op
import spike
from spike.File.BrukerNMR import Import_2D
from spike.Interactive import INTER as I
from spike.Interactive import INTER_2D as I2D
from spike.Interactive.ipyfilechooser import FileChooser
print("\nInteractive module version,",I.__version__)
from datetime import datetime
print('Run date:', datetime.now().isoformat() )
I.initialize()
display(Markdown('## ...program is Ready'))
from importlib import reload  # the two following lines are debugging help

# configurable items
plt.rcParams['figure.figsize'] = (8,4)   # (X,Y) default figure size
I.Activate_Wheel = True                  # True/False wheel control in the graphic cells 

# %% [markdown]
# ---
#
# ### Choose the file
# The `FileChooser()` tool creates a dialog box which allows to choose a file on your disk
#
# - use the `Select` button
# - modify the ( *optional* ) `path` argument, to start the exploration on a given location
# - After the selection, the selected filename is found in `FC.selected`

# %%
FC = FileChooser(path='/home/mad/Documents/',filename='ser')
display(FC)

# %% [markdown]
# ### Import dataset
#
# This is simply done with the `Import_2D()` tool, which returns a `SPIKE` object.
#
# We store the dataset into a variable called d2. 

# %%
print('Reading file ',FC.selected)
d2 = Import_2D(FC.selected)
d2.filename = FC.selected
d2.pulprog = d2.params['acqu']['$PULPROG']
print(d2)                                      # print() of the dataset shows a summary of the parameters
print(d2.pulprog)
display(HTML('<b>title: </b>'+ d2.params['acqu']['title']))    # d1.params is a dictionary which contains the whole 'acqu' and 'proc' Bruker parameters
#d2.set_unit('sec').display(title="%s %s"%(FC.nmrname,d2.pulprog), scale='auto')
plt.figure()
plt.imshow(d2.get_buffer().real, cmap="seismic");

# %% [markdown]
# *Check Proc1DNMR NoteBook for documentation on spectral display*
#
# ## Basic Processing
# We are going to use a basic processing set-up, check the documentation for advanced processing
#
# ### Fourier Transform - modulus mode!

# %%
D2 = d2.copy() # copy the imported data-set to another object for processing
# bk_ftF2 and bk_ftF1 (define in the Bruker plugin) find which FT to apply depending on FnMODE
D2.apod_sin(maxi=0.5,axis='F2').zf(1,2).bk_ftF2()  # chaining  apodisation - zerofill - FT
D2.apod_sin(maxi=0.5,axis='F1').zf(2,1).bk_ftF1()  # chaining  apodisation - zerofill - FT
D2.modulus().set_unit('ppm').rem_ridge()
#D2.display(scale="auto", autoscalethresh=100.0, title="%s %s"%(FC.nmrname,d2.pulprog))  # chain  set to ppm unit - and display
I2D.Show2D(D2)

# %% [markdown]
# ### Advanced Phase sensitive processing
#
#

# %%
D2ph = d2.copy() # copy the imported data-set to another object for processing
# bk_ftF2 and bk_ftF1 (define in the Bruker plugin) find which FT to apply depending on FnMODE
D2ph.apod_sin(maxi=0,axis='F2').zf(1,2).bk_ftF2().bk_pk()  # chaining  apodisation - zerofill - FT - phase
D2ph.apod_sin(maxi=0,axis='F1').zf(2,1).bk_ftF1()  # chaining  apodisation - zerofill - FT
D2ph.set_unit('ppm').rem_ridge()
#D2.display(scale="auto",  autoscalethresh=6.0, title="%s %s"%(FC.nmrname,d2.pulprog))  # chain  set to ppm unit - and display
reload(I2D)
S = I2D.Show2D(D2ph)           # note that you can create the view and store it in a var wiothout displaying it
S.negview.value = True         # and set some parameters (here negative view) before display
S

# %% [markdown]
# ### Rephasing
#
# Use the sliders to adjust the phase parameters,   the pivot can be set with a right click on the spectrum
# Top and Side spectra are taken at the pivot level.
#
#

# %%
reload(I2D)
I2D.Phaser2D(D2ph)

# %% [markdown]
# # An interactive Display
# - show2D +
# - projections
# - color maps
# - title
#
# # todo
# - peak picking

# %%
F1slice = 3.3    # select a F1 (vertical) slice in current unit (here ppm) 
F2slice = 2.05   # select a F2 (horizontal) slice in current unit (here ppm) 

D2ph.col( D2ph.axis2.ctoi(F1slice)).display(title='F1 slice at F2=%.3f ppm'%(F1slice,))
D2ph.row( D2ph.axis1.ctoi(F2slice)).display(title='F2 slice at F1=%.3f ppm'%(F2slice,))

# %%
D2.set_unit('ppm')
b = I2D.Show2D(D2.copy().real(axis='F1').real(axis='F2'), title="%s %s"%(FC.nmrname,d2.pulprog))

# %% [markdown]
# ## Save on disk

# %%
D2.save('example1.gs2')

# %% [markdown]
# # The following entries or not finished yet

# %% [markdown]
# ## Peak-Picker
# - moving the threshold determines the minimum peak intensity
# - peaks are searched only in the selected zoom window

# %% [markdown]
# ## Export a bucket list

# %%
# adapt the parameters below
Zoom = ((0.5,8),(0.5,8))                    # zone to bucket       - in ppm
BucketSize = (0.1,0.1)                 # width of the buckets - in ppm
Output = 'screen'                   # 'screen'  or  'file'  determines output
BucketFileName = 'bucket.csv'     #  the filename if Output (above) is 'file'  - don't forget the .csv extension.

# %%
# the following cell executes the bucketing
if Output == 'file':
    with open(BucketFileName,'w') as F:
        D2.bucket2d(zoom=Zoom, bsize=BucketSize, file=F)
    print('buckets written to %s\n'%op.realpath(BucketFileName))
else:
    D2.bucket2d(zoom=Zoom, bsize=BucketSize)

# %% [markdown]
# *Tools in this page is under intensive development - don't expect things to remain as they are.*

# %%
