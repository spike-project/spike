# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.3.4
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# ** *This notebook is unfinished and still under development* **
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
# - *a complete scientific python environment* ( *tested with python 3.6 - [anaconda](https://www.anaconda.com/) 
#  but it should also work in python 2.7*)
# - [`spike`](https://www.bitbucket.org/delsuc/spike) ( *version 0.99.9 minimum* )
# - [`ipywidgets`](https://ipywidgets.readthedocs.io/en/latest/)  ( *tested with version 7.1* )
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
display(Markdown('## ...program is Ready'))
from importlib import reload  # the two following lines are debugging help
reload(I)                   # and can be removed safely when in production
I.hidecode(message="")
I.hidedoc()

# %% [markdown]
# ### Choose the file
# The `FileChooser()` tool creates a dialog box which allows to choose a file on your disk
#
# - use the `Select` button
# - modify the ( *optional* ) `path` argument, to start the exploration on a given location
# - After the selection, the selected filename is found in `FC.selected`

# %%
FC = FileChooser(path='/DATA/',filename='ser')
display(FC)

# %% [markdown]
# ### Import dataset
#
# This is simply done with the `Import_2D()` tool, which returns a `SPIKE` object.
#
# We store the dataset into a variable, typing the variable name shows a summary of the dataset. 

# %%
print('Reading file ',FC.selected)
d2 = Import_2D(FC.selected)
d2.filename = FC.selected
d2.pulprog = d2.params['acqu']['$PULPROG']
print (d2.params['acqu']['title'])
#d2.set_unit('sec').display(title="%s %s"%(FC.nmrname,d2.pulprog), scale='auto')
plt.imshow(d2.get_buffer().imag, cmap="seismic"); #  "Accent")

# %% [markdown]
# In the current set-up, the figure can be explored *(zoom, shift, resize, etc)* with the jupyter tools displayed  below the dataset.
# The figure can also be saved as a `png` graphic file.
#
# At anytime, the figure can be frozen by clicking on the blue button on the upper right corner, just rerun the cell to make it interactive again.
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
D2.display(scale="auto", autoscalethresh=100.0, title="%s %s"%(FC.nmrname,d2.pulprog))  # chain  set to ppm unit - and display

# %% [markdown]
# `I.Show2D` is a convinient tool to explore 2D's
#
# you can 

# %% [markdown]
# ### Advanced Phase sensitive processing
#
#

# %%
D2 = d2.copy() # copy the imported data-set to another object for processing
# bk_ftF2 and bk_ftF1 (define in the Bruker plugin) find which FT to apply depending on FnMODE
D2.apod_sin(maxi=0,axis='F2').zf(1,2).bk_ftF2().bk_pk()  # chaining  apodisation - zerofill - FT - phase
D2.apod_sin(maxi=0,axis='F1').zf(2,1).bk_ftF1()  # chaining  apodisation - zerofill - FT
D2.set_unit('ppm').rem_ridge()
D2.display(scale="auto",  autoscalethresh=6.0, title="%s %s"%(FC.nmrname,d2.pulprog))  # chain  set to ppm unit - and display

# %% [markdown]
# ### Rephasing
# ( *This is a temporary tool* )

# %%
reload(I2D)
I2D.Phaser2D(D2)

# %% [markdown]
# # An interactive Display

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
# *Tools in this page is under intensive development - things are going to change rapidly.*
