# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 1.2.4
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# ** *This notebook is unfinished and still under development (but not too bad though ;-)* **
#
# # 1D NMR Processing and Display
#
# a simplified environment for processing 1D Bruker NMR datasets with `SPIKE`.
#
# Run each python cell in sequence by using the ⇥Run button above (or typing *shift* Enter).
#
# Cells are meant to be used in order, taking you to the complete analysis, but you can go back at any time.
#
# The SPIKE code used for processing is visible in the cells, and can be used as a minimal tutorial.
#
# `reverse_scroll` inverse the direction of the mouse wheel, whether it is `True` or `False`
#
# ***Remark*** *to use this program, you should have installed the following packages:*
#
# - *a complete scientific python environment* ( *tested with python 3.6 - [anaconda](https://www.anaconda.com/) does probably not work in python 2.7 *)
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
# %matplotlib widget
import os.path as op
import spike
from spike.File.BrukerNMR import Import_1D
from spike.Interactive import INTER as I
from spike.Interactive.ipyfilechooser import FileChooser
display(Markdown('## ...program is Ready'))
from importlib import reload  # the two following lines are debugging help
reload(I)                   # and can be removed safely when in production
I.hidecode()

# %% [markdown]
# ### Choose the file
# The `FileChooser()` tool creates a dialog box which allows to choose a file on your disk
#
# - use the `Select` button
# - modify the ( *optional* ) `path` argument, to start the exploration on a given location
# - After the selection, the selected filename is found in `FC.selected`

# %%
FC = FileChooser(path='/DATA/',filename='fid')
display(FC)

# %% [markdown]
# ### Import dataset
#
# This is simply done with the `Import_1D()` tool, which returns a `SPIKE` object.
#
# We store the dataset into a variable, typing the variable name shows a summary of the dataset. 

# %%
print('Reading file ',FC.selected)
d1 = Import_1D(FC.selected)
d1.filename = FC.selected
d1.set_unit('sec').display(title=FC.nmrname+" fid")

# %% [markdown]
# In the current set-up, the figure can be explored *(zoom, shift, resize, etc)* with the jupyter tools displayed  below the dataset.
# The figure can also be saved as a `png` graphic file.
#
# For more interactivity - see below.
#
# ## Basic Processing
# We are going to use a basic processing set-up, check the documentation for advanced processing
#
# ### Fourier Transform

# %%
D1 = d1.copy() # copy the imported data-set to another object for processing
D1.apod_em(0.3).zf(4).ft_sim().bk_corr().apmin()  # chaining  apodisation - zerofill - FT - Bruker correction - autophase
D1.set_unit('ppm').display(title=FC.nmrname)  # chain  set to ppm unit - and display

# %% [markdown]
# <hr/>
#
# **Following steps are optional**
#
# ### rephasing
# If is is required use the interactive phaser
#
# Use `scale` (or mouse wheel) and the zoom box to tune the display; then use `P0, P1, pivot` to optimize the phase.
#
# `pivot` is where the $1^{st}$ order correction is not acting - and can be moved with the slider or by clicking on the spectrum
#
# Once finished, click on `Done`

# %%
reload(I)
I.Phaser1D(D1, reverse_scroll=True);

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
reload(I)
I.baseline1D(D1, reverse_scroll=True);

# %% [markdown]
# ## Peak-Picker
# - first detect peaks (<span style=color:red>in red</span>) in the zoom window by moving the `threshold`
# - `Add` the detected peaks to the permanent peak list (<span style=color:blue>in blue</span>)
# - recalibrate the spectrum in the `calibration` tab by clicking on a peak and setting the correct ppm value
# - get the peak-list in the `Peak Table` tab

# %%
reload(I)
ph = I.NMRPeaker1D(D1, reverse_scroll=True);

# %% [markdown]
# ## Integrate
# Integration zones are computed from the peaks detected with the Peak-Picker above **required**
#

# %%
reload(I)
D1.real()
ii = I.NMRIntegrate(D1, reverse_scroll=True);

# %% [markdown]
# ## Interactive composite display
# Convenient to set-up your own figure
# (spectral superposition is not operational)

# %%
reload(I)
s = I.Show1Dplus(D1, title=FC.nmrname, reverse_scroll=True);

# %% [markdown]
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
