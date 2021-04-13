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
# # FTICR-MS Processing and Display
#
# a simplified environment for processing 1D Bruker FTICR datasets with `SPIKE`
#
# Run each python cell in sequence by using the ⇥Run button above (or typing *shift* Enter).
#
# Cells are meant to be used in order, taking you to the complete analysis, but you can go back at any time.
#

# %% [markdown]
# ### Initialization of the environment
# the following cell should be run only once, at the beginning of the processing

# %%
# load all python and interactive tools
from __future__ import print_function, division
from IPython.display import display, HTML, Markdown, Image
display(Markdown('## STARTING Environment...'))
# %matplotlib notebook
import os.path as op
import numpy as np
import spike
from spike.Interactive import INTER as I
from spike.Interactive import FTICR_INTER as FI
from spike.Interactive.ipyfilechooser import FileChooser
from spike.File import BrukerMS
display(Markdown('## ... program is Ready'))
I.hidecode()

# %% [markdown]
# ### Choose the file
# Use `FileChooser()` to choose a file on your disk - The optional `base` argument, starts the exploration on a given location.
#
# Bruker files are named `fid` and are contained in a `*.d` directory.

# %%
FC = FileChooser(path='/DATA', filetype='fid')
display(FC)

# %% [markdown]
# (After the selection, the selected filename is found in the `FC.selected` variable)
# ### Import dataset
#

# %%
# This is simply done with the `Import_1D()` tool, which returns a `SPIKE` object.
# We store it into a variable, evaluating the variable show a summary of the dataset. 
print('Reading file ',FC.selected)
d1 = BrukerMS.Import_1D(FC.selected)
d1.filename = FC.selected
d1.set_unit('sec').display(title=FC.selected_path+" transient")
d1

# %% [markdown]
# In the current set-up, the figure can be explored *(zoom, shift, resize, etc)* with the jupyter tools displayed  below the dataset.
# The figure can also be saved as a `png` graphic file.
#
# At anytime, the figure can be frozen by clicking on the blue button on the upper right corner, just rerun the cell for changing it.

# %% [markdown]
# ### Compute Spectrum
#
# many processing methods are available, they can be either applied one by one, or piped by chaining them.
#
# Here we are chaining  apodisation - zerofill - FT - modulus
#
# then setting to `m/z` unit (`Hz` and `points` also available) - finally `display()` is used to display the dataset.
#

# %%
D1 = d1.copy() # copy the imported data-set to another object for processing
D1.kaiser(4).zf(4).rfft().modulus() # kaiser(4) is an apodisation well adapted to FTICR, slightly more resolution than hamming(
D1.set_unit('m/z').display(title=FC.selected_path)  # set to ppm unit - and display

# %% [markdown]
# ### Peak Detection
# The following is used to perform an interactive peak picking, and output the result
# Use the cursor to choose the sensitivity of the picker.

# %%
pkname = op.join(op.dirname(FC.selected),"peak_list.txt")
import spike.plugins.Peaks
print("number of on-screen peaks is limited to %d - zoom in with the top tool to be sure of visualizing them all"%(spike.plugins.Peaks.NbMaxDisplayPeaks))
FI.MSPeaker(D1, pkname);

# %% [markdown]
# ### Calibration
# The calibration used by SPIKE is based on a 2 or 3 parameters equation :
# $$
# f = \frac{A}{m/z} - B + \frac{C}{(m/z)^2}
# $$
# You can change them below:
#

# %%
FI.Calib(D1);

# %% [markdown]
# ---
#
# ### Calibration on reference peaks
# #### *To come soon !*
#

# %% [markdown]
# ### Save processed data
# You can save a dataset, two formats are available:
#
# - Native SPIKE format, `*.msh5` where all informations are stored - run the following cell

# %%
msh5name = op.join(op.dirname(FC.selected),"SpikeProcessed.msh5")
D1.save_msh5(msh5name, compressed=True)
print("File is stored as %s"%msh5name)

# %% [markdown]
# - Or `cvs` format, with only the spectrum (for the peak list, see above) - ( *be carefull this file can be very big*)

# %%
csvname = op.join(op.dirname(FC.selected),"SpikeProcessed.csv")
D1.save_csv(csvname)
print("File is stored as %s"%csvname)

# %% [markdown]
# ### superimpose spectra
# you can superimpose several spectra stored as `.msh5` files in order to compare them

# %%
from IPython.display import display
SP = FI.SuperImpose(base='/DATA',N=6).Show()

# %% [markdown]
# *the following cell display the colormap used here*

# %%
import matplotlib.pylab as plt
plt.figure()
for i ,c in enumerate(FI.Colors):
    plt.plot([0,1],[0,i],'-',color=c,label=c)
plt.legend()

# %%
