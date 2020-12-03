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
# # Display utility for 2D FTICR Spectra
#
# *This little utility allows to interactively explore large 2D FTICR-MS datasets.*
#
# You find here a simple interface the reads and displays the multiresolution 2D files created by `SPIKE` when processing 2D data-sets (usually called `xxx_mr.msh5`).
#
# It is based on the capabilities of both the `SPIKE` library and the `jupyter notebook` interface.
# Thanks to the technology behind, these extremely large files can be accessed rapidly, even on a laptop computer.
#
# This program supersedes the old `Visu2D` program, developped in `Qt`, which could no longer be maintained.
#
# *This is a work inprogress - additional utilities should come soon !*

# %% [markdown]
# ## To use it, 
# execute each executable cell (marked with the `In[]`) either by cliking on the Run icon on the top of the window, or by hitting *shift-Return* on the keyboard

# %% [markdown]
# ### Initialization of the environment
# the following cell should be run only once, at the beginning of the processing

# %%
from IPython.display import display, HTML, Markdown, Image
display(Markdown('## STARTING Environment...'))
# %matplotlib notebook
import spike.Interactive.INTER as I
from spike.Interactive.FTICR_INTER import MR, MR_interact
from spike.Interactive.ipyfilechooser import FileChooser
I.hidecode()

# %% [markdown]
# ### Choose the file
# Use `FileChooser()` to choose a file on your disk - The optional `base` argument, starts the exploration on a given location.
#
# 2D processed files are `*.msh5` files.

# %%
FC = FileChooser('/DATA',filetype='y*.msh5', mode='r')
display(FC)

# %% [markdown]
# the `MR` tool simply loads and describe the content of the file

# %%
import spike.Interactive.FTICR_INTER as FI
FI.SIZEMAX = 16*1024*1024

# %%
MR(FC.selected)

# %% [markdown]
# `MR_interact` loads and display the data-set.
#
# the 
#
# It can be called directly

# %%
MR_interact(FC.selected);

# %% [markdown]
# - 🔍 in  : zoom win
# - 🔍out : zoom out
# - ◀︎ : moves the zoom window left
# - ▶︎ : moves the zoom window right
# - ▲ : moves the zoom window up
# - ▼ : moves the zoom window down
# - ⍇ : back in zoom list
# - ⍈ : forward in zoom list
# - you can also directly enter the zoom coordinates, the click on Update
# - × : lower the levels used for the display
# - ÷ : raise the levels used for the display
# - ≡ : reset default levels
# - ⌘ : reset to initial view
#
# Note that the 2D files contain several version of the spectrum at different resolution - zooming and out may modify the look of the region you are looking to.
# Only the closest zoom contains the unbiaised verion of the spectrum.
#
# Some additional options are possible:
#
# - store the view into a python variable (we'll see other usage below)
# - store the created view into a python variable
# - define behaviour at start-up
# - overload the initial view

# %%
from importlib import reload
reload(FI) 
# complete initialisation, and storing the view into a python var
DI = FI.MR_interact(FC.selected,
                report=False,   # inhibits parameter printing
                show=False,     # does not display on start-up
                figsize=(15,15),# Size of initial display (in cm)
                Debug=False     # Enables live debugging if True
               )
DI.SIZEMAX = 32*1024*1024       # this parameter allows higher quality pictures (default is 8M)
DI._zoom = (380, 700, 380, 700)     # set the initial zoom view, in m/z (F1low , F1High , F2low , F2High)
DI.scale = 3.0                      # set the initial scale
DI.show()                           # and finally show

# %% [markdown]
# There is 1D extraction tool which is handy to examine carefully the details
#
# Just use your stored view and append `.I1D()` to it

# %%
DI.I1D()

# %% [markdown]
# ### to come
# - calibration
# - peak detection
# - superimposition
# - exctaction of arbitrary 1D 
# - locate artifacts due to harmonics

# %%
print (iPython.__version)

# %%
