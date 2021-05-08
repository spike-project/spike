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
# # 2D DOSY NMR Processing and Display
#
#
# a simplified environment for processing DOSY experiments from Bruker NMR datasets.
# The analysis is performed using the PALMA algorithm.
#
# ## Principle
# The PALMA algorithm performs the Inverse Laplace Transform of each columns of the 2D DOSY experiment.
# This inversion is performed by optimizing a target function which maximizes the Entropy of the spectrum, and minimizes its $\ell_1$ norm.
#
# More details are found in the original publication [DOI: 10.1039/c6an01902a](http://doi.org/10.1039/c6an01902a)
# or the [WEB server page](http://palma.labo.igbmc.fr/about.html) on which you will also find [advises](http://palma.labo.igbmc.fr/about.html) to obtain the best possible experiments. 
#
# When using this program, please cite the following reference:
#
# A.Cherni, É.Chouzenoux, and M.-A.Delsuc. (2017) "PALMA, an improved algorithm for DOSY signal processing".
# *Analyst*, **142**(5), p772-779 
#

# %% [markdown]
# ## Usage
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
#
# ## Initialization
# the following cell should be run only once, at the beginning of the processing

# %%
# load all python and interactive tools
from __future__ import print_function, division
from importlib import reload
from IPython.display import display, HTML, Markdown, Image, Javascript
display(Markdown('## STARTING Environment...'))
# %matplotlib notebook
import os.path as op
import spike
from spike.File.BrukerNMR import Import_2D
from spike.plugins.NMR.PALMA import Import_DOSY
from spike.Interactive import INTER as I
from spike.Interactive.ipyfilechooser import FileChooser
display(Markdown('## ...program is Ready'))
I.hidecode()

# %% [markdown]
# ### Choose the file
# The `FileChooser()` tool creates a dialog box which allows to choose a file on your disk
#
# - use the `Select` button
# - modify the ( *optional* ) `path` argument, to start the exploration on a given location
# - choose the `ser  file of your DOSY experiment
#     - *note, the* `difflist` *file should be present*
# - After the selection, the selected filename is found in `FC.selected`

# %%
FC = FileChooser(path='/home/mad',filename='ser')
display(FC)

# %% [markdown]
# ### Import dataset
#
# This is simply done with the `Import_2D()` tool, which returns a `SPIKE` object.
#
# We store the dataset into a variable, typing the variable name shows a summary of the dataset. 

# %%
if not op.exists(op.join(FC.selected_path,'difflist')):
    I.jsalert('No difflist file present - is this a DOSY experiment ?')
print('Reading file ',FC.selected)
d2 = Import_DOSY(FC.selected)
d2.filename = FC.selected
d2.pulprog = d2.params['acqu']['$PULPROG']
print (d2.params['acqu']['title'])
d2.axis2.currentunit = 'sec'
d2.axis1.currentunit = 'points'
d2.display(title="%s/ser %s"%(FC.nmrname,d2.pulprog), scale='auto')

# %% [markdown]
# In the current set-up, the figure can be explored *(zoom, shift, resize, etc)* with the jupyter tools displayed  below the dataset.
# The figure can also be saved as a `png` graphic file.
#
# At anytime, the figure can be frozen by clicking on the blue button on the upper right corner, just rerun the cell to make it interactive again.
#
# ## Data Preparation
# We are going to use a basic processing set-up, check the documentation for advanced processing
#
# ### Fourier Transform

# %%
D2 = d2.copy() # copy the imported data-set to another object for processing
# bk_ftF2 and bk_ftF1 (define in the Bruker plugin) find which FT to apply depending on FnMODE
D2.apod_sin(maxi=0.1,axis='F2').zf(1,2).bk_ftF2().bk_pk()  # chaining  apodisation - zerofill - FT
D2.axis2.currentunit = 'ppm'
D2.display(scale="auto", autoscalethresh=6.0, title="%s %s"%(FC.nmrname,d2.pulprog))  # chain  set to ppm unit - and display

# %% [markdown]
# ### Optimize phase
# If required, rephase first rows, and the correction will be applied to the whole dataset.

# %%
R1 = D2.row(0) # copy the first line to another object for processing
for i in range(4):
    R1 += D2.row(i+1)
Ph = I.Phaser1D(R1)  # phase it
def loc(b):    # patch the on_done action o the phaser to propagate to DOSY
#    Ph.on_Apply(b)
    lp0, lp1 = Ph.ppivot() # get centered values
    if lp0 != 0 or lp1 != 0:
        D2.phase(lp0, lp1, axis=2)
        D2.display(scale="auto", autoscalethresh=6.0, title="%s %s"%(FC.nmrname,d2.pulprog))
Ph.button.on_click(loc)
display(Ph)

# %% [markdown]
# ### baseline correction
# A flat baseline centered on 0 is a requisite for a correct analysis.
# The following tool allows to point baseline positions on a extracted 1D, see the effect of the correction, and apply it to the whole 2D experiment.
#
# - use `add` and `rem` buttons to place points at the location of the selector.
# - the `auto` button populates with a first set of 8 positions

# %%
R1 = D2.row(0) # copy the imported data-set to another object for processing
for i in range(4):
    R1 += D2.row(i+1)
D2.projF2 = R1
b = I.baseline2D_F2(D2)

# %% [markdown]
# ### Verify your data

# %%
I.Show2D(D2)

# %% [markdown]
# ## PALMA processing
# First set you parameters
#
# Processing time is proportional to $N$ `x nbiter x finalsize`
#
# where $N$ is the number of processed columns (chosen with `miniSNR`)

# %%
# Diffusion axis
finalsize = 256  # The final of the DOSY on the diffusion axis
Dmin = 10       # the minimum diffusion coefficient (in µm2/sec) typical is 1-100
Dmax = 10000     # the maximum diffusion coefficient (in µm2/sec) typical is 5000-50000

# Processing
nbiter=1000      # number of iterations - the more the better (usually)
lamda=0.1        # weight between "pure" MaxEnt (1.0) and "pure" l1 (0.0)

# Optionnal parameters
miniSNR=16       # minimum SNR in column to do processing - 32 is optimal -
# uncertainty=1.2  # if >1 allows more room for algo when data quality is poor
# precision=1e-08  # stopping criterium

# MultiProcessing
mppool=None         # a multiprocessing pool to do a parallel processing (faster)
# if you want to use the mp capability, uncomment the following 3 lines
#import multiprocessing as mp
#NProc = 2           # here for 2 cores - adapt to your own set-up
#mppool = mp.Pool(processes=NProc)

# %%
D2.prepare_palma(finalsize, Dmin, Dmax)      # set-up diffusion axis
DD2 = D2.do_palma(nbiter=nbiter, lamda=lamda, mppool=mppool)  # launch computation
if mppool: mppool.close()  # clean-up

# %% [markdown]
# ### And display result
# Diffusion coefficients are along the vertical axis. Values are in $\mu m^2 / sec$ which is $10^{-12} m^2 / sec$.
# Typical diffusion of small molecules in a solvent are between 100 and 1000 $\mu m^2 / sec$ ($10^{-10} - 10^{-9} m^2 / sec$). Expect water around $3000 \, \mu m^2 / sec$ ($3 \, 10^{-10} m^2 / sec$)
#

# %%
DD2.axis1.currentunit = 'Diff'
DD2.projF2 = D2.row(0)
DD2.projF1 = DD2.proj(axis=1,)
I.Show2D(DD2, title="%s %s"%(FC.nmrname,d2.pulprog))

# %% [markdown]
# ## Save the data-set
# either as stand alone native SPIKE files, (there are other formats)

# %%
DD2.save('dosy.gs2')

# %% [markdown]
# *Tools in this page is under intensive development - things are going to change rapidly.*
