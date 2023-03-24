# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
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
# This program supersedes the old `Visu2D` program, developped in `Qt`, which could no longer be maintained.\ables\
# In addition to a regular anaconda python 3 distribution, it requires `pytables` and `ipympl` - to be installed with conda.
#
#
# *This is a work inprogress - additional utilities should come soon !*

# %% [markdown]
# ## To use it, 
# - in the "Cell" pulldown above, select "Run all", the program will start.
# - select the file you want look at, and Load, it will show-up as a full width 2D image.
# *ignore eventual warnings about missing attributes*
#     - the F2/horizontal axis is the high resolution, direct axis. You find fragments along this line
#     - the F1/vertical axis is the low resolution, indirect axis. You have parents along this axis. 
# - with the zoom tool (the square below the spectrum) you can select the region you want to display, you can also dial it in on the top box.
# - with the `scale` slider, you can select the display level.
# *Think to this as an archipelago viewed on the map; the slider changes the sea level, raising the scale raises the floor (or lower the sea level)*
# - The data have a hierarchical multiresolution structure. Zooming does not change the resolution, to force it, either click on the `Redraw` button, or change the scale.
# The smaller the zoom box, the better the resolution.

# %%
### Initialization of the environment
### the following cell should be run only once *(but no harm if you run it twice)* .

from IPython.display import display, HTML, Markdown, Image
display(Markdown('## STARTING Environment...'))
# %matplotlib widget
import spike
import spike.Interactive.INTER as I
import spike.Interactive.FTICR_INTER_v2 as IF2
display(Markdown('## ... program is ready'))
from importlib import reload  # the two following lines are debugging help
reload(IF2)                   # and can be removed safely when in production
I.hidecode(initial='hide',message="")
I.hidedoc(message="")
I.initCSS()
I.Logo()
#I.hidecode(initial='hide', message="")
ms = IF2.MS2Dscene(root='/DATA/DATA/FT-ICR/')

# %% [markdown]
# ### to come
# - calibration
# - peak detection
# - superimposition
# - extraction of arbitrary 1D 
# - locate/remove artifacts due to harmonics

# %% [markdown]
# You can also call use the multi-resolution tool alone:

# %%
filename = "/DATA/My/Dataset/filename.msh5"
try:
    M = IF2.MR_interact(filename,  show=True, figsize=(12,12), Debug=False)
except:
    print("ERROR - check filename")

# %%
