

Tutorial
=================================

Few examples of how to use SPIKE.
We begin first with simple import for both FTICR datasets and Orbitrap datasets then we show how to make more elaborated commands involving data treatment algorithms such as RECITAL and urQRd.

First, open in SPIKE directory a terminal and launch a IPython Notebook document writing:

       ipython notebook

We assume that the data are in a directory next to SPIKE directory named DATA_test.

FTICR
-----

* simple import of native dataset
* Show the FID
* Show the half truncated FID and full FID
* Doing FFT with zerofilling

simple import of native dataset
++++++++++++++++++++++++++++++++++++++++++++

from File.Apex import Import_1D

import numpy as np                 

import matplotlib.pyplot as plt     

from FTICR import FTICRData

Import from Apex

f = Import_1D("../DATA_test/angio_ms_000005.d")

Show the FID
++++++++++++++++++++++++++++++++++++++++++++

f = Import_1D("../DATA_test/angio_ms_000005.d")

f.display(label = "FID")

Show the half truncated FID and full FID
++++++++++++++++++++++++++++++++++++++++++++

f = Import_1D("../DATA_test/angio_ms_000005.d")

f.chsize(len(f.buffer)/2)

ff = f.copy()

ff.buffer = ff.buffer[:len(f.buffer)/2]/2

f.display(label = "FID")

f.display(label = "FID cut", new_fig = False)

Doing FFT with zerofilling
++++++++++++++++++++++++++++++++++++++++++++

Classical FFT with apodisation and zerofilling.

FFT with zerofilling, processing cutting the pipes.

Here instead of writing a single long command with pipelines, the command is cut in many chunks.
This can be used for performing intermediate operations not present in SPIKE.


Orbitrap
---------------

Some examples on how to use NPKv2.

* simple import of native dataset.
* simple FID handling, processing and display
* FFT with zerofilling

We begin first with simple import then we show how to make more elaborated commands involving data treatment algorithms such as RECITAL and urQRd.

simple import of native dataset
++++++++++++++++++++++++++++++++++++++++++++

o = Import_1D("C:/Users/Egor/NPK_V2/DATA_test/ubiquitin_5_scan_res_30000_1.dat")

Show FID
++++++++++++++++++++++++++++++++++++++++++++

o.display(label = "FID")

FFT with zerofilling
++++++++++++++++++++++++++++++++++++++++++++

print o

o.apod_sin(maxi = 0.5).chsize(o.buffer.size*2).rfft().modulus()

o.units = 'm/z'

o.display(label = "zerofill x2")

FFT with zerofilling, processing cutting the pipes.
+++++++++++++++++++++++++++++++++++++++++++++++++++++

o = Import_1D(filename)

o.units = 'm/z'

o.apod_sin(maxi = 0.5)

o.chsize(o.buffer.size*4)

o.rfft()

o.modulus().display(label = "zerofill x4")


urQRd
-------------
is a preprocessing technique used for reducing the noise.
The parameter $k$ given to urQRd is related to the number of expected lines in the spectrum. 
It should be chosen 2 to 3 times larger than this expected number.
Be carefull than the processing time **and** the memory footprint are both proportionnal to this value.

data.units = 'm/z'

data.urqrd(k = 300).rfft().modulus().display(label = "urQRd, rank = 300")

Additional tricks
++++++++++++++++++++++++++++++++++++++++++++

IPython shortcuts.

there are *many* shortcuts and tricks in IPython, read the doc !

a couple of them are really helpfull for MS processing

* you can execute a cell by hitting `shift-return`
* you can get the documentation of any function by adding a ? at the end of its name, eg `o.rfft?`
* you can get all possible values by hitting the `<TAB>` key. Try for instance typing `o. <TAB>`

SPIKE arcanes
++++++++++++++++++++++++++++++++++++++++++++

If needed, you can directly manipulate the numeric data held into the SPIKE dataset:

* the `.get_buffer()` method returns the underlying `numpy` array.
* The `.set_buffer()` method sets it, data can be real or complex.
* Do `.adapt_size()` afterwards if you changed the number of points.

It is also possible to use this sheet as a simple calculator, can be handy some time, for instance for checking charge state.
