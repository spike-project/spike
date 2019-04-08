*******************
First contact SPIKE
*******************


General Organization
====================
**SPIKE** is a library and a set of programs and utilities meant to process, display and analyze spectroscopic data-sets.
It is mostly oriented towards **NMR** and **FT-MS** (FT-ICR and Orbitrap) but can probably be rapidly extended to any new domain.

History of SPIKE
----------------

A Working examples
------------------
	
Plotting
--------
    ytbw

Data File
=========
    ytbw


Advance Usage
=============

SPIKE organizes internally its data in two pieces:

 - a large n-dimensional real array that contains the actual spectroscopic data
 - an Axis() object, which contains all the characteristics along the given axis (calibration, data-type, size, etc…)

This organization is the same regardless we are handling NMR Orbitrap or FTICR Data, the Axis() object is just overloaded with the corresponding set of parameters

Caveat on type
--------------
The main buffer is always a numpy array of type float64. However, some cases requires the spectroscopic data-set to be handled by complexes (for instance phase sensitive NMR of MS spectra). In this case the buffer remains a numpy array of type float64, and the attribute `.type` of the Axis() is set to 1, indicating the buffer should be considered as complex, with the real part located in even entries and imaginary parts located in odd entries (and the size forced to an even number).
All internal SPIKE routine handle this and will consider complexes a data-set with its `.type` set to 1.
For a n-Dimensional data-set, each axes has a type value, allowing hypercomplex arithmetics (see below)

Caveat on size
--------------
	ytbw

Moving data around
------------------
The  main buffer is accessible through the `.get_buffer()` and `.set_buffer()` methods.
`.get_buffer()` will actually return a plain numpy array typed relative to the axis `.type` so that for a complex data-set, this command returns a complex array, whereas for a real data-set the array is real.
So if you want to 

`.get_buffer()` has an optional argument `copy` which is `False` by default. When `False` , only an *Image* of the data is returned, with no memory location being actually copied, so it does not take memory space. 

Because of this mechanism, operations like

.. code:: python

    data.set_buffer( 2*data.get_buffer() )

is actually performed in-place in the main buffer, and does not use any additional memory space.
However, with data being typed real, doing:

    data.set_buffer( 2j*data.get_buffer() )

would create a new buffer, as a `2j` if the complex number $0 +2i$ and the complex numbers thus created take twice the room of the initial data, because of the imaginary part which has to be created.

*Remark both examples are actually a bit contrived, as SPIKE allows doing:  (see below)*

.. code:: python

    data *= 2		# to mean multiply by 2 the values of the data set


If `copy=True` however, a new memory space is created and the data are copied into it.

The main buffer is actually stored in the attribute `.buffer` of an NPKData. Thus `.buffer` is a standard numpy array, so you can actually play with it:

.. code:: python

    x = data.buffer[i]

but this is strongly discouraged, because the data-set is not checked. This syntax is used in some parts of the internal code for optimization reasons.

Data-set arithmetics
--------------------
    
    
Hypercomplex arithmetics
------------------------
	ytbw

