![SPIKE](spike/Interactive/Logo.png)


## SPIKE features
- FT analysis of 1D data-sets
    - apodisation, phasing, modulus, ...
- Analysis of 2D data-sets
    - phase or amplitude modulation
    - complex or hyper-complex algebra
- Robust processing
    - no limit in data-set size
    - parallel processing of the heaviest processing
        - on multi-core desktop using standard python ressources
        - on large clusters, using **MPI** library
- High-Level features
    - noise reduction (filtering, Linear-Prediction, Cadzow, urQRd, sane, ...)
    - automatic or manual baseline correction
    - NUS data processing
    - 1D and 2D Peak-Picker
- Plugin architecture
    - allow easy extension of the core program
    - reduces cross dependences
- Complete spectral display using matplotlib or bokeh transparently
    - zoom, available in several units (depending on the spectroscopy : seconds, Hz, ppm, m/z, etc...)
    - store to png or pdf
- interaction with the Jupyter Notebook environment

## Handles the following Spectroscopies / Spectrometries
- **NMR** 
    - 1D and 2D are fully supported
    - DOSY
    - no nD (n>2) yet
    - relaxation data in progress
- **FT-ICR-MS** 
    - 1D and 2DFT are fully supported
    - LC-MS in progress
- **Orbitrap-MS** 
    - 1D only (!)
- _other spectroscopies are being considered_

## Files can be imported from
* NMR:
    * Bruker Topspin
    * NPK - *Gifa*
    * RS2D SpinIt
    * NMRNoteBook 
* FT-ICR:
    * Bruker Apex 
    * Bruker Solarix
* Orbitrap:
    * Thermofisher raw data
* Other
    * csv and txt files
    * any data in memory in a `Numpy` buffer.

# Usage

## As a processing library
SPIKE is primary meant for being used as a library, code can as simple as :

```python
from spike.File import Solarix 

dd = Solarix.Import_1D('FTICR-Files/ESI_pos_Ubiquitin_000006.d')  # Import create a basic SPIKE object

dd.hamming().zf(2).rfft().modulus()    # we have a simple piped processing scheme
  # here doing apodisation - zerofilling (doubling the size) - FT and modulus.
  # calibration is imported from Bruker - advanced calibration is available

dd.unit = "m/z"
dd.display(zoom=(500,2000))     # display the spectrum for m/z ranging from 500 to 2000

dd.pp(threshold=1E7)         # peak-pick the spectrum in this range
dd.centroid()                # compute centroids

dd.display(zoom=(856.5, 858.5))    # and zoom on the isotopic peak
dd.display_peaks(zoom=(856.5, 858.5), peak_label=True)
```

## interactive mode
SPIKE allows to process datasets interactively from an jupyter (IPython) prompt, and is perfectly working in `jupyter notebook` or even `jupyter lab`

* Look at the examples files ( `eg_*.py` and `*.ipynb` ) for examples and some documentation.
  ( * not fully up to data * )
* display is performed using the `Matplotlib` library.
* large 2D-FT-ICR are handled in batch using the `processing.py` batch program, controlled by parameter file called `*.mscf`
* The batch mode supports multiprocessing, both with MPI and natively on multi-core machines (still in-progress)
* large 2D-FT-ICR are stored in a hierarchical format, easyly displayed with an interactive program.
* data-sets are handled in the HDF5 standard file-format, which allows virtually unlimited file size ( _tested up to 500 Gb_ ).



# 2D FTICR-MS

Processing of 2D FTICR-MS experiment can be done in two manners.

#### directly with a python script, as a regular data-set.
However this approach is limited by the size of computer memory, which has to be able to accept the size of the final data. 
Considering the current possibilities it is the case only for small to medium experiment, and with a suffisantly large memory.

For instance a $1024 \times 256k$ experiment processed in modulus, so that the final spectrum is the same size will require 8GB of memory simply to store the computation.
So don't try it if you have less than 12-16 GB of central memory as some room is needed for the program and the system.

*Be careful here*, as if you hit the memory limit, the computer usually become *so slooow* that usually there is no other solution than rebooting.

#### using the processing.py program
It is a stand alone program, written on the top of SPIKE, which allows the efficient processing
of FT-ICR 2D datasets on files, with no other limit on the size of the final file tha the size of the disk.
It also directly produces multi-resolution files well adapted to the Display Notebook.

It works with a parameter file `param_file.mscf`

to use it, do:
```
python -m spike.processing param_file.mscf
```


### Main program files :
a small description of the files:

- NPKData.py
  the main library, allows all processing for any kind of experiments (1D, 2D and 3D)
  to be used as a library, in a stand-alone program or in IPython interactive session
- NMR.py
   The NPKData library adapted to NMR processing
- FTICR.py
   an extension of NPKData for processing FT-ICR datasets (1D and 2D)
- Orbitrap.py
   an extension of NPKData for processing Orbitrap datasets (1D)


```python

```
