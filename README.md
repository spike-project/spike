# README #

This is the beta version of the **SPIKE** program. A collaborative development for a FT-spectroscopy processing program.

# Release -
This is the version 0.8.2 dated 3 Feb 2016

 - corrected a bug in processing when running under MPI parallel 
 - added warnings in set_col() and set_row() if type do not match.
 - starting to work on the documentation

This is the version 0.8.1 dated 24 Jan 2016

- corrected a bug for Orbitrap related to offsetfreq.
- **WARNING for versions 0.8.0  0.8.1 0.8.2**
    * HDF5 files created with theses versions cannot be read with previous versions
    * HDF5 files created with previous versions cannot be read with these versions - this should be fixed later -

This is the version 0.8.0 dated 23 Jan 2016

- first clean version using the new HDF5 file set-up
 *it is however still preliminary, and many features are still missing - even documented ones*
- File now contains acquisition parameters files in the attached hdf5 sub-group
- datasets now carry store and retrieve the parmeters imported from manufacturers file in d.params
- improved FTMS calibration using 1, 2, and 3 parameters calibration : calibA calibB calibC, retrieve by Import from experimental file
- improved FTMS Hz unit, added the d.axis.offset parameter
- corrected fine details of F1 demodulation and added the parameter freq_f1demodu
- unittests extended, in particular in visu2D
- Starting with this version
  - a stable version will be maintained, downloadable as a zip file in the download page
https://bitbucket.org/delsuc/spike/downloads
  - Two developpement branches will be used, the `default` for the stable version - improved for bugs, and the `devel` branche, used for developping the new features.

Complete history in the [release_notes.md](release_notes.md) file.

## What is SPIKE ? ##

**SPIKE** is a program that allows the processing, the display and the analysis of data-sets obtained from various Fourier-Transform spectroscopies. The name stands for **S**pectrometry **P**rocessing **I**nnovative **KE**rnel.

It allows the processing of **1D** and **2D** FT spectroscopies, implementing Real, Complex and HyperComplex n-dimensionnal Fourier Transform, as well as many other functionalities.

To our knowledge, it is the first program freely available allowing the processing, display and analysis of 2D-FT-ICR (Fourier Transform Ion Cyclotron Resonance).

It is still in very active development.
Many features are missing, and many other while present, are not fully fixed.
However, considering the amount of efforts already present in this code, we decided to make it available.
We believe that even in this partial development stage, this program might prove useful for certain usages.

## Documentation
You can find a **Very preliminary** documentation [here](http://spikedoc.bitbucket.org/)

## SPIKE proposes the following features

####FT analysis of 1D data-sets

* apodisation, phasing, modulus, ...

#### Analysis of 2D data-sets

* phase or amplitude modulation
* complex or hyper-complex algebra

####Robust processing

* no limit in data-set size
* parallel processing of the heaviest processing
    * on multi-core desktop using standard python ressources
    * on large clusters, using **MPI** library

####High-Level features

* noise reduction (filtering, Linear-Prediction, Cadzow, urQRd, sane, ...)
* automatic or manual baseline correction
* 1D and 2D Peak-Picker
        
####Plugin architecture

* allow easy extension of the core program
* reduces cross dependances

####Complete spectral display using matplotlib

* zoom, available in several units (depending on the spectroscopy : seconds, Hz, ppm, m/z, etc...)
* store to png or pdf

##For the moment, SPIKE handles the following Spectroscopies

#### **NMR** 

- 1D and 2D are fully supported

#### **FT-ICR** 

- 1D and 2D are fully supported

#### **Orbitrap** 

- 1D only (!)

#### _other spectroscopies are being considered_

## Files can be imported from

* NMR:
    * Bruker Topspin
    * NMRNoteBook
    * NPK - *Gifa*
* FT-ICR:
    * Bruker Apex 
    * Bruker Solarix
* Orbitrap:
    * Thermofisher raw data
* any data in memory in a `Numpy` buffer.

# Usage

#### As a processing library
SPIKE is primary meant for being used as a library, code can as simple as :
```python
import spike           # insure spike is in your PYTHONPATH
from spike.File import Solarix 

dd = Solarix.Import_1D('FTICR-Files/ESI_pos_Ubiquitin_000006.d')  # Import create a basic SPIKE object

dd.hamming().zf(2).rfft().modulus()    # we have a simple piped processing scheme

dd.unit = "m/z"
dd.display(zoom=(500,2000))     # display the spectrum for m/z ranging from 500 to 2000

dd.pp(threshold=1E7)         # peak-pick the spectrum in this range
dd.centroid()                # compute centroids

dd.display(zoom=(856.5, 858.5))    # and zoom on the isotopic peak
dd.display_peaks(zoom=(856.5, 858.5), peak_label=True)
```

For the moment, SPIKE does not provide a installation script.
If you downloaded spike in a given directory, insure this directory is in your PYTHONPATH.
To do so, either modify the `$PYTHONPATH` environment variable, or add the following lines in the scripts that use SPIKE:
```python
import sys
sys.path.append('the_dir_where_you_put_spike_distrib')
```

#### interactive mode
SPIKE allows to process datasets interactively from an IPython prompt, and is perfectly working in `IPython Notebook` 

* Look at the examples files ( `eg_*.py` ) for examples and some documentation.
  ( * not fully up to data * )
* display is performed using the `Matplotlib` library.
* large 2D-FT-ICR are handled in batch using the `processing.py` batch program, controlled by parameter file called `*.mscf`
* The batch mode supports multiprocessing, both with MPI and natively on multi-core machines (still in-progress)
* large 2D-FT-ICR are stored in a hierarchical format, easyly displayed with an interactive program.
* data-sets are handled in the HDF5 standard file-format, which allows virtually unlimited file size ( _tested up to 200 Gb_ ).
* Version : this is 0.7 beta version

#### running stand-alone programs

processing.py and visu2D.py are two stand alone programs, written on the top of SPIKE.
 - processing.py allowing the efficient processing
   of FT-ICR 2D datasets, with no limit on the size of the final file
   Produces multi-resolution files
 - visu2D.py
   is an interactive tool for visualizing 2D FT-ICR multi-resolution files  
   
syntax :
```
python -m spike.processing param_file.mscf
```
or
```
python -m spike.visu2D param_file.mscf
```

typically, you want to add
```python
python
import sys
sys.path.append('the_dir_where_you_put_spike_distrib')
```
to the header of these scripts, and launch them from the directory which contains SPIKE the distribution.

A more complete documentation is available [here](https://spikedoc.bitbucket.org).


## How do I get SPIKE ? ##
SPIKE is written in pure Python 2.7, and relies on several external libraries.

It requires the following non-standard Python libraries :

* [Numpy](http://docs.scipy.org/doc/numpy/reference/)
* [Scipy](http://docs.scipy.org/doc/scipy/reference/)
* [Matplotlib](http://Matplotlib.org/contents.html)
* HDF5 / [Pytables](http://www.pytables.org/moin) 
* Qt / [PySide](http://qt-project.org/wiki/PySide) *optional* used by visu2D
* MPI / [mpi4py](http://www.mpi4py.scipy.org/) *optionnal* used for parallel processing of large FTICR 2D files

It has been successfully tested in the [**Enthought**](https://enthought.com/) and [**Anaconda**](http://continuum.io/) distributions.

To get it, you can simply 
 - insall the above python distributions
 - download the latest stable version here : https://bitbucket.org/delsuc/spike/downloads
 - *or* `hg clone` the devel branch and keep it up-to-date


## History ##

**SPIKE** is originated from the ** _Gifa_ ** program, developed by M-A Delsuc and others in `FORTRAN 77` since the late eighties.
_Gifa_ has known several mutations, and finally ended as a partial rewrite called **NPK**.
[NPK](http://abcis.cbs.cnrs.fr/NPK/) program is based on some of the original `FORTRAN` code, wrapped in Java and Python, which allows to control all the program possibilities from the Python level.
NPK is a pure computing kernel, with no graphical possibilities, and has been used as a kernel embedded in the commercial program NMRNoteBook, commercialized by NMRTEC.

However, NPK was showing many weaknesses, mostly due to the 32bits organization, and a poor file format. So, when a strong scientific environment became available in Python, a rewrite in pure Python was undertaken. To this initial project, called NPK-V2, many new functionalities were added, and mostly the capability to work in other spectroscopies than NMR.

At some point, we chose to fork NPK-V2 to SPIKE, and make it public.

## Citing SPIKE ##
SPIKE is not published yet, if you happen to use it successfully and wish to cite it, please refer to this site, as well as the following references :

  1.	Tramesel, D., Catherinot, V. & Delsuc, M.-A. Modeling of NMR processing, toward efficient unattended processing of NMR experiments. _J Magn Reson_ **188**, 56–67 (2007).
  2.	van Agthoven, M. A., Chiron, L., Coutouly, M.-A., Delsuc, M.-A. & Rolando, C. Two-Dimensional ECD FT-ICR Mass Spectrometry of Peptides and Glycopeptides. _Anal Chem_ **84**, 5589–5595 (2012).

## Organisation of the Code ##

The main program is `NPKData.py`, which defines NPKData object on which everything is built.

Spectroscopies are defined in the `FTICR.py` and `Orbitrap.py` code, which sub class NPKData.
It is prototyped as an NMR data-set. This set-up is temporary.

Many programs contain routines tests (in an object unittest) that also serve as an example of use.
The code goes through extensive tests daily, using the `unittest` Python library. However, many tests rely on a set of tests data-sets which is more than 1Go large, and not distributed here.


### Main programs :
a small description of the files:

- NPKData.py
   the main library, allows all processing for NMR experiments (1D, 2D and 3D)
   to be used as a library, in a stand-alone program or in IPython interactive session
- FTICR.py
   an extension of NPKData for processing FT-ICR datasets (1D and 2D)
- Orbitrap.py
   an extension of NPKData for processing Orbitrap datasets (1D)

- processing.py
   a stand alone program, written on the top of FTICR.py, allowing the efficient processing
   of FT-ICR 2D datasets, with no limit on the size of the final file
   Produces multi-resolution files
   syntax :   
   
   ```
   python -m spike.processing param_file.mscf
   ```
   
- visu2D.py
   an interactive tool for visualizing 2D FT-ICR multi-resolution files  
   
  ```
  python -m spike.visu2D param_file.mscf
  ```

### Directories
- *Algo*   
   contains algorithms to process data-sets
   (MaxEnt, Laplace, etc...) not everything active !
- *Display*    
   a small utility to choose either for regular Matplotlib display of fake no-effect display (for tests)
- *File*   
   Importers for various file format for spectrometry, as well as the HDF5 SPIKE native format.
- *plugins*   
   Tools automatically plugged in NPK kernel : display utilities, urQRd algorithm and various other tools. 
- *Miscellaneous*    
   "en vrac"
- *Visu*   
   utilities for the Visu2D program
- *util*   
   set of low-level tools used all over in the code
- *v1*    
   a library implementing a partial compatibility with the NPKV_V1 program
- *SPIKE_usage_eg*    
   example of Python programs using the various libraries available
- *example of configuration files*    
    - process_eg.mscf
    - test.mscf

- and *various utilities*

    - NPKConfigParser.py    
    	reads .mscf files
    - NPKError.py  
    	generates error msg
    - QC.py  
    	Quality Check
    - Tests.py  
    	runs all tests
    - dev_setup.py   
  		rolls a new version
    - version.py   
    	defines version number
    - \__init\__.py   
    	defines library
    - rcpylint				
    - To_Do_list.txt 
    - QC.txt  
    - Release.txt  

### Authors and Licence ###
Current Active authors for SPIKE are:

- Marc-André Delsuc     `madelsuc -at- unistra.fr`
- Lionel Chiron         `Lionel.Chiron -at- casc4de.eu`
- Petar Markov          `petar.markov -at- igbmc.fr`
- Christian Rolando     `christian.rolando -at- univ-lille1.fr`

Previous authors:
- Marie-Aude Coutouly . `Marie-Aude.Coutouly -at- nmrtec.com`

Covered code is provided under this license on an "as is" basis, without warranty of any kind, either expressed or implied, including, without limitation, warranties that the covered code is free of defects. The entire risk as to the quality and performance of the covered code is with you. Should any covered code prove defective in any respect, you (not the initial developer or any other contributor) assume the cost of any necessary servicing, repair or correction.

Downloading code and datasets from this page signifies acceptance of the hereunder License Agreement. The code distributed here is covered under the CeCILL license : http://www.cecill.info/index.en.html