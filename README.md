# README #

This is the beta version of the SPIKE program. A joined development by NMRTEC and CNRS.

### What is SPIKE ? ###

* SPIKE is a program coming from a first-development oriented one named NPK-V2 that allows the processing, the display and the analysis of data-sets obtained from various Fourier-Transform spectroscopies.
It stands for Spectrometry Processing Innovative KErnel.
* For the moment, il handles the following data-sets
    * NMR - 1D and 2D are fully supported
    * FT-ICR - 1D and 2D are fully supported
    * Orbitrap - 1D only
    * _other spectroscopies are being considered_
    * Files can be imported from
        * NMR : Bruker topspin
        * FT-ICR : Bruker Apex
        * Orbitrap : Thermofisher raw data
* It allows to process datasets interactively from an ipython prompt
or interactively using the processing.py batch program (aimed towad FT-ICR for the moment)
* The batch mode supports multiprocessing, both with MPI and natively on multi-core machines (still in-progress)
* data-sets are handled in the HDF5 standard file-format, which allows virtually unlimited file size.
* Version : this is 0.5 beta version


### How do I get set up? ###
The program is in python 2.7.

look at the examples files (eg_*.py) and at configuration files (*.mscf)
they contain valuable examples and some documentation.
A more complete documentation is available at [extensive documentation]('_build/html/index.html').

SPIKE requires the following libraries :
* numpy
* scipy
* matplotlib
* Qt / PySide
* Pytables
* mpi4py
* ...

It has been successfully tested in the **Enthought** and **anaconda** [link](http://continuum.io/downloads) distributions.


### Organisation of the Code ###

The main program is NPKData.py, which defines NPKData object on which everything is built.

Spectroscopies are defined in the FTICR.py and Orbitrap.py code, which sub class NPKData
It is prototyped as an NMR data-set, but this will change.

Many programs contain routines tests (in an object unittest) that also serve as an example of use.


#### Main programs :
a small description of the files:
- NPKData.py
   the main library, allows all processing for NMR experiments (1D, 2D and 3D)
   to be used as a library, in a stand-alone program or in ipython interactive session
- FTICR.py
   an extension of NPKData for processing FT-ICR datasets (1D and 2D)
- Orbitrap.py
   an extension of NPKData for processing Orbitrap datasets (1D)

- processing.py
   a stand alone program, written on the top of FTICR.py, allowing the efficient processing
   of FT-ICR 2D datasets, with no limit on the size of the final file
   Produces multi-resolution files
   syntax :  python processing.py param_file.mscf 
- visu2D.py
   an interactive tool for visualizing 2D FT-ICR multi-resolution files
   python visu2D.py param_file.mscf

#### Directories
- Algo
   contains algorithms to process data-sets
   (MaxEnt, Laplace, etc...) not everything active !
- Display
   a small utility to choose either for regular matplotlib display of fake no-effect display (for tests)
- File
   Importers for various file format for spectrometry, as well as the HDF5 SPIKE native format.
- Miscellaneous
   "en vrac"
- Visu
   utilities for the Visu2D program
- util
   set of low-level tools used all over in the code
- v1
   a library implementing a compatibility with the NPKV_V1 program
- SPIKE_usage_eg
   example python programs using the various library available

- example of configuration files
    - process_eg.mscf
    - test.mscf

- and various utilities

    - NPKConfigParser.py	reads .mscf files
    - NPKError.py			generates error msg
    - QC.py				Quality Check
    - Tests.py				runs all tests
    - dev_setup.py			rolls a new version
    - version.py			defines version number
    - __init__.py			defines library
    - rcpylint				
    - To_Do_list.txt
    - QC.txt
    - Release.txt


### Authors and Licence ###
Authors for this code are :

- Marc-André Delsuc - CNRS
- Lionel Chiron - CNRS then NMRTEC then Casc4de
- Marie-Aude Coutouly - NMRTEC

Covered code is provided under this license on an "as is" basis, without warranty of any kind, either expressed or implied, including, without limitation, warranties that the covered code is free of defects. The entire risk as to the quality and performance of the covered code is with you. Should any covered code prove defective in any respect, you (not the initial developer or any other contributor) assume the cost of any necessary servicing, repair or correction.

Downloading code and datasets from this page signifies acceptance of the hereunder License Agreement. The code distributed here is covered under the CeCILL licence. : http://www.cecill.info/index.en.html

*[Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)*