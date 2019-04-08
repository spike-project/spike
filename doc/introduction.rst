Introduction
=================================

This is the beta version of the **SPIKE** program. A collaborative development for a FT-spectroscopy processing program.

What is SPIKE ? 
----------------------

**SPIKE** is a program that allows the processing, the display and the analysis of data-sets obtained from various Fourier-Transform spectroscopies. The name stands for **S**\pectrometry **P**\rocessing **I**\nnovative **KE**\rnel.

It allows the processing of **1D** and **2D** FT spectroscopies, implementing Real, Complex and HyperComplex n-dimensionnal Fourier Transform, as well as many other functionalities.

To our knowledge, it is the first program freely available allowing the processing, display and analysis of 2D-FT-ICR (Fourier Transform Ion Cyclotron Resonance).

It is still in very active development.
Many features are missing, and many other while present, are not fully fixed.
However, considering the amount of efforts already present in this code, we decided to make it available.
We believe that even in this partial development stage, this program might prove useful for certain usages.

For the moment, SPIKE handles the following Spectroscopies

* **NMR** - 1D and 2D are fully supported
* **FT-ICR** - 1D and 2D are fully supported
* **Orbitrap** - 1D only (!)
* other spectroscopies are being considered
* Files can be imported from
    * NMR : Bruker Topspin / NPK (NMRNoteBook) program
    * FT-ICR : Bruker Apex
    * Orbitrap : Thermofisher raw data
    * any data in memory in a `Numpy` buffer.


SPIKE allows to process datasets interactively from an IPython prompt, and is perfectly working in `IPython Notebook` .

* Look at the examples files ( `eg_*.py` ) for examples and some documentation.
* display is performed using the `Matplotlib` library.
* large 2D-FT-ICR are handled in batch using the `processing.py` batch program, controlled by parameter file called `*.mscf`
* The batch mode supports multiprocessing, both with MPI and natively on multi-core machines (still in-progress)
* large 2D-FT-ICR are stored in a hierarchical format, easyly displayed with an interactive program.
* data-sets are handled in the HDF5 standard file-format, which allows virtually unlimited file size ( tested up to 200 Gb ).
* Version : this is 0.6 beta version

A more complete documentation is available [here](https://spikydoc.bitbucket.org).  


How do I get SPIKE ?
------------------------------
SPIKE is written in pure Python 2.7, and relies on several external libraries.

It requires the following non-standard Python libraries :

* [Numpy](http://docs.scipy.org/doc/numpy/reference/)
* [Scipy](http://docs.scipy.org/doc/scipy/reference/)
* [Matplotlib](http://Matplotlib.org/contents.html)
* Qt / [PySide](http://qt-project.org/wiki/PySide)
* HDF5 / [Pytables](http://www.pytables.org/moin) 
* MPI / [mpi4py](http://www.mpi4py.scipy.org/)

It has been successfully tested in the [**Enthought**](https://enthought.com/) and [**Anaconda**](http://continuum.io/) distributions.

History 
----------------------

**SPIKE** is originated from the **Gifa** program, developed by M-A Delsuc and others in `FORTRAN 77` since the late eighties.
_Gifa_ has known several mutations, and finally ended as a partial rewrite called **NPK**.
[NPK](http://abcis.cbs.cnrs.fr/NPK/) program is based on some of the original `FORTRAN` code, wrapped in Java and Python, which allows to control all the program possibilities from the Python level.
NPK is a pure computing kernel, with no graphical possibilities, and has been used as a kernel embedded in the commercial program NMRNoteBook, commercialized by NMRTEC.

However, NPK was showing many weaknesses, mostly due to the 32bits organization, and a poor file format. So, when a strong scientific environment became available in Python, a rewrite in pure Python was undertaken. To this initial project, called NPK-V2, many new functionalities were added, and mostly the capability to work in other spectroscopies than NMR.

At some point, we chose to fork NPK-V2 to SPIKE, and make it public.

Citing SPIKE
----------------------

SPIKE is not published yet, if you happen to use it successfully and wish to cite it, please refer to this site, as well as the following references :

1.	Tramesel, D., Catherinot, V. & Delsuc, M.-A. Modeling of NMR processing, toward efficient unattended processing of NMR experiments. _J Magn Reson_ **188**, 56–67 (2007).
2.	van Agthoven, M. A., Chiron, L., Coutouly, M.-A., Delsuc, M.-A. & Rolando, C. Two-Dimensional ECD FT-ICR Mass Spectrometry of Peptides and Glycopeptides. _Anal Chem_ **84**, 5589–5595 (2012).

Organisation of the Code
-------------------------

The main program is `NPKData.py`, which defines NPKData object on which everything is built.

Spectroscopies are defined in the `FTICR.py` and `Orbitrap.py` code, which sub class NPKData.
It is prototyped as an NMR data-set. This set-up is temporary.

Many programs contain routines tests (in an object unittest) that also serve as an example of use.
The code goes through extensive tests daily, using the `unittest` Python library. However, many tests rely on a set of tests data-sets which is more than 1Go large, and not distributed here.


Main programs :
----------------------

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
   Produces multi-resolution files syntax : 
    
   ```python -m processing param_file.mscf
   ```
   
- visu2D.py
   an interactive tool for visualizing 2D FT-ICR multi-resolution files  
   
  ```python -m visu2D param_file.mscf
  ```

Directories
----------------------

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

and *various utilities*

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
- __init__.py   
	defines library
- rcpylint				
- To_Do_list.txt 
- QC.txt  
- Release.txt  

Authors and Licence
---------------------------
Active authors for SPIKE are :

- Marc-André Delsuc  .  `madelsuc -at- unistra.fr`
- Lionel Chiron      .  `Lionel.Chiron -at- casc4de.eu`
- Marie-Aude Coutouly . `Marie-Aude.Coutouly -at- nmrtec.com`

Covered code is provided under this license on an "as is" basis, without warranty of any kind, either expressed or implied, including, without limitation, warranties that the covered code is free of defects. The entire risk as to the quality and performance of the covered code is with you. Should any covered code prove defective in any respect, you (not the initial developer or any other contributor) assume the cost of any necessary servicing, repair or correction.

Downloading code and datasets from this page signifies acceptance of the hereunder License Agreement. The code distributed here is covered under the CeCILL license : http://www.cecill.info/index.en.html