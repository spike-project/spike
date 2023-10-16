<!-- DO NOT MODIFY this file, it is generated automatically from dev_setup ! -->
# SPIKE Relase Notes
*if you need the latest version, get Spike from github*
*then to make available to you system,* do
    > cd spike
    > python setup.py devel

#### 0.99.33 - October 2023 - not released yet
- Interactive:
    - NMR
        - improved 1D integration
        - improved the Show2D() tool, and added local projections
    	- in Interactive, added AQ in the parameter list
- MS: corrected a bug in pk2pandas() in 2D FTICR which was crashing the action.
- NMR: corrected a logic bug in BrukerNMR importer which triggered when importing data with several expno
- NMR: modified the File.BrukerNMR.Export_dif to correctly export fid
- corrected a subtil bug which always removed the right most point in display() and extract()
- corrected a bug in noise computation which was crashing baseline correction sometimes (and elsewhere)
- many small corrections due to warnings and deprecations in particular in numpy

#### 0.99.32 - March 2023
- Interactive:
    - NMR and MS: example interactive notebooks have been extended, cleaned up and fully checked 
        - better 2D display
        - added documentation 
    - NMR:
        - ProcDOSY.ipynb rewritten
        - ReadSMX.ipynb to import TopSpin processed file into Spike
        - SNR evaluation in Proc1D.ipynb
        - added choice of projection in 2D graphics
        - baseline correction is now fully functional
- updated the plugin mechanism - the previous system was incompatible with the coming python 3.12
    - requires python >= 3.5
    - should not change anything to your set-up
- added  Algo.Linepredic.baselinerollrem() to remove baseline roll using Linear Prediction
- extended real2cpx() to 2D operations (was only 1D)
- NMR: improved importing processed DOSY (PALMA plugin)
    - corrected cyclic import in PALMA
- NMR: in BrukerSMX.py corrected for a warning signaling an upcoming obsolence

#### 0.99.31 - Nov 2022
- NMR:improvements in spike Interactive
    - on spectral superposition
    - improved import/export of peak and integral lists, as csv files
- added BrukerNMR.Export_fid (1D)
- improvements in automatic baseline correction
- improved apodisation functions
    - MS: added the kilgour() apodisation, from Kilgour et al DOI:10.1002/rcm.7190 handy for phased FTICR-MS data
    - all: added the `maxi` parameter for hamming(), hanning(), and kaiser() which becomes a versatile and generalized apodisation
      *( was already available for apod_sin() and apod_sq_sin() )*
    - NMR: normalized gaussenh() (which was not !)
- adapted to python 3.10 - cleaned the code (better grade in QC!)
- **a problem with multiprocessing crashing (core dump!) on my system, still investigating**
- corrected a bug in baseline correction, and improved interface
- NMR: Corrected a bug when importing DOSY from TopSpin

#### 0.99.30 - April 2022
- spike Interactive NMR
    - modified SHOW1D to handle the move tool
    - modified selection to right-click (in phaser, peakpicker, baseline)
    - display of some acquisition parameters
    - in 1D, corrected errors in managing peaks and calibration
    - better handling of baseline corrections, phasing and better graphics
    - dated figure dump
    - corrected displaying integrals
    - improvement in object hierarchy
    - added yratio to SHOW1D, and allow a better zoom experience 
- General
    - adapted to matplotlib 3.5 / python 3.9 
    - improved BrukerNMR.import1D_proc to handle different binary options
    - added yratio to npk.display()
    - rowc() and colc() method get rows and columns in 2D using current unit instead of index (I don't know why this was not done earlier...)
    - added the itoix() ixtoi() ctoix() ixtoc() to access datapoint in complex coordinates
    - added a method Peaks.Peaklist1D.pkadd() to add a peak list to an existing one
    - added a method Peaks.Peaklist1D.pkaggreg(distance) to aggregate close peaks  
    - marker "None" in peak display
    - moved config into $(HOME)/.config/Spike
    - corrected a bug when reading old Gifafile where diffusion is not defined in the header

#### 0.99.29 - Sept 2021
- spike Interactive FTICR - not tested yet -
- spike Interactive NMR
    - **Phaser2D does not work yet**
    - improved the quality of 2D display (smoother interaction, better scale selection)
    - corrected a bug introduced in previous version for displaying 2D in notebooks
- General: hilbert() ihilbert() and real2cpx() basic transforms
    *(were already there, but not isolated as now, should have been done long ago)*
    - added a "None" mode in displaying peaks

#### 0.99.28 - Aug 2021
- spike Interactive NMR
    - interactive phasing is now much smoother and faster
    - added superimposition of spectra in Proc1DNMR - other will be adapted
    - adapted and optimized display in the notebook and in the core code
- General:
    - added a disp_artist field in NPKData which contains the last matplotlib artist generated by display()

#### 0.99.27 - July 2021
- spike interactive in FTICR - many improvements
    - isotopic pattern simulation using imported `isotope` module
    (and fine isotopic as well - thanks to the excellent prgm neutronstar !)
    - processing of narrowband acquisition is correct now
    - better calibration - (still in progress)
    - superimposition of spectra should work now

#### 0.99.26 - July 2021
- FT-ICR:
    - added support for NarrowBand FTICR acqusition - first steps, need to be completed
- General
    - new command zero_dsp() and rem_dsp() to zero or remove the FID header created by the DSP.
    - bug corrections in the 2D peak centroid and display

#### 0.99.25 - June 2021
*0.99.24 was short lived, because additional bugs, introduced in.22 were found and corrected*

- correction for plugins not updated correctly in the pip distribution since 0.99.22,  was making errors in interactive tools
- corrected a bug when fitting peaks on a zoomed region
- NMR: corrections in `Proc1DNMR` notebook in integration tool

#### 0.99.24 - June 2021
- added a `maxdist` flag in the peak aggregator in plugins.Peaks  (default to 10*distance)
- NMR: corrected a bug introduced in 0.99.22 in the interactive phasing when moving the pivot

#### 0.99.23 - May 2021
- MS: first interactive notebook for phasing of FTICR-MS spectra - still a bit rough
- NMR: slight improvements in apmin() (automatic phaser) better search algorithm
- slight improvements in bcorr() (baseline correction), when working on complex spectra
- The logo shows up when loading interactive notebooks - and also a nicer set-up  
- new doc is still in development and missing, sorry (phasing is keeping me busy) 

#### 0.99.22 - April 2021 - not released on pypi -
- change in the plugin set up, there can now be distributed in several directories (in that order)
    - `(distrib dir)/spike/plugins`  - basic plugins allways loaded
    - `$HOME/spike/plugins`  - plugins specific to the user allways loaded
    - `(distrib dir)/spike/plugins/NMR`  - plugins specific to NMR, loaded with `import spike.NMR`
    - `(distrib dir)/spike/plugins/MS`  - plugins specific to MS, loaded with `import spike.FTMS`
    - and as before, plugins with a name starting with a `_` are not loaded
- MS: new `PhaseMS` plugin, which implements quadratic phase correction and permits to phase FTICR spectra
- a .tm() apodisation (trapeze) which emulates D.Kilgour apodisation for phase sensitive FTICR-MS
- rewrote and reorganized  README and documentation

#### 0.99.21 - Feb 2021
- MS: 
    - plugin which implements .diagonal() for computing the diagonal of 2D FTICR spectra 
    - changed the logic to generate downsampled 2D FTICR spectra - smaller files, smaller vignettes
    - correction when reading Apex MS dataset for pulse frequency limits - thanks to Maria van Agthoven
    - correction ThermoFisher/Orbitrap import code - thanks to Will Kew
- small correction in urQRd  - thanks to Will Kew
- small corrections when opening files
- many small bugs corrected

#### 0.99.20 - Nov 2020
- MS: corrected a bad bug which corrupted F1 calibration when loading a 2D-FTICRMS experiment
    this bug was introduced in the 0.99.14 release but was not detected at that time

#### 0.99.19 - May 2020 
- NMR: corrected a bug in BrukerNMR importer... 

#### 0.99.17 - April 2020 
- MS: corrected a bug in BrukerMS importer...

#### 0.99.16 - April 2020 
- MS a new Apex0 bruker importer - to access old datasets, with the "NMR" setup (acqus pdata ...)
- MS a global BrukerMS importer - Import1D - Import2D
- a few corrected bugs

#### 0.99.15 - March 2020 
This release introduces a major modification in the organisation of the `NPKData` object -
which is the central object on which everything is organized.

*Previously* `NPKData.NPKData` was the standard class, which created NMR object, and other classes
(such as `FTICR.FTICRData` or `Orbitrap.OrbiData`) inherited from this class and had to overloaded a few things. 
`NPKData` held also all `Axis` definition, both generic and for NMR.

*Now*,

- `NPKData._NPKData` is a generic object - agnostic about the spectroscopy
- `NPKData` holds also the definitions of generic Axes ( `Axis` but also `TimeAxis` and `LaplaceAxis`) 
- `NMR.NMRData` is the new class for NMR data-sets, `NMR` also contains the definitions for all NMR related Axes.
- `FTICR.FTICRData` and `Orbitrap.OrbiData` now inherit from the `_NPKData` class (through `FTMS`).

In consequence, to create an NMR dataset from scratch, now do:
```
        NMR.NMRData(...)
```
where you were using `NPKData.NPKData` previously

and do 
```
        NPKData._NPKData(...)
```
to create an empty dataset not associated to any spectroscopy

*This should have been done long ago - but I'm so lazy...*

**Other modifications**
- This set-up allows to better adapt compound experiments (LC-NMR LC-MS ...)
- jupytext extension was added to jupyter Notebooks
    - means that a python copy is maintained - only this copy is version controlled
- still improvements in notebooks
- 
- NMR: improvement in the SpinIt importer
- MS: added a Bo attribute
- added the NbMaxPeaks flag in Peak display
- added the self.kind attribute in the Axis class - easier to use than self.NMR !
- small bugs corrections
- adding a complex value to a complex datasets was wrong in complex mode.
- tests in python 2.7 are abandoned - but very few python 3 features are really used...
- REMARK, it was always mentionned that version 1.0 would be rolled out
  when interactive functions in notebooks would be really usefull.
  It will be the next big release probably !

#### 0.99.14 - October 2019 - not released on pypi -
- NMR: lots of improvements in Proc1DNMR notebooks
- improvements in NoteBook for mouse interactivity (click and scroll)
    - requires the additional ipympl module

#### 0.99.13 - October 2019 - not released on pypi -
- MS: added the EasyDisplayFTICR2D for non programers !
- improvements in NoteBook interactivity
- added smoothing in spline baselinecorrection

#### 0.99.12 - September 2019
The 0.99.11 had a bug in the display of 1D NMR experiment - the 0.99.12 corrects it.

#### 0.99.11 - September 2019
- NMR: added a Notebook for processing of DOSY
- many improvement in the interactive Notebooks, and in the interactive library (still work to do though)
- added autpoints computation for spline baseline correction
- corrected axis placement in spectral display (you should not have inverted axis anymore)
- corrected a bug when computing projections

#### 0.99.10 - August 2019
- changed calibration in FTICR-MS - now should better correspond to Bruker, both for linear and quadratic
    - be carefull, the definitions are slightly modified, this should be taken into account when reading files,
    however you should verify the calibration stored into previous files
- added Proc2DNMR Notebook - preliminary!
- continued to improve other Notebooks
- added skewness and kurtosis in bucket lists (optional)
- improved Test suite (should mostly work under Windows now)
- corrected peak-picker so that  width is FWMH after centroid
- added an option in Peaks.pk2pandas to output or not the uncertainties

#### 0.99.9 - June 2019
- improved many aspect of the interactive Notebooks
- improved Proc1DNMR Notebook
    - added peak-picker
    - added integration
    - added bucketing
- Integrate: plugin for 1D NMR data integration
- added peak lists export to pandas
- added limit to the number of peak to be displayed on screen (default 1000)

#### 0.99.8 - April 2019
- corrects a BIG BUG which hampers the import of 1D NMR data-sets,

#### 0.99.7 - April 2019
**Please do not use,** see above
- first version of Interactive Notebooks:
    - ProcessFTICR-MS
    - DisplayFTICR2D
    - Process1DNMR
- The new BrukerMS importer tries Solarix importer and falls back to Apex importer if it fails
- improved Interactive tools
- improved FTICR importers to accept `my_expt.d/fid` as well as `my_expt.d`
- improved peak-picker behavior
- improved error messages in FTICR importers
- cured a bug in Apex.Importxx for a special xml format

#### 0.99.6 - April 2019
- extended and improved tests - finalized installation through PYPI
- support for distribution via pip - you can now do `pip instal spike_py` and spike installed globally on your system.
- still struggling with correct calibration routines for FTICR ! - proceed carefully ! -
- phase() speeded-up by a factor 20 !
- added a autothresh scaling to peakpicking (catching peaks "autothresh" times above the noise level - default is 3)
- slight improvement of peak list reporting (additional key words: format="report" and format="full")
- a bug in extract of complex 1D data-sets was corrected
- added the figure keyword to peaklist display

#### 0.99.3 - March 2019
- Development of Interactive tools, to be used within Jupyter - *should be extended in further releases* -
    - a tool for displaying multiresolution  2D FT-ICR-MS data-sets
    - simple interface in Jupyter for 1D NMR
    - (that part is not tested in python 2)
- added the setup.py prgm, SPIKE is now a regular installable program - still working on it ! -
- `scale="auto"` in 2D display, choose a level `autoscalethresh` (default is 3) times above the noise floor.
- added `gaussenh` apodisation plugin for one-command gaussian enhancement.
- improved display of FTMS spectra
- improved `findnoiselevel` and `findnoiselevel_2D`
- modified `absmax` in NPKData - now a property
- added phase parameters to NMRAxis: `.P0` and `.P1` 
- removed the old Visu2D program - use the jupyter notebook rather !

#### 0.99.2 - January 2019
- added number of local peaks in bucketing
- improved Bruker importer and added support for NEO/TopSpin 4.0 files
- improved the `.set_unit()` method
- improved importing DOSY processed with TopSpin
- corrected a bug for min value in bucketing
- changed pylint/QC defaults -> new values (and corrected a bad bug)
- cleaned the code a little

#### 0.99.1 - November 2018
- added the sane algorithm
- added the pg-sane algorithm
- added the `NPKData.set_unit(unit)` method for pipelining
- added the `NPKData.load_sampling(axis)` method for pipelining
- improved spinit support
- corrected a few bugs
    - `NPKData.save_csv()` now works in python 3
    - `NPKData.copy()` is now more robust

#### 0.99 - April 2018 - temp release branch
We have been developping a lot this last year, and published quite a few results.
The program is now quite stable in most of its features.
Additions and improvements were added to the repository in the `devel` branch, however we neglected updating the more official `default` branch.
This release is an effort to bring everything into normal mode, and hopefully, preparing a 1.0 version !

New in 0.99:

- SPIKE is now fully compatible with python 2 AND python 3
- added the SANE noise denoising algorithm and plugin.
    - an improvement to urQRd
    - more faithfull to small signal intensity
    - slightly different optimum parameters (optimal rank slightly smaller, less iterations needed)
- added the handling of NUS 2D FTICR acquisition
- added the PALMA DOSY processing algo and plugin (NMR).
- added a Linear Prediction plugin
- added the first trial for a m/z calibration plugin (MS)
- added import from SpinIt (NMR)
- added a primitive set of interactive tools to be used in Jupyter notebooks ( `INTER.py` )
- added the possibility to pass a complete dictionary to matplotlib in the .display() method
- added the .center() method for NPKData
- added a plugin implementing a subset of Topspin commands: xf1, xf2, xfb.  (NMR)
- added an line fitter, still very exploratory, only 1D Lorentzian for the moment
- added more controls on plots (new_fig and mpldic arguments of `.display()` )
- added a Spinit importer (preliminary) (NMR)
- added a compress mode in Solarix importer (MS)
- added new automatic tests
- improved and extended the Bucketing plugin, with extended features
- improved the baseline correction code
- improved import/export to Topspin/Bruker NMR files
- improved automatic phaser `.apmin()` (NMR)
- improved the plugin mechanism - with added documentation
- corrected the extract() method which was broken
- corrected a bug when importing Topspin/Bruker NMR datasets, where $NC was not used. (NMR)
- corrected a bug and improved 3 parameters FT-ICR calibration (MS)
- corrected the extract function for NPKData
- corrected a bug with contour plots and matplotlib version > 1.5.0
- modified (improved?) plugin loading code, with additional plugin documentation
- modified the way None values are stored into hdf5 files
- modified `.extract()` code to work in current axis unit
- modified `.mean()` to return complex value is axis is complex
- improved python 3 compatibility. It is not finished yet, but most of the program is python 2/python 3 independent, some parts are still missing, 

- known bugs
    - `NPKData.extract()` method not fully tested
    - `NPKData.save_csv()` is buggy in python 3

#### 0.9 - 8 sept 2016
*never reached the normal distribution - doc partly redundant with 0.8.3*

- added a baseline correction plugins, already quite developed, with 3 different methods
- added an automatic phasing plugin, `.apmin()` still exploratory (NMR)
- added a wavelet filtering plugin (requires the PyWavelet library)
- added a 3D zoom plugin (requires the Mayavi library) 
- added export to Topspin/Bruker files, and added import of processed Topspin files (NMR)
- added the upgrade of files from previous version
- added the `d.axis?.cpxsize` : the size of an axis expressed in spectroscopic points (real of complex)
 different from `d.axis?.size` which is the size of an axis expressed in data points so
   - `d.axis?.cpxsize == d.axis?.size`     is axis is real
   - `d.axis?.cpxsize == d.axis?.size/2`   is axis is complex
- improved the Peak-Picker (mostly the output capabilities)
- improved processing.py for nicer spectra, and possibly faster processing (MS)
- improved visu2D.py, for a greater stability and improved selection syntax
- corrected a bug in `.conv_n_p()` (NMR)
- and many small bugs as well

#### 0.8.3 - April 2016
- ALL spectro.
    - added a new `cpxsize` property, associated to axes and dataset, which counts complex and real entries
    - added: display and peak display now accept a color and markersize arguments
    - improved plugins, plugins with a filename starting with _ do not load
    - improved: automatic baseline correction algorithms have been improved ( `Algo/BC.py` )
    - `finnoiselevel()` set of functions has been rewritten ( `util/signal_tools.py` )
    - standard test now includes testing for `multiprocessing` - *DOES NOT WORK ON ALL DISTRIBUTION* if it is your case,
      set `use_multiprocessing = False` in test.mscf
- NMR
    - added: BrukerNMR now imports TopSpin processed dataset (1r, 2rr)
    - improved: and corrected Laplace axes - for a new DOSY module to come...
    - corrected: conv_n_p() was wrong and has been corrected
    - corrected: `gm_apod()` was wrong and has been corrected_
    - corrected: an error in GifaFile access under Windows
- MS
   - processing.py (2D FTMS) now includes parallel processing in F2 (helping in certain cases)
   - and gives sharper lineshape thanks to kaiser() apodisation
   - files from the previous program version (0.7.x) can now be upgraded and read. just do
      ```    python -m spike.File.HDF5File update your_file.msh5  ```
   - improved `.report()` for FTMS datasets

#### 0.8.2 - 2 Feb 2016
 - corrected a bug in processing when running under MPI parallel 
 - added warning in set_col() and set_row() if type do not match.
 - starting to work on the documentation

#### 0.8.1 - 24 Jan 2016
 - corrected a bug for Orbitrap related to offsetfreq.

#### 0.8.0 - 23 Jan 2016
- first clean version using the new HDF5 file set-up  **WARNING** 
    - HDF5 files created with this version cannot be read with previous versions
    - HDF5 files created with previous versions cannot be read with this version - this should be fixed later -
      File now contains acquisition parameters files in the attached hdf5 sub-group
- datasets now carry store and retrieve the parmeters imported from manufacturers file in d.params
- improved FTMS calibration using 1, 2, and 3 parameters calibration : calibA calibB calibC, retrieve by Import from experimental file
- improved FTMS Hz unit, added the d.axis.offsetfreq parameter
- corrected fine details of F1 demodulation and added the parameter freq_f1demodu
- unittests extended, in particular in visu2D
- Starting with this version
    - a stable version will be maintained, downloadable as a zip file in the download page
  https://bitbucket.org/delsuc/spike/downloads
    - Two developpement branches will be used, the `default` for the stable version - improved for bugs, and the `devel` branche, used for developping the new features.

#### 0.7.1 - 5 Jan 2016
- greatly improved internal compression of msh5 files and speed of processing.py
- many small corrections and bug fixes.

#### 0.7.0 - November 2015
- a plugin mechanism has been created which allows to add very simply new features to the program
    - most new features are implemented through this mechanism
- the organisation of the spectral axes has been complete modified, with the introduction of a Unit class
    - each axis holds its own series of possible units (called .units)
    - and the current unit used for display and selection
    - many commands now have a zoom= ketword that works in the current unit
    - additionnaly, there are itoc and ctoi unit converters
- Thanks to this, NMR data-sets are now correctly handled, DOSY are still in progress and should come soon
    - addtionnaly a plugin for Bruker NMR processing is now implemented
- A complete 1D and 2D peak-picker is now implemented, with many controls and features
- New baseline correction algo have been implemented
- the sane algorithm, which is an evolution from urQrd has been separated from urQRd, so both algo can now be used independently
- Tests have been reorganized and improved
- Importers have been extended - parameters are now brought back to the user
- many others

#### 0.6.4 - march 2015
- added Bruker NMR import
- clean-up of the module, still going on
- Tests improved

#### 0.6.3 - march 2015
- first installeable release

#### 0.6.0 - dec 2014
- Fork to SPIKE
- Large improvements of the display program, renamed visu2D
- Corrected a bug in the hypercomplex modulus, resulting in splitting in 2D-FT-ICR
- many improvements everywhere

#### 0.5.1 - 26 mar 2014
- processing2.py renamed to processing.py   with added features
    - urQRd
- source code reorganized by folders

#### 0.5.0 - 24 mar 2014
- starting new devl effort
- published ! version of urQRd

#### 0.4.1 - 27 Sep 2012
- final (?) version of urQRd
- added data arithmetic
- many other optimisation

#### 0.4.0 - 20 apr 2012 -
new version processing2.py (temporary name)  this one

- processing is performed in steps, F2 from infile to interfile (intermediate file) and F1 from interfile to outfile
- steps are optionnal, F2 or F1 can be performed alone - allowing denoising on the interfile
- processing is faster and mpi enabled - speed-up are better for very large files
- has a better way of computing the smaller spectra - done by downsampling - faster and nicer
- vignette is now 1024x1024 - can be changed using SIZEMIN in config file

#### 0.3.11 - 29 mar 2012
Small tools have been added to modify configuration files and to mix processing.py and ipython visualisation

#### 0.3.10 - 22 jan  2012
processing is now (hopefully) bug free and RAPID !

#### 0.3.9 - 18 jan  2012
fticrvisu.py, processing working, getting all parameters correctly from FTICRData and Apex

#### 0.3.8 - 13 jan  2012
fticrvisu.py, processing working, corrected after Marie came

#### 0.3.7 - 12 dec  2011
correction of Gifa file bug, bug in Apex for narrow band data-sets, changes in msh5 file format

#### 0.3.6
3 Oct  2011 - added  HDF5 file format (.msh5), multiresolution files, configuration files (.mscf), fticrvisu

#### 0.3.5
5 Sept 2011 - added  cadzow in MPI / savitsky-golay / HDF5 still in progress

#### 0.3.4
26 July 2011 - added  autotests / savehdf5 first version

#### 0.3.3
12 July 2011 - first reliable/taged  FTICR version
