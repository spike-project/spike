#!/usr/bin/env python 
# encoding: utf-8

"""
dev_setup.py

To be called any time a new version is rolled out !

Created by Marc-Andre' on 2010-07-20.
"""

from __future__ import print_function
from codecs import encode, decode

ProgramName = "SPIKE"
VersionName = "Development version"
VersionInfo = ["0", "99", "8"]   # Major - Minor - Micro

# Major.minor.micro (int) + name eg SPIKE_2_19_5
# N.M.L
#
# standard way of doing :
# backward incompatible change N -> N+1,  M=0, L=0
# new features    M+1, L=0
# bug fixes L+1

# Release Notes in md syntax !
release_notes="""
# SPIKE Relase Notes
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
"""

# Version control used - only one should be true
UsingSVN = False # subversion
UsingHG = True # mercurial
# projet moved from svn to hg - git not defined yet
if UsingSVN and UsingHG:
    raise Exception("Please define only one flag UsingHG or UsingHG !")

####################### End of configuration ###################
from subprocess import Popen, PIPE

import re
from datetime import date

def generate_notes(fname):
    "write the release notes file"
    with open(fname,'w') as F:
        F.write("<!-- DO NOT MODIFY this file, it is generated automatically from dev_setup ! -->")
        F.write(release_notes)
def generate_version():
    """generates version string, revision id and data"""
    #version string
    version = ".".join(VersionInfo)
    # revision
    if UsingSVN:
        svn = Popen(["svn", "info"], stdout=PIPE).communicate()[0]
        revision = re.search(r'^Revision: (\d*)', svn, re.MULTILINE).group(1)
    elif UsingHG:
        hg = Popen(["hg", "summary"], stdout=PIPE).communicate()[0]
        revision = re.search(r'^parent: ([\d]*):', decode(hg), re.MULTILINE).group(1)
    else:
        revision = "Undefined"
    # date   dd-mm-yyyy
    today = date.today().strftime("%d-%m-%Y")
    return (version, revision, today)

def generate_file(fname):
    """
    write version to the file "name", usually "version.py", used later on
    then version.py is imported at SPIKE initialization.
    No revision version is included
    """
    f = open(fname,"w")
    f.write(r"""
# This file is generated automatically by dev_setup.py 
# Do not edit
""")
    #version = ".".join(VersionInfo)
    #today = date.today().strftime("%d-%b-%Y")
    
    (version, revision, today) = generate_version()
    
    f.write("ProgramName = '%s'\n"%ProgramName)
    f.write("VersionName = '%s'\n"%VersionName)
    f.write("version = '%s'\n"%version)
    f.write("revision = '%s'\n"%revision)
    f.write("rev_date = '%s'\n"%today)
    f.write(r"""
def report():
    "prints version name when SPIKE starts"
    return '''
    ========================
          {0}
    ========================
    Version     : {1}
    Date        : {2}
    Revision Id : {3}
    ========================'''.format(ProgramName, version, rev_date, revision)
print(report())
""")
    f.close()
def generate_file_rev(fname):
    """
    write version to the file "name", usually "version.py", used later on
    then version.py is imported at SPIKE initialization.
        it assumes hg is used on client side ! - no svn version.
    """
    f = open(fname,"w")
    f.write(r"""
# This file is generated automatically by dev_setup.py 
# Do not edit
from subprocess import Popen, PIPE
import re
try:
    hg = Popen(["hg", "summary"], stdout=PIPE).communicate()[0]
    revision = re.search('^parent: ([\d]*):', hg, re.MULTILINE).group(1)
except:
    revision = "--not determined--"
""")
    (version, revision, today) = generate_version()
    f.write("ProgramName = '%s'\n"%ProgramName)
    f.write("VersionName = '%s'\n"%VersionName)
    f.write("version = '%s'\n"%version)
    f.write("rev_date = '%s'\n"%today)
    f.write(r"""
def report():
    "prints version name when SPIKE starts"
    return '''
    ========================
          {0}
    ========================
    Version     : {1}
    Date        : {2}
    Revision Id : {3}
    ========================'''.format(ProgramName, version, rev_date, revision)
print(report())
""")
    #print 'SPIKE version', version, 'date',date "
    f.close()

def do(arg):
    "print and execute"
    print(" ".join(arg))
    retcode = Popen(arg)

def plier():
    "fabrique le zip"
    name = "SPIKE_beta_" + ( "_".join(VersionInfo) )
    dir = "../"+name
    zip = dir+".zip"
    do( ["rm", "-r", dir] )
    do( ["rm", zip] )
    do( ["mkdir", dir] )
    do( ["hg", "clone", ".", dir] )
    do( ["zip", "-r", zip, dir] )
    
    
if __name__ == '__main__':
    # generate version file
    (version, revision, today) = generate_version()
    print("""
=====================
Generating version.py
=====================
Version      :  %s
Revision Id  :  %s
=====================
"""%(version, revision))
    generate_file("version.py")
    generate_file_rev("version_rev.py")
    generate_notes("../release_notes.md")
#    plier()
    # then tests
    test = False
    if test:
        import Tests
        Tests.CLEAN = True  # tells test suite to remove ALL temporary files - will produce them again
        import Display.testplot as testplot
        testplot.PLOT = False   # switches off the display for automatic tests
        Tests.do_Test()
    
