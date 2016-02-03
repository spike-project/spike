#!/usr/bin/env python 
# encoding: utf-8

"""
dev_setup.py

To be called any time a new version is rolled out !

Created by Marc-Andre' on 2010-07-20.
"""

from __future__ import print_function

ProgramName = "SPIKE"
VersionName = "Development version - beta"
VersionInfo = ["0", "8", "2"]   # Major - Minor - Micro

# Major.minor.micro (int) + name eg SPIKE_2_19_5
# N.M.L
#
# standard way of doing :
# backward incompatible change N -> N+1,  M=0, L=0
# new features    M+1, L=0
# bug fixes L+1

# Release Notes in md syntax !
release_notes="""
# Release Notes

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
        revision = re.search(r'^parent: ([\d]*):', hg, re.MULTILINE).group(1)
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
    print '''
    ========================
          %s
    ========================
    Version     : %s
    Date        : %s
    Revision Id : %s
    ========================'''%(ProgramName, version, rev_date, revision)
report()
""")
    #print 'SPIKE version', version, 'date',date "
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
    print '''
    ========================
          %s
    ========================
    Version     : %s
    Date        : %s
    Revision Id : %s
    ========================'''%(ProgramName, version, rev_date, revision)
report()
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
    generate_notes("release_notes.md")
#    plier()
    # then tests
    test = False
    if test:
        import Tests
        Tests.CLEAN = True  # tells test suite to remove ALL temporary files - will produce them again
        import Display.testplot as testplot
        testplot.PLOT = False   # switches off the display for automatic tests
        Tests.do_Test()
    
