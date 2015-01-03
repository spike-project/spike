# encoding: utf-8
"""
dev_setup.py

To be called any time a new version is rolled out !

Created by Marc-André on 2010-07-20.
"""


ProgramName = "SPIKE"
VersionName = "Development version - beta"
VersionInfo = ["0", "6", "0"]   # Major - Minor - Micro

# Major.minor.micro (int) + name eg NPK_2_19_5
# N.M.L
#
# standard way of doing :
# backward incompatible change N -> N+1,  M=0, L=0
# new features    M+1, L=0
# bug fixes L+1

"""
0.6.0 - dec 2014
    - Fork to SPIKE
    - Large improvements of the display program, renamed visu2D
    - Corrected a bug in the hypercomplex modulus, resulting in splitting in 2D-FT-ICR
    - many improvements everywhere
0.5.1 - 26 mar 2014
    - processing2.py renamed to processing.py   with added features
        - urQRd
        - FISTA
    - source code reorganized by folders
0.5.0 - 24 mar 2014
    - starting new devl effort
    - published ! version of urQRd
0.4.1 - 27 Sep 2012
    - final (?) version of urQRd
    - added data arithmetic
    - many other optimisation
0.4.0 - 20 apr 2012 -
new version processing2.py (temporary name)  this one
    - processing is performed in steps, F2 from infile to interfile (intermediate file) and F1 from interfile to outfile
    - steps are optionnal, F2 or F1 can be performed alone - allowing denoising on the interfile
    - processing is faster and mpi enabled - speed-up are better for very large files
    - has a better way of computing the smaller spectra - done by downsampling - faster and nicer
    - vignette is now 1024x1024 - can be changed using SIZEMIN in config file
0.3.11 - 29 mar 2012 - Small tools have been added to modify configuration files and to mix processing.py and ipython visualisation
0.3.10 - 22 jan  2012 - processing is now (hopefully) bug free and RAPID !
0.3.9 - 18 jan  2012 - fticrvisu.py, processing working, getting all parameters correctly from FTICRData and Apex
0.3.8 - 13 jan  2012 - fticrvisu.py, processing working, corrected after Marie came
0.3.7 - 12 dec  2011 - correction of Gifa file bug, bug in Apex for narrow band data-sets, changes in msh5 file format
0.3.6 -  3 Oct  2011 - added  HDF5 file format (.msh5), multiresolution files, configuration files (.mscf), fticrvisu
0.3.5 -  5 Sept 2011 - added  cadzow in MPI / savitsky-golay / HDF5 still in progress
0.3.4 - 26 July 2011 - added  autotests / savehdf5 first version
0.3.3 - 12 July 2011 - first reliable/taged  FTICR version
"""

# Version control used - only one should be true
UsingSVN = False # subversion
UsingHG = True # mercurial
if UsingSVN and UsingHG:
    raise Exception("Please define only one flag UsingHG or UsingHG !")

####################### End of configuration ###################
from subprocess import Popen, PIPE
import re
from datetime import date

def generate_version():
    """generates version string, revision id and data"""
    #version string
    version = ".".join(VersionInfo)
    # revision
    if UsingSVN:
        svn = Popen(["svn", "info"], stdout=PIPE).communicate()[0]
        revision = re.search('^Revision: (\d*)', svn, re.MULTILINE).group(1)
    elif UsingHG:
        hg = Popen(["hg", "summary"], stdout=PIPE).communicate()[0]
        revision = re.search('^parent: ([\d]*):', hg, re.MULTILINE).group(1)
    else:
        revision = "Undefined"
    # date   dd-mm-yyyy
    today = date.today().strftime("%d-%m-%Y")
    return (version, revision, today)

def generate_file(fname):
    """
    write version to the file "name", usually "version.py", used later on
    then version.py is imported at NPK initialization.
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
    "prints version name when NPK starts"
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
    #print 'NPK version', version, 'date',date "
    f.close()

def do(arg):
    "print and execute"
    print " ".join(arg)
    retcode = Popen(arg)

def plier():
    "fabrique le zip"
    name = "NPKV2_beta_" + ( "_".join(VersionInfo) )
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
    print """
=====================
Generating version.py
=====================
Version      :  %s
Revision Id  :  %s
=====================
"""%(version, revision)
    generate_file("version.py")
#    plier()
    # then tests
    test = False
    if test:
        import Tests
        Tests.CLEAN = True  # tels test suite to remove ALL temporary files - will produce them again
        import Display.testplot as testplot
        testplot.PLOT = False   # switches off the display for automatic tests
        Tests.do_Test()
    
