
# This file is generated automatically by dev_setup.py 
# Do not edit
from subprocess import Popen, PIPE
import re
try:
    hg = Popen(["hg", "summary"], stdout=PIPE).communicate()[0]
    revision = re.search('^parent: ([\d]*):', hg, re.MULTILINE).group(1)
except:
    revision = "--not determined--"
ProgramName = 'NPKv2'
VersionName = 'Development version - beta'
version = '0.5.1'
rev_date = '26-03-2014'

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
