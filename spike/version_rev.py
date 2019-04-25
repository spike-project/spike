
# This file is generated automatically by dev_setup.py 
# Do not edit
from subprocess import Popen, PIPE
import re
try:
    hg = Popen(["hg", "summary"], stdout=PIPE).communicate()[0]
    revision = re.search('^parent: ([\d]*):', hg, re.MULTILINE).group(1)
except:
    revision = "--not determined--"
ProgramName = 'SPIKE'
VersionName = 'Development version'
version = '0.99.8'
rev_date = '25-04-2019'

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
