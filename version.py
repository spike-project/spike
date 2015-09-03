
# This file is generated automatically by dev_setup.py 
# Do not edit

from __future__ import print_function

ProgramName = 'SPIKE'
VersionName = 'Development version - beta'
version = '0.6.4'
rev_date = '12-03-2015'

def report():
    "prints version name when SPIKE starts"
    print('''
    ========================
          %s
    ========================
    Version     : %s
    Date        : %s
    ========================'''%(ProgramName, version, rev_date))
report()
