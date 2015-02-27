
# This file is generated automatically by dev_setup.py 
# Do not edit

ProgramName = 'SPIKE'
VersionName = 'Development version - beta'
version = '0.6.0'
rev_date = '03-01-2015'

def report():
    "prints version name when NPK starts"
    print '''
    ========================
          %s
    ========================
    Version     : %s
    Date        : %s
 
    ========================'''%(ProgramName, version, rev_date)
report()
