
# This file is generated automatically by dev_setup.py 
# Do not edit
ProgramName = 'SPIKE'
VersionName = 'Development version - beta'
version = '0.8.2'
revision = '264'
rev_date = '03-02-2016'

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
