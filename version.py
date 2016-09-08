
# This file is generated automatically by dev_setup.py 
# Do not edit
ProgramName = 'SPIKE'
VersionName = 'Development version - beta'
version = '0.9.0'
revision = '331'
rev_date = '08-09-2016'

def report():
    "prints version name when SPIKE starts"
    print( '''
    ========================
          {0}
    ========================
    Version     : {1}
    Date        : {2}
    Revision Id : {3}
    ========================'''.format(ProgramName, version, rev_date, revision))
report()
