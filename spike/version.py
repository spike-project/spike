
# This file is generated automatically by dev_setup.py 
# Do not edit
ProgramName = 'SPIKE'
VersionName = 'Development version'
version = '0.99.5'
revision = '403'
rev_date = '10-04-2019'

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
