import os
from pathlib import Path
# This little program copies the current Notebooks to the local directory

print('--------------------------------------------------------')
nbdir = Path(__file__).resolve().parent/'Notebooks'
lfiles = [i.name for i in nbdir.glob('*.ipynb')]+['README.md']

print("""
Spike comes with a set of interactive tools, and notebooks, collectively called SpikeInteractive or SpI

This program will copy the corresponding files from {0} to the current location:

{1}

files with the same name in the current location will be overriden.

""".format( nbdir, '\n'.join(lfiles)) )

s = 'Y'
s = input("Ok to proceed (Y/n)")
if s not in ('',None,'Y','y'):
     print('Aborted')
     exit(0)

print()
for f in lfiles:
    print('copying',  f, end='...')
    with open(nbdir/f,'r') as ff:
        with open(f,'w') as FF:
            FF.write(ff.read())
    print('Done')

print("""
These programs are jupyter notebooks programs
so the jupyter tool, with the ipympl add-on should be installed
To lanch them, either type
 
 > jupyter notebook

in a terminal, or use the anaconda launcher

This will start a browser with this list,  just click on a file to used it.

The files can be modified, duplicated and copied anywhere 
""")

