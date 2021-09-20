#!/usr/bin/python3
# commande à lancer avant toute release
# les fichiers ici doivent être copiés dans spike/Notebooks,
# Je n'ai pas trouvé comment les distribuer via pip sinon

from pathlib import Path
import shutil

# text files
lfiles = [
"DisplayFTICR2D.ipynb",
"EasyDisplayFTICR2D.ipynb",
"Logo.ipynb",
"Proc1DNMR.ipynb",
"Proc2DNMR.ipynb",
"ProcDOSY.ipynb",
"ProcessFTICR-MS.ipynb",
"ProcessFTICR-MS-PH.ipynb",
"README.md"
]

# binary files
lfilesb = [
"Logo.png",
"Artemisinine 10mg_ml  - CDCl3.pdf"
]

here = Path(__file__).resolve().parent
nbdir = here.parent/'spike'/'Notebooks'

def cpy(source, dest, mode=''):
    with open(source,'r'+mode) as ff:
        with open(dest,'w'+mode) as FF:
            FF.write(ff.read())

for f in lfiles:
    cpy(here/f, nbdir/f)
for f in lfilesb:
    cpy(here/f, nbdir/f, mode='b')
