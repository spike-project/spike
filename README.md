![SPIKE](spike/Interactive/Logo.png)
# 1/ What is SPIKE ? #

**SPIKE** a collaborative development for a FT-spectroscopy processing library.

*This is the version 0.99.32 - March 2023*

**WARNING (31th March 2023) - It seems that the various Jupyter Notebooks which use the ipympl library have a problem if you download a fresh anaconda environment (based on python 3.11)**<br>
*upgrading ipympl to the 0.9.3 version seems to solve the difficulty:*

    conda config --env --add channels conda-forge
    conda install ipympl=0.9.3

**SPIKE**  allows the processing, the display and the analysis of data-sets obtained from various Fourier-Transform spectroscopies. The name stands for **S**pectrometry **P**rocessing **I**nnovative **KE**rnel.

It allows the processing of **1D** and **2D** FT spectroscopies, mostly **NMR** and **FTICR-MS**, but also **Orbitrap**, and other to come.

It implements all needed tools for Data analysis (Fourier Transform, Baseline correction, noise reduction, peak-picking, etc...)

To our knowledge, it is the first program freely available allowing the processing, display and analysis of 2D-FT-ICR (Fourier Transform Ion Cyclotron Resonance), as well as Orbitrap time domain data processing.

It is written in python (currently 3.7 and in python 2.7 up to version 0.99.10).

It is mostly a library for processing, analysis and display of spectral data-sets.
The natural interface is either a python program or a `jupyter notebook` as an interactive front-end.


### SpI *(Spike Interactive)*
An additional program, called **SP**ike **I**nteractive (or **SpI** ) is also developed, it consists in a set of stand alone jupyter Notebooks which allow to use SPIKE process and analyse NMR and MS datasets without any programming knowledge.

It provides simple interactive environment which implement basic analysis can be found in the [SpI](Notebooks) directory.
A part of SpI is also contained in the [`Interactive`](spike/Interactive) directory in the source tree.

It is still in very active development, but also the basis for our research.

Many features are missing, and many other while present, are not fully fixed.
However, considering the amount of efforts already present in this code, we decided to make it available.
We believe that even in this current development stage, this program might prove useful for certain usages.

# 2/ Documentation
Basic presentation of functionalities can be found [Here](Presentation.md)

A complete documentation, with the full API, can be found [Here](http://softwares.casc4de.eu/spike/spikedoc/).
We try to keep it up to date, but it may lag a little behind the current program state.

Most of the documentation is in the code itself, available directly while programming (*try* `function_name?` *in* `jupyter notebook`).

Another source of documentation are [Notebooks](Notebooks) directory.
There are two type of Notebooks, fully developed ones, with interactive tools, which can be used as is,
and example Notebooks, meant for pedagogic purposes.
This is an ongoing effort, and always short of what we'd like to do.

### Release Notes
Can be found [Here](release_notes.md)


# 3/ How to get SPIKE ? ##

### simple installation
You might have python installed on your system, however, SPIKE relies on a whole set of scientific libraries which may dont be installed (see **Dependencies** below).
We advise to install a scientific python environment such as [Anaconda](https://www.anaconda.com/) or [Enthough](https://www.enthought.com/).
Then install SPIKE with the following command in the terminal:

    pip install spike-py

This will install the latest version on your machine. *(be sure to use the pip command which came with the scientific python environment)*.

If you already have SPIKE installed and want to upgrade to the last version do:

    pip install -U spike-py

### Installing additional tools
Several additional tools may be of interest.

**a/ Interactive jupyter notebooks from SpI**
can be found in the [Notebooks](Notebooks) directory or on [github](https://github.com/spike-project/spike/tree/master/Notebooks). 
You will find there fully developed interactive notebooks allowing data analysis without any knowledge in python programming, as well as examples to starting writing small pyhton programs for specific needs.

First the tool `ipympl` needs to be installed (it may not be installed along the scientific environment)
To install it, simply do:

    pip install ipympl

Then copy the related NotBooks to a directory by typing `python -m spike.installSpI` in a terminal in this directory.


**b/ additional plugins** are available in the [plugins/special](spike/plugins/special) directory in the source tree. [check on github](https://github.com/spike-project/spike/tree/master/spike/plugins/specials).

If you want to use them, create a directory called `Spike` in your HOME directory, and a directory called `plugins` in this one (hence `$HOME/Spike/plugins`), and put then there, they will be directly activated. 
You can also develop your own plugins and put them here.

Note that some plugins or extension require additional libraries (`ipympl`, `MPI`, `bokeh`, `mayavi`, ...)

### Licenses
The Spike program is distributed under the [`CeCILL 2.1 FREE SOFTWARE LICENSE AGREEMENT`](http://www.cecill.info)
Which is the GPL license adapted to the French law.

It includes the following additional programs/algorithms:

- [urQRd](http://urqrd.igbmc.fr/) L.Chiron and M-A Delsuc (denoising of large harmonic signals) under CeCiLL licence
- [PALMA](https://github.com/delsuc/PALMA) A.Cherni & M-A.Delsuc (Processing of DOSY based on ILT) under CeCiLL licence
- [progressbar](https://github.com/niltonvolpato/python-progressbar) N.Volpato - under LGPL licence

### source and dependencies
The SPIKE source is available at https://github.com/spike-project/spike

SPIKE is written in pure Python, and relies on several external libraries:

- matplotlib
- numpy
- scipy
- tables
- pandas

It is compatible and fully tested with python 3.7 and 3.9

### developping for SPIKE
check [Here](DevelopmentGuide.md)



# 4/ Citing SPIKE
If you happen to use SPIKE successfully for your research, please cite it, and refer to this site, as well as the following possible references :
- **Main Reference**: first publication of the program itself - *rejected from Anal. Chem. with no real critics except that Reviewer 1 said "too much NMR", Reviewer 2 said "too much MS", !! so I decided to let it on ArXiV*)
    1.    Chiron L., Coutouly M-A., Starck J-P., Rolando C., Delsuc M-A. SPIKE a Processing Software dedicated to Fourier Spectroscopies   https://arxiv.org/abs/1608.06777 (2016)

Other references are also related
- The renewal of Gifa:
    2. Delsuc M-A. "Gifa V.4: A complete package for NMR data set processing" (2020) https://doi.org/10.5281/zenodo.3904595
- presentation of the automation possibilities in NMR
    3. Margueritte, L., Markov, P., Chiron, L., Starck, J.-P., Vonthron S&eacute;n&eacute;cheau, C., Bourjot, M., & Delsuc, M.-A. (2018). "Automatic differential analysis of NMR experiments in complex samples." Magn. Reson. Chem., 80(5), 1387. http://doi.org/10.1002/mrc.4683
- first version of the 2D FT-ICR-MS processing
    4.    van Agthoven, M. A., Chiron, L., Coutouly, M.-A., Delsuc, M.-A. & Rolando, C. "Two-Dimensional ECD FT-ICR" 
Mass Spectrometry of Peptides and Glycopeptides." _Anal Chem_ **84**, 5589-95 (2012).
- first version of the python set-up on which the current SPIKE is loosely based
    5.    Tramesel, D., Catherinot, V. & Delsuc, M.-A. "Modeling of NMR processing, toward efficient unattended processing of NMR experiments. _J Magn Reson_ **188**, 56-67 (2007).

ref 1) is a general purpose reference, the other ones are more specific.


# 5/ Contents of this repository

- [`spike/`](spike) - the python code of the program itself
- [`doc/`](doc) - a directory containing various files used to generate the documentation
- [`Notebooks/`](Notebooks) - the repository of the **SpI** notebooks
- miscelaneous:
    - `QC.py` `QC.txt` a small utility to check the quality of the code using pylint. 



# 6/ Origin of the program

**SPIKE** is originated from the _Gifa_ program, developed by M-A Delsuc and others in `FORTRAN 77` since the late eighties, it still shares similarities, such as command names, or internal data organisation.

_Gifa_ was abandoned around 2005, but it actually experienced a revival in 2020 thanks to the " *ten years challenge* " organised by the [ReScience C Journal](https://rescience.github.io/read/). You can find an installable and runnable version here: https://github.com/delsuc/Gifa

_Gifa_ has known several mutations, and finally ended as a partial rewrite called **NPK**.
The [NPK](http://abcis.cbs.cnrs.fr/NPK/) program is based on some of the original `FORTRAN` code, wrapped in Java and Python, which allows to control all the program possibilities from the Python level.
NPK is purely a computing kernel, with no graphical possibilities, and has been used as a kernel embedded in the commercial program NMRNoteBook, commercialized by NMRTEC.

However, NPK was showing many weaknesses, mostly due to the 32bits organization, and a poor file format. So, when a strong scientific environment became available in Python, a rewrite in pure Python was undertaken. To this initial project, called NPK-V2, many new functionalities were added, and mostly the capability to work in other spectroscopies than NMR.

At some point in 2014, we chose to fork NPK-V2 to SPIKE, and make it public.


# 6/ Authors and Licence
Current Active authors for SPIKE are:

- Marc-Andr&eacute; Delsuc     `madelsuc -at- unistra.fr`
- Laura Duciel          `laura.duciel -at- casc4de.eu`

Previous authors:

- Christian Rolando     `christian.rolando -at- univ-lille1.fr`
- Lionel Chiron         `Lionel.Chiron -at- casc4de.eu`
- Petar Markov          `petar.markov -at- igbmc.fr`
- Marie-Aude Coutouly . `Marie-Aude.COUTOULY - at- datastorm.fr`

Covered code is provided under this license on an "as is" basis, without warranty of any kind, either expressed or implied, including, without limitation, warranties that the covered code is free of defects. The entire risk as to the quality and performance of the covered code is with you. Should any covered code prove defective in any respect, you (not the initial developer or any other contributor) assume the cost of any necessary servicing, repair or correction.

Downloading code and datasets from this page signifies acceptance of the hereunder License Agreement. The code distributed here is covered under the CeCILL license : http://www.cecill.info/index.en.html

