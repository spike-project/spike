#!/usr/bin/env python 
 # encoding: utf-8
"""
The NPK Package

"""
__version__ = "1.03"
__author__ = "Marc A. Delsuc <delsuc@igbmc.fr>, marie-Aude Coutouly, Vincent Catherinot <v.catherinot@nmrtec.com>, Dominique Tramesel <d.tramesel@nmrtec.com>"
__date__ = "june 2010"
NPK_version = __version__


#### Header to set-up the whole NPK environment

# Import NPK java environment
import sys

#__all__ = ["Generic", "Param", "Process1D", "Process2D", "Process3D", "GenericDosy", "ProcessDosy",  "Nucleus", "Bruker", "NPK_version" ] # "Peak","Varian",
from . import Kore
#from . import Generic
#from . import Param
#from . import Process1D
#from . import Process2D
#from . import Process3D
#from . import GenericDosy
#from . import ProcessDosy
Kore.compatibility(locals())

