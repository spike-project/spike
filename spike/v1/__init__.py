#!/usr/bin/env python 
 # encoding: utf-8
"""
The NPK Package


"""
__version__ = "1.2"
# version 1.2 is a reactivation of NPK v1 to use in spike
__author__ = "Marc A. Delsuc <delsuc@igbmc.fr>"
# previous participants marie-Aude Coutouly, Vincent Catherinot <v.catherinot@nmrtec.com>, Dominique Tramesel <d.tramesel@nmrtec.com>" have left long ago
__date__ = "Jan 2020"
NPK_version = __version__


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

