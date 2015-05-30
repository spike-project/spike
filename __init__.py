"""
The Spike Package

"""
#from __future__ import absolute_import
import os, sys

# every thing is in version.py
from version import version as __version__
from version import VersionName as __version_info__
from version import ProgramName as __program_name__
#from version import revision as __revision__
from version import rev_date as __date__

__author__ = "Marc A. Delsuc <delsuc@igbmc.fr>, Marie-Aude Coutouly, Lionel Chiron"
global SPIKE_version
SPIKE_version = __version__


#### Header to set-up the whole SPIKE environment
import spike.NPKData

### plugins to the spike.NPKData class.
# simply put a xxx.py in the plugins folder - and define the interface as described in the doc
from spike.plugins import load

load()
