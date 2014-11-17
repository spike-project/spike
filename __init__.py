"""
The NPK Package

"""
# every thing is in version.py
from version import version as __version__
from version import VersionName as __version_info__
from version import ProgramName as __program_name__
from version import revision as __revision__
from version import rev_date as __date__

__author__ = "Marc A. Delsuc <delsuc@igbmc.fr>, Marie-Aude Coutouly, Lionel Chiron"
global NPK_version
NPK_version = __version__


#### Header to set-up the whole NPK environment

import NPKData
# print dir()