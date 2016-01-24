#!/usr/bin/env python 
# encoding: utf-8

"""
untitled.py

Created by Marc-André on 2010-07-20.
Copyright (c) 2010 IGBMC. All rights reserved.
"""

from __future__ import print_function
        
class NPKError(Exception):
    """ implements NPK generic exception
        adds the named argument data, which can be used to describe the NPKData involved
    """
    def __init__(self, msg="", data=None):
#        super(NPKError, self).__init__()
        self.msg = msg
        self.data = data
    def __str__(self):
        st = self.msg
        if self.data is not None:
            st = """
%s

Data-set involved :
===================
%s
"""%(st,self.data.report())
        return st

