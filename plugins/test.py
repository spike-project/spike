"""
Test procedure for plug-ins
"""
from __future__ import print_function
from spike.NPKData import NPKData_plugin

def fake(dd, title):
    "fake method"
    dd.test_title = title
    return dd

NPKData_plugin("fake", fake)

