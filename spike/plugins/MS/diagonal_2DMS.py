#!/usr/bin/env python
# encoding: utf-8
"""computes diagonal of 2D-MS spectra

Created by DELSUC Marc-André on 2020-12-10.
"""
import numpy as np

from spike.NPKData import NPKData_plugin

def diagonal(self):
    """allows to extract the diagonal of a 2D FTMS spectrum"""
    self.check2D()
    ddiag = self.row(0)               # container
    diag =  np.zeros(self.size2)      # data
    high = int(self.axis2.mztoi(self.axis1.lowmass))  # borders in index
    low = int(self.axis2.mztoi(self.axis1.highmass))
    i = np.arange(low,high)     # convert indices to mz
    z = self.axis2.itomz(i)
    iz = np.int_(np.round(self.axis1.mztoi(z)))
    jz = np.int_(np.round(self.axis2.mztoi(z)))
    diag[jz] = self[iz,jz]      # and copy
    ddiag.set_buffer(diag)
    return ddiag

NPKData_plugin("diagonal", diagonal)