#!/usr/bin/env python 
# encoding: utf-8

"""
This plugins implement the set of Fourier Transform used for NMR

MAD September 2015

NOT FULLY TESTED !!!
"""

from __future__ import print_function
import unittest

from spike import NPKError
from spike.NPKData import NPKData_plugin

########################################################################
def check_real(data, axis):
    if data.axes(axis).itype != 0 : 
        raise NPKError("Axis %d should be real"%axis, data=data)
def check_cpx(data, axis):
    if data.axes(axis).itype != 1 : 
        raise NPKError("Axis %d should be complex"%axis, data=data)
#---------------------------------------------------------------------------
def ft_sim(data):
    """performs the fourier transform of a data-set acquired on a Bruker in
    simultaneous mode
    Processing is performed only along the F2 (F3) axis if in 2D (3D)

    (Bruker QSIM mode)"""
    todo = data.dim
    data.revf().fft(axis=todo)
    return data
NPKData_plugin("ft_sim", ft_sim)

#---------------------------------------------------------------------------
def ft_seq(data):
    """performs the fourier transform of a data-set acquired on a Bruker in simultaneous mode
    Processing is performed only along the F2 (F3) axis if in 2D (3D)

    (Bruker QSIM mode)"""
    todo = data.dim
    data.revf().rfft(axis=todo)
    return data
NPKData_plugin("ft_seq", ft_seq)

#---------------------------------------------------------------------------
def ft_tppi(data, axis = "F1" ):
    """TPPI F1 Fourier transform"""
    if data.dim == 1:
        raise NPKError("Not implemented in 1D", data=data)
    todo = data.test_axis(axis)
    data.revf(axis=todo).rfft(axis=todo)
    return data
NPKData_plugin("ft_tppi", ft_tppi)


#---------------------------------------------------------------------------
def ft_sh(data, axis = "F1" ):
    """ States-Haberkorn F1 Fourier transform"""
    if data.dim == 1:
        raise NPKError("Not implemented in 1D", data=data)
    todo = data.test_axis(axis)
    data.revf(axis=todo).fft(axis=todo)
    return data
NPKData_plugin("ft_sh", ft_sh)

#---------------------------------------------------------------------------
def ft_sh_tppi(data, axis = "F1" ):
    """ States-Haberkorn / TPPI F1 Fourier Transform """
    if data.dim == 1:
        raise NPKError("Not implemented in 1D", data=data)
    todo = data.test_axis(axis)
    data.fft(axis=todo)
    return data
NPKData_plugin("ft_sh_tppi", ft_sh_tppi)


#---------------------------------------------------------------------------
def ft_phase_modu(data, axis = "F1" ):
    """F1-Fourier transform for phase-modulated 2D"""
    if data.dim != 2 :
        raise NPKError("implemented only in 2D",data=data)
    data.flip().revf("F1").fft("F1").flop().reverse("F1")
NPKData_plugin("ft_phase_modu", ft_phase_modu)

#---------------------------------------------------------------------------
def ft_n_p(data, axis ="F1"):
    """F1-Fourier transform for N+P (echo/antiecho) 2D"""
    data.conv_n_p().ft_sh()

NPKData_plugin("ft_n_p", ft_n_p)

class Bruker_NMR_FT(unittest.TestCase):
    def test_1D(self):
        '''Test mostly syntax'''
        from ..NPKData import NPKData
        d1 = NPKData(dim=1)
        d1.axis1.itype = 0
        d1.ft_seq()
        d1.axis1.itype = 1
        d1.ft_sim()
    def test_2D(self):
        '''Test mostly syntax'''
        from ..NPKData import NPKData
        d1 = NPKData(dim=2)
        d1.axis2.itype = 0
        d1.ft_seq()
        d1.axis2.itype = 1
        d1.ft_sim()
        d1.axis1.itype = 0
        d1.ft_tppi()
        d1.axis1.itype = 1
        d1.ft_sh()
        d1.axis1.itype = 1
        d1.ft_sh_tppi()
