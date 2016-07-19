#!/usr/bin/env python 
# encoding: utf-8

"""This plugin implement the set of Fourier Transform used for NMR

and some other related utilities

MAD September 2015

"""

from __future__ import print_function
import unittest

from spike import NPKError
from spike.NPKData import NPKData_plugin, as_float, as_cpx

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
    return data
NPKData_plugin("ft_phase_modu", ft_phase_modu)

#---------------------------------------------------------------------------
def ft_n_p(data, axis ="F1"):
    """F1-Fourier transform for N+P (echo/antiecho) 2D"""
    data.conv_n_p().ft_sh()
    return data
NPKData_plugin("ft_n_p", ft_n_p)
#-------------------------------------------------------------------------------
def bruker_corr(self):
    """
    applies a correction on the spectrum for the time offset in the FID.
    time offset is stored in the axis property zerotime
    """
    delay = self.axes(self.dim).zerotime
    self.phase(0, -360.0*delay, axis=0) # apply phase correction
    return self
NPKData_plugin("bruker_corr", bruker_corr)
#-------------------------------------------------------------------------------
def bruker_proc_phase(self):
    """
    applies a correction on the spectrum for the time offset in the FID and from proc file parameters.
    """
    ph1 = -float(self.params['proc']['$PHC1'])
    ph0 = -float(self.params['proc']['$PHC0'])+ph1/2
    zero = -360*self.axes(self.dim).zerotime
    self.phase(ph0, ph1+zero) #Performs the phase correction from proc file
    return self
NPKData_plugin("bruker_proc_phase", bruker_proc_phase)
#-------------------------------------------------------------------------------  
def conv_n_p(self):
    """
    realises the n+p to SH conversion
    """
    self.check2D()
    for k in xrange(0,self.size1,2):
        a = self.row(k)
        b = self.row(k+1)
#        self.set_row(k,a+b )    # put a+b
#        self.set_row(k+1,(a-b).phase(90,0) )    # put i( a-b )
        # this is about 60% faster !
        self.set_row(k,a.copy().add(b) )    # put sum
        a.buffer += -b.buffer             # compute diff
        a.buffer = as_float( 1j*as_cpx(a.buffer) )  # dephase by 90°
        self.set_row(k+1,a)                 # and put back
    self.axis1.itype = 1
    return self
NPKData_plugin("conv_n_p", conv_n_p)

#########################################################
class Bruker_NMR_FT(unittest.TestCase):
    def test_1D(self):
        '''Test mostly syntax'''
        from ..NPKData import NPKData
        d1 = NPKData(dim=1)
        d1.axis1.itype = 0
        d1.ft_seq()
        d1.axis1.itype = 1
        d1.ft_sim().bruker_corr()
        #d1.bruker_proc_phase()
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

        d1.axis1.itype = 1
        d1.ft_n_p()
        
        d1.axis2.itype = 1
        d1.axis1.itype = 0
        d1.ft_phase_modu()
