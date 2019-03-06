#!/usr/bin/env python 
# encoding: utf-8

"""
FTMS.py

This file defines generic classes for FT-MS Spectroscopy (FT-ICR and Orbitrap)

not meant to be used directly, but rather called from either Orbitrap.py or FTICR.py

Created by Marc-André on 2014-08
Copyright (c) 2014 IGBMC. All rights reserved.
"""

from __future__ import print_function
import math
import unittest
import numpy as np
from . import NPKData
from .NPKError import NPKError


HighestMass = 10000000.0  # the highest m/z of interest ever

class FTMSAxis(NPKData.Axis):
    """
    hold information for one FT-MS axis
    used internally
    """
    def __init__(self, itype=0, currentunit="points", size=1024, specwidth=1E6,  offsetfreq=0.0, left_point = 0.0, highmass=10000.0, calibA=0.0, calibB=0.0, calibC=0.0 ):
#  lowmass_disp=222.0, highmass_disp=222222.0,):
        """
        all parameters from Axis, plus

        specwidth   highest frequency,
        offsetfreq      carrier frequency in heterodyn or lowest frequency if acquisition does notcontains 0.0,

        calibA, calibB, calibC : calibration constant, allowing 1 2 or 3 parameters calibration.
            set to zero if unused
            correspond to Bruker parameter ML1 ML2 ML3 for FTICR
            correspond to Thermo parameter  'Source Coeff1', 'Source Coeff2', 'Source Coeff3' for Orbitrap
        highmass    highest physical m/z of interest
        left_point  coordinates of first data point; usually 0.0 after Fourier Transform; may be different after extraction        
        currentunit default unit used for display and zoom,
            possible values for unit are "points" "m/z"

        """
        super(FTMSAxis, self).__init__(size=size, itype=itype)
        self.specwidth = specwidth
        self.offsetfreq = offsetfreq

        self.calibA = calibA
        self.calibB = calibB
        self.calibC = calibC
        self.left_point = left_point
        self.highmass = highmass      

        self.units["m/z"] = NPKData.Unit(name="m/z", converter=self.itomz, bconverter=self.mztoi)
        self.units["Hz"] = NPKData.Unit(name="Hz", converter=self.itoh, bconverter=self.htoi)
        self.units["sec"]= NPKData.Unit(name="sec", converter=self.itos, bconverter=self.stoi)  # for transients
        self.currentunit = currentunit
        self.attributes = ["specwidth", "highmass", "offsetfreq", "left_point", "calibA", "calibB", "calibC"] + self.attributes

    #-------------------------------------------------------------------------------
    @property
    def borders(self):
        """the (min, max) available windows, used typically for display"""
        if self.currentunit == 'm/z':
            return (int(self.mztoi(self.highmass)), self.size-1)
        else:
            return (0, self.size-1)
    @property
    def lowmass(self):
        """highest mass of interest - defined by the Nyquist frequency limit"""
        return self.itomz(self.size-1)
    #-------------------------------------------------------------------------------
    def extract(self, zoom):
        """
        redefines the axis parameters so that the new axe is extracted for the points [start:end]
        zoom is defined in current axis unit
        """
        #   before index       leftpoint----------------------------size
        #   before hz     off---------------------------------------off+SW
        #
        #   after                    start---------------end
        #   index                    leftpoint'----------size'
        #   Hz             off---------------------------off+SW'
        start, end = self.getslice(zoom)
        new_specwidth = self.itoh(end-1) - self.offsetfreq   # check itoh() for the -1
        self.specwidth = new_specwidth
        self.size = end-start
        self.left_point += start
        return (start, end)

    # The 2 htomz() and mztoh() are used to build all other transfoms
    # they are different in Orbitrap and FTICR and should be redefined there
    def htomz(self, value):
        """
        return m/z (mz) from Hz value
        """
        from warning import warn
        warn("htomz must be overridden by subclasser")
        return value
    def mztoh(self, value):
        """
        return Hz value from m/z (mz)
        """
        from warning import warn
        warn("mztoh must be overridden by subclasser")
        return value
    # the following are generic
    #-------------------------------------------------------------------------------
    def itos(self,value):
        """
        returns time value (s) from point value  - valid for transients
        """
        return 0.5*value/self.specwidth
    def stoi(self,value):
        """
        returns point value (i) from time value (s)
        """
        return 2.0*value*self.specwidth
    #-------------------------------------------------------------------------------
    def htoi(self,value):
        """
        returns point value (i) from Hz value (h)
        """
#        pt_value = value * (self.size+self.left_point-1)/self.specwidth - self.left_point
        pt_value = (value - self.offsetfreq) * (self.size+self.left_point-1)/self.specwidth  - self.left_point
        return pt_value
    #-------------------------------------------------------------------------------
    def itoh(self, value):
        """
        returns Hz value (h) from point value (i)
        """      
#        hz_value = (value + self.left_point)*self.specwidth / (self.size+self.left_point-1)
        hz_value = (value  + self.left_point)*self.specwidth/(self.size+self.left_point-1) + self.offsetfreq
        return hz_value
    #-------------------------------------------------------------------------------
    def itomz(self, value):
        """
        return m/z (mz) from point value (i)
        """
        return self.htomz(self.itoh(value))
    def mztoi(self, value):
        """
        return point value (i) from  m/z (mz) 
        """
        return self.htoi(self.mztoh(value))
    #-------------------------------------------------------------------------------
    def deltamz(self,mz_value):
        "computes the theorical maximum resolution in m/z at m/z location"

        mz = np.fmin(mz_value, HighestMass)
        return 0.5*(self.itomz( self.mztoi(mz)-1) - self.itomz( self.mztoi(mz)+1) )
    def mass_axis(self):
        """return axis containing m/z values, can be used for display"""
        ax = self.points_axis()
        ax[0] = ax[1]       # m/z is undefined in [0], letting null makes a run time error
        return self.itomz(ax)
    mz_axis = mass_axis     # two names
    def freq_axis(self):
        """return axis containing Hz values, can be used for display"""
        return self.itoh(self.points_axis())
    Hz_axis = freq_axis     # two names for this function
#-------------------------------------------------------------------------------

class FTMSData(NPKData.NPKData):
    """
    subclass of NPKData, meant for handling FT-MS data
    allows 1D and 2D data-sets
    
    """
    def __init__(self, dim=1, shape=None, mode="memory", buffer=None, name=None, debug=0):
        """
        dim : dimension of dataset
        shape : shape of the buffer (size1,size2)
        mode : memory : data-set is kept in-memory    /  onfile : data-set is read from file when needed
        buffer : if is not None; used as data
        name : if is not None, data is read from file
        
        Use as a template, this method should be overwriten by subclasser
        """
        from .File.HDF5File import HDF5File, determine_chunkshape
        self.axis1 = FTMSAxis()    # this creates an FTMSAxis so that pylint does not complain - will be overwritten
        if dim == 2:
            self.axis2 = FTMSAxis()
        if name:
            if name.endswith(".msh5"):  # try loading .msh5 file
                if debug>0: print("reading msh5")
                H = HDF5File(name,"r")
                H.load(mode=mode)      # load into memory by default !
                super(FTMSData, self).__init__(buffer=H.data.buffer, debug=debug)
                NPKData.copyaxes(H.data, self)  # and deep copy all axes from file
                self.name = name
                self.hdf5file = H
            else:
                raise Exception("Filename should have a .msh5 extension")
        else:
            if debug>0: print("calling super")
            super(FTMSData, self).__init__(dim=dim, shape=shape, buffer=buffer, name=name, debug=debug)
            for i in range(self.dim):
                axis = self.axes(i+1)
                setattr(self, "axis%d"%(i+1), FTMSAxis(size=axis.size, specwidth=axis.specwidth, itype=axis.itype) )
        if debug>1: print(self.report())
    #------------------------------------------------
    def _gspecwidth(self):
        "copy specwidth to all the axes"
        return [self.axes(i+1).specwidth for i in range(self.dim)]
    def _sspecwidth(self,sw):
        for i in range(self.dim):
            ax = self.axes(i+1)
            ax.specwidth = sw
    specwidth = property(_gspecwidth, _sspecwidth)
    #------------------------------------------------
    def _grmass(self):
        "copy ref_mass to all the axes"
        return [self.axes(i+1).ref_mass for i in range(self.dim)]
    def _srmass(self,mass):
        for i in range(self.dim):
            ax = self.axes(i+1)
            ax.ref_mass = mass
    ref_mass = property(_grmass, _srmass)
    #------------------------------------------------
    def _grfreq(self):
        "copy ref_freq to all the axes"
        return [self.axes(i+1).ref_freq for i in range(self.dim)]
    def _srfreq(self,ref_freq):
        for i in range(self.dim):
            ax = self.axes(i+1)
            ax.ref_freq = ref_freq
    ref_freq = property(_grfreq, _srfreq)
    #------------------------------------------------
    def _ghm(self):
        "copy highmass to all the axes"
        return [self.axes(i+1).highmass for i in range(self.dim)]
    def _shm(self,highmass):
        for i in range(self.dim):
            ax = self.axes(i+1)
            ax.highmass = highmass
    highmass = property(_ghm, _shm)
    #------------------------------------------------
    def trimz(self, axis=0):
        """
        extract the data so as to keep only lowmass-highmass range
        axis determines which axis to trim, axis=0 (default) indicates all axes
        """
        if axis == 0:       # if 0 --> do all axes
            for i in range(self.dim):
                self.trimz(axis=i+1)
        else:               # do only one
            ext = []
            unitaxis = []
            for i in range(self.dim):
                unitaxis.append( self.axes(i+1).currentunit )  # save
                self.axes(i+1).currentunit = 'points'
                if i+1 == axis:
                    ext = ext + [int(self.axes(i+1).mztoi(self.axes(i+1).highmass)), self.axes(i+1).size]
                else:
                    ext = ext + [0,self.axes(i+1).size]
            print("extracting :", ext)
            self.extract(ext)
            for i in range(self.dim):
                self.axes(i+1).currentunit = unitaxis[i]  # restore
        return self

    #----------------------------------------------
    def save_msh5(self, name, compressed=False):
        """
        save data to a HDF5 file
        
        if compressed is True, file is internally compressed using HDF5 compressors (currently zlib)
        not final version !!
        """
        import tables
        from .File.HDF5File import HDF5File, determine_chunkshape
        hf_file = HDF5File(name,"w")
        if compressed:
            hf_file.set_compression(On=True)
        if self.dim == 2:
            chunks = determine_chunkshape(self.size1, self.size2)
        else:
            chunks = None
        hf_file.create_tables()
        hf_file.create_group("/", "resol1")
        table_axes = hf_file.create_table("/resol1", "axes",hf_file.axis_table)
        
        hf_file.create_tables()
        infos = []
        for i in range(self.dim):
            ii = i + 1
            infos.append( self.axes(ii).__dict__)
        hf_file.fill_table(table_axes, infos)
        if self.dim == 2:
            Hbuffer = hf_file.create_carray("/resol1", 'data', tables.Float64Atom(), (self.size1, self.size2), chunk = chunks)
        else:
            Hbuffer = hf_file.create_carray("/resol1", 'data', tables.Float64Atom(), (self.size1,), chunk =  None)
        Hbuffer[...] = self.buffer[...]
        hf_file.close()
        del(hf_file)
        self.name = name
        return self
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# if __name__ == '__main__':
#     unittest.main()
# no unittest for the moment - tested by subclassers
