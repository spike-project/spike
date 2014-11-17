
# encoding: utf-8
"""
Orbitrap.py

Created by Marc-André and Lionel on 10 april 2014
Copyright (c) 2014 IGBMC. All rights reserved.
"""

import math
import unittest
import numpy as np
import File.HDF5File 
import NPKData
from NPKError import NPKError
try:
    from util.ValErr import ValErr
except:
    print "No matplot"

HighestMass = 100000.0  # the highest m/z of interest ever
FREQ0 = 1E7
REF_FREQ = 1887533.975611561     # Ad Hoc values for   REF_FREQ / REF_MASS
REF_MASS = 715.3122

class OrbiAxis(NPKData.Axis):
    """
    hold information for one Orbitrap axis
    used internally
    """
    #A FAIRE
    
    def __init__(self, size=1024, specwidth=FREQ0, itype=0, units="point", ref_mass=REF_MASS, ref_freq=REF_FREQ, highmass=10000.0, left_point = 0.0):
        """
        all parameters from Axis, plus
        specwidth   highest frequency,
        ref_mass    m/z used for calibration
        ref_freq    frequency of the calibration m/z
        highmass    highest m/z of interest - usually defined by the excitation pulse low frequency limit
        left_point  coordinates of first data point; usually 0.0 after Fourier Transform; may be different after extraction
        
        possible values for unit are "points" "m/z"
        
        conversion methods work on numpy arrays as well
        """
        super(OrbiAxis, self).__init__(size=size, itype=itype, units=units)
        self.specwidth = specwidth
        self.ref_mass = ref_mass
        self.ref_freq = ref_freq
        self.highmass = highmass
        self.left_point = left_point
        self.unit_types.append("m/z")
        self.unit_types.append("Hz")
        self.Orbitrap = "Orbitrap"
        for i in ("specwidth", "ref_mass", "ref_freq", "highmass", "left_point", "Orbitrap"):  # updates storable attributes
            self.attributes.insert(0,i)
    #-------------------------------------------------------------------------------
    def _report(self):
        "low level report"
        return "size : %d   sw %f  ref_mass %f ref_freq %f  left_point %f  highmass %f  itype %d units %s"%(self.size, self.specwidth, self.ref_mass, self.ref_freq, self.left_point, self.highmass,self.itype, self.units)
    def report(self):
        "high level reporting"
        
        if self.itype == 0: # Real
            return "Orbitrap axis at %f kHz,  %d real points,  from mz = %8.3f   to m/z = %8.3f  M/DeltaM (M=400) = %.0f"%  \
            (self.specwidth/1000, self.size, self.highmass, self.lowmass, 400.0/self.deltamz(400.))
        else: # Complex
            return "Orbitrap axis at %f kHz,  %d complex pairs,  from mz = %8.3f   to m/z = %8.3f  M/DeltaM (M=400) = %.0f"%  \
            (self.specwidth/1000, self.size/2, self.highmass, self.lowmass, 400.0/self.deltamz(400.))

    #-------------------------------------------------------------------------------
    @property
    def lowmass(self):
        """highest mass of interest - defined by the Nyquist frequency limit"""
        return self.itomz(self.size-1)
    #-------------------------------------------------------------------------------
    def extract(self, (start, end)):
        """redefines the axis parameters so that the new axe is extracted for the points [start:end]"""
        #   before      left----------------------------size
        #   before      -------------------------------SW
        #   after             start---------------end
        #                     left'---------------SW'
        if end <0: end = self.size+end+1
        if start<0 or end>self.size or start>=end:
            raise NPKError("The current axis contains %d points, it cannot be extracted from %d to %d"%(self.size, start, end))
        self.specwidth = self.itoh(end-1)
        self.left_point += int(start)
        self.size = int(end)-int(start)
    #-------------------------------------------------------------------------------
    def htoi(self,value):
        """
        returns point value (i) from Hz value (h)
        """
        pt_value = value * (self.size+self.left_point-1)/self.specwidth - self.left_point
        return pt_value
    #-------------------------------------------------------------------------------
    def itoh(self, value):
        """
        returns Hz value (h) from point value (i)
        """      
        hz_value =   (value + self.left_point)*self.specwidth / (self.size+self.left_point-1)
        
        return hz_value
    def itomz(self, value):
        """
        return m/z (mz) from point value (i)
        """
        return self.ref_mass / (self.itoh(value)/self.ref_freq)**2
 
    def mztoi(self, value):
        """
        return point value (i) from  m/z (mz) 
        """
        return self.htoi(np.sqrt(self.ref_mass/value)*self.ref_freq)
    def deltamz(self,mz_value):
        "computes the theorical resolution in m/z at m/z location"
        mz = np.fmin(mz_value, HighestMass)
        return 0.5*(self.itomz( self.mztoi(mz)-1) - self.itomz( self.mztoi(mz)+1) )
    def mass_axis(self):
        ax = self.points_axis()
        ax[0] = ax[1]       # m/z is undefined in [0], letting null makes a run time error
        return self.itomz(ax)
    mz_axis = mass_axis     # two names
    def freq_axis(self):
        """return axis containing Hz values, can be used for display"""
        return self.itoh(self.points_axis())
    Hz_axis = freq_axis     # two names for this function
#-------------------------------------------------------------------------------
class OrbiData(NPKData.NPKData):
    """
    subclass of NPKData, meant for handling Orbitrap data
    doc to be written ...
    """
#    print "in Orbitrap"
    def __init__(self, dim=1, shape=None, mode="memory", buffer=None, name=None, debug=0):
        self.axis1 = OrbiAxis()    # this creates an OrbiAxis so that pylint does not complain - will be overwritten
        if dim == 2:
            self.axis2 = OrbiAxis()
        if name:
            if name.endswith(".msh5"):  # try loading .msh5 file
                if debug>0: print "reading msh5"
                H = File.HDF5File.HDF5File(name,"r")
                H.load(mode=mode)      # load into memory by default !
                super(OrbiData, self).__init__(buffer=H.data.buffer, debug=debug)
                NPKData.copyaxes(H.data, self)  # and deep copy all axes from file
                self.name = name
                self.hdf5file = H
            else:
                if debug>0: print "reading gifafile"
                super(OrbiData, self).__init__(name=name, debug=debug)
                for i in range(self.dim):
                    axis = self.axes(i+1)
                    setattr(self, "axis%d"%(i+1), OrbiAxis(size=axis.size, specwidth=axis.specwidth, itype=axis.itype) )
        else:
            if debug>0: print "calling super"
            super(OrbiData, self).__init__(dim=dim, shape=shape, buffer = buffer, name = name, debug=debug)
            for i in range(self.dim):
                axis = self.axes(i+1)
                setattr(self, "axis%d"%(i+1), OrbiAxis(size=axis.size, itype=0) )
        if debug>1: print self.report()
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
        return [self.axes(i+1).ref_mass for i in range(self.dim)]
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
    def _gunits(self):
        "copy units to all the axes"
        return [self.axes(i+1).units for i in range(self.dim)]
    def _sunits(self, units):
        for i in range(self.dim):
            ax = self.axes(i+1)
            ax.units = units
    units = property(_gunits, _sunits)
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
            for i in range(self.dim):
                if i+1 == axis:
                    ext = ext + [int(self.axes(i+1).mztoi(self.axes(i+1).highmass)), self.axes(i+1).size]
                else:
                    ext = ext + [0,self.axes(i+1).size]
            print "extracting :", ext
            self.extract(ext)
        return self
    #------------------------------------------------
    def display(self, scale = 1.0, absmax = 0.0, show = False, label = None, new_fig = True, axis = None, mode3D = False, zoom = None, xlabel = "_def_", ylabel = "_def_", figure = None):

        """display the Orbitrap data using NPKDATA display method
        check parameters in NPKDATA

        - copied here - (might be out of date)
        scale   allows to increase the vertical scale of display
        absmax  overwrite the value for the largest point, which will not be computed
            display is scaled so that the largest point is first computed (and stored in absmax),
            and then the value at absmax/scale is set full screen 
        show    will call plot.show() at the end, allowing every declared display to be shown on-screen
                useless in ipython
        label   add a label text to plot
        xlabel, ylabel : axes label (default is self.units - use None to remove)
        axis    used as axis if present, axis length should match experiment length
                in 2D, should be a pair (xaxis,yaxis)
        new_fig will create a new window if set to True (default) (active only is figure==None)
        mode3D  use malb 3D display instead of matplotlib contour for 2D display
        zoom    is a tuple defining the zomm window (left,right) or   ((F1_limits),(F2_limits))
        figure  if not None, will be used directly to display instead of using its own
        
        can actually be called without harm, even if no graphic is available, it will just do nothing.
        
        """
        if self.debug>0: print "zoom at the beginning of display ",zoom,axis
        if self.dim == 1:
            if zoom is None:
                zm = [0,self.size1]
            else:
                zm = list(zoom)     # zoom can be passed as a tuple
            if not self.check_zoom(zm):
                raise NPKError("wrong zoom window : %s "%(str(zm)), data=self)
            if self.axis1.units == "m/z" :
                # if using m/z units, we have to build a zoom window and an axis
                axis = self.axis1.mass_axis()
                zm[0] = int(max(zm[0], self.axis1.mztoi(self.axis1.highmass)))    # restrict to highmass
                if not self.check_zoom(zm):
                    raise NPKError("zoom window is below high-mass limit: %s "%(str(zm)), data=self)
        elif self.dim == 2:
            if zoom is None:
                zm = [ [0,self.size1], [0,self.size2] ]
            else:
                zm = list(zoom)     # zoom can be passed as a tuple
            #print "liste zoom", zm
            if not self.check_zoom(zm):
                raise NPKError("wrong zoom window : %s "%(str(zm)), data=self)
            if self.axis1.units == "m/z" or  self.axis2.units == "m/z":

                # if using m/z units, we have to build a zoom window and an axis
                if axis is None:
                    axis = [0,1]
                for (i,ax) in ( (0,self.axis1), (1,self.axis2) ):
                    if ax.units == "m/z":
                        axis[i] = ax.mass_axis()
                        left = int(max(zm[i][0], ax.mztoi(ax.highmass)))    # restrict to highmass
                        #print "truncated to high-mass :",left
                        zm[i] = [left, zm[i][1]]
                    else:
                        ax[i] = np.arange(ax.size)
            if not self.check_zoom(zm):
                raise NPKError("zoom window is below high-mass limit: %s "%(str(zm)), data=self)
        if self.debug>0: print "in Orbitrap zoom if m/z ",zoom, zm
        # then super.display()
        #print "in Orbitrap zoom ",zoom
        
        super(OrbiData, self).display(scale = scale, absmax = absmax, axis = axis, zoom = zm, show = show, label = label, new_fig = new_fig, mode3D = mode3D, xlabel = xlabel, ylabel = ylabel, figure = figure)
        return self

    #----------------------------------------------
    def save_msh5(self, name):
        """
        save data to a HDF5 file
        
        experimental !
        """
        import tables
        hf_file = File.HDF5File.HDF5File(name,"w")
        if self.dim == 2:
            chunks = File.HDF5File.determine_chunkshape(self.size1, self.size2)
        else:
            chunks = None
        hf_file.create_tables()
        hf_file.createGroup("/", "resol1")
        table_axes = hf_file.createTable("/resol1", "axes", hf_file.axis_table)
        
        hf_file.create_tables()
        infos = []
        for i in range(self.dim):
            ii = i + 1
            infos.append( self.axes(ii).__dict__)
        hf_file.fill_table(table_axes, infos)
        if self.dim == 2:
            Hbuffer = hf_file.createCArray("/resol1", 'data', tables.Float64Atom(), (self.size1, self.size2), chunk = chunks)
        else:
            Hbuffer = hf_file.createCArray("/resol1", 'data', tables.Float64Atom(), (self.size1,), chunk =  None)
        Hbuffer[...] = self.buffer[...]
        hf_file.close()
        del(hf_file)
        self.name = name
        return self
#-------------------------------------------------------------------------------
class Orbi_Tests(unittest.TestCase):
    def setUp(self):
        self.verbose = 1    # verbose > 0 switches messages on
    def announce(self):
        if self.verbose >0:
            print "\n========",self.shortDescription(),'==============='
    def test_atob(self):
        "testing unit conversion functions"
        self.announce()
        Oaxis = OrbiAxis(size = 1000, specwidth = 1667000, itype = 0, units = "point", ref_mass = 344.0974, ref_freq = 419620.0, highmass = 2000.0)
        self.assertAlmostEqual(Oaxis.itoh(0), 0)
        self.assertAlmostEqual(Oaxis.itoh(Oaxis.size-1), Oaxis.specwidth)   # last point is size-1 !!!
        self.assertAlmostEqual(Oaxis.itomz(1023)/Oaxis.deltamz(Oaxis.itomz(1023)), 1023./2, places = 2)    # delta_m / m is 1/size at highest mass - before extraction
        self.assertAlmostEqual(123**2*Oaxis.itomz(123), 321**2*Oaxis.itomz(321))    # verify that f*m/z is constant !
        self.assertAlmostEqual(123*Oaxis.mztoi(123)**2, 321*Oaxis.mztoi(321)**2)    # verify that f*m/z is constant !
        for i in (1,2):
            print Oaxis.report()
            print Oaxis._report()   # low level report
            for x in (1.0, 301.0, Oaxis.size-20.0, Oaxis.size-1.0):    # last point is size-1 !!!
                print "point at index %d is at freq %f, m/z %f"%(x, Oaxis.itoh(x), Oaxis.itomz(x))
                self.assertAlmostEqual(Oaxis.mztoi(Oaxis.itomz(x)), x)
                self.assertAlmostEqual(Oaxis.itoh(Oaxis.htoi(x)), x)
            Oaxis.extract([300,-20])
            
    def test_trim(self):
        """
        Test trimz 
        """
        O = OrbiData(buffer=np.zeros((500, 10000)))
        O.specwidth = 1667000
        O.ref_mass = 344.0974
        O.ref_freq = 419620.0
        O.highmass = 1000.0
        print O.report()
        l1 = int(O.axis1.mztoi(O.axis1.highmass))
        l2 = int(O.axis2.mztoi(O.axis2.highmass))
        O.trimz()
        self.assertEqual(l1, O.axis1.left_point)
        self.assertEqual(l2, O.axis2.left_point)
        self.assertEqual(O.size1, 500-l1)
        self.assertEqual(O.size2, 10000-l2)
        print "2D trimz gain is : %d %%" % (100*(1-(O.size1*O.size2/(500.*10000))))

#-------------------------------------------------------------------------------
if __name__ == '__main__':
    unittest.main()

