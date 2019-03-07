#!/usr/bin/env python 
# encoding: utf-8

"""
The library of processing functions to be used with the NPK program.

This library implement the specific functions needed for MaxEnt processing.
Most of these functions require that the NPK mathematical kernel is loaded.

"""

from __future__ import print_function

__author__ = "Marc A. Delsuc <delsuc@igbmc.fr>"
__date__ = "Oct 2009"

import os.path
import sys
import tempfile
import copy
from Generic import *

#---------------------------------------------------------------------------
def setup_filter2D(si1data,si2data,(lb1,lb2),pos_neg=0):
        """set-up NPK iternal buffers for deconvolution

        filter is created as a dummy 1D, with F2 shape first, and F1 shape afterwards
        """
        dim(2); itype(0)
        szf=max(si1data,si2data)
        if pos_neg:
            chsize(4,szf)
            one()
            row(3);dim(1);mult(-1);dim(2);put("row",3)  # inverse second set
            nchannel(2)
        else:
            chsize(2,szf)
            one()
            nchannel(1)
        for i in range(get_nchannel()):     # usually nchannel==1; so this loop is done once
                row(2*i+1); dim(1)
                em(lb2)
                dim(2); put("row",2*i+1)
        for i in range(get_nchannel()):     # usually nchannel==1; so this loop is done once
                row(2*i+2); dim(1)
                em(lb1)
                dim(2); put("row",2*i+2)
        modifysize(1,get_si1_2d()*get_si2_2d())
        row(1); dim(1); put("filter"); com_filter(2)
        dim(2)
    
#---------------------------------------------------------------------------
def sampling( fname, axis = "F1" ):
    """set-up NPK iternal buffers for applying the MaxEnt analysis on partially sampled data-sets

    fname is the name of the file which contains the sampling function.
    when in 2D, the  axis  parameter is the axis which is currently described.
    takes the value : F1 F2 or F12

    should be called with the data-set to analyse in the main buffer - will use the data buffer

    sampling file format :
    ====================
    in 1D
        one entry per line, with the value of the sampled point

    in 2D
    with axis value  F1 or F2
        one entry per line, with the value of the sampled point along the given axis
    with axis value  F12
        one entry per line, with the two values of the sampled point along the F1 axis first, and F2 axis second
 
    in 2D, calling sampling.g twice in a row with F1 and then with F2 results in a factorized sampling.

    indexes start at 1 and goes up to the size of the sampled axis

    in all cases, lines beginning with a # or a ; are ignored.
eg.
# 1D
1
3
4
...

# 2D F12
1 1
3 2
2 3
...

    returns the number of sampled points in the buffer.

    MAD-VC July 2005 - corrected april 2007
    """
    if ( get_dim() == 1 ):
        axis = "F2"
    elif ( get_dim() == 2 ):
        axis = axis.upper()
    else:
        raise ("sampling is not supported in dim" + str(get_dim()))

    if ( get_dim()== 1 ):
        dimtodo = 1
        window_mode(1)
        itype(0)
        simax=get_si1_1d()
        zero()
    elif ( get_dim() == 2 ):
        dimtodo = 2
        simax = max(get_si2_2d(),get_si1_2d())
        if ( axis == "F1" or axis == "F2"):
            window_mode(1)
            dim(2)
            get("window")
            dim(1)
            chsize(2*simax)
            if ( axis == "F1" ):
                for i in range(simax, 2*simax):
                    setval(i+1, 0.0)
            elif ( axis == "F2" ):
                for i in range(simax):
                    setval(i+1, 0.0)
        elif ( axis == "F12" or axis == "F21" ):
            window_mode(2)
            dim(2)
            zero()

    fin = open(fname)
    f=fin.read()
    lines= f.split("\n")
    fin.close()
    n = 0
    if ( axis == "F2" ):
        for v in lines:
            try:
                setval( int(v), 1.0)
                n = n+1
            except:
                pass
    elif ( axis == "F1" ):
        for v in lines:
            try:
                setval( int(v)+simax, 1.0)
                n = n+1
            except:
                pass
    elif ( axis == "F12" or axis == "F21" ):
        for v in lines:
            try:
                (i,j) = v.split()
                setval( int(i), int(j), 1.0)
                n = n+1
            except:
                pass
    if ( get_dim() == 1 ):
        window_mode(1)
        put("window")
    else:
        dim(2)
        window_mode(2)
        put("window")
    dim( dimtodo )
    return(n)

#---------------------------------------------------------------------------
# reference deconvolution
#---------------------------------------------------------------------------
def ref_deconv(left, right, width ):
    """
    reference deconvolution, uses a given spectral window, and apply deconvolution
    assuming there is a single line in this window.
    See 
    """
    if ( get_dim() != 1 ):
        raise("Still to be done")
    si = get_si1_1d()
    extract(left,right) # extract zone of interest
    sin(0.5)
    dm = int(get_si1_1d()/2)
    reverse()       # and rebuild full size spectrum, with the zone in the middle
    chsize(si/2+dm)
    reverse()
    chsize(si)
    toreal()
    iftbis()

#---------------------------------------------------------------------------
# posneg_corr
#---------------------------------------------------------------------------
def posneg_corr():
    """
    """
    if (get_dim() == 1):
        raise "not written yet"
    elif (get_dim() == 2):
        put("data")
        extract(get_si1_2d()/2, 1, get_si1_2d(), get_si2_2d())
        mult(-1)
        chsize(2*get_si1_2d(), get_si2_2d())
        adddata()
        chsize(get_si1_2d()/2, get_si2_2d())

#------------------------------------------------------------
# Je mets tout dans cet objet, il faudra le couper en morceau
#   - maxent
#   - DOSY (ILTalgo, iltsize)
#   etc..
class maxent_state:
    """ This class maintains a permanent state for MaxEnt processing
    it is a wrapper over the NPK kernel
    it also maintains additional parameters, such as samplings...
    """
    __iltalgo = "MaxEnt"
    __ilttype = "Tabulated"
    __iltsize = 128
    __iter = 20000
    __ndisp = 200
    __lambsp = 30.0
    __miniter = 6
    __lambcont = 2
    __algo = 1
    __linemini = 1
    __wucorrec = 1

    def set_iltalgo(self,st):
        """ string which determines which algo is used for ILT analysis
        MaxEnt or Fit
        """
        if st.upper() == "MAXENT":
            self.__iltalgo = "MaxEnt"
        elif st.upper() == "FIT":
            self.__iltalgo = "Fit"
        else:
            raise "wrong iltalgo value"
    def get_iltalgo(self):
        return self.__iltalgo

    def set_ilttype(self,st):
        """ string which determines which sampling is used for ILT analysis
        REGULAR or TABULATED
        """
        if st.upper() == "REGULAR":
            self.__ilttype = "Regular"
        elif st.upper() == "TABULATED":
            self.__ilttype = "Tabulated"
        else:
            raise "Wrong data type, valid types are REGULAR, TABULATED"
    def get_ilttype(self):
        return self.__ilttype
        

    def set_iltsize(self,x):
        """ size for computing ILT processing
        """
        if (x<1):
            raise "wrong iltsize value"
        self.__iltsize = int(x)
    def get_iltsize(self):
        return self.__iltsize

    def set_iter(self,x):
        if (x<1):
            raise "wrong iter value"
        self.__iter = int(x)
    def get_iter(self):
        return self.__iter

    def set_miniter(self,x):
        if (x<1):
            raise "wrong miniter value"
        self.__miniter = int(x)
    def get_miniter(self):
        return self.__miniter

    def set_ndisp(self,x):
        if (x<1):
            raise "wrong ndisp value"
        self.__ndisp = int(x)
    def get_ndisp(self):
        return self.__ndisp

    def set_lambsp(self,x):
#        if (x<1.0):
#            raise "wrong lambsp value"
        self.__lambsp = float(x)
    def get_lambsp(self):
        return self.__lambsp

    def set_lambcont(self,x):
        if (x not in (0,1,2)):
            raise "wrong lambcont value"
        self.__lambcont = int(x)
    def get_lambcont(self):
        return self.__lambcont

    def set_algo(self,x):
        """ algo determines the global MaxEnt algorithm,
            0 fixed-point with Gull and Daniel equation for Entropy 
            1 fixed-point with Gifa equation for Entropy 
            2 steepest descent equation is used
            3 conjugate gradient is used
        """
        if (x not in (0,1,2,3)):
            raise "wrong algo value"
        self.__algo = int(x)

    def set_linemini(self,x):
        """ xlinemini determines ho line minimisation is performed
            0 single step iteration
            1 line-maximization using parabolic fit.
            2 : line-maximization using bracketing before parabolic
        """
        if (x not in (0,1,2)):
            raise "wrong linemini value"
        self.__linemini = int(x)

    def set_wucorrec(self,x):
        """ determines wether Wu correction is used or not"""
        if (x not in (0,1)):
            raise "wrong wucorrec value"
        self.__wucorrec = int(x)

    def report(self):
        format = """
iltalgo :  %s
iter :     %i
iltsize :  %i
iter :     %i
ndisp :    %i
lambsp :   %f
miniter :  %i
lambcont : %i
algo :     %i
linemini : %i
wucorrec : %i
            """
        return format%(self.__iltalgo,self.__iter,self.__iltsize,self.__iter,self.__ndisp,self.__lambsp,self.__miniter,self.__lambcont,self.__algo,self.__linemini,self.__wucorrec)

    def copyto_kernel(self):
        """ sets the kernel with private values """
        iter(self.__iter)
        ndisp(self.__ndisp)
        lambsp(self.__lambsp)
        miniter(self.__miniter)
        lambcont(self.__lambcont)
        algo(100*self.__wucorrec + 10*self.__linemini + self.__algo)

    def dosy_me_preset(self,  preset = 3 ):
        """
        Presets MaxEnt parameters for DOSY processing
        preset value ranges from 0 to 5
        sets the parameters for a balance between speed (1) and quality (5), 0 is for fit
        """
        print("MaxEnt Preset: "+str(preset))
        if (preset == 1) :
            self.__iltalgo = "MaxEnt"
            self.__iltsize = 64
            self.__iter = 500
            self.__ndisp = 500
            self.__lambsp = 40.0
            self.__miniter = 6
            self.__lambcont = 2 # 1
        elif (preset == 2) :
            self.__iltalgo = "MaxEnt"
            self.__iltsize = 96
            self.__iter = 5000
            self.__ndisp = 500
            self.__lambsp = 40.0
            self.__miniter = 6
            self.__lambcont = 2 # 1
        elif (preset == 3) :
            self.__iltalgo = "MaxEnt"
            self.__iltsize = 128
            self.__iter = 20000
            self.__ndisp = 200
            self.__lambsp = 30.0
            self.__miniter = 6
            self.__lambcont = 2
        elif (preset == 4) :
            self.__iltalgo = "MaxEnt"
            self.__iltsize = 192
            self.__iter = 40000
            self.__ndisp = 200
            self.__lambsp = 10.0
            self.__miniter = 5
            self.__lambcont = 2
        elif (preset == 5) :
            self.__iltalgo = "MaxEnt"
            self.__iltsize = 256
            self.__iter = 100000
            self.__ndisp = 200
            self.__lambsp = 2.0
            self.__miniter = 4
            self.__lambcont = 2
        elif (preset == 0) :
            self.__iltalgo = "Fit"
            self.__iltsize = 256
            self.__iter = 50
            self.__miniter = 6
        else:
            raise "wrong me_preset value"
        
        self.__algo=1   # a verifier
        self.__linemini=1
        self.__wucorrec=0
#        print self.report()

    def me_preset(self,  preset_value=3 ):
        """
        Presets MaxEnt parameters for Fourier processing
        preset value ranges from 1 to 5
        sets the parameters for a balance between speed (1) and quality (5)
        """
        if ( preset_value == 1 ):
            self.__iter = 5
            self.__ndisp = 5
            self.__lambsp = 3.0         # lowering lambsp makes smoother avancement (and slower convergence)
        elif ( preset_value == 2):
            self.__iter = 10
            self.__ndisp = 5
            self.__lambsp = 3.0
        elif ( preset_value == 3):
            self.__iter = 50
            self.__ndisp = 5
            self.__lambsp = 2.0
        elif ( preset_value == 4):
            self.__iter = 200
            self.__ndisp = 10
            self.__lambsp = 2.0
        elif ( preset_value == 5):
            self.__iter = 1000
            self.__ndisp = 10
            self.__lambsp = 1.5
        else:
            raise "wrong me_preset value"
#         stepmax  1.2
        tolerance(1.0e-2)
#         step  1.0
        self.__miniter = 4
        self.__lambcont = 1
        self.__algo = 1
        
