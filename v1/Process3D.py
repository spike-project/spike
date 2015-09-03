#!/usr/bin/env python 
# encoding: utf-8

"""
    This module contains all the routines needed to process 3D NMR spectra
    
    The following functions realize a set of operations, needed for NMR 3D processing,
    
    Processing are divided in pre - FT - Post phases.
    For instance a typical 3D procesing, would be

    pre_3d()
    process3DF2F3()
    process3DF1()
    post_3d()
   
    or

    pre_3d()
    process3DF3()
    process3DF1F2()
    post_3d()

    or

    pre_3d()
    process3DF3()
    process3DF2()
    process3DF1()
    post_3d()

    3D nomenclature is :
    F1 is the slowest varying axis of the 3D
    F2 is the intermediate axis
    F3 is the acquisition axis.
    Note that this is not equivalent to the "spectroscopic" view of 3F, were F1 is associated to the first incremental time of the excitation sequence,
    which may, or may not, be the slowest varying axis.
    
    Additionally,
    F1 planes are planes orthogonal to the F1 axis, thus containing the F2-F3 axes of the 3D.

    3D data-set are not loaded in memory but rather planes extracted from the 3D are loaded planewise,
    and processed as regular in-memory 2D.
    
    Proces3Dxx phases are thus compound actions built upon the 2D couterpart.
    For instance Process3DF2F3 realizes the 2D processing of each F2-F3 planes (aka orthogonal to F1) of the 3D

    all parameters and actions are thus equivalent to there 2D counter part
    (F1 and F2 3D-axes are equivalent to the F1 axis of the 2D, and the F3 axis of the 3D is equivalent to the F2 axis of the 2D)
    They are not cited here but should be checked in the Process2D documentation.
    
    Standard 2D default parameters are used so far, additionnal parameters can be given using the F1- F2- F3 nomenclature.



    Most functions take the following arguments :
    arguments :
    audit
        the opened audit file, if empty or set to none, audit will go to the stdout
    filein      the file name of the input file, 
        will loaded before operation, if name is "memory", then processing takes place on the current 3D kernel buffer
    fileout     the file name of the input file, can be "memory"
        will written after operation, if name is "memory", then processing results is left in 3D kernel buffer
    p_in        the dictionary containing all the processing parameters
        most entries are optionnal as entries are protected with try / except
    f_in        the dictionary containing all the details on the spectrum "filein"
    f_out       details on the spectrum "fileout" will be put in this dictionary
    location    the directory where the output files will be written
    
"""

from __future__ import print_function

__author__ = "Marc A. Delsuc <delsuc@igbmc.fr>"
__date__ = "Oct 2009"


import math
import os
import time
from Generic import *
from Process2D import *

#---------------------------------------------------------------------------
def FT3D(audit,p_in_arg,f_in,f_out,inputfilename,outputfilename):
    """
    FT processing of a 3D FID

    uses  :
        pre_3d()
        process3DF2F3()
        process3DF1()
        post_3d()
    """
    import tempfile
    import os

    p_in = build_dict(("pre_3d","process3DF2F3","process3DF1","post_3d"),p_in_arg)
# this version does not use pre phase, as it is empty for the moment!

    # Initial phase
#    audittrail( audit, "3D-phase", "3D: Initial phase")
#    pre_3d( audit, tempdata3, outputfilename, p_in, f_in, f_out, "." )

    # F2-F3 axes
    tempdata2 = tempfile.mktemp('.npk3D2')
    print(tempdata2)
    audittrail( audit, "3D-phase", "3D: F2-F3 planes processing - based on the following 2D phases")
    process3DF2F3( audit, inputfilename, tempdata2, p_in, f_in, f_out )
    audittrail( audit, "phase", "All F2-F3 axes processed")
    audittrail( audit, "text", "details",
        "from",filein,
        "to",fileout,
        "number of planes processed",n)

    # F1 axis
#    tempdata3 = tempfile.mktemp('.npk3D3')
#    print tempdata3
    audittrail( audit, "3D-phase", "3D: F1 axis processing  - based on the following 2D phases")
    process3DF1( audit, tempdata2, outputfilename, p_in, f_in, f_out )
    audittrail( audit, "phase", "F1 axis processed")
    audittrail( audit, "text", "details",
        "from",filein,
        "to",fileout,
        "number of planes processed",n)
    try:
        os.unlink(tempdata2)
        audittrail( audit, "text", "temporary file removed", "file", tempdata2)
    except:
        pass

    # Final phase
    audittrail( audit, "3D-phase", "3D: Final phase")
    post_3d( audit, outputfilename, p_in, f_in, f_out, "." )
    try:
        os.unlink(tempdata3)
        audittrail( audit, "text", "temporary file removed", "file", tempdata3)
    except:
        pass

#---------------------------------------------------------------------------
def pre_3d( audit, filein, fileout, p_in_arg, f_in, f_out ):
    """
    Pre processing of 3D
    Empty for the moment
    """
    pass

#---------------------------------------------------------------------------
def process3DF2F3( audit, filein, fileout, p_in_arg3D, f_in_arg, f_out ):

    """
    This macro realizes the 2D processing of each F2-F3 planes (aka orthogonal to F1) of the 3D
    """ 
    import Process2D
    L_p_in = p_in_arg3D.change_key('f0','f1').change_key('f1','f2').change_key('f2','f3')
#    L_p_in = change_key_dict('f2', 'f3',  change_key_dict('f1', 'f2', p_in_arg3D))     # in THAT order !
    L_f_in = change_key_dict('f2', 'f3',  change_key_dict('f1', 'f2', f_in_arg))     # in THAT order !
    L_f_out={}
    
    join(filein)
    dim(2)
    (si1,si2) = plane_size("F1")
    n = get_c_sizef1()
    getc("F1",1,1,1,si1,si2)
    audittrail( audit, "text", "first plane loaded in memory for test")
    FT2D(audit,L_p_in,L_f_in,f_out,"memory",fileout+"_trial.gs2")
    # print L_p_in
    # raise "stop"
    
    command = 'Process2D.FT2D("mute", L_p_in, L_f_in, L_f_out, "memory", "memory")'  # action to do
    dict_ = globals()
    for i in locals().keys():
         dict_[i]=locals()[i]
    proc3d( filein, fileout, 'F1', command, dict_)

#---------------------------------------------------------------------------
def process3DF1F2( audit, filein, fileout, p_in_arg, f_in_arg, f_out ):

    """
    This macro realizes the 2D processing of each F1-F2 planes (aka orthogonal to F3) of the 3D
    """ 
    
    print(filein)
    join(filein)
    dim(2)
    (si1,si2)=plane_size("F3")
    n=get_c_sizef1()
    getc("F3",1,1,1,si1,si2)
    audittrail( audit, "text", "first plane loaded in memory for test")
    FT2D(audit,p_in,f_in,f_out,"memory","memory")
    
    command = 'FT2D("mute", p_in, f_in, f_out, "memory", "memory")'  # action to do
    dict_ = globals()
    for i in locals().keys():
        dict_[i]=locals()[i]
    proc3d( filein, fileout, 'F3', command, dict_)

#---------------------------------------------------------------------------
def process3DF1( audit, filein, fileout, p_in_arg, f_in_arg, f_out ):

    """
    This macro realizes the 1D processing of each F1-FID along the F1 (slowest varying) axis of the 3D
    """
    import Process2D
    p_in = change_key_dict('f2', 'f3', p_in_arg)
    p_in["f2_baseline"] = 0 # no F2 baseline
    f_in = change_key_dict('f2', 'f3', f_in_arg)

    print(filein)
    join(filein)
    dim(2)
    (si1,si2)=plane_size("F2")
    n=get_c_sizef1()
    getc("F2",1,1,1,si1,si2)
    audittrail( audit, "text", "first plane loaded in memory for test")
    FT_F1_2D(audit,p_in,f_in,f_out,"memory","memory")

    command = 'FT_F1_2D("mute",p_in,f_in,f_out,"memory","memory")'   # action to do
    dict_ = globals()
    for i in locals().keys():
         dict_[i]=locals()[i]
    proc3d( filein, fileout, 'F2', command, dict_)

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
def process3DF2( audit, filein, fileout, p_in_arg, f_in_arg, f_out ):

    """
    This macro realizes the 1D processing of each F2-FID along the F2 (intermediate varying) axis of the 3D
    """

    p_in = change_key_dict('f2', 'f3', p_in_arg)
    f_in = change_key_dict('f2', 'f3', f_in_arg)

    join(filein)
    dim(2)
    (si1,si2)=plane_size("F2")
    n=get_c_sizef1()
    getc("F2",1,1,1,si1,si2)
    audittrail( audit, "text", "first plane loaded in memory for test")
    FT_F1_2D(audit,p_in,f_in,f_out,"memory","memory")

    command = 'FT_F1_2D("mute",p_in,f_in,f_out,"memory","memory")'   # action to do
    dict = globals()
    for i in locals().keys():
        dict[i]=locals()[i]
    proc3d( filein, fileout, 'F2', command,dict)

#---------------------------------------------------------------------------
def process3DF3( audit, filein, fileout, p_in_arg, f_in_arg, f_out ):

    """
    This macro realizes the 1D processing of each F3-FID along the F3 (classical) axis of the 3D
    """

    p_in = change_key_dict('f2', 'f3', p_in_arg)
    f_in = change_key_dict('f2', 'f3', f_in_arg)

    join(filein)
    dim(2)
    (si1,si2)=plane_size("F1")
    n=get_c_sizef1()
    getc("F2",1,1,1,si1,si2)
    audittrail( audit, "text", "first plane loaded in memory for test")
    FT_F2_2D(audit,p_in,f_in,f_out,"memory","memory")

    command = 'FT_F1_2D("mute",p_in,f_in,f_out,"memory","memory")'   # action to do
    dict = globals()
    for i in locals().keys():
        dict[i]=locals()[i]
    proc3d( filein, fileout, 'F2', command,dict)

#---------------------------------------------------------------------------
def post_3d( audit, file, p_in_arg, f_in, f_out, location ):
    """
    Post processing of 3D
    - projection
    
    every thing is performed in-place

    will be added :
    - calibration
    """
    p_in = p_in_arg
    # %action% 3D_projections
    #  creates and store in files the projections of the 3D
    #  2 files are created per axis :
    #  using a mean algorithm (_M suffix) and using a skyline algorithm (_S suffix)
    # %param% 3D_projections boolean / default 1
    # %return% F1_projections
    # %return% F2_projections
    # %return% F3_projections
    if (key_is_true(p_in,"3D_projections")):
        set_task("Projection computation")
        join(file)
        for axis in ("F1","F2","F3"):
            f_out[axis+"_projections"] = ""
            for algo in ("S","SA"):
                local_proj_3d(axis, algo)
                dim(2)
                projname = os.path.join(location, 'proj%s_%s.gs2'%(axis,algo))
                writec(projname)
                f_out[axis+"_projections"] += projname+" "
        audittrail( audit,  "text", "Projections stored",
                        "F1 projection files", f_out["F1_projections"],
                        "F2 projection files", f_out["F2_projections"],
                        "F3 projection files", f_out["F3_projections"])
        disjoin()


