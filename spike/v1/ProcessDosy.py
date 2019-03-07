#!/usr/bin/env python 
# encoding: utf-8

"""
    This module contains all the routines needed to process DOSY spectra
    
    The following functions realize a set of operations, needed for DOSY processing,
    It is based on an implementation of the Inverse Fourier Transform, by Maximum Entropy
    
    Processing are divided in pre - Transform - Post phases.
    For instance a typical DOSY procesing, would be

    pre_2d()
    pre_ft_2df2()
    ft_2df2()
    post_ft_2df2()
    pre_ilt_2df1()
    ilt_2df1()
    post_ilt_2df1()
    post_2d()

    but other combinations are possible.

    Most functions take the following arguments :
    arguments :
    audit
        the opened audit file, if empty or set to none, audit will go to the stdout
    filein      the file name of the input file, 
        will be loaded before operation, if name is "memory", then processing takes place on the current 2D kernel buffer
    fileout     the file name of the input file, can be "memory"
        will be written after operation, if name is "memory", then processing results is left in 2D kernel buffer
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
from GenericDosy import *
from GenericMaxEnt import *
from Process2D import read_file_2d, write_file_2d, pre_2d, pre_ft_2df2, ft_2df2, post_ft_2df2, post_2d


#---------------------------------------------------------------------------
def PreDosy2D(audit,p_in_arg,f_in,f_out,inputfilename,outputfilename):
    """
    prepare a 2D for processing of a 2D DOSY experiment

    ends with a half processed file
    uses  :
        pre_2d()        ( from Process2D )
        pre_ft_2df2()   ( from Process2D )
        ft_2df2()       ( from Process2D )
        post_ft_2df2()  ( from Process2D )
        pre_dosy()
    """
    p_in = build_dict(("pre_2d","pre_ft_2df2","ft_2df2","post_ft_2df2","post_2d","pre_dosy"),p_in_arg)
    
    # Preparation phase
    audittrail( audit, "phase", "FID analysis phase")
    pre_2d( audit, inputfilename, "memory", p_in, f_in, f_out)

    # Pre FT F2 phase
    audittrail( audit, "phase", "F2_axis - Pre FT phase")
    pre_ft_2df2( audit, "memory", "memory", p_in, f_in, f_out)

    # FT F2 phase
    audittrail( audit, "phase", "F2 axis - FT phase")
    ft_2df2( audit, "memory", "memory", p_in, f_in, f_out )

    # Post FT F2 phase
    audittrail( audit, "phase", "F2 axis - Post FT phase")
    post_ft_2df2( audit, "memory", "memory", p_in, f_in, f_out )

    # DOSY baseline correction
    audittrail( audit, "phase", "F2 axis - Baseline correction")
    pre_dosy( audit, "memory", outputfilename, p_in, f_in, f_out )

#---------------------------------------------------------------------------
def PreDosy3D(audit,p_in_arg,f_in,f_out,inputfilename,outputfilename):
    """
    prepare a 3D for processing of a 3D DOSY experiment

    experiment should be organized as a decaying series along F1 of regular 2D along F2-F3.
    ends with a half processed file
    uses  :
        process3DF2F3()        ( from Process3D )
        pre_dosy()
    """
    p_in = build_dict(("pre_2d","pre_ft_2df2","ft_2df2","post_ft_2df2","post_2d","pre_dosy"),p_in_arg)

    # Preparation phase
    audittrail( audit, "phase", "FID analysis phase")
    pre_2d( audit, inputfilename, "memory", p_in, f_in, f_out)

    # Pre FT F2 phase
    audittrail( audit, "phase", "F2_axis - Pre FT phase")
    pre_ft_2df2( audit, "memory", "memory", p_in, f_in, f_out)

    # FT F2 phase
    audittrail( audit, "phase", "F2 axis - FT phase")
    ft_2df2( audit, "memory", "memory", p_in, f_in, f_out )

    # Post FT F2 phase
    audittrail( audit, "phase", "F2 axis - Post FT phase")
    post_ft_2df2( audit, "memory", "memory", p_in, f_in, f_out )

    # DOSY baseline correction
    audittrail( audit, "phase", "F2 axis - Baseline correction")
    pre_dosy( audit, "memory", outputfilename, p_in, f_in, f_out )

#---------------------------------------------------------------------------
def Dosy2D(audit,p_in_arg,f_in,f_out,inputfilename,outputfilename):
    """
    processing of a 2D DOSY experiment

    starts with a half processed file
    """

    p_in = build_dict(("pre_ilt_2df1","ilt_2df1","post_ilt_2df1"),p_in_arg)

    # Pre ILT phase
    audittrail( audit, "phase", "F1_axis - Pre ILT phase")
    pre_ilt_2df1( audit, inputfilename, "memory", p_in, f_in, f_out )

    # ILT phase
    audittrail( audit, "phase", "2D ILT phase")
    ilt_2df1( audit, "memory", "memory", p_in, f_in, f_out )

    # Post ILT F1 phase
    audittrail( audit, "phase", "2D Post ILT phase")
    post_ilt_2df1( audit, "memory", outputfilename, p_in, f_in, f_out )

#---------------------------------------------------------------------------
def Dosy3D(audit,p_in_arg,f_in,f_out,inputfilename,outputfilename):
    """
    processing of a small 3D DOSY experiment
    larger DOSY3D should be processed with Process3D.DOSY3D()
    
    starts with a half processed file
    """
# this is just a terrible wrapper over DOSY2D
    audittrail( audit, "phase", "3D flattening")
    if ( inputfilename != "memory" ):
        try:
            com_read( inputfilename )
        except:
            raise "Error while reading file ",inputfilename
        if ( get_dim() != 3 ):
            raise "This is not a 3D dataset"
    (si1,si2,si3)=(get_si1_3d(),get_si2_3d(),get_si3_3d())
    modifysize( 1, si1, si2*si3 )
    plane("F1",1)
    dim(2)
    audittrail( audit, "text", "resulting 2D","F1",get_si1_2d(),"F2",get_si2_2d())
    
    p_in = build_dict(("pre_ilt_2df1","ilt_2df1","post_ilt_2df1"),p_in_arg)

    # Pre ILT phase
    audittrail( audit, "phase", "F1_axis - Pre ILT phase")
    pre_ilt_2df1( audit, "memory", "memory", p_in, f_in, f_out )

    # ILT phase
    audittrail( audit, "phase", "2D ILT phase")
    ilt_2df1( audit, "memory", "memory", p_in, f_in, f_out )

    # Post ILT F1 phase
    audittrail( audit, "phase", "2D Post ILT phase")
    post_ilt_2df1( audit, "memory", "memory", p_in, f_in, f_out )

    # tail of the terrible hack
    audittrail( audit, "phase", "3D unflattening")
    dim(3)
    si1 = get_si1_2d()
    chsize(1, si1, si2*si3)
    put("plane","F1",1)
    modifysize(si1, si2, si3 )

    writec(outputfilename)
    audittrail( audit, "text", "resulting 3D",
        "F1",si1,"F2",si2,"F3",si3,
        "size in F1 x F2", si1*si2*si3,
        "file name", outputfilename)

#---------------------------------------------------------------------------
def pre_dosy(audit, filein, fileout, p_in, f_in, f_out):
    """
This macro is simply a wrapper around proc_2d() to realise a f2 baseline correction at the end of the PreDosy2D processing
This is somewhat ad-hoc, and may be changed sometimes

it implements
- f2_baseline   : baseline correction
- select_state  : enforce the final state of the data

%F2-input-domain%  frequency
%F2-output-domain% frequency

%author% Marc-Andre Delsuc
%version% 6.0

    """
    # %action% baseline correction
    #   baseline correction of the hald-processed DOSY spectrum
    # %param% f2_baseline boolean / default 1
    # %param% f2_bcorr_algo enum  offset linear spline / default linear
    #    offset : removes a automatically determined offset
    #    linear - spline : uses the list of user determined pivot points to define the baseline,
    #       : fit with a straight line or a spline
    # %action% change other post_2d() default actions
    # %param% autophase_2d boolean / default 0
    # %param% f1_baseline boolean / default 0
    # %param% projection boolean / default 0

    # %action% select_state
    #   permit to choose the state complex / real of the output file
    #   complex data are changed to real by dropping the imaginary part
    #   real data are changed to complex by computing the Hilbert transform with tocomplex()
    #   this is usually not required as all processing commands prepare the state themselves 
    # %param% select_state boolean / default 1
    #   actuallly does the selection
    # %param% f1_state enum ignore complex real / default ignore
    #   force the f1 axis to real or complex. ignore will let the axis unchanged
    # %param% f2_state enum ignore complex real / default real
    #   force the f2 axis to real or complex. ignore will let the axis unchanged

    post_2d( audit, filein, fileout, p_in, f_in, f_out)
#---------------------------------------------------------------------------
def pre_ilt_2df1( audit, filein, fileout, p_in, f_in, f_out ):

    """
This macro realizes the pre operation on each F1 decay of the 2D

it implements data massaging before ILT:
- ilt_set-up
    dosy_q2 dosy_tab tosy_q2 tosy_tab
- load tabulated file
    file name
- calibration
- F1 noise evaluation
- wrong points filtering

%F1-input-domain%  tabulated
%F1-output-domain% tabulated
%dimensionality%   2

%author% Marc-Andre Delsuc
%version% 6.0

    """

    read_file_2d( audit, filein)

	# if complex in t1, assumes it is an error
    if (get_itype_2d()>=2):
        itype(get_itype_2d()-2)
        audittrail( audit,  "text", "assuming real data-set along F1 !")

	# data should be real in F2!
	toreal("f2")

    # %action% ilt_set-up
    #   determines the kind of analysis to perform
    # %param% ilt_type enum dosy tosy / default dosy
    # %param% sampling_type enum regular tabulated / default regular
    # %param% sampling_tab_file file_name
    #   sampling_type == tabulated implies sampling_tab_file
    set_task("ILT set-up and calibration")

    if (not key_is_true(p_in,"ilt_type")):
        p_in["ilt_type"] = "dosy"

    if (not key_is_true(p_in,"sampling_type")):
        p_in["sampling_type"] = "regular"

    if (p_in["sampling_type"] == "tabulated"):
        load_sq_tab(p_in["sampling_tab_file"])
        audittrail( audit, "text", "ILT set-up",
                       "type of analysis", p_in["ilt_type"],
                       "type of sampling", p_in["sampling_type"],
                       "sampling file", p_in["sampling_tab_file"])
    elif (p_in["sampling_type"] == "regular"):
        audittrail( audit, "text", "ILT set-up",
                       "type of analysis", p_in["ilt_type"],
                       "type of sampling", p_in["sampling_type"])
    else:
        raise ("wrong sampling_type :" + str(p_in["sampling_type"]))

    # %action% ilt_tosy_calibration
    #   Calibrate ilt axis on physical parameters for Tosy
    # %param% ilt_tosy_unit enum usec msec sec / default sec 
    #   value of the Laplace axis unity
    if ( p_in["ilt_type"] == "tosy"):
        if (key_is_true(p_in,"ilt_tosy_unit")):
            if (p_in["ilt_tosy_unit"] not in ("usec","msec","sec") ):
                raise "wrong value for ilt_tosy_unit, accepted values are usec, msec, sec "
        else:
            p_in["ilt_tosy_unit"] ="sec"

    # %action% ilt_dosy_calibration
    #   Calibrate ilt axis on physical parameters for Dosy
    # %param% dosy_big_delta float
    #   "Big Delta"  : diffusion delay in msec
    # %param% dosy_little_delta float
    #   "little Delta"  : gradient duration in msec
    # %param% dosy_seq_typ enum pgse ste bpp_ste ste_2echoes bpp_ste_2echoes oneshot / default ste
    #   the type of DOSY sequence used
    # %param% dosy_nucleus enum 1H 2H 13C 15N 17O 19F 31P / default 1H
    #   the observed nucleus
    # %param% recovery_gradient_delay float / default 0.0
    # %param_cond%  (recovery_gradient_delay >= 0 ) 
    #   Gradient recovery delay
    # %param% max_grad float nonnegative / default 50.0
    #   Maximum Amplificator Gradient Intensity, in G/cm    
    # %param_cond%  (max_grad > 0 )
    # %param% max_tab float / default 100.0
    #   Maximum Tabulated Gradient Value in the tabulated file. 
    #   Bruker users use 100 here
    #   Varian users use 32768 here
    # %param_cond%  (max_tab > 0 )
    # %param% gradient_shape float / default 1.0
    #   integral factor depending on the gradient shape used 
    #   typical values are :
    #       1.0 for rectangular gradients
    #       0.6366 = 2/pi for sine bell gradients
    #       0.4839496 for 4% truncated gaussian (Bruker gauss.100 file)
    # %param_cond%  [ (gradient_shape > 0.0 ) and (gradient_shape <= 1.0 )]
    
    if (key_is_true(p_in,"ilt_dosy_calibration")):
        set_task("automatic calibration")
        bigdelta = p_in["dosy_big_delta"]
        litdelta = p_in["dosy_little_delta"]
        if (key_is_true(p_in,"dosy_seq_typ")):
            seq_type = p_in["dosy_seq_typ"]
            if seq_type not in ("pgse","ste","bpp_ste","ste_2echoes","bpp_ste_2echoes","oneshot"):
                raise "wrong DOSY sequence type"
        else:
            seq_type="ste"
        if (key_is_true(p_in,"dosy_nucleus")):
            nucleus = p_in["dosy_nucleus"]
            if nucleus not in ("1H","2H","13C","15N","17O","19F","31P"):
                raise "unknown nucleus type"
        else:
            nucleus = "1H"
        if (key_is_true(p_in,"recovery_gradient_delay")):
            recovery = p_in["recovery_gradient_delay"]
        else:
            recovery = 0.0
        if (key_is_true(p_in,"max_grad")):
            maxgrad = p_in["max_grad"]
        else:
            maxgrad = 50.0
        if (key_is_true(p_in,"max_tab")):
            maxtab = p_in["max_tab"]
        else:
            maxtab = 100.0
        if (key_is_true(p_in,"gradient_shape")):
            gradshape = p_in["gradient_shape"]
        else:
            gradshape = 1.0
        dfactor = calibdosy(litdelta, bigdelta, recovery, seq_type, nucleus, maxgrad, maxtab, gradshape)
        audittrail( audit, "text", "DOSY Calibration",
                       "sequence type", seq_type,
                       "DFACTOR value", dfactor)
        p_in["dfactor"]=dfactor                       
        f_out["dfactor"]=dfactor
        print("dfactor: %f"%(dfactor))
    # F1 noise evaluation
    # %action% f1_noise
    #   evaluate noise, estimated by finding an empty zone on a given 1D spectrum
    # %param% f1_noise boolean / default 1
    # %param% f1_noise_n integer / default 30
    #    number of different zones where noise is evaluated
    # %param% f1_noise_row integer / default get_si1_2d()
    #    the index of the row on which the noise is evaluated
    # %param% f1_noise_weight float / default 1.0
    #   the noise used during the MaxEnt iteration is weigthed with this scalar,
    #   this permits to compensate for over- or under-estimate of the noise level
    # %param_cond%   [ noise_weight > 0.0 ]
    # %return% noise_in_f1_domain
    if (key_is_not_false(p_in,"f1_noise")):
        set_task("Noise evaluation")
        if (not key_is_true(p_in,"f1_noise_n")):
            p_in["f1_noise_n"] = 30
        if (not key_is_true(p_in,"f1_noise_row")):
            p_in["f1_noise_row"] = get_si1_2d()
        dim(2)
        row(p_in["f1_noise_row"])
        dim(1)
        spec_noise( p_in["f1_noise_n"] )
        dim(2)
        if key_is_true(p_in,"f1_noise_weight"):
            weight=p_in["f1_noise_weight"]
        else:
            weight=1.0
        if weight!=1.0:
            noise(weight*get_noise())
        f_out["noise_in_f1_domain"] = get_noise()
        audittrail( audit, "text", "Estimated F1 Noise",
                       "Estimated on row #", p_in["f1_noise_row"],
                       "Noise weighting", weight,
                       "Noise estimate", f_out["noise_in_f1_domain"])
    else:
        if key_is_true(p_in,"f1_noise_weight"):
            noise(p_in["f1_noise_weight"]*get_noise())
        f_out["noise_in_f1_domain"] = get_noise()
        audittrail( audit, "text", "Input F1 Noise",
                       "Noise weighting", weight,
                       "Noise estimate", f_out["noise_in_f1_domain"])
    # wrong points filtering
    # %action% f1_ignore_point
    #   enter the list of the point to ignore along the (gradient) axis
    # %param% f1_ignore_point boolean / default 0
    # %param% f1_ignore_point_list (list of integer) / default ()
	#   list of points to ignore

# MAD    recopier ce code dans un objet "window"
    dim(2)
#    window("F1",0)   # first resets window  MAD
    if (key_is_true(p_in,"f1_ignore_point")):
        if key_is_true(p_in,"f1_ignore_point_list"):
            t=p_in.raw("f1_ignore_point_list")
            list = t[1:len(t)-1].split(",")
            RemoveRows(list)
            writec("RR.gs2")
            dim(1)
            get("tab")
            RemovePoints(list)
            put("tab")
            writec("RP.gs1")
            audittrail( audit, "text", "Point ignored along F1 axis",
                       "point list :", str(list))
        else:
            raise "F1_ignore_point_list undefined"
    else:
        pass
#        com_window_mode(0)

    write_file_2d( audit, fileout)
# MAD : AD HOC, To be CHANGED !
    dim(1)
    window_reset()
    window_mode(0)

#---------------------------------------------------------------------------
def ilt_2df1( audit, filein_arg, fileout, p_in, f_in, f_out):
    """
This macro realizes the ILT operation on the F1 decay of the 2D

it implements the spectral analysis step :


- damp_width        dmin dmax
- minimum S/N
- col_selection
- me_preset
- me_details
    size
    iteration
    ndisp
    lambcont
    lambsp
    miniter
- ILT
- reverse           : reverses laplace axis


%F1-input-domain%  tabulated
%F1-output-domain% damping
%dimensionality%   2

%author% Marc-Andre Delsuc
%version% 6.0

    """
    read_file_2d( audit, filein_arg)

	# if complex in t1, assumes it is an error
    if (get_itype_2d()>=2):
        itype(get_itype_2d()-2)
        audittrail( audit,  "text", "assuming real data-set along F1 !")

	# data should be real in F2!
    toreal("F2")

    tempofile = NPKtempfile('.gs2')
    writec(tempofile)
    filein=tempofile
    
    # %action% damp_width
    #   sets-up the border for the laplace spectrum
    # %param% damp_width enum diff_standard from_tab defined / default diff_standard
    #   defines how the window will be defined
    #       diff_standard is the set-up for a standard solution : 10-5000 um^2 sec-1
    #       from_tab : the borders are computed from the tabulated values
    #       defined : the values 
    # %param% dmin float > 0  / default 10
    # %param% dmax float > 0  / default 5000
    # %param_cond% [ dmin < dmax ]

    set_task("Setting damping window")
    if(not key_is_true(p_in,"damp_width")):
        p_in["damp_width"] = "diff_standard"
    if (p_in["damp_width"] == "diff_standard"):
        p_in["dmin"] = 10
        p_in["dmax"] = 5000
    elif (p_in["damp_width"] == "from_tab"):
        p_in["dmin"],p_in["dmax"] = auto_damp_width()
    elif (p_in["damp_width"] == "defined"): 
        if(not key_is_true(p_in,"dmin")):
            p_in["dmin"] = 10
        if(not key_is_true(p_in,"dmax")):
            p_in["dmax"] = 5000
    else:
        raise "wrong value for damp_width"
    if (p_in["dmin"] >= p_in["dmax"]):
        raise "wrong values for Dmin and Dmax in parameter file"
    
    audittrail( audit, "text", "damping axis window",
                       "Dmin", p_in["dmin"],
                       "Dmax", p_in["dmax"])

    # %action% minimum_sn
    #   defines the minimum Signal/Noise ratio required for ILT processing
    #   No data will ever be processed if the S/N is below this value, independently of the col_selection set-up
    #   It is considered unsafe to go below 10 for this parameter.
    # %param% minimum_sn float / default 32.0
    if(not key_is_true(p_in,"minimum_sn")):
        p_in["minimum_sn"] = 32
    audittrail( audit, "text", "minimum Signal/Noise ratio",
                       "minimum_SN", p_in["minimum_sn"])

    # %action% col_selection
    #     set-up the conditions for selecting the column to be processed by ILT for DOSY
    # %param% col_selection boolean / default 0
    # %param% col_selection_mode enum valthreshold zone list none / default none
    #   none : no selection (only minimum_sn stands)
    #   valthreshold : absolute value theshold, will not be used if below minimum_sn
    #   zone : only points within that spectral zone will be considered
    #   list : explicit list of the columns to consider
    #
    # %param% col_valthreshold float
    #
    # %param% col_list (integer list)
    #   list og the columns to consider, given in index (real part only)
    #
    # %param% col_selection_left float / default 10.0
    #   the left border of the selection zone, in unit
    # %param% col_selection_left_unit enum ppm hz index / default ppm
    #   the unit in which selection_left is given
    # %param% col_selection_right float / default 0.0
    #   the right border of the selection zone, in unit
    # %param% col_selection_right_unit enum ppm hz index / default ppm
    #   the unit in which selection_right is given
    # %param_cond% (col_selection_left{index} < col_selection_right{index})
    # %return% col_selection_left
    #     the left coordinate of the selection zone in index
    # %return% col_selection_right
    #     the right coordinate of the selection zone in index
    #

    # set-up default
    set_task("Columns checking")
    to_check = range(1,get_si2_2d()+1)  # inital list of col to check
    thresh = get_noise()*p_in["minimum_sn"]
    start = 1
    end = get_si1_2d()
    if(key_is_true(p_in,"col_selection")):
    # first setup selection
        if ("zone" == p_in["col_selection_mode"]):
                if (p_in["col_selection_left_unit"] == 'ppm'):
                    start = (ptoi(p_in["col_selection_left"],1,1))
                elif (p_in["col_selection_left_unit"] == 'hz'):
                    start = (htoi(p_in["col_selection_left"],1,1))
                elif (p_in["col_selection_left_unit"] == 'index'):
                    start = p_in["col_selection_left"]
                else:
                    raise "Error with col_selection_left_unit"

                if (p_in["col_selection_right_unit"] == 'ppm'):
                    end = (ptoi(p_in["col_selection_right"],1,1))
                elif (p_in["col_selection_right_unit"] == 'hz'):
                    end = (htoi(p_in["col_selection_right"],1,1))
                elif (p_in["col_selection_right_unit"] == 'index'):
                    end = p_in["col_selection_right"]
                else:
                    raise "Error with col_selection_right_unit"

                start = (max(1,round(start)))
                end = (min(get_si2_2d(),round(end)))
                if (start > end):
                    raise "Wrong selection zone coordinates"
                to_check = range(start,end+1)
        # set-up threshold
        elif ("valthreshold" == p_in["col_selection_mode"]):
            thresh = max(thresh, p_in["col_valthreshold"])
        elif ("list" == p_in["col_selection_mode"]):
#            to_check = p_in["col_list"]
            t = p_in.raw("col_list")
            to_check = t[1:len(t)-1].split(",")     # remove ( ) ad split with ","
                #eval() should work here, but we have seen a bug for very long strings with jython 2.1
        else:
            raise "wrong selection mode"

    # then do selection
    dim (2) # load in dim(1) the first row of the DOSY file
    row(1)
    dim(1)
    l=[]
    for i in to_check:
        if val1d(int(i)) > thresh:
            l.append(i)

    f_out["col_selection_left"] = start
    f_out["col_selection_right"] = end
    audittrail( audit, "text", "column selection",
                    "number of selected columns", len(l))
    print("Number of selected colums: %i"%(len(l)))

    # %action% me_preset
    #   presets the MaxEnt parameters to default values
    # %param% me_preset_value enum 0 1 2 3 4 5 / default 3
    #   sets the parameter for a balance between speed (1) and quality (5), 0 is for fit
    set_task("Setting Maximum Entropy parameters")
    me = maxent_state()
    if (key_is_true(p_in,"me_preset_value")):
        preset = p_in["me_preset_value"]
    else:
        preset = 3
    me.dosy_me_preset(preset)
    
    # %action% me_details
    #   if this flag is on, default parameters can be set.
    # %param% me_details boolean / default 0
    # %param% me_size integer 
    # %param% me_iteration integer
    # %param% me_int_iteration
    # %param% me_ncheck integer
    # %param% me_lambda_control integer
    # %param% me_lambda_speed
    # %param% me_algo enum Fit MaxEnt
    if (key_is_true(p_in,"me_details")):
        if (key_is_true(p_in,"me_size")):
           me.set_iltsize(p_in["me_size"])
        if (key_is_true(p_in,"me_iteration")):
           me.set_iter(p_in["me_iteration"])
        if (key_is_true(p_in,"me_int_iteration")):
           me.set_miniter(p_in["me_int_iteration"])
        if (key_is_true(p_in,"me_ncheck")):
           me.set_ndisp(p_in["me_ncheck"])
        if (key_is_true(p_in,"me_lambda_control")):
           me.set_lambcont(p_in["me_lambda_control"])
        if (key_is_true(p_in,"me_lambda_speed")):
           me.set_lambsp(p_in["me_lambda_speed"])
        if (key_is_true(p_in,"me_algo")):
           me.set_iltalgo(p_in["me_algo"])

    # set it in kernel
    me.copyto_kernel()
    audittrail(audit, "text","Parameters for ILT processing","MaxEnt parameters",me.report())

    # %action% dump_preset
    #   dump all the parameters to a gtb file
    #   usually used to check the parameter or to prepare for a parallel run
    # %param% dump_preset boolean / default 0
    # %param% dump_preset_file file / default dosy_preset.gtb
    if (key_is_true(p_in,"dump_preset")):
        p=p_in.copy()  # first copy everything, then overwrite specific values
        p["damp_width"]="defined"
        p["noise"]=f_out["noise_in_f1_domain"]
        p["col_selection"]=1
        p["col_selection_mode"]="list"
        p["col_list"]=tuple(l)
        p["me_details"]=1
        p["me_size"]=me.get_iltsize() 
        p["me_iteration"]=me.get_iter()
        p["me_int_iteration"]=me.get_miniter()
        p["me_ncheck"]=me.get_ndisp()
        p["me_lambda_control"]=me.get_lambcont()
        p["me_lambda_speed"]=me.get_lambsp()
        p["me_algo"]=me.get_iltalgo()
        p["execute_ilt"]=1
        p["dump_preset"]=0  # remove dump_preset
        if (key_is_true(p_in,"dump_preset_file")):
            preset_file = p_in["dump_preset_file"]
        else:
            preset_file = "dosy_preset.gtb"
        dict_dump(p,preset_file)
        audittrail(audit, "text", "parameter dump",
            "Parameters for ILT processing have been dumped to",preset_file)

    # %action% ILT
    #   perform the ILT computation itself
    # %param% execute_ilt boolean / default 1
    if (key_is_not_false(p_in,"execute_ilt")):
        set_task("Starting ilt processing %d"%(len(l)))
        t = time.clock()
        dosy2d(filein,p_in["dmin"],p_in["dmax"],p_in["dfactor"],me,l)   # hum, a changer ! MAD
        t = time.clock()-t
        audittrail( audit, "text", "2D ILT computed",
            "number of column computed", len(l),
            "processing time", str(t)+" sec")

    # %action% f1_reverse
    # %param% f1_reverse boolean / default 0
    if(key_is_true(p_in,"f1_reverse")):
        set_task("Reverse F1 axis")
        reverse("f1")
        audittrail( audit, "text", "Reverse F1 spectral axis")

    (sw,fq,off)=damp_to_ppm(p_in["dmin"],p_in["dmax"])	# compute and set pseudo ppm values
    dim(2)
    specw(sw,get_specw_2_2d())
    freq(get_freq(), fq, get_freq_2_2d())
    offset(off,get_offset_2_2d())

    write_file_2d( audit, fileout)
    os.unlink(tempofile)


#---------------------------------------------------------------------------
def post_ilt_2df1( audit, filein, fileout, p_in, f_in, f_out ):
    """
This macro realizes the post F1-ILT operation on a 2D spectrum

- smoothing     : F2 smoothing of DOSY

%F2-input-domain%  time
%F2-output-domain% time
%dimensionality%   2

%author% Marc-Andre Delsuc
%version% 6.0

    """
    read_file_2d( audit, filein)


    write_file_2d( audit, fileout)
