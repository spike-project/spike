#!/usr/bin/env python 
# encoding: utf-8

"""
    This module contains all the routines needed to process 2D NMR spectra
    
    The following functions realize a set of operations, needed for NMR 2D processing,
    You will find operation for FT analysis, MaxEnt Analysis,
    Inverse Fourier, etc..
    
    Processing are divided in pre - FT - Post phases.
    For instance a typical 2D procesing, would be

    pre_2d()
    pre_ft_2df2()
    ft_2df2()
    post_ft_2df2()
    pre_ft_2df1()
    ft_2df1()
    post_ft_2df1()
    post_2d()

    but other combinations are possible.

    Most functions take the following arguments :
    arguments :
    audit
        the opened audit file, if empty or set to none, audit will go to the stdout
    filein      the file name of the input file, 
        will loaded before operation, if name is "memory", then processing takes place on the current 2D kernel buffer
    fileout     the file name of the input file, can be "memory"
        will written after operation, if name is "memory", then processing results is left in 2D kernel buffer
    p_in        the dictionary containing all the processing parameters
        most entries are optionnal as entries are protected with try / except
    f_in        the dictionary containing all the details on the spectrum "filein"
    f_out       details on the spectrum "fileout" will be put in this dictionary
    location    the directory where the output files will be written

    Most operation are performed in-memory, and filein and fileout are not a requisite.
    It is possible to just realize the processing from the 2D data-set in memory, leaving the result in memory.
    
    To do this, simple use "memory" instead of a real filename.
    
"""

from __future__ import print_function

__author__ = "Marc A. Delsuc <delsuc@igbmc.fr> and Vincent Catherinot <v.catherinot@nmrtec.com>"
__date__ = "Oct 2009"

import math
import os
import os.path
import time
from Generic import *
from GenericMaxEnt import *

#---------------------------------------------------------------------------
def FT2D(audit,p_in_arg_,f_in,f_out,inputfilename,outputfilename):
    """
    FT processing of a 2D FID

    uses  :
        pre_2d()
        pre_ft_2df2()
        ft_2df2()
        post_ft_2df2()
        pre_ft_2df1()
        ft_2df1()
        post_ft_2df1()
        post_2d()
    """
    p_in = build_dict(("pre_2d","pre_ft_2df2","ft_2df2","post_ft_2df2","pre_ft_2df1","ft_2df1","post_ft_2df1","post_2d"),p_in_arg_)
        
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

    # Pre FT F1 phase
    audittrail( audit, "phase", "F1_axis - Pre FT phase")
    pre_ft_2df1( audit, "memory", "memory", p_in, f_in, f_out )

    # FT F1 phase
    audittrail( audit, "phase", "F1 axis - FT phase")
    ft_2df1( audit, "memory", "memory", p_in, f_in, f_out )

    # Post FT F1 phase
    audittrail( audit, "phase", "F1 axis - Post FT phase")
    post_ft_2df1( audit, "memory", "memory", p_in, f_in, f_out )

    # Final phase
    audittrail( audit, "phase", "Final phase")
    post_2d( audit, "memory", outputfilename, p_in, f_in, f_out, "." )

#---------------------------------------------------------------------------
def FT_F2_2D(audit,p_in_arg,f_in,f_out,inputfilename,outputfilename):
    """
    processing of the F2 axis of a 2D FID
    used by 3D processing

    uses  :
        pre_2d()
        pre_ft_2df2()
        ft_2df2()
        post_ft_2df2()
    """
    p_in = build_dict(("pre_2d","pre_ft_2df2","ft_2df2","post_ft_2df2"),p_in_arg)

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
    post_ft_2df2( audit, "memory", outputfilename, p_in, f_in, f_out )

#    # Final phase
#    audittrail( audit, "phase", "Final phase")
#    post_2d( audit, "memory", outputfilename, p_in, f_in, f_out, "." )

#---------------------------------------------------------------------------
def FT_F1_2D(audit,p_in_arg,f_in,f_out,inputfilename,outputfilename):
    """
    processing of the F1 axis of a 2D FID
    used by 3D processing

    uses  :
        pre_2d()
        pre_ft_2df1()
        ft_2df1()
        post_ft_2df1()
        post_2d()
    """
    p_in = build_dict(("pre_2d","pre_ft_2df1","pre_ft_2df2","ft_2df1","post_ft_2df1","post_2d"),p_in_arg)
    # Preparation phase
    audittrail( audit, "phase", "FID analysis phase")
    pre_2d( audit, inputfilename, "memory", p_in, f_in, f_out)

    # Pre FT F1 phase
    audittrail( audit, "phase", "F1_axis - Pre FT phase")
    pre_ft_2df1( audit, "memory", "memory", p_in, f_in, f_out )

    # FT F1 phase
    audittrail( audit, "phase", "F1 axis - FT phase")
    ft_2df1( audit, "memory", "memory", p_in, f_in, f_out )

    # Post FT F1 phase
    audittrail( audit, "phase", "F1 axis - Post FT phase")
    post_ft_2df1( audit, "memory", "memory", p_in, f_in, f_out )

    # Final phase
    audittrail( audit, "phase", "Final phase")
    post_2d( audit, "memory", outputfilename, p_in, f_in, f_out, "." )

#---------------------------------------------------------------------------
def MaxEnt2D(audit,p_in_arg,f_in,f_out,inputfilename,outputfilename):
    """
    MaxEnt processing of a 2D FID

    uses  :
        pre_2d()
        pre_ft_2df2()
        pre_ft_2df1()
        maxent_2d()
        post_maxent_2d()
        post_2d()
    """
    p_in = build_dict(("pre_2d","pre_ft_2df2","pre_ft_2df1","maxent_2d","post_maxent_2d","post_2d"),p_in_arg)

    # Preparation phase
    audittrail( audit, "phase", "FID analysis phase")
    pre_2d( audit, inputfilename, "memory", p_in, f_in, f_out)

    # Pre FT F2 phase
    audittrail( audit, "phase", "F2_axis - Pre FT phase")
    pre_ft_2df2( audit, "memory", "memory", p_in, f_in, f_out)

    # Pre FT F1 phase
    audittrail( audit, "phase", "F1_axis - Pre FT phase")
    pre_ft_2df1( audit, "memory", "memory", p_in, f_in, f_out )


    # MaxEnt phase
    audittrail( audit, "phase", "2D MaxEnt phase")
    maxent_2d( audit, "memory", "memory", p_in, f_in, f_out )

    # Post FT F1 phase
    audittrail( audit, "phase", "2D Post MaxEnt phase")
    post_maxent_2d( audit, "memory", outputfilename, p_in, f_in, f_out )

#    # Final phase
#    audittrail( audit, "phase", "Final phase")
#    post_2d( audit, "memory", outputfilename, p_in, f_in, f_out, "." )

#---------------------------------------------------------------------------
def ft_2df1( audit, filein, fileout, p_in, f_in, f_out):
    """
This macro realizes the FT operation on the F1 FID of the 2D

it implements the spectral analysis step :
- truncate          : remove points at the end of the FID
- lp_extend         : extend the FID by LP analysis
- apodisation       : multiply the FID by some windowing function 
- Fourier_transform : performs FFT
- reverse           : reverses spectral axis

%F1-input-domain%  time
%F1-output-domain% frequency
%dimensionality%   2

%author% Marc-Andre Delsuc
%version% 6.0

    """
    debug = 0
    read_file_2d( audit, filein)
    if debug: writec("TTT.gs2")
    
    # %action% truncate
    #   truncates along the F1 dimension by removing the last points
    # %param% f1_truncate boolean / default 0
    # %param% f1_trunc_size integer > 0 <= get_si1_2d() / default get_si1_2d()
    #   FID size after truncation
    # %param_cond% [ (f1_trunc_size > 0 ) and (f1_trunc_size <= get_si1_2d() ) ]
    if(key_is_true(p_in,"f1_truncate")):
        if debug: print("trucation")
        set_task("F1 truncation")
        s = p_in["f1_trunc_size"]
        if (s <= 0):
            raise "Cannot truncate to zero"

        if (s > get_si1_2d()):
            raise "Error with truncation size (larger than actual size)"

        initial = get_si1_2d()
        chsize(min(int(p_in["f1_trunc_size"]),get_si1_2d()), get_si2_2d())
        audittrail( audit, "text", "F1 truncation before FT",
                   "initial size", initial,
                   "final size", p_in["f1_trunc_size"])

    # %action% f1_lp_extend
    #    extend FID with LP algorithm
    # %param% f1_lp_extend boolean / default 0
    # %param% f1_lp_ext_size integer > si1_2d / default 2*get_si1_2d()
    #   final size of FID
    # %param_cond% [ f1_lp_ext_size > get_si1_2d() ]
    # %param% f1_lp_ext_algo enum burg mirror lpsvd lpsvd_stable / default burg
    #   algorithm used
    #   burg and miror are much faster, svd_stable is much slower
    #   miror is to be used when the phase of the spectrum is known before hand (see lp_ext_off)
    # %param% f1_lp_ext_off integer / default 0
    #   offset determines the position of the t=0 point, used by mirror algo
    #     0  no shift : in-phase data set.
    #    -1  acquisition started exactly half a dwell after t=0 - (will need phase 0 180)
    #    n>0 acquisition started exactly n dwell before t=0
    # %param% f1_lp_ext_order integer / default 20
    #    the size of the prediction polynomial used for LP, a rough estimate of the complexity
    # %param_cond% [ f1_lp_ext_order < get_si1_2d()/2 ]
    # %param% f1_lp_ext_apod boolean / default 1
    #   apply a sine bell apodisation after LP extenstion.
    if (key_is_true(p_in,"f1_lp_extend")):
        set_task("F1 linear prediction extension")
        order(int(p_in["f1_lp_ext_order"]))

        if (p_in["f1_lp_ext_algo"] == "burg"):
            burg2d("f1", p_in["f1_lp_ext_size"])

        elif (p_in["f1_lp_ext_algo"] == "mirror"):
            burg2d_mirror("f1", p_in["lp_ext_off"], p_in["lp_ext_size"])   # MAD A VERIF

        elif (p_in["f1_lp_ext_algo"] == "lpsvd"):
            dt2svd(get_si1_1d())
            svd2ar(1)
            ar2dt(p_in["f1_lp_ext_size"], 1)  # MAD A VERIF
            raise "LPSVD Not debugged!"

        elif (p_in["f1_lp_ext_algo"] == "lpsvd_stable"):
            dt2svd(get_si1_1d())
            svd2ar(1)
            ar2rt(1)
            rtreflect(1)
            rt2ar(1)
            ar2dt(p_in["f1_lp_ext_size"], 1)
            raise "LPSVD Not debugged!"

        audittrail( audit, "text", "F1 extension before FT",
                   "algorithm", p_in["f1_lp_ext_algo"],
                   "LP order", p_in["f1_lp_ext_order"])

        if (key_is_not_false(p_in,"f1_lp_ext_apod")):
            sin(0, "f1")
            audittrail( audit, "text", "Apodisation after extension", "function", "sin(0,f1)")
        if debug: writec("TTT.gs2")
        
    # %action% apodize
    #   standard apodisation by a window
    # %param% f1_apodize boolean / default 1
    # %param% f1_apodisation string / default "sin(0.1)"
    #   the string describing the apodisation to apply
    #    "apodisation" should be a suite of window apodisation commands
    #     ie : 'sin(0.5)' or 'exbroad(0.2);sin(0.5)'
    #     the following predefined functions are implemnted
    #     sin (sine bell) sqsin (squared sine bell) expbroad (exponential)
    #     gaussbroad (gaussian) gaussenh (gaussian enhancement)

    if (key_is_true(p_in,"f1_apodize") and (p_in.raw("f1_apodisation") != "none")):
        set_task("F1 apodisation")
        apodise(p_in.raw("f1_apodisation"), "F1")
        audittrail( audit, "text", "F1 apodisation before FT",
                   "function", p_in.raw("f1_apodisation"))

    # %action% Fourier_transform
    #  performs the Fourier transform
    # %param% f1_fourier_transform boolean / default 1
    # %param% f1_ask_for_ftsize boolean / default 0
    #  if set, the final size is determined by f1_ft_size, automatically determined otherwise
    #   size for FT is determined by  size or lp_ext_size or trunc_size
    #   it will be extended to the next 2^n
    # %param% f1_ft_size integer  / default 2*power2(get_si1_2d())
    #   the size for FT if f1_ask_for_ftsize is true
    # %param_cond% (f1_ft_size == power2(n)
    # %param% f1_ft_type enum none ft_sh ft_tppi ft_sh_tppi ft_n_p ft_phase_modu ft rft / default ft_sh
    #    the Fourier transform algorithm, depends on the acquisition scheme
    # %return% f1_size_for_ft
    #   contains the data size right after FT (will usually be different from final size)

    if (key_is_true(p_in,"f1_fourier_transform")):
        set_task("F1 Fourier Transform")
        if (key_is_true(p_in,"f1_ask_for_ftsize")):
            final = int(p_in["f1_ft_size"])
        else:
            # try   size >= lp_ext_size and >= 2*initial size and == 2^n
            if (p_in["f1_lp_extend"]):
                final = 2*power2(get_si1_2d() - 2)
            else:
                final = 4*power2(get_si1_2d() - 2)


        initial = get_si1_2d()

        chsize(final, get_si2_2d())
        if (p_in["f1_ft_type"] != "none"):
            if (p_in["f1_ft_type"] == "ft"):
                # if ($itype_2d == 0 | $itype_2d == 2) itype (%+1)
                ft("f1")
            elif (p_in["f1_ft_type"] == "rft"):
                #  if ($itype_2d == 1 | $itype_2d == 3) itype (%-1)
                rft("f1")
            elif (p_in["f1_ft_type"] == "ft_sh"):
                ft_sh("f1")
            elif (p_in["f1_ft_type"] == "ft_tppi"):
                ft_tppi("f1")
            elif (p_in["f1_ft_type"] == "ft_sh_tppi"):
                ft_sh_tppi("f1")
            elif (p_in["f1_ft_type"] == "ft_n_p"):
                ft_n_p("f1")
            elif (p_in["f1_ft_type"] == "ft_phase_modu"):
                ft_phase_modu("f1")
            else:
                # executer la macro @($p_in[f1_ft_type])
                raise "FT Mode not implemented yet"

        f_out["f1_size_for_ft"] = final
        audittrail( audit,  "text", "F1 Fourier transformation,",
                    "FT type", p_in["f1_ft_type"],
                    "initial size ", initial,
                    "final size ", final)

    # %action% f1_reverse
    # %param% f1_reverse boolean / default 0
    #   value depends on spectrometer
    if(key_is_true(p_in,"f1_reverse")):
        set_task("F1 Spectrum reversal")
        reverse("f1")
        if (get_itype_2d()>1):
            invf("f1")       # if complex, inverse imaginary part
        audittrail( audit, "text", "Reverse F1 spectral axis")

    write_file_2d( audit, fileout)

#---------------------------------------------------------------------------
def ft_2df2( audit, filein, fileout, p_in, f_in, f_out):
    """
This macro realizes the FT operation on the F2 FID of the 2D

it implements the spectral analysis step :
- truncate          : remove points at the end of the FID
- lp_extend         : extend the FID by LP analysis
- apodisation       : multiply the FID by some windowing function 
- Fourier_transform : performs FFT
- causal_corr       : perform causal correction
- reverse           : reverses spectral axis

%F2-input-domain%  time
%F2-output-domain% frequency
%dimensionality%   2

%author% Marc-Andre Delsuc
%version% 6.0

    """
    debug = 0
    read_file_2d( audit, filein)
    
    # %action% truncate
    #   truncates the FID by removing the last points
    # %param% f2_truncate boolean / default 0
    # %param% f2_trunc_size integer / default get_si2_2d()
    #   new FID size after truncation
    # %param_cond% [ (f2_trunc_size > 0 ) and (f2_trunc_size <= get_si2_2d() ) ]
    if (key_is_true(p_in,"f2_truncate")):
            set_task("F2 truncation")
            s = p_in["f2_trunc_size"]
            if (s <= 0):
                raise "Cannot truncate to zero"

            if (s > get_si2_2d()):
                raise "Error with truncation size (larger than actual size)"

            initial = get_si2_2d()
            chsize(get_si1_2d(), min(int(p_in["f2_trunc_size"]),get_si2_2d()))
            audittrail( audit, "text", "FID truncation before FT",
                       "initial size", initial, "final size", get_si2_2d())

    # %action% lp_extend
    #    extend FID with Linear Prediction algorithm
    # %param% f2_lp_extend boolean / default 0
    # %param% f2_lp_ext_size integer / default 2*get_si2_2d()
    #   final size of FID 
    # %param_cond% [ f2_lp_ext_size > get_si2_2d() ]
    # %param% f2_lp_ext_algo enum burg mirror lpsvd lpsvd_stable / default burg
    #   algorithm used
    #   burg and mirror are much faster, svd_stable is much slower
    #   mirror is to be used when the phase of the spectrum is known before hand (see f2_lp_ext_off)
    # %param% f2_lp_ext_off integer / default 0
    #   offset determines the position of the t=0 point, used by mirror algo
    #     0               no shift : in-phase data set.
    #    -1               acquisition started exactly half a dwell after t=0 - (will need phase 0 180)
    #    f2_lp_ext_off>0  acquisition started exactly n dwell before t=0
    # %param_cond% [ f2_lp_ext_off > -1 ]
    # %param% f2_lp_ext_order integer / default 40
    #    the size of the prediction polynomial used for LP, a rough estimate of the complexity
    # %param_cond% [ (f2_lp_ext_order < get_si2_2d()/4 ) and (f2_lp_ext_order > 2) ]
    # %param% f2_lp_ext_apod boolean / default 1
    #   apply a sine bell apodisation after LP extenstion.
    try:
        keytest = p_in["f2_lp_extend"]
    except:
        p_in["f2_lp_extend"]=0  # will be used by FT
    else:
        if(keytest):
            set_task("F2 linear prediction extension")
            order(int(p_in["f2_lp_ext_order"]))
            if (p_in["f2_lp_ext_algo"] == "burg"):
                burg2d( "f2", p_in["f2_lp_ext_size"])

            elif (p_in["f2_lp_ext_algo"] == "mirror"):
                burg2d_mirror("f2", p_in["f2_lp_ext_off"],p_in["f2_lp_ext_size"])   # MAD A VERIF

            elif (p_in["f2_lp_ext_algo"] == "lpsvd"):
                dt2svd(get_si1_1d())
                svd2ar(1)
                ar2dt(p_in["lp_ext_size"], 1)  # MAD A VERIF
                raise "LPSVD Not debugged!"

            elif (p_in["f2_lp_ext_algo"] == "lpsvd_stable"):
                dt2svd(get_si1_1d())
                svd2ar(1)
                ar2rt(1)
                rtreflect(1)
                rt2ar(1)
                ar2dt(p_in["lp_ext_size"],1)
                raise "LPSVD Not debugged!"

            audittrail( audit, "text", "FID extension before FT",
                        "algorithm", p_in["f2_lp_ext_algo"],
                        "LP order", p_in["f2_lp_ext_order"])

            if (key_is_not_false(p_in,"f2_lp_ext_apod")):
                sin(0, "f2")
                audittrail( audit, "text", "Apodisation after extension", "function", "sin(0,f2)")

    # %action% apodize
    #   standard apodisation
    # %param% f2_apodize boolean / default 1
    # %param% f2_apodisation string / default "sin(0.1)"
    #   the string describing the apodisation to apply
    #    "apodisation" should be a suite of window apodisation commands
    #     ie : 'sin(0.5)' or 'exbroad(0.2); sin(0.5)'
    #     the following predefined functions are implemnted
    #     sin (sine bell) sqsin (squared sine bell) expbroad (exponential) gaussbroad (gaussian)
    #     gaussenh (gaussian enhancement)
    if (key_is_true(p_in,"f2_apodize") and (p_in.raw("f2_apodisation") != "none")):
        if debug: writec("TTT.gs2")
        
        #set_task("F2 apodisation")
        if debug: print("F2 Apodisation "+p_in.raw("f2_apodisation"))
        apodise( p_in.raw("f2_apodisation"), "F2")
        audittrail( audit,  "text", "FID apodisation before FT",
                    "function", p_in.raw("f2_apodisation"))
        

    # %action% Fourier_transform
    #  performs the Fourier transform
    # %param% f2_fourier_transform boolean / default 1
    # %param% f2_ask_for_ftsize boolean / default 0
    #  if set, the final size is determined by f2_ft_size, automatically determined otherwise
    #   size for FT is determined by  size or lp_ext_size or trunc_size
    #   it will be extended to the next 2^n
    # %param% f2_ft_size integer  / default 2*power2(get_si2_2d())
    #   the size for FT if f2_ask_for_ftsize is true
    # %param_cond% (f2_ft_size == power2(n)
    # %param% f2_ft_type enum none ft rft ft_seq ft_sim / default ft_sim
    #    the Fourier transform algorithm, depends on spectrometer and acquisition scheme
    # %return% f2_size_for_ft
    #    contains the data size right after FT (will usually be different from final size)
    if (key_is_true(p_in,"f2_fourier_transform")):
        set_task("F2 Fourier Transform")
        if (key_is_true(p_in,"f2_ask_for_ftsize")):
            final = int(p_in["f2_ft_size"])
        else:
            # try   size >= lp_ext_size and >= 2*initial size and == 2^n
            if (p_in["f2_lp_extend"]):
                final = (2*power2(get_si2_2d() - 2))
            else:
                final = (4*power2(get_si2_2d() - 2))
        initial = get_si2_2d()
        chsize( get_si1_2d(), final )
        if (p_in["f2_ft_type"] != "none"):
            if (p_in["f2_ft_type"] == "ft"):
                # if ($itype_2d == 0 | $itype_2d == 2) itype (%+1)
                ft("f2")
            elif (p_in["f2_ft_type"] == "rft"):
                # if ($itype_2d == 1 | $itype_2d == 3) itype (%-1)
                rft("f2")
            elif (p_in["f2_ft_type"] == "ft_sim"):
                ft_sim()
            elif (p_in["f2_ft_type"] == "ft_seq"):
                ft_seq()
            else:
                # executer la macro pointee par $p_in[f2_ft_type]
                pass

        f_out["f2_size_for_ft"] = final
        audittrail( audit,  "text", "F2 Fourier transformation,",
                    "FT type", p_in["f2_ft_type"],
                    "initial size ", initial,
                    "final size ", final )

    # %param_exclusive% causal [ causalize causal_corr ]
    if (key_is_true(p_in,"causalize") and key_is_true(p_in,"causal_corr")):
        raise """
        causalize and causal_corr are incompatible
        causality must be corrected only once"""

    # %action% causal_corr
    #   performs causal correction of the spectrum if not done on the FID
    # %param% causal_corr boolean / default 1
    # %use% AXISF2_ZEROTIMEPOSITION
    # %return% causalize
    #   contains the correction applied by causalize
    if (key_is_true(p_in,"causal_corr")):
        set_task("F2 causal correction")
        delay = float(f_in["axisf2_zerotimeposition"])
        phase( 0.0, -360.0*delay, "f2")
        f_out["causalize"] = -360.0*delay
        audittrail( audit,  "text", "Causal correction",
                    "1st order phase correction", f_out["causalize"]) 

    # %action% f2_reverse
    #   reverses the spectral axis after FT
    # %param% f2_reverse boolean / default 0
    #   depends on spectrometer
    if(key_is_true(p_in,"f2_reverse")):
        reverse("f2")
        if (get_itype_2d()%2 == 1):
            invf("f2")       # if complex, inverse imaginary part
        audittrail( audit, "text", "Reverse F2 spectral axis")

    write_file_2d( audit, fileout )

#---------------------------------------------------------------------------
def pre_ft_2df1( audit, filein, fileout, p_in, f_in, f_out ):

    """
This macro realizes the pre FT operation on each F1 FID of the 2D

it implements FID massaging before FT:
- left_shift       : drops first points of the FID 
- right_shift      : adds empty points on the beginning of the FID
- back_extend      : reconstructs missing points in the beginning of the FID by LP analysis

%F1-input-domain%  time
%F1-output-domain% time
%dimensionality%   2

%author% Marc-Andre Delsuc
%version% 6.0

    """

    read_file_2d( audit, filein)

    # %action% f1_left_shift
    #   shifts the FID to the left by dropping data points, probably always useless !
    # %param% f1_left_shift boolean  / default 0
    # %param% f1_left_shift_size integer  / default 2
    #   number of points to shift
    # %param_cond%  [ (f1_left_shift_size > 0) and (f1_left_shift_size < get_si1_2d()) ]

    if(key_is_true(p_in,"f1_left_shift")):
            set_task("F1 left shift")
            left_shift(p_in["f1_left_shift_size"],'F1')
            audittrail( audit, "text", "F1 left shift",
                       "number of points", p_in["f1_left_shift_size"])

    # %action% f1_right_shift
    #   shifts the FID to the right by adding zero data points
    # %param% f1_right_shift boolean  / default 0
    # %param% f1_right_shift_size integer  / default 2
    #   number of points to shift
    # %param_cond%  [ (f1_right_shift_size > 0) ]
    if(key_is_true(p_in,"f1_right_shift")):
            set_task("F1 right shift")
            right_shift(p_in["f1_right_shift_size"],'F1')
            audittrail( audit, "text", "F1 right shift",
                       "number of points", p_in["f1_right_shift_size"])

    # %action-group% exclusive [ f1_right_shift f1_back_extend ]

    # %action% f1_back_extend
    #   extends the FID backward for reconstruction of initial missing points
    # %param% f1_b_extend boolean  / default 0
    # %param% f1_b_extend_size integer > 0  / default 2*get_si1_2d()
    #   number of missing points to reconstruct
    # %param_cond% [ f1_b_extend_size > 0 ]
    # %param% f1_b_extend_algo enum burg / default burg
    #   algorithm used
    # %return% f1_b_extend_order
    #   the order used during back_extend
    # %return% f1_b_extend_size
    #   the number of points added by back_extend
    if(key_is_true(p_in,"f1_b_extend")):
            set_task("F1 backward extension")
            n = p_in["f1_b_extend_size"]
            order(10)
            if (p_in["f1_b_extend_algo"] == "burg" ):
                burg2d_back("f1", get_si1_2d()+n )

            f_out["f1_b_extend_order"] = get_order()
            f_out["f1_b_extend_size"] = n
            audittrail( audit, "text", "F1 back extension",
                       "number of points", p_in["f1_b_extend_size"],
                       "algorithm", p_in["f1_b_extend_algo"])

    write_file_2d( audit, fileout)

#---------------------------------------------------------------------------
def pre_ft_2df2( audit, filein, fileout, p_in, f_in, f_out ):
    """
This macro realizes the pre FT operation on each F2 FID of the 2D

it implements FID massaging before FT:
- dc_offset        : corrects for constant offset in FID
- causalize        : changes DSP processed FID (Bruker) to causal FID's by Hilbert transform
- flatten_solvent  : removes solvent signal by FID analysis
- left_shift       : drops first points of the FID 
- right_shift      : adds empty points on the beginning of the FID
- back_extend      : reconstructs missing points in the beginning of the FID by LP analysis

%F2-input-domain%  time
%F2-output-domain% time
%dimensionality%   2

%author% Marc-Andre Delsuc
%version% 6.0

    """
    read_file_2d( audit, filein)

    # %action% dc_offset
    #   removes a constant level in the fid, estimated from the last points
    # %param% dc_offset boolean  / default 0
    # %param%  fid_noise_f2zone float / default 0.2
    #   constant offset of the FID, estimated on the last points of the dataset
    if key_is_true(p_in,"dc_offset"):
        set_task("DC-offset correction")
        f_out["pre_ft_mean_shift"]=dc_offset(p_in["fid_noise_f2zone"])
        audittrail( audit, "text", "mean vertical shift of F2 FIDs","value",f_out["pre_ft_mean_shift"])

    # %param_exclusive% causal [ causalize causal_corr ]
    if (key_is_true(p_in,"causalize") and key_is_true(p_in,"causal_corr")):
        raise """
        causalize and causal_corr are incompatible
        causality must be corrected only once"""

    # %action% causalize
    #   removes the non causal header found in certain (Bruker) implementation of the digital filter.
    # %param% causalize boolean / default 0
    # %use% axisf2_zerotimeposition
    # %return% causalize
    #   contains the correction applied by causalize
    if (key_is_true(p_in,"causalize")):
        set_task("F2 causal correction")
        causalize( f_in["axisf2_zerotimeposition"])
        f_out["causalize"] = -360*f_in["axisf2_zerotimeposition"]
        audittrail( audit, "text", "F2 causalisation correction",
                    "1st order phase correction", f_out["causalize"])

    # %action% flatten_solvent
    #   removes water from fit of the FID in the time domain
    # %param% flatten_solvent boolean  / default 0
    # %param% flat_solv_mode enum polynomial moving_average polynomial+moving_average moving_average+polynomial / default polynomial
    #   algorithm used
    if (key_is_true(p_in,"flatten_solvent")):
            set_task("Solvent Flattening")
            flat_solvent( p_in["flat_solv_mode"] )
            audittrail( audit, "text", "apply solvent flattening",
                       "algorithm", p_in["flat_solv_mode"])

    # %action% f2_left_shift
    #   shifts the FID to the left by dropping data points
    # %param% f2_left_shift boolean  / default 0
    # %param% f2_left_shift_size integer  / default 2
    #   number of points to shift
    # %param_cond%  [ (f2_left_shift_size > 0) and (f2_left_shift_size < get_si2_2d()) ]
    if(key_is_true(p_in,"f2_left_shift")):
            set_task("F2 left shift")
            left_shift(p_in["f2_left_shift_size"],'F2')
            audittrail( audit, "text", "F2 left shift",
                       "number of points", p_in["f2_left_shift_size"])

    # %action% f2_right_shift
    #   shifts the FID to the right by adding zero data points
    # %param% f2_right_shift boolean  / default 0
    # %param% f2_right_shift_size integer  / default 2
    #   number of points to shift
    # %param_cond%  [ (f2_right_shift_size > 0) ]
    if (key_is_true(p_in,"f2_right_shift")):
            set_task("F2 right shift")
            right_shift(p_in["f2_right_shift_size"],'F2')
            audittrail( audit, "text", "F2 right shift",
                       "number of points", p_in["f2_right_shift_size"])

    # %action-group% exclusive [ f2_right_shift f2_back_extend ]

    # %action% back_extend 
    #   extends the FID backward for reconstruction of initial missing points
    # %param% f2_b_extend boolean  / default 0
    # %param% f2_b_extend_size integer / default 2
    #   number of missing points to reconstruct
    # %param_cond% [ f2_b_extend_size > 0 ]
    # %param% f2_b_extend_algo enum burg / default burg
    #   algorithm used
    # %return% f2_b_extend_order
    #   the order used during back_extend
    # %return% f2_b_extend_size
    #    the number of points added by back_extend
    if (key_is_true(p_in,"f2_b_extend")):
            set_task("F2 backward extension")
            n = p_in["f2_b_extend_size"]
            order (min(30,n*2))
            if (p_in["f2_b_extend_algo"] == "burg" ):
                burg2d_back( "f2", get_si2_2d() + n)

            f_out["f2_b_extend_order"] = get_order()
            f_out["f2_b_extend_size"] = n
            audittrail( audit,  "text", "F2 back extension",
                        "number of points", p_in["f2_b_extend_size"],
                        "algorithm", p_in["f2_b_extend_algo"])

    write_file_2d( audit, fileout)

#---------------------------------------------------------------------------
def post_ft_2df1( audit, filein, fileout, p_in, f_in, f_out ):
    """
This macro realizes the post FT operation on each F1 spectrum of the 2D

it implements the massaging of the spectrum after FT :
- phase         : apply phase correction on the complex spectrum
- autophase     : compute and apply phase correction 
- invHibert     : apply inverse Hilbert transform
- calibration   : apply ppm calibration
- autocalibrate : autocalibration of the F1 axis based on the F2 calibration
- spectral_zone : extract a spectral zone

%F2-input-domain%  frequency
%F2-output-domain% frequency
%dimensionality%   2

%author% Marc-Andre Delsuc
%version% 6.0

    """
    read_file_2d( audit, filein)

    # we need complex data
    tocomplex("F1")

    
    # %param_exclusive% f1_phase exclusive [ modulus f1_phase f1_autophase ] 

    # %action% phase
    #   applies a F1 phase correction to the spectrum
    # %param% f1_phase boolean / default 0
    # %param% f1_phase_0 float / default 0.0
    #  global order phase correction
    # %param_cond% ( (f1_phase_0 >= -180 ) and (f1_phase_0 <= 180) )
    # %param% f1_phase_1 float / default 0.0
    #  1st order phase correction

    if (key_is_true(p_in,"f1_phase")):
            set_task("F1 phase correction")
            phase(p_in["f1_phase_0"], p_in["f1_phase_1"], "f1")
            audittrail( audit,  "text", "F1 phase correction applied",
                    "0th order", p_in["f1_phase_0"],
                    "1th order", p_in["f1_phase_1"])

    # %action% autophase
    #   automatic phase correction of the spectrum in F1
    #   phase is computed from a set 
    # %param% f1_autophase boolean / default 0
    # %param% f1_phase_algo enum apsl apmin / default apmin
    #  algorithm used
    # %return% f1_autophase_0
    #   the value of the global order phase correction applied by autophase 
    # %return% f1_autophase_1
    #   the value of the 1st order phase correction applied by autophase 
    if (key_is_true(p_in,"f1_autophase")):
            set_task("F1 automatic phase correction")
            col(1)
            dim(1)
            # read("col1.gifa4")
            if ( p_in["f1_phase_algo"] == "apsl" ):
                print("********* running APSL")
                apsl()
            if ( p_in["f1_phase_algo"] == "apmin" ):
                apmin()
            dim(2)
            phase( get_ph0(), get_ph1(), "f1")
            f_out["f1_autophase_0"] = get_ph0()
            f_out["f1_autophase_1"] = get_ph1()
            audittrail( audit, "text", "F1 automatic phase correction applied",
                       "algorithm", p_in["f1_phase_algo"],
                       "0th order", get_ph0(),
                       "1th order", get_ph1())

    # %param_exclusive% real [ modulus f1_invHilbert ]

    # %action% invHilbert
    #   apply an inverse Hilbert transform.
    #   This operation takes a complex spectrm, and generate a real spectrum in-place
    #   i.e. with twice the number of real points.
    # %param% f1_invHilbert boolean / default 0
    if (key_is_true(p_in,"f1_invhilbert")):
            set_task("F1 inverse Hilbert transform")
            invhilbert("F1")
            audittrail( audit, "text", "F1 inverse Hilbert transform applied")

    # %param_exclusive%  [ f1_calibration  f1_autocalibrate ]

    # %action% autocalibrate
    #   calibrates the F1 axis based on the F2 calibration, using IUPAC / IUPAB conventions
    # %param% f1_autocalibrate boolean / default 1
    # %param% autocalibrate_mode   enum standard biomolecule / default standard
    #   if autocalibrate_mode is standard, IUPAC rules are followed,
    #   if autocalibrate_mode is biomolecule then
    #           DSS is used instead of TMS for 2H and 13C
    #           NH3 is used instead of CH3NO2 for 15N
    #           (CH3O)3PO is used instead of H3PO4 for 31P
    try:
        f1_autocalibrate = p_in["f1_autocalibrate"]
    except:
        f1_autocalibrate = 0
    if(f1_autocalibrate):
        set_task("F1 automatic calibration")
        spin = autocalib(p_in["autocalibrate_mode"])
        audittrail( audit, "text", "F1 autocalibration  applied",
                "mode",p_in["autocalibrate_mode"],
                "spin",spin,
                "calibration", get_offset_1_2d()/get_freq_1_2d())

    # %action% calibration
    #   calibrate the ppm scale on the spectrum
    # %param% f1_calibration float / default 0
    #   the coordinate in ppm of the lower-most point in the spectrum
    if (f1_autocalibrate == 0 ):
        try:
            keytest = p_in["f1_calibration"]
        except:
            p_in["f1_calibration"] = 0
        set_task("F1 calibration")
        offset(p_in["f1_calibration"]*get_freq_1_2d(), get_offset_2_2d())
        audittrail( audit, "text", "F1 offset applied",
                   "calibration", get_offset_1_2d()/get_freq_1_2d())

    # %action% spectral_zone
    #   extract one spectral zone of the spectrum
    # %param% f1_spectral_zone boolean / default 0
    # %param% f1_spec_zone_left float / default 10.0
    #   the left border of the extract zone, in unit
    # %param% f1_spec_zone_left_unit enum ppm hz index / default ppm
    #   the unit in which spec_zone_left is given
    # %param% f1_spec_zone_right float / default 0.0
    #   the right border of the extract zone, in unit
    # %param% f1_spec_zone_right_unit enum ppm hz index / default ppm
    #   the unit in which spec_zone_right is given
    # %return% f1_spec_zone_left
    #  the left coordinate of the extracted spectral zone in index
    # %return% f1_spec_zone_right
    #  the right coordinate of the extracted spectral zone in index
    if key_is_true(p_in,"f1_spectral_zone"):
        set_task("F1 spectral zone extraction")
        [l,r] = spectral_zone(p_in["f1_spec_zone_left"], p_in["f1_spec_zone_right"], "F1", p_in["f1_spec_zone_left_unit"], p_in["f1_spec_zone_right_unit"])

        f_out["f1_spec_zone_left"] = l
        f_out["f1_spec_zone_right"] = r
        audittrail( audit, "text", "F1 spectral extraction",
                    "left point", f_out["f1_spec_zone_left"],
                    "right point", f_out["f1_spec_zone_right"],
                    "final size", get_si1_2d())
    else:
        f_out["f1_spec_zone_left"] = 1
        f_out["f1_spec_zone_right"] = get_si1_2d()

    write_file_2d( audit, fileout)

#---------------------------------------------------------------------------
def post_ft_2df2( audit, filein, fileout, p_in, f_in, f_out ):
    """
This macro realizes the post FT operation on each F2 spectrum of the 2D

it implements the massaging of the spectrum after FT :
- phase         : apply phase correction on the complex spectrum
- autophase     : compute and apply phase correction 
- invHibert     : apply inverse Hilbert transform
- calibration   : apply ppm calibration
- spectral_zone : extract a spectral zone

%F2-input-domain%  frequency
%F2-output-domain% frequency
%dimensionality%   2

%author% Marc-Andre Delsuc
%version% 6.0

    """
    read_file_2d( audit, filein)

    tocomplex("F2")

    # %param_exclusive% f2_phase exclusive [ modulus f2_phase f2_autophase ] 

    # %action% phase
    #   applies a F2 phase correction to the spectrum
    # %param% f2_phase boolean / default 0
    # %param% f2_phase_0 float / default 0.0
    #  global order phase correction
    # %param_cond% ( (f2_phase_0 >= -180 ) and (f2_phase_0 <= 180) )
    # %param% f2_phase_1 float / default 0.0
    #  1st order phase correction

    if (key_is_true(p_in,"f2_phase")):
            set_task("F2 phase correction")
            phase(p_in["f2_phase_0"], p_in["f2_phase_1"], "f2")
            audittrail( audit,  "text", "F2 phase correction applied",
                    "0th order", p_in["f2_phase_0"],
                    "1th order", p_in["f2_phase_1"])

    # %action% autophase
    #   automatic phase correction of the spectrum in F2
    #   phase is computed on the first experiment
    # %param% f2_autophase boolean / default 0
    # %param% f2_phase_algo enum apsl apmin / default apmin
    #  algorithm used
    # %return% f2_autophase_0
    #   the value of the global order phase correction applied by autophase 
    # %return% f2_autophase_1
    #   the value of the 1st order phase correction applied by autophase 
    if (key_is_true(p_in,"f2_autophase")):
            set_task("F2 automatic phase correction")
            row(1)
            dim(1)
            if ( p_in["f2_phase_algo"] == "apsl" ):
                apsl()
            if ( p_in["f2_phase_algo"] == "apmin" ):
                apmin()
            dim(2)
            phase( get_ph0(), get_ph1(), "F2")
            f_out["f2_autophase_0"] = get_ph0()
            f_out["f2_autophase_1"] = get_ph1()
            audittrail( audit, "text", "F2 automatic phase correction applied",
                       "algorithm", p_in["f2_phase_algo"],
                       "0th order", get_ph0(),
                       "1th order", get_ph1())

    # %param_exclusive% real [ modulus f2_invHilbert ]

    # %action% invHilbert
    #   apply an inverse Hilbert transform.
    #   This operation takes a complex spectrum, and generates a real spectrum in-place
    #   i.e. with twice the number of real points.
    # %param% f2_invHilbert boolean / default 0
    if (key_is_true(p_in,"f2_invhilbert")):
            set_task("F2 inverse hilbert transform")
            invhilbert("F2")
            audittrail( audit, "text", "F2 inverse Hilbert transform applied")

    # %action% calibration
    #   calibrate the ppm scale on the spectrum
    # %param% f2_calibration float / default 0
    #   the coordinate in ppm of the right-most point in the spectrum
    try:
        keytest = p_in["f2_calibration"]
    except:
        p_in["f2_calibration"]=0
    set_task("F2 calibration")
    offset( get_offset_1_2d(), p_in["f2_calibration"]*get_freq_2_2d())
    audittrail( audit, "text", "F2 offset applied",
               "calibration", get_offset_2_2d()/get_freq_2_2d())

    # %action% spectral_zone
    #   extract one spectral zone of the spectrum
    # %param% f2_spectral_zone boolean / default 0
    # %param% f2_spec_zone_left float / default 10.0
    #   the left border of the extract zone, in unit
    # %param% f2_spec_zone_left_unit enum ppm hz index / default ppm
    #   the unit in which spec_zone_left is given
    # %param% f2_spec_zone_right float / default 0.0
    #   the right border of the extract zone, in unit
    # %param% f2_spec_zone_right_unit enum ppm hz index / default ppm
    #   the unit in which spec_zone_right is given
    # %return% f2_spec_zone_left
    #  the left coordinate of the extracted spectral zone in index
    # %return% f2_spec_zone_right
    #  the right coordinate of the extracted spectral zone in index
    if key_is_true(p_in,"f2_spectral_zone"):
        set_task("F2 spectral zone extraction")
        (l,r) = spectral_zone(p_in["f2_spec_zone_left"], p_in["f2_spec_zone_right"], "F2", p_in["f2_spec_zone_left_unit"], p_in["f2_spec_zone_right_unit"])

        f_out["f2_spec_zone_left"] = l
        f_out["f2_spec_zone_right"] = r
        audittrail( audit, "text", "F2 spectral extraction",
                    "left point", f_out["f2_spec_zone_left"],
                    "right point", f_out["f2_spec_zone_right"],
                    "final size", get_si2_2d())
    else:
        f_out["f2_spec_zone_left"] = 1
        f_out["f2_spec_zone_right"] = get_si2_2d()

    write_file_2d( audit, fileout)

#---------------------------------------------------------------------------
def pp_2d( audit, filein, filepeak, p_in, f_in, f_out ):
    """
This macro realizes the peak picking of a 2D dataset

it implements
- spec_noise        : estimates noise
- prefilter               : apply a filtering to the data-set to remove features below a given separation in Hz
- restrict          : defines the zone in whch the peak-picking is to be done
- peakpick                : do the peak picking
- postfilter              : apply a post filtering on the peak list

%F1-input-domain%  frequency
%F2-input-domain% frequency
%dimensionality%   2

%author% Marc-Andre Delsuc
%version% 6.0

    """
    read_file_2d( audit, filein)

    if ( get_itype_2d() != 0 ):
        toreal("F12")
        audittrail( audit, "text", "complex data-set, removing the imaginary part")

    size_ini1 = get_si1_2d()
    size_ini2 = get_si2_2d()
    put("data")

    # %action% spec_noise
    #   evaluate noise, estimated by finding an empty zone
    # %param% spec_noise boolean / default 0
    # %param% spec_noise_n integer / default 30
    #    number of different zones where noise is evaluated
    # %return% NoiseInFrequencyDomain
    try:
        keytest = p_in["spec_noise"]
    except:
        p_in["spec_noise"] = 0

    #  estimate of the noise in the spectrum
    set_task("Noise evaluation")
    if (p_in["spec_noise"]):
        spec_noise( p_in["spec_noise_n"] )
        f_out["noise"] = get_noise()
    else:
        noise( f_in["noise"]+0.001 )
        f_out["noise"] = get_noise()

    # %action% prefilter
    #    smooths the spectrum before peak-picking, modification is not stored permanently
    # %param% prefilter boolean / default 0
    # %param% f1_prefilter_value float / default 10.0
    #    pre-smoothing parameter for F1 axis in Hz
    # %param% f2_prefilter_value float / default 10.0
    #    pre-smoothing parameter for F2 axis in Hz
    if (key_is_true(p_in,"prefilter")):
        set_task("Spectrum pre-filtering")
        # zerofill to next power of 2
        chsize (2*power2(get_si1_2d()-1), 2*power2(get_si2_2d()-1))
        size_ft1 = get_si1_2d()
        size_ft2 = get_si2_2d()
        # fourier transform
        iftbis("f12")
        if  (p_in["peak_sign"] == "both"):
            contrast = 0.7
        else:
            contrast = 1.4  # peut-etre un peu fort

        gaussenh( p_in["f1_prefilter_value"], contrast, "f1")
        gaussenh( p_in["f2_prefilter_value"], contrast, "f2")
        # reduce trucation artifacts
        sqsin(0, "f12")
        #chsize (max(2k,$size_ft))
        # restore smoothed data
        ftbis("f12")
        chsize( size_ini1, size_ini2 )
        audittrail( audit,  "text", "Pre Filter applied to data-set",
                   "smoothing value F1 in Hz (using gaussenh)", p_in["f1_prefilter_value"],
                   "smoothing value F2 in Hz (using gaussenh)", p_in["f2_prefilter_value"]) 

    # %action% restrict
    #     restricts the peak picking to a certain spectral zone 
    # %param% restrict boolean / default 0
    # %param% f1_restrict_left float / default 10.0
    #   the left border of the F1 restrict zone, in unit
    # %param% f1_restrict_left_unit enum ppm hz index / default ppm
    #   the unit in which F1 restrict_left is given
    # %param% f1_restrict_right float / default 0.0
    #   the right border of the F1 restrict zone, in unit
    # %param% f1_restrict_right_unit enum ppm hz index / default ppm
    #   the unit in which F1 restrict_right is given
    # %param_cond% (f1_restrict_left{index} < f1_restrict_right{index})
    # %return% f1_restrict_left
    #     the left coordinate of the F1 restrict zone in index
    # %return% f1_restrict_right
    #     the right coordinate of the F1 restrict zone in index
    # %param% f2_restrict_left float / default 10.0
    #   the left border of the F2 restrict zone, in unit
    # %param% f2_restrict_left_unit enum ppm hz index / default ppm
    #   the unit in which F2 restrict_left is given
    # %param% f2_restrict_right float / default 0.0
    #   the right border of the F2 restrict zone, in unit
    # %param% f2_restrict_right_unit enum ppm hz index / default ppm
    #   the unit in which f2_restrict_right is given
    # %param_cond% (f2_restrict_left{index} < f2_restrict_right{index})
    # %return% f2_restrict_left
    #     the left coordinate of the F2 restrict zone in index
    # %return% f2_restrict_right
    #     the right coordinate of the F2 restrict zone in index
    if (key_is_true(p_in,"restrict")):
        set_task("Defining peak-picking zone")
        if (p_in["f1_restrict_left_unit"] == 'ppm'):
            l1 = (ptoi(p_in["f1_restrict_left"],2,1))
        elif (p_in["f1_restrict_left_unit"] == 'hz'):
            l1 = (htoi(p_in["f1_restrict_left"],2,1))
        elif (p_in["f1_restrict_left_unit"] == 'index'):
            l1 = p_in["f1_restrict_left"]
        else:
            raise "Error with f1_restrict_left_unit"

        if (p_in["f1_restrict_right_unit"] == 'ppm'):
            r1 = (ptoi(p_in["f1_restrict_right"],2,1))
        elif (p_in["f1_restrict_right_unit"] == 'hz'):
            r1 = (htoi(p_in["f1_restrict_right"],2,1))
        elif (p_in["f1_restrict_right_unit"] == 'index'):
            r1 = p_in["f1_restrict_right"]
        else:
            raise "Error with f1_restrict_right_unit"

        l1 = (max(1,round(l1)))
        r1 = (min(get_si1_2d(),round(r1)))
        if (l1 > r1):
            raise "Wrong f1 restrict zone coordinates"

        if ((r1-l1) < 8):
            raise "f1 restrict zone too small"

        if (p_in["f2_restrict_left_unit"] == 'ppm'):
            l2 = (ptoi(p_in["f2_restrict_left"],2,2))
        elif (p_in["f2_restrict_left_unit"] == 'hz'):
            l2 = (htoi(p_in["f2_restrict_left"],2,2))
        elif (p_in["f2_restrict_left_unit"] == 'index'):
            l2 = p_in["f2_restrict_left"]
        else:
            raise "Error with f2_restrict_left_unit"

        if (p_in["f2_restrict_right_unit"] == 'ppm'):
            r2 = (ptoi(p_in["f2_restrict_right"],2,2))
        elif (p_in["f2_restrict_right_unit"] == 'hz'):
            r2 = (htoi(p_in["f2_restrict_right"],2,2))
        elif (p_in["f2_restrict_right_unit"] == 'index'):
            r2 = p_in["f2_restrict_right"]
        else:
            raise "Error with f2_restrict_right_unit"

        l2 = (max(1,round(l2)))
        r2 = (min(get_si2_2d(),round(r2)))
        if (l2 > r2):
            raise "Wrong f2 restrict zone coordinates"

        if ((r2-l2) < 8):
            raise "f2 restrict zone too small"


        zoom(1, l1, l2, r1, r2)
        f_out["f1_restrict_left"] = l1
        f_out["f1_restrict_right"] = r1
        f_out["f2_restrict_left"] = l2
        f_out["f2_restrict_right"] = r2
        audittrail( audit, "text", "spectral extraction",
                    "left point F1xF2",  f_out["f1_restrict_left"],  " ", f_out["f2_restrict_left"],
                    "right point F1xF2", f_out["f1_restrict_right"], " ", f_out["f2_restrict_right"])
    else:
        zoom(0)
        f_out["f1_restrict_left"] = 1
        f_out["f1_restrict_right"] = get_si1_2d()
        f_out["f2_restrict_left"] = 1
        f_out["f2_restrict_right"] = get_si2_2d()

    # %action% peakpick
    #     do the peak picking, by detecting local extrema 
    # %param% peakpick boolean / default 1
    # %param% noise_thresh float / default 10.0
    #   minimum peak/noise ratio for peak detection
    # %param% ratio_thresh float / default 30.0
    #    maximum tallest peak / peak ratio for peak detection
    # %param% peak_sign enum  positive negative both / default positive
    #   sign of the peaks to detect
    # %return% nb_detect_peaks
    #   the number of detected peaks
    # %return% low_limit
    #   low limit for peak detection
    if (key_is_true(p_in,"peakpick")):
            set_task("Peak Picking")
            com_max()
            mn = (max( geta_max(1)/p_in["ratio_thresh"],p_in["noise_thresh"]*get_noise()))
            minimax(mn, geta_max(1)+1)
            f_out["low_limit"] = mn
            if (p_in["peak_sign"] == "negative"):
                mult(-1.0)
            elif (p_in["peak_sign"] == "both"):
                com_abs()

            pkclear()
            peak(0)
            f_out["nb_detect_peaks"] = get_npk2d()
            get("data")
            pkreset()
            audittrail( audit, "text", "Peak Picking applied",
                   "low limit for peak detection", mn,
                   "number of peak detected", get_npk2d())

    # %action% aggregate
    #   sorts peak list to aggregate peaks close from each other
    # %param% aggregate boolean / default 0
    # %param% aggregate_value float / default 10.0
    #    distance parameter in Hz
    if (key_is_true(p_in,"aggregate")):
            audittrail( audit, "text", "aggregate - Reste a faire")

    # %action% symetrize_peak
    #   symetrizes the peak list to  peaks relative to the diagonal
    # %param% symetrize_peak boolean / default 0
    # %param% symetrize_peak_tol float / default 10.0
    #    tolerance for symetrisation, in Hz
    # %param% symetrize_peak_algo enum remove add / default add
    #   correction to apply,
    #   either add missing peaks, are remove non symetrical peaks
    if (key_is_true(p_in,"symetrize_peak")):
            set_task("Peak symetrisation")
            if (p_in["symetrize_peak_algo"] == "remove"):
                pksym( 0, p_in["symetrize_peak_tol"])
            elif (p_in["symetrize_peak_algo"] == "add"):
                pksym( 1, p_in["symetrize_peak_tol"])

            audittrail( audit,  "text", "Peak symetrisation",
                        "Number of peaks after correction", get_npk2d())

    # ecriture du fichier de peak
    pkwrite_p(filepeak)

    audittrail( audit, "text", "Peak Picking file stored",
               "filename", filepeak)

#---------------------------------------------------------------------------
def maxent_2d( audit, filein, fileout, p_in, f_in, f_out ):
    """
This macro realizes the MaxEnt spectral analysis of a 2D FID

 it implements the following steps :
 - freq_massage      : do some preparation in the frequency domain
    - causal_corr    : perform causal correction
    - phase                      : apply phase correction on the complex spectrum before analysis
    - spectral_zone  : extract a spectral zone before analysis
    - reverse        : reverses spectral axis
 - truncate          : remove points at the end of the FID
 - preconvoluate     : apply a preconcolution function of the FID
 - partialsampling   : implements analysis on partially sampled FID
 - positive_negative : implements two channel processing for positive and negative lines
 - deconvoluate      : deconvoluate spectrum from a given time domain response
 - maxent            : apply Maximum Entropy analysis of the dataset

arguments :
    audit       the opened audit file
    filein      the file name of the input file, can be "memory"
    fileout     the file name of the input file, can be "memory"
    p_in        the dictionary containing all the processing parameters
    f_in        the dictionary containing all the details on the spectrum "filein"
    f_out       details on the spectrum "fileout" will be put in this dictionary

 %F1-input-domain%  time
 %F1-output-domain% frequency
 %F2-input-domain%  time
 %F2-output-domain% frequency
 %dimensionality%   2

 %author% Marc-Andre Delsuc
 %version% 1.0
    """
    read_file_2d( audit, filein)

    # actions in the frequency domain
    # %action% f2_freq_massage
    #  computes a temporary Fourier transform in F2
    # %param% f2_freq_massage boolean / default 1
    # %param% f2_ft_type enum none ft_seq ft_sim ft rft / default ft
    #    the Fourier transform algorithm, depends on spectrometer and acquisition scheme
    # %param% f2_reverse boolean / default 0
    #   reverses the spectral axis after FT
    # %param% causal_corr boolean / default 1
    #   performs causal correction of the spectrum if not done on the FID
    # %param_exclusive% causal [ causalize causal_corr ]
    # %param% f2_phase boolean / default 0
    # %param% f2_phase_0 float / default 0.0
    #    global order phase correction
    # %param% f2_phase_1 float / default 0.0
    #    1st order phase correction
    # %use% AXISF2_ZEROTIMEPOSITION
    # %return% causalize contains the correction applied by causalize
    # %param%
    if (key_is_true(p_in,"f2_freq_massage")):
        set_task("F2 frequency massage")

        # first ft
        final = (4*power2(get_si2_2d() - 2)) # one zerofilling is required !
        initial = get_si2_2d()
        chsize( get_si1_2d(), final )
        if (p_in["f2_ft_type"] != "none"):
            if (p_in["f2_ft_type"] == "ft"):
                # if ($itype_2d == 0 | $itype_2d == 2) itype (%+1)
                ft("f2")
            elif (p_in["f2_ft_type"] == "rft"):
                # if ($itype_2d == 1 | $itype_2d == 3) itype (%-1)
                rft("f2")
            elif (p_in["f2_ft_type"] == "ft_sim"):
                ft_sim()
            elif (p_in["f2_ft_type"] == "ft_seq"):
                ft_seq()
            else:
                # executer la bonne macro @($p_in[f2_ft_type])
                raise "FT Mode not implemented yet"

        f_out["f2_size_for_ft"] = final
        audittrail( audit,  "text", "F2 temporary Fourier transformation,",
                    "FT type", p_in["f2_ft_type"],
                    "initial size", initial,
                    "final size ", f_out["f2_size_for_ft"])
        if(key_is_true(p_in,"f2_reverse")):
            reverse("f2")
            audittrail( audit, "text", "Reverse F2 spectral axis")

        if (key_is_true(p_in,"causal_corr")):
            delay = float(f_in["axisf2_zerotimeposition"])
            phase( 0.0,  (-360.0*delay), "f2")
            initial=initial-2*int(delay)
            f_out["causalize"] = (-360.0*delay)
            audittrail( audit,  "text", "Causal correction",
                        "1st order phase correction", f_out["causalize"],
                        "sized reduced to",initial)

        if (key_is_true(p_in,"f2_phase")):
            phase( p_in["f2_phase_0"], p_in["f2_phase_1"], "F2")
            audittrail( audit,  "text", "F2 phase correction applied",
                        "0th order", p_in["f2_phase_0"],
                        "1th order", p_in["f2_phase_1"])

        # %action% spectral_zone
        #   extract one spectral zone of the spectrum
        # %param% f2_spectral_zone boolean / default 0
        # %param% f2_spec_zone_left float / default 10.0
        #   the left border of the extract zone, in unit
        # %param% f2_spec_zone_left_unit enum ppm hz index / default ppm
        #   the unit in which spec_zone_left is given
        # %param% f2_spec_zone_right float / default 0.0
        #   the right border of the extract zone, in unit
        # %param% f2_spec_zone_right_unit enum ppm hz index / default ppm
        #   the unit in which spec_zone_right is given
        # %return% f2_spec_zone_left
        #  the left coordinate of the extracted spectral zone in index
        # %return% f2_spec_zone_right
        #  the right coordinate of the extracted spectral zone in index
        if key_is_true(p_in,"f2_spectral_zone"):
            set_task("F2 spectral zone extraction")
            (l,r) = spectral_zone(p_in["f2_spec_zone_left"], p_in["f2_spec_zone_right"], "F2", p_in["f2_spec_zone_left_unit"], p_in["f2_spec_zone_right_unit"])

            f_out["f2_spec_zone_left"] = l
            f_out["f2_spec_zone_right"] = r
            audittrail( audit, "text", "F2 spectral extraction",
                        "left point", f_out["f2_spec_zone_left"],
                        "right point", f_out["f2_spec_zone_right"],
                        "final size", get_si2_2d())
            # now modify initial in consequence
            initial = (initial * get_si2_2d())/final    # zero filling is homotetically after inverse FT
            initial = 2*(initial/2)     # has to be even
        else:
            f_out["f2_spec_zone_left"] = 1
            f_out["f2_spec_zone_right"] = get_si2_2d()

        #then back to FID
        real("F2")
        iftbis("F2")
        chsize( get_si1_2d(), initial)
        f_out["f2_mock_fid_size"] = initial
        audittrail( audit,  "text", "Back to F2 time-domain" ,
                    "final size ", initial)

    # actions in the frequency domain
    # %action% f1_freq_massage
    #  computes a temporary Fourier transform in F1
    # %param% f1_freq_massage boolean / default 0
    # %param% f1_ft_type enum none ft_sh ft_tppi ft_sh_tppi ft_n_p ft_phase_modu ft rft / default ft_sh
    #    the Fourier transform algorithm, depends on spectrometer and acquisition scheme
    # %param% f1_reverse boolean / default 0
    #   reverses the spectral axis after FT
    # %param% f1_phase boolean / default 0
    # %param% f1_phase_0 float / default 0.0
    #    global order phase correction
    # %param% f1_phase_1 float / default 0.0
    #    1st order phase correction
    if (key_is_true(p_in,"f1_freq_massage")):
        set_task("F1 frequency massage")
        # first ft
        final = (4*power2(get_si1_2d() - 2)) # one zerofilling is required !
        initial = get_si1_2d()
        chsize( final, get_si2_2d())
        if (p_in["f1_ft_type"] != "none"):
            if (p_in["f1_ft_type"] == "ft"):
                #  if ($itype_2d == 0 | $itype_2d == 2) itype (%+1)
                ft("f1")
            elif (p_in["f1_ft_type"] == "rft"):
                #  if ($itype_2d == 1 | $itype_2d == 3) itype (%-1)
                rft("f1")
            elif (p_in["f1_ft_type"] == "ft_sh"):
                ft_sh("f1")
            elif (p_in["f1_ft_type"] == "ft_tppi"):
                ft_tppi("f1")
            elif (p_in["f1_ft_type"] == "ft_sh_tppi"):
                ft_sh_tppi("f1")
            elif (p_in["f1_ft_type"] == "ft_n_p"):
                ft_n_p("f1")
            elif (p_in["f1_ft_type"] == "ft_phase_modu"):
                ft_phase_modu("f1")
            else:
                # executer @($p_in["f1_ft_type"])
                raise "FT Mode not implemented yet"

        f_out["f1_size_for_ft"] = final
        audittrail( audit,  "text", "F1 temporary Fourier transformation,",
                    "FT type", p_in["f1_ft_type"],
                    "initial size", initial,
                    "final size ", f_out["f1_size_for_ft"])
        if(key_is_true(p_in,"f1_reverse")):
            reverse("f1")
            audittrail( audit, "text", "Reverse F1 spectral axis")

        if (key_is_true(p_in,"f1_phase")):
            phase( p_in["f1_phase_0"], p_in["f1_phase_1"], "F1")
            audittrail( audit,  "text", "F1 phase correction applied",
                        "0th order", p_in["f1_phase_0"],
                        "1th order", p_in["f1_phase_1"])

        # then back to FID
        real("F1")
        iftbis("F1")
        chsize( initial, get_si2_2d())
        f_out["f1_mock_fid_size"] = initial
        audittrail( audit,  "text", "Back to F1 time-domain",
                    "final size ", initial)

    # %action% F2 truncate
    #   removes last points of FID
    # %param% f2_truncate boolean / default 0
    # %param% f2_trunc_size integer / default get_si2_2d()
    #   new FID size after truncation
    # %param_cond% [ (f2_trunc_size > 0 ) and (f2_trunc_size <= get_si2_2d() ) ]
    if(key_is_true(p_in,"f2_truncate")):
        s = p_in["f2_trunc_size"]
        if (s <= 0):
            raise "Cannot truncate to zero"

        if (s > get_si2_2d()):
            raise "Error with truncation size (larger than actual size)"

        initial = get_si2_2d()
        chsize( get_si1_2d(), p_in["f2_trunc_size"])
        audittrail( audit,  "text", "FID truncation before MaxEnt",
                    "initial size", initial,
                    "final size", p_in["f2_trunc_size"])

    # %action% F1 truncate
    # %param% f1_truncate boolean / default 0
    # %param% f1_trunc_size integer > 0 <= get_si1_2d() / default get_si1_2d()
    #   new FID size after truncation
    if(key_is_true(p_in,"f1_truncate")):
        s = p_in["f1_trunc_size"]
        if (s <= 0):
            raise "Cannot truncate to zero"

        if (s > get_si1_2d()):
            raise "Error with truncation size (larger than actual size)"

        initial = get_si1_2d()
        chsize( p_in["f1_trunc_size"], get_si2_2d())
        audittrail( audit,  "text", "F1 truncation before MaxEnt",
                    "initial size", initial,
                    "final size", p_in["f1_trunc_size"])

    # %action% preconvoluate
    #   apply a preconvolution before analysis, this may help stability of the algorithm and enhance noise rejection
    # %param% preconvoluate boolean / default 0
    # %param% f1_exponential_preconv float / default 1.0
    #   exponential preconvolution along F1 axis
    # %param% f2_exponential_preconv float / default 1.0
    #   exponential preconvolution along F2 axis

    if (key_is_true(p_in,"preconvoluate")):
            em(p_in["f1_exponential_preconv"],p_in["f2_exponential_preconv"])
            audittrail( audit, "text", "FID preconvolution before MaxEnt analysis",
                       "F1 exponetial preconvolution", p_in["f1_exponential_preconv"],
                       "F2 exponetial preconvolution", p_in["f2_exponential_preconv"])

    # now create the data for MaxEnt
    put("data") # copy the FID to the data buffer
    noise( f_in["noise"])              # and get noise value
    window_mode(0)          # assuming no windowing
    dim(2)
    window_reset("F1")
    window_reset("F2")

    audittrail( audit,  "text", "Mock FID created",
                "fid size, F1", get_si1_2d(),
                "fid size, F2", get_si2_2d(),
                "noise value", get_noise())

    si1data = get_si1_2d()  # remember values
    si2data = get_si2_2d()

    # %action% partialsampling
    #   set-up for processing data partially sampled in the time domain
    # %param% partialsampling boolean / default 0
    # %param% partialsampling_mode enum F1 F2 F1+F2 F12 / default F1+F2
    #   defines the sampling mod used on the 2D plane
    #   F1 : sampling along the F1 spectral axis only, -- use f1_samplingfile definition
    #   F2 : sampling along one F2 spectral axis only, -- use f2_samplingfile definition
    #   F1+F2 : sampling along both spectral axis, defined independentely
    #   F12 : global sampling of the F1xF2 plane -- use F12_samplingfile definition
    # %param% f1_samplingfile string / default "F1vplist"
    #   the filename of the file containing the list of sample to use along the F1 axis
    # %param% f2_samplingfile string / default "F2vplist"
    #   the filename of the file containing the list of sample to use along the F2 axis
    #   see sampling.g for details
    # %param% f12_samplingfile string / default "F12vplist"
    #   the filename of the file containing the list of sample to use on the F1xF2plane

    if (key_is_true(p_in,"partialsampling")):
        window_mode(1)
        window_reset("F1")
        window_reset("F2")
        if (p_in["partialsampling_mode"] == "F12"):
            window_mode(2)
            sign2 = get_si2_2d()
            sign1 = sampling( p_in["f12_samplingfile"],"F12")
            partial = float(sign1) / (get_si1_2d()*get_si2_2d())
            audittrail( audit,  "text", "Set-up for partial sampling during MaxEnt analysis",
                        "F1xF2 sampling file", p_in["f12_samplingfile"],
                        "number of sampled points", sign1,
                        "mean density of sampling", repr(100 * partial) + " %")
        else:
            if (p_in["partialsampling_mode"] == "F1" or p_in["partialsampling_mode"] == "F1+F2"):
                sign1 = sampling(p_in["f1_samplingfile"],"F1")
                f1 = p_in["f1_samplingfile"]
            else:
                f1 = "not used"
                sign1 = get_si1_2d()

            if (p_in["partialsampling_mode"] == "F2" or p_in["partialsampling_mode"] == "F1+F2"):
                sign2 = sampling(p_in["f2_samplingfile"],"F2")
                f2 = p_in["f2_samplingfile"]
            else:
                f2 = "not used"
                sign2 = get_si2_2d()

            partial = float(sign1*sign2) / (get_si1_2d()*get_si2_2d())
            audittrail( audit, "text", "Set-up for partial sampling during MaxEnt analysis",
                   "F1 sampling file", f1,
                   "number of sampled points", sign1,
                   "F2 sampling file", f2,
                   "number of sampled points", sign2,
                   "mean density of sampling", repr(100 * partial)+" %")
    else:
        partial = 1.0

    # %action% positive_negative
    #   set-up for positive_negative analysis
    # %param% positive_negative boolean / default 0

    if (key_is_true(p_in,"positive_negative")):
        p_in["deconvoluate"]=1      # positive_negative implies deconvolution

    # %action% deconvoluate
    #   describes the deconvolution to apply during analysis,
    # %param% deconvoluate boolean / default 1
    #   this part is not finished yet, only suport for exponential deconvolution so far
    # %param% f1_exponential_deconv float / default 1.0
    #   exponential deconvolution along F1 axis
    # %param% f2_exponential_deconv float / default 1.0
    #   exponential deconvolution along F2 axis

    if (key_is_true(p_in,"deconvoluate")):
        if p_in["f2_exponential_deconv"]:
            lb2=p_in["f2_exponential_deconv"]
        else:
            lb2=0.0
        if p_in["f1_exponential_deconv"]:
            lb1=p_in["f1_exponential_deconv"]
        else:
            lb1=0.0
        setup_filter2D(si1data,si2data,(lb1,lb2),p_in["positive_negative"])
        
        audittrail( audit, "text", "deconvolution set-up for MaxEnt analysis",
                       "F1 exponential deconvolution", p_in["f1_exponential_deconv"],
                       "F2 exponential deconvolution", p_in["f2_exponential_deconv"],
                       "Pos/Neg two channel",p_in["positive_negative"])
    else:
        com_filter(0)
        nchannel(1)
        audittrail( audit, "text", "No deconvolution set-up for MaxEnt analysis")

    # %action% me_preset
    #   presets the MaxEnt parameters to default values
    # %param% me_preset_value enum 0 1 2 3 4 5 / default 3
    #   sets the parameter for a balance between speed (1) and quality (5), 0 is for fit

    me = maxent_state()
    if (key_is_true(p_in,"me_preset_value")):
        preset = p_in["me_preset_value"]
    else:
        preset = 3
    me.me_preset(preset)
    
    # %action% me_details
    #   if this flag is on, default parameters can be set.
    # %param% me_details boolean / default 0
    # %param% me_size integer 
    # %param% me_iteration integer
    # %param% me_int_iteration
    # %param% me_ncheck integer
    # %param% me_lambda_control integer
    # %param% me_lambda_speed 
    # %param% me_algo enum 0 1 2 / default 1
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
           me.set_algo(p_in["me_algo"])

    # set it in kernel
    me.copyto_kernel()
    audittrail(audit, "text","Parameters for MaxEnt processing","MaxEnt parameters",me.report())

    # %action% maxent
    #   apply MaxEnt analysis,
    # %param% maxent boolean / default 1
    # %param% me_preset boolean / default 1
    #   determines convergence parameters from a given preset or gives the details
    # %param% preset_value enum 1 2 3 4 5 / default 3
    #   preset for maxent# 1 is faster, 5 is more accurate
    # %param% iteration integer / default 100
    #   global number of iterations - set by preset
    # %param% control integer / default 10
    #   number of iterations between control step - set by preset
    # %param% lambsp float / default 5.0
    #   the value of the lambda controlling parameter LAMBSP in the Gifa algorithm
    # %param% lambcont enum 0 1 2 / default 1
    #   the value of the lambda controlling parameter LAMBCONT in the Gifa algorithm
    # %param% f1_me_size integer / default 2*power2(get_si1_2d())
    #   F1 size for MaxEnt reconstruction
    # %param% f2_me_size integer / default 2*power2(get_si2_2d())
    #   F2 size for MaxEnt reconstruction
    # %param% noise_weight float / default 1.0
    #   the noise used during the MaxEnt iteration is weigthed with this scalar,
    #   this permits to compensate for over- or under-estimate of the noise level
    # %param_cond%   [ noise_weight > 0.0 ]
    # %param_cond%   [ ispower2(f1_me_size) ]
    # %param_cond%   [ ispower2(f2_me_size) ]
    
    if (p_in["maxent"]):
            get("data")
            if (get_window_mode() != 0):
                apply("window")
            noise (get_noise()*p_in["noise_weight"]*partial)
            audittrail( audit, "text", "Parameters for MaxEnt Analysis",
                       "size of reconstruction in F1", p_in["f1_me_size"],
                       "size of reconstruction in F2", p_in["f2_me_size"],
                       "weighted noise", get_noise(),
                       "algorithm flags", get_algo(),
                       "lambda control", get_lambcont(),
                       "lambda multiplier", get_lambsp(),
                       "global iteration", get_iter(),
                       "control", get_ndisp(),
                       "line minimisation iteration", get_miniter())
            t = time.clock()
            sym_maxent(1)   # MAD 
            maxent(int(p_in["f1_me_size"]), int(p_in["f2_me_size"]))
#            iter(get_iter()+20)
#            maxentcont()
            t = time.clock()-t
            audittrail( audit, "text", "MaxEnt Analysis applied",
                       "number of iteration performed", get_iterdone(),
                       "final chi square", get_chi2(),
                       "final entropy", get_entropy(),
                       "final lambda", get_lambda(),
                       "final convergence criterium", 1-get_convergence(),
                        "MaxEnt processing time", str(t)+" sec")

    write_file_2d( audit, fileout)

#---------------------------------------------------------------------------
def post_2d( audit, filein, fileout, p_in, f_in, f_out, location="."):
    """
This macro realizes the post FT operation on a 2D spectrum

This macro implements the massaging of the spectrum after complete Spectral analysis :
- autophase2D   : computes and apply global phase correction
- modulus       : take the modulus of the complex spectrum
- baseline      : baseline correction
- shear         : apply the shear correction
- tilt          : apply the tilt correction
- symetrize     : apply symetrisation algorithms
- smoothing     : apply a smoothing filter
- median        : apply a media filter
- projection    : compute F1 and F2 projections
- select_state  : enforce the final state of the data
- spec_noise    : estimates noise in final spectrum

%F1-input-domain%  frequency
%F1-output-domain% frequency
%F2-input-domain%  frequency
%F2-output-domain% frequency

%dimensionality%   2

%author% Marc-Andre Delsuc
%version% 6.0
    """
    read_file_2d( audit, filein)

    # %param_exclusive% f1_phase exclusive [ modulus autophase2D] 
    #
    # %action% autophase_2d
    #   automatic phase correction of the spectrum in 2D
    #   phase is computed from a set of lines extracted fro the dataset
    # %param% autophase_2d boolean / default 1
    # %param% phase_2d_algo enum apsl2d apmin2d / default apmin2d
    #  algorithm used for the phase correction
    # %param% phase_2d_axis enum F1 F2 F12 / default F12 
    #  Axis along which the autophasing is to be applied.
    # %return% f1_autophase_0
    #   the value of the global order phase correction applied by autophase 
    # %return% f1_autophase_1
    #   the value of the 1st order phase correction applied by autophase 
    # %return% f2_autophase_0
    #   the value of the global order phase correction applied by autophase 
    # %return% f2_autophase_1
    #   the value of the 1st order phase correction applied by autophase 
    if (key_is_true(p_in,"autophase_2d")):
        if(key_is_true(p_in,"modulus")):
            pass
        else:
            set_task("Automatic 2D phase correction")
            tocomplex("F12") # insure we are cmplx
            if ( p_in["phase_2d_algo"] == "apsl2d" ):
                (p10, p11, p20, p21) = apsl2d(p_in["phase_2d_axis"])
            if ( p_in["phase_2d_algo"] == "apmin2d" ):
                (p10, p11, p20, p21) = apmin2d(p_in["phase_2d_axis"])
#            print p10,p11,p20,p21
            f_out["f1_autophase_0"] = p10
            f_out["f1_autophase_1"] = p11
            f_out["f2_autophase_0"] = p20
            f_out["f2_autophase_1"] = p21
            audittrail( audit, "text", "automatic 2D phase correction applied",
                           "algorithm", p_in["phase_2d_algo"],
                           "F1 0th order", p10,
                           "F1 1th order", p11,
                           "F2 0th order", p20,
                           "F2 1th order", p21)

    # %action% modulus
    #   modulus of the spectrum
    # %param% modulus boolean / default 0
    if(key_is_true(p_in,"modulus")):
        if (get_itype_2d() == 1):
            set_task("F2 modulus applied")
            modulus()
            audittrail( audit,  "text", "F2 modulus applied")
        elif (get_itype_2d() == 2):
            set_task("F1 modulus applied")
            osz = get_si2_2d()
            # flop requires a 2^n buffer !
            chsize( get_si1_2d(), 2*power2(get_si2_2d()-2))
            flop()
            chsize( get_si1_2d(), 2*osz)
            modulus()
            audittrail( audit, "text", "F1 modulus applied")
        elif (get_itype_2d() == 3):
            set_task("Hypercomplex modulus applied")
            modulus()
            audittrail( audit, "text", "Hypercomplex modulus applied")
        else:
            set_task("Taking absolute value")
            absolute()
            audittrail( audit, "text", "Taking Absolute Value")

    # %action% throw_imaginary
    #     at this stage, data should be real
    # %param% throw_imaginary boolean / default 1
    #   throw_imaginary set to false will keep hypercomplex data
    #   Note that most of the following procesing will fail
    if (key_is_not_false(p_in,"throw_imaginary")):
        toreal( "f12" )

    # %action% baseline correction
    #   baseline correction of the spectrum
    # %param% f2_baseline boolean / default 1
    # %param% f2_bcorr_algo enum  offset linear spline quest polynomial moving_average / default moving_average
    #    offset : removes a automatically determined offset
    #    linear - spline : uses the list of user determined pivot points to define the baseline,
    #       : fit with a straight line or a spline
    #    quest : reconstruct the beginning of the FID using Linear Prediction technics,
    #        should usually be followed by offset
    #    polynomial - moving_average : uses statistics to separate signal from baseline,
    #       : apply a polynomial or moveing average correction.
    #    in more than one term is given (blank separated string!), the correction are applied in sequence.
    # %param% spec_noise_n integer / default 10
    #   used by the offset algorithm to determine the offset to correct
    # %param% f2_bcorr_pivots Integerlist / default (int(0.02*get_si2_2d()),int(0.98*get_si2_2d()))
    #   pivot points used by the linear or spline algorithm
    # %param_cond% (f2_bcorr_pivots > 0 and bcorr_pivots <= get_si2_2d())
    # %param% f2_bcorr_radius integer / default 1
    #   radius around pivot points used by the linear or spline algorithm
    # %param% f2_bcorr_order integer  / default 10
    # %param_cond% (f2_bcorr_order > 0 and f2_bcorr_order < get_si2_2d() and f2_bcorr_order == 2*n)
    #   order (number of points corrected in the time domain) of the quest algorithm
    if(key_is_true(p_in,"f2_baseline")):
        set_task("F2 baseline correction")
        toreal("F2")
        pivots = "not used" 
        bcorder = "not used"
        bc = p_in["f2_bcorr_algo"]
        for bcalgo in bc.split():
            if (bcalgo == "offset"):
                bcorr_offset( p_in["spec_noise_n"], "f2")
            elif (bcalgo == "linear"):
                pivots=p_in["f2_bcorr_pivots"]
                bcorr(1,p_in["f2_bcorr_radius"],"F2",pivots)
            elif (bcalgo == "spline"):
                pivots=p_in["f2_bcorr_pivots"]
                bcorr(2,p_in["f2_bcorr_radius"],"F2",pivots)
            elif (bcalgo == "quest"):
                bcorder = p_in["f2_bcorr_order"]        
                bcorr_quest( bcorder, "f2")
            elif (bcalgo == "polynomial"):  # bcorrpx determines which algo is used
                if (get_debug()):
                    print("Baseline correction not available in DEBUG mode")
                else:
                    bcorrp0()
                    bcorr(3, "f2")
            elif (bcalgo == "moving_average"):
                if (get_debug()):
                    print("Baseline correction not available in DEBUG mode")
                else:
                    bcorrp1()
                    bcorr(3, "f2")
                    bcorrp0()
            else:
                raise "Error with baseline correction algorithm"
        audittrail( audit, "text", ("F2 Baseline correction -" + bc),
                         "algorithm", bcalgo,
                         "pivots", pivots,
                         "order", bcorder)

    # %action% baseline correction
    #   baseline correction of the spectrum
    # %param% f1_baseline boolean / default 1
    # %param% f1_bcorr_algo enum  offset linear spline quest polynomial moving_average / default moving_average
    #    offset : removes a automatically determined offset
    #    linear - spline : uses the list of user determined pivot points to define the baseline,
    #       : fit with a straight line or a spline
    #    quest : reconstruct the beginning of the FID using Linear Prediction technics,
    #        should usually be followed by offset
    #    polynomial - moving_average : uses statistics to separate signal from baseline,
    #       : apply a polynomial or moveing average correction.
    #    in more than one term is given (blank separated string!), the correction are applied in sequence.
    # %param% spec_noise_n integer / default 10
    #   used by the offset algorithm to determine the offset to correct
    # %param% f1_bcorr_pivots Integerlist / default (int(0.05*get_si1_2d()),int(0.95*get_si2_2d()))
    #   pivot points used by the linear or spline algorithm
    # %param_cond% (f1_bcorr_pivots > 0 and bcorr_pivots <= get_si1_2d())
    # %param% f1_bcorr_radius integer / default 1
    #   radius around pivot points used by the linear or spline algorithm
    # %param% f1_bcorr_order integer  / default 10
    # %param_cond% (f1_bcorr_order > 0 and f1_bcorr_order < get_si1_2d() and f1_bcorr_order == 2*n)
    #   order (number of points corrected in the time domain) of the quest algorithm
    if(key_is_true(p_in,"f1_baseline")):
        set_task("F1 baseline correction")
        toreal("f1")
        pivots = "not used" 
        bcorder = "not used"
        bc = p_in["f1_bcorr_algo"]
        for bcalgo in bc.split():
            if (bcalgo == "offset"):
                bcorr_offset( p_in["spec_noise_n"], "f1")
            elif (bcalgo == "linear"):
                pivots=p_in["f1_bcorr_pivots"]
                bcorr(1,p_in["f1_bcorr_radius"],"f1",pivots)
            elif (bcalgo == "spline"):
                pivots=p_in["f1_bcorr_pivots"]
                bcorr(2,p_in["f1_bcorr_radius"],"f1",pivots)
            elif (bcalgo == "quest"):
                bcorder = p_in["f1_bcorr_order"]        
                bcorr_quest( bcorder, "f1")
            elif (bcalgo == "polynomial"):  # bcorrpx determines which algo is used
                if (get_debug()):
                    print("Baseline correction not available in DEBUG mode")
                else:
                    bcorrp0()
                    bcorr(3, "f1")
            elif (bcalgo == "moving_average"):
                if (get_debug()):
                    print("Baseline correction not available in DEBUG mode")
                else:
                    bcorrp1()
                    bcorr(3, "f1")
                    bcorrp0()
            else:
                raise "Error with baseline correction algorithm"
        audittrail( audit, "text", ("F1 Baseline correction -" + bc),
                         "algorithm", bcalgo,
                         "pivots", pivots,
                         "order", bcorder)

    # %action% shear
    #  apply the shear correction
    # %param% shear boolean / default 0
    # %param% shear_slope float / default -1
    #   the slope of the shear operation, -1 corresponds to the INADEQUATE to COSY correction
    # %param% shear_pivot float / default 0.5
    #   the location of the invariant point
    # %param_cond% [ (shear_pivot>=0) and (shear_pivot<=1) ] 
    if(key_is_true(p_in,"shear")):
            set_task("Shearing")
            shear( p_in["shear_slope"], p_in["shear_pivot"])
            audittrail( audit, "text", "shear correction applied",
                        "slope", p_in["shear_slope"],
                        "pivot", p_in["shear_pivot"])

    # %action% tilt
    #  apply the tilt correction
    # %param% tilt boolean / default 0
    # %param% tilt_slope float / default 1
    #   the slope of the tilt operation, 1 corresponds to the J-Res correction
    # %param% tilt_pivot float / default 0.5
    #   the location of the invariant point
    # %param_cond% [ (tilt_pivot>=0) and (tilt_pivot<=1) ] 
    if (key_is_true(p_in,"tilt")):
            set_task("Tilting")
            tilt( p_in["tilt_slope"], p_in["tilt_pivot"])
            audittrail( audit, "text", "tilt correction applied",
                        "slope", p_in["tilt_slope"],
                        "pivot", p_in["tilt_pivot"])

    # %action% symetrize
    #  requires the data-set to be square, will reduce to the smallest of both axes
    # %param% symetrize boolean / default 0
    # %param% symetrize_algo enum mean smallest / default smallest
    if(key_is_true(p_in,"symetrize")):
            set_task("Spectrum symetrisation")
            if (get_si1_2d() > get_si2_2d()):
                chsize( 2*power2(get_si1_2d()-2), get_si2_2d())
                iftbis("f1")
                chsize(get_si2_2d(), get_si2_2d())
                ftbis("F1")
            elif (get_si1_2d() < get_si2_2d()):
                chsize( get_si1_2d(), 2*power2(get_si2_2d()-2))
                iftbis("f2")
                chsize( get_si1_2d(), get_si1_2d())
                ftbis("F2")

            if   (p_in["symetrize_algo"].lower() == "mean"):
                sym(1)
            elif (p_in["symetrize_algo"].lower() == "smallest"):
                sym(2)

            audittrail( audit, "text", "Symetrisation applied",
                       "algorithm", p_in["symetrize_algo"],
                       "final size  (F1=F2)", get_si1_2d())

    # %action% smoothing
    #   apply a smoothing filter to the data-set
    # %param% smoothing boolean / default 0
    # %param% smooth_f1_w integer / default 2
    #    size of the smoothing window in f1
    # %param% smooth_f2_w integer / default 2
    #    size of the smoothing window in f2
    # %param% smooth_iteration / default 1
    #   number of loop
    # %param_cond% [ smooth_f1_w > 0 ]
    # %param_cond% [ smooth_f2_w > 0 ]
    if (key_is_true(p_in,"smoothing")):
        set_task("Spectrum smoothing")
        if (key_is_true(p_in,"smooth_iteration")):
            it = p_in["smooth_iteration"]
        else:
            it = 1
        for i in range(it):
            smooth( p_in["smooth_f1_w"], p_in["smooth_f2_w"])
        audittrail( audit,  "text", "Smoothing filter",
                    "F1 number of points", p_in["smooth_f1_w"],
                    "F2 number of points", p_in["smooth_f2_w"],
                    "Number of iterations", p_in["smooth_iteration"])

    # %action% median
    #   apply a median filter to the data-set
    # %param_cond% [ [ median == ((f1_extract_real and f2_extract_real) | modulus) ] 
    # %param% median boolean / default 0
    # %param% median_f1_w integer / default 2
    #    size of the median window in f1
    # %param% median_f2_w integer / default 2
    #    size of the median window in f2
    # %param% median_i integer <= median_w  / default 2
    #    index of point to keep in the median filtering window
    # %param_cond% [ median_f2_w > 0 ]
    # %param_cond% [ median_f1_w > 0 ]
    # %param_cond% [ (median_i > 0) and (median_i <= (median_f1_w*median_f2_w)) ]
    if (key_is_true(p_in,"median")):
        set_task("Median filter applied")
        median( p_in["median_f1_w"], p_in["median_f2_w"], p_in["median_i"]) 
        audittrail( audit,  "text", "Median filter",
                    "F1 number of points", p_in["median_f1_w"],
                    "F2 number of points", p_in["median_f2_w"])

    # %action% projection
    #  creates and store in files the projections of the 2D
    #  2 files are created per axis :
    #  using a mean algorithm (_M suffix) and using a skyline algorithm (_S suffix)
    # %param% projection boolean / default 1
    # %return% f1_projection_M
    # %return% f2_projection_M
    # %return% f1_projection_S
    # %return% f2_projection_S
    if(key_is_true(p_in,"projection")):
            set_task("Projection computation")
            proj("f1", "M")
            dim(1)
            projname = os.path.join(location, 'projf1_M.gifa4')
            writec(projname)
            f_out["f1_projection_M"] = projname
            dim(2)
            proj("f1", "S")
            dim(1)
            projname = os.path.join(location, 'projf1_S.gifa4')
            writec(projname)
            f_out["f1_projection_S"] = projname
            dim(2)
            proj("f2", "M")
            dim(1)
            projname = os.path.join(location, 'projf2_M.gifa4')
            writec(projname)
            f_out["f2_projection_M"] = projname
            dim(2)
            proj("f2", "S")
            dim(1)
            projname = os.path.join(location, 'projf2_S.gifa4')
            writec(projname)
            f_out["f2_projection_S"] = projname
            dim(2)

            audittrail( audit,  "text", "Projections stored",
                        "F1 Mean algorithm file", f_out["f1_projection_M"],
                        "F1 Skyline algorithm file", f_out["f1_projection_S"],
                        "F2 Mean algorithm file", f_out["f2_projection_M"],
                        "F2 Skyline algorithm file", f_out["f2_projection_S"])

    # %action% spec_noise
    #   evaluate noise, estimated by finding an empty zone
    # %param% spec_noise_n integer / default 10
    #    number of different zones where noise is evaluated
    # %return% spec_std_noise
    #  estimate of the noise in the spectrum
    set_task("Noise evaluation")
    spec_noise( p_in["spec_noise_n"])
    f_out["noise"] = get_noise()

    # %action% select_state
    #   permit to choose the state complex / real of the output file
    #   complex data are changed to real by dropping the imaginary part
    #   real data are changed to complex by computing the Hilbert transform with tocomplex()
    #   this is usually not required as all processing commands prepare the state themselves 
    # %param% select_state boolean / default 0
    #   actuallly does the selection
    # %param% f1_state enum ignore complex real / default ignore
    #   force the f1 axis to real or complex. ignore will let the axis unchanged
    # %param% f2_state enum ignore complex real / default ignore
    #   force the f2 axis to real or complex. ignore will let the axis unchanged
    if(key_is_true(p_in,"select_state")):
        (t1,t2)=get_itype(2)
        p1s="ignore"
        p2s="ignore"
        if (p_in["f1_state"] == "complex"):
            if (t1==0):
                set_task("F1 complex reconstruction")
            tocomplex("f1")
            p1s="complex"
        elif (p_in["f1_state"] == "real"):
            toreal("f1")
            p1s="real"
        if (p_in["f2_state"] == "complex"):
            if (t2==0):
                set_task("F2 complex reconstruction")
            tocomplex("f2")
            p2s="complex"
        elif (p_in["f2_state"] == "real"):
            toreal("f2")
            p2s="real"

        if t1 == 1:
            t1s="complex"
        else:
            t1s="real"
        if t2 == 1:
            t2s="complex"
        else:
            t2s="real"
        audittrail( audit,  "text", "Complex state of the data",
                        "F1 axis", "was "+t1s+" now forced to "+p1s,
                        "F2 axis", "was "+t2s+" now forced to "+p2s)

    write_file_2d( audit, fileout)

#---------------------------------------------------------------------------
def post_maxent_2d( audit, filein, fileout, p_in, f_in, f_out,location="." ):
    """
This macro realizes the post MaxEnt operation on a 2D spectrum
 
it implements the massaging of the spectrum after FT :
- level_correction      : correct the slight level offseting by doing a simple baseline correction
- calibration           : apply ppm calibration
- positive_negative_cor : reconstruct a composite 2D from the two channel processing
- smoothing             : apply a smoothing filter
- median                : apply a media filter
- projection            : compute F1 and F2 projections
- select_state          : enforce the final state of the data

arguments :
    audit       the opened audit file
    filein      the file name of the input file, can be "memory"
    fileout     the file name of the input file, can be "memory"
    p_in        the dictionary containing all the processing parameters
    f_in        the dictionary containing all the details on the spectrum "filein"
    f_out       details on the spectrum "fileout" will be put in this dictionary

%F1-input-domain%  frequency
%F1-output-domain% frequency
%F2-input-domain%  frequency
%F2-output-domain% frequency
%dimensionality%  2

%author% Marc-Andre Delsuc
%version% 6.0

    """
    read_file_2d( audit, filein)

    # %action% level_correction
    #   correct the slight level offseting by doing a simple baseline correction
    # %param% level_correction boolean / default 1
    # %param% f2_level_correction_axis boolean / default 1
    #   apply level on F2 axis
    # %param% f2_level_correction_pivots list / default (int(0.01*get_si2_2d()),int(0.99*get_si2_2d()))
    #   where level will be computed
    # %param% f1_level_correction_axis boolean / default 0
    #   apply level on F1 axis
    # %param% f1_level_correction_pivots list / default (int(0.01*get_si1_2d()),int(0.99*get_si1_2d()))
    #   where level will be computed
    print(p_in["f2_level_correction_pivots"])
    if key_is_true(p_in,"level_correction"):
        if key_is_true(p_in,"f2_level_correction_axis"):
            bcorr(1,1,"F2",p_in["f2_level_correction_pivots"])
            audittrail( audit, "text", "level correction applied applied",
               "along F2 axis:", p_in["f2_level_correction_pivots"])
        if key_is_true(p_in,"f1_level_correction_axis"):
            bcorr(1,1,"F1",p_in["f1_level_correction_pivots"])
            audittrail( audit, "text", "level correction applied applied",
               "along F1 axis:", p_in["f1_level_correction_pivots"])

    # %action% positive_negative_cor
    #   reconstruct a composite 2D from the two channel processing
    # %param% positive_negative_cor boolean / default 0
    if key_is_true(p_in,"positive_negative_cor"):
        posneg_corr()
        audittrail( audit, "text", "positive_negative correction applied applied")

    # %action% calibration
    #   calibrate 0ppm on the spectrum
    # %param% f2_calibration float / default 0
    #   the coordinate in ppm of the right-most point in the spectrum
    # %param% f1_autocalibrate boolean / default 1
    #   calibrate the F1 axis based on the F2 calibration.
    #   overwrite the F1 calibration.
    # %param% f1_calibration float / default 0
    #   the coordinate in ppm of the lower-most point in the spectrum
    offset( get_offset_1_2d(), p_in["f2_calibration"]*get_freq_2_2d())
    audittrail( audit, "text", "F2 offset applied",
               "calibration", get_offset_2_2d()/get_freq_2_2d())

    if (key_is_true(p_in,"f1_autocalibrate")):
        spin=autocalib(p_in["autocalibrate_mode"])
        audittrail( audit, "text", "F1 autocalibration applied",
                       "mode",p_in["autocalibrate_mode"],
                       "spin",spin,
                       "calibration", get_offset_1_2d()/get_freq_1_2d())
    else:
        offset(p_in["f1_calibration"]*get_freq_1_2d(), get_offset_2_2d())
        audittrail( audit,  "text", "F1 offset applied",
                        "calibration", get_offset_1_2d()/get_freq_1_2d())

    # %action% smoothing
    #   apply a smoothing filter to the data-set
    # %param% smoothing boolean / default 0
    # %param% smooth_f1_w integer / default 2
    #    size of the smoothing window in f1
    # %param% smooth_f2_w integer / default 2
    #    size of the smoothing window in f2
    # %param% smooth_iteration / default 1
    #   number of loop
    # %param_cond% [ smooth_f1_w > 0 ]
    # %param_cond% [ smooth_f2_w > 0 ]
    if (key_is_true(p_in,"smoothing")):
        set_task("Spectrum smoothing")
        if (key_is_true(p_in,"smooth_iteration")):
            it = p_in["smooth_iteration"]
        else:
            it = 1
        for i in range(it):
            smooth( p_in["smooth_f1_w"], p_in["smooth_f2_w"])
        audittrail( audit,  "text", "Smoothing filter",
                    "F1 number of points", p_in["smooth_f1_w"],
                    "F2 number of points", p_in["smooth_f2_w"],
                    "Number of iterations", p_in["smooth_iteration"])

    # %action% median
    #   apply a median filter to the data-set
    # %param_cond% [ [ median == ((f1_extract_real and f2_extract_real) | modulus) ] 
    # %param% median boolean / default 0
    # %param% median_f1_w integer / default 2
    #    size of the median window in f1
    # %param% median_f2_w integer / default 2
    #    size of the median window in f2
    # %param% median_i integer <= median_w  / default 2
    #    index of point to keep in the median filtering window
    # %param_cond% [ median_f2_w > 0 ]
    # %param_cond% [ median_f1_w > 0 ]
    # %param_cond% [ (median_i > 0) and (median_i <= (median_f1_w*median_f2_w)) ]
    if (key_is_true(p_in,"median")):
        set_task("Median filter applied")
        median( p_in["median_f1_w"], p_in["median_f2_w"], p_in["median_i"]) 
        audittrail( audit,  "text", "Median filter",
                    "F1 number of points", p_in["median_f1_w"],
                    "F2 number of points", p_in["median_f2_w"])

    # %action% projection
    #  creates and store in files the projections of the 2D
    #  2 files are created per axis :
    #  using a mean algorithm (_M suffix) and using a skyline algorithm (_S suffix)
    # %param% projection boolean / default 1
    # %return% f1_projection_M
    # %return% f2_projection_M
    # %return% f1_projection_S
    # %return% f2_projection_S
    if(key_is_true(p_in,"projection")):
            set_task("Projection computation")
            proj("f1", "M")
            dim(1)
            projname = os.path.join(location, 'projf1_M.gifa4')
            writec(projname)
            f_out["f1_projection_M"] = projname
            dim(2)
            proj("f1", "S")
            dim(1)
            projname = os.path.join(location, 'projf1_S.gifa4')
            writec(projname)
            f_out["f1_projection_S"] = projname
            dim(2)
            proj("f2", "M")
            dim(1)
            projname = os.path.join(location, 'projf2_M.gifa4')
            writec(projname)
            f_out["f2_projection_M"] = projname
            dim(2)
            proj("f2", "S")
            dim(1)
            projname = os.path.join(location, 'projf2_S.gifa4')
            writec(projname)
            f_out["f2_projection_S"] = projname
            dim(2)

            audittrail( audit,  "text", "Projections stored",
                        "F1 Mean algorithm file", f_out["f1_projection_M"],
                        "F1 Skyline algorithm file", f_out["f1_projection_S"],
                        "F2 Mean algorithm file", f_out["f2_projection_M"],
                        "F2 Skyline algorithm file", f_out["f2_projection_S"])

    # %action% spec_noise
    #   evaluate noise, estimated by finding an empty zone
    # %param% spec_noise_n integer / default 10
    #    number of different zones where noise is evaluated
    # %return% spec_std_noise
    #  estimate of the noise in the spectrum
    set_task("Noise evaluation")
    spec_noise( p_in["spec_noise_n"])
    f_out["noise"] = get_noise()

    # %action% select_state
    #   permit to choose the state complex / real of the output file
    #   complex data are changed to real by dropping the imaginary part
    #   real data are changed to complex by computing the Hilbert transform with tocomplex()
    #   this is usually not required as all processing commands prepare the state themselves 
    # %param% select_state boolean / default 0
    #   actuallly does the selection
    # %param% f1_state enum ignore complex real / default ignore
    #   force the f1 axis to real or complex. ignore will let the axis unchanged
    # %param% f2_state enum ignore complex real / default ignore
    #   force the f2 axis to real or complex. ignore will let the axis unchanged
    if(key_is_true(p_in,"select_state")):
        (t1,t2)=get_itype(2)
        p1s="ignore"
        p2s="ignore"
        if (p_in["f1_state"] == "complex"):
            if (t1==0):
                set_task("F1 complex reconstruction")
            tocomplex("f1")
            p1s="complex"
        elif (p_in["f1_state"] == "real"):
            toreal("f1")
            p1s="real"
        if (p_in["f2_state"] == "complex"):
            if (t2==0):
                set_task("F2 complex reconstruction")
            tocomplex("f2")
            p2s="complex"
        elif (p_in["f2_state"] == "real"):
            toreal("f2")
            p2s="real"

        if t1 == 1:
            t1s="complex"
        else:
            t1s="real"
        if t2 == 1:
            t2s="complex"
        else:
            t2s="real"
        audittrail( audit,  "text", "Complex state of the data",
                        "F1 axis", "was "+t1s+" now forced to "+p1s,
                        "F2 axis", "was "+t2s+" now forced to "+p2s)
    write_file_2d( audit, fileout)

#---------------------------------------------------------------------------
def pre_2d( audit, filein, fileout, p_in, f_in, f_out ):
    """
This macro realizes the preparation operation on a 2D FID

it implements global FID evaluation before 2D FT:
- data_check       : performs minimum data checking 
- fid_noise        : measure of noise and offset levels in FID of noise and offset levels in FID
- conv_n_p         : conversion of N+P acquisition (echo / antiecho)

arguments :
    audit       the opened audit file
    filein      the file name of the input file, can be "memory"
    fileout     the file name of the input file, can be "memory"
    p_in        the dictionary containing all the processing parameters
    f_in        the dictionary containing all the details on the spectrum "filein"
    f_out       details on the spectrum "fileout" will be put in this dictionary

%F1-input-domain%  time
%F2-input-domain%  time
%F1-output-domain%  time
%F2-output-domain%  time
%dimensionality%   2

%author% Marc-Andre Delsuc
%version% 6.0

    """
    read_file_2d( audit, filein)
    # adhoc ! MAD

    # %action% data_check
    #   performs minimum data integrity
    # %param% data_check boolean / default 1
    #   checks points which can not be checked by the kernel
    #   for the moment, only check data-sets with an odd number of points
    if (key_is_not_false(p_in,"data_check")):
        # check size parity and complex type
        si1=get_si1_2d()
        si2=get_si2_2d()
        if key_is_true(p_in,"f1_fourier_transform"):
            t1=1   # assumes complex, then check for cases which require real
        else:
            t1=0
        if key_is_true(p_in,"f2_fourier_transform"):
            t2=1   # assumes complex, then check for cases which require real
        else:
            t2=0
        # case of odd number of points
        if (si2 % 2 == 1):                   # if SI2 odd
            t2=0
            if key_is_true(p_in,"f2_fourier_transform"):
                if p_in["f2_ft_type"] != 'rft':     # only real case in F2
                    si2 = si2-1        
                    t2=1
        if (si1 % 2 == 1):                   # if SI1 odd
            t1=0
            if key_is_true(p_in,"f1_fourier_transform"):
                if key_is_true(p_in,"conv_n_p") or key_is_true(p_in,"f1_conv_n_p"):
                    si1 = si1-1        
                    t1=1
            if key_is_true(p_in,"f1_fourier_transform"):
                if (p_in["f1_ft_type"] not in ('rft', 'ft_tppi', 'ft_phase_modu')):
                    si1 = si1-1        
                    t1=1
        itype(2*t1+t2)
        chsize(si1,si2)
        audittrail( audit, "text", "Data-set:",
            "size (F1,F2)",str(si1)+" x "+str(si2),
            "Complex type",(t1,t2))

    # %action% fid_noise 
    #   evaluate noise and offset, estimated from the last points
    # %param% fid_noise_F2zone float / default 0.2
    #   zone of the 2D FID (expressed in % from the end) where to compute noise
    # %param_cond% [ (fid_noise_F2zone > 0.0) and (fid_noise_F2zone<=1)]
    # %param% fid_noise_f1zone float / default 0.1
    #   zone of the 2D FID (expressed in % from the end) where to compute noise
    # %param_cond% [ (fid_noise_f1zone > 0.0) and (fid_noise_f1zone<=1)]
    # %return% fid_offset
    #  constant offset of the FID, estimated on the last points of the dataset
    # %return% NOISEINTIMEDOMAIN
    #  estimate of the fid noise, from the standard deviation of the last points of the dataset
    #set_task("FID noise evaluation")
    try:
        lastf1 = p_in["fid_noise_f1zone"]
    except:
        lastf1 = 0.1
    try:
        lastf2 = p_in["fid_noise_f2zone"]
    except:
        lastf2 = 0.2
    evaln ((1.0-lastf1)*get_si1_2d(), (1-lastf2)*get_si2_2d(), get_si1_2d(), get_si2_2d())

    # loads $noise and $shift

    f_out["fid_offset"] = get_shift()
    f_in["noise"] = get_noise()
    audittrail( audit, "text", "FID Noise evaluated","Noise:",get_noise())
    
    # %action% conv_n_p
    #   convert experiment acquired in n+p (echo/antiecho mode) to hypercomplex
    # %param% conv_n_p boolean / default 0
    # %param_exclusive%  [ conv_n_p ft_phase_modu ]
    if key_is_true(p_in,"conv_n_p") or key_is_true(p_in,"f1_conv_n_p"):
            set_task("Conversion from N+P acquisition mode")
            itype(1)
            conv_n_p()
            audittrail( audit, "text", "Conversion from N+P acquisition mode applied")

    write_file_2d( audit, fileout)

#---------------------------------------------------------------------------
def read_file_2d( audit, filein):
    """
    used to read the input binary file, while maintaining the audittrail
    """
    dim(2)
    if ( filein != "memory" ):
        try:
            com_read( filein )
        except:
            raise "Error while reading file ",filein
        if ( get_dim() != 2 ):
            raise "This is not a 2D dataset"

    if ( get_itype_2d() == 2 or get_itype_2d() == 3 ):
        sz = repr(get_si1_2d()/2) + " Complex x"
    else:
        sz = repr(get_si1_2d()) + " Real x"

    if ( get_itype_2d() == 1 or get_itype_2d() == 3 ):
        sz = sz + repr(get_si2_2d()/2) + " Complex"
    else:
        sz = sz + repr(get_si2_2d()) + " Real"

    if ( filein != "memory" ):
            audittrail( audit, "text", "loading in memory", "filename", filein, "size in F1 x F2", sz)
    else:
            audittrail( audit, "text", "using data in memory", "size in F1 x F2", sz)
            
    

#---------------------------------------------------------------------------
def write_file_2d( audit, fileout):
    """
    used to write the output binary file, while maintaining the audittrail
    """
    if ( fileout != "memory" ):
        dim(2)
        writec(fileout)
        if ( get_itype_2d() == 2 or get_itype_2d() == 3 ):
            sz = repr(get_si1_2d()/2) + " Complex x"
        else:
            sz = repr(get_si1_2d()) + " Real x"
            
        if ( get_itype_2d() == 1 or get_itype_2d() == 3 ):
            sz = sz + repr(get_si2_2d()/2) + " Complex"
        else:
            sz = sz + repr(get_si2_2d()) + " Real"
        audittrail( audit, "text", "output file created :",
                   "file name", fileout, "size in F1 x F2", sz)


