#!/usr/bin/env python 
# encoding: utf-8

"""
    This module contains all the routines needed to process 1D NMR spectra
    
    The following functions realize a set of operations, needed for NMR 1D processing,
    You will find operation for FT analysis, MaxEnt Analysis,
    Inverse Fourier, etc..
    
    Processing are divided in pre - FT - Post phases.
    For instance a typical 1D procesing, would be

    pre()
    ft()
    post()

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
    
"""
from __future__ import print_function
__author__ = "Marc A. Delsuc <delsuc@igbmc.fr>"
__date__ = "Oct 2009"

import math
from Generic import *


#---------------------------------------------------------------------------
def FT1D(audit,p_in_arg,f_in,f_out,inputfilename,outputfilename):
    """
    FT processing of a 1D FID
    
    based on pre_ft_1d() ft_1d() post_ft_1d()
    """
    p_in = build_dict(("pre_ft_1d","ft_1d","post_ft_1d"),p_in_arg)

    # Pre ft phase
    audittrail( audit, "phase", "FID preparation phase")
    pre_ft_1d( audit, inputfilename, "memory", p_in, f_in, f_out)
    
    # ft phase
    audittrail( audit, "phase", "Spectral analysis phase")
    ft_1d( audit, "memory", "memory", p_in, f_in, f_out)

    # Post ft phase
    audittrail( audit, "phase", "Post processing phase")
    post_ft_1d( audit, "memory", outputfilename, p_in, f_in, f_out)



#---------------------------------------------------------------------------
def MaxEnt1D(audit,p_in_arg,f_in,f_out,inputfilename,outputfilename):
    """
    MaxEnt processing of a 1D FID

    based on pre_ft_1d() maxent_1d() post_maxent_1d()
    """
    p_in = build_dict(("pre_ft_1d","maxent_1d","post_maxent_1d"),p_in_arg)

    # Pre ft phase
    audittrail( audit, "phase", "FID preparation phase")
    pre_ft_1d( audit, inputfilename, "memory", p_in, f_in, f_out)
        
    # maxent phase
    audittrail( audit, "phase", "Spectral analysis phase by Maximum Entropy")
    maxent_1d( audit, "memory", "memory", p_in, f_in, f_out)

    # Post ft phase
    audittrail( audit, "phase", "Post processing phase")
    post_maxent_1d( audit, "memory", outputfilename, p_in, f_in, f_out)




#---------------------------------------------------------------------------
def pre_ft_1d(  audit, filein, fileout, p_in, f_in, f_out):
    """
This macro realizes the pre FT operation on a 1D FID

- fid_noise        : evaluate of noise and offset levels in FID
- dc_offset        : corrects for constant offset in FID
- causalize        : changes DSP processed FID (Bruker) to causal FID's by Hilbert transform
- flatten_solvent  : removes solvent signal by FID analysis
- left_shift       : drops first points of the FID 
- right_shift      : adds empty points on the beginning of the FID
- back_extend      : reconstructs missing points in the beginning of the FID by LP analysis

%F1-input-domain%  time
%F1-output-domain% time
%dimensionality%   1

%author% Marc-Andre Delsuc
%version% 6.0
    """
    dim(1)
    if ( filein != "memory" ):
        read( filein )
        if ( get_dim() != 1 ):
            raise "This is not a 1D dataset"
        audittrail( audit, "text", "loading in memory", "filename", filein)
    else:
        audittrail( audit, "text", "using data in memory")

    # %action% fid_noise
    #  evaluate noise and offset, estimated from the last points
    # %param% fid_noise_zone float / default 0.2
    #   zone of the FID (expressed in % from the end) where to compute noise
    # %param_cond% [ (fid_noise_zone > 0.0) and (fid_noise_zone<=1)]
    # %return% fid_offset
    #   constant offset of the FID, estimated on the last points of the dataset (defined by fid_noise_zone)
    # %return% NoiseInTimeDomain
    #   estimate of the fid noise, from the standard deviation of the last points of the dataset (defined by fid_noise_zone)
    set_task("FID noise evaluation")
    last = p_in["fid_noise_zone"]        # usually 0.2 which means last 20 % of data set
    evaln((1-last)*get_si1_1d(), get_si1_1d())       # loads get_noise and get_shift

    f_out["fid_offset"] = get_shift()
    f_in["Noise"] = get_noise()
    f_out["Noise"] = get_noise() # a verifer

    # %action% dc_offset
    #   removes a constant level in the fid, estimated from the last points
    # default has been set to off (0), as it seems the 
    # %param% dc_offset boolean  / default 0
    if (key_is_true(p_in,"dc_offset")):
        set_task("DC-offset correction")
        s = get_shift()
        addbase(s)
        audittrail( audit, "text", "vertical shift of FID","offset",s)

    # %param_exclusive% causal [ causalize causal_corr ]
    if (key_is_true(p_in,"causalize") and key_is_true(p_in,"causal_corr")):
        raise """
        causalize and causal_corr are incompatible
        causality must be corrected only once"""

    # %action% causalize
    #   removes the non causal header found in certain (Bruker) implementation of the digital filter.
    # %param% causalize boolean / default 0
    # %use% zerotimeposition
    # %return% causalize
    #   contains the correction applied by causalize
    if (key_is_true(p_in,"causalize")):
        set_task("Causal correction (Bruker digital filter correction")
        causalize(f_in["axisf1_zerotimeposition"])
        f_out["causalize"] = (-360*f_in["axisf1_zerotimeposition"])
        audittrail( audit, "text", "cauzalisation correction",
                   "1st order phase correction", f_out["causalize"])

    # %action% flatten_solvent
    #   removes water from fit of the FID in the time domain
    # %param% flatten_solvent boolean  / default 0
    # %param% flat_solv_mode enum polynomial moving_average polynomial+moving_average moving_average+polynomial / default polynomial
    #    algorithm used
    if (key_is_true(p_in,"flatten_solvent")):
        set_task("Solvent flattening")
        flat_solvent(p_in["flat_solv_mode"])
        audittrail( audit, "text", "apply solvent flattening", "algorithm", p_in["flat_solv_mode"])

    # %action% left_shift
    #   shifts the FID to the left by dropping data points
    # %param% left_shift boolean  / default 0
    # %param% left_shift_size integer  / default 2
    #   number of points to shift
    # %param_cond%  [ (left_shift_size > 0) and (left_shift_size < get_si1_1d()) ]
    if (key_is_true(p_in,"left_shift")):
        set_task("Left shifting of FID")
        left_shift(p_in["left_shift_size"])
        audittrail( audit, "text", "FID left shift", "number of points", p_in["left_shift_size"])

    # %action-group% exclusive [ right_shift b_extend ]
    if (key_is_true(p_in,"right_shift") and key_is_true(p_in,"b_extend")):
        raise """
        right_shift and b_extend are incompatible
        """

    # %action% right_shift
    #   shifts the FID to the right by adding zero data points
    # %param% right_shift boolean  / default 0
    # %param% right_shift_size integer  / default 2
    #   number of points to shift
    # %param_cond%  [ (right_shift_size > 0) ]
    if (key_is_true(p_in,"right_shift")):
        set_task("Right shifting of FID")
        right_shift(p_in["right_shift_size"])
        audittrail( audit, "text", "FID right shift", "number of points", p_in["right_shift_size"])

    # %action% b_extend 
    #   extends the FID backward for reconstruction of initial missing points
    # %param% b_extend boolean  / default 0
    # %param% b_extend_size integer / default 2
    #   number of missing points to reconstruct
    # %param_cond% [ b_extend_size > 0 ]
    # %param% b_extend_algo enum burg svd  / default burg
    #   algorithm used
    # %return% b_extend_order 
    #    the order used during back_extend
    # %return% b_extend_size 
    #    the number of points added by back_extend
    if (key_is_true(p_in,"b_extend")):
        set_task("Backward extension")
        n = p_in["b_extend_size"]
        order(10)  # (min(30,($n*2)))
        if (p_in["b_extend_algo"] == "burg" ):
            burg_back(get_si1_1d() + n)

        elif (p_in["b_extend_algo"] == "svd" ):
            dt2svd( get_si1_1d() )
            svd2ar(2)
            ar2dt(get_si1_1d() + n,  2)


        f_out["b_extend_order"] = get_order()
        f_out["b_extend_size"] = n
        audittrail( audit, "text", "FID back extension",
                   "number of points", p_in["b_extend_size"],
                   "algorithm", p_in["b_extend_algo"])


    write_file_1d( audit, fileout)


#---------------------------------------------------------------------------
# ft_1d
#---------------------------------------------------------------------------
def ft_1d(  audit, filein, fileout, p_in, f_in, f_out ):
    """
This macro realizes the FT operation on a 1D FID

- truncate      : truncates the FID by removing the last points
- lp_extend     : extend FID with a Linear Prediction algorithm
- apodize       : standard apodisation by a window
- fourier_transform : performs the Fourier transform
- causal_corr   : performs causal correction of the spectrum if not done on the FID
- reverse       : reverses the spectral axis after FT


%F1-input-domain%  time
%F1-output-domain% frequency
%dimensionality%   1

%author% Marc-Andre Delsuc
%version% 6.0

    """
    dim(1)
    if ( filein != 'memory' ):
        read( filein )
        if ( get_dim() != 1 ):
            raise "This is not a 1D dataset"
        audittrail( audit, "text", "loading in memory", "filename", filein)
    else:
        audittrail( audit, "text", "using data in memory")

    # %action% truncate
    #   truncates the FID by removing the last points
    # %param% truncate boolean / default 0
    # %param% trunc_size integer / default get_si1_1d()
    #   FID size after truncation
    # %param_cond% [ (trunc_size > 0 ) and (trunc_size <= get_si1_1d() ) ]
    if ( key_is_true(p_in,"truncate")):
        set_task("FID truncation")
        s = p_in["trunc_size"]
        if ( s <= 0 ):
            raise "Cannot truncate to 0"
        if ( s > get_si1_1d() ):
            raise "Error with truncation size (larger than actual size)"
        initial = get_si1_1d()
        s = int(min(get_si1_1d(), s))
        chsize( s )
        audittrail( audit, "text", "FID truncation before FT", "initial_size", initial, "final size", get_si1_1d())

    # %action% lp_extend
    #    extend FID with a Linear Prediction algorithm
    # %param% lp_extend boolean / default 0
    # %param% lp_ext_size integer / default 2*get_si1_1d()
    #   final size of FID
    # %param_cond% [ lp_ext_size > get_si1_1d() ]
    # %param% lp_ext_algo enum burg mirror lpsvd lpsvd_stable / default burg
    #   algorithm used
    #   burg and mirror are much faster, svd_stable is much slower
    #   mirror is to be used when the phase of the spectrum is known before hand (see lp_ext_off)
    # %param% lp_ext_off integer / default 0
    #   offset determines the position of the t=0 point, used by mirror algo
    #     0           no shift : in-phase data set.
    #    -1           acquisition started exactly half a dwell after t=0 - (will need phase 0 180)
    #    lp_ext_off>0 acquisition started exactly n dwell before t=0
    # %param_cond% [ lp_ext_off > -1 ]
    # %param% lp_ext_order integer / default (min(get_si1_1d()/2,20))
    #    the size of the prediction polynomial used for LP, a rough estimate of the complexity
    #  %param_cond% [ (lp_ext_order < get_si1_1d()/2 ) and (lp_ext_order > 2) ]
    #  %param% lp_ext_apod boolean / default 1
    #   apply a sine bell apodisation after LP extenstion.
    if ( key_is_true(p_in,"lp_extend")):
        set_task("Linear prediction extension")
        order( int(p_in["lp_ext_order"]) )

        if ( p_in["lp_ext_algo"] == "burg" ):
            t=get_itype_1d() ; itype(1)     # assume to be complex, burg is not defined on real!
            burg( p_in["lp_ext_size"] )
            print("BURG",p_in["lp_ext_size"],get_si1_1d(), get_order(), t)
            itype(t)

        elif ( p_in["lp_ext_algo"] == "mirror" ):
            t=get_itype_1d() ; itype(1)     # assume to be complex, burg is not defined on real!
            burg( p_in["lp_ext_size"] )
            burg_mirror( p_in["lp_ext_off"], p_in["lp_ext_size"])
            itype(t)

        elif ( p_in["lp_ext_algo"] == "lpsvd" ):
            t=get_itype_1d() ; itype(1)     # assume to be complex, burg is not defined on real!
            burg( p_in["lp_ext_size"] )
            dt2svd( get_si1_1d() )
            svd2ar(1)
            ar2dt( p_in["lp_ext_size"], 1)
            itype(t)

        elif ( p_in["lp_ext_algo"] == "lpsvd_stable" ):
            t=get_itype_1d() ; itype(1)     # assume to be complex, burg is not defined on real!
            burg( p_in["lp_ext_size"] )
            dt2svd( get_si1_1d() )
            svd2ar(1)
            ar2rt(1)
            rtreflect(1)
            rt2ar(1)
            ar2dt( p_in["lp_ext_size"], 1)
            itype(t)

        audittrail(audit, "text", "FID extension before FT",
                  "algorithm", p_in["lp_ext_algo"],
                  "LP order", p_in["lp_ext_order"],
                  "Final size", p_in["lp_ext_size"])

        if ( key_is_not_false(p_in,"lp_ext_apod")):
            sin(0)
            audittrail( audit, "text", "Apodisation after extension", "function", "sin(0)")

    # %action% apodize
    #   standard apodisation by a window
    # %param% apodize boolean / default 1
    # %param% apodisation string / default "expbroad(0.5)"
    #   the string describing the apodisation to apply
    #    "apodisation" should be a suite of window apodisation commands
    #     ie : 'sin(0.5)' or 'exbroad(0.2), sin(0.5)'
    #     the following predefined functions are implemnted
    #     sin (sine bell) sqsin (squared sine bell) expbroad (exponential)
    #     gaussbroad (gaussian) gaussenh (gaussian enhancement)
    if ( key_is_true(p_in,"apodize") and (p_in["apodisation"] != "none")):
        set_task("Apodisation")
        apodise(p_in.raw("apodisation"))
        audittrail( audit, "text", "FID apodisation before FT",
                   "function", p_in["apodisation"])


    # %action% fourier_transform
    #  performs the Fourier transform
    #  size after FT is determined automatically from FID size and apodisation.
    #  the optimum size is determined heuristically from 
    #  - the size of the FID (after truncation) : default is to zerofill once
    #  - the apodisation : resolution is not more than twice the filtered linewidth
    # %param% fourier_transform boolean / default 1
    # %param% ask_for_ftsize boolean / default 0
    #  if set, the final size is determined by ft_size, automatically determined otherwise
    #  it will be extended to the next 2^n
    # %param% ft_size integer  / default 2*power2(get_si1_1d())
    # %param% ft_type enum none ft_seq ft_sim ft rft / default ft_sim
    #    the Fourier transform algorithm, depends on spectrometer and acquisition scheme
    # %param_cond%   [ ispower2(ft_size) ]
    # %return% size_for_ft
    #   contains the data size right after FT (will usually be different from final size)
    if ( key_is_true(p_in,"fourier_transform")):
        set_task("Fourier transformation")
        if ( key_is_true(p_in,"ask_for_ftsize")):
            final = int (p_in["ft_size"])
        else:
            if ( key_is_true(p_in,"lp_extend")):
                final = 2*power2(get_si1_1d()-2)
            else:
                final = 4*power2(get_si1_1d()-2)
        initial = get_si1_1d()
        chsize(final)
        if ( p_in["ft_type"] != "none" ):
            if ( p_in["ft_type"] == "ft" ):
                ft()
            elif ( p_in["ft_type"] == "rft" ):
                rft()
            elif ( p_in["ft_type"] == "ft_sim" ):
                ft_sim()
            elif ( p_in["ft_type"] == "ft_seq" ):
                ft_seq()
            else:
                raise "Unknown Fourier transform"
                # utiliser la macro de nom p_in["ft_type"]
        f_out["size_for_ft"] = final
        audittrail( audit, "text", "Fourier transformation",
                   "FT type", p_in["ft_type"],
                   "initial_size", initial,
                   "final size", final)

    # %action% causal_corr
    #   performs causal correction of the spectrum if not done on the FID
    # %param% causal_corr boolean / default 1
    # %param_exclusive% causal [ causalize causal_corr ]
    # %use% axisF1_zerotimeposition 
    # %return% causalize
    #   contains the correction applied by causalize
    if ( key_is_true(p_in,"causal_corr")):
        set_task("Causal correction")
        try:
            delay = float(f_in["axisf1_zerotimeposition"])
        except:
            raise
            delay = 0.0
        com_phase(0.0, -360.0*delay )
        f_out["causalize"] = -360.0*delay
        audittrail( audit, "text", "Causal correction",
                   "1st order phase correction", f_out["causalize"])


    # %action% reverse
    #   reverses the spectral axis after FT
    # %param% reverse boolean / default 0
    if ( key_is_true(p_in,"reverse")):
        set_task("Spectrum reverse")
        reverse()
        if (get_itype_1d()==1):
            invf()       # if complex, inverse imaginary part
        audittrail( audit, "text", "Reverse spectral axis")

    # ecriture fichier de sortie
    write_file_1d( audit, fileout)



#---------------------------------------------------------------------------
# post_ft_1d
#---------------------------------------------------------------------------
def post_ft_1d( audit, filein, fileout, p_in, f_in, f_out):
    """
This macro realizes the Post processing of a 1D spectrum

- modulus       : takes the complex modulus of the spectrum
- phase         : applies a phase correction to the spectrum
- autophase     : automatically computes the phase correction of the spectrum
- invHilbert    : apply an inverse Hilbert transform
- calibration   : calibrate the ppm scale on the spectrum
- spectral_zone : extract one spectral zone of the spectrum
- baseline_correction : applies a baseline correction to the spectrum
- smoothing     : apply a smoothing filter to the data-set
- median        : apply a median filter to the data-set
- derivative    : compute the nth derivative of the data-set
- spec_noise    : evaluate noise, estimated by finding an empty zone



%F1-input-domain%  frequency
%F1-output-domain% frequency
%dimensionality%   1

%author% Marc-Andre Delsuc
%version% 6.0
    """
    dim(1)
    if ( filein != "memory" ):
        read( filein )
        if ( get_dim() != 1 ):
            raise "This is not a 1D dataset"
        audittrail( audit, "text", "loading in memory", "filename", filein)
    else:
        audittrail( audit, "text", "using data in memory")

    # need complex data
    tocomplex()


    # %param_exclusive% phase  [ modulus phase autophase ]
    # %action% modulus
    #   takes the complex modulus of the spectrum
    # %param% modulus boolean / default 0
    if (key_is_true(p_in,"modulus")):
        set_task("modulus computation")
        modulus()
        audittrail( audit, "text", "modulus applied")

    # %action% phase
    #   applies a phase correction to the spectrum
    # %param% phase boolean / default 0
    # %param% phase_0 float / default 0.0
    #  global order phase correction
    # %param% phase_1 float / default 0.0
    #  1st order phase correction
    if (key_is_true(p_in,"phase")):
        set_task("Phase correction")
        phase( p_in["phase_0"], p_in["phase_1"])
        audittrail( audit,  "text", "phase correction applied",
                     "0th order", p_in["phase_0"],
                     "1th order", p_in["phase_1"])


    # %action% autophase
    #   automatically computes the phase correction of the spectrum
    # %param% autophase boolean / default 1
    # %param% phase_algo enum apsl apmin / default apmin
    #   algorithm used
    # %return% autophase_0
    #    contains the value of the global order phase correction applied by autophase 
    # %return% autophase_1
    #    contains the value of the 1st order phase correction applied by autophase 
    if (key_is_true(p_in,"autophase")):
        set_task("Automatic phase correction")
        # evaluer la bonne macro @($p_in[phase_algo])
        if ( p_in["phase_algo"] == "apmin" ):
            apmin()                          
        if ( p_in["phase_algo"] == "apsl" ):
            apsl()
        f_out["autophase_0"] = get_ph0()
        f_out["autophase_1"] = get_ph1()
        audittrail( audit, "text", "automatic phase correction applied",
                   "algorithm", p_in["phase_algo"],
                   "0th order", get_ph0(),
                   "1th order", get_ph1())

    # %param_exclusive% real [ modulus invHilbert ]
    #
    # %action% invHilbert
    #   apply an inverse Hilbert transform.
    #   This operation takes a complex spectrum, and generates a real spectrum in-place
    #   i.e. with twice the number of real points.
    # %param% invHilbert boolean / default 0
    if (key_is_true(p_in,"invhilbert")):
        set_task("Inverse Hilbert Transform")
        invhilbert()
        audittrail( audit, "text", "inverse Hilbert transform applied")

    # %action% calibration
    #   calibrate the ppm scale on the spectrum
    # %param% calibration float / default 0.0
    #   the coordinate in ppm of the right-most point in the spectrum
    set_task("Calibration")
    offset (p_in["calibration"]*get_freq_1d())
    audittrail( audit,  "text", "offset applied", "calibration", p_in["calibration"])

    # %action% spectral_zone
    #   extract one spectral zone of the spectrum
    # %param% spectral_zone boolean / default 0
    # %param% spec_zone_left float / default 10.0
    #   the left border of the extract zone, in unit
    # %param% spec_zone_left_unit enum ppm hz index / default ppm
    #   the unit in which spec_zone_left is given
    # %param% spec_zone_right float / default 0.0
    #   the right border of the extract zone, in unit
    # %param% spec_zone_right_unit enum ppm hz index / default ppm
    #   the unit in which spec_zone_right is given
    # %param_cond% (spec_zone_left{index} < spec_zone_right{index})
    # %return% spec_zone_left
    #     the left coordinate of the extracted spectral zone in index
    # %return% spec_zone_right
    #     the right coordinate of the extracted spectral zone in index
    if (key_is_true(p_in,"spectral_zone")):
        set_task("Spectral zone extraction")
        if (p_in["spec_zone_left_unit"] == 'ppm'):
            l = ptoi(p_in["spec_zone_left"],1,1)
        elif (p_in["spec_zone_left_unit"] == 'hz'):
            l = htoi(p_in["spec_zone_left"],1,1)
        elif (p_in["spec_zone_left_unit"] == 'index'):
            l = p_in["spec_zone_left"]
        else:
            raise "Error with spec_zone_left_unit"

        if (p_in["spec_zone_right_unit"] == 'ppm'):
            r = ptoi(p_in["spec_zone_right"],1,1)
        elif (p_in["spec_zone_right_unit"] == 'hz'):
            r = htoi(p_in["spec_zone_right"],1,1)
        elif (p_in["spec_zone_right_unit"] == 'index'):
            r = p_in["spec_zone_right"]
        else:
            raise "Error with spec_zone_right_unit"

        l = max(1,round(l))
        r = min(get_si1_1d(),round(r))
        if (l > r):
            raise "Wrong spectral zone coordinates"

        if ((r-l) < 8):
            raise "spectral zone too small"

        if (l != 1 | r != get_si1_1d()):
            extract(int(l), int(r))

        f_out["spec_zone_left"] = l
        f_out["spec_zone_right"] = r
        audittrail( audit, "text", "spectral extraction",
                   "left point", f_out["spec_zone_left"],
                   "right point", f_out["spec_zone_right"],
                   "final size", get_si1_1d() )
    else:
        f_out["spec_zone_left"] = 1
        f_out["spec_zone_right"] = get_si1_1d()

    # %action% baseline_correction
    #   applies a baseline correction to the spectrum
    # %param% baseline boolean / default 1
    # %param% bcorr_algo enum offset linear spline quest polynomial moving_average / default "offset"
    #    offset : removes a automatically determined offset
    #    linear - spline : uses the list of user determined pivot points
    #    to define the baseline, then fit with a straight line or a spline
    #    quest : reconstruct the beginning of the FID using Linear Prediction technics, should usually be followed by offset
    #    polynomial - moving_average : uses statistics to separate signal from baseline,
    #       then apply a polynomial or moveing average correction.
    #    in more than one term is given, the correction are applied in sequence.
    # %param% spec_noise_n integer / default 30
    #   used by the offset algorithm to determine the offset to correct
    # %param% bcorr_pivots Integerlist / default "10 (get_si1_1d()-10)"
    #   pivot points used by the linear or spline algorithm
    # %param_cond% (bcorr_pivots > 0 and bcorr_pivots <= get_si1_1d())
    # %param% bcorr_radius integer / default 1
    #   radius around pivot points used by the linear or spline algorithm
    # %param% bcorr_order integer  / default 10
    # %param_cond% (bcorr_order > 0 and bcorr_order < get_si1_1d() and bcorr_order == 2*n)
    #   order (number of points corrected in the time domain) of the quest algorithm
    # %return% spec_offset
    #     offset computed when using the offset algo
    if (key_is_true(p_in,"baseline")):
        set_task("Baseline correction")
        toreal()
        pivots = "not used"
        bcorder = "not used"
        bcc = 1
        bc = p_in["bcorr_algo"]

        if ((bc.find("linear") + bc.find("spline")) != -2):    # if spline or linear is to be used
            l = p_in["bcorr_pivots"]
            n = 0
            pivots = l
            off = 0

            lst = l.split()

            for l in lst:
                n = n+1
                point_input(float(l))
                point_push()    # store pivots on the point stack
                off = off + val1d(int(l))

            off = off/n

        lst = bc.split()
        for bcalgo in lst:

            if (bcalgo == "offset"):
                spec_noise( p_in["spec_noise_n"] )
                addbase( get_shift())
                f_out["spec_offset"] = get_shift()

            elif (bcalgo == "linear"):
                # il y avait un %% en dernier argument
                bcorr(1, p_in["bcorr_radius"] )

            elif (bcalgo == "spline"):
                # il y avait un %% en dernier argument
                bcorr(2, p_in["bcorr_radius"] )

            elif (bcalgo == "quest"):
                bcorder = p_in["bcorr_order"]   
                bcorr_quest(bcorder)

            elif (bcalgo == "polynomial"):  # bcorrpx determines which algo is used
                if (get_debug()):
                    print("Baseline correction not available in DEBUG mode")
                else:
                    bcorrp0()
                    bcorr(3)

            elif (bcalgo == "moving_average"):
                if (get_debug()):
                    print("Baseline correction not available in DEBUG mode")
                else:
                    bcorrp1()
                    bcorr(3)
                    bcorrp0()

            else:
                raise "Error with baseline correction algorithm"

            audittrail( audit, "text", "Baseline correction -" + repr(bcc),
                       "algorithm", bcalgo,
                       "pivots", pivots,
                       "order", bcorder)
            bcc = bcc+1

    # %action% smoothing
    #   apply a smoothing filter to the data-set
    # %param% smoothing boolean / default 0
    # %param% smooth_w integer / default 5
    #    size of the smoothing window
    # %param% smooth_iteration / default 1
    #   number of loop
    # %param_cond% [ smooth_w > 0 ]
    if (key_is_true(p_in,"smoothing")):
        set_task("Smoothing filter applied")
        toreal()
        if (key_is_true(p_in,"smooth_iteration")):
            it = p_in["smooth_iteration"]
        else:
            it = 1
        for i in range(it):
            smooth(p_in["smooth_w"])
        audittrail( audit, "text", "Smoothing filter",
            "Number of points", p_in["smooth_w"],
            "Number of iterations", p_in["smooth_iteration"])


    # %action% median
    #   apply a median filter to the data-set
    # %param% median boolean / default 0
    # %param% median_w integer  / default 6
    #    size of the median filtering window
    # %param% median_i integer <= median_w  / default 3
    #    index of point to keep in the median filtering window
    # %param_cond% [ median_w > 0 ]
    # %param_cond% [ (median_i > 0) and (median_i <= median_w) ]
    if (key_is_true(p_in,"median")):
        set_task("Median filter applied")
        toreal()
        median(p_in["median_w"], p_in["median_i"]) 
        audittrail( audit, "text", "Median filter", "number of points",
                   p_in["median_w"], "index", p_in["median_i"])

    # %action% derivative
    #   compute the nth derivative of the data-set
    # %param% derivative boolean / default 0
    # %param% deriv_nth integer > 0 / default 1
    # %param% deriv_smooth integer > 0 / default 2
    #    binomial smoothing applied before derivative
    if (key_is_true(p_in,"derivative")):
        set_task("Dataset derivative computed")
        toreal()
        derivative( p_in["deriv_nth"], p_in["deriv_smooth"]) 
        audittrail( audit,  "text", "Derivative filter", "number of points",
                    p_in["deriv_nth"], "index", p_in["deriv_smooth"])

    # %action% select_state
    #   permit to choose the state complex / real of the output file
    #   complex data are changed to real by dropping the imaginary part
    #   real data are changed to complex by computing the Hilbert transform with tocomplex()
    #   this is usually not required as all processing commands prepare the state themselves 
    # %param% select_state boolean / default 0
    #   actuallly does the selection
    # %param% f1_state enum ignore complex real / default ignore
    #   force the f1 axis to real or complex. ignore will let the axis unchanged
    if(key_is_true(p_in,"select_state")):
        (t1,t2)=get_itype(2)
        if (p_in["f1_state"] == "complex"):
            if (t1==0):
                set_task("complex reconstruction")
            tocomplex("f1")
            p1s="complex"
        elif (p_in["f1_state"] == "real"):
            toreal("f1")
            p1s="real"

        if t1 == 1:
            t1s="complex"
        else:
            t1s="real"
        audittrail( audit,  "text", "Complex state of the data",
                        "axis", "was "+t1s+" now forced to "+p1s)

	# %action% spec_noise
    #   evaluate noise, estimated by finding an empty zone
    # %param% spec_noise_n integer / default 30
    #    number of different zones where noise is evaluated
    # %return% spec_std_noise
    #    estimate of the noise in the spectrum
    set_task("Noise evaluation")
    spec_noise( p_in["spec_noise_n"] )
    f_out["Noise"] = get_noise()


    # forces output data to be complex
    if (not key_is_true(p_in,"modulus")):   # HACK - compute imaginary parts to permit rephasing
        tocomplex()

    write_file_1d( audit, fileout )



#---------------------------------------------------------------------------
# write_file_1d
#---------------------------------------------------------------------------
def write_file_1d(  audit, fileout ):
    """
    write a 1D file and update the audittrail
    """
    if ( fileout != "memory" ):
        dim(1)
        writec(fileout)
        if ( get_itype_1d() == 1 ):
            sz = "Complex " + repr(get_si1_1d()/2)
        else:
            sz = "Real " + repr(get_si1_1d())
        audittrail( audit, "text", "output file created :", "file name", fileout,
                   "size", sz)

#---------------------------------------------------------------------------
# pp_1d
#---------------------------------------------------------------------------
def pp_1d(  audit, filein, filepeak, p_in, f_in, f_out ):
    """
This macro realizes the peak picking of a 1D spectrum
- spec_noise : evaluate noise, estimated by finding an empty zone
- prefilter  : smooths the spectrum before peak-picking, modification is not stored permanently
- restrict   : restricts the peak picking to a certain spectral zone 
- peakpick   : do the peak picking, by detecting local extrema 
- aggregate  : sorts peak list to aggregate peaks close from each other

    
%F1-input-domain%  frequency
%F1-output-domain% frequency
%dimensionality%   1

%author% Marc-Andre Delsuc
%version% 6.0
    """
    dim(1)
    if ( filein != "memory" ):
        read( filein )
        if ( get_dim() != 1 ):
            raise "This is not a 1D dataset"
        audittrail( audit, "text", "loading in memory", "filename", filein)
    else:
        audittrail( audit, "text", "using data in memory")

    if ( get_itype_1d() == 1 ):
        real()
        audittrail( audit, "text", "complex data-set, removing the imaginary part")

    size_ini = get_si1_1d()
    put("data")

    # %action% spec_noise
    #    evaluate noise, estimated by finding an empty zone
    # %param% spec_noise boolean / default 0
    # %param% spec_noise_n integer / default 30
    #    number of different zones where noise is evaluated
    # %return% NoiseInFrequencyDomain
    #    estimate of the noise in the spectrum
    set_task("Noise evaluation")
    if (p_in["spec_noise"]):
        spec_noise(p_in["spec_noise_n"])
        f_out["Noise"] = get_noise()
    else:
        noise(f_in["Noise"])
        f_out["Noise"] = get_noise()

    # %action% prefilter
    #    smooths the spectrum before peak-picking, modification is not stored permanently
    # %param% prefilter boolean / default 1
    # %param% prefilter_value float / default 3.0
    #    pre-smoothing parameter in Hz

    if (key_is_true(p_in,"prefilter")):
        set_task("Filtering before Peak-Picking")
        chsize (2*power2(get_si1_1d()-1))   # zerofill to next power of 2
        size_ft = get_si1_1d()
        iftbis()                            # fourier transform
        if  (p_in["peak_sign"] == "both"):
            contrast = 0.7
        else:
            contrast = 1.4             # peut-etre un peu fort
        gaussenh(p_in["prefilter_value"], contrast)
        sqsin(0)                # reduce trucation artifacts
        ftbis()         # restore smoothed data
        chsize(size_ini)
        audittrail( audit, "text", "Pre Filter applied to data-set",
                   "smoothing value in Hz (using gaussenh)",
                   p_in["prefilter_value"])

    # %action% restrict
    #     restricts the peak picking to a certain spectral zone 
    # %param% restrict boolean / default 0
    # %param% restrict_left float / default 10.0
    #   the left border of the restrict zone, in unit
    # %param% restrict_left_unit enum ppm hz index / default ppm
    #   the unit in which restrict_left is given
    # %param% restrict_right float / default 0.0
    #   the right border of the restrict zone, in unit
    # %param% restrict_right_unit enum ppm hz index / default ppm
    #   the unit in which restrict_right is given
    # %param_cond% (restrict_left{index} < restrict_right{index})
    # %return% restrict_left
    #     the left coordinate of the restrict zone in index
    # %return% restrict_right
    #     the right coordinate of the restrict zone in index
    if (key_is_true(p_in,"restrict")):
        set_task("Restricting peak-picking zone")
        if (p_in["restrict_left_unit"] == 'ppm'):
            l = (ptoi(p_in["restrict_left"],1,1))
        elif (p_in["restrict_left_unit"] == 'hz'):
            l = (htoi(p_in["restrict_left"],1,1))
        elif (p_in["restrict_left_unit"] == 'index'):
            l = p_in["restrict_left"]
        else:
            raise "Error with restrict_left_unit"

        if (p_in["restrict_right_unit"] == 'ppm'):
            r = (ptoi(p_in["restrict_right"],1,1))
        elif (p_in["restrict_right_unit"] == 'hz'):
            r = (htoi(p_in["restrict_right"],1,1))
        elif (p_in["restrict_right_unit"] == 'index'):
            r = p_in["restrict_right"]
        else:
            raise "Error with restrict_right_unit"

        l = max(1,round(l))
        r = min(get_si1_1d(),round(r))
        if (l > r):
            raise "Wrong restrict zone coordinates"

        if ((r-l) < 8):
            raise "restrict zone too small"

        zoom(1, l, r)
        f_out["restrict_left"] = l
        f_out["restrict_right"] = r
        audittrail( audit, "text", "spectral extraction",
                   "left point", f_out["restrict_left"],
                   "right point", f_out["restrict_right"])
    else:
        zoom(0)
        f_out["restrict_left"] = 1
        f_out["restrict_right"] = get_si1_1d()

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
        set_task("Peak-picking")
        com_max()
        mn = max(geta_max(1)/p_in["ratio_thresh"],p_in["noise_thresh"]*get_noise())
        minimax(mn, geta_max(1) + 1)
        print("+++ minimax " + repr(mn) + " " + repr(geta_max(1)+1))
        
        f_out["low_limit"] = mn
        if (p_in["peak_sign"] == "negative"):
            mult(-1)
        elif (p_in["peak_sign"] == "both"):
            com_abs()

        pkclear()
        peak()
        f_out["nb_detect_peaks"] = get_npk1d()
        print("+++ detected " + repr(get_npk1d()))
        get("data")
        pkreset()
        audittrail( audit, "text", "Peak Picking applied",
                   "low limit for peak detection", repr(mn), 
                   "number of peak detected", repr(get_npk1d()))

    # %action% aggregate
    #   sorts peak list to aggregate peaks close from each other
    # %param% aggregate boolean / default 0
    # %param% aggregate_value float / default 3.0
    #    distance parameter in Hz
    if (key_is_true(p_in,"aggregate")):
        audittrail( audit, "text", "aggregate - Reste a faire")

    # tab2pk_index(filepeak, "new")
    pkwrite_p(filepeak)

    audittrail( audit, "text", "Peak Picking file stored",
               "filename", filepeak)



#---------------------------------------------------------------------------
# ift_1d
#---------------------------------------------------------------------------
def ift_1d( audit, filein, filout, p_in, f_in, f_out):
    """
This macro Computes the Inverse Fourier transform of a 1D spectrum

- apodize       : standard apodisation
- inverseFourier : performs the inverse Fourier transform
- reverse       : reverses the spectral axis after FT
    
%F1-input-domain%  frequency
%F1-output-domain% time
%dimensionality%   1

%author% Marc-Andre Delsuc
%version% 6.0
    """
    dim(1)
    if ( filein != memory ):
        com_read( filein )
        if ( get_dim() != 1 ):
            raise "This is not a 1D dataset"
        audittrail( audit, "text", "loading in memory", "filename", filein)
    else:
        audittrail( audit, "text", "using data in memory")


    # %action% apodize
    #   standard apodisation
    # %param% apodize boolean / default 0
    # %param% apodisation string / default "sin 0.5"
    #   the string describing the apodisation to apply
    #   see apodise.g for details
    if ( key_is_true(p_in,"apodize") and (p_in["apodisation"] != "none")):
        apodise( p_in.raw("apodisation"))
        audittrail( audit, "text", "spectrum apodisation before iFT",
                   "function", p_in["apodisation"])

    # %action% inverseFourier
    #  performs the inverse Fourier transform
    # %param% inverseFourier boolean / default 1
    # %param% ask_for_iftsize boolean / default 0
    #  if set, the final size is determined by ift_size, automatically determined otherwise
    #   size for FT is determined by current size
    #   it will be extended to the next 2^n
    # %param% ift_size integer  / default 2*power2(get_si1_1d())
    # %param% ft_type enum none ift ifthilbert / default ifthilbert
    #    two iFT algorithms, ifthilbert discard the imaginary part if present and computes it using Hilbert transform
    # %param_cond%   [ ispower2(ift_size) ]
    # %return% size_for_ift contains the data size right after iFT
    if ( key_is_true(p_in,"inverseFourier")):
        initial = get_si1_1d()
        if ( p_in["ask_for_iftsize"] ):
            final = p_in["ift_size"]
        else:
            final = 2*power2(get_si1_1d()-2)
            if ( p_in["ift_type"] == "ifthilbert" and get_itype_1d() == 1 ):
                final = final / 2

        if ( p_in["ift_type"] == "ift" ):
            chsize(final)
            if ( get_itype_1d() == 0 ):
                iftbis()
            else:
                ift()
        elif( p_in["ift_type"] == "ifthilbert" ):
            if ( get_itype_1d() == 0 ):
                chsize(final)
                iftbis()
            else:
                chsize(2*final)
                real()
                iftbis()
        else:
            raise "error in ift type"

        f_out["size_for_ift"] = get_si1_1d()
        audittrail( audit, "text", "inverse Fourier transformation",
                   "iFT type", p_in["ift_type"],
                   "initial size", initial,
                   "final size", f_out["size_for_ift"])

    # %action% reverse
    #   reverses the spectral axis after FT
    # %param% reverse boolean / default 0
    if ( key_is_true(p_in,"reverse") ):
        reverse()
        if (get_itype_1d()==1):
            invf()       # if complex, inverse imaginary part
        auditrail("text", "Reverse time axis")

    write_file_1d( audit, fileout )


#---------------------------------------------------------------------------
# integ_1d
#---------------------------------------------------------------------------
def integ_1d(  audit, filein, integratout, p_in, f_in, f_out):
    """
This macro Computes integrales of a 1D spectrum
- integral    : omputes integral positions from peaks

    
%F1-input-domain%  frequency
%F1-output-domain% frequency
%dimensionality%   1

%author% Marc-Andre Delsuc
%version% 6.0
    """
    dim(1)
    if ( filein != "memory" ):
        read( filein )
        if ( get_dim() != 1 ):
            raise "This is not a 1D dataset"
        audittrail( audit, "text", "loading in memory", "filename", filein)
    else:
        audittrail( audit, "text", "using data in memory")

    print("+++ ITYPE = " + repr(get_itype_1d()))
    if ( get_itype_1d() == 1):
        real()
        audittrail( audit, "text", "complex data-set, removing the imaginary part")

    # %action% integral
    #  computes integral positions from peaks
    # %param% integ_aggregate float / default 30
    #   peak aggregatian width, in Hz
    # %param% integ_extention float / default 15
    #    integral extension, in Hz
    # %return% nb_integral
    #    number of integral zones determined
    set_task("integration")
    inte = 0

    exi = p_in["integ_extention"]*get_si1_1d() / get_specw_1d()
    agi = p_in["integ_aggregate"]*get_si1_1d() / get_specw_1d()

    left = {}
    right = {}
    if ( get_npk1d() == 1 ):
        inte = 1
        left[1] = max(1, geta_pk1d_f(1) - exi)
        right[1] = min( get_si1_1d(), geta_pk1d_f(1) + exi)
    elif ( get_npk1d() > 1 ):
        inte = 0
        left[1] = max(1, geta_pk1d_f(1) - exi)
        i = 2
        while ( i <= get_npk1d() ):
            if ( (geta_pk1d_f(i) - geta_pk1d_f(i-1)) > agi):
                inte = inte+1
                right[inte] = min( get_si1_1d(), geta_pk1d_f(i-1) + exi)
                left[inte+1] =  max(1, geta_pk1d_f(i) - exi)
            i = i+1
        inte = inte+1
        right[inte] = min( get_si1_1d(), geta_pk1d_f(get_npk1d()) + exi)

    f_out["nb_integral"] = inte
    audittrail( audit, "text", "Standard aggregation applied",
               "aggregate width", p_in["integ_aggregate"],
               "extension width", p_in["integ_extention"],
               "number of integral detected", inte)


    # sortie des integrales
    # dbopen $fileinte dbint
    # dbint["description"] = "beginning ending offset value calibration"
    off = {}
    for i in range(1,inte+1):
        left[i] = int(round(left[i]))
        right[i] = int(round(right[i]))
        off[i] = 0.5*(val1d(left[i]) + val1d(right[i]))
    put("data")
    int1d()

    fout = open(integratout, "w")
    for i in range(1,inte+1):
        integral = val1d(right[i]) - val1d(left[i]) - off[i]*float(right[i]-left[i])
        if ( i == 1 ):
            calibration = integral
        # set dbint[$i] = (round($left[$i]) ; round($right[$i]) ; $off[$i] ; $integral ; $integral/$calibration)
        fout.write( repr(i) + "=" +
                    repr(left[i]) + " " + repr(right[i]) + " " +
                    repr(off[i]) + " " + repr(integral) + " " + repr(integral/calibration))
        fout.write("\n")

    fout.close()
    # dbclose dbint
    get("data")

    audittrail( audit, "text", "Integral file stored", "filename", integratout)


#---------------------------------------------------------------------------
# maxent_1d
#---------------------------------------------------------------------------
def maxent_1d( audit, filein, fileout, p_in, f_in, f_out):
    """
This macro realizes the MaxEnt analysis of a 1D FID
- freq_massage  : computes a temporary Fourier transform
- truncate   : truncates the FID by removing the last points
- preconvoluate  :apply a preconvolution before analysis, this may help stability of the algorithm and enhance noise rejection
- partialsampling   : set-up for processing data partially sampled in the time domain
- deconvoluate : apply a deconvolution during analysis,
- maxent : apply MaxEnt analysis,

    
%F1-input-domain%  frequency
%F1-output-domain% frequency
%dimensionality%   1

%author% Marc-Andre Delsuc
%version% 6.0
    """
    dim(1)
    if ( filein != "memory"):
        read( filein )
        if ( get_dim() != 1 ):
            raise "This is not a 1D dataset"
        audittrail( audit, "text", "loading in memory", "filename", filein)
    else:
        audittrail( audit, "text", "using data in memory")

    # actions in the frequency domain
    # %action% freq_massage
    #  computes a temporary Fourier transform
    # %param% freq_massage boolean / default 1
    # %param% ft_type enum none ft_seq ft_sim ft rft / default ft
    #    the Fourier transform algorithm, depends on spectrometer and acquisition scheme
    # %param% phase boolean / default 0
    # %param% phase_0 float / default 0.0
    #  global order phase correction
    # %param% phase_1 float / default 0.0
    #  1st order phase correction
    if ( key_is_true(p_in,"freq_massage") ):
        final = 4*power2(get_si1_1d()-2)
        initial = get_si1_1d()
        chsize(final)
        if ( p_in["ft_type"] != "none" ):
            if ( p_in["ft_type"] == "ft" ):
                ft()
            elif ( p_in["ft_type"] == "rft" ):
                rft()
            elif ( p_in["ft_type"] == "ft_sim" ):
                ft_sim()
            elif ( p_in["ft_type"] == "ft_seq" ):
                ft_seq()
            else:
                pass
                # faut evaluer la macro p_in["ft_type"]
            f_out["size_for_ft"] = final
            audittrail( audit, "text", "temporary Fourier transformation,",
                       "FT type", p_in["ft_type"],
                       "initial size", initial,
                       "final size", f_out["size_for_ft"])
            if ( key_is_true(p_in,"phase") ):
                phase(p_in["phase_0"], p_in["phase_1"])
                audittrail( audit,  "text", "phase correction applied",
                            "0th order", p_in["phase_0"],
                            "1th order", p_in["phase_1"])

            real()
            iftbis()
            chsize(initial)
            f_out["mock_fid_size"] = initial
            audittrail( audit, "text", "Back to time domain",
                       "final size ", initial)

    # %action% truncate
    #   truncates the FID by removing the last points
    # %param% truncate boolean / default 0
    # %param% trunc_size integer / default get_si1_1d()
    #   FID size after truncation
    # %param_cond% [ (trunc_size > 0 ) and (trunc_size <= get_si1_1d() ) ]
    if ( key_is_true(p_in,"truncate")):
        s = int(p_in["trunc_size"])
        if ( s <= 0 ):
            raise "Cannot truncate to 0"
        if ( s > get_si1_1d() ):
            raise "Error with truncation size (larger than actual size)"
        initial = get_si1_1d()
        chsize( min(get_si1_1d(), s))
        audittrail( audit, "text", "FID truncation before FT", "initial_size", initial, "final size", get_si1_1d())

    # %action% preconvoluate
    #   apply a preconvolution before analysis, this may help stability of the algorithm and enhance noise rejection
    # %param% preconvoluate boolean / default 0
    # %param% preconvolution string / default "expbroad(1)"
    #   the string describing the apodisation to apply
    #   see apodise.g for details
    if (key_is_true(p_in,"preconvoluate") and (p_in["preconvolution"] != "none")):
        apodise(p_in.raw("preconvolution"))
        audittrail( audit, "text", "FID preconvolution before MaxEnt analysis",
                   "function", p_in["preconvolution"])


    put("data")
    noise(f_out["Noise"])
#    window_mode(1)
#    window("F1", 0)
#    window("F2", 0)

    audittrail( audit, "text", "Mock FID stored",
               "fid size", get_si1_1d(),
               "noise value", get_noise())


    # %action% partialsampling
    #   set-up for processing data partially sampled in the time domain
    # %param% partialsampling boolean / default 0
    # %param% samplingfile string / default vplist
    #   the filename of the file containing the list of sample to use,
    #   file is relative to the directory containg the input file 
    #   see sampling.g for details

    if (key_is_true(p_in,"partialsampling")):
        sampf = location + os.sep + p_in["samplingfile"]
        sampling(sampf)
        sign = get_returned()
        partial = (sign / get_si1_1d())
        audittrail( audit, "text", "Set-up for partial sampling during MaxEnt analysis",
                   "sampling file", sampf,
                   "number of sampled points", sign,
                   "mean density of sampling", (str(100 * partiall) + "%"))
    else:
        partial = 1.0
#        window_mode(1)
#        partial = 1.0
#        window(0)       # reset window if no partial sampling

    # %action% deconvoluate
    #   apply a deconvolution during analysis,
    # %param% deconvoluate boolean / default 1
    # %param% deconvolution string / default "expbroad(1)"
    #   the string describing the apodisation to apply
    #   see apodise.g for details

    if (key_is_true(p_in,"deconvoluate") and (p_in["deconvolution"] != "none")):
        itype(0)
        one()
        itype(1)
        apodise(p_in.raw("deconvolution"))
        put("filter")  # populate filter buffer with the deconvolution function
        audittrail( audit,  "text", "deconvolution set-up for MaxEnt analysis",
                    "function", p_in["deconvolution"])
        com_filter(1)
        sumcons(0)      # MAD - devrait etre filter 2 - mais ca marche pas ?????
        nchannel(1)
    else:
        com_filter(0)
        nchannel(1)
        audittrail( audit, "text", "No deconvolution set-up for MaxEnt analysis")

    # %action% maxent
    #   apply MaxEnt analysis,
    # %param% maxent boolean / default 1
    # %param% me_preset boolean / default 0
    #   determines convergence parameters from a given preset or gives the details
    # %param% preset_value enum 1 2 3 4 5 / default 3
    #   preset for maxent; 1 is faster, 5 is more accurate
    # %param% iteration integer / default 100
    #   global number of iterations - set by preset
    # %param% control integer / default 10
    #   number of iterations between control step - set by preset
    # %param% lambsp float / default 5.0
    #   the value of the lambda controlling parameter LAMBSP in the Gifa algorithm
    # %param% lambcont enum 0 1 2 / default 1
    #   the value of the lambda controlling parameter LAMBCONT in the Gifa algorithm
    # %param% me_size integer / default 2*get_si1_1d()
    #   size for MaxEnt reconstruction
    # %param% noise_weight float / default 1.0
    #   the noise used during the MaxEnt iteration is weigthed with this scalar,
    #   this permits to compensate for over- or under-estimate of the noise level
    # %param_cond%   [ noise_weight > 0.0 ]
    # %param_cond%   [ ispower2(me_size) ]

    if (key_is_true(p_in,"maxent")):
        if (key_is_true(p_in,"me_preset")):
            me_preset( p_in["preset_value"])
            audittrail( audit,  "text", "Preset for MaxEnt Analysis"
                        "using preset level", p_in["preset_value"])
        else:
            me_preset(3)           
            com_iter(int(p_in["iteration"]))
            ndisp(int(p_in["control"]))
            lambsp(float(p_in["lambsp"]))
            lambcont(int(p_in["lambcont"]))

        noise (get_noise()*p_in["noise_weight"]*partial)
        audittrail( audit, "text", "Parameters for MaxEnt Analysis",
                   "size of reconstruction", p_in["me_size"],
                   "weighted noise", get_noise(),
                   "Algorithm flags", get_algo(),
                   "lambda control", get_lambcont(),
                   "lambda multiplier", get_lambsp(),
                   "global iteration", get_iter(),
                   "control", get_ndisp(),
                   "line minimisation iteration", get_miniter())
        get("data")
        sz=int(p_in["me_size"])
        put("data")
        maxent(sz)
        audittrail.g("text", "MaxEnt Analysis applied",
                     "number of iteration performed", get_iterdone(),
                     "final chi square", getchi2(),
                     "final entropy", get_entropy(),
                     "final lambda", get_lambda(),
                     "final convergence criterium", get_convergence())
        write_file_1d( audit, fileout)
#---------------------------------------------------------------------------
# post_maxent_1d
#---------------------------------------------------------------------------
def post_maxent_1d(  audit, filein, fileout, p_in, f_in, f_out):
    """
    post processing of a 1D spectrum processed by MaxEnt
    """
    dim(1)
    if ( filein != "memory" ):
        read( filein )
        if ( get_dim() != 1 ):
            raise "This is not a 1D dataset"
        audittrail( audit, "text", "loading in memory", "filename", filein)
    else:
        audittrail( audit, "text", "using data in memory")

    # %action% calibration
    #   calibrate 0ppm on the spectrum
    # %param% calibration float / default 0.0
    #   the coordinate in ppm of the right-most point in the spectrum
    offset (p_in["calibration"]*get_freq_1d())
    audittrail( audit, "text", "calibration offset applied", "calibration", p_in["calibration"])


    # %action% spectral_zone
    #   extract one spectral zone of the spectrum
    # %param% spectral_zone boolean / default 0
    # %param% spec_zone_left float / default 10.0
    #   the left border of the extract zone, in unit
    # %param% spec_zone_left_unit enum ppm hz index / default ppm
    #   the unit in which spec_zone_left is given
    # %param% spec_zone_right float / default 0.0
    #   the right border of the extract zone, in unit
    # %param% spec_zone_right_unit enum ppm hz index / default ppm
    #   the unit in which spec_zone_right is given
    # %param_cond% (spec_zone_left{index} < spec_zone_right{index})
    # %return% spec_zone_left  the left coordinate of the extracted spectral zone in index
    # %return% spec_zone_right  the right coordinate of the extracted spectral zone in index
    if (key_is_true(p_in,"spectral_zone")):
        if (p_in["spec_zone_left_unit"] == 'ppm'):
            l = ptoi(p_in["spec_zone_left"],1,1)
        elif (p_in["spec_zone_left_unit"] == 'hz'):
            l = htoi(p_in["spec_zone_left"],1,1)
        elif (p_in["spec_zone_left_unit"] == 'index'):
            l = p_in["spec_zone_left"]
        else:
            raise "Error with spec_zone_left_unit"

        if (p_in["spec_zone_right_unit"] == 'ppm'):
            r = ptoi(p_in["spec_zone_right"],1,1)
        elif (p_in["spec_zone_right_unit"] == 'hz'):
            r = htoi(p_in["spec_zone_right"],1,1)
        elif (p_in["spec_zone_right_unit"] == 'index'):
            r = p_in["spec_zone_right"]
        else:
            raise "Error with spec_zone_right_unit"

        l = (max(1,round(l)))
        r = (min(get_si1_1d(),round(r)))
        if (l > r):
            raise "Wrong spectral zone coordinates"

        if ((r-l) < 8):
            raise "spectral zone too small"

        if (l != 1 or r != get_si1_1d()):
            extract(l,r)

        f_out["spec_zone_left"] = l
        f_out["spec_zone_right"] = r
        audittrail( audit, "text", "spectral extraction",
                   "left point", f_out["spec_zone_left"],
                   "right point", f_out["spec_zone_right"],
                   "final size", get_si1_1d ())
    else:
        f_out["spec_zone_left"] = 1
        f_out["spec_zone_right"] = get_si1_1d()


    # %action% baseline correction
    #   baseline correction of the spectrum
    # %param% baseline boolean / default 1
    # %param% bcorr_algo enum offset linear spline polynomial moving_average polynomial+moving_average moving_average+polynomial / default linear
    # %param% bcorr_pivots Integerlist / default 10.0
    #   pivot points used by the linear or spline algorithm
    # %param_cond% (index > 0 and index <= get_si1_1d())
    # %param% bcorr_radius integer > 0 / default 3
    #   radius around pivot points used by the linear or spline algorithm
    # %return% spec_offset
    #     offset computed from the pivot points when using the offset algo
    if (key_is_true(p_in,"baseline")):
        if (get_debug()):
            print("Baseline correction not available in DEBUG mode")
        else:
            toreal()
            pivots = "not used"
            if (p_in["bcorr_algo"] == "linear" or p_in["bcorr_algo"] == "spline"):
                l = p_in["bcorr_pivots"]
                n = 0
                pivots = l
                off = 0
                lst = l.split()
                for l in lst:
                    n = n+1
                    point_input(int(l))
                    point_push()    # store pivots on the point stack
                    off = off + val1d(int(l))

                off = off/n

            if (p_in["bcorr_algo"] == "offset"):
                addbase(off)
                f_out["spec_offset"] = off
            elif (p_in["bcorr_algo"] == "linear"):
                bcorr( 1, p_in["bcorr_radius"],  "%%")
            elif (p_in["bcorr_algo"] == "spline"):
                bcorr( 2, p_in["bcorr_radius"],  "%%")
            elif (p_in["bcorr_algo"] == "polynomial"):  # bcorrpx determines which algo is used
                bcorrp0()
                bcorr(3)
            elif (p_in["bcorr_algo"] == "moving_average"):
                bcorrp1()
                bcorr(3)
                bcorrp0()
            elif (p_in["bcorr_algo"] == "polynomial+moving_average"):   # these 2 concatenate 2 processing
                bcorrp0()
                bcorr(3)
                bcorrp1()
                bcorr(3)
                bcorrp0()
            elif (p_in["bcorr_algo"] == "moving_average+polynomial"):
                bcorrp1()
                bcorr(3)
                bcorrp0()
                bcorr(3)
            else:
                raise "Error with baseline correction algorithm"

            audittrail( audit, "text", "Baseline correction",
                       "algorithm", p_in["bcorr_algo"],
                       "pivots", pivots)

    write_file_1d( audit, fileout)
    


