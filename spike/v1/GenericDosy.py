#!/usr/bin/env python 
# encoding: utf-8

"""
The library of DOSY processing functions to be used with the NPK program.

This library implement the specific functions needed for NMR-DOSY processing.
Most of these functions require that the NPK mathematical kernel is loaded.

"""

from __future__ import print_function

__author__ = "Marc A. Delsuc <delsuc@igbmc.fr>"
__date__ = "Oct 2009"

import os.path
import sys
import copy
from Generic import *

##################################################################################
def dosybidon(datafile, Dmin, Dmax, Dfactor, me, list_to_do):
    if (get_noise() == 0.0):
        raise 'NOISE should not be 0.0'

# argument check

# MAD
    (here,file) = os.path.split(datafile)
    join(datafile)
    if (me.get_ilttype() == "Tabulated" and get_si_tab() != get_c_sizef1() ) :
        disjoin()
        raise """
SIZE MISMATCH
F1 dimension of input datafile
does not match gradient table size.
Please verify that F1 size = gradient table size"""

# compute sizes to be used before and after (szbf - szaf)
#MAD A REFAIRE szaf = dosysize suffit !
    szbf = get_c_sizef1()
    if (me.get_iltsize() > szbf):
        szaf = me.get_iltsize()
    else:
        szaf = me.get_iltsize()
    dim(2)
    chsize(szaf,get_si2_2d())
    disjoin()


# Count the number of columns to process
    countcol = len(list_to_do)
    if (countcol == 0):
        raise "No columns selected"

# setup the inprogress
    print(('%i columns are to process' % countcol))
    peakcount = 1

#?MAD?
    result = open(os.path.join(here,'dosy_processing.log'),"w")
    result.write("# processing bidon\n")
# Do the processing
    for i in list_to_do:
        print(("Processing col %i ... "%(i)))
        st = ("col %i iter. count %i chi2 %f\\n"%(i, 1000, 12.34))
        result.write(st)

    dim(2)
    one()
    sqsin(0.5,"F1")
    sqsin(0.5,"F1")
    sqsin(0.5,"F2")
    result.close()

	# save dosy parameters
    dim(2)
    dfactor(Dfactor)
    dmin(Dmin)
    dmax(Dmax)
    # fake log(D) in ppm
    freq(get_c_freq(), 1, get_c_freq2())
    specw( math.log(Dmax)-math.log(Dmin), get_c_specwf2()) 
    offset(math.log(Dmax), get_c_offsf2())


################################################################################
def dosy2d(datafile, Dmin, Dmax, Dfactor, me, list_to_do):
    """Realize the actual processing of a 2D DOSY

 It implements the complete DOSY processing for 2D.
 DOSY data are in-memory, diffusion dimension is in F1, NMR spectra in F2

 datafile : name and path of the specified input F2 processed data file.
 Dmin, Dmax, Dfactor determines the window border along the diffusion axis
 me : is a maxent_state object, holding all the details for the processing parameters
 list_to_do : list of the indexes of the column to process.

 %author% Marc-Andre Delsuc
 %version% 1.0
 mar 2006 - first version, rewritte from dosy2d in Gifa

    """
    if (get_noise() == 0.0):
        raise 'NOISE should not be 0.0'
    dim(1)
#    window(0)   # just in case......... MAD
    dim(2)

# MAD
    (here,file) = os.path.split(datafile)
    print(datafile) 
    join(datafile)
    if (me.get_ilttype() == "Tabulated" and get_si_tab() != get_c_sizef1() ) :
        disjoin()
        print("file size: ",get_c_sizef1(),"  tabulated value size : ",get_si_tab())
        raise """
SIZE MISMATCH
F1 dimension of input datafile
does not match gradient table size.
Please verify that F1 size = gradient table size"""

# compute sizes to be used before and after (szbf - szaf)
#MAD A REFAIRE szaf = dosysize suffit !
    szbf = get_c_sizef1()
    if (me.get_iltsize() > szbf):
        szaf = me.get_iltsize()
    else:
        szaf = me.get_iltsize()
    dim(2)
    chsize(szaf,get_si2_2d())
    zero()
    disjoin()


# Count the number of columns to process
    countcol = len(list_to_do)
    if (countcol == 0):
        raise "No columns selected"

# set up avancement file (text file that contains the current number of processed columns during the process)
    avncmt=open(os.path.join(here,'avancement'),"w")
    avncmt.write("1 / %i"%countcol)
    avncmt.close()

#MAD
#message in DOSYm mode
#    if (exist('is_dosym')):
#        topspin_out("starting processing for 1 / %i columns"%countcol)


# setup the inprogress
    print(('%i columns are to process' % countcol))
    peakcount = 0

# fake a regular mode in FIT
#if (($PROC s= 'FIT' | $proc s= 'FIT_2_COMP') & $me.get_ilttype() s= 'REGULAR') :
#  dim 2 col 1 dim 1
#  one tm $si1_1d $si1_1d addbase (1/$si1_1d)
#  mult (0.5/$SPECW_1D)
#  put tab
#endif

# do the processing of selected columns found in the dosydbfile
    dim(2)
    zoom(0)
    #scale 1
    #unit(i)  # for some reason, it does not work unless unit i is chosen
    #unit_x(i)

#?MAD?
    result = open(os.path.join(here,'dosy_processing.log'),"w")

# Do the processing
    join(datafile)
    for il in list_to_do:
        i=int(il)
#        print("Processing col %i   %i / %i"%(i,peakcount+1,len(list_to_do)))
        set_task("Processing col %i   %i / %i"%(i,peakcount+1,len(list_to_do)))
        dim(1)
        chsize(szbf)
        getc("F1",i,1,szbf)
        dfactor(Dfactor)
        dmin(Dmin)
        dmax(Dmax)

        if (me.get_iltalgo() == 'MaxEnt'):
            st=dosymaxent(me)  # this one is a macro
        elif (me.get_iltalgo() == 'Fit'):
            st=dosyfit()
        else:
            raise "error in iltalgo"

        dim(2)
        put("col",i)
        print((("column : %i  "+st)%(i)))
        result.write(("column : %i  "+st+"\n")%(i))
        peakcount=peakcount+1
        avncmt=open(os.path.join(here,'avancement'),"w")
        avncmt.write("%i / %i"%(peakcount,countcol))
        avncmt.close()

    result.close()

	# save dosy parameters
    dim(2)
    dfactor(Dfactor)
    dmin(Dmin)
    dmax(Dmax)
    # fake log(D) in ppm
    freq(get_c_freq(), 1, get_c_freq2())
    specw( math.log(Dmax)-math.log(Dmin), get_c_specwf2()) 
    offset(math.log(Dmax), get_c_offsf2())

    print(("Finally process %i columns :"%(peakcount)))

    if os.path.exists('avancement'):
        os.remove('avancement')
#    if (exist('is_dosym')) topspin_out ("processing finished")


#---------------------------------------
def dosyfit(me):
    raise "Reste a faire ! - helas - "
"""
					dosyfit()
					result.write("col ",i," chi2_1_comp ",get_chi2())
# a teste
					set p2_1 = $p2
					set p1_1 = $p1
					set dp2_1 = $dp2
					set dp1_1 = $dp1
					set chi2_1 = get_chi2()
					if ($proc s= 'FIT_2_COMP' & get_chi2() > ($si1_1d)) :
						dosyfit_2()
					else
						set p3 = 0.0
					endif
# MAD ICI
					if ( get_chi2() < (0.8*$chi2_1) & $p2>0 & $p3>0 ) :
						chsize $szaf
						result.write("chi2_2_comp"#get_chi2())
						result.write('damping 1'#$p2#'+/-'#$dp2#'amplitude 1'#$p1#'+/-'#$dp1)
						result.write('damping 2'#$p4#'+/-'#$dp4#'amplitude 2'#$p3#'+/-'#$dp3)
						set left1 = (dtoi(max(0,$p2-$dp2),1,1))
						set right1 = (dtoi($p2+$dp2,1,1))
						set left2 = (dtoi(max(0,$p4-$dp4),1,1))
						set right2 = (dtoi($p4+$dp4,1,1))
						chsize (%*2)
						zero
						if ($p2 < $dmin | $p2 > $dmax) \
							print "Warning, outside damping limits"
						if ($p4 < $dmin | $p4 > $dmax) \
							print "Warning, outside damping limits"
						unit i   # for some reason, it does not work unless unit i is chosen
						simun $p1 (dtoi($p2,1,1)) 0 0
						gm (2*($right1-$left1)*$SPECW_1D/$SI1_1D)
						put("data")
						zero()
						simun $p3 (dtoi($p4,1,1)) 0 0
						gm (2*($right2-$left2)*$SPECW_1D/$SI1_1D)
						adddata()
						tm(1,1)				#making an apodisation to set somes wiggles positive
						ft()
						#refmacro 1
						real()
						#refmacro 0
					else:
						result.write('damping'#$p2_1#'+/-'#$dp2_1#'amplitu'#$p1_1#'+/-'#$dp1_1)
						if ($p2_1 <0) :
							fprint $result 'negative damping, setting to zero'
							dim 1 chsize $szaf zero
						else
							if ($p2_1 < $dmin | $p2_1 > $dmax) \
							print "Warning, outside damping limits"
							chsize $szaf
							set left = (dtoi(max(0,$p2_1-$dp2_1),1,1))
							set right = (dtoi($p2_1+$dp2_1,1,1))
							chsize (%*2)
							zero simun $p1_1 (dtoi($p2_1,1,1)) 0 0
							gm (2*($right-$left)*$SPECW_1D/$SI1_1D)
							tm 1 1				#making an apodisation to set somes wiggles positive
							ft
							#refmacro 1
							real
							#refmacro 0
						endif
					endif

				chsize $szaf
				set peakcount = (%+1)

				# now feed the DOSY GTB file
				if ($proc s= 'MAXENT') :
					set thenoise = $noise
					set thechi2 = $chi2min
					set thentropy = get_entropy()
					set theiterdone = $nbreiter
					set thelambsp = $lambsp
					set theminiter = $miniter
					set thelambcont = $lambcont
					set thedfactor = $dfactor
					set thedmin = $dmin
					set thedmax = $dmax
				else
					set thenoise = $noise
					set thechi2 = get_chi2()
					set thentropy = 'NoValue'
					set theiterdone = $iter
					set thelambsp = 'NoValue'
					set theminiter = $miniter
					set thelambcont = 'NoValue'
					set thedfactor = $dfactor
					set thedmin = $dmin
					set thedmax = $dmax
				endif
				if ($pflag == 1) :	# for not already processed columns
					#print ($thenoise#$thechi2#$thentropy#$theiterdone#$thelambsp#$theminiter#$thelambcont#$thedfactor#$thedmin#$thedmax)
					set dosy[$i] = (2#$thenoise#$thechi2#$thentropy#$theiterdone#$thelambsp#$theminiter#$thelambcont#$thedfactor#$thedmin#$thedmax)
					dim 2 put col $i
					set pcol = ($pcol+1)
				elsif ($pflag == 3) :	# for already processed columns
					set tmp = (tail($dosy[$i]))
					set oldchi2 = (head(tail($tmp)))
					print ($oldchi2#$thechi2)
					if ($oldchi2 > $thechi2) : # save only if it's better
						set dosy[$i] = (2#$thenoise#$thechi2#$thentropy#$theiterdone#$thelambsp#$theminiter#$thelambcont#$thedfactor#$thedmin#$thedmax)
						dim 2 put col $i
						set srpcol = ($srpcol+1)
					else
						set dosy[$i] = (2#(tail($dosy[$i])))
						dim 2
						set rpcol = ($rpcol+1)
					endif
				endif
				inprogress $peakcount
				open $avncmt
				fprint $avncmt ($peakcount#'/'#$countcol)
				close $avncmt
#				if (exist('is_dosym')) topspin_out ("processing"# $peakcount# '/'# $countcol)

			endif
		endif
=loop
	endfor
"""

def dosymaxent(mearg):
    """
Realise the 1D ILT transformation involved in DOSY processing.
uses the 1D buffer and the tab buffer as found in the kernel.

mearg is a MaxEnt object, which contains all the parameters
It try to optimise lambsp in order to get the fastest iteration possible

It stores temporary results in a file, and get back to it
if something went wrong.

see also : dosy2d INVTLAP INVTLAPCONT

 %author% M-A.D. from Thierry Gostan
    """
    temp = NPKtempfile('.gs1')
    bavard=(1==0)		# mettre 1==1 pour avoir des messages

    me=copy.copy(mearg)

    itermax = me.get_iter()			#target nmbr of iterations
    me.set_iter(me.get_ndisp())               #modifie la valeur du nombre d'iter avec celle du premier affichage

# first iteration is special, try to have at least something working
    do_maxent_ilt(me)
    iterfait = get_iterdone()
    diverged = check_diverged()    # check id ok (entropy and chi2 != NaN)
    nbreiter= get_iterdone()

    while ((diverged) and (iterfait < itermax) ): # searching as much as possible
        me.set_lambsp (me.get_lambsp() / 5)
        lambsp(me.get_lambsp())
        if bavard:
           print("initial divergence - new lambsp : ", me.get_lambsp())
        get("data")
        cont_maxent_ilt(me)
        iterfait = (get_iterdone())

    if check_diverged():    # then have no iteration left, and it never converged !
        chi2min = 0.0
        entropymin = 1
        zero()  # and erase
        st = "iter. count %i chi2 -failed-"%(iterfait)
    else:   # use this as a starting point
        chi2min = get_chi2()
        entropymin = get_entropy()
        st = "iter. count %i chi2 %f"%(iterfait, chi2min)

    writec(temp)

# following iterations  - will be executed only if 1st try converged and me.get_ndisp() < me.get_iter()
    compteur = 0    # count how many time it stalled
    if bavard:
           print("chi2min", chi2min)
    while ((iterfait < itermax) ):   # iterations left to do, and not stalled
        me.set_iter(me.get_iter() + me.get_ndisp())
        iter(me.get_iter())
        chi2old = get_chi2()
        cont_maxent_ilt(me)
        iterfait = (get_iterdone())      # MAD iterdone

        while (check_diverged() and (iterfait<itermax)):      # Pas un nombre et plus d'iter dispo
            read(temp)  # restart from previous
            me.set_lambsp(me.get_lambsp() / 1.2)
            lambsp(me.get_lambsp())
            if bavard:
              print("Diverged - new lambsp : ",me.get_lambsp())
            me.set_iter(me.get_iter() + me.get_ndisp())
            iter(me.get_iter())
            cont_maxent_ilt(me)
            iterfait = (get_iterdone())      # MAD iterdone

            
        if (not check_diverged()):       # if converging
            if (get_chi2() < chi2min):  # if better, (as it should), keep it
                writec(temp)
                chi2min = get_chi2()
                entropymin = get_entropy()
                st = ("iter. count %i chi2 %f"%(iterfait, chi2min))
                nbreiter = me.get_iter()
#                compteur=0
            if ((math.fabs(get_chi2() - chi2old)) < (0.0001 * get_chi2())):   # if not changing by more than 0.0001 == stalled
                compteur = (compteur + 1)
                if bavard:
                   print("Stalled")
            else:   # going up
                compteur = 0

        if (compteur == 2):     #if starting to stall, try to move it 
                me.set_lambsp((me.get_lambsp() * 2))
                lambsp(me.get_lambsp())
                if (me.get_miniter() > 2.5):
                    me.set_miniter(me.get_miniter() / 1.2)
                if bavard:
                    print("Stalled - new me_int_iteration : ",me.get_miniter()," - new labsp : ",me.get_lambsp())

        if (compteur == 4):
                iterfait = itermax
    read(temp)
    try:                    # try removing the temp file, as it may not diseappear
            os.unlink(temp)
    except:
            pass
    return st
##########################################################################################
def dosymaxent_new(mearg):
    """
This version is a rewrite, it does not really work for the moment

Realise the 1D ILT transformation inv-olved in DOSY processing.
uses the 1D buffer and the tab buffer as found in the kernel.

mearg is a MaxEnt object, which contains all the parameters
It try to optimise lambsp in order to get the fastest iteration possible

It stores temporary results in a file, and get back to it
if something went wrong.

see also : dosy2d INVTLAP INVTLAPCONT

 %author% M-A.D. from Thierry Gostan
    """
    stalled = 4 # stalled if 4 major iterations found the prgm not changing
    temp = NPKtempfile('.gs1')

    me=copy.copy(mearg)

    itermax = me.get_iter()			#target nmbr of iterations
    me.set_iter(me.get_ndisp())               #modifie la valeur du nombre d'iter avec celle du premier affichage

# first iteration is special, try to have at least something working
    do_maxent_ilt(me)
    iterfait = get_iterdone()
    diverged = check_diverged()    # check id ok (entropy and chi2 != NaN)

    while ((diverged) and (iterfait < itermax) ): # searching as much as possible
        me.set_lambsp (me.get_lambsp() / 5)
        lambsp(me.get_lambsp())
        print("initial divergence - new lambsp : ", me.get_lambsp())
        get("data")
        cont_maxent_ilt(me)
        iterfait = (get_iterdone())

    if check_diverged():    # then have no iteration left, and it never converged !
        chi2min = 0.0
        entropymin = 1
        zero()  # and erase
        st = "iter. count %i chi2 -failed-"%(iterfait)
    else:   # use this as a starting point
        chi2min = get_chi2()
        entropymin = get_entropy()
        st = "iter. count %i chi2 %f"%(iterfait, chi2min)

    writec(temp)

# following iterations  - will be executed only if 1st try converged and me.get_ndisp() < me.get_iter()
    compteur = 0    # cout how many time it stalled
    compt_moved = 0    # count how time lambda was modified for stalling
    print("chi2min", chi2min)
    while ((iterfait < itermax) and (compteur <stalled) and (chi2min>1.0 and chi2min!=0.0)):   # iterations left to do, and not stalled
        me.set_iter(me.get_iter() + me.get_ndisp())
        iter(me.get_iter())
        chi2old = get_chi2()
        cont_maxent_ilt(me)
        iterfait = (get_iterdone())      # MAD iterdone

        if (check_diverged()):      # if diverged, restart from previous, and reduce lambsp
            read(temp)  # restart from previous
            me.set_lambsp(me.get_lambsp() / 1.2)
            lambsp(me.get_lambsp())
            print("Diverged - new lambsp : ",me.get_lambsp())
            
        else:       # if converging
            if (get_chi2() < chi2min):  # if better, (as it should), keep it
                writec(temp)
                chi2min = get_chi2()
                entropymin = get_entropy()
                st = ("iter. count %i chi2 %f"%(iterfait, chi2min))
                compteur=0
                compt_moved=0
            elif ((math.fabs(get_chi2() - chi2old)) < (0.0001 * get_chi2())):   # if not changing by more than 0.0001 == stalled
                compteur = (compteur + 1)
                print("Stalled")
            else:   # going up
                if (get_chi2()/chi2min < 4 ):
                    me.set_lambsp((me.get_lambsp() * 1.2))
                    lambsp(me.get_lambsp())
                    print("Going up - new lambsp : ",me.get_lambsp())
                else: 
                    me.set_lambsp((me.get_lambsp() / 2))
                    lambsp(me.get_lambsp())
                    print("Diverging - new lambsp : ",me.get_lambsp())
                compt_moved =compt_moved+1
                compteur = 0
        if (compteur >= 2):     #if starting to stall, try to move it 
            if (get_chi2()/chi2min < 4 ):
                me.set_lambsp((me.get_lambsp() * 1.2))
                lambsp(me.get_lambsp())
                print("Stalled - new lambsp : ",me.get_lambsp())
            else: 
                me.set_lambsp((me.get_lambsp() / 2))
                lambsp(me.get_lambsp())
                print("Diverging - new lambsp : ",me.get_lambsp())
#            me.set_miniter(me.get_miniter() / 1.2)
            compt_moved =compt_moved+1
            compteur = 0
        if (compt_moved >= 4):     #if helpless, try more miniiteration (costly !)
            compt_moved = 0
            me.set_miniter(me.get_miniter()+1)
            miniter(me.get_miniter())
            me.set_lambsp((me.get_lambsp() / 2))
            lambsp(me.get_lambsp())
            print("Stalled - new me_int_iteration : ",me.get_miniter()," - new labsp : ",me.get_lambsp())
    return st

def do_maxent_ilt(me):
    """ copy me to kernel and run maxent accordingly"""
    me.copyto_kernel()
    if (me.get_ilttype() == "Tabulated") :
        invtlap(me.get_iltsize())
    elif (me.get_ilttype() == "Regular"):
        invlap(me.get_iltsize())    # this one is a direct - not really debuged so far !
    else:
        raise "error with ilttype"

def cont_maxent_ilt(me):
    """Continue from a do_maxent_ilt run"""
    if (me.get_ilttype() == "Tabulated") :
        invtlapcont()
    elif (me.get_ilttype() == "Regular"):
        invlapcont()    # this one is a direct - not really debuged so far !
    else:
        raise "error with ilttype"

def check_diverged():    
    """ check if MaxEnt is ok (entropy and chi2 != NaN)
        true if something diverged
    """
    diverged=(1==0)
    try:
        ent = 1.0*get_entropy()
    except:
        print("diverged by entropy")
        diverged = (1==1)
        ent=1.0
    try:
        chi2 = 1.0*get_chi2()
    except:
        print("diverged by chi2")
        chi2=1.0
        diverged = (1==1)
    return diverged

#------------------------------------------------------------
def load_sq_tab(filename):
    """
    load in TAB buffer a simple 1D series;
    takes the square of the values, usefull for loading DOSY buffers
    skip # and ; comments
    
    destroy the content of the 1D buffer!
    see also : load
    """
    try:
        fin = open(filename)
    except:
        raise filename," cannot be accessed"
    # read file
    list = []
    f=fin.read()
    ls= f.split("\n")
    for v in ls:
        v = v.lstrip()
        if ( v and (v[0] != '#') and (v[0] != ';')):    # skip empty lines and comments
            list.append(float(v)**2)
    fin.close()
    # copy to 1D buffer
    dim(1)
    chsize(len(list))
    for i in range(len(list)):
            setval( i+1,list[i])
    put("tab")

#------------------------------------------------------------
def calibdosy(litdelta, bigdelta, recovery=0.0, seq_type='ste', nucleus='1H', maxgrad=50.0, maxtab=50.0, gradshape=1.0, unbalancing=0.0):
    """
    returns the DOSY calibrating factor from the parameters
    
      bigdelta float
       "Big Delta"  : diffusion delay in msec

      litdelta float
       "little Delta"  : gradient duration in msec

      seq_type enum "pgse","ste","bpp_ste","ste_2echoes","bpp_ste_2echoes","oneshot" / default ste
       the type of DOSY sequence used
        pgse : the standard hahn echoe sequence
        ste : the standard stimulated echoe sequence
        bpp_ste : ste with bipolar gradient pulses
        ste_2echoes : ste compensated for convection
        bpp_ste_2echoes : bpp_ste compensated for convection
        oneshot : the oneshot sequence from Pelta, Morris, Stchedroff, Hammond, 2002, Magn.Reson.Chem. 40, p147
        
      nucleus enum "1H","2H","13C","15N","17O","19F","31P" / default 1H
       the observed nucleus

     recovery float
       Gradient recovery delay

     maxgrad float
       Maximum Amplificator Gradient Intensity, in G/cm    / default 50.0

     maxtab float
       Maximum Tabulated Gradient Value in the tabulated file. / default 100.0
       Bruker users with gradient list in % use 100 here
       Bruker users with gradient list in G/cm (difflist) use maxgrad here
       Varian users use 32768 here

     gradshape float
       integral factor depending on the gradient shape used / default 1.0
       typical values are :
           1.0 for rectangular gradients
           0.6366 = 2/pi for sine bell gradients
           0.4839496 for 4% truncated gaussian (Bruker gauss.100 file)
       Bruker users using difflist use 1.0 here, as it is already included in difflist

    """
# MAD : modified august-sept 2007 - corrected ste; added oneshot; added PGSE

    g = (maxgrad / maxtab)*1E-4 # now in Tesla/cm
    aire = g*gradshape*litdelta
    if nucleus == "1H":
        gama = 2.675E8                           # rad.s-1.T-1
    elif nucleus == '2H':
        gama = 0.411E8                           # rad.s-1.T-1
    elif nucleus =='13C':
        gama = 0.673E8                           # rad.s-1.T-1
    elif nucleus == '15N':
        gama = -0.271E8                          # rad.s-1.T-1
    elif nucleus == '17O':
        gama = -0.363E8                          # rad.s-1.T-1
    elif nucleus =='19F':
        gama = 2.517E8                           # rad.s-1.T-1
    elif nucleus =='31P':
        gama = 1.083E8                           # rad.s-1.T-1
    else:
        raise'Unknown nucleus'

    K = ((gama * aire)**2)			  # Calcul de q^2

# equation references are in 	Jerschow,A.;Muller,N.;JMR;125;1997;Suppresion of convection artifacts
    if seq_type == 'ste':
        K = (K * (bigdelta + ((2 * litdelta)/3) + recovery))                 # cm 2 sec-1 pour Q
    elif seq_type == 'bpp_ste':
        K = (K * (bigdelta + ((2 * litdelta)/3) + ((3 * recovery)/4))) # cm 2 sec-1 pour Q
    elif seq_type == 'ste_2echoes':
        K = (K * (bigdelta + ((4 * litdelta)/3) + (2 * recovery)))     # cm 2 sec-1 pour Q
    elif seq_type == 'bpp_ste_2echoes':
        K = (K * (bigdelta + ((4 * litdelta)/3) + ((3 * recovery)/2))) # cm 2 sec-1 pour Q
    elif seq_type == 'oneshot':
        K = (K * (bigdelta + litdelta * (unbalancing * unbalancing - 2) / 6 + recovery * (unbalancing * unbalancing - 1) / 2))
    elif seq_type == 'pgse':
        K = (K * bigdelta + (2 * litdelta)/3)
    else:
        raise 'Unknown sequence'

    K = (K * 1e-8)
    return(1/K)

#------------------------------------------------------------
def damp_to_ppm(dmin,dmax,reverse=0):
    """compute pseudo ppm axis parameters

    if reverse = 1, means axis was reversed with respect to standard MaxEnt processing (as in TopSpin/Dosym

    """
    if (reverse == 0):
        fq=-1.0
        off=(math.log(dmax)/math.log(10))
        sw=(math.log(dmax)-math.log(dmin))/math.log(10)
    else:
        fq=1.0
        off=(math.log(dmin)/math.log(10))
        sw=(math.log(dmax)-math.log(dmin))/math.log(10)
    return (sw,fq,off)
#------------------------------------------------------------
def auto_damp_width():
    """
    uses the tab buffer to determine the optimum dmin and dmax for ILT processing
    """
    mn = geta_tab(2)
    mx = mn
    for i in range(3,get_si_tab()+1):
        x =  geta_tab(i)
        mn = min(mn,x)
        mx = max(mx,x)
    print(mn,mx)
    idmax = (get_dfactor()*2/(mn*mn))
    idmin = (get_dfactor()*0.5/(mx*mx))
    logidmin = (int(math.log(idmin)/math.log(10)))
    edmin = (math.exp(math.log(10)*logidmin))
    dmin = (int((int(idmin/edmin))*edmin))

# calculate dmax from idmax by taking it's upper 10th decimal
    logidmax = (int(math.log(idmax)/math.log(10)))
    edmax = (math.exp(math.log(10)*logidmax))

    dmax = (int(((int(idmax/edmax))+1)*edmax))

    return dmin, dmax

#------------------------------------------------------------
def RemoveRows(list):
    """ remove the rows of the current 2D  which index is given in list """
    print(get_si1_2d(), get_si2_2d())
    l=list
    l.sort(lambda x,y :cmp(int(x),int(y)) )   # sort the list
    for i in (l[::-1]): # then go through backward
         print("***",i)
         dim(2)
         for j in range(int(i),get_si1_2d()):
             row(j+1)
             put("row",j)
         chsize(get_si1_2d()-1,get_si2_2d())
    writec("RR.gs2")

#------------------------------------------------------------
def RemovePoints(list):
    """ remove the points of the current 1D which index is given in list """
    l=list
    l.sort(lambda x,y :cmp(int(x),int(y)) )   # sort the list
    for i in (l[::-1]): # then go through backward
        print("***",i)
        dim(1)
        for j in range(int(i),get_si1_1d()):
            setval(j,val1d(j+1))
        chsize(get_si1_1d()-1)

#------------------------------------------------------------
class maxent_state:
# Je mets tout dans cet objet, il faudra le couper en morceau
#   - maxent
#   - DOSY (ILTalgo, iltsize)
#   etc..
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
lambsp :   %i
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
        elif ( preset_value == 2):
            self.__iter = 10
            self.__ndisp = 5
        elif ( preset_value == 3):
            self.__iter = 20
            self.__ndisp = 5
        elif ( preset_value == 4):
            self.__iter = 40
            self.__ndisp = 10
        elif ( preset_value == 5):
            self.__iter = 100
            self.__ndisp = 10
        else:
            raise "wrong me_preset value"
#         stepmax  1.2
        tolerance(1.0e-2)
#         step  1.0
        self.__lambsp = 5.0
        self.__miniter = 4
        self.__lambcont = 1
        self.__algo = 1
        
