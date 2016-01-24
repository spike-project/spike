#!/usr/bin/env python 
# encoding: utf-8

"""
    Utility for handling Bruker files
"""

from __future__ import print_function

__author__ = "Marc A. Delsuc <delsuc@igbmc.fr>"
__date__ = "Oct 2009"

import struct
import re
import os
import time

zerotimeposition=0.0

from Kore import *
import unittest
import ConfigParser

#compatibility(locals())
#import Kore
#K = Kore.kore

################################################################
def find_acqu_gene(dir=".",acqulist=('acqus','acqu','ACQUS','ACQU')):
    """
    find a Bruker acqu file associated to the directory dir and return its name
    
    """
    import os.path as op

# search acqu param file
    found = 0
    l = acqulist
    for a in l:
        filename=(op.join(dir,a))
        try:    # check all versions of acqu file
            f = open(filename)
            f.close()
            found = 1
            break
        except:
            print(a+" not found")
    if (found == 0):
        raise "No acqu file of type "+l[0]+" found in "+dir
    return(filename)
################################################################
def find_acqu(dir="."):
    """
    find a Bruker acqu file associated to the directory dir and return its name
    
    """
    return( find_acqu_gene(dir,('acqus','acqu','ACQUS','ACQU')) )

################################################################
def find_acqu2(dir="."):
    """
    find a Bruker acqu2 file associated to the directory dir and return its name
    
    """
    return( find_acqu_gene(dir,('acqu2s','acqu2','ACQU2S','ACQU2')) )

################################################################
def find_acqu3(dir="."):
    """
    find a Bruker acqu3 file associated to the directory dir and return its name
    
    """
    return( find_acqu_gene(dir,('acqu3s','acqu3','ACQU3S','ACQU3')) )

################################################################
def find_proc_gene(dir=".",proclist=('procs','proc','PROCS','PROC')):
    """
    find a Bruker proc file associated to the directory dir and return its name
    
    """
    import os.path as op
    import glob

# search acqu param file
    found = 0
    l = proclist
    for a in l:
        ll = glob.glob( op.join(dir,"pdata","*",a) )
        for filename in ll:
            try:    # check all versions of acqu file
                f = open(filename)
                f.close()
                found = 1
                break
            except:
                print(filename+" not found")
        if (found == 1):
            break
    if (found == 0):
        raise "No proc file of type "+l[0]+" found in "+dir
    return(filename)

################################################################
def find_proc(dir="."):
    """
    find a Bruker proc file associated to the directory dir and return its name
    
    """
    return( find_proc_gene(dir,('procs','proc','PROCS','PROC')) )

################################################################
def find_proc2(dir="."):
    """
    find a Bruker proc file associated to the directory dir and return its name
    
    """
    return( find_proc_gene(dir,('proc2s','proc2','PROC2S','PROC2')) )

################################################################
def find_proc3(dir="."):
    """
    find a Bruker proc file associated to the directory dir and return its name
    
    """
    return( find_proc_gene(dir,('proc3s','proc3','PROC3S','PROC3')) )

################################################################
def read_param(filename="acqus"):
    """ 
    load a Bruker acqu or proc file as a dictionnary
    
    arrayed values are stored in python array
    
    comments (lines starting with $$) are stored in the special entrey [comments]
    
    M-A Delsuc jan 2006
    oct 2006 : added support for array
    """
    debug = 0
    try:
        fin = open(filename)
    except:
        raise filename," cannot be accessed"
    # read file
    dict = {}
    dict['comments']=""
    f=fin.read()
    fin.close()
    ls= f.split("\n")

#    for v in ls:
    while ls:
        v=ls.pop(0)
        v = v.lstrip()
        if debug: print("-",v,"-")
        if (re.search(r"^\$\$",v)):  # print comments
            dict['comments']=dict['comments']+"\n"+v
        else:
            m=re.match(r"##(.*)= *\(0\.\.([0-9]*)\)(.*)$",v )   # match arrays
            if (m is not None):
                if debug: print("ARRAY",v,m.group(1,2,3))
                (key,numb,line)=m.group(1,2,3)
                v=ls.pop(0)
                v = v.lstrip()
                while (not re.match(r"##",v)):    # concatenate all array lines
                    line = line+" "+v
                    v=ls.pop(0)
                    if debug: v = v.lstrip()
                ls.insert(0,v)
                array=line.split()
                if debug: print(key,numb,len(array),array)
                if ((int(numb)+1) != len(array)):   # (0..9) is 10 entries !
                    raise "size mismatch in array"
                dict[key] = array
                continue
            m=re.match(r"##(.*)= *<(.*)>",v )   #match string
            if (m is not None): 
                if debug: print("STRING",v)
                (key,val) = m.group(1,2)
                dict[key] = val
                continue
            m=re.match(r"##(.*)= *(.*)$",v )   #match value
            if (m is not None):
                if debug: print("VAL",v)
                (key,val) = m.group(1,2)
                dict[key] = val
                continue
# debug code
    if debug:
        for i in dict.keys():
            print(i+" = "+str(dict[i]))
    return dict

################################################################
def read_1D(size,filename="fid",bytorda=0):
    """
    Reads in a Bruker 1D fid
    
    size is the number of data-points in the fid
    uses struct
    """
    import struct

# read binary
    dim(1)
    chsize(int(size))

    f=open(filename,"rb")
    buf = f.read(4*size)
    ibuf = struct.unpack(">%di"%(size),buf)  # upack the string as integers > is used to keep the normal endianess

    for i in range(size):
        setval(i+1,ibuf[i])      # copy to 1D buffer
    f.close()
    
################################################################
def read_2D(sizeF1,sizeF2,filename="ser"):
    """
    Reads in a Bruker 2D fid

    sizeF1 is the number of fid
    sizeF2 is the number of data-points in the fid
    uses struct
    """
    import struct
# read binary
    dim(2)
    chsize(int(sizeF1),int(sizeF2))
    f=open(filename,"rb")
    for i1 in range(sizeF1):
        for i2 in range(0,sizeF2,64):   # read by 64 steps
            buf = f.read(256)
            ibuf = struct.unpack(">64i",buf)   # > is used to keep the normal endianess
            bufsz = min(64, sizeF2-i2)
            for i3 in range(bufsz):
                setval(i1+1,i2+i3+1,ibuf[i3])      # copy to 2D buffer
    f.close()
################################################################
def read_3D(sizeF1,sizeF2,sizeF3,filename="ser"):
    """
    Reads in a Bruker 3D fid

    sizeF1 x sizeF2 is the number of fid
    sizeF3 is the number of data-points in the fid
    uses struct
    """
    import struct
# read binary
    dim(3)
    chsize(int(sizeF1),int(sizeF2),int(sizeF3))
    f=open(filename,"rb")
    for i1 in range(sizeF1):
        print(i1, end=' ')
        for i2 in range(sizeF2):
            for i3 in range(0,sizeF3,64): # read by 64 steps
                buf = f.read(256)
                ibuf = struct.unpack(">64i",buf)   # > is used to keep the normal endianess
                bufsz = min(64, sizeF3-i3)
                for i4 in range(bufsz):
                    setval(i1+1,i2+1,i3+i4+1,ibuf[i4])      # copy to 3D buffer
    f.close()
################################################################
def read_1D_java(size,filename="fid"):
    """
    Reads in a Bruker 1D fid

    size is the number of data-points in the fid
    uses java low level file access
    """
    from java.io import RandomAccessFile
#    from nmrtec.util.IO import SwapByte
    import jarray

# read binary
    dim(1)
    chsize(int(size))
    f=RandomAccessFile(filename,"r")
    for i in range(size):
        x = f.readInt()
#        x = SwapByte.swap(x)       # uncoment if swapbyte is required
#        print i,x
        setval(i+1,float(x))      # copy to 1D buffer

################################################################
def read_2D_java(sizeF1,sizeF2,filename="ser"):
    """
    Reads in a Bruker 2D fid

    sizeF1 is the number of fid
    sizeF2 is the number of data-points in the fid
    uses java low level file access
    """
    from java.io import RandomAccessFile
#    from nmrtec.util.IO import SwapByte
    import jarray

# read binary
    dim(2)
    chsize(int(sizeF1),int(sizeF2))
    f=RandomAccessFile(filename,"r")
    for i1 in range(sizeF1):
        for i2 in range(sizeF2):
            x = f.readInt()
#           x = SwapByte.swap(x)       # uncoment if swapbyte is required
            setval(i1+1,i2+1,float(x))      # copy to 2D buffer
        if ((sizeF2 % 256) != 0):  # pb with Bruker files, where size is always a 256 multiple
            for i2 in range(((sizeF2 / 256)+1)*256 - sizeF2):
                x = f.readInt()
################################################################
def read_3D_java(sizeF1,sizeF2,sizeF3,filename="ser"):
    """
    Reads in a Bruker 3D fid
    
    sizeF1 x sizeF2 is the number of fid
    sizeF3 is the number of data-points in the fid
    uses java low level file access
    """
    from java.io import RandomAccessFile
#    from nmrtec.util.IO import SwapByte
    import jarray

# read binary
    dim(3)
    chsize(int(sizeF1),int(sizeF2),int(sizeF3))
    f=RandomAccessFile(filename,"r")
    for i1 in range(sizeF1):
        print(i1, end=' ')
        for i2 in range(sizeF2):
            for i3 in range(sizeF3):
                x = f.readInt()
#               x = SwapByte.swap(x)       # uncoment if swapbyte is required
                setval(i1+1,i2+1,i3+1,float(x))      # copy to 3D buffer
            if ((sizeF3 % 256) != 0):  # pb with Bruker files, where size is always a 256 multiple
                for i3 in range(((sizeF3 / 256)+1)*256 - sizeF3):
                    x = f.readInt()

################################################################
def zerotime(acqu):
    """get digital filter parameters, if any
    
  The zerotime function computes the correction for the Bruker figital filter.
  the phase correction to apply is computed given the 3 parameters :
  DSPFIRM DSPFVS DECIM
  as found in the acqus parameter file in XwinNMR

  correction is then -360*zerotime in firstorder correction

  dspfvs is not used so far
  oct 2006
  """
    tabdelay=   [
       [ 179,201,533,709,1097,1449,2225,2929,4481,5889,8993,11809,18017,23649,36065,47329,72161, 94689,144353,189409,288737 ],
       [ 184,219,384,602, 852,1668,2312,3368,4656,6768,9344,13568,18560,27392,36992,50040,73856,110336,147584,220928,295040 ],
       [ 184,219,384,602, 852,1668,2292,3369,4616,6768,9264,13568,18560,27392,36992,50040,73856,110336,147584,220928,295040 ],
       [  11, 17, 23, 35,  47,  71,  95, 143, 191, 287, 383,  575,   -1,   -1,   -1,   -1,   -1,    -1,    -1,    -1,    -1 ],
       [  60, 90,118,179, 244, 360, 492, 724, 980,1444,1958, 2886, 3912, 5768, 7820,11532,   -1,    -1,    -1,    -1,    -1 ],
       [  -1, -1, 58,152, 202, 318, 418, 642, 842,1290,1690, 2586, 3386,   -1,   -1,   -1,   -1,    -1,    -1,    -1,    -1 ]
       ]
    decim_offset = [
       2, 3,  4,  6,   8,  12,  16,  24,  32,  48,  64,   96,   128,  192,  256,  384, 512,   768,  1024,  1536,  2048 ]

    try:
        decim = int(acqu['$DECIM'])
#        print "DECIM=",decim
    except:
        decim=1
    if (decim!=1):
        dspfvs = int(acqu['$DSPFVS'])
#        print "DSPFVS=",dspfvs
        dspfirm = int(acqu['$DSPFIRM'])
#        print "DSPFIRM=",dspfirm
# first special cases
        if (dspfvs >= 20  and  dspfvs <= 23):
            zerotimeposition= float(acqu['$GRPDLY'])     # new (aka 2005) DRU
            return(zerotimeposition)

        if (dspfvs == 15 and decim == 3 and float(acqu['$SW_h']) >= 104000):
            z = ( -110.0) / (decim)
            zerotimeposition = z/2
            return(zerotimeposition)

        if (dspfvs == 0 or decim == 1):
            return(0.0)

        try:
            j = decim_offset.index(decim)
        except:
            raise "*** wrong value for DECIM"
        try:
            d = tabdelay[(dspfvs) - 10][j];
#            print "d=",d
        except:
             raise "*** wrong value for DSPFVS " + dspfirm
        if (d == -1):
             raise "*** wrong DECIM/DSPFVS parameter combination"
        z = (float(d)/float(decim))
        zerotimeposition = z/2
        return(zerotimeposition)

    else:
        zerotimeposition = 0.0
        return(zerotimeposition)

    return(zerotimeposition)
    
################################################################
def offset(acqu,proc):
    """computes the offset """
    try:
        ppmPointUn = float(proc['$OFFSET'])
        ppmWidth = float(acqu['$SW_h']) / float(acqu['$SFO1'])
        calibrationOffset = float(ppmPointUn - ppmWidth)*  float(acqu['$SFO1'])
    except:
        calibrationOffset=0.0
    return(calibrationOffset)


################################################################
def Import_1D(filename="fid",outfile=""):
    """
    Imports a 1D Bruker fid
    
    """
    import os.path as op
    import Generic
    global zerotimeposition

    if (not op.exists(filename)):
        raise filename+" : file not found"
    dire=op.dirname(filename)
    acqu = read_param(find_acqu(dire))
    size= int(acqu['$TD'])  # get size
    print("importing 1D FID, size =",size)
    read_1D(size,filename,bytorda=int(acqu['$BYTORDA']))

# then set parameters
    specw(float(acqu['$SW_h']))
    freq(float(acqu['$SFO1']),float(acqu['$SFO1']))

    proc = read_param(find_proc(dire))
    com_offset( offset(acqu, proc) )
    
    zerotimeposition = zerotime(acqu)
    if (outfile != ""):
        dim(1)
        writec(outfile)
#        join(outfile)
#        putheader("zerotimeposition", repr(zerotimeposition))
#        disjoin()
        param={}
        param['axisf1_zerotimeposition'] = zerotimeposition
        Generic.dict_dump(param,outfile+'.gtb')
    return (get_si1_1d(),)
    
################################################################
def Import_2D(filename="ser",outfile=""):
    """
    Imports a 2D Bruker ser
    
    """
    import os.path as op
    import Generic

    if (not op.exists(filename)):
        raise filename+" : file not found"
    dir=op.dirname(filename)
    acqu = read_param(find_acqu(dir))
    acqu2 = read_param(find_acqu2(dir))
    sizeF1= int(acqu2['$TD'])  # get size
    sizeF2= int(acqu['$TD'])  # get size
    print("importing 2D FID, size =",sizeF1,"x",sizeF2)
    read_2D(sizeF1,sizeF2,filename)
# then set parameters
    freq(float(acqu['$SFO1']), float(acqu2['$SFO1']), float(acqu['$SFO1']))
    specw(float(acqu2['$SFO1'])*float(acqu2['$SW']), float(acqu['$SW_h']))

    proc = read_param(find_proc(dir))
    proc2 = read_param(find_proc2(dir))
    com_offset( offset(acqu2, proc2), offset(acqu, proc))
    zerotimeposition = zerotime(acqu)
    if (outfile != ""):
        dim(2)
        writec(outfile)
#        K.join(outfile)
#        K.putheader("zerotimeposition", repr(zerotimeposition))
#        K.disjoin()
        param={}
        param['axisf2_zerotimeposition'] = zerotimeposition
        Generic.dict_dump(param,outfile+'.gtb')
    return (get_si1_2d(),get_si2_2d())
    

################################################################
def Import_3D(filename="ser",outfile=""):
    """
    Imports a 3D Bruker ser
    
    """
    import os.path as op
    import NPK.Generic

    if (not op.exists(filename)):
        raise filename+" : file not found"
    dir=op.dirname(filename)
    acqu = read_param(find_acqu(dir))
    acqu2 = read_param(find_acqu2(dir))
    acqu3 = read_param(find_acqu3(dir))
    proc = read_param(find_proc(dir))
    proc2 = read_param(find_proc2(dir))
    proc3 = read_param(find_proc3(dir))
    AQORDER = int(proc['$AQORDER']) # t'was wrongly reading in proc3 before jan 2009 - MAD
    if AQORDER != 0:
        (acqu3,acqu2) = (acqu2,acqu3) # exchange acqu2 and acqu3
    sizeF1= int(acqu3['$TD'])  # get size
    sizeF2= int(acqu2['$TD'])  # get size
    sizeF3= int(acqu['$TD'])  # get size
    print("importing 3D FID, size =",sizeF1,"x",sizeF2,"x",sizeF3)
    read_3D(sizeF1,sizeF2,sizeF3,filename)

# then set parameters
    freq(float(acqu['$SFO1']), float(acqu3['$SFO1']), float(acqu2['$SFO1']),float(acqu['$SFO1']))
    specw(float(acqu3['$SFO1'])*float(acqu3['$SW']),float(acqu2['$SFO1'])*float(acqu2['$SW']), float(acqu['$SW_h']))
    com_offset( offset(acqu3, proc3), offset(acqu2, proc2), offset(acqu, proc))

    zerotimeposition = zerotime(acqu)
    if (outfile != ""):
        dim(3)
        writec(outfile)
#        K.join(outfile)
#        K.putheader("zerotimeposition", repr(zerotimeposition))
#        K.disjoin()
        param={}
        param['axisf3_zerotimeposition'] = zerotimeposition
        NPK.Generic.dict_dump(param,outfile+'.gtb')
    return (get_si1_3d(),get_si2_3d(),get_si3_3d())

################################################################
def Import(audit,filename=".",outfile="-"):
    """
    Imports a file or a set of files located in a given directory
    
    """
    import os.path as op
    import glob
    import Generic
    if op.isfile(filename):
        if op.basename(filename) in ("fid","FID"):
            if outfile == "-":
                outfile=op.join(op.dirname(filename),"data.gf1")
            size = Import_1D(filename,outfile)
            Generic.audittrail( audit, "text", "imported 1D file", filename+"  converted to   ", outfile, "Size", size )

        elif op.basename(filename) in ("ser","SER"):
            try:
                t=find_acqu3(op.dirname(filename))
                # then we are in 3D
                state = "3D"
            except:
                state = "2D"
            if (state == "3D"):
                if outfile == "-":
                    outfile=op.join(op.dirname(filename),"data.gf3")
                size = Import_3D(filename,outfile)
                Generic.audittrail( audit, "text", "imported 3D file", filename+"  converted to   ", outfile, "Size", size )
            else:
                if outfile == "-":
                    outfile=op.join(op.dirname(filename),"data.gf2")
                size = Import_2D(filename,outfile)
                Generic.audittrail( audit, "text", "imported 2D file", filename+"  converted to   ", outfile, "Size", size )

    elif op.isdir(filename):
        print("**DIR**"+filename)
        for i in glob.glob(op.join(filename,"*")):
            if (i != "."):
                Import(audit,i, "-")


################################################################
def calibdosy(file="acqus"):
    """
    get the parameters from the acqus file and compute the calibration (dfactor)

    returns : (BigDelta, litteldelta, recovery, sequence, nucleus)

    assumes that you are using the standard Bruker set-up, and have started the dosy acquisition ith the dosy macro.

    from calibdosy.g in Gifa 5.2006
    """
    import re

    param=read_param(file)  # load values
    nuc1 = param["$NUC1"]
    p30 =  float(param["$P"][30])*1e-6
    p1 =  float(param["$P"][1])*1e-6
    p19 =  float(param["$P"][19])*1e-6
    d16 =  float(param["$D"][16])
    d20 =  float(param["$D"][20])
    pulprog = param["$PULPROG"]
# STEBP_2echos Bruker avance sequences
    if (re.search('dstebp',pulprog)):
        sequence = 'bpp_ste_2echoes'	
    	if nuc1  in ('1H','15N','13C','31P','19F','17O'):
            nucleus = nuc1
        delta = (2*p30)
        bdelta = (d20-(10*p1)-(8*p30)-(8*d16)-(8*d17)-(2*p19))
        tau = d16
# STE_2echos Bruker avance sequences
    elif re.search('dstegp',pulprog):
        sequence = 'ste_2echoes'
    	if nuc1  in ('1H','15N','13C','31P','19F','17O'):
            nucleus = nuc1
       	delta = p30
        bdelta = (2*(d20-(2*p1)-(p30)-(2*d16)-(p19)))
        tau = d16
# BPP_LED NMRtec and Bruker Avance sequences
    elif re.search('stegpbp|ledbp',pulprog):
        sequence = 'bpp_ste'
    	if nuc1  in ('1H','15N','13C','31P','19F','17O'):
            nucleus = nuc1
        delta = 2*p30
        bdelta = d20-(4*p1)-(2*p30)-(3*d16)-(p19)
        tau = 2*d16
# LEDgp/STEgp Bruker Avance sequence
    elif re.search('stegp|led',pulprog):
    	sequence = 'ste'
    	if nuc1  in ('1H','15N','13C','31P','19F','17O'):
            nucleus = nuc1
        delta = p30
        bdelta = d20-(2*p1)-(p30)-(2*d16)-(p19)
       	tau = d16
    else:
        print("Unsupported pulse program.")
    return (bdelta,delta,tau,sequence, nucleus)
