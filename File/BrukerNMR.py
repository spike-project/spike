#!/usr/bin/env python 
# encoding: utf-8

"""
    Utility to Handle NMR Bruker files

partitly based on NPK v1 code
"""

from __future__ import print_function

__author__ = "Marc Andre' Delsuc"
__date__ = "november 2014"

import os
import os.path as op
import glob
import re
import struct
import time
import unittest
import numpy as np

from ..NPKData import NPKData

verbose = False   # change this for vebose importers
################################################################
def find_acqu_proc_gene(dir, acqulist):
    """
    find a Bruker acqu or proc file associated to the directory dir and return its name
    """
# search acqu param file
    found = 0
    for a in acqulist:
        filename=(op.join(dir,a))
        if op.exists(filename):
            break
    else:
        raise Exception("no parameter file found in "+dir)
    return(filename)
################################################################
def find_acqu(dir="."):
    """
    find a Bruker acqu file associated to the directory dir and return its name
    
    """
    return( find_acqu_proc_gene(dir,('acqus','acqu','ACQUS','ACQU')) )

################################################################
def find_acqu2(dir="."):
    """
    find a Bruker acqu2 file associated to the directory dir and return its name
    
    """
    return( find_acqu_proc_gene(dir,('acqu2s','acqu2','ACQU2S','ACQU2')) )

################################################################
def find_acqu3(dir="."):
    """
    find a Bruker acqu3 file associated to the directory dir and return its name
    
    """
    return( find_acqu_proc_gene(dir,('acqu3s','acqu3','ACQU3S','ACQU3')) )

################################################################
def find_proc_down(dir, proclist):
    """
    find a Bruker proc file associated to the directory dir and return its name
    
    search in pdada/PROCNO and returns the first one
    """
# search acqu param file
    for fname in glob.glob( op.join(dir,"pdata","*","*") ):
        if op.basename(fname) in proclist:  # found
            break
    else:
        raise Exception("No proc file found in "+dir)
    return(fname)

################################################################
def find_proc(dir=".", down=True):
    """
    find a Bruker proc file associated to the directory dir and return its name
    if  down is True - search in dir/pdata/* other searches here
    """
    if down:
        return( find_proc_down(dir,('procs','proc','PROCS','PROC')) )
    else:
        return( find_acqu_proc_gene(dir,('procs','proc','PROCS','PROC')) )
################################################################
def find_proc2(dir=".", down=True):
    """
    find a Bruker proc file associated to the directory dir and return its name
    
    """
    if down:
        return( find_proc_down(dir,('proc2s','proc2','PROC2S','PROC2')) )
    else:
        return( find_acqu_proc_gene(dir,('proc2s','proc2','PROC2S','PROC2')) )

################################################################
def find_proc3(dir=".", down=True):
    """
    find a Bruker proc file associated to the directory dir and return its name
    
    """
    if down:
        return( find_proc_down(dir,('proc3s','proc3','PROC3S','PROC3')) )
    else:
        return( find_acqu_proc_gene(dir,('proc3s','proc3','PROC3S','PROC3')) )

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
    with open(filename) as fin:
        # read file
        dict = {}
        dict['comments']=""
        f=fin.read()
        fin.close()
        ls= f.split("\n")

    #    for v in ls:
        while ls:
            v=ls.pop(0)
            v = v.strip()
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
def read_1D(size, filename="fid", bytorda=1):
    """
    Reads in a Bruker 1D fid as a numpy float array
    
    size is the number of data-points in the fid
    uses struct
    does not check endianess
    """
# read binary
    if bytorda == 0:
        fmt = "<%di"
    else:
        fmt = ">%di"   # > is used to keep the normal endianess
    with open(filename,"rb") as F:
        buf = F.read(4*size)
        ibuf = struct.unpack(fmt%(size), buf)  # upack the string as integers
    npkbuf = np.empty(size, dtype=np.float64)
    npkbuf[:] = ibuf[:]
    return npkbuf
    
################################################################
def read_2D(sizeF1, sizeF2, filename="ser", bytorda=1):
    """
    Reads in a Bruker 2D fid as a numpy float array

    sizeF1 is the number of fid
    sizeF2 is the number of data-points in the fid
    """
    npkbuf = np.empty((sizeF1, sizeF2), dtype=np.float64)
# read binary
    if bytorda == 0:
        fmt = "<64i"
    else:
        fmt = ">64i"   # > is used to keep the normal endianess
    with open(filename,"rb") as F:
        for i1 in range(sizeF1):
            for i2 in range(0, sizeF2, 64):   # read by 64 steps - MAndatory !
                buf = F.read(256)
                data = struct.unpack(fmt, buf)
                bufsz = min(64, sizeF2-i2)
                npkbuf[i1,i2:i2+bufsz] = data[:bufsz]
                # for i3 in range(bufsz):
                #     setval(i1+1,i2+i3+1,ibuf[i3])      # copy to 2D buffer
    return npkbuf
################################################################
def read_3D(sizeF1, sizeF2, sizeF3, filename="ser", bytorda=1):
    """
    Reads in a Bruker 3D fid

    sizeF1 x sizeF2 is the number of fid
    sizeF3 is the number of data-points in the fid

    should be on file.
    """
    raise Exception("Not performed yet")
# read binary
    if bytorda == 0:
        fmt = "<64i"
    else:
        fmt = ">64i"   # > is used to keep the normal endianess
    f=open(filename,"rb")
    for i1 in range(sizeF1):
        for i2 in range(sizeF2):
            for i3 in range(0,sizeF3,64): # read by 64 steps
                buf = f.read(256)
                ibuf = struct.unpack(fmt,buf)
                bufsz = min(64, sizeF3-i3)
                for i4 in range(bufsz):
                    setval(i1+1,i2+1,i3+i4+1,ibuf[i4])      # copy to 3D buffer
    f.close()

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
        decim = int(float(acqu['$DECIM']))
#        print("DECIM=",decim)
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
def Import_1D(filename="fid", outfile=None):
    """
    Imports a 1D Bruker fid as a NPKData
    
    """
    if (not op.exists(filename)):
        raise Exception(filename+" : file not found")
    dire=op.dirname(filename)
    acqu = read_param(find_acqu(dire))
    size= int(acqu['$TD'])  # get size
    if verbose: print("importing 1D FID, size =",size)
    data = read_1D(size, filename, bytorda=int(acqu['$BYTORDA']))
    d = NPKData(buffer=data)
# then set parameters
    d.axis1.specwidth = float(acqu['$SW_h'])
    d.axis1.frequency = float(acqu['$SFO1'])
    d.frequency = d.axis1.frequency
    if acqu['$AQ_mod'] == '0':  # QF
        d.axis1.itype = 0
    else:                       # complex mode
        d.axis1.itype = 1
    proc = read_param(find_proc(dire))
    d.axis1.offset = offset(acqu, proc)

    d.axis1.zerotime = zerotime(acqu)
    if outfile is not None:
        raise Exception("Not implemented yet")
    pardic = {"acqu": acqu, "proc": proc} # create ad-hoc parameters
    d.params = pardic   # add the parameters to the data-set
    return d
################################################################
def Import_1D_proc(filename="1r"):
    """
    Imports a 1D Bruker 1r processed file as a NPKData
    if 1i exists imports the complex spectrum

    """
    if (not op.exists(filename)):
        raise Exception(filename+" : file not found")
    dire=op.dirname(filename)
    proc = read_param(find_proc(dire, down=False))
    diracq = op.dirname(op.dirname(op.dirname(filename)))
    acqu = read_param(find_acqu(diracq))
    if verbose:     print("importing 1D spectrum")
    data = np.fromfile(filename, 'i4').astype(float)
    if op.exists(op.join(dire,'1i')):  # reads imaginary part
        data = data + 1j*np.fromfile(filename, 'i4')
    d = NPKData(buffer=data)
# then set parameters
    d.axis1.specwidth = float(acqu['$SW_h'])
    d.axis1.frequency = float(acqu['$SFO1'])
    d.frequency = d.axis1.frequency
    d.axis1.offset = offset(acqu, proc)
    d.axis1.zerotime = zerotime(acqu)
    pardic = {"acqu": acqu, "proc": proc} # create ad-hoc parameters
    d.params = pardic   # add the parameters to the data-set
    return d
    
################################################################
def FnMODE(acqu, proc):
    """
    complex type along F1 for a 2D or 3D
    search FnMODE in acqu
    and if absent, search MC2 in proc
    None    0
    QF      1
    QSEQ    2
    TPPI    3
    States  4
    States-TPPI 5
    Echo-AntiEcho   6
    
    returns either 0 (real) or 1 (complex)
    """
    try:    # find acquisition mode
        mode = acqu['$FnMODE']  # new school
    except KeyError:
        try:
            mode = proc['$MC2'] # old school
        except KeyError:
            mode = '0'          # desperate
    if mode in ('4', '5'):
        r = 1
    else:
        r = 0
    return r
    
def Import_2D(filename="ser", outfile=None):
    """
    Imports a 2D Bruker ser
    
    """
    if (not op.exists(filename)):
        raise Exception(filename+" : file not found")
    dire=op.dirname(filename)
    acqu = read_param(find_acqu(dire))
    acqu2 = read_param(find_acqu2(dire))
    proc = read_param(find_proc(dire))
    proc2 = read_param(find_proc2(dire))
    sizeF1= int(acqu2['$TD'])  # get size
    sizeF2= int(acqu['$TD'])  # get size

    data = np.fromfile(filename, 'i4').astype(float)
    if op.exists(op.join(dire,'1i')):  # reads imaginary part
        data = data + 1j*np.fromfile(filename, 'i4')
    d = NPKData(buffer=data)
    
    data = read_2D(sizeF1, sizeF2, filename,  bytorda=int(acqu['$BYTORDA']))
    d = NPKData(buffer=data)
# then set parameters
    d.axis1.specwidth = float(acqu2['$SW_h'])
    d.axis1.frequency = float(acqu2['$SFO1'])
    d.axis2.specwidth = float(acqu['$SW_h'])
    d.axis2.frequency = float(acqu['$SFO1'])
    d.frequency = d.axis2.frequency
    if acqu['$AQ_mod'] == '0':  # QF
        d.axis2.itype = 0
    else:                       # complex mode
        d.axis2.itype = 1
    d.axis1.itype = FnMODE(acqu2, proc2)
    d.axis1.offset = offset(acqu2, proc2)
    d.axis2.offset = offset(acqu, proc)

    d.axis2.zerotime = zerotime(acqu)
    if outfile is not None:
        raise Exception("Not implemented yet")

    pardic = {"acqu": acqu, \
        "acqu2": acqu2, \
        "proc": proc, \
        "proc2": proc2} # create ad-hoc parameters
    d.params = pardic   # add the parameters to the data-set
    return d

def Import_2D_proc(filename="2rr", outfile=None):
    """
    Imports a 2D Bruker ser
    
    """
    from BrukerSMX import BrukerSMXHandler
    if (not op.exists(filename)):
        raise Exception(filename+" : file not found")
    if verbose:     print("importing 2D spectrum")
    SMX = BrukerSMXHandler(op.dirname(filename))
    SMX.read_2D()
    datar = SMX.data_2d_2rr.astype('float')     # loads 2rr anyhow

    if SMX.data_2d_2ir is not None:             # complex in F2
        datar = datar +  1j*SMX.data_2d_2ir.astype('float')

    if SMX.data_2d_2ri is not None:             # complex in F1 stored independtly
        datai = SMX.data_2d_2ri.astype('float')
        if SMX.data_2d_2ii is not None:
            datai =  datai + 1j*SMX.data_2d_2ii.astype('float')
        data = np.concatenate((datar, datai)) # then concat
    else:
        data = datar

    d = NPKData(buffer=data)
    if SMX.data_2d_2ri is not None:     # swap if was concatenated
        d.swap(axis='F1')
# then set parameters
    d.axis1.specwidth = float(SMX.acqu2['$SW_h'])
    d.axis1.frequency = float(SMX.acqu2['$SFO1'])
    d.axis2.specwidth = float(SMX.acqu['$SW_h'])
    d.axis2.frequency = float(SMX.acqu['$SFO1'])
    d.frequency = d.axis2.frequency
    d.axis1.offset = offset(SMX.acqu2, SMX.proc2)
    d.axis2.offset = offset(SMX.acqu, SMX.proc)

    d.axis2.zerotime = zerotime(SMX.acqu)
    if outfile is not None:
        raise Exception("Not implemented yet")

    pardic = {"acqu": SMX.acqu, \
        "acqu2": SMX.acqu2, \
        "proc": SMX.proc, \
        "proc2": SMX.proc2} # create ad-hoc parameters
    d.params = pardic   # add the parameters to the data-set
    return d

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
    if verbose:     print("importing 3D FID, size =",sizeF1,"x",sizeF2,"x",sizeF3)
    read_3D(sizeF1, sizeF2, sizeF3, filename,  bytorda=int(acqu['$BYTORDA']))

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

    pardic = {"acqu": acqu, \
        "acqu2": acqu2, \
        "acqu3": acqu3, \
        "proc": proc, \
        "proc2": proc2, \
        "proc3": proc3} # create ad-hoc parameters
    d.params = pardic   # add the parameters to the data-set
    return d

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
    d17 =  float(param["$D"][17])
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

#----------------------------------------------
class Bruker_Tests(unittest.TestCase):
    """ A FAIRE"""
    def setUp(self):
        self.verbose = 1    # verbose > 0 switches messages on
    def announce(self):
        if self.verbose >0:
            print("\n========",self.shortDescription(),'===============')
    def test_import(self):
        from ..Tests import filename
        name = filename("Lasalocid-Tocsy/dataset/ser")
        d = Import_2D(name)
        self.assertEqual(d.axis1.itype, 1)
        self.assertEqual(d.axis2.itype, 1)
        self.assertAlmostEqual(d[234,567], 8729.0)
        self.assertAlmostEqual(d.axis2.frequency, 400.131880611)
        self.assertAlmostEqual(d.axis1.specwidth, 5201.560468140)
        self.assertAlmostEqual(d.axis2.zerotime, 70.1875)
        d.apod_sin(axis=2,maxi=0.5).revf().fft(axis=2).apod_sin(axis=1,maxi=0.5).fft(axis=1).modulus()
        self.assertAlmostEqual(d[65,43], 4469.4687352772253)
    #-------------------------------------------------

if __name__ == '__main__':
    unittest.main()
