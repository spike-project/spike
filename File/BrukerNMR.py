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

################################################################
def find_acqu_gene(dir=".",acqulist=('acqus','acqu','ACQUS','ACQU')):
    """
    find a Bruker acqu file associated to the directory dir and return its name
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
    print(bytorda)
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
                npkbuf[i1,i2:i2+bufsz] = data[:]
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
        print(i1, end=' ')
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
    print("importing 1D FID, size =",size)
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
    print("importing 2D FID, size =",sizeF1,"x",sizeF2)
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
    d.axis1.offset = offset(acqu, proc)
    d.axis2.offset = offset(acqu2, proc2)

    d.axis2.zerotime = zerotime(acqu)
    if outfile is not None:
        raise Exception("Not implemented yet")

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
    print("importing 3D FID, size =",sizeF1,"x",sizeF2,"x",sizeF3)
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
    return (get_si1_3d(),get_si2_3d(),get_si3_3d())

################################################################
class Exporter(object):
    '''
    the following code is for exprt to Bruker file
    only the binary data files are created (ser, 1rr, 2rr, etc..)
    they should inserted into full-fledge existing data directory
    
    original code from Lionel Chiron
    
    NOT DEBUGGED - In development
    '''
    def __init__(self, direc="."):
        self.dir = direc
        self.param_acq = read_param( find_acqu(direc) )
        self.param_proc = read_param( find_proc(direc) )
        
    def write_file(self, data, filename):
        '''
        data written as integers. 
        '''
        with open(filename, 'wb') as f:
            if self.param_acq['$BYTORDA'] == '0': 
                f.write(data.astype('<i4').tostring()) # little endian
            else:
                f.write(data.astype('>i4').tostring()) # big endian

    def save_fid_1d(self, data):
        """
        save data as a 1d fid in self.dir
        """
        fname = op.join(self.dir, 'fid')
        self.write_file(data, fname)

    def save_spec_1d(self, data, procno=1):
        """
        save data as a 1d spectrum in self.dir
        """
        fname = op.join(self.dir, 'pdata', proco, '1r')
        self.write_file(real(data), fname)
        if isinstance(data[0] , complex):
            fname = op.join(self.dir, 'pdata', proco, '1i')
            self.write_file(data, fname.imag)

#### not debugged further down !
    def reorder_bck_subm(self,data):
        """
        Reorder flat matrix back to sbmx Bruker data.
        self.sub_per_dim : [nb of submatrices in t1, nb of submatrices in t2]
        self.nsubs : total number of submatrices
        self.param_proc['$SI'] : shape of the 2D data. 
        self.param_acq['$XDIM'] : size submatrix
        """
        print("reorder matrix back")
        self.prepare_mat()
        interm = data.reshape(self.dim_mat)
        mat = []
        for sub_num, sub_idx in enumerate(np.ndindex(tuple(self.sub_per_dim))):
            zipshape = zip(sub_idx, self.dim_sub_mat)
            sub_slices = [slice(i * j, (i + 1) * j) for i, j in zipshape ]
            mat.append(list(np.ravel(interm[sub_slices])))
        data = np.array(mat).reshape(self.dim_mat)
        return data

    def reorder_subm(self, data):
        """
        Reorder sbmx binary Bruker data to flat matrix.
        self.sub_per_dim : [nb of submatrices in t1, nb of submatrices in t2]
        self.nsubs : total number of submatrices
        self.param_proc['$SI'] : shape of the 2D data. 
        self.param_acq['$XDIM'] : size submatrix
        """
        print("reorder matrix")
        self.prepare_mat()
        #longmat = int(self.param_proc['$SI']), int(self.param_proc2['$XDIM'])*self.sub_per_dim[1]
        print("data.shape ",data.shape)
        print("longmat ",self.long_mat)
        interm = data.reshape(self.long_mat)      
        mat = []
        for sub_num, sub_idx in enumerate(np.ndindex(tuple(self.sub_per_dim))):
            zipshape = zip(sub_idx, self.dim_sub_mat)
            sub_slices = [slice(i * j, (i + 1) * j) for i, j in zipshape ]
            slt2 = slice(sub_num*self.dim_sub_mat[0],(sub_num+1)*self.dim_sub_mat[0]) # dimension t1
            slt1 = slice(0,self.dim_sub_mat[1])# dimension t2
            print("self.rdata[sub_slices].shape ",self.rdata[sub_slices].shape)
            print("interm[slt1, slt2].shape ",interm[slt1, slt2].shape)
            print("self.rdata[sub_slices].shape",self.rdata[sub_slices].shape)
            self.rdata[sub_slices] = interm[slt2, slt1]
        data = self.rdata
        return data


    def save_spec_2d(self, big = False):
        """ 
        Write Bruker binary data to file
        big or little endianess.
        """
        print("save_denoised_2d")
        #filename = op.join(self.addrproc, procfile)
        #print "remove ", filename
        if self.param_acq2['$FnMODE'] == '0':
            list_proc= ['2rr','2ri']
            for name_proc in list_proc:
                if name_proc == '2rr':
                    if os.path.exists(self.addr_data_2rr):
                        os.remove(self.addr_data_2rr)
                    self.write_file(self.data_2d_denoised_2rr, self.addr_data_2rr)
                else:
                    if os.path.exists(self.addr_data_2ri):
                        os.remove(self.addr_data_2ri)
                    self.write_file(self.data_2d_denoised_2ri, self.addr_data_2ri)
                    
        if self.param_acq2['$FnMODE'] == '3' :
            list_proc= ['2rr','2ir']
            for name_proc in list_proc:
                if name_proc == '2rr':
                    if os.path.exists(self.addr_data_2rr):
                        os.remove(self.addr_data_2rr)
                    self.write_file(self.data_2d_denoised_2rr, self.addr_data_2rr)
                else:
                    if os.path.exists(self.addr_data_2ir):
                        os.remove(self.addr_data_2ir)
                    self.write_file(self.data_2d_denoised_2ir, self.addr_data_2ir)
        
        
        if self.param_acq2['$FnMODE'] == '1':
            list_proc= ['2rr','2ii']
            for name_proc in list_proc:
                if name_proc == '2rr':
                    if os.path.exists(self.addr_data_2rr):
                        os.remove(self.addr_data_2rr)
                    self.write_file(self.data_2d_denoised_2rr, self.addr_data_2rr)
                else:
                    if os.path.exists(self.addr_data_2ii):
                        os.remove(self.addr_data_2ii)
                    self.write_file(self.data_2d_denoised_2ii, self.addr_data_2ii)

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
