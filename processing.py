# encoding: utf-8
"""
Processing.py

This program realises the processing of an FTICR data

Created by Marc-Andre on 2011-09-23.
Copyright (c) 2011 IGBMC. All rights reserved.
"""
import sys, os, time
import unittest
import numpy as np
from numpy import fft as npfft
import tables
from scipy.signal import decimate, lfilter, cheby1, medfilt, medfilt2d
from NPKConfigParser import NPKConfigParser
from FTICR import *
from File.Apex import Import_2D as Import_2D_Apex
from File.Solarix import Import_2D as Import_2D_Solarix
Import_2D = {'Apex': Import_2D_Apex,'Solarix': Import_2D_Solarix }
from NPKData import copyaxes
from File.HDF5File import HDF5File, determine_chunkshape
import util.progressbar as pg
import util.mpiutil as mpiutil
from util.signal_tools import findnoiselevel_offset  as findnoiselevel
from NPKData import as_cpx
import itertools
import multiprocessing as mp
from util.simple_logger2 import TeeLogger
import pickle


debug = 1   # 0 means no debugging
interfproc = False


'''
Processing for performing urQRd  on 2D FTICR datasets.
previous version was named processing2-urqrd-superresol
under Linux or MacOsX : mpirun -n nbproc python processing.py (configfile.mscf)
under Windows : mpiexec -n nbproc python processing.py (configfile.mscf)
'''
# some static parameters
# the largest dataset allowed
LARGESTDATA = 8*1024*1024*1024  # 8 Giga points  (that's 64Gb !)
# the smallest axis size allowed
SIZEMIN = 1024
#HIGHMASS = 100000   # kludge to use mztoi !



def intelliround(x):
    "returns a number rounded to the nearest 'round' (easy to FT) integer"
    from math import log
    lx = int(log(x)/log(2))   # say 3
    err = 2*x
    r = 0
    for t in (2**lx, 2**(lx+1), 3*2**(lx-1), 5*2**(lx-2), 3*3*2**(lx-3),  7*2**(lx-2)):
        #      8      16          12          10           9            14  - increassing complexity
        if abs(t-x)<err:
            err = abs(t-x)
            r = t
            if err == 0: break
    return r
def pred_sizes_zf(d0, zf = 0, sizemin = SIZEMIN):
    """
    given an input data set, determines the optimum size s1,s2 to process it with a zerofilling of zf
    zf = +n is doubling n times along each axis
    zf = -n is halving n times along each axis
    zf = 0 is no zerofiling
    however, axes can never get smaller than sizemin
    returns (si1, si2, ...) as the dataset dimension
    """
    def dopow2(s, zf, sizemin):
        "return s * 2^zf rounded to the nearest 2^n or p*2^n"
        sm = min(sizemin, s)    # not smaller than sizemin, unless s is already smaller
        s1 = s*pow(2, zf)        # the zf level
        s = max(sm, s1)         # not smaller than sm
        return intelliround(s)
    def dopow(s, zf, sizemin):
        "do the math  not used"
        sm = min(sizemin, s)    # not smaller than sizemin, unless s is already smaller
        s1 = s*pow(2, zf)        # the zf level
        s = max(sm, s1)         # not smaller than sm
        return int(s)
    r = []      # will build into r
    for i in range(d0.dim):
        r.append( dopow2(d0.axes(i+1).size, zf, sizemin)) # apply l() for each axes
    if debug > 0: print r, reduce(lambda x, y:x*y, r)/1024/1024, 'Mpoint'     # report if debug
    return tuple(r) 
   
def pred_sizes(d0, szmult = (1,1), sizemin = SIZEMIN):
    """
    given an input data set, determines the optimum size s1,s2 to process it
        with a size multiplicant of szmult
    szmult (szm1, szm2)     where szm1 is multiplicant for s1 and szm2 for s2
    szmx = 1 : no change  /  2 : size doubling  /  0.5 : size halving
    any strictly positive value is possible, 0.2 0.33 1.1 2 2.2 5 etc...
    
    however, axes can never get smaller than sizemin
    returns (si1, si2, ...) as the dataset dimension
    """
    def dosize(s, szm, sizemin):
        "return s * 2^zf rounded to the nearest 2^n or 3*2^n"
        sm = min(sizemin, s)    # not smaller than sizemin, unless s is already smaller
        s1 = s*szm              # the szm level
        s = max(sm, s1)         # not smaller than sm
        return intelliround(s)
    r = []      # will build into r
    for i in range(d0.dim):
        r.append( dosize(d0.axes(i+1).size, szmult[i], sizemin)) # apply l() for each axes
    if debug > 0: print r, reduce(lambda x, y:x*y, r)/1024/1024, 'Mpoint'     # report if debug
    return tuple(r) 

def comp_sizes(d0,  zflist=None, szmlist=None, largest = LARGESTDATA, sizemin = SIZEMIN, vignette = True):
    """
    return a list with data-sizes, computed either
        zflist : from zerofilling index    eg : (1,0,-1)
        szmlist : from multiplicant pairs  eg : (2,2)

    largest determines the largest dataset allowed
    sizemini determines the minimum size when downzerofilling
    when vignette == True (default) a minimum size data (defined by sizemini) is appended to the list
    """
    sizes = []
    if (zflist!=None) and (szmlist!=None):
        raise Exception("Please define only one value : zerofilling or sizes multipliers")
    if (zflist==None) and (szmlist==None):
        zflist=[0]  # do at least something
    szres = []
    if (szmlist!=None):         # generate from parameters
        sm1,sm2 = szmlist
        while True:
            (si1,si2) = pred_sizes(d0, (sm1,sm2))
            print (sm1, sm2)
            print (si1, si2)
            if si1*si2 <= 1.5*SIZEMIN*SIZEMIN:
                break
            szres.append( (si1,si2))
            (sm1, sm2) = (0.25*sm1, 0.25*sm2)   # divide by 4, and insure floats
    if (zflist!=None):
        for zf in zflist:
            szres.append( pred_sizes_zf(d0, zf) )
    for (si1,si2) in szres:     # then verify
        while si1*si2 > largest:    # in size
            si2 /= 2
            print "Warning, reducing SI2 to %s"%si2
        sz = (si1,si2)
        if not sz in sizes:     # insure this size is not here yet ! 
            sizes.append(sz)
    if vignette:            # insure that last entry is a vignette, not smaller than sizemin
#        sz = (sizemin, sizemin)
        sz1,sz2 = sizes[-1]     # takes last entry
        while sz1 >= sizemin:
            sz1 /= 2
        while sz2 >= sizemin:
            sz2 /= 2
        if not (sz1,sz2) in sizes:
            sizes.append( (2*sz1, 2*sz2) )
    if debug>0: print "sizes to process", sizes
    return sizes

def apod(d, size, axis = 0):
    "apply sin 0.5 apodisation and change size"
    if d.size1<size:    # zerofilling
        d.apod_sin(maxi = 0.5, axis = axis)
        d.chsize(size)
    elif d.size1>size:
        d.chsize(size)  # truncating
        d.apod_sin(maxi = 0.5, axis = axis)
    else:
        d.apod_sin(maxi = 0.5, axis = axis)
    return d

def do_proc_F2(dinp, doutp):
    "scan all rows of dinp, apply proc() and store into doutp"
    size = doutp.axis2.size
    scan = min(dinp.size1, doutp.size1)      # min() because no need to do extra work !
    #scan = dinp.size1 # when was it done? 
    widgets = ['Processing F2: ', pg.Percentage(), ' ', pg.Bar(marker='-',left='[',right=']'), pg.ETA()]
    pbar= pg.ProgressBar(widgets=widgets, maxval=scan) #, fd=sys.stdout)
    print "############  in do_proc_F2 #########"
    print "dinp.axis1.itype ", dinp.axis1.itype 
    print "dinp.axis2.itype ", dinp.axis2.itype
    print "doutp.axis1.itype ", doutp.axis1.itype 
    print "doutp.axis2.itype ", doutp.axis2.itype
    print "dinp.axis1.size ", dinp.axis1.size 
    print "dinp.axis2.size ", dinp.axis2.size 
    print "doutp.axis1.size ", doutp.axis1.size 
    print "doutp.axis2.size ", doutp.axis2.size 
    print "########################### doutp.report() "
    print doutp.report()
    #print dir(doutp)
    for i in xrange(scan):
        # if i%(scan/16) == 0:                    # print avancement
        #     print "proc row %d / %d"%(i,scan)
        r = dinp.row(i)
        apod(r, size)
        r.rfft()    
        doutp.set_row(i,r)
        pbar.update(i+1)
        if interfproc:
            output = open('InterfProc/progbar.pkl', 'wb')
            pb = ['F2', int((i+1)/float(scan)*100)]
            pickle.dump(pb, output)
            output.close()
    pbar.finish()

def do_proc_F1(dinp, doutp):
    "scan all cols of dinp, apply proc() and store into doutp"
    size = doutp.axis1.size
    scan = min(dinp.size2, doutp.size2)      # min() because no need to do extra work !
    widgets = ['Processing F1: ', pg.Percentage(), ' ', pg.Bar(marker = '-',left='[',right=']'), pg.ETA()]
    pbar= pg.ProgressBar(widgets = widgets, maxval = scan) #, fd=sys.stdout)
    for i in xrange(scan):
        # if i%(scan/16) == 0:                    # print avancement
        #     print "proc col %d / %d"%(i,scan)
        c = dinp.col(i)
        apod(c, size)
        c.rfft()    
        doutp.set_col(i,c)
        pbar.update(i)
    pbar.finish()

def do_proc_F1_modu(dinp, doutp):
    "as do_proc_F1, but applies hypercomplex modulus() at the end"
    size = 2*doutp.axis1.size
    scan =  min(dinp.size2, doutp.size2)
    widgets = ['Processing F1 modu: ', pg.Percentage(), ' ', pg.Bar(marker = '-',left = '[',right=']'), pg.ETA()]
    pbar = pg.ProgressBar(widgets=widgets, maxval=scan) #, fd=sys.stdout)
    d = FTICRData( buffer = np.zeros((2*doutp.size1,2)) )     # 2 columns - used for hypercomplex modulus 
    for i in xrange( scan):
        # if i%(scan/16) == 0:                    # print avancement
        #     print "proc col %d / %d"%(i,scan)
        d.chsize(2*doutp.size1, 2)     # 2 columns - used for hypercomplex modulus 
        for off in (0,1):
            p = dinp.col(2*i+off)
            apod(p, size)
            p.rfft()    
            d.set_col(off,  p )
        d.axis1.itype = 1
        d.axis2.itype = 1
        d.modulus()
        doutp.set_col(i,d.col(0))
        pbar.update(i+1)
    pbar.finish()

def _do_proc_F1_flip_modu(data):
    "given a pair of columns, return the processed flipped FTed modulused column"
    c0, c1, shift, size, parameter = data
    d = FTICRData(buffer = np.zeros((c0.size1, 2)))     # 2 columns - used for hypercomplex modulus 
    d.set_col(0,  c0 )
    d.set_col(1,  c1 )
    d.axis1.itype = 0
    d.axis2.itype = 1
    if parameter.samplingfile is not None:      # NUS ?
        d.axis1.load_sampling(parameter.samplingfile)   # load it
        samp = d.axis1.get_sampling() # and store it aside
        if parameter.samplingfile_fake:    # samplingfile_fake allows to fake NUS on complete data-sets
            d.set_buffer(d.get_buffer()[samp])     # throw points
        d.zf()                          # add missing points by padding with zeros
    d.flipphase(0.0, 180*shift) # equivalent to  d.flip()  d.phase()  d.flop()
    if parameter.do_urqrd:
        d.urqrd(k = parameter.urqrd_rank, axis = 1)#
    if parameter.samplingfile is not None:      # NUS ?
        if not parameter.do_urqrd:
            d.axis1.set_sampling(samp)
            d.set_buffer( d.get_buffer()[samp])     # throw points created before after fliphase is done
        else:    # if NUS but not urQRd, put zeros on unmeasured points
            d.zf()

    apod(d, size, axis = 1)
    d.rfft(axis = 1)        # this rfft() is different from npfft.rfft() one !
    d.modulus()
    return d.col(0).get_buffer()   # return raw data

def iterarg(dinp, rot, size, parameter ):
    "an iterator used by the processing to allow  multiprocessing or MPI set-up"
    for i in range(0, dinp.size2, 2):
        c0 = dinp.col(i)
        c1 = dinp.col(i+1)
#        print i, c0, c1, rot, size
        yield (c0, c1, rot, size, parameter)

def do_proc_F1_flip_modu(dinp, doutp, parameter):
    "as do_proc_F1, but applies flip and then complex modulus() at the end"
    size = 2*doutp.axis1.size
    scan =  min(dinp.size2, doutp.size2)
    widgets = ['Processing F1 flip-modu: ', pg.Percentage(), ' ', pg.Bar(marker='-',left = '[',right = ']'), pg.ETA()]
    pbar= pg.ProgressBar(widgets = widgets, maxval = scan) #, fd=sys.stdout)
    shift = doutp.axis1.mztoi( doutp.axis1.highmass )   # frequency shift in points, computed from location of highmass
    hshift = doutp.axis1.itoh(shift)                    # the same in Hz
    if parameter.samplingfile is not None:                      #    NUS 
        dinp.axis1.load_sampling(parameter.samplingfile)       # load sampling file, and compute rot in non-NUS space
        cdinp = dinp.col(0)
        cdinp.zf()
        rot = cdinp.axis1.mztoi( dinp.axis1.highmass )
        print "11111111", shift, rot
        del(cdinp)
    else:                                                   # plain mode
        rot = dinp.axis1.mztoi( dinp.axis1.highmass )       # rot correction is applied in the starting space
    if debug>0: print "LEFT_POINT", shift
    doutp.axis1.left_point = shift      
    doutp.axis1.specwidth += hshift    # correction of specwidth
    xarg = iterarg(dinp, rot, size, parameter)      # construct iterator for main loop
    if parameter.mp:  # means multiprocessing //
        res = Pool.imap(_do_proc_F1_flip_modu, xarg)
        for i,buf in enumerate(res):
            doutp.buffer[:,i] = buf
#            doutp.set_col(i,p)
            pbar.update(i+1)
    elif mpiutil.MPI_size > 1:      # code for MPI processing //
        mpiutil.mprint('MPI NEW STYLE')
        res = mpiutil.enum_imap(_do_proc_F1_flip_modu, xarg)    # apply it
        for i,buf in res:       # and get results
            doutp.buffer[:,i] = buf
#            doutp.set_col(i, p)
            pbar.update(i+1)
            if interfproc:
                output = open('InterfProc/progbar.pkl', 'wb')
                pb = ['F1', int((i+1)/float(scan)*100)]
                pickle.dump(pb, output) # for Qt progressbar
                output.close()
        if interfproc:
            output = open('InterfProc/progbar.pkl', 'wb')
            pb = ['end']
            pickle.dump(pb, output) # for Qt progressbar
            output.close()
    else:       # plain non //
        res = itertools.imap(_do_proc_F1_flip_modu, xarg)
        for i,buf in enumerate(res):
            doutp.buffer[:,i] = buf
#            doutp.set_col(i,p)
            pbar.update(i+1)
    pbar.finish()

def do_process2D(dinp, datatemp, doutp, parameter):
    """
    apply the processing to an input 2D data set : dinp
    result is found in an output file : doutp
    
    dinp and doutp should have been created before, size of doutp will determine the processing
    will use a temporay file if needed
    """
    if debug>0:
        for f,d in ((parameter.infile, dinp), (parameter.interfile, datatemp), (parameter.outfile, doutp)):
            print "----------", f
            print d
    # in F2
    t00 = time.time()
    if parameter.do_F2:
        do_proc_F2(dinp, datatemp)
        print_time(time.time()-t00, "F2 processing time")
    # in F1
    if parameter.do_F1:
        t0 = time.time()
        if parameter.do_flip and parameter.do_modulus:
            do_proc_F1_flip_modu(datatemp, doutp, parameter)
        elif parameter.do_modulus:
            do_proc_F1_modu(datatemp, doutp)
        else:
            do_proc_F1(datatemp, doutp)
        print_time(time.time()-t0, "F1 processing time")
    # remove ridge computed on the last 10% rows
    if parameter.do_F1 and parameter.do_rem_ridge:      
        from util.rem_ridge import rem_ridge
        rem_ridge(doutp)
    if parameter.compress_outfile :  # fastclean is the trick for compression
        doutp.fastclean(nsigma=parameter.compress_level, axis=1)
    print_time(time.time()-t00, "F1-F2 processing time")

def downsample2D(data, outp, n1, n2):
    """
    takes data (a 2D) and generate a smaller dataset downsampled by factor (n1,n2) on each axis
    then returned data-set is n1*n2 times smaller
    - simply takes the mean
    ** Not fully tested on non powers of 2 **
    """
    for i in xrange(0, data.size1, n1):
        temp = np.zeros(data.size2/n2)
        for j in xrange(n1):
            if n2>1:
                yy = decimate(data.row(i+j).buffer, n2, ftype = "fir")   # filter along F2
            else:
                yy = data.row(i+j).buffer
            temp += yy
        outp.buffer[i/n1,:] = (1.0/n1)*temp
    copyaxes(data, outp)
    outp.axis1.left_point = outp.axis1.left_point/n1
    outp.axis2.left_point = outp.axis2.left_point/n2
    outp.adapt_size()
    return outp

def load_input(name):
    """load input file and returns it, in read-only mode"""
    if debug>0: print "reading", name
    hf = HDF5File(name, "r")    # reading the data 
    d0 = hf.load()
    d0.hdf5file = hf
    return d0

class Proc_Parameters(object):
    """this class is a container for processing parameters"""
    def __init__(self, configfile = None):
        "initialisation, see processe.mscf for comments on parameters"
        # processing
        self.do_F2 = True
        self.do_F1 = True
        self.do_modulus = True
        self.do_rem_ridge = True
        self.do_flip = True
        self.do_urqrd = False
        self.urqrd_rank = 20
        self.zflist = None
        self.szmlist = None
        self.mp = False
        self.nproc = 4
        # files
        self.apex = None
        self.format = None
        self.infile = None
        self.interfile = None
        self.outfile = None
        self.compress_outfile = True
        self.compress_level = 1.0
        self.samplingfile = None
        self.samplingfile_fake = False
        self.tempdir = "/tmp"
        self.largest = LARGESTDATA
        if configfile:
            self.load(configfile)
    def load(self, cp):
        "load from cp config file - should have been opened with ConfigParser() first"
        self.apex =    cp.get( "import", "apex")                                        # input file
        self.format =    cp.get( "import", "format")                                    # used format Apex or Solarix
        self.infile =  cp.get( "processing", "infile")                                  # input file
        self.interfile = cp.get( "processing", "interfile", None)                       # intermediatefile
        self.outfile = cp.get( "processing", "outfile")                                 # output file
        self.compress_outfile = cp.getboolean( "processing", "compress_outfile", str(self.compress_outfile))
        self.compress_level = cp.getfloat( "processing", "compress_level", self.compress_level)
        self.tempdir = cp.get( "processing", "tempdir", ".")                            # dir for temporary file
        self.samplingfile = cp.get( "processing", "samplingfile")
        self.samplingfile_fake = cp.getboolean( "processing", "samplingfile_fake",str(self.samplingfile_fake))
        self.largest = cp.getint( "processing", "largest_file", 8*LARGESTDATA)            # largest allowed file
        self.largest = self.largest/8                                                   # in byte in the configfile, internally in word
        self.do_modulus = cp.getboolean( "processing", "do_modulus", str(self.do_modulus))   # do_modulus
        self.do_flip = cp.getboolean( "processing", "do_flip", str(self.do_flip))       # do_flip
        self.do_urqrd = cp.getboolean( "processing", "do_urqrd", str(self.do_urqrd))    # do_urqrd
        self.urqrd_rank = cp.getint( "processing", "urqrd_rank", self.urqrd_rank)       # do_urqrd
        self.do_rem_ridge = cp.getboolean( "processing", "do_rem_ridge", str(self.do_rem_ridge))
        self.mp = cp.getboolean( "processing", "use_multiprocessing", str(self.mp))
        self.nproc = cp.getint( "processing", "nb_proc", self.nproc)
        self.do_F1 = cp.getboolean( "processing", "do_F1", str(self.do_F1))
        self.do_F2 = cp.getboolean( "processing", "do_F2", str(self.do_F2))
        # load zflist  or  szmlist
        zflist =  cp.get( "processing", "zerofilling", self.zflist)                              # get zf levels
        if zflist:
            self.zflist = [int(i) for i in zflist.split()]                              # got a string, change in a list of int
            self.zflist.sort()
            self.zflist.reverse()                                                       # sort and reverse zflist, so that biggest is first
        else:
            self.zflist = None
        szmlist =  cp.get( "processing", "sizemultipliers", self.szmlist)                          #  get szm levels
        if szmlist:
            self.szmlist = [float(i) for i in szmlist.split()]                          # got a string, change in a list of values
            print self.szmlist
        else:
            self.szmlist = None
        # verifications
        if cp.has_option("processing", "sizemulipliers"):   # that nasty bug was around once.
            raise Exception('Error on the name of sizemultiplier parameter, sizemuliplier instead of sizemultiplier')
        if not self.do_F1 and not self.do_F2:
            raise Exception("no processing !")
        if self.interfile is None and not (self.do_F1 and self.do_F2):
            raise Exception("Partial processing, without intermediate file")
        if self.do_F1 and self.do_flip and not self.do_modulus:
            raise Exception("do_flip but not do_modulus is not implemented !")
        for f1, f2 in ((self.infile, self.interfile), (self.interfile,self.outfile), (self.infile, self.outfile)):
            if f1 == f2:
                raise Exception("input and output files have the same name : %s - this is not possible"%f1)
        if self.samplingfile and self.do_urqrd:
            raise Exception("urQRd cannot be applied on NUS a data-set")
        if (self.zflist!=None) and (self.szmlist!=None):
            raise Exception("Please define only one value : zerofilling or sizes multipliers")
        if self.mp and mpiutil.MPI_size > 1:
            raise Exception("use_multiprocessing is not compatible with MPI")
        
    def report(self):
        "print a formatted report"
        print "------------ processing parameters ------------------"
        for i in dir(self):
            if not i.startswith('_'):
                v = getattr(self,i)
                if not callable(v):
                    print i, ' :', v
        print "-----------------------------------------------------"

########################################################################################
class Test(unittest.TestCase):
    """tests """        
    def test_intelli(self):
        "testing 'intelligent' rounding"
        r = []
        for i in range(16,33):
            r.append(intelliround(i))
        self.assertEqual(r, [16, 16, 18, 20, 20, 20, 24, 24, 24, 24, 24, 28, 28, 28, 32, 32, 32])
                        #   [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32]
    def test_zf(self):
        "testing zerofilling computation"
        print self.test_zf.__doc__
        d = FTICRData(dim = 2)
        d.axis1.size = 1024
        d.axis2.size = 10*1024-10    # 10240 - 10 = 10230
        sizes = comp_sizes(d, zflist=(1,0,-1))
        if SIZEMIN == 1024:
            self.assertEqual( sizes, [(2048, 20480), (1024, 10240), (1024, 5120), (1024, 1280)])
        sizes = comp_sizes(d, szmlist=(3, 1.5) )
        if SIZEMIN == 1024:
            self.assertEqual( sizes, [(3072, 14336), (1024, 3584), (1024, 1792)])
    def test_proc(self):
        "apply a complete processing test"
        main(["prgm", "test.mscf",])


########################################################################################
def Report_Table_Param():
    print "---------------SETTINGS---------------------"
    print "| MAX_COLUMNS ", tables.parameters.MAX_COLUMNS
    print "| MAX_NODE_ATTRS " , tables.parameters.MAX_NODE_ATTRS
    print "| MAX_GROUP_WIDTH ", tables.parameters.MAX_GROUP_WIDTH
    print "| MAX_UNDO_PATH_LENGTH ", tables.parameters.MAX_UNDO_PATH_LENGTH
    print "| CHUNK_CACHE_NELMTS ", tables.parameters.CHUNK_CACHE_NELMTS
    print "| CHUNK_CACHE_PREEMPT " , tables.parameters.CHUNK_CACHE_PREEMPT 
    print "| CHUNK_CACHE_SIZE ", tables.parameters.CHUNK_CACHE_SIZE 
    print "| METADATA_CACHE_SIZE ", tables.parameters.METADATA_CACHE_SIZE 
    print "| NODE_CACHE_SLOTS ", tables.parameters.NODE_CACHE_SLOTS
    print "| IO_BUFFER_SIZE ", tables.parameters.IO_BUFFER_SIZE
    print "| BUFFER_TIMES ", tables.parameters.BUFFER_TIMES 
    print "| EXPECTED_ROWS_EARRAY ", tables.parameters.EXPECTED_ROWS_EARRAY
    print "| EXPECTED_ROWS_TABLE ", tables.parameters.EXPECTED_ROWS_TABLE
    print "| PYTABLES_SYS_ATTRS ", tables.parameters.PYTABLES_SYS_ATTRS
    #print "| MAX_THREADS ",tables.parameters.MAX_THREADS 
    print "---------------SETTINGS---------------------"

def Set_Table_Param():
#    if debug>0: return
    tables.parameters.CHUNK_CACHE_PREEMPT = 1
    tables.parameters.CHUNK_CACHE_SIZE = 100*1024*1024
    tables.parameters.METADATA_CACHE_SIZE  = 100*1024*1024
    tables.parameters.NODE_CACHE_SLOTS = 100*1024*1024
    #tables.parameters.EXPECTED_ROWS_EARRAY = 100
    #tables.parameters.EXPECTED_ROWS_TABLE =100
    #tables.parameters.MAX_THREADS = 8
    #tables.parameters.PYTABLES_SYS_ATTRS = False

########################################################################################
def print_time(t,st="Processing time"):
    "prints processing time"
    h = int(t/3600)
    m = int ((t-3600*h)/60)
    s = int(t-3600*h - 60*m)
    print " %s : %d:%02d:%02d"%(st, h, m, s)


def main(argv = None):
    """
    Does the whole on-file processing, 
    syntax is
    processing.py [ configuration_file.mscf ]
    if no argument is given, the standard file : process.mscf is used.
    """
    TeeLogger(erase=True)
    global Pool     # This global will hold the multiprocessing.Pool if needed
    Pool = None
    t0 = time.time()
    t00 = t0
    ######### read arguments
    if not argv:
        argv = sys.argv
    try:                        # First try to read config file from arg list
        configfile = argv[1]
    except IndexError:          # then assume standard name
        configfile = "process.mscf"
    if interfproc:
        output = open('InterfProc/progbar.pkl', 'wb')
        pb = ['F2', 0]
        pickle.dump(pb, output)
        output.close()
        print "using %s as configuration file" %configfile
    #### get parameters from configuration file - store them in a parameter object
    cp = NPKConfigParser()
    cp.readfp(open(configfile))
    print "reading config file"
    param = Proc_Parameters(cp) # parameters from config file.. 
    # get optionnal parameters
    opt_param = {}
    for p in ("F1_specwidth", "F2_specwidth", "highmass", "ref_mass", "ref_freq"):
        v = cp.getfloat( "import", p, 0.0)
        if v != 0.0:
            opt_param[p] = v
    if param.mp:
        Pool = mp.Pool(param.nproc)     # if multiprocessing, creates slaves early, while memory is empty !
    param.report()
    ######## determine files and load inputfile
    ### input file either raw to be imported or already imported
    imported = False
    if not os.path.exists(param.infile):
        print "importing %s into %s"%(dir, param.infile)
        d0 = Import_2D[param.format](param.apex, param.infile)
        imported = True
        if opt_param != {}: # if some parameters were overloaded in config file
            # hum close, open, close, open ...
            d0.hdf5file.close()
            del(d0)
            hf = HDF5File(param.infile,"rw")
            for item in opt_param:
                if item.startswith('F1_'):
                    fileitem = item[3:]
                    hf.axes_update(axis = 1, infos = {fileitem:opt_param[item]})
                    print "Updating axis F1 %s to %f"%(fileitem, opt_param[item])
                elif item.startswith('F2_'):
                    fileitem = item[3:]
                    hf.axes_update(axis = 2, infos = {fileitem:opt_param[item]})
                    print "Updating axis F2 %s to %f"%(fileitem, opt_param[item])
                else:
                    hf.axes_update(axis = 1, infos = {item:opt_param[item]})
                    hf.axes_update(axis = 2, infos = {item:opt_param[item]})
                    print "Updating all axes %s to %f"%(item, opt_param[item])
            hf.close()
            d0 = load_input(param.infile)
    else:
        d0 = load_input(param.infile)
    d0.check2D()    # raise error if not a 2D
    # if param.samplingfile is not None:
    #     d0.axis1.load_sampling()
    if imported:
        print_time( time.time()-t0, "Import")
    else:
        print_time( time.time()-t0, "Load")
    Set_Table_Param()
    if debug>0:
        Report_Table_Param()
        print d0.report()
    ### compute final sizes
    allsizes = comp_sizes(d0, zflist=param.zflist, szmlist=param.szmlist, largest = param.largest)
    if debug>0: print allsizes
    (sizeF1, sizeF2) = allsizes.pop(0)   # this is the largest, to be processed by FT
    ### prepare intermediate file
    if debug>0: print "preparing intermediate file "
    if param.interfile is None:     # We have to create one !
        interfile = os.path.join(param.tempdir,'tmpfile_for_{}.msh5'.format(os.path.basename(param.outfile[:-4])))  
        print "creating TEMPFILE:",interfile
    else:
        interfile = param.interfile
    if param.do_F2:     # create
        temp =  HDF5File(interfile, "w")
        datatemp = FTICRData(dim = 2)
        copyaxes(d0, datatemp)
        if param.do_modulus:
            datatemp.axis1.size = min(d0.size1, sizeF1)
            datatemp.axis2.size = 2*sizeF2
        else:
            datatemp.axis1.size = min(d0.size1, sizeF1)
            datatemp.axis2.size = sizeF2
        temp.create_from_template(datatemp)
    else:                # already existing
        datatemp = load_input(param.interfile)
    ### prepare output file
    if debug>0: print "preparing output file "
    if param.do_F1:
        hfar =  HDF5File(param.outfile, "w") #, debug=debug)  # OUTFILE for all resolutions
        d1 = FTICRData( dim = 2 )   # create dummy 2D
        copyaxes(d0, d1)        # copy axes from d0 to d1
        d1.axis2.size = sizeF2
        d1.axis1.size = sizeF1
        group = 'resol1'
        if param.compress_outfile:    # file is compressed
            hfar.set_compression(True)
        hfar.create_from_template(d1, group)
        if debug>0:
            print "######################### d1.report() ################"
            print d1.report()
            print "######################### Checked ################"
    else:
        d1 = None
    print """
=============================
processing FT
============================="""
    t0 = time.time()
    do_process2D(d0, datatemp, d1, param) # d0 original, d1 processed
    # close input and temp file
    try:
        d0.hdf5file.close()
    except AttributeError:      # depends on how d0 was loaded
        pass
    datatemp.hdf5file.close()
    if param.do_F1:
        hfar.axes_update(group = group,axis = 1, infos = {'specwidth':d1.axis1.specwidth, 'left_point':int(d1.axis1.left_point)})
    if param.interfile is None:
        temp.close()
        os.unlink(interfile)
    print "==  FT Processing finished  =="
    print_time(time.time()-t0, "FT processing time")
    if param.do_F1:
        downprevious = d1       # used to downsample by step   downprevious -downto-> down
        t0 = time.time()
        for (i, (sizeF1, sizeF2)) in enumerate(allsizes):
            if (downprevious.size1%sizeF1) != 0 or  (downprevious.size2%sizeF2) != 0:
                print "downsampling not available for level %d : %d x %d -> %d x %d"%((i+1), downprevious.size1, downprevious.size2, sizeF1, sizeF2)
                continue
            zflevel = "level %d"%(i+1)
            print """
================
downsampling %s
================""" % zflevel
            group = 'resol%d'%(i+2)     # +2 because we poped the first value
            if debug > 0: print "downsampling", group, (sizeF1, sizeF2)
            down = FTICRData( dim = 2 )   # create dummy 2D
            copyaxes(d1, down)        # copy axes from d0 to d1
            down.axis1.size = sizeF1
            down.axis2.size = sizeF2
            #create_branch(hfar, group, d1)
            hfar.create_from_template(down, group)
            if debug > 0: print down
            downsample2D(downprevious, down, int(downprevious.size1/sizeF1), int(downprevious.size2/sizeF2))
            hfar.axes_update(group = group, axis = 1, infos = {'left_point': down.axis1.left_point})
            downprevious = down
        print_time(time.time()-t0, "Downsampling time")
    # close output file
    hfar.close()
    print "== Processing finished  =="
    print_time(time.time() - t00, "Total processing time")
    if param.mp:
        Pool.close()    # finally closes multiprocessing slaves

# default values
if __name__ == '__main__':
    mp.freeze_support()
    if mpiutil.MPI_size < 2:            # this is single processor
        main()
    else:                       # this is a MPI run
        if mpiutil.MPI_rank == 0:   # master proc
            main()
            mpiutil.shutdown()
        else:               # slave proc
            mpiutil.slave()
