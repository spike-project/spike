#!/usr/bin/env python 
# encoding: utf-8

"""
big.py

First try for off-memory processing, using pytables

Created by MAD on 2011-04-13.
Copyright (c) 2011 IGBMC. All rights reserved.

"""
 
from __future__ import print_function
import Apex as ap
import FTICR
import tables
import array
import time
import numpy as np
import downsample
import Cadzow
import pickle
import File.HDF5File as hf5
tables.parameters.CHUNK_CACHE_PREEMPT = 1
tables.parameters.CHUNK_CACHE_SIZE = 16844808


###########################################################
# This code allows to pickle methods, and thus to use multiprocessing on methods
# This very nice trick comes from :
# Steven Bethard
# http://bytes.com/topic/python/answers/552476-why-cant-you-pickle-instancemethods
#
def _pickle_method(method):
    
   func_name = method.im_func.__name__
   obj = method.im_self
   cls = method.im_class
   if func_name.startswith('__') and not func_name.endswith('__'):
       cls_name = cls.__name__.lstrip('_')
   if cls_name: func_name = '_' + cls_name + func_name
   
   return _unpickle_method, (func_name, obj, cls)

def _unpickle_method(func_name, obj, cls):
   for cls in cls.mro():
       try:
           func = cls.__dict__[func_name]
       except KeyError:
           pass
       else:
           break
   return func.__get__(obj, cls)
import copy_reg
import types
copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)
# end of the trick
###########################################################


def impbig():
    Apex.Ser2D_to_H5f(8*1024, 512*1024, filename="/Echange/big/ser",outfile="/Echange/big/ser.h5",chunks = (2*8,2*512))

def new_tb(outfile, sizeF1, sizeF2, chunkshape):
    """
    creates empty pytable file
    """
    h5f = hf5.HDF5File(outfile,"w")
    group_resol = h5f.createGroup("/", 'resol1')
    h5f.createCArray("/resol1",'data', tables.Float64Atom(), (sizeF1,sizeF2), chunk = chunkshape)
    return h5f.hf.root.resol1.data

def ft_2_mp_V1(inf, outf):
    i = 0
# 1ere version
    print("type inf",type(inf))
    print("type ",type(inf.iter_all_row()))
    for r in inf.iter_all_row():
        if i%64 == 0:
            print("\nproc row %d"%(i))
        i += 1
        r.apod_sin(maxi=0.5).chsize(inf.size2).rfft()
        outf.set_row(i,r)
def func(row):
    return row.apod_sin(maxi=0.5).chsize(row.size1).rfft()
    
from string import ascii_lowercase
ascii=ascii_lowercase
#print ascii.next()

import threading
posth=0
def threadfunc(row):
    ascii[posth]=threading.Thread(None,func,None,(row),{})# thread
    ascii[posth].start()
    posth+=1
    
    
def func2(arg):
    (data, i) = arg
    return data.row(i).apod_sin(maxi=0.5).chsize(row.size1).rfft()
    
def ft_2_mp_V2(inf, outf, N_proc=8):
    #2eme version
    import multiprocessing as mp
    import itertools

    print("type inf ", type(inf))
    print("inf.size1 ",inf.size1)
    pool = mp.Pool(processes=N_proc)
    print("type inter ", type(inf.iter_all_row))
    #print "pickle inf"
    #fpic=open("/Echange/big/picpic",'w')
    #pickle.dump(inf,fpic)
    gen=inf.iter_all_row
    #result = pool.imap(func, gen)
    result = itertools.imap(threadfunc, gen)
    i = 0
    for i in range(3):
        result
    '''
    for d1D in result:
        if i%64 == 0:
            print "\nproc row %d"%(i)
        outf.set_row(i,d1D)
        i += 1
    '''
def ft_2_mp_V3(inf, outf, N_proc=8):
    #3eme version
    import multiprocessing as mp
    import itertools as it
        
    print("type inf ", type(inf))
    print("inf.size1 ",inf.size1)
    pool = mp.Pool(processes=N_proc)
    args = it.izip( it.repeat(inf), xrange(inf.size1) )
    #result = pool.imap(func2, args)
    result = itertools.imap(func, gen)
    i = 0
    for d1D in result:
        if i%64 == 0:
            print("\nproc row %d"%(i))
        outf.set_row(i,d1D)
        i += 1
  
            
def ft_2(inf, outf):
    for i in xrange(inf.size1):
        if i%64 == 0:
            print("\nproc row %d"%(i), end=' ')
        else:
            print(".", end=' ')
        r = inf.row(i).apod_sin(maxi=0.5).chsize(inf.size2).rfft()
        outf.set_row(i,r)

def ft_1_mod(inf,outf):
    print("ft_1-mod")
    print("size hauteur",inf.size2)
    for i in xrange(0,inf.size2,2):
        if i%64 == 0:
            print("\nproc col %d"%(i), end=' ')
        else:
            print(".", end=' ')
        d = FTICR.FTICRData( buffer=np.zeros((2*inf.size1,2)) )     # used for hypercomplex modulus
        d.set_col(0, inf.col(i).apod_sin(maxi=0.5).chsize(2*inf.size1).rfft())
        d.set_col(1, inf.col(i+1).apod_sin(maxi=0.5).chsize(2*inf.size1).rfft())
        d.axis1.itype = 1
        d.axis2.itype = 1
        d.modulus()
        outf.set_col(i/2,d.col(0))

def procbig(infile="/Echange/big/ser.h5", tmpfile="/Echange/big/tmp.h5", outfile="/Echange/big/smx.h5", chunks=(16, 1024)):
    """
    do the 2D FT on file, using tmpfile as temporaty storage
    all files are pytables
    """
    
    t0 = time.time()
    h = hf5.HDF5File(infile,"rw")
    inf = FTICR.FTICRData(buffer=h.hf.root.resol1.data)   # fichier d'entree
    print("inf : ")
    print("end inf")

    tmpb = new_tb(tmpfile, inf.size1, inf.size2, chunkshape=chunks)
    tmp = FTICR.FTICRData(buffer=tmpb)
    print("tmp : ")
    

    outb = new_tb(outfile, inf.size1, inf.size2/2, chunkshape=chunks)
    out = FTICR.FTICRData(buffer=outb)
    print("out :",out)

    print("Avant FT2 ", time.time()-t0)
    print("exec ft_2mp")
    #ft_2_mp_V2(inf, tmp)
    ft_2(inf, tmp)
    #ft_2_mp_V1(inf, tmp)

    print("Avant FT1 ", time.time()-t0)
    ft_1_mod(tmp, out)

    print("Total ", time.time()-t0)

def clean(d1):
    """clean d1 by removing constant columns - averaged on last 20 points"""
    c = d1.row(d1.size1-30)     # compute correction
    for i in xrange(19):
        dd = d1.row(d1.size1-29+i)
        c.buffer -= dd.buffer
    c.buffer *= 1.0/20
    for i in xrange(d1.size1):
        dd = d1.row(i)
        dd.buffer -= c.buffer
        d1.set_row(i,dd)        # then apply
    return d1

def down(inf="/Echange/big/smx.h5", out="downsampled.gs2"):
    import downsample
    h = hf5.HDF5File(inf,"r")
    inf = FTICR.FTICRData(buffer=h.hf.root.resol1.data)   # fichier d'entree
    inf.specwidth = 1000*500.0
    inf.highmass = 2000
    inf.ref_mass = 344.0974          # reference mass
    inf.ref_freq =  1000*419.62      # freq for the reference mass    
    out = downsample.downsample2D(inf, n1=8, n2=16)
    out.specwidth = 1000*500.0
    out.highmass = 2000
    out.ref_mass = 344.0974          # reference mass
    out.ref_freq =  1000*419.62      # freq for the reference mass    
    out.save(out)
    
def extract(infile, outfile, chunkshape):
    h = hf5.HDF5File(infile,"r")
    inf = FTICR.FTICRData(buffer=h.hf.root.resol1.data)   # fichier d'entree
    inf.specwidth = 1000*500.0
    inf.highmass = 2000
    inf.ref_mass = 344.0974          # reference mass
    inf.ref_freq =  1000*419.62      # freq for the reference mass
    x1 = inf.axis1.mztoi(288)
    x2 = inf.axis2.mztoi(288)
    print("je coupe : (%d %d) x (%d %d)"%(x1, inf.size1, x2, inf.size2))
    outb = new_tb(outfile, inf.size1-x1-1, inf.size2-x2-1, chunkshape=chunkshape)
    out = FTICR.FTICRData(buffer=outb)
    for i in xrange(x1, inf.size1):
        if i%64 == 0:   print("\nproc row %d"%(i), end=' ')
        r = inf.row(i)
        r.buffer = r.buffer[x2:]
        r.adapt_size()
        out.set_row(i-x1, r)

def zoomz(data, scale, mz1, mz2):
    """2D zoom in m/z from mz1 to mz2"""
    if mz1 < mz2:
        (mz2, mz1) = (mz1, mz2)
    # le zoom doit -être donné en index croissants - donc en m/z décroissant
    z = [ [ data.axis1.mztoi(mz1), data.axis1.mztoi(mz2) ],
        [ data.axis2.mztoi(mz1), data.axis2.mztoi(mz2) ]]
    data.display(scale=scale, zoom = z)
    

if __name__ == '__main__':
    t0 = time.time()
    #Apex.Ser2D_to_H5f(2*1024, 64*1024, filename="/Echange/FT-ICR-Cyto/cytoC_2D_000001.d/ser", outfile="/Echange/big/ser1.h5",chunks = (2*8,512))
    #procbig(infile="/Echange/big/ser1.h5", tmpfile="/Echange/big/tmp1.h5", outfile="/Echange/big1/smx.h5", chunks=(2*8,512))
#    Apex.Ser2D_to_H5f(8*1024, 512*1024, filename="/Echange/big/ser", outfile="/Echange/big/ser.h5",chunks = (2*8,2*512))
#    print "Import ", time.time()-t0
    d = ap.Import_2D("/Volumes/XeonData/Developement/MS-FTICR/cytoC_2D_000001.d", "/Users/mac/Desktop/Chunk.hf")
    procbig(infile="/Users/mac/Desktop/Chunk.hf", tmpfile="/Users/mac/Desktop/Tmp.hf", outfile="/Users/mac/Desktop/Smx.hf", chunks=(257,8193))
    #down(inf="/Echange/big/smx2.h5", out="downsampled2.gs2")
    print("Total Total ", time.time()-t0)
