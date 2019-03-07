#!/usr/bin/env python 
# encoding: utf-8

"""
    Utility to Handle NMR Spinit files

"""

from __future__ import print_function

__author__ = "Marc Andre' Delsuc"
__date__ = "August 2017"

import os
import os.path as op
import glob
import re
import struct
import time
import unittest
import shutil
import xml.etree.ElementTree as ET

import numpy as np

from ..NPKData import NPKData

debug = False
verbose = True   # change this for verbose importers

################################################################
def load_header(filename="header.xml"):
    """ 
    loads a header.xml as an ET.tree, and keep it in memory

    returns headertree
    xml parameters are found in headertree.getroot()
    
    """
    tree = ET.parse(filename)
    return tree

def read_param(filename="header.xml"):
    """ 
    loads a header.xml as a dictionnary key:value
    
    """
    debug = 0
    headertree = load_header(filename)
    xroot = headertree.getroot()
    xparams = xroot[0]
    params = {}
    for entry in xparams:
        k = entry.find('key').text
        params[k] = ""
        for xv in entry.iter('value'):
            try:
                params[k] += xv.text.strip()
            except AttributeError:
                pass
# debug code
    if debug:
        for i in dict.keys():
            print(i+" = "+str(dict[i]))
    return params

def modify_val(headertree, key, value):
    """
    modify an entry from a headertree loaded with load_header()
    key should be present in headertree
    key and val are both strings

    headertree can then written back to disk with
        headertree.write('header.xml')
    """
    debug = 1
    for entry in headertree.getroot()[0]:
        k = entry.find('key').text
        if k == key:                  # match
            vals = entry.find('value')
            bef = vals.find('value').text
            vals.find('value').text = value
            if debug:
                print (key, 'Before:', bef, 'After:', vals.find('value').text)

################################################################
def read_1D(size, filename="data.dat"):
    """
    Reads in a Spinit 1D fid as a numpy float array
    
    size is the number of data-points in the fid
    uses struct
    does not check endianess
    """
    print("size is ", size)
    fmt = ">%df"
    si = size
    rl = np.empty((si/2))
    im = np.empty((si/2))
    print("****** bits number is  {0} ".format(8*si))
    dt = np.dtype(">f")
    ibuf = np.fromfile(filename, dt)
    # with open(filename,"rb") as F:
    #     buf = F.read(8*si)
    # ibuf = struct.unpack(fmt%(2*si), buf)  # unpack the string as pairs of float
    print("passing in complex")
    print("rl.size ",rl.size)
    print("im.size ",im.size)
    print("ibuf.size ",ibuf.size)
    rl[:] = ibuf[::2]
    im[:] = ibuf[1::2]
    npkbuf = rl + 1j*im
    print("complex is OK")
    return npkbuf
    
################################################################
def read_2D(sizeF1, sizeF2, filename="data.dat"):
    return npkbuf

################################################################
def zerotime(acqu):
    """get digital filter parameters, if any
    """
    try:
        zerotimeposition = float( acqu['DIGITAL_FILTER_SHIFT'] )
    except:
        zerotimeposition = 0.0
    return(zerotimeposition)
################################################################
def offset(acqu, proc):
    """
    computes the offset from spinit to spike
    """
    try:
        ppmPointUn = float(proc['$OFFSET'])
        ppmWidth = float(acqu['$SW_h']) / float(acqu['$SFO1'])
        calibrationOffset = float(ppmPointUn - ppmWidth)*  float(acqu['$SFO1'])
    except:
        calibrationOffset=0.0
    return (calibrationOffset)
################################################################
def Import_1D(filename="data.dat"):
    """
    Imports a 1D spinit fid as a NPKData
    
    """
    if (not op.exists(filename)):
        raise Exception(filename+" : file not found")
    dire=op.dirname(filename)
    acqu = read_param(op.join(dire,'header.xml'))
    size= int(acqu['ACQUISITION_MATRIX_DIMENSION_1D'])  # get size
    if verbose: print("importing 1D FID, size =",size)
    data = read_1D(size, filename)
    print("After read_1D, size is : ", data.size)
    d = NPKData(buffer=data)
#    d.revf()
# then set parameters
    d.axis1.specwidth = float(acqu['SPECTRAL_WIDTH'])           # in Hz
    d.axis1.frequency = float(acqu['OBSERVED_FREQUENCY'])/1E6   # in MHz
    d.axis1.itype = 1      # always complex
    d.axis1.offset = float(acqu['SR']) - d.axis1.specwidth/2

    d.axis1.zerotime = zerotime(acqu)
    pardic = {"header": acqu} # create ad-hoc parameters
    d.params = pardic   # add the parameters to the data-set
    return d

################################################################
def Export_1D(d, filename="data.dat", template="header.xml", kind=None):
    """
    export a 1D NPKData as a spinit
    kind: 1DFID, 1DSPEC
    """
    if d.dim >1:
        raise Exception("not implemented yet for Dim>1")
    if d.itype == 1:     # complex
        array = d.get_buffer()
    else:                   # real
        array = d.get_buffer() + 0.0j   # add a zero imaginary part

    dirname = op.dirname(filename)
    if not op.isdir(dirname):
        os.makedirs(dirname)
    # array.astype(">f").tofile(filename)
    array.view(np.float_).astype('>f').tofile(filename)

# then set parameters
    headertree = load_header(template)
    modify_val(headertree, 'MODALITY', "NMR" )
    modify_val(headertree, 'ACCU_DIM', "1" )
    modify_val(headertree, 'ACQUISITION_MATRIX_DIMENSION_1D', str(d.axis1.cpxsize) )
    modify_val(headertree, 'SPECTRAL_WIDTH', str(d.axis1.specwidth) )
    modify_val(headertree, 'OBSERVED_FREQUENCY', str(d.axis1.frequency*1E6) )
    modify_val(headertree, 'SR', str(d.axis1.offset+d.axis1.specwidth/2) )
    print ('    **SR is wrong ! to-be-changed**')
    try:
        STATE = [str(d.axis1.state),"0", "0", "0"]
    except:
        STATE = d.params['header']['STATE']
    modify_val(headertree, 'STATE',  STATE)
    modify_val(headertree, 'DIGITAL_FILTER_SHIFT', str(d.axis1.zerotime) )
    headerfile = op.join( dirname, "header.xml")
    headertree.write(headerfile)

if __name__ == '__main__':
    unittest.main()
