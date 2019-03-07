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

opd = os.path.dirname
opb = os.path.basename
opj = os.path.join

import numpy as np

from ..NPKData import NPKData, NPKData_plugin

debug = False
verbose = False   # change this for verbose importers

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
    headertree = load_header(filename)
    xroot = headertree.getroot()
    params = {}
    for entry in xroot.findall(".//entry"):
        k = entry.find('key').text
        paramsk = [kk.text for kk in entry.findall('./value/value')]  # yes 'value/value' !
        if len(paramsk) == 1 :   # singleton list  make it string
            params[k] = paramsk[0]
        else:
            params[k] = paramsk
    return params

def has_param(headertree, key):
    '''
    Checks if the parameter (key) exists
    '''
    root = headertree.getroot()
    entry = root.find(".//entry/[key='%s']"%(key))  
    return entry != None
    
def modify_val(headertree, key, values):
    """
    modify an entry from a headertree loaded with load_header()
    key should be present in headertree
    values is either a single value or a list, depending on key type, and lengths should match
    key and val are both strings

    headertree can then written back to disk with
        headertree.write('header.xml')
    """
    root = headertree.getroot()
    entry = root.find(".//entry/[key='%s']"%(key))
    print(entry)
    xmlval = entry.find('./value')
    if type(values) is list:
        _vals = values
    else:
        _vals = [values,]
    if len(xmlval.findall('./value')) != len(_vals):
        print (len(xmlval.findall('./value')) , len(_vals))
        raise Exception('XML value lengths do not match')
    for p,v in zip(xmlval.findall('./value'), _vals):
        p.text = str(v)

def state(val):
    '''
    Produces the entry block for STATE from the given list val. 
    '''
    return ET.fromstring('''
                
        <entry>
            <key>STATE</key>
            <value xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:type="listNumberParam">
                <name>STATE</name>
                <description>STATE.description</description>
                <displayedName>STATE</displayedName>
                <locked>true</locked>
                <lockedToDefault>false</lockedToDefault>
                <group>Miscellaneous</group>
                <category>Miscellaneous</category>
                <rolesEnum>User</rolesEnum>
                <value>{0}</value>
                <value>{1}</value>
                <value>{2}</value>
                <value>{3}</value>
                <defaultValue>0</defaultValue>
                <defaultValue>0</defaultValue>
                <defaultValue>0</defaultValue>
                <defaultValue>0</defaultValue>
                <maxValue>2147483647</maxValue>
                <minValue>-2147483648</minValue>
                <numberEnum>Integer</numberEnum>
            </value>
         </entry>
      
        '''.format(*val))

def add_state(headertree, value):
    """
    add an entry from a headertree loaded with load_header()
    value is a list of strings
    """
    head = headertree.getroot()[0]
    head.append(state(value))
     
def data_representation(val):
    '''
    Produces the entry block for DATA_REPRESENTATION from the given val (type is list). 
    '''
    return ET.fromstring('''
       <entry>
            <key>DATA_REPRESENTATION</key>
            <value xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:type="listTextParam">
                <name>DATA_REPRESENTATION</name>
                <displayedName>DATA_REPRESENTATION</displayedName>
                <locked>true</locked>
                <lockedToDefault>false</lockedToDefault>
                <group>Dimension</group>
                <category>Miscellaneous</category>
                <rolesEnum>User</rolesEnum>
                <value>{0}</value>
                <value>{1}</value>
                <value>{2}</value>
                <value>{3}</value>
                <defaultValue>COMPLEX</defaultValue>
                <defaultValue>REAL</defaultValue>
                <defaultValue>REAL</defaultValue>
                <defaultValue>REAL</defaultValue>
            </value>
        </entry>      
        '''.format(*val))

def add_data_representation(headertree, value):
    """
    add an entry from a headertree loaded with load_header()
    value is a list of strings
    """
    head = headertree.getroot()[0]
    head.append(data_representation(value))

################################################################
def read_1D(size, filename="data.dat", debug=0):
    """
    Reads in a Spinit 1D fid as a numpy float array
    
    size is the number of data-points in the fid
    uses struct
    does not check endianess
    """
    if debug>0: print("size is ", size)
    fmt = ">%df"
    si = size
    if debug>0: print("****** bits number is  {0} ".format(8*si))
    dt = np.dtype(">f")
    ibuf = np.fromfile(filename, dt)
    try:
        rl = np.empty((si/2))
        im = np.empty((si/2))
        rl[:] = ibuf[::2]
        im[:] = ibuf[1::2]
    except:
        rl = np.empty((si))
        im = np.empty((si))
        rl[:] = ibuf[::2]
        im[:] = ibuf[1::2]
    if debug>0:
        print("passing in complex")
        print("rl.size ",rl.size)
        print("im.size ",im.size)
        print("ibuf.size ",ibuf.size)
    npkbuf = rl + 1j*im
    if debug>0: print("complex is OK")
    return npkbuf
    
################################################################
def read_2D(sizeF1, sizeF2, filename="data.dat"):
    '''
    Reads the 2D files and return a buffer
    '''
    import struct
    npkbuf = np.empty((sizeF1, sizeF2), dtype='complex') # , dtype='complex'
    data_file = open(filename, mode='rb')

    for j in range(0, sizeF1):
        for i in range(0, sizeF2):
            re, im = struct.unpack('>ff', data_file.read(8)) # do something here
            npkbuf[j,i] = re+1j*im
                        
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
    Computes the offset from spinit to spike
    """
    try:
        ppmPointUn = float(proc['$OFFSET'])
        ppmWidth = float(acqu['$SW_h']) / float(acqu['$SFO1'])
        calibrationOffset = float(ppmPointUn - ppmWidth)*  float(acqu['$SFO1'])
    except:
        calibrationOffset=0.0
    return (calibrationOffset)
################################################################
def get_acquisition_mode(acqu):
    """
    Retrieves the acquisition mode list from the header file. 
    known values:
    REAL,
    COMPLEX,
    TPPI,
    COMPLEX_TPPI,
    PHASE_MODULATION,
    ECHO_ANTIECHO
    """
    if 'ACQUISITION_MODE' in acqu:
        # multi valued list, keep only F1 mode
        return acqu['ACQUISITION_MODE'][1]
    elif 'Phase Mod' in acqu:
        return acqu['Phase Mod']
    return None
################################################################
def get_data_representation(acqu):
    """
    Retrieves the data representation list from the header file. 
    known values:
    REAL,
    COMPLEX
    """
    if 'DATA_REPRESENTATION' in acqu:
        return acqu['DATA_REPRESENTATION']
    else:
        mode = get_acquisition_mode(acqu)        
        if mode in ('STATES', 'COMPLEX', 'COMPLEX_TPPI'):
            return ['COMPLEX', 'COMPLEX', 'REAL', 'REAL']
    
    return ['COMPLEX', 'REAL', 'REAL', 'REAL']
    
def ftF1_spinit(data, debug=0):
    '''
    Plugin spinit for perfoming the FT along the F1 axis according to the kind of acquisition.
    Cases taken in account : 
    TPPI
    COMPLEX
    PHASE_MODU
    COMPLEX_TPPI
    ECHO_ANTIECHO
    '''
    mode = get_acquisition_mode(data.params['header'])
    if debug>0: print("mode is ", mode)
    if mode == 'TPPI':
        data.revf(axis='F1')
        return data.ft_tppi()
    elif mode in ('STATES', 'COMPLEX'):
        data.revf(axis='F1')
        return data.ft_sh()
    elif mode in ('PHASE_MODU', 'PHASE_MODULATION', 'NONE', None):
        return data.ft_phase_modu().reverse(axis='F1')
    elif mode == 'COMPLEX_TPPI':
        data.revf(axis='F1')
        return data.ft_sh_tppi()
    elif mode == 'ECHO_ANTIECHO':
        data.revf(axis='F1')
        return data.ft_n_p()  
    if debug>0: print("no case encountered return a phase_modu processing")
    # cas par defaut: PHASE_MODULATION
    return data.ft_phase_modu().reverse(axis='F1')
    
NPKData_plugin('ftF1_spinit', ftF1_spinit)    
                
################################################################
def Import_1D(filename="data.dat"):
    """
    Imports a 1D spinit fid as a NPKData
    
    """
    if (not op.exists(filename)):
        raise Exception(filename+" : file not found")
    dire=op.dirname(filename)
    acqu = read_param(op.join(dire,'header.xml'))
    size= int(acqu['MATRIX_DIMENSION_1D'])  # get size
    if verbose: print("importing 1D FID, size =",size)
    data = read_1D(size, filename)
    print("After read_1D, size is : ", data.size)
    d = NPKData(buffer=data)
# then set parameters
    d.axis1.specwidth = float(acqu['SPECTRAL_WIDTH'])           # in Hz
    d.axis1.frequency = float(acqu['OBSERVED_FREQUENCY'])/1E6   # in MHz
    d.axis1.itype = 1      # always complex
    d.axis1.offset = float(acqu['SR']) - d.axis1.specwidth/2
    d.digital_filter_removed = bool(acqu['DIGITAL_FILTER_REMOVED'])

    d.axis1.zerotime = zerotime(acqu)
    pardic = {"header": acqu}          # create ad-hoc parameters
    d.params = pardic                  # add the parameters to the data-set
    return d

def Import_2D(filename="data.dat"):
    """

    Imports a 2D spinit fid as a NPKData
    
    """
    if (not op.exists(filename)):
        raise Exception(filename + " : file not found")
    dire = op.dirname(filename)
    acqu = read_param(op.join(dire,'header.xml'))
    size1 = int(acqu['MATRIX_DIMENSION_2D'])  # get size
    size2 = int(acqu['MATRIX_DIMENSION_1D'])  # get size
    if verbose:
        print("size1", size1)
        print("size2", size2)
    
    if int(acqu['MATRIX_DIMENSION_3D']) != 1:
        print("MATRIX_DIMENSION_3D != 1, ignoring")
    if int(acqu['MATRIX_DIMENSION_4D']) != 1:
        print("MATRIX_DIMENSION_4D != 1, ignoring")
        
    data = read_2D(size1, size2, filename)
    print (data.shape)
    d = NPKData(buffer=data, dim=2)
    d.axis1.specwidth = float(acqu['SPECTRAL_WIDTH_2D'])               # in Hz
    try:
        d.axis1.frequency = float(acqu['OBSERVED_FREQUENCY_2D'])/1E6   # in MHz
    except:
        print('OBSERVED_FREQUENCY_2D is missing, fold back to OBSERVED_FREQUENCY')
        d.axis1.frequency = float(acqu['OBSERVED_FREQUENCY'])/1E6      # in MHz
        
    representation = get_data_representation(acqu)        
    if representation[1] == 'REAL':
        d.axis1.itype = 0      #  real
    elif representation[1] == 'COMPLEX':
        d.axis1.itype = 1      # complex
    else:
        print('DATA_REPRESENTATION unknown: ', mode)
        d.axis1.itype = 0      #  real
    ###
    d.axis1.offset = -d.axis1.specwidth/2 + float(acqu['OFFSET_FREQ_2']) +float(acqu['SR'][1])

    d.axis2.specwidth = float(acqu['SPECTRAL_WIDTH'])           # in Hz
    d.axis2.itype = 1      # always complex
    d.axis2.offset = -d.axis2.specwidth/2 + float(acqu['OFFSET_FREQ_1']) +float(acqu['SR'][0])
#    d.axis1.offset = d.axis2.offset + d.axis2.specwidth/2 - d.axis1.specwidth/2  #### Only true in Homonuclear

    d.axis2.zerotime = zerotime(acqu)
    pardic = {"header": acqu}          # create ad-hoc parameters
    d.params = pardic                  # add the parameters to the data-set
    return d

################################################################
def Export_1D(d, filename="data.dat", template="header.xml", kind=None):
    """
    export a 1D NPKData as a spinit
    kind: 1DFID, 1DSPEC
    """
    if d.dim >1:
        raise Exception("not implemented yet for Dim>1")
    array = d.get_buffer() # + 0.0j   # add a zero imaginary part if lacking

    dirname = op.dirname(filename)
    if not op.isdir(dirname):
        os.makedirs(dirname)
    array.view(np.float_).astype('>f').tofile(filename)
# then set parameters
    headertree = load_header(template)
    modify_val(headertree, 'MODALITY', "NMR" )
    modify_val(headertree, 'ACCU_DIM', "1" )
    modify_val(headertree, 'MATRIX_DIMENSION_1D', str(d.axis1.cpxsize) )
    modify_val(headertree, 'SPECTRAL_WIDTH', str(d.axis1.specwidth) )
    modify_val(headertree, 'OBSERVED_FREQUENCY', str(d.axis1.frequency*1E6) )
    modify_val(headertree, 'SR', str(d.axis1.offset+d.axis1.specwidth/2) )
    print ('    **SR is wrong ! to-be-changed**')
    modify_val(headertree, 'DIGITAL_FILTER_SHIFT', str(d.axis1.zerotime) )
    headerfile = op.join( dirname, "header.xml")
    headertree.write(headerfile)

################################################################
def Export_2D(d,  filename="data.dat", template="header.xml", kind=None, debug=0):
    """
    export a 2D NPKData as a spinit
    """
    array = d.get_buffer().flatten() + 0.0j
    if debug>0:
        print( "type is ",  type(d.get_buffer()))
        print("array.size ", array.size)
    dirname = op.dirname(filename)
    if not op.isdir(dirname):
        os.makedirs(dirname)
    array.view(np.float_).astype('>f').tofile(filename)
    if debug>0: array.view(np.float_).astype('>f').tofile('z:/bordel/data.dat')
    # then set parameters
    headertree = load_header(template)
    modify_val(headertree, 'MATRIX_DIMENSION_1D', str(d.axis2.cpxsize) ) 
    modify_val(headertree, 'MATRIX_DIMENSION_2D', str(d.axis1.size) ) 
    modify_val(headertree, 'MATRIX_DIMENSION_3D', "1" ) 
    modify_val(headertree, 'MATRIX_DIMENSION_4D', "1" ) 
       
    ####
    mode = get_acquisition_mode(d.params['header'])   # retrieving acquisition mode
    if debug>0: print("mode is ", mode)
    d.revf(axis='F1')
    if has_param(headertree, 'STATE'):
        modify_val(headertree, 'STATE', ["1", "1", "0", "0"])
    else: 
        add_state(headertree, ["1", "1", "0", "0"])
        
    drepr = {1:'COMPLEX',0:'REAL'}
    data_representation = [drepr[1], drepr[d.axis1.itype], drepr[0], drepr[0]]
        
    if has_param(headertree, 'DATA_REPRESENTATION'):
        modify_val(headertree, 'DATA_REPRESENTATION', data_representation)
    else: 
        add_data_representation(headertree, data_representation)
 
    headerfile = op.join( dirname, "header.xml")
    headertree.write(headerfile)                   # writes in the headerfile
    if debug>0: headertree.write('spinit2D_header.xml')


if __name__ == '__main__':
    unittest.main()
