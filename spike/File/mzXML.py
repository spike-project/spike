#!/usr/bin/env python 
# encoding: utf-8

from __future__ import print_function, division
import numpy as np
import xml.dom.minidom
import base64
import struct
import unittest

class mzXML():
    '''
    Draft for reading mzXML files
    '''
    def __init__(self, filename):
        self.doc = xml.dom.minidom.parse(filename)
        self.debug = 1
        self.get_params()
        
    def txt2np(self, txt):
        '''
        from mzXML to numpy
        Float 32, big endian.
        Returns numpy format.
        '''
        data = base64.b64decode(txt)
        endian = '>'
        precision = 'f'
        count = len(data) // struct.calcsize(endian + precision)
        nump = np.array(struct.unpack(endian + precision * count, data[0:len(data)]))
        return nump
        
    def np2txt(self, nump):
        '''
        from numpy to mzXML
        Float 32, big endian.
        Returns text format.
        '''
        endian = '>'
        precision = 'f'
        count = len(nump)
        binary = struct.pack(endian + precision * count, *nump)
        txt = base64.b64encode(binary)
        return txt
    
    def scattr(self, scan, param, kind=None):
        '''
        Scan attributes
        kind specifies if we want integer, float etc.. 
        '''
        param_value = scan[0].getAttribute(param)
        if kind:
            return kind(param_value)
        else:
            return param_value
    
    def get_params(self):
        '''
        Extract parameters for limits and size
        '''
        spec = self.doc.getElementsByTagName('peaks')
        scan = self.doc.getElementsByTagName('scan')
        num = self.txt2np(spec[0].childNodes[0].toxml()) # numpy data for spectrum with axis
        self.spec = num[0::2]
        self.mzaxis = num[1::2]
        self.scan_size = len(self.txt2np(spec[0].childNodes[0].toxml()))/2 # size of the scan
        self.dic_param = {'msLevel':int, 'peaksCount':int, 'polarity':None, 'scanType':None,
            'filterLine':None, 'retentionTime':None, 'lowMz':float, 'highMz':float,
             'basePeakMz':float, 'basePeakIntensity':float, 'totIonCurrent':int, 'collisionEnergy':int}
        for p in self.dic_param:
            if self.debug > 2:
                print("param ", p)
            try:
                setattr(self, p, self.scattr(scan, p, self.dic_param[p])) 
            except:
                print("parameter {} not found".format(p))
        
    def save_mzXML(self, x, y, namefile_final):
        '''
        Save data in the new mzXML file.
        '''
        spec = self.doc.getElementsByTagName('peaks')
        txt = self.np2txt(np.array(zip(x,y)).ravel())
        spec[0].childNodes[0].replaceWholeText(txt) # Replacing in the mzXML
        self.doc.writexml(open(namefile_final,'w')) # Saving the mzXML
        print('mzXML saved')

class mzXML_Tests(unittest.TestCase):  
    def test_simple(self):
        from ..Tests import filename, directory
        mzxml = mzXML(directory()+'/testsmzXML/20131129_V_EV_bC_1499_CID30_r7k_25usc_range410-2000.mzXML')
        dic_param = {'msLevel':2, 'peaksCount':35775, 'polarity':"+", 'scanType':"Full",
            'filterLine':"FTMS + p NSI Full ms2 1499.60@cid30.00 [410.00-2000.00]",
            'retentionTime':"PT0.0366S", 'lowMz':410.001, 'highMz':2069.29, 
            'basePeakMz':1151.69, 'basePeakIntensity':13452.2, 'totIonCurrent':681064,
            'collisionEnergy':30}
        for p in dic_param:
            getat = getattr(mzxml, p)
            self.assertEqual(getattr(mzxml,p), dic_param[p])
        
if __name__ == '__main__':
    unittest.main()