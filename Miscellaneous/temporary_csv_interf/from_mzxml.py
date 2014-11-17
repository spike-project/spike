#!/usr/bin/env python

'''
Read scans to a mzXML file
Based on PyMsXML (http://edwardslab.bmcb.georgetown.edu/software/PyMsXML.html)
'''
#########################################################################################
#                                                                                       #
# iQmerge software tool for merging CID and HCD scans for iTRAQ/TMT type experiments    #                                                                            #
# Written by Ashoka D. Polpitiya                                                        #
# for the Translational Genomics Research Institute (TGen), Phoenix, AZ                 #
#                                                                                       #
# Copyright 2011, Translational Genomics Research Institute                             #
# E-mail: ashoka@tgen.org                                                               #
# Website: http://iqmerge.googlecode.com                                                #
# --------------------------------------------------------------------------------------#
#                                                                                       #
# Licensed under the Apache License, Version 2.0 (the "License");                       #
# you may not use this file except in compliance with the License.                      #
# You may obtain a copy of the License at                                               #
#                                                                                       #
#       http://www.apache.org/licenses/LICENSE-2.0                                      #
#                                                                                       #
# Unless required by applicable law or agreed to in writing, software                   #
# distributed under the License is distributed on an "AS IS" BASIS,                     #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.              #
# See the License for the specific language governing permissions and                   #
# limitations under the License.                                                        #
#                                                                                       #
#########################################################################################

import sys,os, os.path
import gzip
from xml.sax import *
from mzXML.mzXMLmdsax import mzXMLmdsax
from mzXML.mzXMLspecsax import mzXMLspecsax
from mzXML.ToMzXML import ToMzXML

class FromMzXML:
    def __init__(self,filename):
        self.filename = filename
        self.parser = None

        if not os.path.exists(self.filename):
            print >>sys.stderr, "Filename %s does not exist."%(self.filename,)
            sys.exit(1)

        if self.filename.endswith('.gz'):
            self.compressfmt = 'gz'
        else:
            self.compressfmt = None
        print >>sys.stderr, "Initialized reading Filename %s ."%(self.filename,)
        sys.stderr.flush()

    def __del__(self):
        self.close()

    def open(self):

        if self.parser is None:
            # print >>sys.stderr, "Parsing out meta-data..."
            # sys.stderr.flush()
            self.parser = make_parser()
            self.handler = mzXMLmdsax()
            self.parser.setContentHandler(self.handler)
            try:
                if self.compressfmt is None:
                    self.parser.parse(self.filename)
                elif self.compressfmt == 'gz':
                    self.parser.parse(gzip.open(self.filename,'r'))
                else:
                    print >>sys.stderr, "Bad compressfmt specification: %s"%(self.compressfmt,)
                    sys.exit(1)
            except SAXException:
                pass

            self.md = self.handler.md

            # print >>sys.stderr, self.md
            # sys.stderr.flush()

            # print >>sys.stderr, "Set up spectrum parser..."
            # sys.stderr.flush()
            self.parser = make_parser()
            self.handler = mzXMLspecsax()
            self.parser.setContentHandler(self.handler)

    def close(self):
        pass
    
    def spectra(self):
        self.open()
        while True:
            # print >>sys.stderr, "Spectrum parser to read spectra from",self.handler.scanstart
            # sys.stderr.flush()
            try:
                if self.compressfmt is None:
                    self.parser.parse(self.filename)
                elif self.compressfmt == 'gz':
                    self.parser.parse(gzip.open(self.filename,'r'))
                else:
                    print >>sys.stderr, "Bad compressfmt specification: %s"%(self.compressfmt,)
                    sys.exit(1)
            except SAXException:
                pass
            for (s,d) in self.handler.spectra:
                # print >>sys.stderr, "Yeild spectrum number:",d['num']
                yield (s,d)
            if self.handler.done:
                break

    def getMsRunMetaData(self):
        self.open()
        d = {}
        if self.md.has_key('startTime'):
            d.update({'startTime:PT%fS':self.md['startTime']})
        if self.md.has_key('endTime'):
            d.update({'endTime:PT%fS':self.md['endTime']})
        return d

    def getFilenames(self):
        return self.md['parentFile']

    def msInstrumentID(self):
        return "1"

    def peakDetect(self,level):
        return False

    def msManufacturer(self):
        return self.md.get('msManufacturer','')

    def msModel(self):
        return self.md.get('msModel','')

    def msIonisation(self):
        return self.md.get('msIonisation','')

    def msMassAnalyzer(self):
        return self.md.get('msMassAnalyzer','')

    def msDetector(self):
        return self.md.get('msDetector','')

    def acquisitionSoftware(self):
        return self.md.get('acquisitionSoftware','')

    def acquisitionSoftwareVersion(self):
        return self.md.get('acquisitionSoftwareVersion','')

    def precursorPrecision(self,mz):
        return 0

    def precursorTune(self):
        return False

    def lc(self):
        return True

    def maldi(self):
        return not self.lc()
    
if __name__ == '__main__':
    #filename = 'C:/Ashoka/Research/IQmerge/Data/CAR_testslice.mzXML'
    #filename = "D:/Research/iQmerge/Data/CAR_testslice.mzXML"
    #r = FromMzXML(filename)
    #meta = r.getMsRunMetaData()
    #x = ToMzXML(r,"C:/Ashoka/Research/IQmerge/Data/testfile2.xml")
    #x = ToMzXML(r,"D:/Research/iQmerge/Data/tstFile.xml")
    #x.write()
    #print "Done..."
    pass 