#!/usr/bin/env python

'''
Created on May 5, 2011

@author: apolpitiya
'''

import sys
import sha
import tempfile
import gzip
import zlib
import os
from base64 import b64encode
from urllib import quote
from xml.sax import *

class NameVersion:
    def __init__(self):
        self.name_ = "iQmerge"
        self.majorVer = 0
        self.minorVer = 5
        self.revVer = 4
    def version(self):
        return "%d.%d.%d"%(self.majorVer,self.minorVer,self.revVer)
    def name(self):
        return self.name_

class ToMzXML:
    def __init__(self,reader,filename,cmpfmt=None,filt=None,peaksCompress=False,version="3.0"):
        self.filename = filename
        self.compressfmt = cmpfmt
        self.reader = reader
        self.initSHA()
        self.writeByteCounter = 0
        self.actualScanCount = 0
        self.scanIndices = {}
        self.filter = filt
        if not version in ('2.1','2.2','3.0'):
            raise "Bad mzXML version \"%s\""%(version,)
        self.version = version
        if peaksCompress and float(self.version) < 3.0:
            raise "Cannot compress peaks until mzXML version 3.0"
        self.peakcompress = peaksCompress

    def __del__(self):
        if hasattr(self,'tmpFileName') and self.tmpFileName and self.tmpFileName != '':
            os.unlink(self.tmpFileName)

    def getFilename(self):
        return self.filename
   
    def initSHA(self):
        self.sha1 = sha.new()

    def updateSHA(self,s):
        self.sha1.update(s)

    def getSHA(self):
        return self.sha1.hexdigest().lower()

    def writeString(self,fh,data):
        self.writeByteCounter += len(data)
        self.sha1.update(data)
        fh.write(data)

    def write(self,debug=False):

        # Make a temporary file for the scan data to count the number of spectra.
        (tmpFileFD,self.tmpFileName) = tempfile.mkstemp(dir='.',prefix='.pymsxml')
        tmpFile = os.fdopen(tmpFileFD,'wb')

        # Has a number of side effects, it fills in scanIndices, and sets self.actualScanCount.
        self.write_scans(tmpFile,debug)

        tmpFile.close()

        # Reset!

        self.initSHA()
        self.writeByteCounter = 0

        if self.compressfmt is None and not self.filename.endswith('.gz'):
            xmlFile = open(self.filename,'wb')
        elif self.compressfmt == 'gz' or self.filename.endswith('.gz'):
            if self.filename.endswith('.gz'):
                xmlFile = gzip.open(self.filename,'wb')
            else:
                xmlFile = gzip.open(self.filename+'.gz','wb')
        else:
            print >>sys.stderr, "Bad compressfmt specification"
            sys.exit(1)

        self.writeString (xmlFile,'<?xml version="1.0" encoding="ISO-8859-1"?>\n')
 
        outStr = "<mzXML "
        outStr += "xmlns=\"http://sashimi.sourceforge.net/schema_revision/mzXML_%s\" "%(self.version,)
        outStr += "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
        outStr += "xsi:schemaLocation=\"http://sashimi.sourceforge.net/schema_revision/mzXML_%s http://sashimi.sourceforge.net/schema_revision/mzXML_%s/mzXML_idx_%s.xsd\""%(self.version,self.version,self.version)
        outStr += ">\n"

        self.writeString (xmlFile,outStr)
 
        outStr = "<msRun"
        outStr += " scanCount=\"%d\"" % (self.actualScanCount,)
        md = self.reader.getMsRunMetaData()
        for (k,v) in md.items():
            k,f = k.split(':')
            if self.version not in ('2.2',) or (k != 'startTime'and k != 'endTime'):
                outStr += (" %%s=\"%s\""%(f,))%(k,v)
        outStr += ">\n"

        self.writeString (xmlFile,outStr)

        filehashes = {}
       
        for f in self.reader.getFilenames():
            if len(f) == 2:
                outStr = "<parentFile "
                outStr += "fileName=\"%s\" " % (quote(f[0]),)
                outStr += "fileType=\"RAWData\" "
                outStr += "fileSha1=\"%s\"" % (f[1],)
                outStr += "/>\n"
            else:
                outStr = "<parentFile "
                outStr += "fileName=\"%s\" " % (f[0],)
                outStr += "fileType=\"%s\" " % (f[1],)
                outStr += "fileSha1=\"%s\""  % (f[2],)
                outStr += "/>\n"                
            self.writeString (xmlFile,outStr)

        if float(self.version) >= 3.0:
            outStr = "<msInstrument msInstrumentID=\"%s\">\n" % (self.reader.msInstrumentID(),)
        else:
            outStr = "<msInstrument>\n"
           
        outStr += "<msManufacturer category=\"msManufacturer\" value=\"%s\"/>\n" % (self.reader.msManufacturer(),)
        outStr += "<msModel category=\"msModel\" value=\"%s\"/>\n" % (self.reader.msModel(),)
        if self.reader.msIonisation():
            outStr += "<msIonisation category=\"msIonisation\" value=\"%s\"/>\n" % (self.reader.msIonisation(),)
        if self.reader.msMassAnalyzer():
            outStr += "<msMassAnalyzer category=\"msMassAnalyzer\" value=\"%s\"/>\n" % (self.reader.msMassAnalyzer(),)
        if self.reader.msDetector():
            outStr += "<msDetector category=\"msDetector\" value=\"%s\"/>\n" % (self.reader.msDetector(),)

        self.writeString (xmlFile,outStr)
 
        outStr = "<software "
        outStr += "type=\"acquisition\" "
        outStr += "name=\"%s\" "%(self.reader.acquisitionSoftware(),)
        outStr += "version=\"%s\""%(self.reader.acquisitionSoftwareVersion(),)
        outStr += "/>\n"

        self.writeString (xmlFile,outStr)
       
        outStr = "</msInstrument>\n"

        self.writeString (xmlFile,outStr)

        outStr = "<dataProcessing>\n"

        self.writeString (xmlFile,outStr)

        nv = NameVersion()

        outStr = "<software "
        outStr += "type=\"conversion\" "
        outStr += "name=\"%s\" "%(nv.name(),)
        outStr += "version=\"%s\"" % (nv.version(),)
        outStr += "/>\n"

        self.writeString (xmlFile,outStr)

        if (self.reader.peakDetect(1)):
            outStr = "<processingOperation "
            outStr += "name=\"peak_detection_level_1\" "
            outStr += "value=\"true\""
            outStr += "/>\n"
            outStr += "<processingOperation "
            outStr += "name=\"peak_detection_level_1_threshold\" "
            outStr += "value=\"1%\""
            outStr += "/>\n"
            outStr += "<processingOperation "
            outStr += "name=\"peak_detection_level_1_software\" "
            outStr += "value=\"%s\""%(self.reader.acquisitionSoftware(),)
            outStr += "/>\n"
            self.writeString (xmlFile,outStr)
        if (self.reader.peakDetect(2)):
            outStr = "<processingOperation "
            outStr += "name=\"peak_detection_level_2\" "
            outStr += "value=\"true\""
            outStr += "/>\n"
            outStr += "<processingOperation "
            outStr += "name=\"peak_detection_level_2_threshold\" "
            outStr += "value=\"1%\""
            outStr += "/>\n"
            outStr += "<processingOperation "
            outStr += "name=\"peak_detection_level_2_software\" "
            outStr += "value=\"%s\""%(self.reader.acquisitionSoftware(),)
            outStr += "/>\n"
            self.writeString (xmlFile,outStr)

        outStr = "<processingOperation "
        outStr += "name=\"min_peaks_per_spectra\" "
        outStr += "value=\"1\""
        outStr += "/>\n"

        self.writeString (xmlFile,outStr)

        outStr = "</dataProcessing>\n"

        self.writeString (xmlFile,outStr)

        if self.reader.maldi():
            outStr = '<spotting>\n'
            for d in self.reader.plateData():
                outStr += '<plate plateID="%s" spotXCount="%s" spotYCount="%s">\n' % \
                          (d['plateID'],d['spotXCount'],d['spotYCount'])
                outStr += '<plateManufacturer category="plateManufacturer" value="%s"/>\n' % \
                          (d['plateManufacturer'],)
                outStr += '<plateModel category="plateModel" value="%s"/>\n' % \
                          (d['plateModel'],)
                for s in self.reader.spotData(d['plateID']):
                    outStr += '<spot spotID="%s" spotXPosition="%s" spotYPosition="%s">\n' % \
                              (s['spotID'],s['spotXPosition'],s['spotYPosition'])
                    outStr += '<maldiMatrix category="maldiMatrix" value="%s"/>\n' % \
                              (s['maldiMatrix'],)
                    outStr += '</spot>\n'
                outStr += '</plate>\n'
            outStr += '</spotting>\n'
            self.writeString(xmlFile,outStr)

        if self.reader.lc():
            pass
       
        scanOffset = self.writeByteCounter

        #Start writing scans now. They are in the temp file
        tmpFile = open(self.tmpFileName,'rb')
        while True:
            tmpStr = tmpFile.read(1024)
            if tmpStr == '':
                break
            self.writeString(xmlFile,tmpStr)
        tmpFile.close()
        os.unlink(self.tmpFileName)
        self.tmpFileName = ''

        outStr = "</msRun>\n"
        self.writeString (xmlFile,outStr)

        indexOffset = self.writeByteCounter
                                                 
        outStr = "<index "
        outStr += "name=\"scan\" "
        outStr += ">\n"

        self.writeString (xmlFile,outStr)

        for i in xrange(1,self.actualScanCount+1):
            outStr = "<offset id=\"%d\">" % (i,)
            outStr += "%d</offset>\n" % (self.scanIndices[i] + scanOffset,)
            self.writeString(xmlFile,outStr)

        outStr = "</index>\n"
        self.writeString (xmlFile,outStr)

        outStr = "<indexOffset>%d</indexOffset>\n" % (indexOffset,)

        self.writeString (xmlFile,outStr)

        self.writeString (xmlFile,"<sha1>")

        outStr = self.getSHA()

        self.writeString (xmlFile,outStr)
       
        self.writeString (xmlFile,"</sha1>\n")
 
        outStr = "</mzXML>\n"
        self.writeString (xmlFile,outStr)

        xmlFile.close()

    def write_scans(self,xmlFile,debug=False):
        '''
        scans are read and written to the tempfile (xml)
        '''
        msLevel = 0
        scanNumber = 0
        ancestors = []
        self.writeByteCounter = 0
       
        for (s,d) in self.reader.spectra():

            if self.filter is not None and not self.filter.test(d):
                continue

            if debug and scanNumber >= 10:
                break
           
            if not d.has_key('msLevel'):
                print >>sys.stderr, "Required scan attributes missing."
                sys.exit(1)

            prevLevel = msLevel
            msLevel = d['msLevel']

            if prevLevel < msLevel and prevLevel > 0:
                # We went "in" a level, push scan number of parent
                ancestors.append((scanNumber,prevSpec))
            elif prevLevel > msLevel and msLevel > 0:
                if len(ancestors) == 0:
                    pass #print >>sys.stderr, "No ancestor for scan %s at level %s"%(scanNumber,msLevel)
                else:
                    ancestors.pop()
           
            outStr = ''
            if prevLevel > msLevel:
                for m in xrange(0,prevLevel-msLevel+1):
                    outStr += "</scan>\n"
            else:
                if prevLevel > 0 and prevLevel == msLevel:
                    outStr += "</scan>\n"

            self.writeString(xmlFile,outStr)

            scanNumber = scanNumber + 1

            self.scanIndices[scanNumber] = self.writeByteCounter

            totIonCurrent = 0
            maxmz = None
            minmz = None
            basePeakMz = None
            basePeakIntensity = None
            peaksCount = 0
            for i in xrange(0,len(s),2):
                x = s[i]; y = s[i+1]
                if minmz is None or minmz > x:
                    minmz = x
                if maxmz is None or maxmz < x:
                    maxmz = x
                totIonCurrent += y
                if basePeakIntensity is None or basePeakIntensity < y:
                    basePeakIntensity = y
                    basePeakMz = x
                peaksCount += 1

            outStr = "<scan num=\"%d\" msLevel=\"%d\"" % (scanNumber,msLevel)
            outStr += " peaksCount=\"%d\"" % (peaksCount,)
            outStr += " lowMz=\"%f\"" % (minmz,)
            outStr += " highMz=\"%f\"" % (maxmz,)
            outStr += " totIonCurrent=\"%f\"" % (totIonCurrent,)
            outStr += " basePeakMz=\"%f\"" % (basePeakMz,)
            outStr += " basePeakIntensity=\"%f\"" % (basePeakIntensity,)

            for (k,v) in d.items():
                if k.startswith('scan.'):
                    k,f = k.split(':')
                    k = k[5:]
                    outStr += (" %%s=\"%s\""%(f,))%(k,v)
            outStr += ">\n"
            self.writeString(xmlFile,outStr)

            if self.version in ('2.2',):
                any = False
                for (k,v) in d.items():
                    if k.startswith('scanOrigin.'):
                        if not any:
                            any = True
                            outStr = "<scanOrigin"
                        k,f = k.split(':')
                        k = k[11:]
                        outStr += (" %%s=\"%s\""%(f,))%(k,v)
                if any:
                    outStr += "/>\n"
                    self.writeString(xmlFile,outStr)

            if msLevel > 1:

                if not d.has_key('precursorMz') :
                    print >>sys.stderr, "Required precursorMz attribute missing."
                    sys.exit(1)

                # We scan the parent spectrum in the region of the
                # precursorMz to establish the intensity. Optionally,
                # we tune the precursorMz itself, on the basis of this

                pMz = d['precursorMz']
                tol = self.reader.precursorPrecision(pMz)
                pMzLB = pMz-tol
                pMzUB = pMz+tol
                pMzIntensity = 0

                if len(ancestors) > 0:
                    for i in xrange(0,len(ancestors[-1][1]),2):
                        x = ancestors[-1][1][i]; y = ancestors[-1][1][i+1];
                        if x < pMzLB:
                            continue
                        if x > pMzUB:
                            break
                        if pMzIntensity < y:
                            pMzIntensity = y
                            pMzTune = x

                    if self.reader.precursorTune():
                        pMz = pMzTune

                outStr = "<precursorMz"
                for (k,v) in d.items():
                    if k.startswith('precursorMz.'):
                        k,f = k.split(':')
                        k = k[12:]
                        outStr += (" %%s=\"%s\""%(f,))%(k,v)
                if len(ancestors)>0:
                    outStr += " precursorScanNum=\"%d\""%(ancestors[-1][0],)
                    outStr += " precursorIntensity=\"%f\">"%(pMzIntensity,)
                else:
                    outStr += ">"
                outStr += "%f</precursorMz>\n" % (pMz,)
                self.writeString(xmlFile,outStr)

            if self.reader.maldi() and self.version in ('2.2',):
                outStr = '<maldi'
                for (k,v) in d.items():
                    if k.startswith('maldi.'):
                        k,f = k.split(':')
                        k = k[6:]
                        outStr += (" %%s=\"%s\""%(f,))%(k,v)
                outStr += " plateID=\"%s\""%(d['plateID'],)
                outStr += " spotID=\"%s\""%(d['spotID'],)
                outStr += "/>\n"
                self.writeString(xmlFile,outStr)

            if sys.byteorder != 'big':
                s.byteswap()

            if debug:
                s = s[:20]

            specstr = s.tostring()
            if self.peakcompress:
                specstr = zlib.compress(specstr,9)
                outStr = "<peaks precision=\"32\" byteOrder=\"network\" contentType=\"m/z-int\" compressionType=\"zlib\" compressedLen=\"%d\">"%(len(specstr),)
            else:
                if float(self.version) >= 3.0:
                    outStr = "<peaks precision=\"32\" byteOrder=\"network\" contentType=\"m/z-int\" compressionType=\"none\" compressedLen=\"%d\">"%(len(specstr),)
                else:
                    outStr = "<peaks precision=\"32\" byteOrder=\"network\" pairOrder=\"m/z-int\">"
                   
            self.writeString(xmlFile,outStr)

            # Spec says the byte order shall be
            # network, which is the same as big endian.

            if sys.byteorder != 'big':
                s.byteswap()

            outStr = b64encode(specstr)
            self.writeString(xmlFile,outStr)

            outStr = "</peaks>\n"
            self.writeString(xmlFile,outStr)

            if self.reader.maldi() and self.version not in ('2.2',):
                outStr = ''
                for (k,v) in d.items():
                    if k.startswith('maldi.'):
                        k,f = k.split(':')
                        k = k[6:]
                        outStr += '<nameValue'
                        outStr += ' name="maldi.%s"'%(k,)
                        outStr += (' value="%s"'%(f,))%(v,)
                        outStr += '/>\n'
                outStr += '<nameValue'
                outStr += ' name="maldi.plateID"'
                outStr += ' value="%s"'%(d['plateID'],)
                outStr += '/>\n'
                outStr += '<nameValue'
                outStr += ' name="maldi.spotID"'
                outStr += ' value="%s"'%(d['spotID'],)
                outStr += '/>\n'
                self.writeString(xmlFile,outStr)

            if self.version not in ('2.2',):
                outStr = ''
                for (k,v) in d.items():
                    if k.startswith('scanOrigin.'):
                        k,f = k.split(':')
                        k = k[11:]
                        outStr += '<nameValue'
                        outStr += ' name="scanOrigin.%s"'%(k,)
                        outStr += (' value="%s"'%(f,))%(v,)
                        outStr += '/>\n'
                self.writeString(xmlFile,outStr)

            outStr = ''
            for (k,v) in d.items():
                if k.startswith('nameValue.'):
                    k,f = k.split(':')
                    k = k[10:]
                    outStr += '<nameValue'
                    outStr += ' name="%s"'%(k,)
                    outStr += (' value="%s"'%(f,))%(v,)
                    outStr += '/>\n'
            if len(outStr) > 0:
                self.writeString(xmlFile,outStr)

            xmlFile.flush()

            prevSpec = s

        outStr = ""
        for m in xrange(0,len(ancestors)+1):
            outStr += "</scan>\n"

        self.writeString(xmlFile,outStr)
        self.actualScanCount = scanNumber

