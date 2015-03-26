import numpy as np
import xml.dom.minidom
import base64
import struct


class mzXML():
    '''
    Draft for reqding mzXML files
    '''
    def __init__(self, filename):
        self.doc = xml.dom.minidom.parse(filename)
        self.get_params()
    
    def txt2np(self, txt):
        '''
        from mzXML to numpy
        '''
        data = base64.b64decode(txt)
        endian = '>'
        precision = 'f'
        count = len(data) / struct.calcsize(endian + precision)
        nump = np.array(struct.unpack(endian + precision * count, data[0:len(data)]))
        return nump
        
    def np2txt(self, nump):
        '''
        from numpy to mzXML
        '''
        endian = '>'
        precision = 'f'
        count = len(nump)
        binary = struct.pack(endian + precision * count, *nump)
        txt = base64.b64encode(binary)
        return txt

    def get_params(self):
        '''
        Extract parqmeters for limits and size
        '''
        spec = self.doc.getElementsByTagName('peaks')
        scan = self.doc.getElementsByTagName('scan')
        num = self.txt2np(spec[0].childNodes[0].toxml())
        print 'type(num) ', type(num)
        self.spec = num[0::2]
        self.mzaxis = num[1::2]
        self.scan_size = len(self.txt2np(spec[0].childNodes[0].toxml()))/2
        self.hmz = float(scan[0].getAttribute('highMz'))
        self.lmz = float(scan[0].getAttribute('lowMz'))
        
    def save_mzXML(self, x, y, namefile_final):
        '''
        Save the processed data in the new mzXML file.
        '''
        spec = self.doc.getElementsByTagName('peaks')
        txt = self.np2txt(np.array(zip(x,y)).ravel())
        spec[0].childNodes[0].replaceWholeText(txt) # Replacing in the mzXML
        self.doc.writexml(open(namefile_final,'w')) # Saving the mzXML
        print 'mzXML saved'
