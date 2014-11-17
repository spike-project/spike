import numpy as np
import re
from matplotlib import pyplot as plt
from time import time

def read_dat(f):
    """
    reads Orbitrap .dat FID files
    returns a numpy array with the FID
    """
    F = open(f,'rb')
    allfile = F.read()
    linelength = re.findall("Data Points:.*$", allfile, re.MULTILINE)[0] #,
    lenghtfid = re.findall(r'\d+', linelength)[0]
    print lenghtfid
    offset = allfile.find("Data:\n")+6

    print "offset ",offset
    F.close()
    F = open(f,'rb')
    F.seek(offset)
    buf = np.fromfile(F, dtype = 'f4')
    return np.array(buf)

if __name__ == '__main__':
    addr = '/media/sdb/Orbitrap/increasing_microscans/20130523_ubi_5_scan_res_100000_10micro_1.dat'
    t0 = time()
    fid = read_dat(addr)
    print "time is ",(time()-t0)*1000, " ms"
    plt.plot(fid)
    plt.show()