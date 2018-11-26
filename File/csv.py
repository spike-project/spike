#!/usr/bin/env python 
# encoding: utf-8

"""
    Utility to import and export data in text and csv files

    all functions compress transparently if the filenames end with .gz
    Marc-André adapted from some Lionel code
"""

from __future__ import print_function

__author__ = "Marc André Delsuc"
__date__ = "april 2014"

import numpy as np
from numpy.compat import asbytes, asstr
import unittest
import os
import gzip

def do_open(filename,flag):
    "opens regular and gz files"
    if filename.endswith('.gz'):
        return gzip.open(filename,flag)
    else:
        return open(filename,flag)

def save(data,filename, delimiter=',', fmt='%.18E'):
    "save 1D data in txt, single column, no unit - with attributes as pseudo comments "
    with do_open(filename, 'wb') as F:
        F.write(asbytes("# %s \n"%data.axis1.report()))  # first description
        for att in data.axis1.attributes:       # then attributes
            F.write(asbytes("#$%s%s%s\n"%(att, delimiter, getattr(data.axis1,att))))
        np.savetxt(F, data.get_buffer(), delimiter=delimiter, fmt=fmt)

def load(filename, column=0, delimiter=','):
    """
    load 1D data from txt or csv file,
    attribute are in pseudo-comments starting with #$
    value are in columns, separated by delimiter - only the column given in arg will be loaded
    column = 0 is fine for text files
    column = 1 is fine for csv files with currentunit
    returns a numpy buffer and an attribute dictionary
    """
    buf = []
    att = {}
    with do_open(filename, 'r') as F:
        for ll in F:
            l = asbytes(ll)
            fields = l.split(asbytes(delimiter))
            if l.startswith(asbytes('#')):    # first comments
                if l.startswith(asbytes('#$')):  #  #$key value  for parameters
                    k = asstr(fields[0][2:])
                    v = asstr(fields[1]).rstrip()
                    try:
                        v = float(v)        # eventually turn it into number
                    except ValueError:
                        pass
                    att[k] = v
                    print(k, v)
            else:
                buf.append( float( fields[column] ) )
    return np.array(buf), att

def save_unit(data, filename, delimiter=',', fmt='%.18E'):
    """save 1D data in csv,
    in 2 columns, with attributes as pseudo comments
    """
    step = data.axis1.itype+1
    ax = data.axis1.unit_axis()[::step]
    y = data.get_buffer()
    with do_open(filename, 'wb') as F:
        F.write(asbytes("# %s \n"%data.axis1.report()))
        for att in data.axis1.attributes:       # then attributes
            F.write(asbytes("#$%s%s%s\n"%(att, delimiter, getattr(data.axis1,att))))
        np.savetxt(F, np.c_[ax,y], fmt=fmt, delimiter=delimiter)
    
def Import_1D(filename, column=0, delimiter=','):
    """
    import a 1D file stored as csv
    header as comments (#)
    parameters in pseudocomments :
        #$key value
    then one value per line
    column and delimiter  as in load()
    """
    from ..NPKData import NPKData
    from ..FTICR import FTICRData
    from ..Orbitrap import OrbiData
    buf, att = load(filename, column=column, delimiter=delimiter)
    if "NMR" in att.keys():
        d = NPKData(buffer=buf)
    elif "FTICR" in att.keys():
        d = FTICRData(buffer=buf)
    elif "Orbitrap" in att.keys():
        d = OrbiData(buffer=buf)
    for kk,v in att.items():
        k = asstr(kk)
        if k in d.axis1.attributes:
            setattr(d.axis1, k, v)
        else:
            print("Warning - wrong attribute: ",k,v)
    d.adapt_size()
    return d


class csvTests(unittest.TestCase):
    """ - Testing NPKData basic behaviour - """
    def test_csv(self):
        from ..Tests import filename, directory
        from ..FTICR import FTICRData
        from ..NPKData import NPKData
        name2D = filename("dosy-cluster2.gs2")
        A = FTICRData(buffer=np.zeros(10000))
        A.specwidth = 1667000
        A.ref_mass = 344.0974
        A.ref_freq = 419620.0
        A.highmass = 1000.0
        A.currentunit="m/z"
        A.save_txt(filename("1D_test.txt"))
        B = FTICRData(dim=2)
        B.load_txt(filename("1D_test.txt"))
        self.assertAlmostEqual((A-B).get_buffer().max(), 0.0)
        os.unlink(filename("1D_test.txt"))
        
        d = NPKData(name=name2D)
        r = d.row(133)
        r.currentunit = "ppm"
        r.save_csv(filename("test2.csv.gz"), fmt="%.14g")
        rr = Import_1D(filename("test2.csv.gz"),column=1)
        self.assertAlmostEqual((r-rr).get_buffer().max(), 0.0)
        self.assertEqual(r.axis1.currentunit, rr.axis1.currentunit)
        os.unlink(filename("test2.csv.gz"))
if __name__ == '__main__':
    unittest.main()
