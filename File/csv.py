#!/usr/bin/env python 
# encoding: utf-8

"""
    Utility to import and export data in text and csv files

    all functions compress transparently if the filenales end with .gz
    Marc-André adapted from some Lionel code
"""

from __future__ import print_function

__author__ = "Marc André Delsuc"
__date__ = "april 2014"

import numpy as np
import unittest
import os
import gzip

def save(data,filename, delimiter=','):
    "save 1D data in txt, single column, no unit - with attributes as pseudo comments "
    if filename.endswith('.gz'):
        do_open = gzip.open
    else:
        do_open = open
    with do_open(filename, 'w') as F:
        F.write("# %s \n"%data.axis1.report())  # first description
        for att in data.axis1.attributes:       # then attributes
            F.write("#$%s%s%s\n"%(att, delimiter, getattr(data.axis1,att)))
        np.savetxt(F, data.get_buffer(), delimiter=delimiter)

def load(filename, column=0, delimiter=','):
    """
    load 1D data from txt or csv file,
    attribute are in pseuo-coments startin with #$
    value are in columns, separated by delimiter - only the columun given in arg will be loaded
    column = 0 is fine for text files
    column = 1 is fine for csv files with currentunit
    returns a numpy buffer and an attribute dictionary
    """
    buf = []
    att = {}
    if filename.endswith('.gz'):
        do_open = gzip.open
    else:
        do_open = open
    with do_open(filename, 'r') as F:
        for l in F:
            fields = l.split(delimiter)
            if l.startswith('#'):    # first comments
                if l.startswith('#$'):  #  #$key value  for parameters
                    k = fields[0][2:]
                    v = fields[1].rstrip()
                    try:
                        v = float(v)        # eventually turn it into number
                    except ValueError:
                        pass
                    att[k] = v
                    print(k, v)
            else:
                buf.append( float( fields[column] ) )
    return np.array(buf), att

def save_unit(data, filename, delimiter=','):
    """save 1D data in csv,
    in 2 columns, with attributes as pseudo comments
    """
    ax = data.axis1.unit_axis()
    y = data.get_buffer()
    if filename.endswith('.gz'):
        do_open = gzip.open
    else:
        do_open = open
    with do_open(filename, 'w') as F:
        F.write("# %s \n"%data.axis1.report())
        for att in data.axis1.attributes:       # then attributes
            F.write("#$%s%s%s\n"%(att, delimiter, getattr(data.axis1,att)))
        np.savetxt(F, np.c_[ax,y], delimiter=delimiter)

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
    if "Orbitrap" in att.keys():
        d = OrbiData(buffer=buf)
    for k,v in att.items():
        if k in d.axis1.attributes:
            setattr(d.axis1, k, v)
        else:
            print("Warning - wrong attributes : ",k,v)
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
        r.save_csv(filename("test2.csv.gz"))
        rr = Import_1D(filename("test2.csv.gz"),column=1)
        self.assertAlmostEqual((r-rr).get_buffer().max(), 0.0)
        self.assertEqual(r.axis1.currentunit, rr.axis1.currentunit)
        os.unlink(filename("test2.csv.gz"))
if __name__ == '__main__':
    unittest.main()
