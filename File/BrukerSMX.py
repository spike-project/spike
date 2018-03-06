#
"""
Code to handle 2rr submatrix Bruker file format

used internally by BrukerNMR, not for final use

Original code from L.Chiron
adapted my M-A Delsuc
"""
from __future__ import print_function

import os.path as op
import numpy as np
from . import BrukerNMR as Bruker

debug = False
class BrukerSMXHandler(object):
    '''
    Reads/writes Bruker 2rr files

    Writing not fully implemented

    '''
    def __init__(self, addrproc):
        """
        Creates an object that holds everything
        addrproc the address of data the directory : eg   /here/expname/12/pdata/1
        """
        self.addrproc = op.abspath(addrproc)
        self.addr_expno = op.dirname(op.dirname(addrproc))
        #####
        self.acqu = Bruker.read_param( Bruker.find_acqu(self.addr_expno) )
        self.proc = Bruker.read_param( Bruker.find_proc(self.addrproc, down=False) )

        self.acqu2 = Bruker.read_param( Bruker.find_acqu2(self.addr_expno) )
        self.proc2 = Bruker.read_param( Bruker.find_proc2(self.addrproc, down=False) )
        self.si1 = int(self.proc2['$SI'])
        self.si2 = int(self.proc['$SI'])
        self.xd1 = int(self.proc2['$XDIM'])
        self.xd2 = int(self.proc['$XDIM'])

        if debug:
            for att in ("addrproc", "addr_expno", 'si1', "si2", "xd1", "xd2"):
                print("%s : %s"%(att,getattr(self,att)) )

    def read_smx(self):
        '''
        Reads the 2D "smx" (2rr 2ri 2ir 2ii) files and keeps it in self.data_2d_2xx
        data taken as integers.
        '''
        if debug:   print ("reading 2D")
        self.prepare_mat()
        for name in ('2rr','2ir','2ri','2ii'):

            fn = op.join(self.addrproc, name)
            dataname = "data_2d_"+name          # result stored in object as data_2d_2xx

            if op.exists(fn):
                if debug:   print ("read", name)
                mat =  self.read_file(fn)
                setattr(self, dataname, self.reorder_subm(mat) )
            else:
                setattr(self, dataname, None)
        if debug:
            print ("2rr.shape ", self.data_2d_2rr.shape)
            print ("read finished")

    def write_smx(self):
        '''
        writes the prepared self.data_2d_2xx (x in {r,i} ) into 2xx files

        should have been prepared elsewhere
        missing data are not written
        '''
        if debug:
            print ("writing 2D")
            print ("2rr.shape ", self.data_2d_2rr.shape)
        self.prepare_mat()
        for name in ('2rr','2ir','2ri','2ii'):
            try:
                dataset = getattr(self,  'data_2d_' + name)
            except:
                continue    # go to next if not found
            else:   # or do the stuff
                if dataset is not None:
                    fn = op.join(self.addrproc, name)
                    data = self.reorder_bck_subm(dataset)
                    self.write_file(data, fn)
                    if debug:   print (name, " written")
        if debug:   print ("read finished")

    def prepare_mat(self):
        '''
        sub_per_dim : [nb of submatrices in t1, nb of submatrices in t2]
        self.si2 == self.proc['$SI'] : dimension 2 of the 2D data.
        self.si1 == self.proc2['$SI'] : dimension 1 of the 2D data.
        self.xd2 == self.acqu['$XDIM'] : size submatrix in F2
        self.xd1 == self.acqu2['$XDIM'] : size submatrix in F1
        '''
        self.dim_mat = (self.si1, self.si2)
        self.dim_sub_mat = (self.xd1, self.xd2)
        zipshape = zip(self.dim_mat, self.dim_sub_mat)
        if debug:
            print ("self.dim_mat, self.dim_sub_mat ",self.dim_mat, self.dim_sub_mat )
            print("zipshape", zipshape)
        self.sub_per_dim = [int(i / j) for i, j in zipshape]
        if debug:
            print ("self.sub_per_dim ", self.sub_per_dim)
            print ("self.dim_mat[0] ",self.dim_mat[0])
        self.long_mat = self.dim_mat[0]*self.sub_per_dim[1], self.dim_sub_mat[1]

    def reorder_bck_subm(self, data):
        """
        Reorder flat matrix back to sbmx Bruker data.
        self.sub_per_dim : [nb of submatrices in t1, nb of submatrices in t2]
        self.nsubs : total number of submatrices
        self.proc['$SI'] : shape of the 2D data.
        self.acqu['$XDIM'] : size submatrix
        """
        if debug:   print ("reorder matrix back")
        interm = data.reshape(self.dim_mat)
        mat = []
        for sub_num, sub_idx in enumerate(np.ndindex(tuple(self.sub_per_dim))):
            zipshape = zip(sub_idx, self.dim_sub_mat)
            sub_slices = [slice(i * j, (i + 1) * j) for i, j in zipshape ]
            mat.append(list(np.ravel(interm[sub_slices])))
        data = np.array(mat).reshape(self.dim_mat)
        return data

    def reorder_subm(self, data):
        """
        Reorder sbmx binary Bruker data to flat matrix.
        self.sub_per_dim : [nb of submatrices in t1, nb of submatrices in t2]
        self.nsubs : total number of submatrices
        self.proc['$SI'] : shape of the 2D data.
        self.acqu['$XDIM'] : size submatrix
        """
        if debug:    print ("reorder matrix")
        rdata = np.empty((self.si1, self.si2)) #, dtype = 'complex128')
        interm = data.reshape(self.long_mat)
        for sub_num, sub_idx in enumerate(np.ndindex(tuple(self.sub_per_dim))):
            zipshape = zip(sub_idx, self.dim_sub_mat)
            sub_slices = [slice(i * j, (i + 1) * j) for i, j in zipshape ]
            slt2 = slice(sub_num*self.dim_sub_mat[0],(sub_num+1)*self.dim_sub_mat[0]) # dimension t1
            slt1 = slice(0,self.dim_sub_mat[1])# dimension t2
            # print "self.rdata[sub_slices].shape ",self.rdata[sub_slices].shape
            # print "interm[slt1, slt2].shape ",interm[slt1, slt2].shape
            # print "self.rdata[sub_slices].shape",self.rdata[sub_slices].shape
            rdata[sub_slices] = interm[slt2, slt1]
        return rdata

    def read_file(self, filename):
        '''
        data read as integers.
        '''
        if self.proc['$BYTORDP'] == '0':
            fmt = '<i4'       # little endian
        else:
            fmt = '>i4'      # big endian
        d = np.fromfile(filename, fmt)
        return d
    def write_file(self, data, filename):
        '''
        data written as integers.
        '''
        if self.proc['$BYTORDP'] == '0':
            fmt = '<i4'       # little endian
        else:
            fmt = '>i4'      # big endian
        with open(filename, 'wb') as f:
            f.write(data.astype(fmt).tostring())

