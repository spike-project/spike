#
"""
Code to handle 2rr submatrix Bruker file format

Original code from L.Chiron
adapted my M-A Delsuc
"""
from __future__ import print_function
import numpy as np
from time import time
import os.path as op
from . import BrukerNMR as Bruker

debug = False
class BrukerSMXHandler(object):
    '''
    Reads/Writes Bruker 2rr files
    
    Writing not fully implemented

    '''    
    def __init__(self, addrproc):
        """
        Creates an object that hold everything
        addrproc the address of data the directory : eg   /here/expname/12/pdata/1 
        """
        self.addrproc = addrproc
        self.addr_expno = op.dirname(op.dirname(addrproc))
        #####
        self.acqu = Bruker.read_param( Bruker.find_acqu(self.addr_expno) )
        self.acqu2 = Bruker.read_param( Bruker.find_acqu2(self.addr_expno) )
        self.proc = Bruker.read_param( Bruker.find_proc(self.addrproc, down=False) )
        self.proc2 = Bruker.read_param( Bruker.find_proc2(self.addrproc, down=False) )
        ####
        self.si1 = int(self.proc2['$SI'])
        self.si2 = int(self.proc['$SI'])
        self.xd1 = int(self.proc2['$XDIM'])
        self.xd2 = int(self.proc['$XDIM'])
        if debug:
            for att in ("addrproc", "addr_expno", 'si1', "si2", "xd1", "xd2"):
                print("%s : %s"%(att,getattr(self,att)) )
    def load(self):
        '''
        Reads the data file (1D or 2D) and keeps it in self.data_1d for 1D
        in self.data_2d_2rr and self.data_2d_2ri for 2D.
        '''
        if os.path.exists(op.join(self.addr_expno,'fid')):
            self.read_1D()
        elif os.path.exists(op.join(self.addr_expno,'ser')):
            self.read_2D()

    def read_2D(self):
        '''
        Reads the 2D "fid" file and keeps it in self.data
        data taken as integers. 
        '''
        if debug:   print ("read 2D")
        self.data_DIM = '2d'

        fn = op.join(self.addrproc,'2rr')
        self.data_2d_2rr = self.reorder_subm(np.fromfile(fn, 'i4'))

        fn = op.join(self.addrproc,'2ri')
        if op.exists(fn):
            if debug:   print ("read 2ri")
            self.data_2d_2ri = self.reorder_subm(np.fromfile(fn, 'i4'))
        else:
            self.data_2d_2ri = None

        fn = op.join(self.addrproc,'2ir')
        if op.exists(fn):
            if debug:   print ("read 2ir")
            self.data_2d_2ir = self.reorder_subm(np.fromfile(fn, 'i4'))
        else:
            self.data_2d_2ir = None
        fn = op.join(self.addrproc,'2ri')
        if op.exists(fn):
            if debug:   print ("read 2ii")
            self.data_2d_2ii = self.reorder_subm(np.fromfile(fn, 'i4'))
        else:
            self.data_2d_2ii = None

        if debug:   print ("read finished")
                    
    def prepare_mat(self):
        '''
        sub_per_dim : [nb of submatrices in t1, nb of submatrices in t2]
        self.proc['$SI'] : dimension 2 of the 2D data. 
        self.proc2['$SI'] : dimension 1 of the 2D data. 
        self.acqu['$XDIM'] : size submatrix
        '''
        if debug:   print ("prepare matrix")
        self.rdata = np.empty((self.si1, self.si2)) #, dtype = 'complex128')
        if debug:
            print ("rdata.shape ",self.rdata.shape)
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

    def reorder_bck_subm(self,data):
        """
        Reorder flat matrix back to sbmx Bruker data.
        self.sub_per_dim : [nb of submatrices in t1, nb of submatrices in t2]
        self.nsubs : total number of submatrices
        self.proc['$SI'] : shape of the 2D data. 
        self.acqu['$XDIM'] : size submatrix
        """
        if debug:   print ("reorder matrix back")
        self.prepare_mat()
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
        self.prepare_mat()
        #longmat = int(self.proc['$SI']), int(self.proc2['$XDIM'])*self.sub_per_dim[1]
        # print "data.shape ",data.shape
        # print "longmat ",self.long_mat
        interm = data.reshape(self.long_mat)      
        mat = []
        for sub_num, sub_idx in enumerate(np.ndindex(tuple(self.sub_per_dim))):
            zipshape = zip(sub_idx, self.dim_sub_mat)
            sub_slices = [slice(i * j, (i + 1) * j) for i, j in zipshape ]
            slt2 = slice(sub_num*self.dim_sub_mat[0],(sub_num+1)*self.dim_sub_mat[0]) # dimension t1
            slt1 = slice(0,self.dim_sub_mat[1])# dimension t2
            # print "self.rdata[sub_slices].shape ",self.rdata[sub_slices].shape
            # print "interm[slt1, slt2].shape ",interm[slt1, slt2].shape
            # print "self.rdata[sub_slices].shape",self.rdata[sub_slices].shape
            self.rdata[sub_slices] = interm[slt2, slt1]
        data = self.rdata
        return data

    def write_file(self, data, filename):
        '''
        data written as integers. 
        '''
        f = open(filename, 'wb')
        if self.acqu['$BYTORDA'] == '0': 
            f.write(data.astype('<i4').tostring()) # little endian
        else:
            f.write(data.astype('>i4').tostring()) # big endian
        f.close()
        print ("rewrote ", filename)



