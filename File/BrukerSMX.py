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
        if op.exists( op.join(self.addrproc, '1r') ) or op.exists( op.join(self.addr_expno, 'fid') ):
            self.dim = 1
        elif op.exists( op.join(self.addrproc, '2rr') )  or op.exists( op.join(self.addr_expno, 'ser') ):
            self.dim = 2
        else:
            raise Exception('Cannot determine data dimension')
        if self.dim == 2:
            self.acqu2 = Bruker.read_param( Bruker.find_acqu2(self.addr_expno) )
            self.proc2 = Bruker.read_param( Bruker.find_proc2(self.addrproc, down=False) )
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

    def reorder_bck_subm(self, data):
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
        if self.proc['$BYTORDP'] == '0':
            f.write(data.astype('<i4').tostring()) # little endian
        else:
            f.write(data.astype('>i4').tostring()) # big endian
        f.close()
        print ("rewrote ", filename)

class BrukerWriter(BrukerSMXHandler):
    """
    Exports a 1D or a 2D npkdata to a  Bruker 1r / 2rr file, using another SMX as a template

    """
    def __init__(self, d, fname, templname):
        """initialize
        d is an NMR NPKData 

        fname and templname are procno :    eg   /here/expname/12/pdata/1 

        and the files are created in the fname directory
        a pdata/procno should already exists as a template

        if d contains metadata parameters from Bruker, there will be used,
        however all files common to fname and templname expno will not be updated

        if fname and templname are exactly the same, (or templname is None)
            the proc and procs files will be overwriten
        """
        self.d = d
        if d.dim>2:
            raise Exception('Not implemented yet')
        # we use templname as basis
        super(BrukerWriter, self).__init__(templname)

        self.filename = op.abspath(fname)
        
        # check and build dir trees from bottom
        self.fexpno = op.dirname(op.dirname(self.filename))

    def create_dir_tree(self):
        texpno = op.dirname(op.dirname(template))
        escratch = False    # indicates expno Bruker dir tree has to be created from scratch
        pscratch = False    # indicates procno Bruker dir tree has to be created from scratch
        if not op.exists(fexpno):
            os.makedirs(fexpno)
            escratch = True
            for f in glob.glob(op.join(texpno,'*')):        # copy metadata from template
                if op.isfile(f) and op.basename(f) not in ('ser','fid'):  # but not raw data
                    if debug:   print('**CP**', f, fexpno)
                    shutil.copy( f, op.join(fexpno, op.basename(f)) )
        if not op.exists(filename):
            os.makedirs(filename)
            fscratch = True
            for f in glob.glob(op.join(template,'*')):        # copy metadat from template
                if op.isfile(f) and op.basename(f) not in ('1r','1i','2rr','2ri','2ir','2ii'):  # but not raw data
                    if debug:   print('**CP**', f, filename)
                    shutil.copy( f, op.join(filename, op.basename(f)))

    def write_params(self):
        # now we create Bruker parameter files if creating from scratch and f contains meta data
        # we provide a user warning if it is not possible
        if escratch:
            warn = False
            if self.d.dim == 1:
                pnamelist = ('acqu',)
            elif self.d.dim == 2:
                pnamelist = ('acqu','acqu2')
            for pname in pnamelist:    # create parameter files
                try:
                    par = self.d.params[pname]
                except AttributeError, KeyError:
                    warn = True
                else:
                    write_param(par, op.join(fexpno, pname) )
                    write_param(par, op.join(fexpno, pname+'s') )
            if warn:
                print("Warning, acqu/acqus files have not been updated")
        if fscratch:    # here we create acqu/acqus files
            warn = False
            if self.d.dim == 1:
                pnamelist = ('proc',)
            elif self.d.dim == 2:
                pnamelist = ('proc','proc2')
            for pname in pnamelist:    # create parameter files
                try:
                    par = self.d.params[pname]
                except AttributeError, KeyError:
                    warn = True
                else:
                    write_param(par, op.join(filename, pname) )
                    write_param(par, op.join(filename, pname+'s') )
            if warn:
                print("Warning, proc/procs files have not been updated")
        # load template params
        proc = read_param(find_proc(template, down=False))
        acqu = read_param(find_acqu(texpno))
        if self.d.dim == 2:
            proc2 = read_param(find_proc2(template, down=False))
            acqu2 = read_param(find_acqu2(texpno))

    def scale_data(self):
        # scale between 2^28 and 2^29
        bufabs = abs(self.d.buffer)
        bmax = bufabs.max()
        NC_proc = 0
        while bmax <2**28:
            bmax *= 2
            NC_proc -= 1
        while bmax >2**29:
            bmax /= 2
            NC_proc += 1
        if debug:   print("NC_proc :", NC_proc)
        buffinal = self.d.buffer * (2**(-NC_proc))
    def update_params(self):
        # update a few parameters and write proc files
        if self.d.dim == 1:
            proc['$SI'] = str(self.d.axis1.cpxsize)
            proc['$SF'] = str(self.d.axis1.frequency)
            proc['$SW_p'] = str(self.d.axis1.specwidth)
            proc['$OFFSET'] = str(revoffset(self.d.axis1.offset, acqu, proc))
            proc['$YMAX_p'] = str(buffinal.max())
            proc['$YMIN_p'] = str(buffinal.min())
            write_param(proc, op.join(filename, 'proc') )
            write_param(proc, op.join(filename, 'procs') )
        if self.d.dim == 2:
            proc['$SI'] = str(self.d.axis2.cpxsize)
            proc2['$SI'] = str(self.d.axis1.cpxsize)
            proc['$SF'] = str(self.d.axis2.frequency)
            proc2['$SF'] = str(self.d.axis1.frequency)
            proc['$SW_p'] = str(self.d.axis2.specwidth)
            proc2['$SW_p'] = str(self.d.axis1.specwidth)
            proc['$OFFSET'] = str(revoffset(self.d.axis2.offset, acqu, proc))
            proc2['$OFFSET'] = str(revoffset(self.d.axis1.offset, acqu2, proc2))
            proc['$YMAX_p'] = str(buffinal.max())
            proc['$YMIN_p'] = str(buffinal.min())
            write_param(proc, op.join(filename, 'proc') )
            write_param(proc, op.join(filename, 'procs') )
            write_param(proc2, op.join(filename, 'proc2') )
            write_param(proc2, op.join(filename, 'proc2s') )
    def write_files(self):
        # create binary files
        if self.d.dim == 1:
            if self.d.axis1.itype == 0:
                writebin(op.join(filename, '1r'), buffinal)
            else:
                writebin(op.join(filename, '1r'), buffinal[::2])
                writebin(op.join(filename, '1r'), buffinal[1::2])
        if self.d.dim == 2:
            if self.d.axis2.itype == 0:
                if self.d.axis1.itype == 0:
                    writebin(op.join(filename, '2rr'), buffinal)
                else:
                    writebin(op.join(filename, '1r'), buffinal[::2])
                    writebin(op.join(filename, '1r'), buffinal[1::2])


