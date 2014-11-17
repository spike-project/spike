
# encoding: utf-8
"""
HDF5File.py

Created by Marc-André Delsuc, Marie-Aude Coutouly on 2011-07-13.
Copyright (c) 2011 __NMRTEC__. All rights reserved.

API dealing with HDF5File. For now it is non surclassing tables, you have to use *.hf. to access all tables functionalities
"""

import sys
import os
import unittest
import numpy as np
import tables
import FTICR
import array
import ConfigParser
import time

# File_version is written into the file and tag the file itself 
# to be changed only if the file format changes
__file_version__ = "0.8"
"""
__file_version__ 0.8 axes are gathered in one table, it's a more HDF5 way to deal with informations
__file_version__ 0.7 has the ability to update the axes
__file_version__ 0.6 is first stabilized multiresolution version
"""
# version is the library version number
__version__ = "0.8"
"""
v 0.8 axes are gathered in one table, it's a more HDF5 way to deal with informations
v 0.6 Uses tables.parameter to optimise the use of the cache
"""

def determine_chunkshape(size1, size2):
    """
    returns optimum size for chuncks for a dataset of file size1, size2
    and update cachesize for accomodating dataset
    """
    c1 = int(size1/64. +1)
    c2 = int(size2/64. +1)
#    update_cache_size()
    return (c1, c2)

class HDF5File(object):
    """
    defines the interface to simply (read/write) access HDF5 files
    standard methods are load() and save()

    standard sequence to read is
    H = HDF5File(filename,"r")
    B = H.get_data()      # B is a FTICRdata
    H.close()

    or
    H = HDF5File(filename,"r")
    H.load()
    B = H.data      # B is a FTICRdata
    H.close()


    and to write
    H = HDF5File(filename,"w")
    H.set_data(B)         # where B is a FTICRdata; do not use    H.data = B
    H.save()
    H.close()

    """
    def __init__(self, fname, access = 'r', info = None, nparray = None, fticrd = None, debug = 0):
        # still have to deal with changing the file_version when file is opened reading only
        self.debug = debug
        self.fname = fname
        self.info = None
        self.nparray = None
        self.chunks = None
        self.filters = None    # for compressing : tables.Filters(complevel=1, complib='zlib')
        
        if access not in ("w", "r", "rw", "a","r+"):
            raise Exception(access + " : acces key not valid")
        if os.path.isfile(self.fname):
            self.checkversion()
        if access == "r":     
            if (self.debug > 0):
                print "open in read mode"
            self.access = access
            self.hf = tables.openFile(self.fname, self.access)
            self.get_file_infos()
        elif access == "w":
            if (self.debug > 0):
                print "Open HDF5 File with writing rights"
            self.access = access
            self.hf = tables.openFile(self.fname, self.access)
            self.create_generic()
            if (info is not None):
                if (self.debug > 0):
                    print "Create HDF5 File from info"
                self.info = info
                self.data = None
                self.create_HDF5_info()
            elif (nparray is not None):
                if (self.debug > 0):
                    print "Create HDF5 File from nparray"
                data = FTICR.FTICRData(buffer = nparray)
                self.create_from_template(data)
            elif (fticrd is not None):
                if (self.debug > 0):
                    print "Create HDF5 File from fticrdata"
                self.create_from_template(fticrd)
            else:
                if (self.debug > 0): print "without any data nor infos"
        elif access == "r+":
            if (self.debug > 0):
                print "open in modifying mode r+"
            self.access = "r+"
            self.hf = tables.openFile(self.fname, self.access)
            self.get_file_infos()
        elif access == "rw":
            if (self.debug > 0):
                print "open in modifying mode rw"
            self.access = "r+"
            self.hf = tables.openFile(self.fname, self.access)
            self.get_file_infos()
        elif access == "a":
            if (self.debug > 0):
                print "open in modifying mode a"
            self.access = "a"
            self.hf = tables.openFile(self.fname, self.access)
            self.get_file_infos()
        else:
            raise " Internal error in HDF5 creation - This should never happen"
            #self.f.close()
    #----------------------------------------------
    def checkversion(self):
        """
        check file version and exit if incompatible
        """
        self.hf = tables.openFile(self.fname,"r")
        try:
            infos = self.hf.root.generic_table.read()
        except:
            raise Exception("The file %s does not seem to be a valid file"%self.fname)
        info = infos[0]     # only first entry is used
        self.hf.close()
        if self.debug > 0 : print "File_Version", info["File_Version"]
        
        if info["File_Version"] != __file_version__:
            msg = """
The file {0} is from version {1} while this program handles file version {2}.
You have to update your msh5 file running update functions.
to do this run the following program :

python HDF5File.py update {0}

""".format(self.fname, info["File_Version"],__file_version__)
            raise Exception(msg)
        return True
    #----------------------------------------------
    def set_compression(self, On=False):
        "sets HDF5 file compression to zlib if On is True; to none otherwise"
        if On:
            self.filters = tables.Filters(complevel=1, complib='zlib')
        else:
            self.filters = None
        return On
    #----------------------------------------------
    def set_data(self, data, group='resol1'):
        """
        Take the ser_file and the params and put all the informations in the HDF5File
        """
        self.dim = data.dim
        self.create_tables()
        group_resol = self.createGroup("/", group)
        group_axe = self.createGroup(group_resol, 'axes')
        self.dims = []
        for i in range(self.dim):
            ii = i + 1
            infos = data.axes(ii).__dict__
            self.fill_table(ii, group_axe, infos)
            # For the moment we retrieve the dimension at this point
            # We put them in the self.dims list
            self.dims.append(infos["size"])
        
        #self.data_valid = True
    #----------------------------------------------
    def set_data_from_fticrd(self, buff, group='resol1'):
        """
        sets the FTICRdata attached to the (to be written) file
        """
        self.dim = buff.dim
        self.create_tables()
        group_resol = self.createGroup("/", group)
        group_axe = self.createGroup(group_resol, 'axes')
        self.dims = []
        for i in range(self.dim):
            ii = i + 1
            infos = buff.axes(ii).__dict__
            self.fill_table(ii, group_axe, infos)
            # For the moment we retrieve the dimension at this point
            # We put them in the self.dims list
            self.dims.append(infos["size"])
        self.data = buff.buffer
        self.data_valid = True
    #----------------------------------------------
    def fill_table(self, table, infos):
        """
        Fill in the given table. Axis is the dimension we are processing
        """

        if (self.debug > 0): print "table",table
        rows = table.row
        
        for i in xrange(len(infos)):
            for key in infos[i].keys():
                rows["sampling"] = "uniform"            # default value for sampling
                if key in self.axis_table:
                    rows[key] = infos[i][key]
            
            rows.append()
        table.flush()
        table.close()
    #----------------------------------------------
    def create_generic(self, owner = "NMRTEC"):
        """
        A table is created with all generic informations about the file : owner, method, HDF5 Release,CreationDate, Last modification
        """
        import time
        # Generic infos that don't depend on the axis 
        self.generic_table = {}
        self.generic_table["Owner"] = tables.StringCol(itemsize=16)
        self.generic_table["Method"] = tables.StringCol(itemsize=16)
        self.generic_table["HDF5_Version"] = tables.StringCol(itemsize=16)
        self.generic_table["Version"] = tables.StringCol(itemsize=16)
        self.generic_table["File_Version"] = tables.StringCol(itemsize=16)
        self.generic_table["Creation_Date"] = tables.StringCol(itemsize=16)
        self.generic_table["Last_Modification_Date"] = tables.StringCol(itemsize=16)
        
        generic = self.createTable("/", "generic_table", self.generic_table)
        rows = generic.row
        rows["Owner"] = owner
        rows["Method"] = "FTICR-MS"
        rows["HDF5_Version"] = tables.hdf5Version
        rows["Version"] = __version__
        rows["File_Version"] = __file_version__
        rows["Creation_Date"] = time.strftime("%m/%d/%Y %I:%M:%S %p", time.localtime(os.path.getctime(self.fname)))
        rows["Last_Modification_Date"] = time.strftime("%m/%d/%Y %I:%M:%S %p", time.localtime(os.path.getmtime(self.fname)))
        
        rows.append()
        generic.flush()
        generic.close()
    #----------------------------------------------
    def create_tables(self):
        """
        Creates the different tables needed in a HDF5File according to the info parameters given 
        If you don't pass any info dictionnary, it will take parameters from the self.header
        """
        # Infos that have to be known on each axis 
        self.axis_table = {}
        self.axis_table["size"] = tables.Int32Col()
        self.axis_table["highmass"] = tables.Float32Col()
        self.axis_table["itype"] = tables.Int32Col()
        self.axis_table["left_point"] = tables.Int32Col()
        self.axis_table["ref_freq"] = tables.Float32Col()
        self.axis_table["ref_mass"] = tables.Float32Col()
        self.axis_table["specwidth"] = tables.Float32Col()
        self.axis_table["units"] = tables.StringCol(itemsize=16)
        self.axis_table["sampling"] = tables.StringCol(itemsize=16)
    #----------------------------------------------
    def position_array(self, group="resol1"):
        """
        Fill in the HDF5 file with the given buffer, HDF5 file is created with the given numpy array and the corresponding tables
        """
        if (self.debug > 1): print type(self)
        if self.dim == 2:
            if (self.info is not None):
                sizeF1 = self.info["Dim1"] 
                sizeF2 = self.info["Dim2"]
            else:
                sizeF1 = self.dims[0]
                sizeF2 = self.dims[1]
            if (self.debug > 0): print sizeF1, sizeF2
            if self.chunks is not None:
                chunks = self.chunks
            else:
                chunks = determine_chunkshape(sizeF1, sizeF2)

            if (self.debug > 0): print "chunks", chunks
            Carray = self.createCArray("/"+group, 'data', tables.Float64Atom(), (sizeF1, sizeF2), chunk = chunks)
            # if you have given a numpy array as buffer, it uses it to fillin the data
            if (self.nparray is not None):
                for i1 in xrange(sizeF1):
                    tbuf = self.nparray[i1, 0:sizeF2]
                    Carray[i1, 0:sizeF2] = tbuf[0:sizeF2]
            elif (self.data is not None):
                if (self.debug > 1): print type(self.data)
                for i1 in xrange(sizeF1):
                    tbuf = self.data.buffer[i1, 0:sizeF2]
                    Carray[i1, 0:sizeF2] = tbuf[0:sizeF2]
        elif self.dim == 1:
            if (self.info is not None):
                sizeF1 = self.info["Dim1"] 
            else:
                sizeF1 = self.dims[0]
            Carray = self.createCArray("/"+group, 'data', tables.Float64Atom(), (sizeF1,))

            if (self.nparray is not None):
                Carray = self.nparray[:]
            elif (self.data is not None):
                Carray = self.nparray[:]
            
    #----------------------------------------------
    def get_data(self, group="resol1", mode="onfile"):
        """
        loads and returns the FTICRdata attached to the self file
        same parameters as load()
        """ 
        self.load(group=group, mode=mode)
        return self.data
    #----------------------------------------------        
    def load(self, group="resol1", mode="onfile"):
        """
        loads the data into memory,
        set self.data as a FTICRData 

        group defines which group is loaded (default is resol1)
        mode defines how it is loaded in memory,
            "onfile" (default ) the data is kept on file and loaded only on demand.
                the capability of modifying the data is determined by the way the file was opened
                the data cannot be modified unless the file was opened with access='w' or 'rw'
            "memory" the data is copied to a memroy buffer and can be freely modified
                warning - may saturate the computer memory, there is no control
            if you want to load data into memory after having opened in 'onfile" mode, then do the following :
            h.load(mode="onfile")
            b = d.data.buffer[...]     # data are now copied into a new memory buffer b using ellipsis syntax
            d.data.buffer = b           # and b is used as the data buffer.
        """
        if (self.debug > 0): print group
        hfgroup = getattr(self.hf.root, group)
        try:
            hbuf = hfgroup.data
        except:
            raise Exception("File seems to be corrupted, wrong resolution or not a genuine FTICR file.")
        if mode == "onfile":
            self.data = FTICR.FTICRData(buffer = hbuf) # no copy
        elif mode == "memory":
            self.data = FTICR.FTICRData(buffer = hbuf[...]) # copy !
        else:
            raise Exception("wrong mode, only 'onfile' and 'memory' allowed")
        if (self.debug > 1): print type(self.data.buffer)
        # try:
        # except:
        #     raise Exception(" You might be trying to load a file with another architecture 'maybe different resolutions'")
            
        
        for table in self.hf.walkNodes("/"+group,"Table"): # get the list of tables in the axes group
           values = table.read()
        dtypes = ("highmass","itype","left_point","ref_freq","ref_mass","sampling","size","specwidth","units")
        
        for i in range(len (values)):
            dico = {}
            for j in range(len(values[i])):
                dico[dtypes[j]] = values[i][j]
            setattr(self.data, "axis%d"%(i+1), FTICR.FTICRAxis(highmass=dico["highmass"],left_point = dico["left_point"],size=dico["size"], specwidth=dico["specwidth"], itype=dico["itype"], ref_mass = dico["ref_mass"], ref_freq =dico["ref_freq"],units=dico["units"]))
           
                #self.axes_update(axis = i+1, infos = dico)
        return self.data
        
    #----------------------------------------------
    def get_info(self):
        """
        Retrieve info from self.nparray
        """
        info = {}
        info["Dim"] = len(self.nparray.shape)
        # info["HDF5release"] = 0
        # info["HDF5version"] = 1
        info["Sampling"] = "Uniform"
        info["Dim1"] = self.nparray.shape[0]
        if info["Dim"] > 1:
            info["Dim2"] = self.nparray.shape[1]
            if info["Dim"] > 2:
                info["Dim3"] = self.nparray.shape[2]
        return info
    #----------------------------------------------
    def axes_update(self, group = "resol1", axis = 2, infos = None):
        """routine called when you want to modify the information on a given axis
        group is the group name, default is resol1
        axis is the dimension we want to adjust
        infos is a dictionnary with al fields we want to adjust
        """
        for item in infos :
            self.table_update(group, axis, item, infos[item])
    #----------------------------------------------
    def table_update(self, group = "resol1", axis = 2, key = 'highmass', value = 4000.0 ):
        """Microchangement in the wanted table"""
        for table in self.hf.walkNodes("/"+group,"Table"): # get the list of tables in the axes group
            #print table
            #if (str(table.name) == str(name)):# get the table with the correct name
            #    try :
                
            table.modifyColumn((axis-1), column = value, colname = key)
            #    except :
            #        raise Exception("You might have problem accessing your MSH5 file")

            
    #----------------------------------------------        
    def save(self, ser_file, group="resol1"):
        """
        save the ser_file to the HDF5 file
        """
        import array
    #    import platform # platform seems to be buggy on MacOs, see http://stackoverflow.com/questions/1842544
        if sys.maxint == 2**31-1:   # the flag used by array depends on architecture - here on 32biy
            flag = 'l'              # Apex files are in int32
        else:                       # here in 64bit
            flag = 'i'              # strange, but works here.
            chunks = self.determine_chunkshape()
        Carray = self.createCArray("/"+group, 'data', tables.Float64Atom(), (self.dims[0], self.dims[1]), chunk = chunks)
        with open(ser_file,"rb") as f:
            for i1 in xrange(int(self.dims[0]) ):
                tbuf = f.read(4*int(self.dims[1]) )
                abuf = np.array(array.array(flag, tbuf))
                Carray[i1, 0:int(self.dims[1])] = abuf[0:int(self.dims[1]) ]
        
    #----------------------------------------------        
    def save_fticrd(self):
        """
        save the FTICRData to the H5F file
        """
        if not self.data_valid:
            raise Exception("buffer data not set")
        self.position_array()
    #----------------------------------------------
    def create_from_template(self, data, group='resol1'):
        """
        Take params from the empty FTICR data and put all the informations in the HDF5File
        creates an empty data, and attach it to data.buffer
        data is created in group, with default value 'resol1'
        """
        self.dim = data.dim
        self.create_tables()
        group_resol = self.createGroup("/", group)
        table_axes = self.createTable(group_resol, 'axes', self.axis_table)
        infos = []
        self.dims = []
        for i in range(self.dim):
            ii = i + 1
            infos.append( data.axes(ii).__dict__)
            if (self.debug > 0): print infos
            
            # For the moment we retrieve the dimension at this point
            # We put them in the self.dims list
            self.dims.append(infos[i]["size"])
        #self.data = data       # when you use synthetic data
        self.fill_table(table_axes, infos)
        self.data = None
        self.position_array("/"+group)
        hfgroup = getattr(self.hf.root, group)
        data.buffer = hfgroup.data
        self.data_valid = True
        data.hdf5file = self
    #----------------------------------------------
    def create_HDF5_info(self):
        """
        Creates a HDF5 file, takes info as parameter
        """
        
        # creates  the tables according to the Dimension
        self.create_tables()
        # fills tables with the values given in info
        self.position_array()
    #----------------------------------------------
    def create_HDF5_nparray(self):
        """
        Creates a HDF5 file, takes nparray as parameters
        """
        # creates self.info from self.nparray 
        # it can just retrieve Dim info for now
        self.info = self.get_info()
        # creates  the tables according to the Dimension
        self.create_tables()
        # fills tables and carray with the values given self.info and self.nparray
        self.position_array()
    #----------------------------------------------
    def close(self):
        """
        Closes HDF5File
        """
        if (self.debug > 0) : print "ABOUT TO CLOSE ", self.fname
        self.hf.close()
        if (self.debug > 0) : print "IT's CLOSED ", self.fname
    #----------------------------------------------
    def createGroup(self, where, name):
        """
        Create a group in the given hf_file
        """ 
        group = self.hf.createGroup(where, name)
        return group
    #----------------------------------------------
    def createTable(self, where, name, description):
        """
        Create a Table in the given hf_file at the given position with the right description
        """
        table = self.hf.createTable( where, name, description)
        return table
    #----------------------------------------------
    def createCArray(self, where, name, data_type, shape, chunk = None):
        """
        Create a CArray in the given hf_file
        """
        if len(shape) == 2:             # If we have a 2D spectrum
            if not chunk:
                chunk = determine_chunkshape(shape[0], shape[1])
        # filters is set from outside, see self.set_compression()
        array = self.hf.createCArray(where, name, data_type, shape,  chunkshape = chunk,  filters = self.filters)
        return array
    
    #----------------------------------------------
    def determine_chunkshape(self, sizeF1 = None, sizeF2= None):    # recopié ici et passé en méthode, juste pour jouer avec la taille du chunck
        """
        Determine a good chunkshape according to the size of each axis
        """
        if (self.info is not None):
            sizeF1 = self.info["Dim1"] 
            sizeF2 = self.info["Dim2"]
        elif (sizeF1 is None):
            sizeF1 = self.dims[0]
            sizeF2 = self.dims[1]
        sz12 = float(sizeF2) / float(sizeF1)
        n2 = 1
        n1 = 1024
        while  ( (float(n2)/ float(n1)) < sz12) and (n1 > 1 ) :
            n1 = n1/2
            n2 = n2*2
            if (self.debug > 0): print "si1 x si2 : %d %d   n1 x n2 : %d %d" % (float(sizeF1), float(sizeF2), n1, n2)
        return(n1, n2)
    #----------------------------------------------
    def get_file_infos(self):
        """
        Read the generic_table and return the informations
        """
        infos = self.hf.root.generic_table.read()
        print "******************* \n %s is a %s file created by %s on %s with code version %s.\n HDF5 version is %s and API version is %s.\n Last modification has been made %s \n********************"% ( self.fname, infos[0]['Method'], infos[0]['Owner'], infos[0]['Creation_Date'], infos[0]['Version'], infos[0]['HDF5_Version'], infos[0]['File_Version'], infos[0]['Last_Modification_Date'])
################################################
# Update functions
def update(fname, debug = 1):
    """update so that the file is up to date"""
    hf = tables.openFile(fname,"r")
    try:
        infos = hf.root.generic_table.read()
    except:
        raise Exception("The file %s does not seem to be a valid file"%fname)
    info = infos[0]     # only first entry is used
    fileversion = float(info["File_Version"])
    if debug > 0 : print "File_Version", fileversion
    hf.close()
    hf.close()
    if fileversion == float(__file_version__):
        print "File is up to date, there is nothing to do"
    if  fileversion < 0.6:
        raise Exception("The file %s is from version %s, which is too old to be updated, sorry"%(fname, fileversion))
    if fileversion < 0.7:
        print "updating from 0.6 to 0.7"
        up0p6_to_0p7(fname, debug=debug)
        fileversion = 0.7
    if fileversion < 0.8:
        print "updating from 0.7 to 0.8"
        up0p7_to_0p8(fname, debug=debug)
        fileversion = 0.8
    print "The file %s has been fully updated to %s"%(fname,__file_version__)
#----------------------------------------------
def up0p6_to_0p7(fname, debug = 1):
    """docstring for up0p6_to_0p7
    Function that deals with changing HDF5 files created with file_version 0.6 to be read with 0.7
    It modifies 
    """
    import time
    print "modifying dtype of highmass"
    hf = tables.openFile(fname,"r+")
    for table in hf.walkNodes("/","Table"): # get the list of tables in the axes group
        try :
            table.description.highmass.type = 'float32'
        except:
            print " no highmass field in this table"
    hf.root.generic_table.modifyColumn(0, column = "0.7", colname = "File_Version")
    hf.root.generic_table.modifyColumn(0, column = time.strftime("%m/%d/%Y %I:%M:%S %p", time.localtime()), colname = "Last_Modification_Date")
    
    hf.close()
#----------------------------------------------
def up0p7_to_0p8(fname, debug = 1):
    """docstring for up0p7_to_0p8
    Function that deals with changing HDF5 files created with file_version 0.7 to be read with 0.8
    """
    hf = tables.openFile(fname,"r+")
    
    description = ""
    for group in hf.iterNodes("/","Group"): # get the list of tables in the file
        infos = []
        for table in hf.iterNodes(group.axes,"Table"): # get the list of tables in the file
            if (table.name != "generic_table"):
                infos.append(table.read())
                description = table.description
        hf.removeNode(group,"axes",True)

        table = hf.createTable (group, "axes" , description)
        for i in xrange(len(infos)):
            infos[i]["sampling"] = "uniform"
            table.append(infos[i])
        table.flush()
        table.close()
    
    hf.root.generic_table.modifyColumn(0, column = "0.8", colname = "File_Version")
    hf.root.generic_table.modifyColumn(0, column = time.strftime("%m/%d/%Y %I:%M:%S %p", time.localtime()), colname = "Last_Modification_Date")
    if (debug > 0):
        print "We have to gather all axes in one table"
    hf.close()

#----------------------------------------------
class HDF5_Tests(unittest.TestCase):

    def setUp(self):
        rootfiles = os.getcwd()
        # make ini file from example
        self.TestFolder = os.path.join('..','DATA_test')
        self.DataFolder = os.path.join(self.TestFolder, 'cytoC_2D_000001.d')
        self.name_write = os.path.join(self.TestFolder,"file_write")
        self.name_fticr = os.path.join(self.TestFolder,"file_fticr")
        self.name_npar = os.path.join(self.TestFolder,"file_npar")
        self.npar_fticr = os.path.join(self.TestFolder,"npar_fticr")
        self.name_chunk = os.path.join(self.TestFolder,"file_fticr")
        self.name_get = os.path.join(self.TestFolder,"file_fticr")
                
        self.verbose = 1    # verbose > 0 switches messages on
    def announce(self):
        if self.verbose > 0:
            print "\n========", self.shortDescription(), '==============='
    def test_get_data(self):
        "Test routine that opens a HDF5 file reading, gets the headers and the buffer"
        self.announce()
        print "testing get_data the type file"
        
        hdf5 = HDF5File(self.name_get,"rw")
        hdf5.load()
        B = hdf5.data
        hdf5.close()
    #----------------------------------------------
    def test_create_from_fticr(self):
        "Test routine that creates a HDF5 file according to a given FTICRData"
        self.announce()
        import Apex  as ap
        fticrdata = ap.Import_2D(self.DataFolder)
        h5f = HDF5File(self.name_fticr, "w", fticrd = fticrdata, debug =1)
        h5f.close()
    #----------------------------------------------
    def test_create_from_nparray(self):
        "Testing the creation of a HDF5 file according to a given nparray"
        self.announce()
        data_init = 10*np.random.random((2048, 65536))  # example of buffer (numpy array) from which you might create a HDF5 file
        h5f = HDF5File(self.name_npar, "w", nparray = data_init, debug =1)
        h5f.close()
    #----------------------------------------------
    def test_nparray_to_fticr(self):
        "Test routine that creates a HDF5 file according to a given nparray and returns the buffer (FTICRData)"
        self.announce()
        data_init = 10*np.random.random((2048, 65536)) 
        B = nparray_to_fticrd(self.npar_fticr, data_init)
        print B.axis1.size
    #----------------------------------------------
    def test_axes_update(self):
        "Test routine that overloads the parameter from an axis"
        import Apex  as ap
        self.announce()
        import time
        
        d = HDF5File(self.name_get,"rw",debug = 1)
        d.axes_update(axis = 2, infos = {'highmass':4200.00, 'units':'m/z'})
        d.close()
    #----------------------------------------------
    def _test_chunkshape(self):
        "Test routine that creates , processes Data according to chunkshapes"
        self.announce()
        from time import time
        from Apex import Import_2D
   
        # Dim1 = 8192
        # Dim2 = 16384
        # 
        # data = 10*np.random.random((Dim1, Dim2))
        # fticr = FTICRData(buffer = data)
        # print fticr.buffer
        # t1 = time()
        # HF = HDF5File("/Users/mac/Desktop/Chunk.hf5","w")
        # HF.get_file_infos()
        # HF.chunks = (25,100)
        # HF.create_from_template(fticr)
        # fticr.hdf5file = HF
        # 
        
        #print HF.hf.root.data[:]
        t1 = time()
        d = Import_2D(self.DataFolder, self.name_chunk)
        print "---------------SETTINGS---------------------"
        print "CHUNK_CACHE_PREEMPT " , tables.parameters.CHUNK_CACHE_PREEMPT 
        print "CHUNK_CACHE_SIZE ", tables.parameters.CHUNK_CACHE_SIZE 
        print "METADATA_CACHE_SIZE ", tables.parameters.METADATA_CACHE_SIZE 
        print "NODE_CACHE_SLOTS ", tables.parameters.NODE_CACHE_SLOTS
        print "---------------SETTINGS---------------------"
       #d.hdf5file.close()
       #d = HF.get_data()
        print "import", time()-t1, "secondes"
        t0 = time()
        t00 = t0
        d.rfft(axis=2)
        print "rfft2", time()-t0, "secondes"
       # HF = HDF5File("/Users/mac/Desktop/Chunk.hf","r")
       # d = HF.get_data()
        t0 = time()
        # d.rfft(axis=2)
        # print "rfft2",time()-t0,"secondes"
        # #print d.buffer.__dict__
        # #HF.close()
        # #print d.buffer[:]
        # t0 = time()
        # #print type(d)
        d.rfft(axis=1)
        print "rfft1", time()-t0, "secondes"
        t0 = time()
        d.modulus()
        print "modulus", time()-t0, "secondes"
        print "calcul", time()-t00, "secondes"


#----------------------------------------------
def nparray_to_fticrd(name, nparray):
    """
    """
    h5f = HDF5File(name, "w", None, nparray)
    h5f.get_file_infos()
    h5f.close()
    h5f2 = HDF5File(name, "r")
    B = h5f2.get_data()
    h5f2.close()
    return B
    
#----------------------------------------------
def syntax(prgm):
    print """
****** ERROR ******
syntax is
{0}
   or 
{0} Tests     (both run self tests)
   or
{0} update file_name    to update a old file to current
""".format(prgm)
    sys.exit(1)
#----------------------------------------------
if __name__ == '__main__':
    """ can be used as a prgm for updating files """
    # get
    argv = sys.argv
    try:
        action = argv[1]
    except:
        action = 'Tests'

    # do
    if action == 'update':
        try:
            fname = argv[2]
        except:
            syntax(argv[0])
        update(fname, debug = 1)
    elif action == "Tests":
        print "running tests"
        argv.pop()
        unittest.main()
    else:
        syntax(argv[0])
