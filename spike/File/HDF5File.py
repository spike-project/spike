#!/usr/bin/env python 
# encoding: utf-8

"""
HDF5File.py

Created by Marc-André Delsuc, Marie-Aude Coutouly on 2011-07-13.

API dealing with HDF5File. For now it is non surclassing tables, you have to use *.hf. to access all tables functionalities
"""

from __future__ import print_function
import sys
import os
import unittest
import numpy as np
import time
import array
import tempfile

import json
import tables
from tables.nodes import filenode

if sys.version_info[0] < 3:
    import ConfigParser
else:
    import configparser as ConfigParser
    xrange = range


# File_version is written into the file and tag the file itself 
# to be changed only if the file format changes
if sys.version_info[0] < 3:
    __file_version__ = "0.9"
else:
    __file_version__ = b"0.9"
"""
__file_version__ 0.9 new list of attributes for FTMS axes
    additional files are now stored along the data in "attached" group.
__file_version__ 0.89 transitionnal between 0.8 and 0.9 
__file_version__ 0.8 axes are gathered in one table, it's a more HDF5 way to deal with informations
__file_version__ 0.7 has the ability to update the axes
__file_version__ 0.6 is first stabilized multiresolution version
"""
# version is the library version number
__version__ = "0.901"
"""
v 0.901  when storing a string 'None' it will be read back as the python object None - no change on the file format
v 0.9 porting from pytables v2 to v3.x / new list of attributes for FTMS axes / storage for files and objects / better compression
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

# The dictionary used to define axes tables
FTICR_AXISvp9 = {
    "itype":tables.Int32Col(),
    "size":tables.Int32Col(),
    "FTICR":tables.StringCol(itemsize=16),
    "sampling":tables.StringCol(itemsize=16),
    "specwidth":tables.Float32Col(),
    "highmass":tables.Float32Col(),
    "offsetfreq":tables.Float32Col(),
    "left_point":tables.Int32Col(),
    "calibA":tables.Float64Col(),
    "calibB":tables.Float64Col(),
    "calibC":tables.Float64Col(),
    "highfreq":tables.Float32Col(),
    "lowfreq":tables.Float32Col(),
    }

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

    
    HDF5File have the capacity to store and retrieve complete files and python objects:
        with
    lis = [any kind of list or tuple]  # works also with dict and nested list and dict
        then
    H.store_internal_object(lis, "name_of_storage")
        will store the object, and
    lis_back = H.retrieve_object("name_of_storage")
        will retrieve it
    data are stored using JSON, so anything compatible will do  
    """
    def __init__(self, fname, access = 'r', info = None, nparray = None, fticrd = None, compress=False, debug = 0, verbose=False):
        """
        access:
            r: Read-only; no data can be modified.
            w: Write; a new file is created (an existing file with the same name would be deleted).
            a: Append; an existing file is opened for reading and writing, and if the file does not exist it is created.
            r+: It is similar to ‘a’, but the file must already exist.
        """
        # still have to deal with changing the file_version when file is opened reading only
        from  .. import FTICR
        #import spike.FTICR as FTICR
        import getpass
        try:
            owner = getpass.getuser()
        except:
            owner = None
        self.debug = debug
        self.fname = fname
        self.info = None
        self.nparray = None
        self.chunks = None
        self.filters = None    # for compressing : tables.Filters(complevel=1, complib='zlib')
        self.set_compression(compress)
        if access not in ("w", "r", "rw", "a", "r+"):
            raise Exception(access + " : acces key not valid")
        if os.path.isfile(self.fname):
            self.checkversion()
        if access == "r":     
            if (self.debug > 0):
                print("open in read mode")
            self.access = access
            self.hf = tables.open_file(self.fname, self.access)
            if verbose: self.get_file_infos()
        elif access == "w":
            if (self.debug > 0):
                print("Open HDF5 File with writing rights")
            self.access = access
            self.hf = tables.open_file(self.fname, self.access)
            # generic_table contains info on the file
            self.create_generic(owner=owner)
            # create a "attached" group for storing files using filenode mode - compressed -
            self.hf.create_group("/",'attached',filters=tables.Filters(complevel=1, complib='zlib'))
            if (info is not None):
                if (self.debug > 0):
                    print("Create HDF5 File from info")
                self.info = info
                self.data = None
                self.create_HDF5_info()
            elif (nparray is not None):
                if (self.debug > 0):
                    print("Create HDF5 File from nparray")
                data = FTICR.FTICRData(buffer = nparray)
                self.create_from_template(data)
            elif (fticrd is not None):
                if (self.debug > 0):
                    print("Create HDF5 File from fticrdata")
                self.create_from_template(fticrd)
            else:
                if (self.debug > 0): print("without any data nor infos")
        elif access == "r+":
            if (self.debug > 0):
                print("open in modifying mode r+")
            self.access = "r+"
            self.hf = tables.open_file(self.fname, self.access)
            if verbose: self.get_file_infos()
        elif access == "rw":
            if (self.debug > 0):
                print("open in modifying mode rw")
            self.access = "r+"
            self.hf = tables.open_file(self.fname, self.access)
            if verbose: self.get_file_infos()
        elif access == "a":
            if (self.debug > 0):
                print("open in modifying mode a")
            self.access = "a"
            self.hf = tables.open_file(self.fname, self.access)
            if verbose: self.get_file_infos()
        else:
            raise " Internal error in HDF5 creation - This should never happen"
            #self.f.close()
    #----------------------------------------------
    def checkversion(self):
        """
        check file version and exit if incompatible
        """
        self.hf = tables.open_file(self.fname,"r")
        try:
            infos = self.hf.root.generic_table.read()
        except:
            raise Exception("The file %s does not seem to be a valid file"%self.fname)
        info = infos[0]     # only first entry is used
        self.hf.close()
        if self.debug > 0 : print("File_Version", info["File_Version"])
        
        if info["File_Version"] != __file_version__:
            msg = """

WARNING
The file {0} is from version {1} while this program handles file version {2}.
You have to upgrade your msh5 file applying the update function.
to do this run the following spike command:

spike.File.HDF5File.update("{0}")

or by running the following program: ( adapt to your local spike setup ) :

python -m spike.File.HDF5File update {0}
""".format(self.fname, info["File_Version"],__file_version__)
            raise Exception(msg)
        return True

    ######### file nodes ################
    # store arbitrary objects in a hdf5 node
    #----------------------------------------------
    def open_internal_file(self, h5name, access='r', where='/attached'):
        """
        opens a node called h5name in the file, which can be accessed as a file.
        returns a file stram which can be used as a classical file.
        
        access is either 
            'r' : for reading an existing node
            'w' : create a node for writing into it
            'a' : for appending in an existing node
        file is stored in a h5 group called h5name

        eg.
        F = h5.open_internal_file('myfile.txt', 'w', where='/files')
        # create a node called '/files/myfile.txt' (node 'myfile.txt' in the group '/files')
        F.writelines(text)
        F.close()
        # and write some text into it

        # then, latter on
        F = h5.open_internal_file('myfile.txt', 'r', where='/files')
        textback = F.read()
        F.close()

        This is used to add parameter files, audit_trail, etc... to spike/hdf5 files
        
        it is based on the filenode module from pytables
        """
        import warnings
        if access == 'r':
            v  = self.hf.get_node(where=where, name=h5name)
            F = filenode.open_node(v, 'r')
        elif access == 'a':
            v  = self.hf.get_node(where=where, name=h5name)
            F = filenode.open_node(v, 'a+')
        elif access == 'w':
            with warnings.catch_warnings():   # remove warnings, as the dot in name activates them
                warnings.simplefilter("ignore")
                F = filenode.new_node(self.hf, where=where, name=h5name)
        return F
        
    def store_internal_file(self, filename, h5name=None, where='/attached'):
        """
        Store a (text) file into the hdf5 file,
            filename: name of the file to be copied
            h5name: is its internal name (more limitation than in regular filesystems)
                copied from os.path.basename(filename) by default
            where: group where the file is copied into the hdf5
        file content will be retrieved using    open_internal_file(h5name,'r)
        """
        if h5name is None:
            h5name = os.path.basename(filename)
        node = self.open_internal_file(h5name, 'w', where=where)
        with open(filename, 'rb') as F:
            node.write( F.read() )
        node.close()
    #----------------------------------------------
    def retrieve_internal_file(self, h5name, where='/attached'):
        """
        returns the content of a n internal file stored with store_internal_file() or written directly
        """
        v  = self.hf.get_node(where=where, name=h5name)
        F = filenode.open_node(v,'r')
        content = F.read()
        F.close()
        return content
    #----------------------------------------------
    def store_internal_object(self, obj, h5name, where='/'):
        """
        store a python object into the hdf5 file
        object are then retrieve with retrieve_object()

        uses JSON to serialize obj
            - so works only on values, lists, dictionary, etc... but not functions or methods
        """
        node = filenode.new_node(self.hf, where=where, name=h5name, filters=tables.Filters(complevel=1, complib='zlib'))
        js = json.dumps(obj, ensure_ascii=True)
        node.write(js.encode())
        node.close()
    #----------------------------------------------
    def retrieve_object(self, h5name, where='/', access='r'):
        """
        retrieve a python object stored with store_internal_object()
        """
        v  = self.hf.get_node(where=where, name=h5name)
        F = filenode.open_node(v,'r')
        js = F.read()
        F.close()
        obj = json.loads(js.decode())
        return obj
    ###############################################
    def set_compression(self, On=False):
        "sets Carray HDF5 file compression to zlib if On is True; to none otherwise"
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
        group_resol = self.create_group("/", group)
        group_axe = self.create_group(group_resol, 'axes')
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
        group_resol = self.create_group("/", group)
        group_axe = self.create_group(group_resol, 'axes')
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

        if (self.debug > 0): print("table",table)
        rows = table.row
        
        for i in xrange(len(infos)):
            for key in infos[i].keys():
#                rows["sampling"] = "uniform"            # default value for sampling
                if key in self.axis_table:
                    rows[key] = infos[i][key]
            
            rows.append()
        table.flush()
#        table.close()
    #----------------------------------------------
    def create_generic(self, owner=None):
        """
        A table is created with all generic informations about the file : owner, method, HDF5 Release, CreationDate, Last modification
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
        
        generic = self.create_table("/", "generic_table", self.generic_table)
        rows = generic.row
        if owner is None:
            rows["Owner"] = "Unknown"
        else:
            rows["Owner"] = owner
        rows["Method"] = "FTICR-MS"
        rows["HDF5_Version"] = tables.hdf5_version
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
        Creates the different tables needed in a HDF5File/FTICR  
        """
        # MAD TO BE ADAPTED !!!
        # Infos that have to be known on each axis 
        
        self.axis_table = FTICR_AXISvp9
#         self.axis_table = {}
#         self.axis_table["itype"] = tables.Int32Col()
#         self.axis_table["size"] = tables.Int32Col()
# #        self.axis_table["units"] = tables.StringCol(itemsize=16)
#         self.axis_table["FTICR"] = tables.StringCol(itemsize=16)
#         self.axis_table["sampling"] = tables.StringCol(itemsize=16)
#         self.axis_table["specwidth"] = tables.Float32Col()
#         self.axis_table["highmass"] = tables.Float32Col()
#         self.axis_table["offsetfreq"] = tables.Float32Col()
#         self.axis_table["left_point"] = tables.Int32Col()
#         self.axis_table["calibA"] = tables.Float64Col()
#         self.axis_table["calibB"] = tables.Float64Col()
#         self.axis_table["calibC"] = tables.Float64Col()
#         self.axis_table["highfreq"] = tables.Float32Col()
#         self.axis_table["lowfreq"] = tables.Float32Col()

    #----------------------------------------------
    def position_array(self, group="resol1"):
        """
        Fill in the HDF5 file with the given buffer, HDF5 file is created with the given numpy array and the corresponding tables
        """
        if (self.debug > 1): print(type(self))
        if self.dim == 2:
            if (self.info is not None):
                sizeF1 = self.info["Dim1"] 
                sizeF2 = self.info["Dim2"]
            else:
                sizeF1 = self.dims[0]
                sizeF2 = self.dims[1]
            if (self.debug > 0): print(sizeF1, sizeF2)
            if self.chunks is not None:
                chunks = self.chunks
            else:
                chunks = determine_chunkshape(sizeF1, sizeF2)

            if (self.debug > 0): print("chunks", chunks)
            Carray = self.create_carray("/"+group, 'data', tables.Float64Atom(), (sizeF1, sizeF2), chunk = chunks)
            # if you have given a numpy array as buffer, it uses it to fillin the data
            if (self.nparray is not None):
                for i1 in xrange(sizeF1):
                    tbuf = self.nparray[i1, 0:sizeF2]
                    Carray[i1, 0:sizeF2] = tbuf[0:sizeF2]
            elif (self.data is not None):
                if (self.debug > 1): print(type(self.data))
                for i1 in xrange(sizeF1):
                    tbuf = self.data.buffer[i1, 0:sizeF2]
                    Carray[i1, 0:sizeF2] = tbuf[0:sizeF2]
        elif self.dim == 1:
            if (self.info is not None):
                sizeF1 = self.info["Dim1"] 
            else:
                sizeF1 = self.dims[0]
            Carray = self.create_carray("/"+group, 'data', tables.Float64Atom(), (sizeF1,))

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
        from  .. import FTICR
        if (self.debug > 0): print(group)
        hfgroup = getattr(self.hf.root, group)
        #print('############ ', dir(hfgroup))
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
        if (self.debug > 1): print(type(self.data.buffer))

#  Nw axes details
#        for table in self.hf.walk_nodes("/"+group,"Table"): # get the list of tables in the axes group
#          values = table.read()
#           print(values)
#           exit()
        table = getattr(hfgroup,"axes")
        values = table.read()   # load all axes
        fields = [i[0] for i in values.dtype.descr]  # rename / rewrite of dtype above
        for i in range(len (values)): # for each axis
            # dico = {}
            # for j in range(len(values[i])):   # for each entry in the axis table load into "dico"
            #     dico[fields[j]] = values[i][j]
            ax = FTICR.FTICRAxis()  # create empty axis
            for j in range(len(values[i])):    # for all entries in axis table
                vv = values[i][j]                       # fetch value
                if vv == 'None':  vv = None             # the string None codes for the object None !!!
                setattr(ax, fields[j], vv)              # copy table entries into axis attribute
            setattr(self.data, "axis%d"%(i+1), ax)      # and set axis
        
        self.data.adapt_size()      # axis sizes have to be updated
        
        #finally load params object
        try:
            self.data.params = self.retrieve_object('params')
        except:
            print("WARNING: %s file does not have the params attribute"%self.fname)
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
        for table in self.hf.walk_nodes("/"+group,"Table"): # get the list of tables in the axes group
            #print table
            #if (str(table.name) == str(name)):# get the table with the correct name
            #    try :
                
            table.modify_column((axis-1), column = value, colname = key)
            #    except :
            #        raise Exception("You might have problem accessing your MSH5 file")

            
    #----------------------------------------------        
    def save(self, ser_file, group="resol1"):
        """
        save the ser_file to the HDF5 file
        """
        import array
    #    import platform # platform seems to be buggy on MacOs, see http://stackoverflow.com/questions/1842544
        if sys.maxsize  == 2**31-1:   # the flag used by array depends on architecture - here on 32biy
            flag = 'l'              # Apex files are in int32
        else:                       # here in 64bit
            flag = 'i'              # strange, but works here.
            chunks = self.determine_chunkshape()
        Carray = self.create_carray("/"+group, 'data', tables.Float64Atom(), (self.dims[0], self.dims[1]), chunk = chunks)
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
        group_resol = self.create_group("/", group)
        table_axes = self.create_table(group_resol, 'axes', self.axis_table, )
        infos = []
        self.dims = []
        for i in range(self.dim):
            ii = i + 1
            infos.append( data.axes(ii).__dict__)
            if (self.debug > 0): print(infos)
            
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
    def flush(self):
        """
        flushes all nodes but does not close
        """
        if (self.debug > 0) : print("flushing ", self.fname)
        for node in self.hf.walk_nodes("/","Table"):
            node.flush()
        for node in self.hf.walk_nodes("/","Array"):
            node.flush()
        
    #----------------------------------------------
    def close(self):
        """
        Closes HDF5File
        """
        if (self.debug > 0) : print("ABOUT TO CLOSE ", self.fname)
        self.hf.close()
        if (self.debug > 0) : print("IT's CLOSED ", self.fname)
    #----------------------------------------------
    def create_group(self, where, name):
        """
        Create a group in the given hf_file
        """ 
        group = self.hf.create_group(where, name)
        return group
    #----------------------------------------------
    def create_table(self, where, name, description):
        """
        Create a Table in the given hf_file at the given position with the right description
        """
        table = self.hf.create_table( where, name, description)
        return table
    #----------------------------------------------
    def create_carray(self, where, name, data_type, shape, chunk = None):
        """
        Create a CArray in the given hf_file
        """
        if len(shape) == 2:             # If we have a 2D spectrum
            if not chunk:
                chunk = determine_chunkshape(shape[0], shape[1])
        # filters is set from outside, see self.set_compression()
        array = self.hf.create_carray(where, name, data_type, shape,  chunkshape = chunk,  filters = self.filters)
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
            if (self.debug > 0): print("si1 x si2 : %d %d   n1 x n2 : %d %d" % (float(sizeF1), float(sizeF2), n1, n2))
        return(n1, n2)
    #----------------------------------------------
    def get_file_infos(self):
        """
        Read the generic_table and return the informations
        """
        infos = self.hf.root.generic_table.read()
        print("******************* \n %s is a %s file created by %s on %s with file version %s.\n HDF5 version is %s and API version is %s.\n Last modification has been made %s \n********************"% ( self.fname, infos[0]['Method'], infos[0]['Owner'], infos[0]['Creation_Date'], infos[0]['File_Version'], infos[0]['HDF5_Version'], infos[0]['Version'], infos[0]['Last_Modification_Date']))

################################################
# Update functions
def update(fname, debug = 1):
    """update so that the file is up to date"""
    hf = tables.open_file(fname,"r")
    try:
        infos = hf.root.generic_table.read()
    except:
        raise Exception("The file %s does not seem to be a valid file"%fname)
    info = infos[0]     # only first entry is used
    fileversion = float(info["File_Version"])
    if debug > 0 : print("File_Version", fileversion)
    hf.close()
    hf.close()
    if fileversion == float(__file_version__):
        print("File is up to date, there is nothing to do")
    if  fileversion < 0.6:
        raise Exception("The file %s is from version %s, which is too old to be updated, sorry"%(fname, fileversion))
    if fileversion < 0.7:
        print("updating from 0.6 to 0.7")
        up0p6_to_0p7(fname, debug=debug)
        fileversion = 0.7
    if fileversion < 0.8:
        print("updating from 0.7 to 0.8")
        up0p7_to_0p8(fname, debug=debug)
        fileversion = 0.8
    if fileversion < 0.9:
        print("updating from 0.8 to 0.9")
        up0p8_to_0p9(fname, debug=debug)
        fileversion = 0.9
        print("The file %s has been fully updated to %s"%(fname,__file_version__))
#----------------------------------------------
def up0p6_to_0p7(fname, debug = 1):
    """docstring for up0p6_to_0p7
    Function that deals with changing HDF5 files created with file_version 0.6 to be read with 0.7
    It modifies 
    """
    import time
    print("modifying dtype of highmass")
    hf = tables.open_file(fname,"r+")
    for table in hf.walk_nodes("/","Table"): # get the list of tables in the axes group
        try :
            table.description.highmass.type = 'float32'
        except:
            print(" no highmass field in this table")
    hf.root.generic_table.modify_column(0, column = "0.7", colname = "File_Version")
    hf.root.generic_table.modify_column(0, column = time.strftime("%m/%d/%Y %I:%M:%S %p", time.localtime()), colname = "Last_Modification_Date")
    
    hf.close()
#----------------------------------------------
def up0p7_to_0p8(fname, debug = 1):
    """docstring for up0p7_to_0p8
    Function that deals with changing HDF5 files created with file_version 0.7 to be read with 0.8
    """
    hf = tables.open_file(fname,"r+")
    
    description = ""
    for group in hf.iter_nodes("/","Group"): # get the list of tables in the file
        infos = []
        for table in hf.iter_nodes(group.axes,"Table"): # get the list of tables in the file
            if (table.name != "generic_table"):
                infos.append(table.read())
                description = table.description
        hf.remove_node(group,"axes",True)

        table = hf.create_table (group, "axes" , description)
        for i in xrange(len(infos)):
            infos[i]["sampling"] = "uniform"
            table.append(infos[i])
        table.flush()
        table.close()
    
    hf.root.generic_table.modify_column(0, column = "0.8", colname = "File_Version")
    hf.root.generic_table.modify_column(0, column = time.strftime("%m/%d/%Y %I:%M:%S %p", time.localtime()), colname = "Last_Modification_Date")
    if (debug > 0):
        print("We have to gather all axes in one table")
    hf.close()
#----------------------------------------------
def up0p8_to_0p9(fname, debug = 1):
    """
    Function that deals with changing HDF5 files created with file_version 0.8 to be read with 0.9 lib
    """
    hf = tables.open_file(fname,"r+")

    # creates /attached/
    hf.create_group("/", 'attached', filters=tables.Filters(complevel=1, complib='zlib'))

    # now modify axes tables

    for group in hf.iter_nodes("/","Group"): # get the list of tables in the file
        print("GROUP",group._v_name)
        if group._v_name.startswith('resol'):
            axesv8 = getattr(group,"axes")
            print(axesv8)
            # do unit conversion
            newaxes = []
            for iax, ax in enumerate(axesv8):
                newax = {}
                # from 0.! to 0.9 FTMS freq coordinate system has changed, 
                # so here we translate from ond to new (where sw is invariant)
                sw = axesv8[iax]['specwidth']
                size = axesv8[iax]['size']
                left_point = axesv8[iax]['left_point']
                newax['left_point'] = 0
                newax['offsetfreq'] = left_point*sw/(size-1+left_point)  # only approximate to a few Hz !
                newax['specwidth'] = sw - newax['offsetfreq']  # and remove this valuee
                newax['calibA'] = axesv8[iax]['ref_mass']*axesv8[iax]['ref_freq']
                newax['calibB'] = 0.0
                newax['calibC'] = 0.0
                newax['FTICR'] = 'FTICR'
                newax['highfreq'] = newax['specwidth']
                newax['lowfreq'] = newax['calibA']/axesv8[iax]["highmass"]
                # the following are common to v08 and v09
                for key in ["size", "highmass", "itype", "sampling"]:
                    newax[key] = axesv8[iax][key]
#                print (iax, newax)
                newaxes.append(newax)
            # remove old axes
            hf.remove_node(group,"axes",True)
            # create a new on
            axtable = hf.create_table (group, "axes" , FTICR_AXISvp9)
            # and fills it
            rows = axtable.row
            for ax in newaxes:
                for key in ax.keys():
                    rows[key] = ax[key]
                rows.append()
            axtable.flush()
    
    # update version
    hf.root.generic_table.modify_column(0, column = "0.9", colname = "File_Version")
    hf.root.generic_table.modify_column(0, column = time.strftime("%m/%d/%Y %I:%M:%S %p", time.localtime()), colname = "Last_Modification_Date")
    print("file update succesful")
    hf.close()

#----------------------------------------------
class HDF5_Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        "This one is called before running all the tests"
        from ..Tests import filename, directory
        rootfiles = os.getcwd()
        print ('Setting up tests')
        cls.TestFolder = directory()
        cls.DataFolder = filename('ubiquitine_2D_000002.d')#'cytoC_2D_000001.d')
        cls.name_file1 = filename("test_file.msh5")
        cls.name_file2 = filename("test_file2.msh5")
        cls.name_chunk = filename("test_chunk.msh5")
        # Testing the creation of a HDF5 file according to a given nparray - and creating for next tests
        data_init = 10*np.random.random((2048, 32*1024))  # example of buffer (numpy array) from which you might create a HDF5 file
        threshold = data_init.std()
        data_init[abs(data_init)<threshold] = 0.0   # but a large heap or zeros to test compression
        h5f = HDF5File(cls.name_file1, "w", nparray = data_init, debug =1, compress=True)
        h5f.close()
        cls.verbose = 1    # verbose > 0 switches messages on
    def announce(self):
        if self.verbose > 0:
            print("\n========", self.shortDescription(), '===============')
    def test_get_data(self):
        "Test routine that opens a HDF5 file reading, gets the headers and the buffer"
        self.announce()        
        hdf5 = HDF5File(self.name_file1,"r")
        hdf5.load()
        B = hdf5.data
        hdf5.close()
    #----------------------------------------------
    def test_create_from_fticr(self):
        "Test routine that creates a HDF5 file according to a given FTICRData"
        self.announce()
        from . import Apex  as ap
        fticrdata = ap.Import_2D(self.DataFolder)
        h5f = HDF5File(self.name_file2, "w", fticrd = fticrdata, debug =1, compress=True)
        h5f.close()
    #----------------------------------------------
    def test_nparray_to_fticr(self):
        "Test routine that creates a HDF5 file according to a given nparray and returns the buffer (FTICRData)"
        self.announce()
        data_init = 10*np.random.random((2048, 65536)) 
        B = nparray_to_fticrd(self.name_file2, data_init)
        print(B.axis1.size)
    #----------------------------------------------
    def test_axes_update(self):
        "Test routine that overloads the parameter from an axis"
        from . import Apex  as ap
        self.announce()
        import time
        
        d = HDF5File(self.name_file1,"rw",debug = 1)
        d.axes_update(axis = 2, infos = {'highmass':4200.00, 'itype':1})
        d.close()
    #----------------------------------------------
    def _test_chunkshape(self):     # just TOO SLOWWW
        "Test routine that creates , processes Data according to chunkshapes"
        self.announce()
        from time import time
        from .Apex import Import_2D
   
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
        print("---------------SETTINGS---------------------")
        print("CHUNK_CACHE_PREEMPT " , tables.parameters.CHUNK_CACHE_PREEMPT) 
        print("CHUNK_CACHE_SIZE ", tables.parameters.CHUNK_CACHE_SIZE) 
        print("METADATA_CACHE_SIZE ", tables.parameters.METADATA_CACHE_SIZE) 
        print("NODE_CACHE_SLOTS ", tables.parameters.NODE_CACHE_SLOTS)
        print("---------------SETTINGS---------------------")
       #d.hdf5file.close()
       #d = HF.get_data()
        print("import", time()-t1, "secondes")
        t0 = time()
        t00 = t0
        d.rfft(axis=2)
        print("rfft2", time()-t0, "secondes")
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
        print("rfft1", time()-t0, "secondes")
        t0 = time()
        d.modulus()
        print("modulus", time()-t0, "secondes")
        print("calcul", time()-t00, "secondes")
    #----------------------------------------------
    def test_filenodes(self):
        "Test routines that work with filenodes"
        self.announce()
        # 1st objects
        h5f = HDF5File(self.name_file1, "r+", debug =1)
        obj =  ['foo', {'bar': ('baz', None, 1.0, 2)}]  # dummy complex object
        h5f.store_internal_object(obj, h5name='test')
        # then files
        name = 'scan.xml'
        fname = os.path.join(self.DataFolder, name)
        h5f.store_internal_file(fname, h5name=name)
        h5f.flush()
        h5f.close()
        # now retrieve
        h5f = HDF5File(self.name_file1, "r", debug =1)
        objt = h5f.retrieve_object('test')
        self.assertEqual(objt[1]['bar'][2],1.0)
        with self.assertRaises(Exception):
            objt2 = h5f.retrieve_object('foo')
        h5f = HDF5File(self.name_file1, "r")
        F = h5f.open_internal_file(h5name=name)
        for i in range(17):
            l = F.readline()
        print('###',l.strip())
        self.assertEqual(l.strip(),
            b'<scan><count>15</count><minutes>0.5997</minutes><tic>2.37E7</tic><maxpeak>3.618E6</maxpeak></scan>')
#            b"<scan><count>15</count><minutes>0.4828</minutes><tic>1.398E7</tic><maxpeak>3.108E5</maxpeak></scan>")
# depends on self.DataFolder first is for 'ubiquitine_2D_000002.d' second is for 'cytoC_2D_000001.d'
        with self.assertRaises(Exception):
            G = h5f.open_internal_file(h5name="foo.bar")
        F.close()
        h5f.close()
        h5f.close()
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
    print("""
****** ERROR ******
Usage:
{0} Tests               runs self tests
   or
{0} update file_name    to update a old file to current
""".format(prgm))
    sys.exit(1)
#----------------------------------------------
if __name__ == '__main__':
    """ can be used as a prgm for updating files """
    # get
    argv = sys.argv
    try:
        action = argv[1]
    except:
        syntax(argv[0])
        sys.exit(0)
    # do
    if action == 'update':
        try:
            fname = argv[2]
        except:
            syntax(argv[0])
        update(fname, debug = 1)
    elif action == "Tests":
        print("running tests")
        argv.pop()
        unittest.main()
    else:
        syntax(argv[0])
        sys.exit(1)
        