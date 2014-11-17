import  NPKConfigParser as npkcfg 
from File.HDF5File import HDF5File
from util.debug_tools import*
import os

@dec_class_pr
@decclassdebugging
class Visu_Parameters(object):
    """
    this class is a container for visualization parameters
    """
    def __init__(self, configfile = None):
        '''
        '''
        self.multresfile = None
        if configfile:
            self.load(configfile)       # Loads data from informations in the configfile
            
    def load(self, cp):
        '''
        Loads in self the information from the configfile.
        '''
        self.multresfile =  cp.get( "visu2d", "multresfile")            # multiresolution file
        print "name of the msh5 to be loaded is  ", self.multresfile
        
    def report(self):
        '''
        Show all the parameters in self for the class Visu_parameters.
        '''
        print "------------ processing parameters ------------------"
        for i in dir(self):
            if not i.startswith('_'):
                v = getattr(self,i)
                if not callable(v):
                    print i, ' :', v
        print "-----------------------------------------------------"

@dec_class_pr
@decclassdebugging
class LOAD(object):
    '''
    Class to load the resolutions from the visu2D.mscf file
    '''
    def __init__(self, configfile = None, msh5file = None):
        '''
        Retrieves addresse of the msh5 from visu2D.mscf and creates
        a list containing resolutions as FTICRData files.
        '''
        print "in LOAD"
        self.d = []         # list that will contain all the resolutions for window C
        self.d_interm = []      #  list that will contain all the resolutions only in point for window D
        self.mode_point = True      # mode point
        self.active = True
        #### with configfile
        if configfile:
            cp = npkcfg.NPKConfigParser()  # Instantiate Config Parser.
            cp.readfp(open(configfile)) 
            self.param = Visu_Parameters(cp) # parameters from config file.. 
            self.param.report()             # check configfile
            
        elif msh5file:
            self.param = Visu_Parameters()
            self.param.multresfile =  msh5file
        self.namefile = os.path.basename(self.param.multresfile)[:-4]
        ###
        self.resolu = "resol1"              # initialising to the resolution 1
        self.loadres()              # Loading all the resolutions
     
    def loadres(self):
        '''
        Load the different resolutions from Hdf5 files
        and put the in a list of FTICRdata objects.
        This list is in self.d
        self.d[0]
        Count also the number of resolutions
        ''' 
        self.NBRES = 0          # number of resolutions
        f = HDF5File(self.param.multresfile, 'r')           # Reads the Hdf5 file
        for gp in f.hf.walkGroups("/"):     # all groups in file
            for nd in gp._f_listNodes():    # all nodes in group
                if nd._v_name == 'data':    # if node is named "data" it is a dataset
                    self.d.append( f.get_data(gp._v_name) )     # resolutions for C window
                    self.d_interm.append( f.get_data(gp._v_name) )     # resolutions for D window
                    self.NBRES += 1
        self.resmin = self.d_interm[len(self.d_interm)-1]           # resolution for D window
        print "all the resolutions are loaded for both windows"

if __name__ == '__main__':
    l = LOAD()
    print "type(l.d[0])", type(l.d[0])
    print dir(l)