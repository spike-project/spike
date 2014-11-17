import  NPKConfigParser as npkcfg 
from File.HDF5File import HDF5File
from util.debug_tools import*   

@dec_class_pr
@decclassdebugging
class Proc_Parameters(object):
    """this class is a container for visualization parameters"""
    def __init__(self, configfile = None):
        "initialisation, see visu2D.mscf for comments on parameters"
        self.multresfile = None
        if configfile:
            self.load(configfile)
    def load(self, cp):
        "load from cp config file - should have been opened with ConfigParser() first"
        self.multresfile =  cp.get( "visu2d", "multresfile") # multiresolution file

    def report(self):
        "print a formatted report"
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
    def __init__(self):
        '''
        Initialize , recuperate addresses from visu2D.mscf and create multiresolution file 
        '''
        print "in LOAD"
        self.d = [] # list that will contain all the resolutions. 
        self.d_interm = [] #  list that will contain all the resolutions only in point
        self.mode_point = True # mode point
        self.active = True
        cp = npkcfg.NPKConfigParser()
        cp.readfp(open("Visu/visu2D.mscf"))
        self.param = Proc_Parameters(cp) # parameters from config file.. 
        self.param.report()
        
        # rootfiles = os.getcwd()# address of the working directory
        # config = ConfigParser.ConfigParser()
        # config.read(rootfiles+"/visu2D.mscf")# reading the config file visu2D.mscf
        # self.param.multresfile = config.get('visu2d', "multresfile")# address of the multiresolution file
        
        self.change_resolution(u = "resol1"# initialising to the resolution 1
        self.loadres()
     
    def loadres(self):
        '''
        Load the different resolutions from Hdf5 files
        and put the in a list of FTICRdata objects.
        This list is named after res.d
        Count also the number of resolutions
        ''' 
        self.NBRES = 0 # number of resolutions
        f = HDF5File(self.param.multresfile, 'r')
        for gp in f.hf.walkGroups("/"):     # all groups in file
            for nd in gp._f_listNodes():    # all nodes in group
                if nd._v_name == 'data':    # if node is named "data" it is a dataset
                    self.d.append( f.get_data(gp._v_name) )     # resolutions for C window
                    self.d_interm.append( f.get_data(gp._v_name) )     # resolutions for D window
                    self.NBRES += 1
        self.resmin = self.d_interm[len(self.d_interm)-1] # resolution for D window
        print "all the resolutions are loaded for both windows"

if __name__ == '__main__':
    l = LOAD()
    print "type(l.d[0])", type(l.d[0])