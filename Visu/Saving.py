from util.debug_tools import* 
import time
import os  

@dec_class_pr
@decclassdebugging
class SAVE(object):
    '''
    Class for savings the 2D, 3D, profiles etc.. 
    Permits to save the object in the same location : a directory named with the date etc.. 
    '''
    def __init__(self, data):
        self.dir_saved = None
        self.data = data
        
    def dir_save(self):
        '''
        Makes a directory for the data if it doesn't exist
        The directory is named with the day, the month, the year and the hour. 
        '''
        data_dir = os.path.dirname(self.data.param.multresfile)                                 # taking root dir from memorized path.
        date = time.strftime('%d-%m-%Y-%Hh', time.localtime())                                  # Makes the date. 
        self.dir_saved = os.path.join(data_dir, "data_saved_" + date)
        if not os.path.exists(self.dir_saved):                                                  # Checks if directory yet exists.
            os.mkdir(self.dir_saved)

    def prep_path_save(self, namefile):
        '''
        prepares address path_save.
        '''
        self.dir_save()                                                                         # Makes self.dir_saved                                                          # retrieves file 's name from QlineEdit
        path_save = os.path.join(self.dir_saved, namefile)                                      # address where to save the data
        print "data saved at ", path_save                                                       # tells where it is saved.
        return path_save
     
if __name__ == '__main__':
    pass