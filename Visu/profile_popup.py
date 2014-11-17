from util.debug_tools import*  
import sys, os, math
from Pyside_PyQt4 import*
from Matplotlib_QtCanvas import*
import numpy as np
from numpy import linspace
from functools import partial

'''
Popup window for managing profiles. 
Permits to save data as csv or pdf. 
'''

@dec_class_pr
@decclassdebugging
class NavigToolbar(NavigationToolbar):
    '''
    Customized navigation toolbar with CSV and PDF added functions.
    '''
    def __init__(self, canvas, parent, data, save, fig, axes, name_profile, namefile):
        self.canvas = canvas
        self.data_profile = data
        self.fig = fig
        self.axes = axes
        self.namefile = namefile 
        self.prep_name_profile = self.namefile + '_' + name_profile
        self.dial = Dialog(data, save, self)                                                        # Instantiates dialog box for saving CSV
        NavigationToolbar.__init__(self, canvas, parent)
        self.clearButtons = []
        self.next = None
        for c in self.findChildren(QToolButton):                                                    # goes through existing buttons.
            if self.next is None:
                self.next = c                                                                                     # Customizing buttons
            if str(c.text()) in ('Save',):                                                          #'Subplots','Customize'
                c.defaultAction().setVisible(False)                                                 # Deactivates visibility of unused buttons
                continue                                                                                 # Need to keep track of pan and zoom buttons
            if str(c.text()) in ('Pan','Zoom'):
                self.clearButtons.append(c)
                self.next = None
        self.custom_button_csv()
        self.custom_button_pdf()
        self.custom_button_fullscale()
    
    def fullscale(self):
        self.axes.set_ylim(0, self.data_profile.buffer.max())
        self.canvas.draw()  
    
    def Button_template(self, add_icon, action):
        '''
        General Button template. 
        Used in custom_button_csv and custom_button_pdf
        '''
        icon = QtGui.QIcon(add_icon)                                                                # embeds icon in QtGui
        act = QAction(action, self)                                                                # creates an action 
        act.setIcon(icon)                                                                          # adds the icon
        act.setToolTip(action)                                                                     # adds the tooltip
        button = QToolButton(self)                                                                  # creates a button
        button.setDefaultAction(act)                                                               # associates an action to the button.
        self.insertWidget(self.next.defaultAction(), button)                                        # insert the button in the toolbar.
        if action.find('save') == 0:
            kind_saved = action[4:]
            print "kind_saved ", kind_saved
            if self.data_profile.axis1.units == 'm/z':
                act.triggered.connect(partial(self.dial.open_file_dialog, kind_saved))                 # Opens window dialog box
            else:
                print "can be saved only in m/z mode."
        else:
            act.triggered.connect(getattr(self, action))
        
    def custom_button_csv(self):
        '''
        Button for saving as CSV. 
        Makes double columns CSV files using NPKv2 csv code.
        '''
        self.Button_template("Visu/iconsUi/CSV.png", 'savecsv')
    
    def custom_button_pdf(self):
        '''
        Button for saving PDFs.
        '''
        self.Button_template("Visu/iconsUi/pdf_profile.png", 'savepdf')
    
    def custom_button_fullscale(self):
        '''
        Button for rescalling at fullscale.
        '''
        self.Button_template("Visu/iconsUi/fullscale.png", 'fullscale')

@dec_class_pr
@decclassdebugging
class PROFILE(QMainWindow): 
    '''
    Popup window for the profiles.
    Creates a QMainWindow in which is made a customized toolbar.
    -Possible to save in PDF, CSV
    -Function for having full scale.
    '''
    def __init__(self, data_profile, save, name_profile, namefile,  ptlx = None):

        QMainWindow.__init__(self)
        self.data_profile = data_profile
        self.name_profile = name_profile
        self.namefile = namefile
        self.ptlx = ptlx                       # used in case of diagonal profile.
        self.save = save                                                         # saving utilities
        self.kind_line = '-'
        self.name_window = 'Profile ' + self.name_profile                                           # name for the popup
        self.window_size = 10                                                                       # size of popup window
        self.prepare_window()                                                                       # prepares canvs, fig etc.. 
        self.make_window()                                                                          # gathers widget, toolbar etc.. 

    def prepare_window(self):
        '''
        Prepares the window with name, size, canvas
        '''
        self.setWindowTitle(self.name_window)
        self.frame = QWidget()
        self.fig = Figure((self.window_size, self.window_size))
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.frame)
        self.axes = self.fig.add_subplot(111)
        
    def make_window(self):
        '''
        Makes the profile window
        '''
        self.toolbar = NavigToolbar(self.canvas, self.frame,\
                    self.data_profile, self.save, self.fig, self.axes,\
                    self.name_profile, self.namefile)                                                 # Instantiates the toolbar.
        vbox = QVBoxLayout()
        vbox.addWidget(self.toolbar)                                                                # adds the toolbar
        vbox.addWidget(self.canvas)                                                                 # adds the canvas
        self.frame.setLayout(vbox)
        self.setCentralWidget(self.frame)
        self.plot_profile()
    
    def axes_mz(self, data):
        '''
        Calculates axes in the case of x or y profile.
        '''
        end, beg = data.mztoi(data.lowmass), data.mztoi(data.highmass)
        return data.itomz( linspace(beg, end, self.data_profile.buffer.size ))

    def plot_profile(self):
        '''
        Plots the profile in the window
        '''
        if debug(self):
            print "self.data_profile.units ", self.data_profile.units
        if self.data_profile.along == 'diag':
            axis_profile = self.data_profile.axes(2).itomz(self.ptlx)                                       # m/z diago                                           
        elif self.data_profile.along == 'x':
            axis_profile = self.axes_mz(self.data_profile.axes(1))                                          # m/z axis x values
        elif self.data_profile.along == 'y':
            axis_profile = self.axes_mz(self.data_profile.axes(2))                                          # m/z axis y values
        if self.data_profile.axis1.units == 'm/z':
            print "makes plot for m/z format "
            self.axes.set_xlabel('m/z')
            self.axes.plot(axis_profile, self.data_profile.buffer, self.kind_line)                          # Plots the data in m/z in the popup window
        elif self.data_profile.axis1.units == 'points':
            self.axes.set_xlabel('points')
            self.axes.plot(self.data_profile.buffer, self.kind_line)                                        # Plots the data in points in the popup window
        self.canvas.draw()                                                                          # print in the canvas.

@dec_class_pr
@decclassdebugging
class Dialog(QtGui.QWidget):
    '''
    Dialog box for saving CSV and PDF files.
    '''
    def __init__(self, data, save, toolbar):
        QtGui.QWidget.__init__(self)
        self.data_profile = data
        self.save = save
        self.toolbar = toolbar
        self.name_dialog = "save profile"                                       # name of the dialog input window
        self.message_input = "enter name profile"                   # message for input.

    def open_file_dialog(self, kind_saved):
        """
        Opens a file dialog 
        Permits to save CSV and PDF.
        """
        qd = QtGui.QInputDialog()
        name_data, result = qd.getText(self, self.name_dialog, self.message_input, QtGui.QLineEdit.Normal, self.toolbar.prep_name_profile)
        name_data = self.save.prep_path_save(str(name_data))
        print "name_data ", name_data
        if kind_saved == 'csv':
            print "saving CSV file"
            self.data_profile.save_csv(name_data + '.csv')
        elif kind_saved == 'pdf':
            self.toolbar.fig.savefig(name_data + '.pdf')

def profile_popup(data):
    '''
    '''
    prof = PROFILE(data)
    prof.show()

if __name__ == "__main__":
    x = [i*0.1 for i in range(1000)]
    y = np.array([math.sin(j) for j in x])
    import FTICR
    data = FTICR.FTICRData(buffer = y)
    data.specwidth = 1667000
    data.ref_mass = 344.0974
    data.ref_freq = 419620.0
    data.highmass = 1000.0
    data.units = "m/z"
    profile_popup(data)
