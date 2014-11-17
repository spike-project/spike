# -*- coding: utf-8 -*-
'''
Created by Lionel Chiron  02/10/2013 
Copyright (c) 2013 __NMRTEC__. All rights reserved.
'''
import os, sys, re, csv
from PySide import QtGui, QtCore
from matplotlib import pyplot as plt
import numpy as np
import pickle as pkl
import glob
import scipy as sp
import scipy.interpolate as spinterp
import scipy.fftpack as fft
from graphic_tools import Subplotnum, Centroid, Label
from mail_error_logger import Logger

'''
Read orbitrap files.
Makes the spectrum in m/z
Methods:
    load
    show
'''

class DATA(object):
    def __init__(self, data, f):
        print "data in class data is  ", data
        self.x, self.y = data
        self.size = os.path.getsize(f)/float(1e6)
        self.max = self.y.max()
        self.type = os.path.splitext(f)[1][1:]

class LOAD_data(object):
    def __init__(self, droplist):
        self.droplist = droplist
        self.data = {} # dictionary containing the data 
        self.data_name = [] # list of the data names
        self.read_dic = {'pkl': self.__read_pkl,
                         'dat':self.__read_dat,
                         'csv':self.__read_csv}
        self.b = 0 # parameter for subplots
        self.load()
    
    def prepare_load(self, i, f):
        base, ext = os.path.splitext(f)
        namef = os.path.basename(base) + ext
       
        return namef, ext
    
    def load(self):
        '''
        Fills the dictionary "self.data" with "self.droplist" contains
        b == 1, only .pkl data
        b == 2, only .dat data
        b == 3, both .pkl and .dat data.. 
        '''
        print "in orbitrap_fticr load droplist"
        print "in orbitrap_fticr reinitialize combobox "
        
        if self.droplist is not None:
            for i, f in enumerate(self.droplist):
                print "file is ", f
                self.extract_data(i, f)
            #self.save_rep_last_file() # save directory of the last file dropped on the interface.
        else:
            self.droplist = None
    
    def extract_data(self, i, f):
        '''
        Extracts data according to extension
        '''
        namef, ext = self.prepare_load(i,f)
        self.kind_subplot(ext)
        #for type_read in self.read_dic:
        self.data[namef] = DATA(self.read_dic[ext[1:]](f), f)
        self.data_name.append(namef)
     
    def kind_subplot(self, ext):
        '''
        Changes parameter self.b to determine the kind of plot
        '''
        if ext == '.pkl':
            if self.b == 0 or self.b == 2 : self.b += 1       
        elif ext == '.dat':
            if self.b == 0 or self.b == 1 : self.b += 2
        elif ext == '.csv':
            print "csv file"
            if self.b == 0 or self.b == 2 : self.b += 1
    
    def __read_dat(self, f):
        """
        reads Orbitrap .dat FID files
        returns a numpy arrays
        """
        print "reading dat data "
        with open(f,'rb') as F:
            pos = 0
            for l in F:
                pos += len(l)
                if l.startswith("Data Points"):
                    print "number of points is",  re.findall(r'\d+', l)[0]
                if l.startswith("Data:"):
                    break
            F.seek(pos)
            y = np.array(np.fromfile(F, dtype = 'f4').tolist()) # ndarray to list to np.array
            x = np.arange(y.size)
        print "x.size, y.size ",x.size, y.size
        return x, y
    
    def __read_csv(self, f):
        '''
        Read csv files
        '''
        print "reading csv data "
        cr = csv.reader(open(f,"rb"))
        absc = []
        ordon = []
        for i, row in enumerate(cr):
            if i > 0 :
                absc.append(np.float(row[0]))
                ordon.append(np.float(row[1]))
        x, y = np.array(absc), np.array(ordon)
        return x, y
    
    def __read_pkl(self, f):
        '''
        Reads .pkl files
        '''
        print "reading pkl data "
        y = pkl.load(open(f,"rb")) # spectrum in point format
        x = np.array(len(y))
        return x, y

class orbitrap_fticr():
    '''
    Orbitrap specific treatments
    '''
    def __init__(self, interface = False,
                 do_centroid = False, do_label = False,
                 xy_lim = None, kind_data = None,
                 selected_spec = 0, data = None,
                 send_error = ''):
        if not interface :
            self.fig = plt.figure(figsize = (16, 4))
        else:
            self.fig = interface.canvas.fig
        self.interface = interface
        self.canvas = interface.canvas
        self.ax = {}
        self.refnorm = 1
        self.refresol = 1
        self.interface.ui.comboBox.clear()
        if data is None :
            self.load()
        else:
            self.use_data(data)
        self.combo_items()
        self.__makesubplot(self.b) # Makes subplot for spectrum and fid
        self.save_rep_last_file() # save directory of the last file dropped on the interface.
        self.selected_spec = selected_spec
        self.keep_select_spec()
        self.mz_ref = 715.3122 # reference mass 
        self.do_centroid = do_centroid # flag for plottig or not the centroids
        self.do_label = do_label # flag for showing the labels after peakpicking.
        #self.do_split = do_split # flag for showing the labels after peakpicking.
        self.trunc_deb = 0.1 # removes low frequencies.
        self.xy_lim = xy_lim # limits for zoom
        self.define_kind_data(kind_data)
        self.thresh_div = 100
        self.threshold = 1e6
        self.send_error = send_error
        if self.send_error == 'sending errors' :
            sys.stderr = Logger()
        
    def load(self):
        l = LOAD_data(self.canvas.droplist)
        self.loading = l
        self.data = l.data
        self.data_name = l.data_name
        self.b = l.b
        if self.canvas.droplist != []:
            self.keep_data = [self.data, self.data_name, self.b]
        else:
            self.keep_data = None
        
    def use_data(self, data):
        print "reusing data "
        self.data = data[0]
        self.data_name = data[1]
        self.b = data[2]
      
    def combo_items(self):
        '''
        Loads the items in the combobox
        '''
        for i, d in enumerate(self.data_name):
            self.load_combobox(i) # add item and icon to combobox.
    
    def max_y_data(self):
        '''
        Find the maximum on y axis in all the data
        '''
        max_y  = 0
        for namef in self.data_name:
            print '###### ', namef
            y_val = self.data[namef].max
            if y_val > max_y:
                max_y = y_val
        return max_y
       
    def keep_select_spec(self):
        '''
        Conserves the selection for the combobox
        '''
        self.interface.ui.comboBox.setCurrentIndex(self.selected_spec)
    
    def define_kind_data(self, kind_data):
        self.kind_data = kind_data
        if self.kind_data is None:
            self.kind_data = 'orbitrap' # orbitrap or fticr
    
    def load_combobox(self, i):
        '''
        Loads the ith icon for combobox
        '''
        self.interface.ui.comboBox.addItem(QtGui.QIcon('images/fig' + str(i) + '.jpg'), '')
    
    def gives_size(self):
        '''
        Gives the size of the file in Mb
        '''
        try:
            self.interface.ui.label_3.setText("size is " +
                str(int(self.current_data.size)) + " Mb")
        except :
            pass
        
    @property
    def current_data(self):
        return self.data[self.data_name[self.selected_spec]]
        
    def set_combo_last_index(self):
        '''
        Position index combo to last element in droplist
        '''
        self.interface.ui.comboBox.setCurrentIndex(len(self.canvas.droplist)-1)
    
    def save_rep_last_file(self):
        '''
        Save the directory of the last file dropped on the interface
        '''
        try:
            print "#### save last file directory "
            if len(self.interface.canvas.droplist) > 0 :
                rep = open('last_rep.txt', 'w')
                rep.write(os.path.dirname(self.interface.canvas.droplist[-1]))
                rep.close()
        except Exception:
            pass
         
    def mixt_subplot(self):
        '''
        Subplot spectrum and data
        '''
        xlabeldat = 'time in s'; ylabeldat = 'Intensity'
        xlabelpkl = 'm/z'; ylabelpkl = ''
        sp = Subplotnum(2, 1, self.fig) # Two subplots
        self.ax['spec'] = sp.next(xlabel = xlabelpkl, ylabel = ylabelpkl)
        
        self.ax['dat'] = sp.next(xlabel = xlabeldat, ylabel = ylabeldat)
    
    def one_kind_subplot(self, b):
        '''
        Subplot data or spectrum only
        '''
        xlabeldat = 'time in s'; ylabeldat = 'Intensity'
        xlabelpkl = 'm/z'; ylabelpkl = ''
        sp = Subplotnum(1, 1, self.fig)# only one kind of data, one subplot
        if b == 1: # only .pkl
            self.ax['spec'] = sp.next(xlabel = xlabelpkl, ylabel = ylabelpkl) # creates only spectrum
            self.gca = self.ax['spec']
            self.toolb_mode = plt.get_current_fig_manager().toolbar.mode
            print "get gca #################", self.gca.get_xlim()
        elif b == 2: # only .dat
            self.ax = {'dat': sp.next(xlabel = xlabeldat, ylabel = ylabeldat)} # creates only fid
    
    def __makesubplot(self, b):
        '''
        Subplot for spectrum and fid
        b == 1, only .pkl data
        b == 2, only .dat data
        b == 3, both .pkl and .dat data.. 
        '''
        if b == 3 : # two kinds of data: Fids and spectra
            self.mixt_subplot()
        else : # case b == 1 or b == 2
            self.one_kind_subplot(b) 

    def __prep_orbit_xpec_spec(self, namef):
        '''
        Prepare data for orbitrap m/z unit. Called in plot_spec
        '''
        if self.kind_data == 'orbitrap' :
            data_size = self.data[namef].y.size
            point_ref = data_size * 377500./(2*999982)
            trunc = self.trunc_deb * data_size
            xaxis = self.mz_ref/(((1 + np.arange(1.0*data_size))/point_ref)**2)
            xspec = xaxis[trunc:]
            spec =  self.data[namef].y[trunc:]
        elif self.kind_data == 'fticr' :
            xspec, spec = self.data[namef].x, self.data[namef].y # spectrum
        return xspec, spec
    
    def makes_centro_labels(self, namef, xspec, spec):
        if self.data_name.index(namef) == self.selected_spec :   # if spectrum is selected with combobox..
            print "name of the spectrum is ", namef
            if self.do_centroid : self.centroid.plot(xspec, spec) # plot centroids
            if self.do_label :
                self.label.plot(xspec, spec) # plot labels
    
    def keep_labels(self):
        self.label_vec = self.label.label_vec
        self.peaks_vec = self.label.peaks_vec
    
    def replot_labels_zoom(self):
        
        self.label = Label(self.ax, 'spec', self.threshold) # instantiates label
        self.label.label_vec = self.label_vec
        self.label.peaks_vec = self.peaks_vec
     
        print "####### ########## #############  in replot_labels_zoom"
        print '################', self.label.label_vec
        print self.label.peaks_vec
        deltax = abs(self.xy_lim[0][0]-self.xy_lim[0][1])/4
        diff = (self.label.label_vec-self.label.peaks_vec)
        new_arr = diff/abs(diff)*deltax
        self.label.label_vec = self.label.peaks_vec + new_arr
        self.label.label_above_arrow() # avoid having labels behind arrows
        self.label.peaklabel()
    
    def current_zoom(self):
        print "in current zoom "
        print "self.xy_lim ", self.xy_lim
        self.ax['spec'].set_xlim(self.xy_lim[0])
        self.ax['spec'].set_ylim(self.xy_lim[1])
    
    def __plot_spec(self, namef):
        '''
        Subplots prepared in makesubplot
        plot the data, the  centroids and the labels
        called by function show
        '''
        xspec, spec = self.__prep_orbit_xpec_spec(namef) # extract x and y from data
        print "namef ", namef
        self.ax['spec'].plot(xspec, spec, label = namef)# plot spectrum
        self.gives_size() # give the size of the file
        self.ax['spec'].set_ylim(0, 1.2*self.max_y_data()) # height of the figure.

        try:
            self.current_zoom()# keep the view in the current zoom values.. 
        except Exception:
            pass

        self.makes_centro_labels(namef, xspec, spec) # makes centroids and labels.
        
    def __plot_fid(self, namef):
        '''
        plot the Fid
        '''
        self.ax['dat'].plot(self.data[namef].y, label = namef)
    
    def defines_thresh(self):
        '''
        Threshold for centroids and labels
        '''
        try:
            self.threshold = int(float(self.interface.ui.lineEdit.text())*1e6)# retrieves threshold
        except Exception:
            pass
    
    def show(self):
        '''
        self.data contains .pkl or/and .dat .
        print spectrum in m/z
        truncation for avoiding high spurious m/z
        '''
        if self.data != {}:
            print "self.data != {}"
            print "self.data"
            self.defines_thresh() # takes thres from lineEdit or from initialization.
            if self.do_centroid : self.centroid = Centroid(self.ax,'spec', self.threshold )# instantiates centroid
            if self.do_label : self.label = Label(self.ax, 'spec', self.threshold) # instantiates label
            for namef in self.data_name: # calls the files in the sorted maneer.
                print "namef ", namef
                f, ext = os.path.splitext(namef)
                if ext in ['.pkl', '.csv'] :
                    print "plot m/z "
                    self.__plot_spec(namef) # plot spectrum
                elif ext == '.dat':
                    print "plot fid "
                    self.__plot_fid(namef) # plot fid
            for key in self.ax.iterkeys():
                self.ax[key].legend()
            print "######### end of orbitrap show"