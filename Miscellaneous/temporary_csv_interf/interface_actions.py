# -*- coding: utf-8 -*-
'''
Created by Lionel Chiron  02/10/2013 
Copyright (c) 2013 __NMRTEC__. All rights reserved.
'''
import sys, os
from matplotlib import pyplot as plt
from PySide import QtGui, QtCore
from PySide.QtCore import QObject
from functools import partial
import os.path as op
from time import sleep

class interact(object):
    '''
    Makes the link between interface and orbitrap library.
    Functions are :
        clear
        copy to clipboard
        make_ctx_menu
        centroid
        label
    '''

    def __init__(self, interface, orbitrap_fticr):
        self.interface = interface
        self.orbitrap_fticr = orbitrap_fticr
        self.xy_lim = None
        self.orig_view = False
        self.dial = Dialog(interface, self)
        self.refresh_plot(fig = False)
    
    def swap(self, a, value1, value2, button):
        print "swapping between orbitrap and fticr"
        if a == value1 :
            a = value2
            print 'passing from ', value1, ' to ', value2
        elif a == value2 :
            a = value1
            print 'passing from ', value2, ' to ', value1
        button.setText(a)
        return a
        
    def keep_orig_view(self):
        '''
        Saving the first view
        '''
        if self.orig_view == False:
            self.save_xy_lim()
            self.xy_lim_orig = self.xy_lim
            print "in keep_orig_view self.xy_lim ", self.xy_lim
            self.orig_view = True
        
    def make_connections(self):
        '''
        Makes connections to buttons, menus, drop signal etc.. 
        '''
        #button_release_event
        self.click = self.interface.canvas.mpl_connect('button_press_event', self.onclick)
        self.release = self.interface.canvas.mpl_connect('button_release_event', self.onrelease)
        self.drop = self.interface.canvas.connect(QtCore.SIGNAL('dropped'), self.drop_refresh)
        #self.interface.ui.pushButton_5.clicked.connect(self.clear)
        #self.interface.ui.toolButton.clicked.connect(self.dial.open_file_dialog)
        self.interface.ui.pushButton.clicked.connect(self.toggle_orbit_fticr)
        self.interface.ui.comboBox.connect(QtCore.SIGNAL('currentIndexChanged(int)'), self.read_combobox)
    
    #def connect_menubar_file(self):
    #    '''
    #    Menu File
    #    '''
    #    self.connect_menubar_export()

    def connect_menubar_edit(self):
        '''
        Menu Edit
        '''
        self.connect_menubar_clipboard()
        self.connect_menubar_clear()
        self.connect_menubar_open_last_rep()
     
    def connect_menubar_tool(self):
        '''
        Menu Tool
        '''
        self.connect_menubar_peakpick()
        self.connect_menubar_centroid()
        self.connect_menubar_change_size()
    
    def connect_menubar_navigation(self):
        '''
        Menu Navigation
        '''
        self.connect_menubar_toolbar()
        self.connect_menubar_origview()
    
    def connect_menubar_help(self):
        self.connect_menubar_send_error()
    ############################################## Connects the menus
        
    def make_menubar(self):
        '''
        Makes the whole menubar
        '''
        #self.connect_menubar_file()
        self.connect_menubar_edit()
        self.connect_menubar_tool()
        self.connect_menubar_navigation()
        self.connect_menubar_help()
    
    ############################################## In each menu
    ##################################### Menu File
    def connect_menubar_export(self):
        '''
        Menubar File, copies to clipboard
        '''
        a = QtGui.QAction('&export centroids to mzML', )
        a.setShortcut(QtGui.QKeySequence('Ctrl+E'))
        a.setStatusTip('Export in xmzML format')
        a.triggered.connect()
        if not hasattr(self,'mzml'):
            self.mzml = self.interface.menubar_file.addAction(a)
    
    ##################################### Menu Edit
    def connect_menubar_clipboard(self):
        '''
        Menubar Edit, copies to clipboard
        '''
        a = QtGui.QAction('&copier', self.interface.menubar_edit)
        a.setShortcut(QtGui.QKeySequence('Ctrl+C'))
        a.setStatusTip('Copy to clipboard')
        a.triggered.connect(self.interface.copy_clipboard)
        if not hasattr(self,'copy'):
            self.copy = self.interface.menubar_edit.addAction(a)
            
    def connect_menubar_clear(self):
        '''
        Menubar Edit , clears last spectrum of fid.
        '''
        a = QtGui.QAction('&clear', self.interface.menubar_edit)
        a.setShortcut(QtGui.QKeySequence('Ctrl+D'))
        a.setStatusTip('Clear')
        a.triggered.connect(self.clear)
        if not hasattr(self,'clear_data'):
            self.clear_data = self.interface.menubar_edit.addAction(a)
    
    def connect_menubar_open_last_rep(self):
        '''
        Menubar Edit, opens the last directory used.
        '''
        a = QtGui.QAction('&open last directory', self.interface.menubar_edit)
        a.setShortcut(QtGui.QKeySequence('Ctrl+R'))
        a.setStatusTip('open last directory')
        a.triggered.connect(self.dial.open_file_dialog)
        if not hasattr(self,'open_last_dir'):
            self.open_last_dir = self.interface.menubar_edit.addAction(a)
    
    ##################################### Menu Tools
    
    def connect_menubar_change_size(self):
        '''
        Menubar Tools, changes the size of the spectrum.
        '''
        a = QtGui.QAction('&Change_size', self.interface.menubar_tools) #
        a.setShortcut(QtGui.QKeySequence('Ctrl+S'))
        a.setStatusTip('Change size')
        a.triggered.connect(self.redefines_max_y)
        if not hasattr(self,'cs'):
            self.cs = self.interface.menubar_tools.addAction(a)
        
    def connect_menubar_peakpick(self):
        '''
        Menubar Tools, makes the peakpicking
        '''
        a = QtGui.QAction('&Peakpicking', self.interface.menubar_tools) #
        a.setIcon(QtGui.QIcon('images/peakpicking.png'))
      
        a.setShortcut(QtGui.QKeySequence('Ctrl+Q'))
        a.setStatusTip('Peakpicking')
        a.triggered.connect(self.label)
        if not hasattr(self,'pp'):
            self.pp = self.interface.menubar_tools.addAction(a)
            
    def connect_menubar_centroid(self):
        '''
        Manubar Tools, makes the centroids
        '''
        a = QtGui.QAction('&resolutions', self.interface.menubar_tools)
        a.setIcon(QtGui.QIcon('images/resolution.png'))
        a.setShortcut(QtGui.QKeySequence('Ctrl+W'))
        a.setStatusTip('Resolutions')
        a.triggered.connect(self.centroid)
        if not hasattr(self,'resol'):
            self.resol = self.interface.menubar_tools.addAction(a)
      
    def connect_menubar_toolbar(self):
        '''
        Show or not the toolbar
        '''
        a = QtGui.QAction('&toolbar', self.interface.menubar_navigation)
        a.setShortcut(QtGui.QKeySequence('Ctrl+Z'))
        a.setStatusTip('show toolbar')
        a.triggered.connect(self.toolbar)
        if not hasattr(self,'tb'):
            self.tb = self.interface.menubar_navigation.addAction(a)
            
    def connect_menubar_origview(self):
        '''
        Show the first view
        '''
        a = QtGui.QAction('&first view', self.interface.menubar_navigation)
        a.setShortcut(QtGui.QKeySequence('Ctrl+F'))
        a.setStatusTip('show first view')
        a.triggered.connect(self.original_view)
        if not hasattr(self,'origview'):
            self.origview = self.interface.menubar_navigation.addAction(a)
    
    
    ##################################### Menu Help  
    def connect_menubar_send_error(self):
        '''
        Menubar Help, selects error message sending or deselects.
        '''
        a = QtGui.QAction('&send mail error message', self.interface.menubar_help)
        a.setShortcut(QtGui.QKeySequence('Ctrl+M'))
        a.setStatusTip('sendind error message')
        a.triggered.connect(self.send_error)
        if not hasattr(self,'send_err'):
            self.send_err = self.interface.menubar_help.addAction(a)
    
    #################################################
    
    def connect_clipboard(self):
        '''
        Menubar copy paste.
        '''
        a = QtGui.QAction('copier', self.interface.menu)
        a.triggered.connect(self.interface.copy_clipboard)
        self.interface.menu.addAction(a)

    def connect_split_spec(self):
        a = QtGui.QAction('split_spec', self.interface.menu)
        a.triggered.connect(self.split_spectra)
        self.interface.menu.addAction(a)
        
    def connect_remove_last_spec(self):
        '''
        Removes last spectrum or raw data
        '''
        a = QtGui.QAction('remove_last_spec', self.interface.menu)
        a.triggered.connect(self.remove_last_spec)
        self.interface.menu.addAction(a)
 
    def __actions_ctx_menu(self):
        '''
        Actions of the context menu
        '''
        self.connect_clipboard()
        self.connect_split_spec()
        self.connect_remove_last_spec()

    def __actions_ctx_submenu(self):
        print "makes submenu"
        self.interface.submenu = QtGui.QMenu(self.interface.menu)
        self.interface.submenu.setTitle("Submenu")
        self.interface.menu.addMenu(self.interface.submenu)   
    
    def make_ctx_menu(self):
        '''
        Makes context menu and submenus.
        '''
        self.__actions_ctx_menu()
        #self.actions_ctx_submenu()
    
    def make_connections_and_ctxmenu(self):
        '''
        Applying connexion and context menu
        '''
        self.make_connections()
        self.make_ctx_menu()
        self.make_menubar()
        #self.define_buttons()
        
    def save_xy_lim(self):
        '''
        saves figure coordinates when releasing the mouse
        '''
        try:
            print "save coordinates after release "
            xlim = self.orbitrap_fticr.gca.get_xlim()
            ylim = self.orbitrap_fticr.gca.get_ylim()
            self.xy_lim = [xlim, ylim]
            print "############  coord are ", self.xy_lim
        except :
            self.xy_lim = None
            
    def init_orbitrap_fticr(self):
        '''
        Reinitializes the part orbitrap_fticr
        called by refresh_plot
        '''
        self.orbitrap_fticr.__init__(interface = self.interface,
                       do_centroid = self.orbitrap_fticr.do_centroid,
                       do_label = self.orbitrap_fticr.do_label,
                       xy_lim = self.xy_lim,
                       kind_data = self.orbitrap_fticr.kind_data,
                       selected_spec = self.orbitrap_fticr.selected_spec,
                       data = self.orbitrap_fticr.keep_data,
                       send_error = self.orbitrap_fticr.send_error
                       )    
   
    def refresh_plot(self, label = False, fig = False):
        '''
        Refresh plot
        reinits interface, reinits orbitrap_fticr and makes the connexions
        '''
        print "refresh_plot"
        droplist = self.interface.canvas.droplist
        self.interface.init_interf() #reinitialize interface
        if fig :
            fig = self.interface.canvas.fig #conserving fig
            self.interface.canvas.droplist = droplist #conserving droplist
        self.init_orbitrap_fticr()
        print "####################### after init_orbitrap self.xy_lim ", self.xy_lim
        self.orbitrap_fticr.show()
        self.keep_orig_view() # keep the first view
        self.make_connections_and_ctxmenu()
        
    def on_motion(self, event):
        '''
     
        '''
        dx = event.xdata
        dy = event.ydata
    
    def onclick(self, event):

        self.xclick = event.xdata
        self.yclick = event.ydata
    
    def sorting_xy(self):
        x_int = [self.xclick, self.xrelease]
        y_int = [self.yclick, self.yrelease]
        x_int.sort()
        y_int.sort()
        self.xy_lim = [x_int, y_int]
        print "after sorting ", self.xy_lim
    
    def onrelease(self, event):
        '''
        At each release of the mouse, saves the limits
        '''
        mode = self.interface.navi_toolbar.mode
        print "mode is ", mode
        if mode == "zoom rect" :
            print "in zoom rect mode "
            self.xrelease = event.xdata
            self.yrelease = event.ydata
            self.sorting_xy()
            
            if self.orbitrap_fticr.do_label:
                print "######################## zooming"
                self.orbitrap_fticr.do_label = False
                self.orbitrap_fticr.keep_labels()
                self.refresh_plot(fig = True)
                self.orbitrap_fticr.replot_labels_zoom()
                self.orbitrap_fticr.do_label = True
    
    def select_format(self):
    
        '''
        Position name of data format.
        '''
        base, ext = op.splitext(self.interface.canvas.droplist[-1])# read extension last file
        cnd1 = self.orbitrap_fticr.kind_data == 'orbitrap' and ext == '.csv'
        cnd2 = self.orbitrap_fticr.kind_data == 'fticr' and ext == '.pkl' or ext == 'dat'
        if cnd1 or cnd2:
            self.toggle_orbit_fticr()
    
    def add_new_data(self):
        droplist = self.interface.canvas.droplist
        last_drop = droplist[-1]
        pos = len(droplist)-1
        self.orbitrap_fticr.loading.extract_data(pos, last_drop)

    def drop_refresh(self):
        '''
        Executed when adding a new file for example.
        '''
        self.select_format() # indicate the data format.
        self.add_new_data() # add last element of droplist.. 
        self.refresh_plot(fig = True)
        self.orbitrap_fticr.gives_size() # gives the size of data
        self.orbitrap_fticr.set_combo_last_index() # position combo to last index
    
    def clear_all(self):
        '''
        Clear plot
        '''
        self.orbitrap_fticr.keep_data = None
        self.refresh_plot(fig = False)
    
    def clear(self): #_current
        '''
        remove the current curve given by self.orbitrap_fticr.selected_spec
        '''
        self.interface.canvas.droplist.pop(self.orbitrap_fticr.selected_spec)
        print "self.orbitrap_fticr.selected_spec ", self.orbitrap_fticr.selected_spec
        print "self.interface.canvas.droplist ", self.interface.canvas.droplist        
        self.orbitrap_fticr.keep_data = None
        self.refresh_plot(fig = True)
        
    def centroid(self):
        '''
        plot centroids in spectra
        '''
        self.interface.ui.label_4.setText('find resolutions')
        self.orbitrap_fticr.do_centroid = not(self.orbitrap_fticr.do_centroid)
        self.refresh_plot(fig = True)
        
    def label(self):
        '''
        plot labels
        '''
        self.interface.ui.label_4.setText('peakpicking')
        self.orbitrap_fticr.do_label = not(self.orbitrap_fticr.do_label)
        self.refresh_plot(fig = True)
        
    def original_view(self):
        self.xy_lim  = self.xy_lim_orig
        self.refresh_plot(fig = True)
    
    def send_error(self):
        a = self.orbitrap_fticr.send_error
        a = self.swap(a, '', 'sending errors', self.interface.ui.label)
        self.orbitrap_fticr.send_error = a  
    
    def split_spectra(self):
        self.xy_lim  = self.xy_lim_orig
        self.refresh_plot(fig = True)
    
    def remove_last_spec(self):
        '''
        Removes last element in droplist and refreshes
        '''
        self.interface.canvas.droplist.pop()
        self.refresh_plot(fig = True)
    
    def toolbar(self):
        '''
        show toolbar or not
        '''
        self.interface.show_toolbar = not(self.interface.show_toolbar)
        self.refresh_plot(fig = True)
    
    def toggle_orbit_fticr(self):
        '''
        Passing from orbitrap to fticr and vice versa
        '''
        a = self.orbitrap_fticr.kind_data
        button = self.interface.ui.pushButton
        a = self.swap(a, 'orbitrap', 'fticr', button)
        self.orbitrap_fticr.kind_data = a
           
    def read_combobox(self):
        '''
        Takes the index of combobox
        '''
        index_combo = self.interface.ui.comboBox.currentIndex()
        print "index combobox is ", index_combo
        self.orbitrap_fticr.selected_spec = index_combo
        print "in read combo self.orbitrap_fticr.selected_spec is ",self.orbitrap_fticr.selected_spec
        #self.interface.ui.label.setText("index combo " + str(index_combo))
        self.orbitrap_fticr.gives_size() # gives size of selected data..
        self.orbitrap_fticr.do_label = False
        self.orbitrap_fticr.do_centroid = False
        
    def redefines_max_y(self):
        '''
        changes the spectrum size 
        '''
        print "changes size of data"
        self.target_max = float(self.interface.ui.lineEdit_2.text())*1e6# retrieves max spectrum.
        print "targeted size is ", self.target_max
        print "self.orbitrap_fticr.current_data.max ", self.orbitrap_fticr.current_data.max
        ratio = self.target_max/float(self.orbitrap_fticr.current_data.max)
        print "changement is done with ratio ", ratio
        self.orbitrap_fticr.current_data.y *= ratio
        self.refresh_plot(fig = True)
        self.orbitrap_fticr.current_data.max = self.target_max
       
class Dialog(QtGui.QWidget):

    def __init__(self, interface, interact):
        QtGui.QWidget.__init__(self)
        self.interface = interface
        self.interact = interact
 
    #----------------------------------------------------------------------
    def open_file_dialog(self):
        """
        Opens a file dialog and sets the label to the chosen path
        """
        qd = QtGui.QFileDialog
        f = open('last_rep.txt', 'r')
        addr = f.readline()
        print "last directory is ", addr
        f.close()
        self.path, _ = qd.getOpenFileName(self, "Open File", addr) #, None, None,  QtGui.QFileDialog.
        if self.path != '':
            self.interface.canvas.droplist.append(self.path)
            self.interact.drop_refresh()
        
    def about(self):
        '''Popup a box with about message.'''
        QMessageBox.about(self, "About PyQt, Platform and the like",
                """<b> About this program </b> v %s
                #<p>Copyright  2013 Lionel Chiron. 
                #All rights reserved in accordance with
                ## - NO WARRANTIES!
                ##<p>This application can be used for
                ##displaying OS and platform details.
                ###<p>Python %s -  PySide version %s - Qt version %s on %s"""  % 
                (__version__, platform.python_version(), PySide.__version__,
                 PySide.QtCore.__version__, platform.system()))
        