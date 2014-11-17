class PARAM_ZOOM(object):
    '''
    Central object for zoom management. 
    '''
    def __init__(self, data):
        self.data = data                                                                    # 
        self.data_resmin = self.data.d[len(self.data.d)-1]                                  # Minimal resolution
        self.zoom_coord = [0, 0, self.data_resmin.size2, self.data_resmin.size1]            # initial zoom
        self.zoom_coord_old = self.zoom_coord                                               #
        self.movezoo = [0, 0, 0, 0]                                                         #
        self.greyzoo = None                                                                 # flag for displaying the zoom with greyzone
        self.zoomready = False                                                             # flag for zooming
        self.listview_index = 0                                                             # number of coordinates registered.
        self.sliderpos = 1                                                                  # initial position of the slider for the scale
        self.scale = 1                                                                      # scale for plotting the 2D views
        self.listview = [{'zoom' : self.zoom_coord, 'resolution' : 'resol1',\
            'scale' : self.scale, 'slider' : self.sliderpos}]                       # first zoom, resolution and scale
        self.mouse_motion = None
        self.area3Dmax = 5e3        # maximum area for 3D zoom
    
    def zoom_diag_vector(self):
        '''
        Vector from diagonal for drag etc.. 
        '''
        vx = self.zoom_coord[0] - self.zoom_coord[2]  # displacement in x
        vy = self.zoom_coord[1] - self.zoom_coord[3] # displacement in y
        return vx, vy
    
    def report(self):
        print "############ report paramz"
        for elem in dir(self):
            if elem.find('_') != 0:
                if not hasattr(getattr(self, elem), '__call__') :
                    print elem, getattr(self, elem)
        print "############ "