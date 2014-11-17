from util.debug_tools import*   

@dec_class_pr
@decclassdebugging
class CONVERT(): 
    """
    Conversion between m/z and points. 
    """ 
    def __init__(self, display, data, paramz):
        self.display = display
        self.data = data
        self.paramz = paramz

    def mztoi(self, d, coorr, ax):
        '''
        mz to point for each coordinate
        '''
        #newc = int(d.axes(ax).mztoi(coorr))
        newc = d.axes(ax).mztoi(coorr)
        return newc

    def mztoi_all(self, d, llx, lly, urx, ury):
        '''
        transforming from "m/z" coordinate to "point" coordinates of coorr (zoom window)
        '''
        llx, lly = d.axes(2).mztoi(llx), d.axes(1).mztoi(lly)
        urx, ury = d.axes(2).mztoi(urx), d.axes(1).mztoi(ury)
        return llx, lly, urx, ury

    def itomz_all(self, d, llx, lly, urx, ury):
        '''
        transforming from "point" coordinates to "mz" coordinates.
        '''
        if debug(self):
            print "d.axes(2).itomz(llx) ", d.axes(2).itomz(llx)
            print "d.axes(2).highmass ", d.axes(2).highmass
        llx, lly = min(d.axes(2).itomz(llx), d.axes(2).highmass), min(d.axes(1).itomz(lly), d.axes(1).highmass)
        urx, ury = max(d.axes(2).itomz(urx), d.axes(2).lowmass),  max(d.axes(1).itomz(ury), d.axes(1).lowmass)
        return llx, lly, urx, ury

    def pass_to_pt(self, llx, lly, urx, ury):
        '''
        transforming from point coordinate to m/z coordinates with mode point or m/z conditon
        '''
        if debug(self):
            print "in convert.pass_to_pt "
            print "llx, lly, urx, ury ", llx, lly, urx, ury
            print " self.data.mode_point ", self.data.mode_point
        dd = self.display.res2dd()                                                             # retrieves fticrdata
        if not self.data.mode_point :                                                       # if m/z
            llx, lly, urx, ury = self.mztoi_all(dd, llx, lly, urx, ury)
            print "after self.mztoi_all llx, lly, urx, ury are ", llx, lly, urx, ury
        return llx, lly, urx, ury
    
        
    def pass_to_curr_mode(self, llx, lly, urx, ury):
        '''
        If in mz/mode pass coordinates in "m/z mode", if in point mode pass the coordinates in "point mode".
        '''
        if debug(self):
            print "self.data.mode_point ", self.data.mode_point
            print "in graphic_tools.convert.pass_to_curr_mode"
            print "coordinates llx, lly, urx, ury = ", llx, lly, urx, ury
        if not self.data.mode_point :                                                       # if m/z
            dd = self.display.res2dd()                                                      # retrieves fticrdata
            llx, lly, urx, ury = self.itomz_all(dd, llx, lly, urx, ury)                    # calculating new zoom window for m/z view
        return llx, lly, urx, ury
    
    def maxres(self,llx, lly, urx, ury):
        '''
        If mode point, converts to coordinates with maximal resolution
        '''
        if self.data.mode_point:
            llx, lly, urx, ury = self.itomz_all(self.display.currentd, llx, lly, urx, ury)
            llx, lly, urx, ury = self.mztoi_all(self.display.RESMAX, llx, lly, urx, ury)
        return llx, lly, urx, ury
            
    def to_npk(self):
        "return npk formated zoom"
        zc = self.paramz.zoom_coord
        self.set(zc[0], zc[1], zc[2], zc[3])
        zc = self.paramz.zoom_coord
        zz = [[zc[1], zc[3]], [zc[0], zc[2]]]
        return zz
    
    def set(self, xa, ya, xb, yb):
        '''
        Setter to put the coordinate in the right order and avoid having values outside.
        '''
        #gives the x0,y0,x1,y1 in good order
        x0 = min(xa, xb)
        y0 = min(ya, yb)
        x1 = max(xa, xb)
        y1 = max(ya, yb)
        #avoid to  have coordinates outside the screen
        x0 = max(0, x0)
        y0 = max(0, y0)
        x1 = min(x1, self.display.currentd.size2 - 1)
        y1 = min(y1, self.display.currentd.size1 - 1)
        self.paramz.zoom_coord = [x0, y0, x1, y1]
        
if __name__ == '__main__':
    print "hello"