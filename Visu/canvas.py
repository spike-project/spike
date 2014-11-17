from util.debug_tools import*
from Matplotlib_QtCanvas import*
from Pyside_PyQt4 import*

@dec_class_pr
@decclassdebugging
class Qt4MplCanvas(FigureCanvas): 
    """Class for integrating Matplotlib in Qt""" 
    def __init__(self, parent, paramz):
        '''
        Canvas
        '''
        #print "in class Qt4MplCanvas"
        self.fig = Figure() 
        self.axes = self.fig.add_subplot(111)
        FigureCanvas.__init__(self, self.fig) # set the parent widget 
        self.setParent(parent)
        FigureCanvas.updateGeometry(self)
        self.paramz = paramz

    def contextMenuEvent(self, event):
        '''
        Context menu
        '''
        menu = QtGui.QMenu(self)
        Action = menu.addAction("peaks")
        Action.triggered.connect(gtools.createpeaks)
        Action = menu.addAction("profile")
        Action.triggered.connect(selectL)
        menu.exec_(event.globalPos())
        
if __name__ == '__main__':
    print "hello"