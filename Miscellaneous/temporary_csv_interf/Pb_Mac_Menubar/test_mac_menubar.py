from PySide import QtCore, QtGui

app = QtGui.QApplication([''])
mw = QtGui.QMainWindow()
mbar = QtGui.QMenuBar() # menu visible

#mbar = QtGui.QMenuBar(mw) # menu not visible
def hello():
    print "hello"
menu = QtGui.QMenu(mbar)
menu.setTitle("Testing")
mbar.addMenu(menu)
if True:
    
    
    #mw.setMenuBar(mbar)
    a = QtGui.QAction('&Save', menu)
    a.triggered.connect(hello)
    menu.addAction(a)
    print dir(menu)
mw.show()
#mw.raise_()
app.exec_()