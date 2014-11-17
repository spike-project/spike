from PySide import QtCore, QtGui

app = QtGui.QApplication([''])
mw = QtGui.QMainWindow()
mbar = QtGui.QMenuBar() # menu visible

#mbar = QtGui.QMenuBar(mw) # menu not visible
def hello():
    print "hello"

def add_Menu(mbar, name):
    '''
    Menu for macOsX.
    '''
    menu = QtGui.QMenu(mbar)
    menu.setTitle(name)
    mbar.addMenu(menu)
    return menu

menu = add_Menu(mbar, "Try")

a = QtGui.QAction('&Save', menu)
a.triggered.connect(hello)
menu.addAction(a)
print dir(menu)


mw.show()
app.exec_()