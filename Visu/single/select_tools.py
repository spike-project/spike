from util.debug_tools import*  
from Visu.Pyside_PyQt4 import*

@dec_class_pr
@decclassdebugging
class SELECT_TOOLS(object):
    def __init__(self, display, paramz, interf):
        self.display = display
        self.paramz = paramz
        self.interface = interf
        self.tools = ['profile', 'manual_profile', 'zoom', 'drag', 'cursor', 'scale']
        for t in self.tools:                                                            # Initializes to False
            setattr(self, t, False)
        
    def change_to(self, tool):
        '''
        Changes to the selected tool.
        '''        
                                  
        for t in self.tools:                                                            # Reinitializes to False
            setattr(self, t, False)
        print "### change tool to ", tool
        setattr(self, tool, True)     # makes the tool flag to True
        if tool == 'zoom':
            self.paramz.zoomready = False

        self.interface.ui.label_14.setText(tool)                                              
        try:
            #cursor = self.list_tool_cursor[tool]
            #self.display.qmc.setCursor(QtGui.QCursor(cursor))
            self.display.setcursor(tool)   #  if cursor tool exists pass to this kind of cursor.. 
            
        except:
            pass


    def report(self):
        '''
        Show values of the tools.
        '''
        for t in self.tools:
            attr = getattr(self, t)
            print "tool {} : {}".format(t, attr)

    