'''
decorate the class with decclassdebugging and in each method makes :
if debug(self) : instructions.. 
'''
import sys, os
import inspect

def pr(self, var, message = ''):
    '''
    print
    '''
    frame = inspect.currentframe()
    loc = frame.f_back.f_locals
    def show(value):
        if message == '':
            print var + ' ', value 
        else:
            print message + ' ' + var + ' ', value 
    try:
        show(eval(var))
    except Exception:
        show(loc[var])

def dec_class_pr(cl):
    setattr(cl, 'pr', pr)
    return cl

def name_up():
    frame = sys._getframe(2)
    return frame.f_code.co_name

def debug(cl):
    return getattr(cl, name_up()+'_debug')

def decclassdebugging(cl):
    '''
    decorator for debugging. Allows to trigger specifically which method in a class we want ot observe.
    for debugging you have to add in the code :   if fista.method_name_debug: 
    Then when instantiating the class add fista.method_name_debug = True to activate the debugging.
    '''
    for m in dir(cl):
        if len(m.split('__')) == 1: 
            setattr(cl, m + '_debug', False)
    return cl

@dec_class_pr
@decclassdebugging
class test(object):
    def __init__(self, ):
        obj_test = False
    
    def hip(self):
        print "passing by hip"
        if debug(self) :
            print "hip is being debugged.. hiphip hourra !"
        
    def hop(self, ):
        print "you are in hop"
        if debug(self) :
            print "hop debug"


if __name__ == '__main__':
    '''
    On this example just hip is activated..
        so debug(self) is True only in hip.
    '''
    t = test()
    t.hip_debug = True
    t.hip()
    t.hop()
    
