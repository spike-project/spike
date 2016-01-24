#!/usr/bin/env python 
# encoding: utf-8

'''
Code permitting to make easily subplot without having to anticipate the number of plot. 
Just the number of columns is needed. 
typical syntax is :
    sub = subpl(nbsub_h = 2) # suplot organized in two columns
    sub.next() # adding new subplot
    sub.plot(np.arange(6), np.arange(6)**2, 'g', label = 'one') # first plot in first subplot
    sub.title('One')
'''

from __future__ import print_function

import numpy as np
import unittest

def add_meth(attr):
    '''
    Evaluates methods created in dec_class. 
    Replaces pyplot method by methods which keep the arguments. 
    for example, subpl.title() will stock its arguments in the list "subpl.ltitle"
    '''
    def wrapped(*args, **kwargs ):#
        '''
        Appends the useful arguments to the list "subpl.lmeth" instead of executing pyplot method "meth"
        '''
        eval("subpl.l" + attr + ".append([args[1:], kwargs])") # 
        if attr == 'plot': # if the method is plot, in this case add a new element to "subpl.l_num"
            eval("subpl.l_num.append(subpl.numplot)") # list for keeping the number of plot/ subplot. 
    return classmethod(wrapped)

def add_class(l): # l is the list of methods to be added for replacing matplotlib pyplot methods.
    '''
    Each method in the list l is associated to a method defined with add_meth.
    This method as the same name as the pyplot method, but it has a different function. 
    It will append the arguments passed in the plotting code to a list "lmeth" which will contains the arguments. 
    Those
    
    '''
    def dec_class(cls):
        '''
        Creates a methods from the list "l" mimicking pyplot methods. 
        and in parallel creates a list for each fake method that will stock the arguments
        passed when writing code for plotting. 
        '''
        setattr(cls, 'numplot', 0) # number of plots set to 0
        for attr in l: # read the name of methods to be mimicked. 
            setattr(cls, attr, add_meth(attr)) # eg. cls.plot equivalent to cls.lplot.append()
            setattr(cls, 'l' + attr, []) # makes the lists, eg. lplot=[]
            if attr == 'plot':
                setattr(cls, 'l_num', []) # class to tregister num of plots
        return cls
    return dec_class

l = ['title', 'xlabel', 'ylabel', 'plot', 'xticklabels', 'ax'] #

@add_class(l) # adding method keeping the arguments before execution.
class subpl():
    '''
    Iterator for returning automatically
    the number of the subplot iteratively.
    nbsub_h is the number of subplot horizontally
    which_plot is the kind of plot, fake plot or matplotlib pyplot plot.. 
    '''
    def __init__(self, nbsub_h = 1, which_plot = None):
        global plt
        if which_plot:
            print("using injected pyplot")
            plt = which_plot # possibility to replace plot by fakeplot with testplot.
        else:
            plt = globals()['plt']
        self.fig = plt.figure()     # Make a new figure
        subpl.numplot = 0           # reinitialize the number of subplots
        self.nbsub_h = nbsub_h          # number of horizontal subplots
        self.nbsub_v = 0            # number of vertical subplots
        for meth in dir(self):
            if meth[0] == "l":
                setattr(subpl, meth, []) #
    
    def next(self):
        '''
        Increments the number of plots "subpl.numplot"
        And calculates the number of vertical lines. 
        '''
        subpl.numplot += 1          # number of plot created.
        self.nbsub_v = (subpl.numplot-1)/self.nbsub_h +1        # vertical number of lines used for the whole plot. 
    
    def read_kargs(self, l):
        '''
        '''
        lkargs = ''
        for name in l:
            n = l[name]
            if type(n) == str: nkargs = '"'+ n + '"'
            else: nkargs = str(n)   
            lkargs += ',' + name + '=' + nkargs
        return lkargs
    
    def make_str_kargs(self, func, i):
        '''
        Extracts the list containing the arguments from sublist [1] of each method,
        in the list "self.lmeth".
        '''
        try:
            str_kargs = self.read_kargs(eval('self.l'+ func + '[i][1]')) 
        except :
            str_kargs = ''
        return str_kargs
    
    def select_pos(self, func):
        '''
        
        '''
        if func == 'plot':
            pos = '[i][0]'
        else:
            pos = '[numpl-1][0]'
        return pos
    
    def show(self):
        '''
        Takes all the list of arguments from the fake methods sub.meth
        And evaluates them one after the other with the correct pyplot corresponding method. 
        '''
        nbticks = None
        for i, numpl in enumerate(self.l_num): # read all subplots
            self.ax = self.fig.add_subplot(self.nbsub_v, self.nbsub_h, numpl )
            for func in l: # l = ['title', 'plot', 'xlabel', 'ylabel' etc.. ]
                str_kargs = self.make_str_kargs(func, i) # makes string for key arguments
                pos = self.select_pos(func)
                arg =  '(*self.l'+ func + pos + str_kargs +')' # unfolds the list of argumetns from self.lmeth
                if func != 'plot':
                    try:
                        expr_eval = 'self.ax.set_'+ func + arg
                        eval(expr_eval)  # set parameters before plotting
                    except:
                        pass
                    if func == 'xticklabels':
                        try:
                            nbticks = len(self.lxticklabels[numpl-1][0][0])-1
                        except :
                            nbticks = None
                #### Making the plots
                elif func == 'plot': # equivalent of plt.plot
                    expr_eval = 'self.ax.' + func + arg
                    eval(expr_eval) # plot
                    if str_kargs != '':
                        plt.legend()
        if nbticks: # number of ticks
            if not hasattr(plt,'FAKE'):     # fakeplot defines FAKE but not MaxNLocator
                plt.gca().xaxis.set_major_locator( MaxNLocator(nbins = nbticks) )
        plt.show()

class Test_dynsubplot(unittest.TestCase):
    '''
    Unittests for dynsubplot with four subplots.
    First plot contains 2 plots. 
    '''      
    def test_dynsub(self):
        # from matplotlib import pyplot as plt
        import spike.Display.testplot as testplot
        global plt
        plt = testplot.plot()
        ## from matplotlib.ticker import MaxNLocator
        sub = subpl(nbsub_h = 2) # suplot organized in two columns
        sub.next() # adding subplot 1
        sub.plot(np.arange(6), np.arange(6)**2, 'g', label = 'one') # first plot in first subplot
        sub.title('One')
        sub.ylabel('y1')
        sub.xlabel('x1')
        sub.plot(np.cos(np.arange(11)), label = 'two') # plotting second plot in first subplot
        sub.next()# adding subplot 2
        sub.title('Second sub')
        sub.plot(np.sin(np.arange(15)), label = 'three') # plotting second plot in second subplot
        sub.ylabel('y3')
        sub.next()# adding subplot 3
        sub.ylabel('y4')
        sub.title('Third')
        sub.plot(np.arange(21), 'g', label = 'four') # plotting first plot in third subplot
        sub.next()# adding subplot 4
        sub.ylabel('y5')
        sub.title('Forth subplot')
        sub.plot(np.cos(np.arange(21))-np.arange(21), 'r', label = 'five') # plotting first plot in fourth subplot
        sub.show()
    
        sub = subpl(nbsub_h = 2) # suplot organized in two columns
        sub.next() # adding subplot 1
        sub.plot(np.arange(10), np.arange(10)**3, 'g', label = 'one') # first plot in first subplot
        sub.title('One')
        sub.ylabel('y1')
        sub.xlabel('x1')
        sub.plot(np.sin(np.arange(20)), label = 'two') # plotting second plot in first subplot
        sub.next()# adding subplot 2
        sub.title('Second sub')
        sub.plot(np.cos(np.arange(41)), label = 'three') # plotting second plot in second subplot
        sub.ylabel('y3')
        sub.show()

if __name__ == '__main__':
    unittest.main()
    

    
