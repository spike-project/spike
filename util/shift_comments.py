#!/usr/bin/env python 
# encoding: utf-8

'''
Usage:  python shift_comment.py file.py nbcolumn
nb column is the column nb at which we want the comment to be placed.
If no nb is given will shift the comments from the end of the line at a regular intervall.
'''

from __future__ import print_function
import sys
import re

class SHIFT_COMMENT():
    '''
    '''
    def __init__(self, namefile, poscomment = None):
        self.shift_comm = 10 # shift the comment to ten column on the right from the end of the line.
        self.namefile = namefile
        self.f = open(namefile, 'r')
        self.g = open(namefile[:-3] + "corr.py", 'w')
        self.poscomment = int(poscomment) # desired position for the comments.
    
    def make_linecorr_fixedpos(self, posin, comment, line):
        '''
        Correction of the line for fixed position
        '''
        if self.poscomment > posin:  # simpler case, the desired position is longer than the line.
            self.linecorr = line[: posin] + (self.poscomment - posin)*' ' + comment # corrected line
        else: # the desired position is shorter than the line.
            self.linecorr = line[: self.poscomment] + comment # corrected line
        print(self.linecorr)
    
    def make_linecorr_shiftedpos(self, lenline, comment, line):
        '''
        Correction of the line for shifted comment position
        '''
        self.linecorr = line[:lenline] + self.shift_comm*' ' + comment # corrected line
 
    def read_correc(self):
        '''
        Read the file and moves the comments to the right place.
        '''
        for line in self.f.readlines():
            if re.findall('# ', line) != []:        # find the existence of comment.. 
                print(line)
                posin = line.find('# ')
                print(posin)
                lenline = len(line[:posin].rstrip())
                comment = line[posin: ]         # take all the line from posin to the end
                if self.poscomment:
                    self.make_linecorr_fixedpos(posin, comment, line)
                else:
                    self.make_linecorr_shiftedpos(lenline, comment, line)
                self.g.write(self.linecorr)
            else:
                self.g.write(line)
    
if __name__== '__main__':
    namefile = sys.argv[1]
    poscomment = sys.argv[2]
    SC = SHIFT_COMMENT(namefile, poscomment)
    SC.read_correc()