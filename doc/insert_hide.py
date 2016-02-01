import glob
import os, sys
import os.path as op
from os.path import join
from bs4 import BeautifulSoup

'''
Utility for hiding Unittest etc in the documentation with javascript.
Find htmls in th folder and insert the javascript 
synthax : python -m spike.util.insert_hide address_html_folder a'_debug' c'_Tests'
a : attribute
c : class
'''

class INSERTJS(object):
    '''
    Insert js code in htmls
    '''
    def __init__(self, folder, listinj):
        self.listinj = listinj
        self.folder = folder
        print self.listinj, self.folder

    def recurswalk(self, somepath): # Walking recursively in th folder
        for root, dirs, files in os.walk(somepath):
            for curfile in files:
                yield join(root, curfile)
                
    def insert_tag(self, inj):
        '''
        Insert a new tag
        '''        
        new_tag = self.soup.new_tag("script")        # make a script tag
        new_tag.insert(1, inj)           # insert the javascript for hiding attributes finishing in "_debug"
        self.soup.body.append(new_tag)               # append the script tag at the end of body
    
    def inject_in_all(self):
        print "inject_in_all"
        for path in self.recurswalk(self.folder):
            if path[-4:] == 'html':
                print "path is ", path
                with open(path,'r') as html_doc:
                    self.soup = BeautifulSoup(html_doc)          # Beautiful soup object
                    for inj in self.listinj:
                        self.insert_tag(inj)
                with open(path,'w') as html_doc_modif:
                    html_doc_modif.write(self.soup.prettify("utf-8")) # Rewrite the file with modification.


if __name__ == "__main__":
    addr = sys.argv[1]
    list_rem = sys.argv[2:]
    print addr, list_rem
    kindpatt = {'a' : '.attribute', 'c' : '.class'}
    listinj = [r'''
        $('{0}').each(function(){{
            if ($(this).children('dt').attr('id').search('{1}') !=-1 ) {{ $(this).hide()}}
        }})
    '''.format(kindpatt[patt[0]], patt[1:]) for patt in list_rem]
    ijs = INSERTJS(folder=addr, listinj=listinj)
    ijs.inject_in_all()
    