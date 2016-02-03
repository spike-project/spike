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

    def recurswalk(self, somepath): # Walking recursively in th folder
        '''
        Walk in all the folders
        '''
        for root, dirs, files in os.walk(somepath):
            for curfile in files:
                yield join(root, curfile)
                
    def insert_tag(self, inj):
        '''
        Insert a new tag
        '''        
        new_tag = self.soup.new_tag("script")        # make a script tag
        new_tag.insert(1, inj)                       # insert the javascript for hiding attributes finishing in "_debug"
        self.soup.body.append(new_tag)               # append the script tag at the end of body
    
    def inject_in_all(self):
        '''
        Insert the javascripts in all the htmls
        '''
        print "### inject_in_all"
        for path in glob.glob(op.join(self.folder, '*.*')):
            if path[-4:] == 'html':
                with open(path,'r') as html_doc:
                    self.soup = BeautifulSoup(html_doc)          # Beautiful soup object
                    for inj in self.listinj:
                        if not inj in self.soup.prettify():      # if the injection is not yet done.. 
                            self.insert_tag(inj)
                        else:
                            print "### injection yet existing "
                with open(path,'w') as html_doc_modif:
                    html_doc_modif.write(self.soup.prettify("utf-8")) # Rewrite the file with modification.
        print "### injection finished "

if __name__ == "__main__":
    addr = sys.argv[1]
    list_rem = sys.argv[2:]
    kindpatt = {'a' : '.attribute', 'c' : '.class'}
    listinj = [r'''
    var tog_{2} = function(){{
            $('{0}').each(function(){{
                    if ($(this).children('dt')
                               .attr('id')
                               .search('{1}') !=-1 )
                     {{ $(this).toggle()}}
                }}) // end each
            }}; // end function
     tog_{2}();

     var def_butt = $('<button/>')
                        .text('{0}')
                        .attr('type','submit')
     var butt = $('<div/>').append(def_butt)
                           .attr('id','tog_{2}')
     butt.insertBefore('.footer')
     var stl = {{'text-align':'center'}}
     $('#tog_{2}').css(stl).click(function(){{tog_{2}()}})

    '''.format(kindpatt[patt[0]], patt[1:], patt[0]) for patt in list_rem]
    ijs = INSERTJS(folder=addr, listinj=listinj)
    ijs.inject_in_all()
    