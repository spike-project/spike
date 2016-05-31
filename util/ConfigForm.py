#!/usr/bin/env python 
# encoding: utf-8
from __future__ import print_function
"""
ConfigForm.py

Handles ConfigParser files via a html interface

Creates a HTML form from the cfg file and store a modified file

usage : python ConfigForm.py configfile.cfg

Entry in the config files should of the form :

# comment on section1
# comment on section1
[section1]

# comment on key1
# comment on key1
# %tag%
key1 = value

# comment on key2
# comment on key2
# %tag% %validators% %formaters%
key2 = value

Remark 
valid %tag% are    and %validators% (optionnal except for %options% )

'%text%':           %min% %max% %extension%     default tag is tag is absent
'%float%':          %min% %max% 
'%integer%':        %min% %max% 
'%radio%':          %options%
'%select%':         %options%
'%file%':           %extension%
'%infile%':         %extension%
'%outfile%':        %extension%
'%boolean%':            -
'%hidden%':

valid %validators% are
%min:val%    %max:val%   : limit on the value or on the length of the string
%options:list:of:valid:values%  : list of options
%extension:.txt% : constraint the entry to finish by a given string

valid %formaters% are
%length:n% : contraint the length of the HTML text field to be n char long
%tooltip:text% : add text to the field tooltip - some html formating can be used

Created by DELSUC Marc-Andre on 2015-04-15
"""
"""
Requires :

Flask
wtforms
jinja2


Still to be done
BUGS
- doc of first section
- booleans
- debug %extension    and make it multiple
- two options with same name in 2 diff sections collide !
improvments
- better layout / css
- find a way to stop the app
- use flash message for validation
"""

import sys
import os
import unittest
from collections import OrderedDict, defaultdict
import ConfigParser
import threading, webbrowser
import random
import socket

import wtforms
from wtforms import Form, validators, ValidationError
from jinja2 import Template
from flask import Flask, Response, render_template, request, url_for, redirect
from PyQt4 import QtGui

###########
HOST = "localhost"      # the name of the web service, if localhost, the service is only local
PORT = 5000             # the port of the service - choose a free one
DEBUG = False
STARTWEB = False         # If true, opens a browser when launched
METHOD = 'POST'          # 'POST' or 'GET'
__version__ = "0.2"     # 0.1 first trial 0.2 performs ok but still rough
########## Graphical variables
COLOR_TOOLTIP = "PapayaWhip"
###########
#mini templates
"""

"""
head = """
<!DOCTYPE html>
<html>
<head>

<script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.2/jquery.min.js"></script>
<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.4/js/bootstrap.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/bootstrap-select/1.7.2/js/bootstrap-select.min.js"></script>
<script>

$(document).ready(function(){$('[data-toggle="tooltip"]').tooltip();});

</script>

<!-- Css -->
<link rel="stylesheet" href="http://maxcdn.bootstrapcdn.com/bootstrap/3.3.4/css/bootstrap.min.css"><!-- Scripts -->
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/bootstrap-select/1.7.2/css/bootstrap-select.css">

<style>
.CF-tooltip + .tooltip > .tooltip-inner { background-color: PapayaWhip; text-align : left}
.CF-tooltip + .tooltip > .tooltip-arrow { border-bottom-color: PapayaWhip; }
#namefile, #namefile_button
{
    display:inline;
}


</style>

<title>{{title}}</title>
</head>

""" 

foot = """
<p><small><i>ConfigForm, version %s  - author : M.A. Delsuc</i></small></p>
"""%(__version__)

# known Field types syntax {%key%: (Field,type)}

cfFieldDef = {
    '%boolean%': (wtforms.StringField, bool),
    '%text%': (wtforms.StringField, str),
    '%float%':   (wtforms.FloatField, float),
    '%integer%':   (wtforms.IntegerField, int),
    '%radio%':   (wtforms.RadioField, str),
    '%file%':   (wtforms.FileField, str),
    '%infile%':   (wtforms.FileField, str),
    '%outfile%':   (wtforms.StringField, str),
    '%select%':   (wtforms.SelectField, str),
    '%hidden%':  (wtforms.HiddenField, str)
}
cfFieldKeys = cfFieldDef.keys()
if DEBUG: print(cfFieldKeys)

# known validators syntax {%key: Validator}  # Validator field not used so-far
cfvalDef = {
    "%options": validators.AnyOf,
    "%extension":  validators.Regexp,
    "%min": validators.NumberRange,
    "%max": validators.NumberRange,
    "%length": None,   # length is not a WTForm validator, but a HTML format
    "%tooltip": None   # tooltip is not a WTForm validator, but a HTML format
}
cfvalDefKeys = cfvalDef.keys()
if DEBUG: print(cfvalDefKeys)

############

def dForm(dico):
    """
    This function creates a wtforms object from the description given in the (Ordered) dict dico 
    dico = {"name_of_field" : Field(), ...}

    This is the key, wtforms requires a class describing the form, we have to build it dynamically.
    """
    class C(Form):
        pass
    for i,j in dico.items():
        setattr(C,i,j)
    return C
def dynatemplate(kw, method, action, class_):
    "builds a simplistic template string, see  dTemplate()  - used for tests -"
    temp = ['<form method="{0}" action="{1}">'.format(method,action)]
    for i in kw:
        temp.append('<div>{{{{ form.{0}.label }}}}: {{{{ form.{0}(class="{1}") }}}}</div>'.format(i,class_) )
    temp.append('</form>')
    return "\n".join(temp)
def dTemplate(dico, method, action, class_):
    """
    This function creates a simple Template object from the description given in the (Ordered) dict dico
    method is either POST or GET
    action is the callback url
    class_ is the css class attached to each field
    
    used for tests
    """
    return Template(dynatemplate(dico, method, action, class_))
    
class FIND_FOLDER_FILE(QtGui.QWidget):
    '''
    PyQt interface for chosing the folder for the processing.
    '''
    def __init__(self, parent = None):
        QtGui.QWidget.__init__(self, parent)
        self.curDir = os.getcwd()

    def browse(self, kind):
        '''
        Search for a folder or a file
        '''
        print("saving config")
        if kind == "folder":
            self.selected = QtGui.QFileDialog.getExistingDirectory(self, "Select Folder",  self.curDir)
        elif kind == "file":
            self.selected = QtGui.QFileDialog.getOpenFileName(self,'Open file', self.curDir)
        print("#################### self.selected ", self.selected)

def search_folder_file(kind):
    '''
    Opens the PyQt interface for searching the folder
    When the folder is chosen the program stops the Qt interface. 
    '''
    app = QtGui.QApplication(sys.argv)
    ff = FIND_FOLDER_FILE()
    ff.browse(kind) 
    return ff.selected
   

############################
class cfField(object):
    """
    This class holds all the parameters of an option in the config file

    subparse_comment() loads it    
    """
    def __init__(self, name, value):
        self.name = name
        self.value = value
        self.class_ = "CF"
        self.type_ =  "%text%"
        self.length = 60
        self.tooltip = []
        self.tooltip_color = COLOR_TOOLTIP
        self.doc = ""
        self.section = ""
        self.lopt = None
        self.meta = None      # meta will be used by produce()
        self.validators = [validators.DataRequired()]   # all field are required, initiailize the list
        
    def subparse_comment(self, line):
        """parse the comment, and search for meta comments"""
        import re
#        words = line.split()
        words = re.findall("%.+?%", line)
        if not words:
            words = [line]
        print(words)
        word = words.pop(0)     # get first one
        if word not in cfFieldKeys:  # we do not have a meta comment
            self.doc.append(line)
            return
        # now we do !
        self.type_ = word
        self.meta = line
        for word in words:      # then process meta
            val = word[:-1].split(":")  # extract sub options, removing the final %
            if val[0] in cfvalDefKeys:  # we have a validator
                if val[0] == "%min":
                    if cfFieldDef[self.type_][1] == str:     # min and max are both numbers and string
                        self.validators.append( validators.Length(min=int(val[1])) )
                        self.tooltip.append("minimum length : %d"%int(val[1]))
                    else:
                        self.validators.append( validators.NumberRange(min=float(val[1])) )
                        self.tooltip.append("minimum value : %d"%float(val[1]))
                elif val[0] == "%max":
                    if cfFieldDef[self.type_][1] == str:     # min and max are both numbers and string
                        self.validators.append( validators.Length(max=int(val[1])) )
                        self.length = int(val[1])
                        self.tooltip.append("maximum length : %d"%int(val[1]))
                    else:
                        self.validators.append( validators.NumberRange(max=float(val[1])) )
                        self.tooltip.append("maximum value : %d"%float(val[1]))
                elif val[0] == "%extension":
                    print("%Etension:xxx% is currently not implemented")
                    print("XX EXT", val[0], '%s$'%(val[1],))
                    self.tooltip.append("file extension : %s"%(val[1],) )
#                    self.validators.append( validators.Regexp(r'%s$'%val[1], message="should end with %s"%val[1]) )
                elif  val[0] == "%options":
                    lopt = val[1:]
                    self.lopt = lopt
                    self.validators.append( validators.AnyOf(values=lopt ) )
                    self.tooltip.append("choose one entry")
                elif val[0] == "%tooltip":
                    self.tooltip.append(val[1])
                elif val[0] == "%length":
                    try:
                        lgt = int(val[1])
                    except:
                        lgt = 60
                        print("WARNING, wrong value in line ",line)
                    self.length = lgt
    def opt_templ(self):
        "build the template for an option stored in cfField"
        if self.type_ != "%hidden%":
            doc = "<br/>\n".join(self.doc)
            
            style_tooltip = "'background-color:{0} ; color:Black' ".format(self.tooltip_color)
            html_tooltip = '<p style=' + style_tooltip + '>' +  '<br/>'.join(self.tooltip) + '</p>'
            
            if self.tooltip:
                ttip = 'data-toggle="tooltip" data-html="true" data-placement="top"  title="{}" class = "CF-tooltip" '.format(html_tooltip)
                print("ttip ", ttip)
            else:
                ttip = ""
            # then produce
            if self.type_ in ("%text%", "%outfile%"):    # if string, add SIZE
                print("self.section, self.name ", self.section, self.name)
                templ = '<div class="{2}"><p>{3}<br/><a style="color:black" href="#" {4}>{{{{ form.{0}_{1}.label}}}}</a>: \
                    {{{{form.{0}_{1}(class="{2}",SIZE="{5}")}}}}</p></div>'.format( \
                                self.section, self.name, self.class_, doc, ttip, self.length)
            
            elif self.type_ in ("%select%"): # Select boxes with Bootstrap class selectpicker
                print("self.section, self.name ", self.section, self.name)
                templ = '<div class="{2}" ><p>{3}<br/><a style="color:black" href="#" {4}>{{{{ form.{0}_{1}.label}}}}</a>: \
                         {{{{form.{0}_{1}(class="selectpicker")}}}} </p></div>'.format( \
                                self.section, self.name, self.class_, doc, ttip)     
                                
            elif self.type_ in ("%infile%"):
                templ = ''' 
                    </br>
                    <button  type="submit" name="submitform" value="chose_file-{0}_{1}"  class="btn btn-default"> Chose File </button>
                      {{% if form.{0}_{1}.infilename %}}
                        <p style="color:blue;"> {{{{form.{0}_{1}.infilename}}}} </p>
                      {{% endif %}}
                '''.format(self.section, self.name)

            elif self.type_ in ("%boolean%"):    # Boolean checkboxes
                 print("self.section, self.name ", self.section, self.name)
                 
                 templ = '<div class="{2}"><p>{3}<br/><a style="color:black" href="#" {4}>{{{{ form.{0}_{1}.label}}}}</a>:\
                   <input type="checkbox" name="{0}_{1}" checked> </p>  </div>  \n'.format(self.section, self.name, self.class_, doc, ttip)
            else:
                print("self.section, self.name ", self.section, self.name)
                templ = '<div class="{2}"><p>{3}<br/><a style="color:black" href="#" {4}>{{{{ form.{0}_{1}.label}}}}</a>: \
                    {{{{form.{0}_{1}(class="{2}")}}}}</p></div>'.format( \
                                self.section, self.name, self.class_, doc, ttip)                
        else:   # special if hidden 
            templ = '{{{{form.{0}_{1}}}}}'.format(self.section, self.name)
        return templ
    def __repr__(self):
        return "%s_%s %s %s %s\n%s"%(self.section, self.name, self.value, self.type_, self.validators, self.doc)
    def __str__(self):
        return self.__repr__()
class ConfigForm(object):
    """
    This class reads a configparser file
    creates a HTML form for it
    starts a broswer on it (using Flask) which allows to modify the values and store them back.

    """
    def __init__(self, cffilename):
        self.cffilename = cffilename
        self.cp = ConfigParser.ConfigParser()
        self.configfile = open(cffilename).readlines()
        self.cp.read(cffilename)
        self.doc_sec = {}   # contains comments for sections
        self.doc_opt = {}   # and options
        self.options = OrderedDict()   # contains all cfFormField()
        self.read_comments()
        self.parse_comments()
        self.method = METHOD
        self.callback = '/callback'
        self.class_ = "CF"      # this will be used in html templates
        self.port = PORT 
        self.host = HOST
        self.template = None
        self.form = None
    def reload(self):
        "reload the content from the file"
        self.configfile = open(self.cffilename).readlines()
        self.cp = ConfigParser.ConfigParser()
        self.cp.read(self.cffilename)
        self.read_comments()
        self.parse_comments()
        self.template = None        # clear cache
        self.form = None            # clear cache
    def read_comments(self):
        """
        this method reads the comments associated to sections and options
        comments should be put BEFORE the commented entry
        ## comments (double #) are silently skipped
        """
        currentsection = "initial"
        currentoption = None
        currentcom = []
        doc_sec = defaultdict(list)
        doc_opt = defaultdict(list)
        for l in self.configfile:
            l.strip()
            if l.startswith('##'):    # skip double comments
                continue
            if l.startswith('#'):     # store comments
                currentcom.append(l[1:].strip())
                continue
            m = self.cp.SECTCRE.match(l)
            if m:  # a new section
                currentsection = m.group('header').strip().lower()  # ConfigParser puts all keys in lower case
                currentoption = None
                doc_sec[currentsection] += currentcom
                currentcom = []
            m = self.cp.OPTCRE.match(l)
            if m:  # a new option
                currentoption = m.group('option').strip().lower()
                doc_opt[currentoption] += currentcom
                currentcom = []
        self.doc_sec = doc_sec
        self.doc_opt = doc_opt
    def parse_comments(self):
        """parse self.doc_opt and generates self.options"""
        for sec in self.cp.sections():
            for opt in self.cp.options(sec):
                cle = "%s_%s"%(sec,opt)
                cff = cfField(opt, self.cp.get(sec, opt))
                cff.doc = []
                cff.section = sec
                cff.opt = opt
                for l in self.doc_opt[opt]:
                    cff.subparse_comment(l)
                if DEBUG:
                    print(cff)
                    print()
                self.options[cle] = cff

    def sect_templ(self, sec):
        "build the template for a section 'sec' "
        doc = "<br/>\n".join(self.doc_sec[sec])
        print("doc ", doc)
        options_list = [cff.opt_templ()   for cff in self.options.values()  if cff.section == sec] 
        
        options_templ = "  \n".join(options_list )
        print("options_templ ", options_templ)
        sec_templ = """
<div class={0}.section>
<hr>
    <div class="jumbotron" style="background-color:Bisque;">
    <h2  class="text-center">{1}</h2>
    </div> 
<div class="container">
<b>{2}</b>
{3}
</div>
</div>

        """.format(self.class_, sec, doc, options_templ)
        return sec_templ

    def buildtemplate(self):
        "builds the jinja2 template string, from the config file"
        templ = [head]
        templ += ['<body>']
        templ += ['<div class="container">']
        templ += ['''
            <h1 style="text-align: center;">Configuration for file :  </h1>
            <h1 style="text-align: center;">{{filename}}</h1>
        ''']
        templ += ['<form method="{0}" action="{1}">'.format(self.method, self.callback)]
        for sec in self.cp.sections():
            # print "#################   self.sect_templ(sec) ", self.sect_templ(sec)
            templ.append( self.sect_templ(sec) )
        templ.append('<hr>\n    <input type="submit" name="submitform" value="Reload Values from file" class="btn btn-default" />')
        templ.append('    <input type="submit" name="submitform" value="Validate Values" class="btn btn-default"/>\n')
        templ.append('</form>')
        templ.append(foot)
        templ += ["</div>"] # Ending container
        templ += ['''
</body>
</html>
        ''']
        return "\n".join(templ)
    def buildforms(self):
        "build the wtform on the fly from the config file"
        dico = OrderedDict()
        values = {}
        for cle,content in self.options.items():        # construct the dict for dForm
            # print "cle, content", cle, content
            Field_type = cfFieldDef[content.type_][0]
            if content.type_ == "%select%":
                dico[cle] = Field_type(content.name, choices=zip(content.lopt,content.lopt), validators=content.validators)
            else:
                dico[cle] = Field_type(content.name, validators=content.validators)
            values[cle] = content.value
        df = dForm(dico)                    # create the form builder
        form = df(**values)                 # create the form and load the values at once
        return form

    def render(self):
        'associate template and wtform to produce html'
        if self.form is None:                   # simple caching method
            self.form = self.buildforms()
        if self.template is None:               # simple caching method
            self.template = Template( self.buildtemplate() )
        
        html = self.template.render(form=self.form, title="ConfigForm", filename=self.cffilename)
        return html
    
    def produce(self):
        """produces the text of the current defined config file, as a list, one line per entry"""
        text = ["## File generated by ConfigForm v %s"%(__version__)]
        for sec in self.cp.sections():
            text.append('##########################################')
            text += ["# %s"%d for d in self.doc_sec[sec] ]
            text.append('[%s]'%sec)
            text.append('')
            for opt in self.cp.options(sec):
                cle = "%s_%s"%(sec,opt)
                text += ["# %s"%d for d in self.options[cle].doc ]
                if self.options[cle].meta is not None:
                    text.append('# %s'%(self.options[cle].meta) )
                text.append('%s = %s'%(opt,self.form[cle].data))
                text.append('')
        # print "######## text produced is ", text
        return text
    def writeback(self, filename):
        """writes back the config file with the current values"""
        with open(filename,'w') as F:
            F.write( "\n".join( self.produce() ) )
        
##########################
def BuildApp(cffile):
    app = Flask(__name__, static_url_path='/static')
    app.config['SECRET_KEY'] = 'kqjsdhf'
    cf = ConfigForm(cffile)
    cf.callback = '/callback'
    # cf.ask_folder = '/ask_folder'
    if DEBUG: print(cf.doc_sec.keys())

    @app.route('/')
    def index():
        '''
        First page when launching the interface.
        '''
        return cf.render()

    @app.route('/callback',  methods = [cf.method])
    def callback():
        '''
        route called by the form - validate entries
        '''
        principal_types = [str, float, int] # All non boolean types used. 
        bool_cles = [] # List for registering boolean keys. 
        not_bool_errors = 0 # counters for non boolean errors after validation.
        ###
        text = [head]
        text += ['<body>']
        text += ['<div class="container">']

        if request.form["submitform"] == "Reload Values from file":
            cf.reload()
            return redirect(url_for('index'))
            
        elif request.form["submitform"].find("chose_file") != -1: # Searching for a folder
            print(" section and name for chose file is : ", request.form["submitform"].split("-")[1])
            cf.infilename = request.form["submitform"].split("-")[1]
            return redirect(url_for('ask_folder'))
            
        else :   # Validate request.form["submitform"] == "Validate Values"
            for cle, content in cf.options.items():
                print("############## cle", cle)
                print("############## content", content)
                type_ = cfFieldDef[content.type_][1]
                if  type_ in principal_types : # 
                    try:
                        cf.form[cle].data = type_((request.form[cle]))
                    except:
                        try:
                            cf.form[cle].data = request.form[cle]
                        except:
                            try:
                                cf.form[cle].data = cf.form[cle].infilename # Case infile, handled with PyQt. 
                            except:
                                print("no file chosen")
                                cf.form[cle].data = None # None needed for test with cf.form.validate()
                else : # checkbox
                    if request.form.getlist(cle) == []:
                        cf.form[cle].data = False
                    else:
                        cf.form[cle].data = True
                    bool_cles.append(cle) # Append the boolean keys for avoiding validation error on booleans. 
                if DEBUG:
                    print("DEBUG")
                    text.append("%s : %s<br/>"%(cle, request.form[cle]))
            text_fine = ["<h1>Everything is fine</h1>\n<p>All entries validated</p>"]
            text_fine.append('<p> <a href="{}">write file {}</a>'.format(url_for('write'), cf.cffilename) )
            if not cf.form.validate() and len(bool_cles) > 0:  # error in validation
                text_error = []
                for name, msges in cf.form.errors.items():
                    nmspl = name.split('_')
                    sec = nmspl[0]
                    opt = "_".join(nmspl[1:])
                    name_cle = sec + '_' + opt
                    if not name_cle in bool_cles:
                        not_bool_errors += 1 # increment the not boolean errors
                        for msg in msges:
                            text_error.append('<li>in section <b>{}</b> entry <b>{}</b> : {}</li>\n'.format(sec, opt, msg))
                if not_bool_errors > 0:       
                    text.append("<h1>Error in config file</h1><p>The following error are detected :</p>")
                    text += text_error        
                    text.append("</ul>")
                else:
                    text += text_fine # Case only boolean errors.
            else:
                text += text_fine # No error at all 
            text.append('<p class="btn btn-default"> <a href="%s">back to form</a> </p> '%(url_for('index'),) )
            text.append(foot)
            text += ["</div>"] # Ending container
            text += ['''
            </body>
            </html>
                    ''']
            valid = "\n".join(text)
            return valid
    @app.route('/write')
    def write():
        tmpf = "_tempfile.cfg"
        while os.path.exists(tmpf):
            tmpf = "_tempfile_%d.cfg"%(int(1E6*random.random()))
        print(tmpf)
        cf.writeback(tmpf)      # write to tempfile
        os.rename(cf.cffilename, cf.cffilename+"~") # move away initial file
        os.rename(tmpf, cf.cffilename)
        return redirect(url_for('bye'))
        
    @app.route('/ask_folder')
    def ask_folder():
        cf.form[cf.infilename].infilename = search_folder_file('file') # opens PyQt interface for searching the folder.
        return redirect(url_for('index'))
        
    @app.route('/show')
    def show():
        return "<pre>" + "\n".join( cf.produce() ) + "</pre>"
    @app.route('/bye')
    def bye():
        'quitting'
        return "<H1>Bye !</H1>"
        sys.exit(0)
    return app

##########################
class ConfigFormTests(unittest.TestCase):
    "unittests"
    def setUp(self):
        self.verbose = 1    # verbose >0 switches on messages
    def announce(self):
        if self.verbose >0:
            print(self.shortDescription())
            print("----------")
    def _test_dForm(self):
        """testing dForm"""
        self.announce()
        otest = OrderedDict()
        otest['urQRd'] = wtforms.SelectField('do urQRd processing', choices=[(True,'yes'),(False,'no')] )
        otest['urank'] = wtforms.IntegerField('choose urQRd rank',[validators.NumberRange(min=4, max=100)])
        otest['File'] = wtforms.FileField('File to process')
        df = dForm(otest)
        form2 = df(urank=123)
        temp = dTemplate(otest, "POST", "/param", "classX")
        html = temp.render(form=form2)
        print(html)
        self.assertTrue(html.startswith('<form method="POST" action="/param">') )
    def test_render(self):
        "test render"
        filename = "test.cfg"
        cf = ConfigForm(filename)

##########################
def main():
    "called at start-up"
    global HOST, PORT, DEBUG, STARTWEB
    
    import argparse
        # print >> sys.stderr, "wrong input"
        # print >> sys.stderr, "Usage: python ConfigForm.py configfile.cfg"
        # return 2

    # Parse and interpret options.
    parser = argparse.ArgumentParser()
    parser.add_argument('configfile', nargs='?', default='config.cfg', help='the configuration file to analyse, default is config.cfg')
    parser.add_argument('--doc', action='store_true', help="print a description of the program")
    parser.add_argument('-w', '--webserver', default=HOST, help="the hostname of the server, default is %s"%HOST)
    parser.add_argument('-p', '--port', type=int, default=PORT, help="the port on which the servers run, default is %d"%PORT)
    parser.add_argument('-s', '--start', help="start the browser on http://WEBSERVER:PORT", action='store_true')
    parser.add_argument('-d', '--debug', default=DEBUG, help="enter debug mode", action='store_true')
    args = parser.parse_args()

    if args.doc:
        print(__doc__)
        sys.exit(0)

    PORT = args.port
    DEBUG = args.debug
    STARTWEB = args.start
    HOST = args.webserver
    filename = args.configfile

    print("Processing ",filename)
    input_file = open(filename).readlines() # read the file
                                            # a bad quick hack for checking the file is ok.
    url = "http://{0}:{1}".format(HOST, PORT)
    print(url)
    app = BuildApp(filename)
    if STARTWEB:
        threading.Timer(1.5, lambda: webbrowser.open(url)).start() # open a page in the browser.
    try:
        app.run(host=HOST, port = PORT, debug = DEBUG)
    except socket.error:
        print("wrong port number", file=sys.stderr)
        print("try using another port number (defined in the application header)", file=sys.stderr)
        return 3
    return 0

if __name__ == '__main__':
#    unittest.main()
    sys.exit(main())
    
