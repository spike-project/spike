#!/usr/bin/env python 
# encoding: utf-8

"""
This class extend the standard dictionnary behaviour.
It is used to read and store all the parameters used by the NPK program.

- it adds the possibility to load from files and dump/store to files the content of the dictionnary
  the file will be of the form :
  key=value, with one entry per line

- when extracted from the dictionnnary, values are interpreted in several manner
    o 1K is replaced by 1024 and all multiples (4k; 24k; etc...)
    o 1M is replaced by 1024*1024 and all multiples (4M; 24M; etc...)
    o all arithmetic is evaluated : math.pi/2, 0.1*get_si1_2d(), etc...
  
"""

from __future__ import print_function

__author__ = "Marc A. Delsuc <delsuc@igbmc.fr>"
__date__ = "Oct 2009"

# UserDict is used for 2.1 compatibility required by jython

import collections
try:
    from UserDict import UserDict
except ImportError:
    from collections import UserDict


#---------------------------------------------------------------------------
def parse(line):
    import re
    """Parse one property entry and prepare it to be entered in the  dictionnary
    """
    l = re.split(r'(?<!\\)=', line,1)   #matches = but not \=
#    dkey = re.sub(r'\\=', '=',l[0]).lower()
    dkey = re.sub(r'\\=', '=',l[0])
    val = re.sub(r'\\=', '=',l[1])      # replaces \= by =
    
#    l = line.strip().split("=",2)
#    dkey=l[0].lower()
#    val = l[1]

#    fval=NPKevaluate(val)
#    print dkey,fval
    return (dkey,val)

#---------------------------------------------------------------------------
def NPKevaluate(val):
    """
    - val is interpreted in several manner
        o 1K is replaced by 1024 and all multiples (4k; 24k; etc...)
        o 1M is replaced by 1024*1024 and all multiples (4M; 24M; etc...)
        o all arithmetic is evaluated : math.pi/2, 0.1*get_si1_2d(), etc...
    """
    import math
    # from NPK import *   # will be used by eval()
    found=0
    try:    # try eval()
        fval=eval(val)
        if  isinstance(fval, collections.Callable) :  # we want to avoid functions !
            pass
#            print " callable"
        else:
#            print " worked"
            found=1
    except:
#        print " failed"
        pass
    if (not found):
        try:    # try to type value
            if ( val.endswith("k")):
                fval = 1024*int(val[0:len(val)-1])
                found=1
            elif ( val.endswith("M")):
                fval = 1024*1024*int(val[0:len(val)-1])
                found=1
            else:
                fval = float(val)
                found=1
        except:
            pass
    if (not found):
        try:    # try to type value
            if ( val.lower() == "false" ):
                fval = 0
                found=1
            elif ( val.lower() == "true" ):
                fval = 1
                found=1
        except:
            fval = val
            found=1
    if (not found):
        fval=val
 #   print val, fval
    return fval


#---------------------------------------------------------------------------
class NPKParam(UserDict):
    "NPKParam class handles the parameter set used by NPK"
    #---------------------------------------------------------------------------
    def __init__(self,dict_=()):
        import copy
        if dict_ == ():
            ldict = {}
        else:
            ldict = copy.copy(dict_)
        self.data=ldict
    #---------------------------------------------------------------------------
    def __str__(self):
        return 'parameters :\n'+repr(self)

    #---------------------------------------------------------------------------
    def __getitem__(self,key):
        """overloading of the dictionnary [] method
        evaluates the content
        """
        i=NPKevaluate(self.data[key])
        return i

    #---------------------------------------------------------------------------
    def raw(self,key):
        """
        acces to the raw data present in the dictionnary, necessary in certain cases
        """
        return self.data[key]

    #---------------------------------------------------------------------------
    def dump(self,fname):
        """
        dump the content of the Parameter dictionnary as a property list file
        entries are not evaluated
        
        one entry per line with the following syntax :
        entry=value
        
        """
        try:
            fout = open(fname, 'w')
        except:
            print("Error while opening file :", sys.exc_info()[0])
            raise

        fout.write("#Property list file, dumped :"+time.strftime("%a, %d %b %Y %H:%M:%S %Z", time.localtime())+"\n")

        for i in self.keys():
            lk = re.sub('=', r'\\=', i)      # replaces = by \=
            lv = re.sub('=', r'\\=', str(self.data[i]))      # replaces = by \=

            fout.write(lk+"="+lv+"\n")
        fout.close()

    #---------------------------------------------------------------------------
    def store(self,fname):
        """
        write the content of the Parameter dictionnary as a property list file
        equivalent to dump, BUT entries are evaluated
        
        one entry per line with the following syntax :
        entry=value
        
        """
        try:
            fout = open(fname, 'w')
        except:
            print("Error while opening file :", sys.exc_info()[0])
            raise

        fout.write("#Property list file, output :"+time.strftime("%a, %d %b %Y %H:%M:%S %Z", time.localtime())+"\n")

        for i in self.keys():
            lk = re.sub('=', r'\\=', i)      # replaces = by \=
            lv = re.sub('=', r'\\=', self[i])      # replaces = by \=

            fout.write(lk+"="+lv+"\n")
        fout.close()

    #---------------------------------------------------------------------------
    def load(self,fname):
        """
        load a property list file as a dictionary
        
        one entry per line with the following syntax :
        entry=value
        
        keys are set to lowercase
        """
        try:
            fin = open(fname)
            f = fin.read()
            ls = f.split("\n")
            for v in ls:
                if ( v == "" or v[0] == "#"):   # skip comments and empty lines
                    pass
                else:
                    (dkey,fval) = parse(v.strip())
                    try:            # check key existence
                        val=self.data[dkey]
                    except:
                        pass
                    # else:
                    #     if val != fval:
                    #         print "WARNING, key -",dkey,"- defined twice in",fname
                    #         print "         previous value :",val," new value",fval
                    self.data[dkey]=fval
            fin.close
        except:
            print("File "+fname+" not found.")
            print("Creating empty parameter list\n")

    #---------------------------------------------------------------------------
    def build_default(self,default_list):
        """ build the default parameter dictionary
        
        used in standard actions,
        returns a dictionary built from the default parmaters (see do_default.py and Param/*)
        and the additional parameters defined in the optionnal p_in_arg overwrite the default values
        
        """
#        from NPK import NPK_PATH
#        base = os.path.join(NPK_PATH,"Param")
        from . import Generic
        base = os.path.join(Generic.get_npk_path(),"Param")
        for l in default_list:
            self.load( os.path.join(base,l) +".gtb" )
#            for i in pp.keys():
#                self[i]=pp[i]
    #---------------------------------------------------------------------------
    def change_key(self,patternOut, patternIn):
        """
        goes though the given dictionnay (which remains unchanged)
        and changes in keys the pattern "patternIn" to "patternOut"
        and returns the modified dictionnary

        typically used in 3D processing :
        
        p_in.change_key('f1', 'f2').change_key('f2', 'f3')

            substitutes F2 (of the 3D) by F1 (of the 2D plane)
            substitutes F3 (of the 3D) by F2 (of the 2D plane)

        """
        import re
        verbose = 0
        if verbose: print("change_key\n", patternIn,patternOut, self)
        p = NPKParam()
        if (len(self)>0):
            for i in self.keys():
                j=re.sub(patternIn,patternOut,i)
                try:
                    p[j]=self.raw(i)        # we want raw() values
                except:
                    p[j]=self[i]
        if verbose:
            print(p)
            raise Exception("stop due to verbose mode")
        return p

# if launched with
# NPK -e Param.py 

if (__name__ == "__main__"):
  import math
  p=NPKParam()
  p["a"]="pa"
  p["b"]="1k"
  p["c"]="math.pi"
  p["d"]="2*3*4"
  for i in p.keys():
      print(i,":",p[i],type(p[i]))
