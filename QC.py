#!/usr/bin/env python 
# encoding: utf-8

"""
run_pylint.py

Created by Marc-André on 2011-10-26.
Copyright (c) 2011 IGBMC. All rights reserved.

python 2 only so far !

"""

from __future__ import print_function, division
from subprocess import Popen, PIPE

import sys
import re
from codecs import encode, decode
#import glob

CMD = "pylint"
ARGS = ["--rcfile=rcfile3"]

def msg(string, sep='='):
    "print string in a box"
    print(sep*(len(string)+4))
    print('|', string, '|')
    print(sep*(len(string)+4))

def tracked_files(excluded=('none',)):
    """return list of hg tracked files in folder as a tuple : (all,modified) 
    excluded is a list of patterns (unix style) which are excluded
    """
    import fnmatch
    hgall = Popen(["hg", "stat", "-A"], stdout=PIPE).communicate()[0]
    lines = [l for l in hgall.split(b"\n")]
    hg = []
    modif = []
    for line in lines:
        sp = line.split()
        cont = False
        if len (sp) < 1:      # exit on conditions
            continue
        if sp[0] not in (b'C', b'M'):
            continue
        fname = decode(sp[1])
        for pat in excluded:
            if fnmatch.fnmatch(fname, pat): # find pattern
                print("excluding", fname)
                cont = True
        if cont:
            continue
        hg.append(fname)
        if sp[0] in (b'M'):
            modif.append(decode(sp[1]))
    return (hg, modif)
class WritableObject(object):
    "dummy input/output stream for pylint"
    def __init__(self):
        self.content = []
    def write(self, string):
        "dummy write"
        self.content.append(string)
    def readlines(self):
        "dummy read"
        return self.content
    def readline(self):
        "dummy read"
        return self.content.pop(0)
    def __iter__(self):
        return self
    def next(self):
        if self.content:
            return self.content.pop(0)
        else:
            raise(StopIteration)
    def __next__(self):
        return self.next()
def find_depedencies(reader):
    "goes thru the pylint report reader, gather dependencies, and return as a set"
    dep = set()
    for l in reader:
        if l == "Duplication":
            break
        m = re.search(r"^\s*(\w+)", l)
        if m:
            print(m.group(1))
            dep.add(m.group(1))
    print(dep)
    return dep

def run_pylint(fich):
    "run pylint on the given file and return synthetic results (note, error, reorder, warning)"
    from pylint import lint
    from pylint.reporters.text import TextReporter
#    pyl = lint.Run(args)
    pylint_output = WritableObject()
    lint.Run([fich]+ARGS, reporter=TextReporter(pylint_output), exit=False)
    ignored_error = ('Exter', 'E1103')  # 5 first char of ignored errors
    ignored_warn = ('none',)            # 5 first char of ignored warnings
    errors = 0
    warn = 0
    reorder = 0
    note = 0.0
    for l in pylint_output:
        if l.startswith('E0001'):
            errors = 10000  # indicates crash
            break
        elif l.startswith('E'):
            if l[0:5] not in ignored_error:
                errors += 1
        elif l.startswith('W'):
            if l[0:5] not in ignored_warn:
                warn += 1
        elif l.startswith('R') and l.split()[0] != 'Raw':
            reorder += 1
        elif l.startswith('Your code has been rated'):
            m = re.match(r'Your code has been rated at (-?\d+.\d.)/10', l)
            note = float(m.group(1))
        # elif l.find('External dependencies') !=-1:
        #     print "##############"
        #     dep = find_depedencies(pylint_output)
        #     print "##############"
    if errors == 10000:    # means crash
        note = -100
#    print(fich, note, errors, warn, reorder)
    return (note, errors, warn, reorder )

def process(files):
    "apply pylint to files"
    # msg("Checking files")
    stats = {}
    for f in files:
        if f.endswith('.py'):
            print("checking %s"%f)
            res = run_pylint(f)
            stats[f] = res
        else:
            print("excluding", f)
    return stats
def processmp(files):
    "apply pylint to files in a multiprocessing manner"
    import multiprocessing as mproc
    pool = mproc.Pool()
    # msg("Checking files")
    stats = {}
    todo = [ f for f in files if f.endswith('.py')]
    res = pool.imap(run_pylint, todo)
    for f in todo:
        print("checked %s"%f)
        stats[f] = res.next()
    return stats
def report(stats, modif):
    "print a tabular report"
    text = []
    title = "  Note/10   File                          Error Warning"
    nchar = len(title)
    text.append("\n"+"-"*nchar)
    text.append( title )
    text.append( "-"*nchar )
    top = 1
    mean = 0.0
    ptop = 0
    pmed = 0
    plow = 0
    errtot = 0
    kstats = sorted(list(stats.keys()), key = lambda x: stats[x][0])
    for k in reversed(kstats):
        note = stats[k][0]
        err = stats[k][1]
        errtot += err
        warn = stats[k][2]
        if k in modif:
            head = 'M'
        else:
            head = ' '

        if err != 10000:    # means crash
            mean += note

        if note >= 5.0:     # count program per class
            ptop += 1
        elif note >= 0.0:
            pmed += 1
        else:
            plow += 1

        if top == 1:    # draw separators
            if note < 5.0:
                text.append( "-"*nchar )
                top = 0
        elif top == 0:
            if note < 0.01:
                text.append( "*"*nchar )
                top = -1

        if err == 10000:
            text.append( "{} ******  {:35s} HAS SYNTAX ERRORS".format(head, k) )
        elif err > 0:
            text.append( "{} {:5.2f}  {:35s} {:3d}   {:3d} ".format(head, note, k, err, warn) )
        else:
            text.append( "{}{:6.2f}  {:35s}       {:3d}".format(head, note, k, warn) )
    ncode = ptop+pmed+plow
#    assert(len(stats)==ncode)
    text.append( "-"*nchar )
    mnote = mean/ncode
    median = stats[kstats[ncode//2]][0]
    aggreg = 0.5*mnote + 0.5*(10*ptop + 5*pmed)/ncode
    text.append( "\nOverall aggregated note value is : {:5.2f}/10\n".format(aggreg) )
    text.append( "score is {} / {} / {}  (top / med / low)".format(ptop, pmed, plow ))
    text.append( "Mean - Median note values are : {:5.2f}/10 - {:5.2f}/10".format(mnote, median ))
    if median <5.0:
        text.append( "  Median is below 5.0 - let's try to improve it")
    text.append( "Total number of errors : {:d} on a total of {:d} files".format(errtot, ncode) )
    text.append( "-"*nchar )
    return "\n".join(text)

def message():
    "print pylint command"
    print("    to run pylint, type :")
    print(CMD, end=' ') 
    for a in ARGS:
        print(a, end=' ')
    print("filename.py")
    print("    to check only errors, type :")
    print(CMD, "-E", end=' ')
    for a in ARGS:
        print(a, end=' ')
    print("filename.py")
    
def main(excluded = ['.hgignore'], files = 'hg'):
    """does the work
    excluded is a list of excluded file pattern
    files is the list of files to analyse, if files == 'hg' (defaullt) then files are taken from a mercurial repository.
    """
    #collect
    if files == 'hg':
        (hg, modif) = tracked_files(excluded=excluded)
    else:
        modif = []
        hg = files
    if modif:
        msg("Warning, the following files are modified but not commited", sep="*")
        for i in modif:
            print(i)
    stats = processmp(hg)
    t = report(stats, modif)
    print(t)
    F = open("QC.txt","w")
    F.write(t)
    F.close()
    message()

if __name__ == '__main__':
    sys.exit(main(excluded=('v1/*', 'util/dynsubplot.py', 'Miscellaneous/*', 'Visu/*', 'SPIKE_usage_eg/previous-to-clean/*')))
    # excluding dynsubplot.py because the add_class method is just too confusing to pylint
#    sys.exit( main( files=glob.glob('*.py')+glob.glob('*/*.py') ) )

