#!/usr/bin/env python 
# encoding: utf-8

'''
Program for controlling if files are not deteriorated during transmission.
Usage : write in command line "python hashcheck.py directory-pathname"
Lionel Chiron, 11 oct 2012
copyright NMRTEC S.A.S.
'''

from __future__ import print_function

import hashlib
import fnmatch
import os,sys
import json
import unittest

def md5sum(filename):
    '''
    from filename make a hash name.
    able to hash big files
    '''
    md5 = hashlib.md5()
    chunksize = 128*md5.block_size
    with open(filename,'rb') as f: 
        for chunk in iter(lambda: f.read(chunksize), b''): 
             md5.update(chunk)
    return md5.hexdigest()

def hashdir(rootpath, pattern = '*' ):
    '''
    from directory rootpath make a dictionary
    whose entries are the hash numbers associated to address
    of the files
    '''
    dicf = {}
    for root, dirs, files in os.walk(rootpath):
        for filename in fnmatch.filter(files, pattern):
            key = os.path.join(root, filename)
            hashnb = md5sum(key)
            dicf[hashnb] = key    
    return dicf

def hash_make_or_control(rootpath = ".", verbose = False):
    '''
    Makes the file hashfile in the directory hash if it doesn't exists.
    If the hash directory exists, check if files are ok or not.
    '''
    if rootpath is None:
        rootpath = sys.argv[1] # taking the argument from the command line.
    dirhash = os.path.join(rootpath,'hash')
    if not os.path.exists(dirhash ): # hash directory created
        os.makedirs(dirhash )
        hashres = hashdir(rootpath)
        f = open(os.path.join(dirhash,'hashfile'),'wb')
        json.dump(hashres,f)
        f.close()
        print("%s created"%os.path.join(dirhash,'hashfile'))
    else : # compare to hashed files
        allok = True
        addhash = os.path.join(dirhash,'hashfile')
        hash1 = json.load(open(addhash,'rb')) #
        hash2 = hashdir(rootpath)
        for key in hash1 :
            if key in hash2 :
                if verbose: print("ok for %s : %s"%(hash2[key], key))
            else :
                print("error in %s : %s"%(hash1[key], key))
                allok = False
        if allok:
            print("hashcheck : Ok, directory not corrupted or modified.")
        else:
            print("This directory has been corrupted or modified")

class Test_hashcheck(unittest.TestCase):
    def test_hash(self):
        from ..Tests import directory
        hash_make_or_control(rootpath = directory() )
        
if __name__ == '__main__':
    hash_make_or_control()
### Test ####

