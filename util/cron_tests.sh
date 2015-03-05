#!/bin/sh

/usr/local/bin/hg clone ssh://hg@bitbucket.org/delsuc/SPIKE $HOME/tested_npkv2
cd $HOME/tested_npkv2
$HOME/anaconda/bin/ipython Tests.py
rm -R $HOME/tested_npkv2