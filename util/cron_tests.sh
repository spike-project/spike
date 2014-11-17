#!/bin/sh

/usr/local/bin/hg clone https://delsuc@bitbucket.org/NPK_v2/draft $HOME/tested_npkv2
cd $HOME/tested_npkv2
$HOME/anaconda/bin/ipython Tests.py
rm -R $HOME/tested_npkv2