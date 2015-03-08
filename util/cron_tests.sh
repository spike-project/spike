#!/bin/sh
rm -fR $HOME/tested_npkv2_prev
mv $HOME/tested_npkv2 $HOME/tested_npkv2_prev
cd $HOME/tested_npkv2
/usr/local/bin/hg clone ssh://hg@bitbucket.org/delsuc/SPIKE
python -c "import spike.Tests as t; t.main()"