#!/bin/sh
rm -fR $HOME/tested_spike_prev
mv $HOME/tested_spike $HOME/tested_spike_prev
cd $HOME/tested_spike
/usr/local/bin/hg clone ssh://hg@bitbucket.org/delsuc/spike
python -c "import spike.Tests as t; t.main()"