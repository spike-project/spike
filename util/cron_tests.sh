#!/bin/sh
# This can be used in a cron job to test spike regularly
rm -fR $HOME/tested_spike_prev
mv $HOME/tested_spike $HOME/tested_spike_prev
mkdir $HOME/tested_spike
cd $HOME/tested_spike
/usr/local/bin/hg clone ssh://hg@bitbucket.org/delsuc/spike
python -m spike.Tests -m youradress@mail.com -D /Volumes/DATA_test 