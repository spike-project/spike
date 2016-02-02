#!/bin/sh
# This can be used in a cron job to test spike regularly

HG=/usr/local/bin/hg
PYTHON=/Volumes/biak_1ToHD/rdc/anaconda/bin/python
rm -fR $HOME/tested_spike_prev
mv $HOME/tested_spike $HOME/tested_spike_prev
mkdir $HOME/tested_spike
cd $HOME/tested_spike
$HG clone ssh://hg@bitbucket.org/delsuc/spike
$PYTHON -m spike.Tests -m youradress@mail.com -D /Volumes/DATA_test 

# create spikedoc.bitbucket.org
# this account must have read/write access 
$HG clone ssh://hg@bitbucket.org/spikedoc/spikedoc.bitbucket.org

#build doc
sh spike/doc/script_doc.sh 
cd spikedoc.bitbucket.org
# add new files
for i in $($HG stat | awk '/\?/ {print $2}'); do $HG add $i; done
$HG commit -m 'automatic nightly build'
$HG push

# update to devel branch
cd spike
$HG pull && $HG update devel
cd ..
# can run tests here !
$PYTHON -m spike.Tests -m youradress@mail.com -D /Volumes/DATA_test 

