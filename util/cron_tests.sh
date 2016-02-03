#!/bin/sh
# This can be used in a cron job to test spike regularly

# adpt those to your own set-up
# the commands
HG=/usr/local/bin/hg
PYTHON=/usr/local/anaconda/bin/python
# your email
EMAIL=youradress@mail.com
#location of DATA files for tests
DATAdir=/Volumes/DATA_test

# then prepare directories
rm -fR $HOME/tested_spike_prev
mv $HOME/tested_spike $HOME/tested_spike_prev
mkdir $HOME/tested_spike
cd $HOME/tested_spike

# clone code from repository
$HG clone ssh://hg@bitbucket.org/delsuc/spike

# and run tess
$PYTHON -m spike.Tests -m $youradress@mail.com -D $DATAdir

# Now build the documentation
# this account must have read/write access to the repository for the doc be push backed
# get the code
$HG clone ssh://hg@bitbucket.org/spikedoc/spikedoc.bitbucket.org

#build doc
sh spike/doc/script_doc.sh 
cd spikedoc.bitbucket.org

# add new files to repository
for i in $($HG stat | awk '/\?/ {print $2}'); do $HG add $i; done
$HG commit -m 'automatic nightly build'
$HG push
cd ..

# Eventually:
# update to devel branch
cd spike
$HG pull && $HG update devel
cd ..
# and run tests here !
$PYTHON -m spike.Tests -m $youradress@mail.com -D $DATAdir

