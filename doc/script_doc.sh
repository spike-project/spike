# Start from directory one level above spike dir.

# get doc from source
sphinx-apidoc -F -o spike_make_sphinx spike
# copy rst files and configuration
cp  spike/doc/*.rst   spike_make_sphinx/
cp  spike/doc/conf.py   spike_make_sphinx/

# build the doc
make -C spike_make_sphinx html

# modify result to hide *_debug attributes and *_Tests classes
python -m spike.doc.insert_hide spike_make_sphinx/_build/html a'_debug' c'_Tests'

# copy for the spikedoc.bitbucket.org web site
mkdir  -p spikedoc.bitbucket.org           # should already be there if ou want to upload
cp -r spike_make_sphinx/_build/html/*  spikedoc.bitbucket.org/

# remove temporary directory
rm -R spike_make_sphinx/*
