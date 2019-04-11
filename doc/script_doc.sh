# Start from directory one level above spike dir.
#cd ..
rm -R spike_make_sphinx/*

# get doc from source
sphinx-apidoc -F -o spike_make_sphinx spike
# copy rst files and configuration
cp  doc/*.rst   spike_make_sphinx/
cp  doc/conf.py   spike_make_sphinx/
pandoc README.md -t rst -o Readme.rst
pandoc README.md -s -t latex -o README.pdf
cp Readme.rst spike_make_sphinx/Readme.rst
pandoc release_notes.md -t rst -o spike_make_sphinx/release_notes.rst

# build the doc
make -C spike_make_sphinx html

# modify result to hide *_debug attributes and *_Tests classes
#python2  doc/insert_hide.py spike_make_sphinx/_build/html a'_debug' c'_Tests'

# copy for the spikedoc.org web site
mkdir  -p spikedoc           # should already be there if ou want to upload
cp -r spike_make_sphinx/_build/html/*  spikedoc/

# remove temporary directory
#rm -R spike_make_sphinx/*

#