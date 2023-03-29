# Start from directory one level above spike dir.
#cd ..
rm -R spike_make_sphinx/* dist/* build/*

# get doc from source
Version=$(grep "version = '.*'" spike/version.py |grep -o [0-9\.]*)

# copy rst files
cp  doc/conf.py spike_make_sphinx/
cp  doc/*.rst   spike_make_sphinx/
cp  doc/sphinxMakefile   spike_make_sphinx/Makefile
# modiy configuration

pandoc README.md -t rst -o spike_make_sphinx/Readme.rst
pandoc release_notes.md -t rst -o spike_make_sphinx/release_notes.rst

# still to do Presentation
# and         DevelopmentGuide

#pandoc release_notes.md -s -t latex -o release_notes.pdf

# build the doc
make -C spike_make_sphinx html

# modify result to hide *_debug attributes and *_Tests classes
#python2  doc/insert_hide.py spike_make_sphinx/_build/html a'_debug' c'_Tests'

# copy for the spikedoc.org web site
mkdir  -p spikedoc           # should already be there if ou want to upload
cp -r spike_make_sphinx/_build/html/*  spikedoc/

# remove temporary directory
#rm -R spike_make_sphinx/*

#sphinx-apidoc --full -V $Version -A M-A.Delsuc -o spike_make_sphinx_sep spike
