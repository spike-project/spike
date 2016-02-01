sphinx-apidoc -F -o spike_make_sphinx spike
cp -r spike/doc/*   spike_make_sphinx/
cd spike_make_sphinx
make html
cd ..
python -m spike.doc.insert_hide spike_make_sphinx/_build/html a'_debug' c'_Tests'
cp -r spike_make_sphinx/_build/html/*  spikedoc.bitbucket.org/
rm -R spike_make_sphinx