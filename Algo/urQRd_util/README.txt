Applic_urqrd
== Description ==
Applic_urqrd.py is used for performing denoising on harmonic signals. 
Stable version: 0.1

== Installation ==
For running the code you have to install Canopy https://www.enthought.com/products/canopy/ 
or Anaconda https://store.continuum.io/cshop/anaconda/
In the config file "Applic_urqrd.mscf" enter the path to the main package "draft_Fista".
i.e : npk_dir = path to draft_urqrd
"draft_Fista" contains the directory urqrd_applic, util, v1 and other dependencies. 
In /draft_urqrd/Applic_urqrd, copy the file "Applic_urqrd_eg.mscf" and change its name to "Applic_urqrd.mscf"  (no "_eg").


== Configuration ==
Before using urqrd you need to change the config file "Applic_urqrd.mscf". 
Four sections need to be changed : config, data, fista and store.
The first one is for specifying the path to the main package (see here above).
The following ones are :
= data =
	
= fista =
	
== Running ==
To run the program once the installation and configuration finished, in Canopy browse to the file /draft_urqrd/urqrd_applic/Applic_urqrd.py and execute. 
In Anaconda launch Python from /bin/python2.7 i.e write : Anaconda/bin/python2.7  /draft_urqrd/urqrd_applic/Applic_urqrd.py
	
== Changelog ==

= 0.1 =



