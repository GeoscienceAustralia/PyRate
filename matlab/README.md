Creating and running Matlab executables
=======================================
It is useful to have a copy of the original Matlab version of Pirate (from which PyRate was ported) for testing purposes.

However, it may not be possible for all users to run Pirate because they don't have a Matlab license. Fortunately, Matlab offers the ability to create native "deployed Matlab" executables that simply require the appropriate "Matlab Runtime" (http://au.mathworks.com/products/compiler/mcr/) and this has been done for Pirate and all the necessary files are contained in this folder.

What's in this folder?
----------------------
* `pirate_v3.2beta1_modified.zip` is the original source code that was used to create PyRate
	* this is the version that came from Sudipta Basak's dropbox
	* the code has been slighly modified to
		* work better when being called from the command line
		* store all output as a single `.mat` file as opposed to individual files
			* the structure of the `.mat` file can be seen in PyRate/pyrate/utils/cmp_out.py
	* code changes have been tracked with git
* pirate folder
	* this folder contains the deployed Matlab executables created from running the command `>> mcc -mv pirate.m` from the 64 bit linux version of Matlab R2015a
	* the user must install the 64 bit Linux Matlab Runtimes from http://au.mathworks.com/products/compiler/mcr/
	* the executable can then be run with `$./pirate pirate_config_file.cfg`


