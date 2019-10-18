.. :changelog:

Release History
===============
0.4.0 (2019-09-24)
-----------------------
Added
+++++
- Python 3.7 support.
- Optional 'conv2tif' step.
- Building of docs integrated with Travis CI.
- Coherence masking, view coherence masking section in 'input_parameters.conf'
  for options.
- Parameter validation.
- SLC and coherence file lists for file discovery.

Changed
+++++++
- 'linrate' command has been renamed to 'process'.
- CLI structure, config files now need to provided with '-f' flag.
- Reduced console output, default setting is now INFO.
- Restructure of code layout, src modules now in 'core' directory and scripts
  at top level.

Removed
+++++++
- Unused luigi code.
- References to Matlab.
- Python 3.6 support.

0.3.0 (2019-07-26)
-----------------------
Added
+++++
- `apt_install.sh` in utils/ that lists Ubuntu/apt package requirements.
- `load_modules.sh` in utils/ that sets up Raijin environment.

Fixed
+++++
- Errors being caused by newer version of networkx. Nx 2.3 now supported.

Removed
+++++++
- Unused Python and OS packages.
- environment.yml - conda env will now be installed using 'requirements.txt'.
- HPC directory - hpc README.rst moved to docs.
- setup.cfg - no longer needed.
- Luigi functionality - hasn't been operational and is reported as vulnerable.
  Single machine parallelism is achieved with joblib. 

Changed
+++++++
- Requirements now managed by 'requirements.txt' file, parsed by setup.py.
- Requirements now split across base 'requirements.txt' and files for dev 
  ('requirements-dev.txt') and testing ('requirements-test.txt')
- Moved default config files to top level source directory.
- Pinned Python dependencies to specific versions.
- Travis build now installs GDAL from apt.
- Travis only builds on master, develop and \*-travis branches.
- Consolidated documentation into /docs.
- Updated install instructions for Ubuntu and NCI.

0.2.0 (2017-05-22)
------------------
- Stable beta release.

0.1.0 (2017-01-31)
------------------
- First release on PyPI.
