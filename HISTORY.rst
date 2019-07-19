.. :changelog:

Release History
---------------
Unreleased (2019-07-xx)
+++++++++++++++++++++++
Added
-----
* `apt-install.sh` in utils/ that lists Ubuntu/apt package requirements.
* `load_modules.sh` in utils/ that sets up Raijin environment.

Fixed
-----
* Errors being caused by newer version of networkx. Nx 2.3 now supported.

Removed
-------
* Unused Python and OS packages.
* environment.yml - conda env will now be installed using 'requirements.txt'.
* HPC directory - hpc README.rst moved to docs.
* setup.cfg - no longer needed.

Changed
-------
* Requirements now managed by 'requirements.txt' file, parsed by setup.py.
* Requirements now split across base 'requirments.txt' and files for dev 
  ('requirements-dev.txt') and testing ('requirements-test.txt')
* Moved default config files to top level source directory.
* Pinned Python dependencies to specific versions.
* Travis build now installs GDAL from apt.
* Travis only builds on master, develop and \*-travis branches.


0.2.0 (2017-05-22)
++++++++++++++++++

* Stable beta release.

0.1.0 (2017-01-31)
++++++++++++++++++

* First release on PyPI.


