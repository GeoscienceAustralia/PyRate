===================
Releasing a version
===================

These are instructions on how to make a `PyRate` Release.

Prepare the release
-------------------

What to merge
^^^^^^^^^^^^^

Decide whether develop will be merged with master, or some other release branch forking off develop will be merged with master. 
For a small release with not much active development by multiple people, merging develop directly in to master is a reasonable approach. 

Level of release
^^^^^^^^^^^^^^^^

Decide on level of release (major, semi-major, minor) and increment version number accordingly.

Read up on `semantic versioning`_ if you are unfamiliar.

.. _semantic versioning: https://semver.org/ 

Update the variable __version__ accordingly in this file_ and commit change.

.. _file: https://github.com/GeoscienceAustralia/PyRate/blob/develop/setup.py

Update Documentation
^^^^^^^^^^^^^^^^^^^^

Review, update and commit changes to the documentation based on changes in the release including the authors_ and changelog_: 

.. _authors: https://github.com/GeoscienceAustralia/PyRate/blob/master/docs/authors.rst
.. _changelog: https://github.com/GeoscienceAustralia/PyRate/blob/master/docs/history.rst

* Edit .rst files in docs/ dir.  

* pip install -r requirements-dev.txt 

* cd into pyrate/docs and run make html to build the html pages,  

* copy the new _build dir to Windows 

* review html files in a browser. 

If the year has changed: update the year in all of the code files that contain this copyright line, and commit changes:

#   Copyright 2021 Geoscience Australia 

This can be achieved using this command line code from the main PyRate folder, replace the years appropriately: 

   `find . -type f -exec sed -i 's/2021\(.*Geoscience*\)/2022\1/' {} +`


Review and publish on GitHub 
----------------------------

* Lodge a `Pull Request`_ (PR) that will merge to master (draft). Add the changelog in to the details of the PR. Seek peer review from appropriate team members, does the software function as expected? Has anything broken?

* `Draft a new release`_ in GitHub. Add the changelog content to the release notes.

* After addressing any review comments in the PR and after all CI tests are passing, merge the PR. 

* Publish the draft release in GitHub and point it to the tip of master (this will tag that commit with the release number). 

* Check that the html link to tar.gz Source Code file from github release page matches download_url in setup.py. Commit if you change it. 

* Merging to master should trigger Actions CI to build and deploy the web documentation. Monitor CI for its progress, then check the web pages have been updated_. 

.. _Pull Request: https://github.com/GeoscienceAustralia/PyRate/pulls

.. _updated: https://geoscienceaustralia.github.io/PyRate/

.. _Draft a new release: https://github.com/GeoscienceAustralia/PyRate/releases

Publish release to PyPi
-----------------------

The following steps will prepare the package and upload it to the Python Package Index. The final upload step requires log-in credentials, contact the core development team for assistance with this.

`pip install -r requirements-dev.txt – install dependencies`

`python3 setup.py sdist bdist_wheel - prepare package`

`python3 -m twine upload dist/*`

To check it worked visit the packages PyPi page_

.. _page: https://pypi.org/project/Py-Rate/


Following this `PyRate`can be installed using the command `pip install Py-Rate`.
