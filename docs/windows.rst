Windows
-------

These instructions have been tested on Windows 10. PyRate setup requires MPI, Git and QGIS as prerequisites.

Download and install MPI from:

::

    https://www.microsoft.com/en-us/download/details.aspx?id=100593

Download and install Git from:

::

    https://git-scm.com/

Download and install QGIS version 3.10 from:

::

    https://www.qgis.org/en/site/forusers/download.html

Open QGIS terminal from "C:\Program Files (x86)\QGIS 3.10\OSGeo4W". Execute the following commands to clone PyRate
and install required python modules:

::

    cd %UserProfile%\Desktop
    git clone -b test --single-branch https://github.com/GeoscienceAustralia/PyRate
    cd PyRate
    pip install -r requirements-dev.txt
    pip install -r requirements-test.txt
    pip install pyproj==2.2.1
    python setup.py install


QGIS terminal will be required to run PyRate workflow.