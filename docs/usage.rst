Usage
=====

Configuration
-------------

Various parameters_ can be set in the configration file. Example parameters for running PyRate with GAMMA
format interferograms are contained in the example configuration file *input_parameters.conf*.

.. _parameters: https://geoscienceaustralia.github.io/PyRate/config.html


Workflow
--------

After following the steps in `Installation <installation.html>`__, an
executable program ``pyrate`` is created.

Use ``help`` for the different command line options:

::

    >> pyrate --help
    Usage: pyrate [OPTIONS] COMMAND [ARGS]...

    Options:
      -v, --verbosity [DEBUG|INFO|WARNING|ERROR]
                                      Level of logging
      --help                          Show this message and exit.

    Commands:
      linrate
      postprocess
      prepifg

The ``pyrate`` program has three command line options corresponding to
different parts of the PyRate workflow:

1. prepifg
2. linrate
3. postprocess

Below we discuss these options.

prepifg: Preparing input interferograms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first step of PyRate is to convert the GAMMA or ROI\_PAC format
unwrapped interferograms into geotiff format, followed by applying
multi-looking and cropping operations. These procedures are all
performed by the ``pyrate prepifg`` command:

::

    >> pyrate prepifg --help
    Usage: pyrate prepifg [OPTIONS] CONFIG_FILE

    Options:
      --help  Show this message and exit.

The ``prepifg`` command is used as follows:

::

    pyrate prepifg /path/to/config_file

The two major steps during the ``prepifg`` operation are described
below.

Data formatting: convert to geotiff
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``prepifg`` command will determine the input format from the value
specified at the *processor:* keyword in the config file (0: ROI\_PAC;
1: GAMMA)

Each GAMMA geocoded unwrapped interferogram requires three header files
to extract metadata required for data formatting: a geocoded DEM header
file (*\*.dem.par*), and the master and slave epoch SLC parameter files
(*\*.slc.par*).

The path and name of the DEM header file are specified in the config
file under the *demHeaderFile:* keyword.

The SLC parameter files should be in the directory specified in the
config file under the *slcFileDir:* keyword. SLC parameter files for a
particular interferogram are found automatically by date-string pattern
matching.

Each ROI\_PAC geocoded unwrapped interferogram requires its own
header/resource file (*\*.unw.rsc*). These header files need to be
stored in the same directory as the interferograms.

In addition, the geocoded DEM header file (*\*.dem.rsc*) is required and
its path and name are specified in the config file under the
*demHeaderFile:* keyword. The geographic projection in the parameter
*DATUM:* is extracted from the DEM header file.

Image transformations: multilooking and cropping
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``prepifg`` command will also perform multi-looking (image
sub-sampling) and cropping of the input interferograms.

Two example configuration files are provided in the *configs/*
directory, one each for ROI\_PAC and GAMMA prepifg configuration. Either
configuration file can be used with ``prepifg``.

linrate: Main workflow and linear rate and time series analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    >> pyrate linrate --help
    Usage: pyrate linrate [OPTIONS] CONFIG_FILE

    Options:
      -r, --rows INTEGER  divide ifgs into this many rows
      -c, --cols INTEGER  divide ifgs into this many columns
      --help              Show this message and exit

This is the core of the PyRate processing workflow, handled by the
``linrate`` command:

::

    pyrate linrate path/to/config_file -c 3 -r 4

This command will perform the time series and linear rate analysis and
has the option to break the interferograms into a number of tiles in
``r`` rows and ``c`` columns. For example, the above command will break
the interferograms into 12 tiles and will produce 12 linear rate and
time series products corresponding to each tile.

The optional rows and columns arguments can be used to create smaller
tiles of the full size interferograms. This enables large interferograms
to be more easily be accommodated in memory. The number of tiles chosen
should be as small as possible that fits in the system memory.

Optionally, an orbital error correction and a spatio-temporal filter
operation to estimate and remove atmospheric phase screen signals is
applied to the interferograms prior to time series and linear rate
analysis. The corrected interferograms are updated on disk and the
corrections are not re-applied upon subsequent runs. This functionality
is controlled by the *orbfit:* and *apsest:* options in the
configuration file.

Non-optional pre-processing steps include: - Minimum Spanning Tree
matrix calculation - Identification of a suitable reference pixel -
Removal of reference phase from interferograms - Calculation of
interferogram covariance - Assembly of the variance-covariance matrix

postprocess: Putting the tiles back together
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The last step of the PyRate workflow is to re-assemble the tiles and
save geotiff files of the final time series and linear rate products.

::

    >> pyrate postprocess --help
    Usage: pyrate postprocess [OPTIONS] CONFIG_FILE

    Options:
      -r, --rows INTEGER  divide ifgs into this many rows
      -c, --cols INTEGER  divide ifgs into this many columns
      --help              Show this message and exit.

Make sure to use the same number of rows and columns that was used in
the previous ``linrate`` step:

::

    pyrate postprocess path/to/config_file -c 3 -r 4

Multiprocessing
---------------

PyRate can use standard multi-threading simply by turning
``parallel:  1`` in the configuration
file to take advantage of multiple cores on a single PC.
