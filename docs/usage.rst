Usage
=====

Configuration
-------------

Various parameters_ can be set in the configration file. Example parameters for running PyRate with GAMMA
format interferograms are contained in the example configuration file *input_parameters.conf*.

.. _parameters: https://geoscienceaustralia.github.io/PyRate/config.html

File Discovery
~~~~~~~~~~~~~~

To allow flexibility in the file types that can be processed, PyRate requires
file lists to be provided. This allows PyRate to identify files of each
type without relying on file extensions. The file path to these lists are 
provided under the following keywords in the configuration file:

.. note::

    Filenames should be provided with their absolute system path, with each
    name on a separate line.

``ifgfilelist``: this is the list of interferograms to be processed.

.. note::

    Interferogram filenames must contain a pair of epochs. Any naming convention
    is appropriate so long as an epoch pair of format ``YYYYMMDD-YYYYMMDD``
    exists in the filename.

    Example of an interferogram file list:
    ::

        /absolute/path/to/20150702-20150920_interferogram
        /absolute/path/to/20151219-20160109_interferogram
        /absolute/path/to/20160202-20160415_interferogram

``slcfilelist``: this is the list which contains the pool of available
GAMMA SLC header files.

.. note::

    Header filenames must contain an epoch. The epoch must be
    in the format ``YYYYMMDD``.

    Example of a GAMMA SLC header file list:
    ::

        /absolute/path/to/20150702_header
        /absolute/path/to/20150920_header
        /absolute/path/to/20151219_header
        /absolute/path/to/20160109_header
        /absolute/path/to/20160202_header
        /absolute/path/to/20160415_header

``cohfilelist``: this is the list which contains the pool of available
coherence files (used during optional coherence masking).

.. note::

    Coherence filenames must contain a pair of epochs. The epoch pair must be
    in the format ``YYYYMMDD-YYYYMMDD``.

    Example of a coherence file list:
    ::

        /absolute/path/to/20150702-20150920_coherence
        /absolute/path/to/20151219-20160109_coherence
        /absolute/path/to/20160202-20160415_coherence

The epochs in filenames are used to match the corresponding header or coherence
files to each interferogram. It is recommended to provide all available headers/coherence
files in their respective lists, as only the necessary files will be
used. This allows you to process a subset of interferograms by removing
entries in ``ifgfilelist`` without needing to modify anything else.


Workflow
--------

After following the steps in `Installation <installation.html>`__, an
executable program ``pyrate`` is created.

Use ``--help`` for the different command line options:

::

    Usage: pyrate [OPTIONS] COMMAND [ARGS]...

        CLI for carrying out PyRate workflow. Typical workflow:

            Step 1: conv2tif

            Step 2: prepifg

            Step 3: process

            Step 4: merge

        Refer to https://geoscienceaustralia.github.io/PyRate/usage.html for more
        details.

    Options:
      -V, --version                   Show the version and exit.
      -v, --verbosity [DEBUG | INFO | WARNING | ERROR]
                                      Level of logging
      --help                          Show this message and exit.

    Commands:
      conv2tif    Convert interferograms to geotiff.
      merge             Reassemble computed tiles and save as geotiffs.
      prepifg           Perform multilooking and cropping on geotiffs.
      process           Time series and Stacking computation.

The ``pyrate`` program has four command line options corresponding to
different steps in the PyRate workflow:

1. conv2tif
2. prepifg
3. process
4. merge

Below we discuss these steps.

conv2tif: Converting input interferograms to Geotiff format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before PyRate can process GAMMA or ROI\_PAC interferograms, they need to be
converted into geotiff format using the ``conv2tif`` command.

::

    >> pyrate conv2tif --help
    Usage: pyrate conv2tif -f CONFIG_FILE

      Convert interferograms to geotiff.

    Options:
      --help  Show this message and exit.

The ``conv2tif`` command will determine the input format from the value
specified at the *processor:* keyword in the config file (0: ROI\_PAC;
1: GAMMA)

Each GAMMA geocoded unwrapped interferogram requires three header files
to extract metadata required for data formatting: a geocoded DEM header
file (``demHeaderFile`` in config) and the master and slave epoch SLC
header files (supplied by ``slcfilelist`` in config).

SLC files for a particular interferogram are found automatically by
date-string pattern matching based on epochs. 

Each ROI\_PAC geocoded unwrapped interferogram requires its own
header/resource file. These header files need to be
stored in the same directory as the interferograms.

In addition, the geocoded DEM header file is required and
its path and name are specified in the config file under ``demHeaderFile``.
The geographic projection in the parameter *DATUM:* is extracted from the DEM
header file.

Upon completion, geotiff formatted copies of the input files will be placed
in the directory the input files are located in. Note that ``conv2tif``
will not perform the conversion if geotiffs for the provided input files
already exist.

prepifg: Preparing input interferograms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The second step of PyRate is applying multi-looking and cropping
operations to the geotiff interferograms.
These procedures are all performed by the ``prepifg`` command:

::

    >> pyrate prepifg --help
    Usage: pyrate prepifg  -f CONFIG_FILE [OPTIONS]

    Options:
      --help  Show this message and exit.

The ``prepifg`` command is used as follows:

::

    pyrate prepifg -f /path/to/config_file

Coherence masking
^^^^^^^^^^^^^^^^^

If specified, ``prepifg`` will perform coherence masking on the
interferograms before multilooking and cropping is performed. This requires
corresponding coherence images for each interferogram. The purpose
of coherence masking is to remove poor quality phase observations and leave a set of
high-quality pixels for analysis. Pixels with coherence values below a certain threshold
will be set to Not-a-Number (NaN). Note that the number of pixels with numeric phase values
(i.e. pixels not equal to NaN) in each interferogram will be different after coherence masking.

Coherence masking is enabled by setting the ``cohmask`` parameter to ``1`` in
the configuration file. A threshold, ``cohthresh`` needs to be provided. If
``cohfiledir`` is provided, this is where PyRate will look for coherence
images. If not provided it will look in the observations directory where the
interferograms exist (``obsdir`` in config). The available coherence
filenames need to be specified in a file list and provided as the
``cohfilelist`` parameter.

Image transformations: multilooking and cropping
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``prepifg`` step will perform multi-looking (image
sub-sampling) and cropping of the input interferograms in geotiff format.
The purpose of this is to reduce the resolution of the interferograms to
reduce the computational complexity of performing the time series and
stacking analysis.

An example configuration file is provided in the root source directory
as ``input_parameters.conf``.

process: Main workflow, including stacking and time series analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    >> pyrate process --help
    Usage: pyrate process [OPTIONS] CONFIG_FILE

    Options:
      -f,   --config_file   STRING   Path to configuration file
      -r,   --rows          INTEGER  Divide interferograms into this many rows
      -c,   --cols          INTEGER  Divide interferograms into this many columns
      --help                         Show this message and exit

This is the core of the PyRate processing workflow, handled by the
``process`` command:

::

    pyrate process -f path/to/config_file -c 3 -r 4

This command will perform the time series and stacking analysis and
has the option to break the interferograms into a number of tiles in
``r`` rows and ``c`` columns. For example, the above command will break
the interferograms into 12 tiles and will produce 12 stacking and
time series products corresponding to each tile.

The optional rows and columns arguments can be used to split the full-size
interferograms into smaller tiles. This enables large interferograms
to be more easily accommodated in memory. The number of tiles chosen
should be as small as possible that fits within the available system memory.

Optionally, an orbital error correction and a spatio-temporal filter
operation to estimate and remove atmospheric phase screen (APS) signals is
applied to the interferograms prior to time series and stacking
analysis. The corrected interferograms are updated on disk and the
corrections are not re-applied upon subsequent runs. This functionality
is controlled by the ``orbfit`` and ``apsest`` options in the
configuration file.

Non-optional pre-processing steps include: - Minimum Spanning Tree
matrix calculation - Identification of a suitable reference pixel -
Removal of reference phase from interferograms - Calculation of
interferogram covariance - Assembly of the variance-covariance matrix

merge: Putting the tiles back together
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The last step of the PyRate workflow is to re-assemble the tiles and
save geotiff files of the final time series and stacking products.

::

    >> pyrate merge --help
    Usage: pyrate merge -f CONFIG_FILE [OPTIONS]

    Options:
      -f,   --config_file   STRING   Path to configuration file
      -r,   --rows          INTEGER  Divide interferograms into this many rows
      -c,   --cols          INTEGER  Divide interferograms into this many columns
      --help                         Show this message and exit

Make sure to use the same number of rows and columns that was used in
the previous ``process`` step:

::

    pyrate merge -f path/to/config_file -c 3 -r 4

Parallel processing
-------------------

PyRate can be run in parallel using standard multi-threading simply by turning
``parallel:  1`` in the configuration file to take advantage of multiple cores
on a single machine. The parameter ``processes`` sets the number of threads.

Alternatively, PyRate can be parallelised using MPI by using ``mpirun``:

::

    # Modify '-n' based on the number of processors available.
    mpirun -n 4 pyrate conv2tif -f input_parameters.conf
    mpirun -n 4 pyrate prepifg -f input_parameters.conf
    mpirun -n 4 pyrate process -f input_parameters.conf -c 2 -r 2
    mpirun -n 4 pyrate merge -f input_parameters.conf -c 2 -r 2

In the case that PyRate is run using MPI, standard multi-threading is automatically
disabled (i.e. equivalent to setting ``parallel:  0``).

