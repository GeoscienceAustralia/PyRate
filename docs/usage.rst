Usage
=====

Configuration
-------------

Various parameters_ can be set in the configration file. Example parameters for running PyRate with GAMMA
format interferograms are contained in the example configuration file *input_parameters.conf*.

.. _parameters: https://geoscienceaustralia.github.io/PyRate/config.html

File Discovery
~~~~~~~~~~~~~~

To allow flexibility in the file types the can be processed, PyRate requires
file lists to be provided. This allow PyRate to identify what files are of
which type without relying on file extensions. The path to 
these lists are provided under the following keywords in the configuration 
file:

.. note::
    
    Filenames should be provided without the preceding path, wtih each
    name on a separate line.
    
``ifgfilelist``: this is the list of interferograms to be processed.

.. note::

    Interferogram filenames must contain an epoch pair. Any naming convention 
    is appropriate so long as an epoch pair of format ``XXXXXXXX-YYYYYYYY`` 
    exists in the filename. 

    Example of an interferogram file list:
    ::

        20150702-20150920_interferogram
        20151219-20160109_interferogram
        20160202-20160415_interferogram
    
``slcfilelist``: this is the list which contains the pool of available
GAMMA headers. 

.. note::

    Header filenames must contain an epoch. The epoch must be
    in the format ``XXXXXXXXX``.

    Example of a GAMMA header file list:
    ::

        20150702_header
        20150920_header
        20151219_header
        20160109_header
        20160202_header
        20160415_header
       
``cohfilelist``: this is the list which contains the pool of available
coherence files (used in optional coherence masking). 

.. note::

    Coherence filenames must contain an epoch pair. The epoch pair must be 
    in the format ``XXXXXXX-YYYYYYYY``.

    Example of a coherence file list:
    ::

        20150702-20150920_coherence
        20151219-20160109_coherence
        20160202-20160415_coherence
    
The epochs in filenames are used to match the corresponding header or coherence
files to each interferogram. It is recommended to provide all available headers/coherence 
files in their respective lists, as only the necessary files will be
used. This allows you to process a subset of interferograms by reducing
the names in ``ifgfilelist`` without needing to modify anything else.
    

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
      conv2tif  Convert interferograms to geotiff.
      merge       Reassemble computed tiles and save as geotiffs.
      prepifg           Perform multilooking and cropping on geotiffs.
      process           Time series and linear rate computation.

The ``pyrate`` program has four command line options corresponding to
different parts of the PyRate workflow:

1. conv2tif
2. prepifg
3. process
4. merge

Below we discuss these options.

conv2tif: Converting input intergerograms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before PyRate can process GAMMA or ROI\_PAC intergerograms, they need to be
converted into geotiff format by the ``conv2tif`` command.

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
parameter files (supplied by ``slcfilelist`` in config).

The SLC parameter files should be in the directory specified in the
config file under ``slcFileDir``. SLC files for a
particular interferogram are found automatically by date-string pattern
matching based on epochs. If ``slcFileDir`` is not provided, PyRate will
fallback to looking in the observations direcotry (``obsdir`` in config). 

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
operations to the converted interferograms. 
These procedures are all performed by the ``pyrate prepifg`` command:

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

If specified, ``prepifg`` will perform coherence masking of the unwrapped
interferogram before multilooking and cropping is performed. This requires
corresponding coherence images for each unwrapped interferogram. 

Coherence masking is enabled  by setting the ``cohmask`` argument to ``1`` in
the configuration file. A threshold, ``cohthresh`` needs to be provided. If
``cohfiledir`` is provided, this is where PyRate will look for coherence 
images. If not provided it will look in observations directory where the
unwrapped interferograms exist (``obsdir`` in config). The available coherence 
filenames need tospecified in a file list and provided as the 
``cohfilelist`` parameter.

Image transformations: multilooking and cropping
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``prepifg`` command will perform multi-looking (image
sub-sampling) and cropping of the input interferograms in geotiff format.
The purpose of this is to reduce the resolution of the interferograms to 
reduce the computational complexity of performing the time series and 
linear rate analysis.

An example configuration file is provided in the root source directory
as ``input_parameters.conf``. 

process: Main workflow and linear rate and time series analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    >> pyrate process --help
    Usage: pyrate process [OPTIONS] CONFIG_FILE

    Options:
      -f, --config_file STRING path to configuration file
      -r, --rows INTEGER  divide ifgs into this many rows
      -c, --cols INTEGER  divide ifgs into this many columns
      --help              Show this message and exit

This is the core of the PyRate processing workflow, handled by the
``process`` command:

::

    pyrate process -f path/to/config_file -c 3 -r 4

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
is controlled by the ``orbfit`` and ``apsest`` options in the
configuration file.

Non-optional pre-processing steps include: - Minimum Spanning Tree
matrix calculation - Identification of a suitable reference pixel -
Removal of reference phase from interferograms - Calculation of
interferogram covariance - Assembly of the variance-covariance matrix

merge: Putting the tiles back together
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The last step of the PyRate workflow is to re-assemble the tiles and
save geotiff files of the final time series and linear rate products.

::

    >> pyrate merge --help
    Usage: pyrate merge -f CONFIG_FILE [OPTIONS]

    Options:
      -f, --config_file STRING path to configuration file
      -r, --rows INTEGER  divide ifgs into this many rows
      -c, --cols INTEGER  divide ifgs into this many columns
      --help              Show this message and exit.

Make sure to use the same number of rows and columns that was used in
the previous ``process`` step:

::

    pyrate merge -f path/to/config_file -c 3 -r 4

Multiprocessing
---------------

PyRate can use standard multi-threading simply by turning
``parallel:  1`` in the configuration
file to take advantage of multiple cores on a single PC.
