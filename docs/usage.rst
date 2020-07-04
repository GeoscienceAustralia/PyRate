Usage
=====

Configuration
-------------

Various parameters_ are controlled by values given in a text-based configuration file.
The user can choose any filename for this configuration file.
Example parameters for running `PyRate` with `GAMMA`-format interferograms are
contained in the example configuration file ``PyRate/input_parameters.conf``.

.. _parameters: https://geoscienceaustralia.github.io/PyRate/config.html


Workflow
--------

File Discovery
^^^^^^^^^^^^^^

To allow flexibility in the file types that can be processed, `PyRate` requires
file lists to be provided. This allows `PyRate` to identify files of each
type without relying on file extensions. The file path to these lists are 
provided under ``ifgfilelist``, ``hdrfilelist`` and ``cohfilelist`` keywords
in the configuration file.

.. note::

    - Filenames should be provided in these lists with their absolute system path,
      with each name on a separate line.
    - Filenames should only have a single period (".") delimiting the file extension.
      Multiple periods in a filename will raise an error during runtime.

The ``ifgfilelist`` gives the list of interferograms to be processed.
Example of an interferogram file list for `GAMMA` flat-binary files:

::

    /absolute/path/to/20150702-20150920.unw
    /absolute/path/to/20151219-20160109.unw
    /absolute/path/to/20160202-20160415.unw


.. note::

    - Interferogram filenames must contain a pair of date epochs.
      Any naming convention is appropriate as long as an epoch pair of format
      ``YYYYMMDD-YYYYMMDD`` exists in the filename.

The ``hdrfilelist`` gives the list of `GAMMA` SLC header text files.
Example of a `GAMMA` SLC header file list:

::

    /absolute/path/to/20150702_slc.par
    /absolute/path/to/20150920_slc.par
    /absolute/path/to/20151219_slc.par
    /absolute/path/to/20160109_slc.par
    /absolute/path/to/20160202_slc.par
    /absolute/path/to/20160415_slc.par


.. note::

    - SLC header filenames must contain a single date epoch in the format ``YYYYMMDD``.

The ``cohfilelist`` is an optional list which contains the pool of all available
coherence files to be used for optional coherence masking.
Example of a coherence file list for `GAMMA` flat-binary files:

::

    /absolute/path/to/20150702-20150920.coh
    /absolute/path/to/20151219-20160109.coh
    /absolute/path/to/20160202-20160415.coh

.. note::

    - Like interferograms, coherence filenames must contain a pair of epochs.
      The epoch pair must be in the format ``YYYYMMDD-YYYYMMDD``. Otherwise, any
      naming convention is appropriate.

The date epochs in filenames are used to match the corresponding SLC header
or coherence files to each interferogram. It is recommended to provide a complete
list of available headers/coherence files for a stack in their respective lists,
since only the necessary files will be used. This allows you to process a subset
of interferograms by removing entries in ``ifgfilelist`` without needing to modify
all three lists.

`PyRate` Workflow
^^^^^^^^^^^^^^^^^

After `Installation <installation.html>`__, an
executable program ``pyrate`` is created in the system path:

::

    >> pyrate --help
    usage: pyrate [-h] [-v {DEBUG,INFO,WARNING,ERROR}]
                  {conv2tif,prepifg,process,merge,workflow} ...

    PyRate workflow:

        Step 1: conv2tif
        Step 2: prepifg
        Step 3: process
        Step 4: merge

    Refer to https://geoscienceaustralia.github.io/PyRate/usage.html for
    more details.

    positional arguments:
      {conv2tif,prepifg,process,merge,workflow}
        conv2tif            Convert interferograms to geotiff.
        prepifg             Perform multilooking and cropping on geotiffs.
        process             Main processing workflow including corrections, time series and stacking computation.
        merge               Reassemble computed tiles and save as geotiffs.
        workflow            Run all the PyRate processes

    optional arguments:
      -h, --help            show this help message and exit
      -v {DEBUG,INFO,WARNING,ERROR}, --verbosity {DEBUG,INFO,WARNING,ERROR}
                            Increase output verbosity

 .. note::

    - If running on NCI, be sure to first load the correct modules and virtual environment:
      ``source ~/PyRate/scripts/nci_load_modules.sh`` 

The ``pyrate`` program has four command line options corresponding to
different steps in the `PyRate` workflow:

1. ``conv2tif``
2. ``prepifg``
3. ``process``
4. ``merge``

All four steps are required. 
A fifth option, ``workflow``, is available that will run all four steps in order.
In the following sub-sections we discuss each of the available steps.



``conv2tif``: Converting flat-binary interferograms to Geotiff format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before `PyRate` can process interferograms in flat-binary file format, they need to be
converted into geotiff format using the ``conv2tif`` step.

::

    >> pyrate conv2tif --help
    usage: pyrate conv2tif [-h] -f CONFIG_FILE

    optional arguments:
      -h, --help            show this help message and exit
      -f CONFIG_FILE, --config_file CONFIG_FILE
                            Pass configuration file

The ``conv2tif`` step is used as follows:

::

    >> pyrate conv2tif -f /path/to/config_file

.. note::

    - Users of the `GAMMA` software can skip the ``conv2tif`` step if they have generated
      geotiffs using the `GAMMA` program ``data2geotiff``, which is included in all
      `GAMMA` software distributions.
    - In this case, ``ifgfilelist`` and ``cohfilelist`` would contain the absolute
      paths to these geotiff files. Even when using geotiff files, the SLC header files
      are still required by ``prepifg``.
    - If a DEM is to be processed by ``prepifg``, it's file format should match the
      input interferograms (e.g. geotiff or flat-binary files).

Upon completion of ``conv2tif`` geotiff formatted copies of the input files will be placed
in the ``outdir`` directory defined in the configuration file.

.. note::

     - ``conv2tif`` will not perform the conversion if geotiffs for the provided
       input files already exist.


``prepifg``: Preparing input interferograms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``prepifg`` is the second step of `PyRate`, which applys multi-looking, cropping
and coherence masking operations to the geotiff interferograms.

::

    >> pyrate prepifg --help
    usage: pyrate prepifg [-h] -f CONFIG_FILE

    optional arguments:
      -h, --help            show this help message and exit
      -f CONFIG_FILE, --config_file CONFIG_FILE
                            Pass configuration file

The ``prepifg`` step is used as follows:

::

    >> pyrate prepifg -f /path/to/config_file


**Coherence masking**

If specified via the ``cohmask`` parameter, ``prepifg`` will perform coherence masking
on the interferograms before multi-looking and cropping is performed. This requires
corresponding coherence images for each interferogram. The purpose
of coherence masking is to remove poor quality phase observations and leave a set of
high-quality pixels for analysis. Pixels with coherence values below a certain threshold
(defined with ``cohthresh`` parameter) will be set to Not-a-Number (NaN). 

.. note::

    - The number of pixels with numeric phase values (i.e. pixels not equal to NaN)
      in each interferogram will be different after coherence masking.

Coherence masking is enabled by setting the ``cohmask`` parameter to ``1`` in
the configuration file. A threshold, ``cohthresh`` needs to be provided. 
For every pixel where the coherence is lower than ``cohthresh`` the phase will be
changed to a NaN.
The available coherence files need to be specified in a list file as described above
and defined in the ``cohfilelist`` parameter.


**Multi-looking**

The ``prepifg`` step will perform optional multi-looking (image sub-sampling) 
of the input interferograms in geotiff format. The purpose of multi-looking is twofold:

- Reduce the spatial resolution of the interferograms in order to improve the
computational efficiency of `PyRate` analysis.
- Reduce the general phase noise in the interferograms.

To multi-look, set ``ifglksx`` and ``ifglksy`` to the integer subsampling factor
in the x (easting) and y (northing) dimensions respectively. Separate parameters for x
and y gives flexibility for users in case they want to achieve different spatial
resolution in in each dimension.

.. note::

    - For example, a value of ``2`` will reduce the resolution by half.
      A value of ``1`` will keep the resolution the same as the input interferograms
      (i.e. no multi-looking).
    - It is recommended to try a large multi-look factor to start with (e.g. ``10``
      or greater), and subsequently reduce the multi-looking factor once the user
      has experience with processing a certain dataset.


**Cropping**

The ``prepifg`` step will perform optional spatial cropping of the input interferograms.
This is useful if you are focussing on a specific area of interest within the full
extent of the input interferograms. The advantage of cropping is that `PyRate`
analysis will be computationally more efficient.

To crop, set ``ifgcropopt`` to ``3`` and provide the geographic latitude and longitude
bounds in the ``ifgxfirst`` (west), ``ifgxlast`` (east), ``ifgyfirst`` (north), and
``ifgylast`` (south) parameters.


``process``: Main workflow, including stacking and time series analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``process`` is the core of the `PyRate` processing workflow. This step will perform
a series of corrections to the interferogram phase data before running the time series
and stacking analysis.

::

    >> pyrate process --help
    usage: pyrate process [-h] -f CONFIG_FILE

    optional arguments:
      -h, --help            show this help message and exit
      -f CONFIG_FILE, --config_file CONFIG_FILE
                            Pass configuration file

The ``process`` step is used as follows:

::

    >> pyrate process -f path/to/config_file


Optionally, an orbital error correction and a spatio-temporal filter
operation to estimate and remove atmospheric phase screen (APS) signals is
applied to the interferograms prior to time series and stacking
analysis. The corrected interferograms are updated on disk and the
corrections are not re-applied upon subsequent runs. This functionality
is controlled by the ``orbfit`` and ``apsest`` options in the
configuration file.

Non-optional pre-processing steps include: 
- Minimum Spanning Tree matrix calculation,
- Identification of a suitable reference phase area,
- Removal of reference phase from interferograms,
- Calculation of interferogram covariance,
- Assembly of the variance-covariance matrix.

Following the above processing steps the time series and stacking
calculations are run. Time series is optional, controlled by the 
``tscal`` parameter. Stacking is not optional.


``merge``: Reassemble the tiles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``merge`` is the last step of the `PyRate` workflow, which produces geotiff
files containing the final time series and stacking products.
``merge`` will also re-assemble tiles that were generated during the ``process``
step when run with MPI. Tiling is discussed in the :ref:`parallel_label` section
below.

::

    >> pyrate merge --help
    usage: pyrate merge [-h] [-f CONFIG_FILE]

    optional arguments:
      -h, --help            show this help message and exit
      -f CONFIG_FILE, --config_file CONFIG_FILE
                            Pass configuration file

The ``merge`` step is used as follows:

::

    >> pyrate merge -f path/to/config_file

After running the ``merge`` step, several geotiff products will appear in the
directory defined by the ``outdir`` parameter.


``workflow``: Run the full PyRate workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``workflow`` is a fifth option that will run all four steps (``conv2tif``, ``prepifg``,
``process``, and ``merge``) in order as a single job.

::

    >> pyrate workflow --help
    usage: pyrate workflow [-h] [-f CONFIG_FILE]

    optional arguments:
      -h, --help            show this help message and exit
      -f CONFIG_FILE, --config_file CONFIG_FILE
                            Pass configuration file

The ``workflow`` is used as follows:

::

    pyrate workflow -f path/to/config_file

.. note::

    - ``workflow`` will only be useful for users starting with flat-binary input files,
      since ``conv2tif`` is the first step to be run as part of this full workflow.


Input Files
-----------

`PyRate` currently supports input files generated by the `GAMMA` and `ROI\_PAC`
interferometry softwares. `PyRate` will determine the input format from the 
``processor:`` parameter in the configuration file (``0``: `ROI\_PAC`;
``1``: `GAMMA`).

.. note::

    - Support and development of `ROI\_PAC` has been discontinued.
    - `ROI\_PAC` support in `PyRate` will be deprecated in a future release.

`GAMMA`
^^^^^^^

Each `GAMMA` geocoded unwrapped interferogram requires three header files
to extract metadata required for data formatting: a geocoded DEM header
file (``demHeaderFile`` keyword in the configuration file) and the relevant
SLC image header files (``*slc.par``) found in the ``hdrfilelist``.
The header files for the first and second SLC images used in the formation
of a particular interferogram are found automatically by date-string pattern
matching based on date epochs given in the filenames.
A DEM with matching size and geometry to the interferograms can also be processed.
The DEM absolute path and filename are set with the ``demfile`` parameter.

`ROI\_PAC`
^^^^^^^^^^

Each `ROI\_PAC` geocoded unwrapped interferogram requires its own
header/resource file (``*.rsc``). These header files need to be
listed in the defined ``hdrfilelist``. In addition, the geocoded DEM
header file is required and its path and name are specified in the config file under
``demHeaderFile``. The geographic projection in the parameter ``DATUM:`` is extracted
from the DEM header file.
A DEM with matching size and geometry to the interferograms can also be processed.
The DEM absolute path and filename are set with the ``demfile`` parameter.

.. _parallel_label:

Parallel Processing
-------------------

By their very nature, interferograms are large files. This is particularly the case
for `Sentinel-1`_, which has an image swath of 250 km and a pixel resolution on the order
of tens of metres in IW-mode.
Consequently, InSAR processing can be computationally expensive and time consuming.
It therefore makes sense to parallelise processing operations wherever possible.

.. _`Sentinel-1`: https://sentinel.esa.int/web/sentinel/user-guides/sentinel-1-sar

`PyRate` can be run in parallel using standard multi-threading simply by turning
``parallel:  1`` in the configuration file to take advantage of multiple cores
on a single machine. The parameter ``processes`` sets the number of threads.

Alternatively, `PyRate` can be parallelised on a system with an installed MPI library
by using ``mpirun``:

::

    # Modify '-n' based on the number of processors available.
    mpirun -n 4 pyrate conv2tif -f path/to/config_file
    mpirun -n 4 pyrate prepifg -f path/to/config_file
    mpirun -n 4 pyrate process -f path/to/config_file
    mpirun -n 4 pyrate merge -f path/to/config_file

.. note::

    - In the case that `PyRate` is run using ``mpirun``, standard multi-threading is automatically
      disabled (i.e. equivalent to setting ``parallel:  0``).

- During ``conv2tif`` and ``prepifg``, parallelism is achieved by sending sub-lists of input
  files to each process.
- In the multi-threading case, parallelism in the ``process`` step is achieved by splitting
  individual pixels across the available processes in a pixel-by-pixel nature.
- In the MPI case, parallelism in the ``process`` step is achieved by splitting the images
  in to a grid of tiles, where the number of tiles equals the number of processes passed with
  the ``-n`` option. The number of tiles in x and y dimension are automatically calculated by
  `PyRate`, ensuring a roughly equivalent number in both dimensions. One of the functions of the 
  ``merge`` step is to reassemble these tiles in to the full image for each output product.
