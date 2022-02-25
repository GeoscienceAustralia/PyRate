.. :changelog:

Release History
===============

0.6.1 (2022-02-18)
------------------
Added
+++++
- List generator in ``utils/listcreator.sh`` for easier input generation from GAMMA output to PyRate.

Fixed
+++++
- Fix wrong sign in Y-intercept output file.
- Fix and simplify how user supplied reference pixel is validated and cropping issue.
- Add metadata for reference pixel latitude and longitude to output files.

0.6.0 (2021-10-18)
------------------
Added
+++++
- Geometry and baseline calculations, making use of GAMMA software MLI metadata and baseline files.
- DEM error estimation and correction functionality (untested prototype).
- Unwrapping error detection and masking functionality, making use of phase closure loops.
- Tests to check that the independent and network orbital methods are both able to recover the input parameters in synthetic examples.
- Tests for temporal and spatial gaussian filter options. Compare PyRate code against ``scipy.ndimage`` equivalents.
- More output visualisation scripts in the ``utils/`` directory, including an interferogram plotting script.
- New "Crop A" unit test dataset derived from Sentinel-1 data over Mexico City.
- Scaling factor parameter ``velerror_nsig`` for uncertainty/error output products.
- Calculate and output coherence statistics files during ``prepifg``.
- Add line-of-sight projection functionality for output products.
- Add signal polarity switching functionality for output products.
- Error handling for more than two MLI.par files in the header list matching an interferograms date-pair.

Fixed
+++++
- Fix bugs in the application of NaN masking and radian-to-millimetre conversion of interferogram data in ``correct`` step.
- Fix bugs in the implementation of the orbital error network and independent methods.
- Fix bugs in the implementation of spatial and temporal filters in the APS module.
- Fix bug in the application of reference phase subtraction that resulted in NaNs in the reference window.
- Fix file handling behaviour in ``prepifg`` so that a large number of files are not left open.

Changed
+++++++
- Enabled multi-looking and parallel processing for the independent orbital method.
- Simplify output filenames produced by PyRate.
- Move output files to dedicated named sub-directories. 
- Use of the ``offset`` parameter in orbital module; now renamed to ``intercept``.
- Input rasters to ``conv2tif`` can now have un-equal X and Y pixel resolutions.
- Change units of temporal filter from years to days.
- Moved Continuous Integration from Travis service to GitHub Actions.
- Made MPI an optional system dependency for PyRate (previously required).

Removed
+++++++
- Remove unused cython, glob2 and pillow dependencies.
- Remove support for Python 3.6.
- Remove support for filter types other than Gaussian in the APS module, and associated config options ``slpfmethod``, ``slpforder``, ``tlpfmethod``.
- ``config`` module deprecated; functionality moved to ``configuration``, ``shared`` and ``constants`` modules.
- Deprecated ``obs_dir``, ``slc_dir`` and ``coh_file_dir`` configuration parameters.

0.5.0 (2020-09-08)
------------------
Added
+++++
- New functionality "``linear_rate``" to calculate linear regression of
  cumulative displacement time series for every pixel as part of the ``timeseries`` step.
- Script for plotting ``timeseries`` and ``linear_rate`` output geotiff products using Matplotlib.
  To use, additional dependencies listed in ``requirements-plot.txt`` are required.
- Correction data (except ``maxvar`` and ``vcmt``) applied to the ifg data is saved to disk
  and re-used on subsequent repeat runs. Corrections are only re-calculated if config
  parameters change between runs.
- MPI parallelisation of APS spatio-temporal filter correction.
- Unit test coverage for refpixel lat/lon to x/y conversion and ``aps`` module.

Fixed
+++++
- Re-enable ``ifglksx`` and ``ifglksy`` to be different values, resulting in different
  resolutions in x and y dimensions in multi-looked interferograms.
- Re-enable ``orbfitlksx`` and ``orbfitlksy`` to be different values, resulting in different
  resolutions in x and y dimensions during network orbit correction.
- Screen messages from main process only during MPI runs.

Changed
+++++++
- ``process`` step has been renamed ``correct``. Stacking and timeseries have been removed from
  this step and are now invoked by separate ``timeseries`` and ``stack`` command line options.
- Processing of coherence files by ``conv2tif`` and ``prepifg`` is now triggered by the presence
  of ``cohfilelist`` in the config file. If the list is present, multilooked/cropped coherence
  files are saved to disk, regardless of whether ``cohmask`` is 0 or 1.
- Parallelisation capability is refactored - MPI and multiprocessing both now use a common
  tiling framework for ``stack``, ``timeseries`` and ``mst`` algorithms.
- Introduced a simplified and standardised file naming format for files produced by the
  ``prepifg`` step. Information previously saved in the filename (e.g. ``ifglksx``, ``ifglksy``,
  and ``ifgcropopt`` values) is now saved to the geotiff header instead.

Removed
+++++++
- Redundant ``tscal`` config parameter was deprecated - not needed now there is a ``timeseries``
  step invokable on the command line.
- Unused ``Pillow``, ``cython`` and ``glob2`` dependencies.
- Deprecated function ``pyrate.prepifg_helper.prepare_ifgs``, which is no longer needed.

0.4.3 (2020-08-04)
------------------
Added
+++++
- Ability to define the order of steps in the ``process`` workflow
  (default order unchanged).
  
Fixed
+++++
- Nil

Changed
+++++++
- ``prepifg`` output interferograms are saved as read-only files.
- ``process`` makes a writable copy of the ``prepifg`` output data
  at the beginning of each run.
- The selected reference pixel is saved to disk and re-used on subsequent
  ``process`` runs.  
- Saving of incremental time series (``tsincr``) products is optional,
  controlled by the ``savetsincr`` configuration parameter (default is on).

Removed
+++++++
- Removed obsolete InSAR terminology from code, docs and test data
  (changed to `first` and `second` images).
- Stopped using ``unittest`` unit test framework in favour of exclusively
  using ``pytest``.

0.4.2 (2020-06-26)
------------------
Added
+++++
- Save full-res coherence files to disk in ``conv2tif`` step if ``cohmask = 1``.
- Save multi-looked coherence files to disk in ``prepifg`` step if ``cohmask = 1``.
- Additional ``DATA_TYPE`` geotiff header metadata for above coherence files.
- ``conv2tif`` and ``prepifg`` output files have a tag applied to filename dependent
  on data type, i.e. ``_ifg.tif``, ``_coh.tif``, ``_dem.tif``.
- Metadata about used reference pixel is added to interferogram geotiff headers:
  lat/lon and x/y values; mean and standard deviation of reference window samples.
- Quicklook PNG and KML files are generated for the ``Stack Rate`` error map by default.

Fixed
+++++
- Ensure ``prepifg`` treats input data files as `read only`.
- Fix the way that the reference phase is subtracted from interferograms
  during ``process`` step.
- Manual entry of ``refx/y`` converted to type ``int``.

Changed
+++++++
- User supplies latitude and longitude values when specifying a reference pixel in
  the config file. Pixel x/y values are calculated and used internally.
- Move ``Stack Rate`` masking to a standalone function ``pyrate.core.stack.mask_rate``,
  applied during the ``merge`` step and add unit tests.
- Skip ``Stack Rate`` masking if threshold parameter ``maxsig = 0``.
- Provide log message indicating the percentage of pixels masked by 
  ``pyrate.core.stack.mask_rate``.
- Refactor ``pyrate.core.stack`` module; expose two functions in documentation:
  i) single pixel stacking algorithm, and
  ii) loop function for processing full ifg array.
- Refactor ``pyrate.merge`` script; remove duplicated code and create reusable
  generic functions.
- Colourmap used to render quicklook PNG images is calculated from min/max values of
  the geotiff band.
- Updated ``test`` and ``dev`` requirements.

Removed
+++++++
- Deprecate unused functions in ``pyrate.core.config`` and corresponding tests.
- Static colourmap ``utils/colourmap.txt`` that was previously used to render
  quicklook PNG images is removed. 

0.4.1 (2020-05-19)
------------------
Added
+++++
- Python 3.8 support.
- Algorithm to automatically calculate rows and columns for tiling.
  User no longer specifies these as part of the CLI, but can optionally
  specify ``rows`` and ``cols`` in the configuration file.
- Improvements to the test suite, including systems-wide tests.
- Improved logging.

Fixed
+++++
- Fixed bug in resampling/multi-looking when coherence masking is used.
  This bugfix will result in significantly fewer ``nan`` pixels in the outputs.
- Fixed a bug in how NaNs are handled during coherence masking and multi-looking.
  Output rasters will contain ``nan`` as the nodata value.

Changed
+++++++
- ``Linear Rate`` algorithm has been renamed ``Stack Rate``.
- User supplies full paths to input files in respective file lists.
- All files generated by `PyRate` saved to user-defined ``outdir`` directory.
- Renamed ``slcfilelist`` parameter to ``hdrfilelist``.
- Log files are generated in the ``outdir`` and every `PyRate` step produces independent log files.

Removed
+++++++
- Deprecate the use of ``obsdir``, ``slcfiledir`` and ``cohdir`` configuration variables.
- Deprecate ``parallel = 2`` option; splitting image via rows for parallelisation.

0.4.0 (2019-10-31)
------------------
Added
+++++
- Python 3.7 support.
- Optional ``conv2tif`` step.
- Building of docs integrated with Travis CI.
- Coherence masking, view coherence masking section in ``input_parameters.conf``
  for options.
- Input parameter validation.
- SLC and coherence file lists for file discovery.
- Create quick view png for rate map product.
- Add support for reading interferogram in Geotiff format.
- Add detailed validation and hints for configuration parameters
- Add system tests for all 3 types of input formats

Changed
+++++++
- ``linrate`` step has been renamed to ``process``.
- ``postprocess`` step has been renamed to ``merge``.
- ``converttogeotiff`` step has been renamed to ``conv2tif``.
- CLI structure: config files now need to be provided with ``-f`` flag.
- Reduced console output, default verbosity setting is now ``INFO``.
- Restructure of code layout, src modules now in ``PyRate/pyrate/core`` directory
  and scripts at ``PyRate/scripts``.
- Reference pixel values are expected to be in latitude and longitude values.

Removed
+++++++
- Unused luigi code.
- References to Matlab.
- Unused tests for legacy api.

0.3.0 (2019-07-26)
------------------
Added
+++++
- ``utils/apt_install.sh`` script that lists Ubuntu/apt package requirements.
- ``utils/load_modules.sh`` script that sets up NCI Raijin HPC environment.

Fixed
+++++
- Errors being caused by newer version of ``networkx``; v2.3 now supported.

Removed
+++++++
- Unused Python and OS packages.
- environment.yml - conda env will now be installed using ``requirements.txt``.
- HPC directory - hpc README.rst moved to docs.
- setup.cfg - no longer needed.
- Luigi functionality - hasn't been operational and is reported as vulnerable.
  Single machine parallelism is achieved with joblib. 

Changed
+++++++
- Requirements now managed by ``requirements.txt`` file, parsed by ``setup.py``.
- Requirements now split across base ``requirements.txt`` and separate files
  for dev (``requirements-dev.txt``) and testing (``requirements-test.txt``).
- Moved default config files to top level source directory.
- Pinned Python dependencies to specific versions.
- Travis build now installs GDAL from apt.
- Travis only builds on master, develop and \*-travis branches.
- Consolidated documentation into ``PyRate/docs``.
- Updated install instructions for Ubuntu and NCI.

0.2.0 (2017-05-22)
------------------
- Stable beta release.

0.1.0 (2017-01-31)
------------------
- First release on PyPI.
