Troubleshooting
===============

Some common/known issues with `PyRate` that users may encounter are described below.

If your issue is not covered below, please contact the `Geoscience Australia InSAR Team`_ for help.

.. _`Geoscience Australia InSAR Team`: mailto:insar@ga.gov.au

Corrections being skipped
-------------------------
**Problem**: When running the ``process`` step, many corrections are reported as ``Skipped: interferograms already corrected``, but I want to try different processing parameters!

::

    >> pyrate correct -f input_parameters.conf
    16:43:16 main 97 24732 INFO 0/0 Verbosity set to INFO.
    16:43:16 shared 1294 24732 INFO 0/0 Running process serially
    16:43:17 process 86 24732 INFO 0/0 Found 13 unique epochs in the 17 interferogram network
    16:43:17 process 134 24732 INFO 0/0 Searching for best reference pixel location
    16:43:18 process 155 24732 INFO 0/0 Selected reference pixel coordinate: (38, 58)
    16:43:18 process 170 24732 INFO 0/0 Calculating orbital correction
    16:43:18 shared 1255 24732 INFO 0/0 Skipped: interferograms already corrected
    16:43:18 process 198 24732 INFO 0/0 Calculating reference phase
    16:43:19 shared 1255 24732 INFO 0/0 Skipped: interferograms already corrected
    16:43:19 process 105 24732 INFO 0/0 Calculating minimum spanning tree matrix
    16:43:19 process 342 24732 INFO 0/0 Calculating the temporal variance-covariance matrix
    16:43:20 process 391 24732 INFO 0/0 Calculating time series using SVD method
    16:43:20 timeseries 152 24732 INFO 0/0 Calculating timeseries in serial
    16:43:21 process 323 24732 INFO 0/0 Calculating rate map from stacking
    16:43:21 process 326 24732 INFO 0/0 Stacking of tile 0
    16:43:21 stack 64 24732 INFO 0/0 Calculating stack rate in serial
    16:43:21 process 314 24732 INFO 0/0 PyRate workflow completed

**Reason**: `PyRate` updates the phase values in the input interferogram geotiff files as corrections are applied during the ``process`` step. Metadata is then added to the geotiff header to indicate the correction has been applied. This metadata is then checked upon subsequent runs to see if the correction should be applied.

**Solution**: Start again from ``prepifg`` step, creating new cropped/multi-looked interferograms that have not been corrected.

.. note::

    We plan to change this workflow behaviour in a future `PyRate` release, recognising that
    it would be convenient to be able to quickly test the impact of parameter changes.


ValueError: too many values to unpack (expected 2)
--------------------------------------------------
**Problem**: During ``prepifg`` step, the following error is encountered:

::

    Traceback (most recent call last):
      File "pyrate/main.py", line 131, in <module>
        main()
      File "pyrate/main.py", line 103, in main
        prepifg.main(params)
      File "~/PyRateVenv/lib/python3.7/site-packages/Py_Rate-0.4.0-py3.7.egg/pyrate/prepifg.py", line 57, in main
        do_prepifg(gtiff_paths, params)
      File "~/PyRateVenv/lib/python3.7/site-packages/Py_Rate-0.4.0-py3.7.egg/pyrate/prepifg.py", line 88, in do_prepifg
        _prepifg_multiprocessing(gtiff_path, xlooks, ylooks, exts, thresh, crop, params)
      File "~/PyRateVenv/lib/python3.7/site-packages/Py_Rate-0.4.0-py3.7.egg/pyrate/prepifg.py", line 112, in _prepifg_multiprocessing
        coherence_path=coherence_path, coherence_thresh=coherence_thresh)
      File "~/PyRateVenv/lib/python3.7/site-packages/Py_Rate-0.4.0-py3.7.egg/pyrate/core/prepifg_helper.py", line 187, in prepare_ifg
        op = output_tiff_filename(raster.data_path, out_path)
      File "~/PyRateVenv/lib/python3.7/site-packages/Py_Rate-0.4.0-py3.7.egg/pyrate/core/shared.py", line 1222, in output_tiff_filename
        fname, ext = os.path.basename(inpath).split('.')
    ValueError: too many values to unpack (expected 2)

**Reason**: This is caused by multiple double dots being used in the file names.

**Solution**: Ensure only a single dot is given in the file names, delimiting the file extension.
For example, rename ``20151219-20160112.unw.tif`` to ``20151219-20160112_unw.tif``.


Stack Rate map appears to be blank/empty
----------------------------------------
**Problem**: Output of Stack Rate algorithm contains NaN values causing “merge“ to fail:

::

    Traceback (most recent call last):
    File "~/PyRateVenv/bin/pyrate", line 11, in <module> load_entry_point('Py-Rate==0.4.0', 'console_scripts', 'pyrate')()
    File "~/PyRateVenv/lib/python3.7/site-packages/Py_Rate-0.4.0-py3.7.egg/main.py", line 154, in main merge_handler(args.config_file)
    File "~/PyRateVenv/lib/python3.7/site-packages/Py_Rate-0.4.0-py3.7.egg/main.py", line 87, in merge_handler merge.main(config.__dict__)
    File "~/PyRateVenv/lib/python3.7/site-packages/Py_Rate-0.4.0-py3.7.egg/merge.py", line 64, in main create_png_from_tif(output_folder_path)
    File "~/PyRateVenv/lib/python3.7/site-packages/Py_Rate-0.4.0-py3.7.egg/merge.py", line 130, in create_png_from_tif minimum, maximum, mean, stddev = srcband.GetStatistics(True, True)
    File "/apps/gdal/3.0.2/lib/python3.7/site-packages/osgeo/gdal.py", line 2610, in GetStatistics return _gdal.Band_GetStatistics(self, *args)
    RuntimeError: ~/out/stack_rate.tif, band 1: Failed to compute statistics, no valid pixels found in sampling.
    Traceback (most recent call last):
    File "~/PyRateVenv/bin/pyrate", line 11, in <module> load_entry_point('Py-Rate==0.4.0', 'console_scripts', 'pyrate')()
    File "~/PyRateVenv/lib/python3.7/site-packages/Py_Rate-0.4.0-py3.7.egg/main.py", line 154, in main merge_handler(args.config_file)
    File "~/PyRateVenv/lib/python3.7/site-packages/Py_Rate-0.4.0-py3.7.egg/main.py", line 87, in merge_handler merge.main(config.__dict__)
    File "~/PyRateVenv/lib/python3.7/site-packages/Py_Rate-0.4.0-py3.7.egg/merge.py", line 64, in main create_png_from_tif(output_folder_path)
    File "~/PyRateVenv/lib/python3.7/site-packages/Py_Rate-0.4.0-py3.7.egg/merge.py", line 130, in create_png_from_tif minimum, maximum, mean, stddev = srcband.GetStatistics(True, True)
    File "/apps/gdal/3.0.2/lib/python3.7/site-packages/osgeo/gdal.py", line 2610, in GetStatistics return _gdal.Band_GetStatistics(self, *args)
    RuntimeError: ~/out/stack_rate.tif, band 1: Failed to compute statistics, no valid pixels found in sampling.
    Traceback (most recent call last):
    File "~/PyRateVenv/bin/pyrate", line 11, in <module> load_entry_point('Py-Rate==0.4.0', 'console_scripts', 'pyrate')()
    File "~/PyRateVenv/lib/python3.7/site-packages/Py_Rate-0.4.0-py3.7.egg/main.py", line 154, in main merge_handler(args.config_file)
    File "~/PyRateVenv/lib/python3.7/site-packages/Py_Rate-0.4.0-py3.7.egg/main.py", line 87, in merge_handler merge.main(config.__dict__)
    File "~/PyRateVenv/lib/python3.7/site-packages/Py_Rate-0.4.0-py3.7.egg/merge.py", line 64, in main create_png_from_tif(output_folder_path)
    File "~/PyRateVenv/lib/python3.7/site-packages/Py_Rate-0.4.0-py3.7.egg/merge.py", line 130, in create_png_from_tif minimum, maximum, mean, stddev = srcband.GetStatistics(True, True)
    File "/apps/gdal/3.0.2/lib/python3.7/site-packages/osgeo/gdal.py", line 2610, in GetStatistics return _gdal.Band_GetStatistics(self, *args)
    RuntimeError: ~/out/stack_rate.tif, band 1: Failed to compute statistics, no valid pixels found in sampling.

**Reason**: The ``maxsig`` parameter is too low, resulting in stack rate values being replaced by NaNs. ``maxsig`` is a threshold for masking stack rate pixels according to the corresponding stack error estimate saved in ``out/tmpdir/stack_error_*.npy``.

**Solution**: Increase ``maxsig``, then re-run ``process`` and ``merge`` steps. Maximum permittable value for ``maxsig`` is 1000 mm.


Failure of APS spatial low pass filter
---------------------------------------
**Problem**: Atmospheric corrections during “process“ fails due to the interpolated grid used for correction being empty:

::

    +8s pyrate.aps:INFO Applying spatial low pass filter
    Traceback (most recent call last):
    File "~/PyRateVenv/bin/pyrate", line 11, in <module>
    load_entry_point('Py-Rate==0.3.0.post3', 'console_scripts', 'pyrate')()
    File "~/PyRateVenv/lib/python3.6/site-packages/Click-7.0-py3.6.egg/click/core.py", line 764, in __call__
    return self.main(*args, **kwargs)
    File "~/PyRateVenv/lib/python3.6/site-packages/Click-7.0-py3.6.egg/click/core.py", line 717, in main
    rv = self.invoke(ctx)
    File "~/PyRateVenv/lib/python3.6/site-packages/Click-7.0-py3.6.egg/click/core.py", line 1137, in invoke
    return _process_result(sub_ctx.command.invoke(sub_ctx))
    File "~/PyRateVenv/lib/python3.6/site-packages/Click-7.0-py3.6.egg/click/core.py", line 956, in invoke
    return ctx.invoke(self.callback, **ctx.params)
    File "~/PyRateVenv/lib/python3.6/site-packages/Click-7.0-py3.6.egg/click/core.py", line 555, in invoke
    return callback(*args, **kwargs)
    File "~/PyRateVenv/lib/python3.6/site-packages/Py_Rate-0.3.0.post3-py3.6.egg/pyrate/scripts/main.py", line 69, in linrate
    run_pyrate.process_ifgs(sorted(dest_paths), params, rows, cols)
    File "~/PyRateVenv/lib/python3.6/site-packages/Py_Rate-0.3.0.post3-py3.6.egg/pyrate/scripts/run_pyrate.py", line 391, in process_ifgs
    _wrap_spatio_temporal_filter(ifg_paths, params, tiles, preread_ifgs)
    File "~/PyRateVenv/lib/python3.6/site-packages/Py_Rate-0.3.0.post3-py3.6.egg/pyrate/aps.py", line 63, in _wrap_spatio_temporal_filter
    spatio_temporal_filter(tsincr, ifg, params, preread_ifgs)
    File "~/PyRateVenv/lib/python3.6/site-packages/Py_Rate-0.3.0.post3-py3.6.egg/pyrate/aps.py", line 86, in spatio_temporal_filter
    ts_aps = mpiops.run_once(spatial_low_pass_filter, ts_hp, ifg, params)
    File "~/PyRateVenv/lib/python3.6/site-packages/Py_Rate-0.3.0.post3-py3.6.egg/pyrate/mpiops.py", line 54, in run_once
    f_result = f(*args, **kwargs)
    File "~/PyRateVenv/lib/python3.6/site-packages/Py_Rate-0.3.0.post3-py3.6.egg/pyrate/aps.py", line 192, in spatial_low_pass_filter
    _interpolate_nans(ts_lp, params[cf.SLPF_NANFILL_METHOD])
    File "~/PyRateVenv/lib/python3.6/site-packages/Py_Rate-0.3.0.post3-py3.6.egg/pyrate/aps.py", line 208, in _interpolate_nans
    _interpolate_nans_2d(a, rows, cols, method)
    File "~/PyRateVenv/lib/python3.6/site-packages/Py_Rate-0.3.0.post3-py3.6.egg/pyrate/aps.py", line 224, in _interpolate_nans_2d
    method=method
    File "~/PyRateVenv/lib/python3.6/site-packages/scipy-1.3.0-py3.6-linux-x86_64.egg/scipy/interpolate/ndgriddata.py", line 226, in griddata
    rescale=rescale)
    File "interpnd.pyx", line 846, in scipy.interpolate.interpnd.CloughTocher2DInterpolator.__init__
    File "qhull.pyx", line 1836, in scipy.spatial.qhull.Delaunay.__init__
    File "qhull.pyx", line 276, in scipy.spatial.qhull._Qhull.__init__
    ValueError: No points given

**Solution**: Use more interferograms as input and/or reduce the threshold parameters ``ts_pthr``, ``pthr``, ``tlpfpthr`` in the configuration file.

In general, users are advised to input a network of small-baseline interferograms
that has at least 2 interferometric connections per SAR image epoch. Furthermore,
make sure that ``ts_pthr``, ``pthr`` and ``tlpfpthr`` are smaller than the number
of image epochs. To check that ``process`` worked correctly, users can check that
the ``tsincr_*.npy`` and ``tscuml*.npy`` arrays in the ``/<outdir>/tmpdir`` contain numeric values and not NaNs.


Out of memory errors
--------------------
**Problem**: `PyRate` is memory intensive. You may receive various out of memory errors if there is not enough memory to accommodate the images being processed.

::

    joblib.externals.loky.process_executor.TerminatedWorkerError: A worker process managed by the executor was unexpectedly terminated. This could be caused by a segmentation fault while calling the function or by an excessive memory usage causing the Operating System to kill the worker. The exit codes of the workers are {EXIT(1), EXIT(1), EXIT(1)}

**Solution**: Increase the amount of memory available. On HPC systems this can be done by increasing the value provided to the ``mem`` argument when submitting a PBS job, e.g.:

::

    mem=32Gb

Incorrect modules loaded on Gadi
----------------------------------
**Problem**: `PyRate` requires certain versions of Python, GDAL and Open MPI to be loaded on Gadi and other HPC systems. While sourcing the `PyRate/scripts/nci_load_modules.sh` script will load the correct modules, you may need to unload previously unloaded modules.

Example of errors caused by module conflicts::

    ERROR:150: Module 'python3/3.7.2' conflicts with the currently loaded module(s) 'python3/3.4.3-matplotlib'
    ERROR:150: Module 'gdal/2.2.2' conflicts with the currently loaded module(s) 'gdal/2.0.0'

**Solution**: Purge the loaded modules and source the ``nci_load_modules.sh`` script:

::

    module purge
    source ~/PyRate/scripts/nci_load_modules.sh
