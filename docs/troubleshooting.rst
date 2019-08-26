Troubleshooting
===============

Out of memory errors
--------------------
PyRate is memory intensive. You may recieve various out of memory errors if 
there is not enough memory to accomdate the images being processed.

Error::

    joblib.externals.loky.process_executor.TerminatedWorkerError: A worker process managed by the executor was unexpectedly terminated. This could be caused by a segmentation fault while calling the function or by an excessive memory usage causing the Operating System to kill the worker. The exit codes of the workers are {EXIT(1), EXIT(1), EXIT(1)}

**Solution**: increase the amount of memory available. On HPC systems this can
be done by increasing the value provided to the ``mem`` argument when 
launching a PBS job.

::

    mem=32gb

Incorrect modules loaded on Raijin
==================================
PyRate requires certain versions of Python, GDAL and OpenMPI to be loaded
on Raijin and other HPC systems. While sourcing the ``PyRate/utils/load_modules.sh``
script will load the correct modules, you may need to unload previously unloaded modules.

Example of errors caused by module conflicts::

    ERROR:150: Module 'python3/3.7.2' conflicts with the currently loaded module(s) 'python3/3.4.3-matplotlib'
    ERROR:150: Module 'gdal/2.2.2' conflicts with the currently loaded module(s) 'gdal/2.0.0'

**Solution**: Unload the conflicting modules and re-source the ``load_modules.sh`` script.

::

    module unload python3/3.4.3-matplotlib
    module unload gdal/2.0.0
    source PyRate/utils/load_modules.sh
