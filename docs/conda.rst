Instructions for setting up an Anaconda virtual environment on a typical
Red Hat linux server.

Anaconda
--------

Install
^^^^^^^

For 64-bit linux download this version of the conda installer:
https://repo.continuum.io/archive/Anaconda2-4.1.1-Linux-x86\_64.sh

Run the installer using the linux bash interpreter:

::

    bash Anaconda2-4.1.1-Linux-x86_64.sh

Accept all default options. This will install anaconda in the
``~/anaconda2`` directory.

Set up proxys for pip (only applicable if working behind proxy):
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

::

    export http_proxy=<proxy_server>:<port>
    export https_proxy=<proxy_server>:<port>

Set up ``~/.condarc`` with proxy pass information (only applicable if working behind proxy):
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For conda installer to work you need to create a ``~/.condarc`` file
with the following proxy-pass information:

::

    proxy_servers:
        http: <proxy_server>:<port>
        https: <proxy_server>:<port>
        

Clone the PyRate repo:
""""""""""""""""""""""

The rest of the instructions assume that you cloned the ``PyRate`` Git
repository in your home directory:

::

    cd ~
    git clone git@github.com:GeoscienceAustralia/PyRate.git

This will prompt for your rsa key passphrase. If you have access to
``PyRate`` repo and once you have entered the passphrase, ``PyRate``
will clone in your current (home) directory.

Set up environment variables
""""""""""""""""""""""""""""

You need to add the directory containing the ``pyrate/`` directory to
the ``PYTHONPATH`` environment variable, and the environment variable
``PYRATEPATH`` needs to point to the root directory:

::

    export PYRATEPATH="~/PyRate"
    export PYTHONPATH=$PYRATEPATH/:$PYTHONPATH

Activate ``anaconda`` and install ``conda-env``
"""""""""""""""""""""""""""""""""""""""""""""""

::

    source ~/anaconda2/bin/activate ~/anaconda2/
    conda install -c conda conda-env        

The ``conda-env`` package enables YAML based installation.

Note that using the ``source activate`` command only works cleanly when
using the bash shell

Create the ``pyrate`` environment
"""""""""""""""""""""""""""""""""

The required packages for PyRate are listed in the environment.yml (YAML
format) file:

::

    conda env create -f $PYRATEPATH/environment.yml

Activate the ``pyrate`` environment
"""""""""""""""""""""""""""""""""""

::

    source ~/anaconda2/bin/activate pyrate

\

Deactivate
""""""""""

Once you are done using ``PyRate`` you could deactivate the conda
environment using:

::

    source deactivate

Back to main Anaconda environment if you need to:
"""""""""""""""""""""""""""""""""""""""""""""""""

::

    source activate /home/user/anaconda2

