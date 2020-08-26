Docker
^^^^^^

Docker can be used to run `PyRate` on Linux, Windows and MacOS systems.

.. note::

    - Docker is the recommended method to use `PyRate` under Windows.

The only system pre-requisites for using docker are Git and Docker Desktop
(e.g. for Windows https://docs.docker.com/docker-for-windows/ ).
The other system dependencies required by PyRate are installed in to the
docker container itself during the build operation.
Once Git and Docker Desktop are installed, and Docker Desktop is running,
execute the following at a command line prompt (e.g. PowerShell in Windows 10)::

    git clone git@github.com:GeoscienceAustralia/PyRate.git
    cd PyRate
    docker build -t pyrate-image .
    docker run -it pyrate-image

.. note::

    - The image name “pyrate-image” is not mandatory.
      You are free to name the docker image whatever you choose.

The ``docker build`` command builds the docker container described in the
``Dockerfile`` in the `PyRate` repo. It first starts by setting up an Ubuntu linux
system image, and then installs the required system and python dependencies.

Once the docker container is running (you will see a different-looking command
prompt in the PowerShell), execute the following commands::

    source /usr/local/bin/virtualenvwrapper.sh
    workon pyrate
    cd PyRate
    python3 setup.py install
    pyrate --help
    pyrate workflow -f input_parameters.conf

The ``pyrate`` executable program is built as part of the ``docker build`` step.
If the ``pyrate`` executable is reported as "not found", re-run the compilation::

    cd PyRate
    python3 setup.py install

