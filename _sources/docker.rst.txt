Docker
------

Docker can be used to run PyRate on Linux, Windows and MacOS. Please follow the instruction at https://hub.docker.com/ to
download and install Docker on your machine. Then run:

::

    git clone git@github.com:GeoscienceAustralia/PyRate.git
    cd PyRate
    docker build -t pyrate-image:1.0 .
    docker run -it --rm -v %cd%:/PyRate pyrate-image:1.0 /bin/bash
