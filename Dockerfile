# Build: docker build -t pyrate .
FROM osgeo/gdal:ubuntu-small-3.4.3

SHELL ["/bin/bash", "-c"]

# update system and install pthon tooling + thirdparty deps
RUN apt-get update
RUN apt-get install -y python3-pip python3-venv python3-wheel python3-dev
RUN apt-get install -y openmpi-bin libopenmpi-dev sqlite3

# update setuptools
RUN pip3 install --upgrade setuptools

# Copy PyRate source into image
COPY . /usr/src/PyRate
WORKDIR /usr/src/PyRate

# This is a hack the CI scripst also use, unit tests see, to assume it's read-only?
RUN chmod 444 tests/test_data/small_test/tif/geo_070709-070813_unw.tif  

# Create virtual env & install PyRate into it...
RUN python3 -m venv /usr/local/PyRate

# Dance around numpy/GDAL installation ordering requirements first...
RUN source /usr/local/PyRate/bin/activate && pip3 install $(grep numpy requirements.txt)
RUN source /usr/local/PyRate/bin/activate && pip3 install GDAL==$(gdal-config --version)

# Then install our requirements verbatim
RUN source /usr/local/PyRate/bin/activate && pip3 install -r requirements.txt
RUN source /usr/local/PyRate/bin/activate && pip3 install -r requirements-test.txt
RUN source /usr/local/PyRate/bin/activate && pip3 install -r requirements-dev.txt

# For running like a program
RUN python3 setup.py install --prefix=/usr/local/PyRate
ENTRYPOINT [ "/bin/bash", "--init-file", "/usr/local/PyRate/bin/activate", "-i", "-c" ]

# For unit testing
ENV PYTHONPATH "/usr/src/PyRate:${PYTHONPATH}"
ENV OMPI_ALLOW_RUN_AS_ROOT 1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM 1
CMD [ "pytest", "-m", "'not slow and not mpi'" ]
