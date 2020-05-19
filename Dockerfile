FROM ubuntu:16.04

ENV LANG C.UTF-8
ENV LANG=en_AU.UTF-8
ENV WORKON_HOME=$HOME/venvs
ENV GDALINST=$HOME/gdalinstall
ENV GDALBUILD=$HOME/gdalbuild
ENV PROJINST=$HOME/gdalinstall
ENV GDALVERSION="3.0.2"
ENV PROJVERSION="6.1.1"
ENV PROJOPT="--with-proj=$GDALINST/gdal-$GDALVERSION"
ENV PATH=$GDALINST/gdal-$GDALVERSION/bin:$PATH
ENV LD_LIBRARY_PATH=$GDALINST/gdal-$GDALVERSION/lib:$LD_LIBRARY_PATH
ENV PROJ_LIB=$GDALINST/gdal-$GDALVERSION/share/proj
ENV GDAL_DATA=$GDALINST/gdal-$GDALVERSION/share/gdal


RUN apt-get update \
    && apt-get install -y build-essential checkinstall libreadline-gplv2-dev \
    libncursesw5-dev libssl-dev libsqlite3-dev tk-dev libgdbm-dev libc6-dev \
    libbz2-dev openssl curl libffi-dev

RUN mkdir -p $HOME/opt

RUN cd $HOME/opt \
    && curl -O https://www.python.org/ftp/python/3.7.7/Python-3.7.7.tgz \
    && tar -xzf Python-3.7.7.tgz \
    && cd Python-3.7.7 \
    && ./configure --enable-shared --enable-optimizations --prefix=/usr/local LDFLAGS="-Wl,--rpath=/usr/local/lib" \
    && make altinstall

RUN apt-get install -y build-essential python3-pip \
    apt-utils git libgdal-dev libatlas-base-dev openmpi-bin \
    libopenmpi-dev gfortran wget libhdf5-serial-dev sqlite3 vim

# update pip
RUN python3.7 -m pip install pip --upgrade
RUN python3.7 -m pip install wheel
RUN pip3 install --upgrade setuptools


RUN mkdir -p $GDALBUILD $GDALINST $PROJBUILD $PROJINST

ENV GDALOPTS="  --with-geos \
            --with-expat \
            --without-libtool \
            --with-libz=internal \
            --with-libtiff=internal \
            --with-geotiff=internal \
            --without-gif \
            --without-pg \
            --without-grass \
            --without-libgrass \
            --without-cfitsio \
            --without-pcraster \
            --with-netcdf \
            --with-png=internal \
            --with-jpeg=internal \
            --without-gif \
            --without-ogdi \
            --without-fme \
            --without-hdf4 \
            --with-hdf5 \
            --without-jasper \
            --without-ecw \
            --without-kakadu \
            --without-mrsid \
            --without-jp2mrsid \
            --without-mysql \
            --without-ingres \
            --without-xerces \
            --without-odbc \
            --with-curl \
            --without-sqlite3 \
            --without-idb \
            --without-sde \
            --without-perl \
            --without-python"

RUN cd $PROJBUILD && wget -q https://download.osgeo.org/proj/proj-$PROJVERSION.tar.gz \
    && tar -xzf proj-$PROJVERSION.tar.gz \
    && cd proj-$PROJVERSION \
    && ./configure --prefix=$PROJINST/gdal-$GDALVERSION \
    && make -s -j 2 \
    && make install

RUN  cd $GDALBUILD \
    && wget -q http://download.osgeo.org/gdal/$GDALVERSION/gdal-$GDALVERSION.tar.gz \
    && tar -xzf gdal-$GDALVERSION.tar.gz

RUN cd $GDALBUILD/gdal-$GDALVERSION \
    && ./configure --prefix=$GDALINST/gdal-$GDALVERSION $GDALOPTS $PROJOPT

RUN cd $GDALBUILD/gdal-$GDALVERSION && make

RUN cd $GDALBUILD/gdal-$GDALVERSION  && make install

RUN pip install virtualenv virtualenvwrapper
ENV VIRTUALENVWRAPPER_PYTHON=/usr/local/bin/python3.7

ADD . / PyRate/
SHELL ["/bin/bash", "-c"]
RUN source /usr/local/bin/virtualenvwrapper.sh \
    && mkvirtualenv -p python3.7 pyrate \
    && cd PyRate \
    && sed -i '/^GDAL/d' requirements.txt \
    && workon pyrate \
    && pip install -r requirements.txt -r requirements-dev.txt -r requirements-test.txt \
    && pip install GDAL==$(gdal-config --version) \
    && python setup.py install
