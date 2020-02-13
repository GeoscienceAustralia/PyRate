#!/bin/sh
#   This script is part of the PyRate software package.
#
#   Copyright 2020 Geoscience Australia
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
# originally contributed by @rbuffat to Toblerity/Fiona
set -e

GDALOPTS="  --with-ogr \
            --with-geos \
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
            --without-netcdf \
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
            --without-bsb \
            --without-grib \
            --without-mysql \
            --without-ingres \
            --without-xerces \
            --without-odbc \
            --with-curl \
            --without-sqlite3 \
            --without-dwgdirect \
            --without-idb \
            --without-sde \
            --without-perl \
            --without-php \
            --without-ruby \
            --with-python \
            --with-static-proj4=/usr/lib"

# Create build dir if not exists
if [ ! -d "$GDALBUILD" ]; then
  mkdir $GDALBUILD;
fi

if [ ! -d "$GDALINST" ]; then
  mkdir $GDALINST;
fi

ls -l $GDALINST

if [ "$GDALVERSION" = "1.9.2" -a ! -d "$GDALINST/gdal-$GDALVERSION" ]; then
  cd $GDALBUILD
  wget http://download.osgeo.org/gdal/gdal-$GDALVERSION.tar.gz
  tar -xzf gdal-$GDALVERSION.tar.gz
  cd gdal-$GDALVERSION
  ./configure --prefix=$GDALINST/gdal-$GDALVERSION $GDALOPTS
  make -s -j 2
  make install
fi

if [ "$GDALVERSION" != "1.9.2" -a ! -d "$GDALINST/gdal-$GDALVERSION" ]; then
  cd $GDALBUILD
  wget http://download.osgeo.org/gdal/$GDALVERSION/gdal-$GDALVERSION.tar.gz
  tar -xzf gdal-$GDALVERSION.tar.gz
  cd gdal-$GDALVERSION
  ./configure --prefix=$GDALINST/gdal-$GDALVERSION $GDALOPTS
  make -s -j 2
  make install
fi

# change back to travis build dir
cd $TRAVIS_BUILD_DIR
