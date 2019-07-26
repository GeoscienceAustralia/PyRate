FROM ubuntu:18.04

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update -y && apt-get -y install git gdal-bin=2.2.3* libgdal-dev libatlas-base-dev openmpi-bin libopenmpi-dev libnetcdf13

ENV C_INCLUDE_PATH=/usr/include/gdal
ENV CPLUS_INCLUDE_PATH=/usr/include/gdal
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

RUN apt install -y python3 python3-dev python3-pip
RUN pip3 install --upgrade setuptools

ADD . / PyRate/
RUN cd PyRate && sed -i 's/GDAL//g' requirements.txt
RUN cd PyRate && pip3 install -r requirements.txt
RUN pip3 install GDAL==$(gdal-config --version) --global-option=build_ext --global-option="-I/usr/include/gdal"
RUN cd PyRate && pip3 install -r requirements-dev.txt
RUN cd PyRate && python3 setup.py install