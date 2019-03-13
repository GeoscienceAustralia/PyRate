#   This Python module is part of the PyRate software package.
#
#   Copyright 2017 Geoscience Australia
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
"""
This Python module contains the base clase that defines
conversion routines for handling metadata, reading of
various formats and output to GeoTIFF (which can be extended to other formats).
"""

import gdal
import rasterio
from pyrate import ifgconstants as ifc


#def metadata(header):
#    """
#    Standardise the metadata to enable simpler output for the
#    write_geotif function.
#
#    :param header:
#    :type dict:
#        Dictionary containing interferogram metadata.
#
#    :return:
#        A dict containing specific flags for either ifg, dem or
#        incidence.
#    """
#    # similiar to what is written in the Conversion base class' _pass_header
#
#    return


def rio_dataset(out_fname, columns, rows, driver="GTiff", bands=1,
                dtype="float32", metadata=None, crs=None, geotransform=None,
                creation_opts=None):
    """
    Initialises a rasterio dataset object for writing image data.
    The reason being that the write_function shouldn't need to deal with
    reading data. But we shouldn't be required to read the whole file in
    order to write data either.
    So a generic write_geotif is fine, but works better when you have
    read all the data into memory.
    """
    kwargs = {
        "width": columns,
        "height": rows,
        "driver": driver,
        "count": bands,
        "crs": crs,
        "dtype": dtype,
        "transform": geotransform
    }

    # driver creation options
    if creation_opts is not None:
        for k, v in creation_opts.items():
            kwargs[k] = v

    # create output dataset
    outds = rasterio.open(out_fname, 'w', **kwargs)

    # global metadata
    if metadata is not None:
        outds.update_tags(**metadata)

    return outds


#def gdal_dataset(out_fname, columns, rows, driver="GTiff", bands=1,
#                 dtype='float32', metadata=None, crs=None,
#                 geotransform=None, creation_opts=None):
#    """
#    Initialises a py-GDAL dataset object for writing image data.
#    """
#    gdal_dtype = gdal.GDT_Float32 if dtype == 'float32' else gdal.GDT_Int16
#
#    # create output dataset
#    driver = gdal.GetDriverByName(driver)
#    outds = driver.Create(out_fname, columns, rows, bands, dtype,
#                          options=creation_opts)
#
#    # geospatial info
#    outds.SetGeoTransform(geotransform)
#    outds.SetProjection(crs)
#
#    # global metadata
#    if metadata is not None:
#        outds.SetMetadataItem(k, str(v)) # just following PyRate's example, but what if value is numeric???
#
#    return outds


def write_geotiff(data, out_fname, dtype='float32', metadata=None, crs=None,
                  geotransform=None, creation_opts={'compress': 'packbits'}):
    """
    A generic routine for writing a NumPy array to a geotiff.
    """
    # only support "2 <= dims <= 3"
    if data.ndim == 3:
        count, height, width = data.shape
    elif data.ndim == 2:
        height, width = data.shape
    else:
        msg = "Only support dimensions of '2 <= dims <= 3'."
        raise TypeError(msg)

    kwargs = {
        "width": width,
        "height": height,
        "driver": "GTiff",
        "count": count,
        "crs": crs,
        "dtype": dtype,
        "transform": geotransform
    }

    # driver creation options
    if creation_opts is not None:
        for k, v in creation_opts.items():
            kwargs[k] = v

    # create output dataset
    with rasterio.open(out_fname, 'w', **kwargs) as outds:

        # global metadata
        if metadata is not None:
            outds.update_tags(**metadata)

        # write data
        if data.ndim = 3:
            for band_id, blob in enumerate(data):
                outds.write(blob, band_id+1)
        else:
            # single 2D band only (1D or > 3D should have failed earlier)
            outds.write(data, 1)


class Conversion(object):

    """
    *****************************************************************
    NOTE!
    This presents a concept for a future direction when there is a
    greater need to support more formats.
    *****************************************************************

    The base conversion class.
    Ideally this is subclassed for each format type that requires
    conversion.
    eg GAMMA, ROIPC, SNAP-HDF5, SNAP-BEAM-DIMAP etc

    And the read methods are unique to each format type.
    That way the workflow code just calls {conversion_class}.write_geo_tiff()

    Still need to handle multibands though.
    """

    def __init__(self, filename, header, out_fname):
        # we need to specify the native block read size
        # eg 1 row by all columns for BSQ raw binary files (gamma raw output file)
        # so the _read routine can define an iterator to read native blocks
        # and write them out, before reading the next block (save memory)

        self._filename = filename
        self._header = header
        self._out_fname = out_fname
        self._ifg = ifc.PYRATE_WAVELENGTH_METRES in header
        self_incidence = 'FILE_TYPE' in header
        self._ifg_proc = header[ifc.PYRATE_INSAR_PROCESSOR]
        self._ncols = header[ifc.PYRATE_NCOLS]
        self._nrows = header[ifc.PYRATE_NROWS]

        self._pass_header()

    def _pass_header(self):
        """
        Extract relevant info from the header into dict containing
        metadata for direct dumping into gdal file creation.
        """
        hdr = {}
        if self.ifg:
            for k in [ifc.PYRATE_WAVELENGTH_METRES, ifc.PYRATE_TIME_SPAN,
                      ifc.PYRATE_INSAR_PROCESSOR,
                      ifc.MASTER_DATE, ifc.SLAVE_DATE,
                      ifc.DATA_UNITS, ifc.DATA_TYPE]:
                # is it always str or could they be numeric?
                hdr[k] = str(self._header[k])
            if self.ifg_proc == GAMMA:
                for k in [ifc.MASTER_TIME, ifc.SLAVE_TIME, ifc.PYRATE_INCIDENCE_DEGREES]:
                # is it always str or could they be numeric?
                    hdr[k] = str(self._header[k])
        elif self.incidence:
            hdr[ifc.DATA_TYPE] = ifc.INCIDENCE
        else: # must be dem
            hdr[ifc.DATA_TYPE] = ifc.DEM

        self._geotransform = [
            header[ifc.PYRATE_LONG],
            header[ifc.PYRATE_X_STEP],
            0,
            header[ifc.PYRATE_LAT],
            0,
            header[ifc.PYRATE_Y_STEP]
        ]

        srs = osr.SpatialReference()
        res = srs.SetWellKnownGeogCS(header[ifc.PYRATE_DATUM])
        self._crs = res.ExportToWkt()

        self._metadata = hdr

    def _read_block(self):
        # this method needs to be overridded with the specific read routine
        # i.e. raw GAMMA, raw ROIPC, tif GAMMA, SNAP HDF5 etc
        # this should return a read iterator
        # where each call will read a portion
        # eg for data, window in {conversion_class}._read_block()
        raise NotImplementedError

    def write_geotiff(self):
        # output a geotiff according to PyRate's format spec

        # the method will operate as
        # create output dataset (via rio_dataset or gdal_dataset)
        # read a blob of data
        # write that blob of data to the appropriate block in the tif
        # read next blob of data

        # eg
        # outds = rio_dataset(...)
        # for data, window in {conversion_class}._read_block():
        #     outds.write(data, 1, window=window)
        raise NotImplementedError

    @property
    def filename(self):
        return self._filename

    @property
    def out_fname(self):
        return self._out_fname

    @property
    def ncols(self):
        return self._ncols

    @property
    def nrows(self):
        return self._nrows

    @property
    def ifg(self):
        return self._ifg

    @property
    def incidence(self):
        return self._incidence

    @property
    def ifg_proc(self):
        return self._ifg_proc

    @property
    def metadata(self):
        return self._metadata

    @property
    def crs(self):
        return self._crs

    @property
    def geotransform(self):
        return self._geotransform


class ConvertGamma(Conversion):

    """GAMMA file format converter."""

    def __init__(self, filename, header, out_fname):
        super(ConvertGamma, self).__init__(filename, header, out_fname)

    def write_geotiff(self):
        """
        Overwrite the base classes write_geotiff method.
        Also ignoring the need for a read iterator,
        which could be scrapped altogether.
        """
        # these two variables could be properties (unique to raw binary files)
        bytes_per_col, fmtstr = _data_format(self.ifg_proc, self.ifg,
                                             self.ncols)
        # this could be a property (unique to raw binary files)
        row_bytes = self.ncols * bytes_per_col

        # initialise the gdal dataset
        outds = gdal_dataset(self.out_fname, self.ncols, self.nrows,
                             metadata=self.metadata, crs=self.crs,
                             geotransform=self.geotransform,
                             creation_opts=['compress=packbits'])

        band = outds.GetRasterBand(1)

        with open(self.filename, 'rb') as src:
            for y in range(selfnrows):
                data = struct.unpack(fmtstr, src.read(row_bytes))
                band.WriteArray(np.array(data).reshape(1, self.ncols), yoff=y)

        # can't remember if there is a Close() method
        band = None
        outds = None
