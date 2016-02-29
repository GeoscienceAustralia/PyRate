import os, datetime, numpy
from osgeo import gdal
import numpy as np
import numpy.ma as ma
import netCDF4 as nc
from functools import partial
from glob import glob
try:
    import Image
except:
    from PIL import Image

#: The name of the data variable in netcdf files used herein.
DATA_VARIABLE_NAME = 'data'

#: The name of the time dimension and variable in netcdf files used herein.
TIME_VARIABLE_NAME = 'date'

#: The name of the *Y* dimension and variable in netcdf files used herein.
LAT_VARIABLE_NAME = 'lat'

#: The name of the *X* dimension and variable in netcdf files used herein.
LON_VARIABLE_NAME = 'lon'

#: The name of the file that the frequency of *NaN* values is written to.
NA_FREQENCY_IMAGE_NAME = 'na_frequency.png'

#: The name of the file that holds the time series stack of images in the
#: working folder.
IMAGE_STACK_FILE_NAME = 'test.ncdf'

#: The projection of the data. This is returned in HTTP responses and must
#: be suitable for use with openlayers.
PYRATE_DATA_PROJECTION = 'EPSG:4326'

class Tile(object):
    """
    A lightweight wrapper around a GDAL dataset providing comparison based on
    the date specified in the metatdata (retrieved with
    ``GetMetadata()['DATE']``)  and access to the date and data via the
    properties *date* and *data* respectively.

    :param path: The path to the GDAL dataset. Checked for existence upon
        construction.
    """

    def __init__(self, path):
        if not os.path.exists(path):
            raise Exception("path '{}' does not exist".format(path))
        self.path = path
        self._date = None
        self._tileData = None

    def __str__(self):
        return self.path

    def __cmp__(self, other):
        # should I assert identical product, location, ...?
        return cmp(self.date, other.date)

    @property
    def date(self):
        """
        The date of the tile. Note that for interferograms this is the first of the two dates in the
        metdatdata (extracted useing ``Dataset.GetMetadata()['DATE'])``).
        """

        if self._date is None:
            try:
                gdat = gdal.Open(self.path)
                dateStr = gdat.GetMetadata()['DATE']
                self._date = datetime.datetime.strptime(dateStr, '%Y-%m-%d').date().strftime('%Y%m%d')
            finally:
                gdat = None

        if self._date is None:
            raise Exception("could not get date from '{}'".format(self.path))

        return self._date

    @property
    def data(self):
        """
        The data from (band 1 of) the dataset.
        """

        if self._tileData is None:
            try:
                gdat = gdal.Open(self.path)
                self._tileData = gdat.GetRasterBand(1).ReadAsArray()
            finally:
                gdat = None

        return self._tileData





class TimeSeries(object):
    """
    Represents a time series. This has two attributes *times* (a list of
    :py:class:`datetime.date` objects) and *data* (a list of values of
    unspecified type - but probably :py:class:`float`\s).

    :param list times: The times for the time series.
    :param list data: The data for the time series.
    """

    def __init__(self, times, data):
        self.times = [datetime.datetime.strptime(str(t), '%Y%m%d') for t in times]
        self.data = [float(d) for d in data]

    def __str__(self):
        return '\n'.join(['(%s, %s)' % (time, data) for time, data in zip(self.times, self.data)])





def getTiles(rootDir, ext):
    """
    Return a list of :py:class:`Tile`\s, for all files ending with *ext* in
    *rootDir*.
    """

    tileNames = glob(os.path.join(rootDir, '*.{}'.format(ext)))
    return [Tile(p) for p in tileNames]





#: The default color scheme to use for images.
DEFAULT_COLOR_TABLE = [
    (  0, (0,     0, 255)),
    ( 30, (0,   255, 255)),
    ( 50, (0,   255,   0)),
    ( 70, (255, 255,   0)),
    (100, (255,   0,   0))]





def generateColorSchemeForNans(nLayers):
    """
    Creates a grey scale color scheme with *nLayers* shades of grey. The color
    scheme has the same structure as :py:data:`DEFAULT_COLOR_TABLE` except that
    the first element in each of the tuples is a propoprtion, not a percentage.
    """

    mx = 255.
    assert(nLayers <= mx)
    return [(i/nLayers, (i*mx/nLayers,)*3) for i in range(1, nLayers+1)]





def generateImage(imageData, discrete, colorTable=DEFAULT_COLOR_TABLE):
    """
    Create a png image from *imageData* using the color scheme specified by
    *colorTable*. If *discrete* is *True*, then the color scheme will be
    discrete, otherwise it will be continuous (interpolated). It the case that
    it is discrete, the color for a pixel is based on the upper bound specified
    in the color table.

    :param imageData: A numpy `masked array`_ of the image data.

    :param bool discrete: Should the color scheme be discrete (*True*) or
        continuous (*False*).

    :param colorTable: If *colorTable* is iterable, its entries are tuples (or
        at least tuple like) with two elements;

        1. the first is a float scalar which specifies the lower percentile the
        color applies to, and

        2. the second is an rgb triple the values of which should be between 0
        and 255 inclusive (though this is not checked). See
        :py:data:`DEFAULT_COLOR_TABLE` for an example.

        If colorTable is callable, then it should return a structure similar to
        that described above, but the first entry in each tuple is should be
        the actual value for the lower bound (not the percentile). See
        :py:data:`generateColorSchemeForNans` for example.

    .. todo:: At the moment we only use :py:data:`DEFAULT_COLOR_TABLE`. If we
        allow the user to provide one in the future, should the RGBs in
        it (*colorTable*) be scaled, or do we specify that to the user? In the
        latter case, we may still want to check the values.

    .. _masked array: http://docs.scipy.org/doc/numpy/reference/maskedarray.html
    """

    maxImageValue = 255

    def applyDiscreteColor(b1, b2):
        inds = np.where(np.logical_and(~imageData.mask, imageData > b1[0], imageData <= b2[0]))
        red[inds] = b2[1][0]
        green[inds] = b2[1][1]
        blue[inds] = b2[1][2]

    def applyContinuousColor(b1, b2):
        inds = np.where(np.logical_and(~imageData.mask, imageData > b1[0], imageData <= b2[0]))
        scaled = (imageData[inds] - b1[0]) / (b2[0] - b1[0])
        oneminusscaled = 1 - scaled
        red[inds]   = (scaled * b2[1][0] + oneminusscaled * b1[1][0])
        green[inds] = (scaled * b2[1][1] + oneminusscaled * b1[1][1])
        blue[inds]  = (scaled * b2[1][2] + oneminusscaled * b1[1][2])

    if hasattr(colorTable, '__call__'):
        colorBands = colorTable()
    else:
        percentiles = np.percentile(imageData[~imageData.mask], [c[0] for c in colorTable])
        colorBands = [(p, rgb[1]) for p, rgb in zip(percentiles, colorTable)]

    red   = np.zeros(imageData.shape, np.uint8)
    green = np.zeros(imageData.shape, np.uint8)
    blue  = np.zeros(imageData.shape, np.uint8)
    alpha = np.zeros(imageData.shape, np.uint8)
    alpha.fill(maxImageValue)
    alpha[imageData.mask] = 0

    colorizer = applyDiscreteColor if discrete else applyContinuousColor

    for b1, b2 in zip(colorBands[:-1], colorBands[1:]):
        colorizer(b1, b2)

    return (
        Image.merge("RGBA", map(Image.fromarray, (red, green, blue, alpha))),
        colorBands)




def makePube(rootDir, workingDir, ext):
    """
    Build the netcdf 'image stack' with name *outFile* for all files in
    *rootDir* with file extension *ext*. The dataset has as many time slices
    as there are files.

    .. note:: The projection of maps returned from the API is set to
       :py:data:`PYRATE_DATA_PROJECTION`. If the raster data from which the
       pube is built does not return a GeoTransform (see
       :py:meth:`gdal.Dataset.GetGeoTransform`) which is compatible with this
       projection, then funny things will happen when one clicks on the map,
       draws transects or otherwise. Note that the default projection returned
       by :py:meth:`gdal.Dataset.GetGeoTransform` is probably not compatible.

    .. todo:: At present, the projection of the map in the browser is set to be
       the same as that of the data (see :py:data:`pyrate.viewer.web.BROWSER_PROJECTION`),
       so there is no need to transform locations being sent through the API.
       If one had a map that was in a different projection, this would cease
       to work. Reprojection could be done either on the browser or in Python.
       If done in Python, the appropriate place is probably
       :py:func:`pyrate.viewer.web.pixelScaler`, which would need to be provided
       with the projection for the map in the browser.

       The main reason I have provided this functionality (which is quite trivial
       for points) is that if one is serving to a map with a different projection
       to that of pyrate data, the images returned by the API would also need
       reprojection. The place to do this is not so clear cut as the functions
       :py:func:`generateImage` and :py:func:`pyrate.viewer.web.getImageData`
       need to do the same thing, so ideally some refactoring should be done.
    """

    outFile = os.path.join(workingDir, IMAGE_STACK_FILE_NAME)

    tiles = getTiles(rootDir, ext)
    tiles.sort()

    nTiles = len(tiles)
    if not nTiles:
        raise Exception('no tiles found')

    # get the
    inputProt = gdal.Open(tiles[0].path)
    nx = inputProt.RasterXSize
    ny = inputProt.RasterYSize
    (xul, xsz, xrot, yul, yrot, ysz) = inputProt.GetGeoTransform()
    inputProt = None

    ndat = nc.Dataset(outFile, 'w', 'NETCDF4')

    tdim = ndat.createDimension(TIME_VARIABLE_NAME, None)
    ydim = ndat.createDimension(LAT_VARIABLE_NAME, ny)
    xdim = ndat.createDimension(LON_VARIABLE_NAME, nx)

    lats  = ndat.createVariable( LAT_VARIABLE_NAME, 'f8', ( LAT_VARIABLE_NAME,))
    lons  = ndat.createVariable( LON_VARIABLE_NAME, 'f8', ( LON_VARIABLE_NAME,))
    times = ndat.createVariable(TIME_VARIABLE_NAME, 'i8', (TIME_VARIABLE_NAME,))
    data  = ndat.createVariable(DATA_VARIABLE_NAME, 'f4', (TIME_VARIABLE_NAME, LAT_VARIABLE_NAME, LON_VARIABLE_NAME))

    lons[:] = xul + xsz*numpy.arange(nx, dtype='float64')
    lats[:] = yul + ysz*numpy.arange(ny, dtype='float64')

    ndat.minx = xul
    ndat.maxx = xul + nx * xsz
    ndat.miny = yul + ny * ysz
    ndat.maxy = yul

    shape = (ny, nx)
    at = 0
    grid = np.zeros((len(ydim), len(xdim)), dtype=np.int)
    for tile in tiles:
        print 'loading tile {} with date {}'.format(tile, tile.date)
        tileData = tile.data
        if tile.data.shape != shape:
            # should we also check spatial extent here?
            raise Exception('tile {} has shape {}, but should have {}'.format(
                tile.data.shape, shape))
        data[at,:,:] = tileData
        grid += (~np.isnan(tileData)).astype(np.int)
        times[at] = tile.date
        at += 1

    ndat.close()

    imageData = ma.masked_array(grid, mask=False)

    image, colorTable = generateImage(imageData, False,
        colorTable=partial(generateColorSchemeForNans, at))

    outFileName = os.path.join(workingDir, NA_FREQENCY_IMAGE_NAME)
    with open(outFileName, 'wb') as buf:
        image.save(buf, format="PNG")





def getDimensions(stackFile):
    """
    Get the geographic and image space dimensions from a netCDF file created by
    :py:func:`makePube`.

    :return: Two lists containing the geographic and image space dimensions
        of the data in *stackFile*. The geographic dimensions are of the form
        ``[minx, miny, maxx, maxy]`` and the image space dimensions are
        ``[0, 0, nx, ny]`` where ``nx`` and ``ny`` are the width and height
        of the image in pixels respectively.
    """

    ndat = None
    try:
        ndat = nc.Dataset(stackFile)

        pixelDimensions = [0, 0,
            len(ndat.dimensions[LON_VARIABLE_NAME]),
            len(ndat.dimensions[LAT_VARIABLE_NAME])]

        geoDimensions = [
            getattr(ndat, 'minx'),
            getattr(ndat, 'miny'),
            getattr(ndat, 'maxx'),
            getattr(ndat, 'maxy')]

        return geoDimensions, pixelDimensions

    finally:
        if ndat is not None:
            ndat.close()





def extractTS(ncdfFile, row, col):
    """
    Extract a time series for a pixel at *(row, col)* of the results.

    :param str ncdfFile: The name of the netcdf file to extract the time series
        from.

    :param int row: The row of the data (where the bottom row is zero) to
        extract data for.

    :param int col: The column of the data (where the bottom row is zero) to
        extract data for.

    :return: An instance of :py:class:`TimeSeries`.
    """

    if not os.path.exists(ncdfFile):
        raise Exception("file '{}' does not exist".format(ncdfFile))

    ndat = nc.Dataset(ncdfFile) or None

    if row < 0 or col < 0 or \
        len(ndat.dimensions[LAT_VARIABLE_NAME]) <= row or \
        len(ndat.dimensions[LON_VARIABLE_NAME]) <= col:
            raise Exception('point (row:{}, col:{}) is out of bonds'.format(row, col))

    try:
        row = len(ndat.dimensions[LAT_VARIABLE_NAME]) - int(row) - 1
        col = int(col)
        timeDat = ndat.variables[TIME_VARIABLE_NAME][:]
        dataDat = ndat.variables[DATA_VARIABLE_NAME][:, row, col]
        return TimeSeries(timeDat, dataDat)
    finally:
        if ndat:
            ndat.close()
