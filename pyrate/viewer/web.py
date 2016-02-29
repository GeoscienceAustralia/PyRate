# -*- coding: utf-8 -*-
"""
Script for viewing the inputs/outputs from pyrate. Uses
"""

import os
import sys
import json
import urllib
import gc
import datetime
import types
import math
import StringIO as io
import webbrowser as wb

import numpy as np
import numpy.ma as ma

from osgeo import gdal

try:
    import Image
except:
    from PIL import Image
from pyrate.viewer.rbf import Rbf
#from scipy.interpolate import Rbf
from functools import wraps
from glob import glob
from flask import Flask, request, Response
from jinja2 import Environment, PackageLoader
from pyrate.viewer.datamanager import (
    PYRATE_DATA_PROJECTION,
    NA_FREQENCY_IMAGE_NAME,
    IMAGE_STACK_FILE_NAME,
    makePube,
    extractTS,
    getDimensions,
    generateImage)



gc.enable()

DEBUG = False
MISSING_VALUE = -99999999
MAX_REGION_RADIUS = 20
REGION_RADIUS_STOP_FRAC = 0.01
DEFAULT_KERNEL = 'gaussian'

#: The dimensions of the images we are working with.
#:
#: This is set to something more meaningful if the script is run.
DATASET_GEO_DIMENSIONS = None

if len(sys.argv) < 2:
    print 'usage: web.py <working_dir> [port_no]'
    sys.exit(0)

#: The name of the host the server will be running on.
SERVER_HOST = "localhost"

#: The port the server will be running on.
SERVER_PORT = sys.argv[2] if len(sys.argv) > 2 else 5000

#: The URL for the server.
CLIENT_URL = 'http://{}:{}'.format(SERVER_HOST, SERVER_PORT)

#: The projection used for the browser
BROWSER_PROJECTION = PYRATE_DATA_PROJECTION

#: The directory we are 'running against'.
WORKING_DIR = sys.argv[1]

#: The NetCDF file that contains the stack.
TIME_SERIES_FILE = os.path.join(WORKING_DIR, IMAGE_STACK_FILE_NAME)



#: The Flask app. This must be defined before it is used below.
app = Flask(__name__, static_folder='./static', static_url_path='/static')

#: The Fask environment.
env = Environment(loader=PackageLoader('pyrate.viewer', 'static'))



def getImageData(fileName):
    """
    Get the raster data from the image *fileName*.

    :param fileName: The name of the file containing the image data. This is
        opened using GDAL.

    :return: A dictionary with items:

        - *data*: The image data as a numpy array.
        - *date*: The date extracted from the dataset using::

            Dataset.GetMetadata()['DATE']
    """

    gdat = None
    try:
        gdat = gdal.Open(fileName)
        return {
            'data':gdat.GetRasterBand(1).ReadAsArray(),
            'date':str(datetime.datetime.strptime(str(gdat.GetMetadata()['DATE']), '%Y-%m-%d').date())}

    finally:
        gdat = None



class _FauxRbf(Rbf):
    """
    A 'wrapper' for :py:class:`pyrate.viewer.rbf.Rbf` used in calculating
    weights given to pixels a given distance from a point.
    """

    @staticmethod
    def _getKernel(kernel):
        if isinstance(kernel, types.StringTypes):
            # this is using the private API of Rbf... which may change without
            # notice.
            return getattr(Rbf, "_h_" + kernel)

        return kernel

    def __call__(self, r):
        return self.kernel(r)

    def __init__(self, epsilon, kernel):
        self.epsilon = epsilon
        import new
        self.kernel = new.instancemethod(
            _FauxRbf._getKernel(kernel),
            self,
            _FauxRbf)





def extractSubsetForSmoothing(data, row, col, epsilon, kernel, maxNeighborhood):
    """
    Extract the row and columns extents describing the region of *data* to be
    used for smoothing for a value at the pixel at *(row, col)*.

    :param data: The numpy array from which data is to be extracted. It is only
        required here to determine the shape of itself (that is only *data.shape*
        is used).

    :param float row: The row at which an estimate is to be made. It may be a
        point number but is converted to int before use.

    :param float col: The column at which an estimate is to be made. It may be a
        point number but is converted to int before use.

    :param float epsilon: The epsilon to be used in the kernel smoother
        (see :py:func:`scipy.interpolate.Rbf`).

    :param float kernel: The kernel to be used in the kernel smoother
        (see :py:func:`scipy.interpolate.Rbf`).

    :param int maxNeighborhood: The maximum radius (minus 1) of the neighbourhood
        to use when smoothing.
    """

    kernel = _FauxRbf(epsilon, kernel)
    tot = kernel(0.)
    at = 1
    while True:
        curr = kernel(at)
        tot += curr
        if (curr / tot) < REGION_RADIUS_STOP_FRAC or at > maxNeighborhood:
            break
        at += 1

    col = int(col)
    row = int(row)
    at = int(max(at, 5))
    shape = data.shape
    return [
        max(col-at, 0),        # xmin
        min(col+at, shape[1]), # xmax
        max(row-at, 0),        # ymin
        min(row+at, shape[0])] # ymax

    return slice(max(row-at, 0), min(row+at, shape[0])), slice(max(col-at, 0), min(col+at, shape[1]))



def smoother(xs, ys, zin, epsilon, kernel, maxNeighborhood):
    """
    This is adapted from http://wiki.scipy.org/Cookbook/RadialBasisFunctions.
    """

    def doSmooth(y, x):
        shape = zin.shape
        y = shape[0] - y
        [xmin, xmax, ymin, ymax] = extractSubsetForSmoothing(zin, y, x, epsilon, kernel, maxNeighborhood)
        z = zin[ymin:ymax, xmin:xmax]
        try:
            notnans = ~np.isnan(z)
            yg, xg = np.mgrid[ymin:(ymax+1), xmin:(xmax+1)] + 0.5
            xnn = xg[notnans].flatten()
            ynn = yg[notnans].flatten()
            znn =  z[notnans].flatten()
            if len(znn) == 0:
                return MISSING_VALUE
            else:
                rbf = Rbf(xnn, ynn, znn, epsilon=epsilon, function=kernel)
                res = rbf([x], [y])
                rbf = None
                gc.collect()
                return res
        except:
            return MISSING_VALUE

    assert(len(xs) == len(ys))
    return [doSmooth(y, x) for y, x in zip(ys, xs)]



def pixelScaler(x, y):
    """
    Convert *x* and *y* to pixel (column, row) coordinates. This is done using
    :py:data:`DATASET_GEO_DIMENSIONS` and :py:data:`DATASET_PIXEL_DIMENSIONS`.
    """

    lonDif = DATASET_GEO_DIMENSIONS[2] - DATASET_GEO_DIMENSIONS[0]
    latDif = DATASET_GEO_DIMENSIONS[3] - DATASET_GEO_DIMENSIONS[1]
    x = math.floor(DATASET_PIXEL_DIMENSIONS[2] * (x - DATASET_GEO_DIMENSIONS[0]) / lonDif)
    y = math.floor(DATASET_PIXEL_DIMENSIONS[3] * (y - DATASET_GEO_DIMENSIONS[1]) / latDif)
    return x, y



def jsonify(func):
    """
    Decorator that converts the outputs from *func* to JSON. if
    :py:data:`pyrate.viewer.web.DEBUG` is *False*, exceptions thrown from *func*
    (*e*) are JSON encoded as ``{'exception':str(e)}``. Otherwise they are left
    unchecked (which makes dubugging easier).
    """

    if DEBUG:
        # then we don't want to catch exceptions.
        def wrapper(*args, **kwargs):
            res = func(*args, **kwargs)
            return Response(json.dumps(res, indent=4), mimetype="application/json")
    else:
        def wrapper(*args, **kwargs):
            try:
                res = func(*args, **kwargs)
            except Exception, e:
                res = {'exception':str(e)}
            return Response(json.dumps(res, indent=4), mimetype="application/json")

    return wraps(func)(wrapper)



@app.route("/")
def index():
    """
    Get the index page for the viewer.
    """

    template = env.get_template('index.html')
    if DATASET_GEO_DIMENSIONS is None:
        return 'Could not get dataset dimensions'
    return template.render(
        base_url = CLIENT_URL,
        extent = DATASET_GEO_DIMENSIONS,
        missingValue = MISSING_VALUE,
        initialFilename = NA_FREQENCY_IMAGE_NAME,
        projection = BROWSER_PROJECTION)



@app.route("/files")
@jsonify
def files():
    """
    Get a JSON encoded list of tif files that are in *WORKING_DIR*.
    """

    files = [os.path.basename(p) for p in glob(os.path.join(WORKING_DIR, '*.tif'))]
    files.append(NA_FREQENCY_IMAGE_NAME)
    return files



@app.route('/image/<fileName>')
@jsonify
def getRegionImage(fileName):
    """
    Return a JSON encoded PNG for creating an OpenLayers layer. *fileName* must
    be a GDAL compatible dataset in :py:data:`WORKING_DIR`.

    This expects the HTTP parameter *discrete* to be in the request. If the
    value of this parameter is *true* (case insensitive) then the color scheme
    for the image will be discrete, otherwise it will be continuous.

    :return: A dictionary with items:

        - *url*: The URL for the encoded image. At present this is a `data URL`_,
        - *imageExtent*: The extent of the image (in pixels),
        - *imageSize*: The size of the image in pixels.

    .. todo: The image Extent is used to determine where the image is in space
        in the the viewer (in our case `OpenLayers`_). This only works for us
        because the projection we use there is pixel based, but we want to move
        to using some geographic projection, this will need to be changed.

    .. _data URL: http://tools.ietf.org/html/rfc2397

    .. _OpenLayers: http://openlayers.org/
    """

    if fileName == NA_FREQENCY_IMAGE_NAME:
        png = Image.open(os.path.join(WORKING_DIR, NA_FREQENCY_IMAGE_NAME))
        buf = io.StringIO()
        png.save(buf, format="PNG")
        url = 'data:image/png;base64,{}'.format(urllib.quote(buf.getvalue().encode("base64")))
        buf.close()
        return {
            'url':url,
            'imageExtent': DATASET_GEO_DIMENSIONS,
            'imageSize': DATASET_PIXEL_DIMENSIONS[2:],
            'projection': PYRATE_DATA_PROJECTION}

    discrete = request.args['discrete'].lower() == 'true'
    filePath = os.path.join(WORKING_DIR, fileName)
    data = gdat = None
    try:
        gdat = gdal.Open(filePath)
        data = gdat.GetRasterBand(1).ReadAsArray()
        noDataValue = gdat.GetRasterBand(1).GetNoDataValue()
        if noDataValue != noDataValue:
            data = ma.masked_invalid(data)
        else:
            data = ma.masked_equal(data, noDataValue)
        nx = gdat.RasterXSize
        ny = gdat.RasterYSize
        (xul, xsz, xrot, yul, yrot, ysz) = gdat.GetGeoTransform()
    finally:
        gdat = None

    if data is None:
        raise Exception("data not extracted from '{}'".format(fileName))

    image, colorTable = generateImage(data, discrete)
    buf = io.StringIO()
    image.save(buf, format="PNG")

    content = {
        'url':'data:image/png;base64,{}'.format(urllib.quote(buf.getvalue().encode("base64"))), #'/static/result.png'
        'imageExtent':DATASET_GEO_DIMENSIONS,#[0, 0, nx, ny],
        'imageSize':[nx, ny],
        'projection': PYRATE_DATA_PROJECTION,
        'colorTable':colorTable}
    buf.close()

    return content



@app.route("/ts")
@jsonify
def getTimeSeries():
    """
    Extract a time series (pixel drill) through the tif files in
    :py:data:`WORKING_DIR`.

    This expects the HTTP parameters *x* and *y* to be in the request. These
    are converted to :py:class:`float` and passed to
    :py:func:`pyrate.viewer.datapube.extractTS`. Please see that function for
    how they are interpreted.

    :return: A dictionary containing items:

        - *times*: A list of times as strings (formatted by calling *str* on
            instances of :py:class:`datetime.date`).

        - *data*: A list of (floating point) data values.
    """

    y = float(request.args['y'])
    x = float(request.args['x'])
    x, y = pixelScaler(x, y)
    ts = extractTS(
        TIME_SERIES_FILE,
        int(y),
        int(x))
    data = [float(d) if d==d else MISSING_VALUE for d in ts.data]
    times = [str(t.date()) for t in ts.times]

    return {
        'times':times,
        'data':data}



@app.route("/transect")
@jsonify
def transect():
    """
    Get estimates of the the value at a set of points along a line for each
    image found in :py:data:`WORKING_DIR`.

    This expects the HTTP parameter *geom* to be in the request. This must be
    a JSON encoded dictionary with elements *start* and *end* (case sensitive)
    specifying the start and end of the line respectively. An example of the
    value of the parameter is::

        {"start":[17.24,52.83],"end":[20.82,50.73]}

    :return: A list of dictionaries with items:

        - *data*: The estimates at the points.
        - *date*: The date from the image (as returned by
            :py:func:`pyrate.viewer.web.getImageData`).

        An example of which is::

            [
                {
                    "date": "2006-08-28",
                    "data": [
                        -11.900258597604543,
                        -12.415132600332,
                        -12.178865800445514,
                        -12.697776081923053,
                        -12.039423493666977
                    ]
                },
                {
                    "date": "2006-10-02",
                    "data": [
                        -2.314936565531913,
                        -2.1408805229372234,
                        -2.147135196349211,
                        -2.2752048465306034,
                        -1.8302705710257476
                    ]
                }
            ]
    """

    def evalPoints(p1, p2, n):
        """
        Get an evenly spaced set of locations between the points *p1* and *p2*,
        which are assumed to be in the projection of the data and converted to
        pixels using :py:func:`pixelScaler`.
        """

        p1 = pixelScaler(*p1)
        p2 = pixelScaler(*p2)
        return (
            np.linspace(p1[0], p2[0], n),
            np.linspace(p1[1], p2[1], n))

    def getSmoothed(fileName, xs, ys, epsilon, kernel, maxNeighborhood):
        """
        use :py:func:`pyrate.viewer.web.smoother` to get estimates of the
        surface at the locations defined by *xs* and *ys*.

        :param str fileName: The name of the file containing the image data.
        :param iterable xs: The 'X' locations at which to estimate the surface.
        :param iterable xs: The 'Y' locations at which to estimate the surface.
        """

        data = getImageData(os.path.join(WORKING_DIR, fileName))
        return {
            'data':[float(d) for d in smoother(
                xs, ys, data['data'], epsilon, kernel, maxNeighborhood)],
            'date':data['date']}

    def exampleSmooth(epsilon, kernel, maxNeighborhood):
        r = np.array(range(0, maxNeighborhood+1), dtype=np.float)
        k = _FauxRbf(epsilon, kernel)
        return {
            'x':[v for v in r],
            'y':[k(v) for v in r]}

    kernel = DEFAULT_KERNEL # could get this from the user in future
    points = json.loads(request.args['geom'])
    epsilon = float(request.args['epsilon'])
    maxNeigborhood = int(request.args['maxNeighborhood'])
    xat, yat = evalPoints(points['start'], points['end'], 5)
    files = [os.path.basename(p) for p in glob(os.path.join(WORKING_DIR, '*.tif'))]

    return {
        'transects':[getSmoothed(fn, xat, yat, epsilon, kernel, maxNeigborhood) for fn in files],
        'weights':exampleSmooth(epsilon, kernel, maxNeigborhood)}



def shutdown_server():
    """
    Shut down the web server. This is based on
    `this snippet <http://flask.pocoo.org/snippets/67/>`_.
    """

    func = request.environ.get('werkzeug.server.shutdown')
    if func is None:
        raise RuntimeError('Not running with the Werkzeug Server')
    func()

@app.route("/shutdown", methods=["POST"])
def shutdown():
    shutdown_server()
    return "shutting down"



if __name__ == "__main__":
    # make the stack if it does not already exist. This must be done before the
    # definition of :py:data:`DATASET_GEO_DIMENSIONS`.
    if not os.path.exists(TIME_SERIES_FILE):
        makePube(WORKING_DIR, TIME_SERIES_FILE, 'tif')

    DATASET_GEO_DIMENSIONS, DATASET_PIXEL_DIMENSIONS = getDimensions(TIME_SERIES_FILE)

    print 'DATASET_GEO_DIMENSIONS', DATASET_GEO_DIMENSIONS
    print 'DATASET_PIXEL_DIMENSIONS', DATASET_PIXEL_DIMENSIONS

    # open the web browser. we have to do this before we start the server, because
    # the call the starts the server blocks... doing this seems to work though, but
    # is perhaps because the time out of the browser is long enough, though this is
    # just a guess.
    if os.environ.get('WERKZEUG_RUN_MAIN') == 'true' or not DEBUG:
        wb.open(CLIENT_URL)

    # and start the server. this blocks, so we need to do it last.
    app.run(host=SERVER_HOST, port=SERVER_PORT, debug=DEBUG)
