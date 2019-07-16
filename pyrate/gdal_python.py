#   This Python module is part of the PyRate software package.
#
#   Copyright 2017 Geoscience Australia
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
"""
This Python module contains bindings for the GDAL library
"""
# pylint: disable=too-many-arguments,R0914
import gc
from osgeo import gdal, gdalnumeric, gdalconst
from PIL import Image, ImageDraw
import numpy as np
from pyrate import ifgconstants as ifc

gdal.SetCacheMax(2**15)
GDAL_WARP_MEMORY_LIMIT = 2**10
LOW_FLOAT32 = np.finfo(np.float32).min*1e-10


def world_to_pixel(geo_transform, x, y):
    """
    Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate
    the pixel location of a geospatial coordinate;
    see: http://pcjericks.github.io/py-gdalogr-cookbook/raster_layers.html

    :param list geo_transform: Affine transformation coefficients
    :param float x: longitude coordinate
    :param float y: latitude coordinate

    :return: col: pixel column number
    :rtype: int
    :return: line: pixel line number
    :rtype: int
    """
    ul_x = geo_transform[0]
    ul_y = geo_transform[3]
    xres = geo_transform[1]
    yres = geo_transform[5]
    col = int(np.round((x - ul_x) / xres))
    line = int(np.round((ul_y - y) / abs(yres))) # yres has negative size

    return col, line


def crop(input_file, extents, geo_trans=None, nodata=np.nan):
    """
    Adapted from http://karthur.org/2015/clipping-rasters-in-python.html

    Clips a raster (given as either a gdal.Dataset or as a numpy.array
    instance) to a polygon layer provided by a Shapefile (or other vector
    layer). If a numpy.array is given, a "GeoTransform" must be provided
    (via dataset.GetGeoTransform() in GDAL). Returns an array. Clip features
    must be a dissolved, single-part geometry (not multi-part). Modified from:

    http://pcjericks.github.io/py-gdalogr-cookbook/raster_layers.html
    #clip-a-geotiff-with-shapefile

    Arguments:
        rast            A gdal.Dataset or a NumPy array
        features_path   The path to the clipping features
        geo_trans              An optional GDAL GeoTransform to use instead
        nodata          The NoData value; defaults to -9999.

    :param str input_file: input image file path
    :param list extents: extents for the cropped area
    :param list geo_trans: An optional GDAL GeoTransform to use instead
    :param int nodata: The NoData value; defaults to -9999

    :return: clip: cropped part of the image
    :rtype: ndarray
    :return: gt2: geotransform parameters for the cropped image
    :rtype: list
    """

    def image_to_array(i):
        """
        Converts a Python Imaging Library (PIL) array to a gdalnumeric image.
        """
        arr = gdalnumeric.fromstring(i.tobytes(), 'b')
        arr.shape = i.im.size[1], i.im.size[0]
        return arr

    raster = gdal.Open(input_file)
    # Can accept either a gdal.Dataset or numpy.array instance
    if not isinstance(raster, np.ndarray):
        if not geo_trans:
            geo_trans = raster.GetGeoTransform()
        raster = raster.ReadAsArray()
    else:
        if not geo_trans:
            raise ValueError('geo transform must be supplied')

    # Convert the layer extent to image pixel coordinates
    min_x, min_y, max_x, max_y = extents
    ul_x, ul_y = world_to_pixel(geo_trans, min_x, max_y)
    lr_x, lr_y = world_to_pixel(geo_trans, max_x, min_y)

    # Calculate the pixel size of the new image
    px_width = int(lr_x - ul_x)
    px_height = int(lr_y - ul_y)

    # If the clipping features extend out-of-bounds and ABOVE the raster...
    if geo_trans[3] < max_y:
        # In such a case... ul_y ends up being negative--can't have that!
        # iY = ul_y
        ul_y = 0

    # Multi-band image?
    try:
        clip = raster[:, ul_y:lr_y, ul_x:lr_x]

    except IndexError:
        clip = raster[ul_y:lr_y, ul_x:lr_x]

    # Create a new geomatrix for the image
    gt2 = list(geo_trans)
    gt2[0] = min_x
    gt2[3] = max_y

    # Map points to pixels for drawing the boundary on a blank 8-bit,
    #   black and white, mask image.
    points = [(min_x, min_y), (max_x, min_y), (max_x, max_y), (min_y, max_y)]
    pixels = []


    for point in points:
        pixels.append(world_to_pixel(gt2, point[0], point[1]))

    raster_poly = Image.new('L', size=(px_width, px_height), color=1)
    rasterize = ImageDraw.Draw(raster_poly)
    rasterize.polygon(pixels, 0)  # Fill with zeroes


    mask = image_to_array(raster_poly)

    # Clip the image using the mask
    try:
        clip = gdalnumeric.choose(mask, (clip, nodata))

    # If the clipping features extend out-of-bounds and BELOW the raster...
    except ValueError:
        # We have to cut the clipping features to the raster!
        rshp = list(mask.shape)
        if mask.shape[-2] != clip.shape[-2]:
            rshp[0] = clip.shape[-2]

        if mask.shape[-1] != clip.shape[-1]:
            rshp[1] = clip.shape[-1]

        mask.resize(*rshp, refcheck=False)

        clip = gdalnumeric.choose(mask, (clip, nodata))

    # AttributeError: 'numpy.ndarray' object has no attribute 'close'
    # raster.close()
    raster = None

    return clip, gt2


def resample_nearest_neighbour(input_tif, extents, new_res, output_file):
    """
    Nearest neighbor resampling and cropping of an image.

    :param str input_tif: input geotiff file path
    :param list extents: new extents for cropping
    :param float new_res: new resolution for resampling
    :param str output_file: output geotiff file path

    :return: dst: resampled image
    :rtype: ndarray
    """
    dst, resampled_proj, src, _ = _crop_resample_setup(extents, input_tif,
                                                       new_res, output_file)
    # Do the work
    gdal.ReprojectImage(src, dst, '', resampled_proj,
                        gdalconst.GRA_NearestNeighbour)
    return dst.ReadAsArray()


def _crop_resample_setup(extents, input_tif, new_res, output_file,
                         dst_driver_type='GTiff', out_bands=2):
    """
    Convenience function for crop/resample setup
    """
    # Source
    src_ds = gdal.Open(input_tif, gdalconst.GA_ReadOnly)
    src_proj = src_ds.GetProjection()

    # source metadata to be copied into the output
    meta_data = src_ds.GetMetadata()

    # get the image extents
    min_x, min_y, max_x, max_y = extents
    geo_transform = src_ds.GetGeoTransform()  # tuple of 6 numbers

    # Create a new geotransform for the image
    gt2 = list(geo_transform)
    gt2[0] = min_x
    gt2[3] = max_y
    # We want a section of source that matches this:
    resampled_proj = src_proj
    if new_res[0]:  # if new_res is not None, it can't be zero either
        resampled_geotrans = gt2[:1] + [new_res[0]] + gt2[2:-1] + [new_res[1]]
    else:
        resampled_geotrans = gt2

    px_height, px_width = _gdalwarp_width_and_height(max_x, max_y, min_x,
                                                     min_y, resampled_geotrans)

    # Output / destination
    dst = gdal.GetDriverByName(dst_driver_type).Create(
        output_file, px_width, px_height, out_bands, gdalconst.GDT_Float32)
    dst.SetGeoTransform(resampled_geotrans)
    dst.SetProjection(resampled_proj)

    for k, v in meta_data.items():
        dst.SetMetadataItem(k, v)

    return dst, resampled_proj, src_ds, src_proj


def _gdalwarp_width_and_height(max_x, max_y, min_x, min_y, geo_trans):
    """
    Modify pixel height and width
    """
    # modified image extents
    ul_x, ul_y = world_to_pixel(geo_trans, min_x, max_y)
    lr_x, lr_y = world_to_pixel(geo_trans, max_x, min_y)
    # Calculate the pixel size of the new image
    px_width = int(lr_x - ul_x)
    px_height = int(lr_y - ul_y)
    return px_height, px_width  # this is the same as `gdalwarp`


def crop_resample_average(
        input_tif, extents, new_res, output_file, thresh,
        out_driver_type='GTiff',
        match_pirate=False):
    """
    Crop, resample, and average a geotiff image.

    :param str input_tif: Path to input geotiff to resample/crop
    :param tuple extents: Cropping extents (xfirst, yfirst, xlast, ylast)
    :param list new_res: [xres, yres] resolution of output image
    :param str output_file: Path to output resampled/cropped geotiff
    :param float thresh: NaN fraction threshold
    :param str out_driver_type: The output driver; `MEM` or `GTiff` (optional)
    :param bool match_pirate: Match Matlab Pirate output (optional)

    :return: resampled_average: output cropped and resampled image
    :rtype: ndarray
    :return: out_ds: destination gdal dataset object
    :rtype: gdal.Dataset
    """
    dst_ds, _, _, _ = _crop_resample_setup(
        extents, input_tif, new_res, output_file,
        out_bands=2, dst_driver_type='MEM')

    # make a temporary copy of the dst_ds for Pirate style prepifg
    tmp_ds = gdal.GetDriverByName('MEM').CreateCopy('', dst_ds) \
        if (match_pirate and new_res[0]) else None

    resampled_average, src_ds_mem = \
        gdal_average(dst_ds, input_tif, thresh)
    src_dtype = src_ds_mem.GetRasterBand(1).DataType
    src_gt = src_ds_mem.GetGeoTransform()

    # write out to output geotif file
    driver = gdal.GetDriverByName(out_driver_type)

    # required to match Matlab Pirate output
    if tmp_ds:
        _matlab_alignment(input_tif, new_res, resampled_average, src_ds_mem,
                          src_gt, tmp_ds)

    # write final pyrate GTiff
    out_ds = driver.Create(output_file, dst_ds.RasterXSize, dst_ds.RasterYSize,
                           1, src_dtype, options=['compress=packbits'])

    out_ds.GetRasterBand(1).SetNoDataValue(np.nan)
    out_ds.GetRasterBand(1).WriteArray(resampled_average)
    out_ds.SetGeoTransform(dst_ds.GetGeoTransform())
    out_ds.SetProjection(dst_ds.GetProjection())
    # copy metadata
    for k, v in dst_ds.GetMetadata().items():
        if k == ifc.DATA_TYPE:
            # update data type metadata
            if v == ifc.ORIG:
                out_ds.SetMetadataItem(ifc.DATA_TYPE, ifc.MULTILOOKED)
            elif v == ifc.DEM:
                out_ds.SetMetadataItem(ifc.DATA_TYPE, ifc.MLOOKED_DEM)
            elif v == ifc.INCIDENCE:
                out_ds.SetMetadataItem(ifc.DATA_TYPE, ifc.MLOOKED_INC)
            elif v == ifc.MULTILOOKED:
                pass
            else:
                raise TypeError('Data Type metadata not recognised')
        else:
            out_ds.SetMetadataItem(k, v)

    return resampled_average, out_ds


def _matlab_alignment(input_tif, new_res, resampled_average, src_ds_mem,
                      src_gt, tmp_ds):
    """
    Correction step to match python multi-look/crop output to match that of
    Matlab Pirate code. Modifies the resampled_average array in place.
    """
    src_ds = gdal.Open(input_tif)
    data = src_ds.GetRasterBand(1).ReadAsArray()
    xlooks = ylooks = int(new_res[0] / src_gt[1])
    xres, yres = _get_matlab_resampled_data_size(xlooks, ylooks, data)
    nrows, ncols = resampled_average.shape
    # Matlab Pirate does nearest neighbor resampling for the last
    # [yres:nrows, xres:ncols] cells without nan_conversion
    # turn off nan-conversion
    src_ds_mem.GetRasterBand(1).SetNoDataValue(LOW_FLOAT32)
    # nearest neighbor resapling
    gdal.ReprojectImage(src_ds_mem, tmp_ds, '', '',
                        gdal.GRA_NearestNeighbour)
    # only take the [yres:nrows, xres:ncols] slice
    if nrows > yres or ncols > xres:
        resampled_nearest_neighbor = tmp_ds.GetRasterBand(1).ReadAsArray()
        resampled_average[yres - nrows:, xres - ncols:] = \
            resampled_nearest_neighbor[yres - nrows:, xres - ncols:]


def gdal_average(dst_ds, input_tif, thresh):
    """
    Perform subsampling of an image by averaging values

    :param gdal.Dataset dst_ds: Destination gdal dataset object
    :param str input_tif: Input geotif
    :param float thresh: NaN fraction threshold

    :return resampled_average: resampled image data
    :rtype: ndarray
    :return src_ds_mem: Modified in memory src_ds with nan_fraction in Band2. The nan_fraction
        is computed efficiently here in gdal in the same step as the that of
        the resampled average (band 1). This results is huge memory and
        computational efficiency
    :rtype: gdal.Dataset
    """
    src_ds, src_ds_mem = _setup_source(input_tif)
    src_ds_mem.GetRasterBand(2).SetNoDataValue(-100000)
    src_gt = src_ds.GetGeoTransform()
    src_ds_mem.SetGeoTransform(src_gt)
    gdal.ReprojectImage(src_ds_mem, dst_ds, '', '', gdal.GRA_Average)
    # dst_ds band2 average is our nan_fraction matrix
    nan_frac = dst_ds.GetRasterBand(2).ReadAsArray()
    resampled_average = dst_ds.GetRasterBand(1).ReadAsArray()
    resampled_average[nan_frac >= thresh] = np.nan
    return resampled_average, src_ds_mem


def _setup_source(input_tif):
    """convenience setup function for gdal_average"""
    src_ds = gdal.Open(input_tif)
    data = src_ds.GetRasterBand(1).ReadAsArray()
    src_dtype = src_ds.GetRasterBand(1).DataType
    mem_driver = gdal.GetDriverByName('MEM')
    src_ds_mem = mem_driver.Create('',
                                   src_ds.RasterXSize, src_ds.RasterYSize,
                                   2, src_dtype)
    src_ds_mem.GetRasterBand(1).WriteArray(data)
    src_ds_mem.GetRasterBand(1).SetNoDataValue(0)
    # if data==0, then 1, else 0
    nan_matrix = np.isclose(data, 0, atol=1e-6)
    src_ds_mem.GetRasterBand(2).WriteArray(nan_matrix)
    src_ds_mem.SetGeoTransform(src_ds.GetGeoTransform())
    return src_ds, src_ds_mem


def _get_matlab_resampled_data_size(xscale, yscale, data):
    """convenience function mimicking the Matlab Pirate output size"""
    xscale = int(xscale)
    yscale = int(yscale)
    ysize, xsize = data.shape
    xres, yres = int(xsize / xscale), int(ysize / yscale)
    return xres, yres
