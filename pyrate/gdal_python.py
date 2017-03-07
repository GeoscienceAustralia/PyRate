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
    the pixel location of a geospatial coordinate; from:
    http://pcjericks.github.io/py-gdalogr-cookbook/raster_layers.html#clip-a-geotiff-with-shapefile
    """
    ul_x = geo_transform[0]
    ul_y = geo_transform[3]
    res = geo_transform[1]
    pixel = int(np.round((x - ul_x) / res))
    line = int(np.round((ul_y - y) / res))
    return pixel, line


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

    # for point in range(pts.GetPointCount()):
    #     points.append((pts.GetX(point), pts.GetY(point)))

    for point in points:
        pixels.append(world_to_pixel(gt2, point[0], point[1]))

    raster_poly = Image.new('L', size=(px_width, px_height), color=1)
    rasterize = ImageDraw.Draw(raster_poly)
    rasterize.polygon(pixels, 0)  # Fill with zeroes

    # If the clipping features extend out-of-bounds and ABOVE the raster...
    # SB: we don't implement clipping for out-of-bounds and ABOVE the raster
    # We might need this down the line when we have looked at `maximum crop` in
    # detail
    # if geo_trans[3] < max_y:
    #     # The clip features were "pushed down" to match the bounds of the
    #     #   raster; this step "pulls" them back up
    #     premask = image_to_array(raster_poly)
    #     # We slice out the piece of our clip features that are "off the map"
    #     mask = np.ndarray((premask.shape[-2] - abs(iY),
    #                        premask.shape[-1]), premask.dtype)
    #     mask[:] = premask[abs(iY):, :]
    #     mask.resize(premask.shape)  # Then fill in from the bottom
    #
    #     # Most importantly, push the clipped piece down
    #     gt2[3] = max_y - (max_y - geo_trans[3])
    #
    # else:
    #     mask = image_to_array(raster_poly)
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
    raster = None  # manual close
    return clip, gt2


def resample_nearest_neighbour(input_tif, extents, new_res, output_file):
    """
    nearest neighbor resampling via GDAL
    """
    dst, resampled_proj, src, _ = crop_rasample_setup(extents,
                                                      input_tif,
                                                      new_res,
                                                      output_file)
    # Do the work
    gdal.ReprojectImage(src, dst, '', resampled_proj,
                        gdalconst.GRA_NearestNeighbour)
    return dst.ReadAsArray()


def crop_rasample_setup(extents, input_tif, new_res, output_file,
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

    px_height, px_width = gdalwarp_width_and_height(max_x, max_y, min_x, min_y,
                                                    resampled_geotrans)

    # Output / destination
    dst = gdal.GetDriverByName(dst_driver_type).Create(
        output_file, px_width, px_height, out_bands, gdalconst.GDT_Float32)
    dst.SetGeoTransform(resampled_geotrans)
    dst.SetProjection(resampled_proj)

    for k, v in meta_data.items():
        dst.SetMetadataItem(k, v)

    return dst, resampled_proj, src_ds, src_proj


def gdalwarp_width_and_height(max_x, max_y, min_x, min_y, geo_trans):
    """
    Modify pixel height and width based on world_to_pixel
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
    Crop, resample, and average a geotif
    Parameters
    ----------
    input_tif: str
        path to input geotif to resample/crop
    extents: tuple
        georeferenced extents for new file: (xfirst, yfirst, xlast, ylast)
    output_file: str
        output resampled/cropped file name
    new_res: list
        [xres, yres] Sets resolution output Ifg metadata.
    thresh: float
        nan fraction threshold
    out_driver_type: str, optional
        the output driver type. `MEM` or `GTiff`.
    match_pirate: bool, optional
        whether to match matlab pirate style resampled/croppped output
    """
    dst_ds, _, _, _ = crop_rasample_setup(
        extents, input_tif, new_res, output_file,
        out_bands=2, dst_driver_type='MEM')

    # make a temporary copy of the dst_ds for pirate style prepifg
    tmp_ds = gdal.GetDriverByName('MEM').CreateCopy('', dst_ds) \
        if (match_pirate and new_res[0]) else None

    resampled_average, src_ds_mem = \
        gdal_average(dst_ds, input_tif, thresh)
    src_dtype = src_ds_mem.GetRasterBand(1).DataType
    src_gt = src_ds_mem.GetGeoTransform()

    # write out to output geotif file
    driver = gdal.GetDriverByName(out_driver_type)

    # required to match matlab output
    if tmp_ds:
        matlab_alignment(input_tif, new_res, resampled_average, src_ds_mem,
                         src_gt, tmp_ds)

    # write final pi/pyrate GTiff
    out_ds = driver.Create(output_file, dst_ds.RasterXSize, dst_ds.RasterYSize,
                           1, src_dtype)

    out_ds.GetRasterBand(1).SetNoDataValue(np.nan)
    out_ds.GetRasterBand(1).WriteArray(resampled_average)
    out_ds.SetGeoTransform(dst_ds.GetGeoTransform())
    out_ds.SetProjection(dst_ds.GetProjection())
    # copy metadata
    for k, v in dst_ds.GetMetadata().items():
        if k == ifc.DATA_TYPE:
            out_ds.SetMetadataItem(ifc.DATA_TYPE, ifc.MULTILOOKED)
        else:
            out_ds.SetMetadataItem(k, v)

    return resampled_average, out_ds


def matlab_alignment(input_tif, new_res, resampled_average, src_ds_mem, src_gt,
                     tmp_ds):
    """
    Correction step to match python multilook/crop ouput to match that of
    matlab pirate code.

    Parameters
    ----------
    input_tif: str
        path to input geotif to resample/crop
    new_res: list
        [xres, yres] Sets resolution output Ifg metadata.
    resampled_average: ndarray
        ndarray from previous step with average resampling
        applied to phase data
    src_ds_mem: gdal.Dataset
        gdal memory dataset object
    src_gt: tuple
        geotransform tuple
    tmp_ds: gdal.Dataset
        gdal memory dataset object

    Modifies the resampled_average array in place.
    """
    src_ds = gdal.Open(input_tif)
    data = src_ds.GetRasterBand(1).ReadAsArray()
    xlooks = ylooks = int(new_res[0] / src_gt[1])
    xres, yres = get_matlab_resampled_data_size(xlooks, ylooks, data)
    nrows, ncols = resampled_average.shape
    # pirate does nearest neighbor resampling for the last
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
    Parameters
    ----------
    dst_ds: gdal.Dataset
        destination gdal dataset object
    input_tif: str
        input geotif
    thresh: float
        nan fraction threshold

    Returns
    -------
    resampled_average: ndarray
        ndarray of of ifg phase data
    src_ds_mem: gdal.Dataset
        modified in memory src_ds with nan_fraction in Band2. The nan_fraction
        is computed efficiently here in gdal in the same step as the that of
        the resampled average (band 1). This results is huge memory and
        computational efficiency.

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
    """convenience setup function"""
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


def get_matlab_resampled_data_size(xscale, yscale, data):
    """
    convenience function mimicking the matlab output size
    """
    xscale = int(xscale)
    yscale = int(yscale)
    ysize, xsize = data.shape
    xres, yres = int(xsize / xscale), int(ysize / yscale)
    return xres, yres
