#   This Python module is part of the PyRate software package.
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
"""
This Python module contains bindings for the GDAL library
"""

import numexpr as ne
import numpy as np
from PIL import Image, ImageDraw
from osgeo import gdal
from osgeo import gdal_array
from osgeo import gdalconst

from core import shared, ifgconstants as ifc, prepifg_helper

gdal.SetCacheMax(2 ** 15)
GDAL_WARP_MEMORY_LIMIT = 2 ** 10
LOW_FLOAT32 = np.finfo(np.float32).min * 1e-10


def coherence_masking(src_ds, coherence_ds, coherence_thresh):
    """Perform coherence masking on raster in-place.
    
    Based on gdal_calc formula provided by Nahidul:
    gdal_calc.py -A 20151127-20151209_VV_8rlks_flat_eqa.cc.tif
     -B 20151127-20151209_VV_8rlks_eqa.unw.tif
     --outfile=test_v1.tif --calc="B*(A>=0.8)-999*(A<0.8)"
     --NoDataValue=-999

    Args:
      ds: The interferogram to mask as GDAL dataset.
      coherence_ds: The coherence GDAL dataset.
      coherence_thresh: The coherence threshold.
      src_ds: 

    Returns:

    """
    coherence_band = coherence_ds.GetRasterBand(1)
    src_band = src_ds.GetRasterBand(1)
    # ndv = src_band.GetNoDataValue()
    ndv = np.nan
    coherence = coherence_band.ReadAsArray()
    src = src_band.ReadAsArray()
    var = {"coh": coherence, "src": src, "t": coherence_thresh, "ndv": ndv}
    formula = "where(coh>=t, src, ndv)"
    res = ne.evaluate(formula, local_dict=var)
    src_band.WriteArray(res)


def world_to_pixel(geo_transform, x, y):
    """Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate
    the pixel location of a geospatial coordinate;
    see: http://pcjericks.github.io/py-gdalogr-cookbook/raster_layers.html

    Args:
      list: geo_transform: Affine transformation coefficients
      float: x: longitude coordinate
      float: y: latitude coordinate
      geo_transform: param x:
      y: returns: col: pixel column number
      x: 

    Returns:
      int: col: pixel column number

    """
    ul_x = geo_transform[0]
    ul_y = geo_transform[3]
    xres = geo_transform[1]
    yres = geo_transform[5]
    col = int(np.round((x - ul_x) / xres))
    line = int(np.round((ul_y - y) / abs(yres)))  # yres has negative size

    return col, line


def crop(input_file, extents, geo_trans=None, nodata=np.nan):
    """Adapted from http://karthur.org/2015/clipping-rasters-in-python.html
    
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

    Args:
      str: input_file
      list: extents
      list: geo_trans
      int: nodata
      input_file: param extents
      geo_trans: Default value
      nodata: Default value
      extents: 

    Returns:
      ndarray: clip: cropped part of the image

    """

    def image_to_array(i):
        """Converts a Python Imaging Library (PIL) array to a gdalnumeric image.

        Args:
          i: 

        Returns:

        """
        arr = gdal_array.fromstring(i.tobytes(), "b")
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
            raise ValueError("geo transform must be supplied")

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

    raster_poly = Image.new("L", size=(px_width, px_height), color=1)
    rasterize = ImageDraw.Draw(raster_poly)
    rasterize.polygon(pixels, 0)  # Fill with zeroes

    mask = image_to_array(raster_poly)

    # Clip the image using the mask
    try:
        clip = gdal_array.choose(mask, (clip, nodata))

    # If the clipping features extend out-of-bounds and BELOW the raster...
    except ValueError:
        # We have to cut the clipping features to the raster!
        rshp = list(mask.shape)
        if mask.shape[-2] != clip.shape[-2]:
            rshp[0] = clip.shape[-2]

        if mask.shape[-1] != clip.shape[-1]:
            rshp[1] = clip.shape[-1]

        mask.resize(*rshp, refcheck=False)

        clip = gdal_array.choose(mask, (clip, nodata))

    # AttributeError: 'numpy.ndarray' object has no attribute 'close'
    # raster.close()
    raster = None

    return clip, gt2


def resample_nearest_neighbour(input_tif, extents, new_res, output_file):
    """Nearest neighbor resampling and cropping of an image.

    Args:
      str: input_tif: input geotiff file path
      list: extents: new extents for cropping
      float: new_res: new resolution for resampling
      str: output_file: output geotiff file path
      input_tif: param extents:
      new_res: param output_file:
      extents: 
      output_file: 

    Returns:
      ndarray: dst: resampled image

    """
    dst, resampled_proj, src, _ = _crop_resample_setup(extents, input_tif, new_res, output_file)
    # Do the work
    gdal.ReprojectImage(src, dst, "", resampled_proj, gdalconst.GRA_NearestNeighbour)
    return dst.ReadAsArray()


def _crop_resample_setup(extents, input_tif, new_res, output_file, dst_driver_type="GTiff", out_bands=2):
    """Convenience function for crop/resample setup

    Args:
      extents: param input_tif:
      new_res: param output_file:
      dst_driver_type: Default value = "GTiff")
      out_bands: Default value = 2)
      input_tif: 
      output_file: 

    Returns:

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

    px_height, px_width = _gdalwarp_width_and_height(max_x, max_y, min_x, min_y, resampled_geotrans)

    # Output / destination
    dst = gdal.GetDriverByName(dst_driver_type).Create(output_file, px_width, px_height, out_bands,
                                                       gdalconst.GDT_Float32)
    dst.SetGeoTransform(resampled_geotrans)
    dst.SetProjection(resampled_proj)

    for k, v in meta_data.items():
        dst.SetMetadataItem(k, v)

    return dst, resampled_proj, src_ds, src_proj


def _gdalwarp_width_and_height(max_x, max_y, min_x, min_y, geo_trans):
    """Modify pixel height and width

    Args:
      max_x: param max_y:
      min_x: param min_y:
      geo_trans: 
      max_y: 
      min_y: 

    Returns:

    """
    # modified image extents
    ul_x, ul_y = world_to_pixel(geo_trans, min_x, max_y)
    lr_x, lr_y = world_to_pixel(geo_trans, max_x, min_y)
    # Calculate the pixel size of the new image
    px_width = int(lr_x - ul_x)
    px_height = int(lr_y - ul_y)
    return px_height, px_width  # this is the same as `gdalwarp`


def crop_resample_average(input_raster, extents, resolution, output_file, thresh, header=None, coherence_path=None,
                          coherence_thresh=None):
    """Crop, resample, and average a geotiff image.

    Args:
      str: input_raster: Path to input geotiff to resample/crop
      tuple: extents: Cropping extents (xfirst, yfirst, xlast, ylast)
      list: resolution: [xres, yres] resolution of output image
      str: output_file: Path to output resampled/cropped geotiff
      float: thresh: NaN fraction threshold
      str: out_driver_type: The output driver; `MEM` or `GTiff` (optional)
      bool: match_pyrate: Match Legacy output (optional)
      dict: header: dictionary of metadata
      input_raster: param extents:
      resolution: param output_file:
      thresh: param header:  (Default value = None)
      coherence_path: Default value = None)
      coherence_thresh: Default value = None)
      extents: 
      output_file: 
      header: (Default value = None)

    Returns:
      ndarray: resampled_average: output cropped and resampled image

    """
    out_driver_type = "GTiff"
    match_pyrate = False
    dst_ds, _, _, _ = _crop_resample_setup(extents, input_raster, resolution, output_file, out_bands=2,
                                           dst_driver_type="MEM")

    # make a temporary copy of the dst_ds for PyRate style prepifg
    tmp_ds = gdal.GetDriverByName("MEM").CreateCopy("", dst_ds) if (match_pyrate and resolution[0]) else None

    src_ds, src_ds_mem = _setup_source(input_raster)

    if coherence_path and coherence_thresh:
        coherence_raster = prepifg_helper.dem_or_ifg(coherence_path)
        coherence_raster.open()
        coherence_ds = coherence_raster.dataset
        coherence_masking(src_ds_mem, coherence_ds, coherence_thresh)
    elif coherence_path and not coherence_thresh:
        raise ValueError(
            f"Coherence file provided without a coherence "
            f"threshold. Please ensure you provide 'cohthresh' "
            f"in your config if coherence masking is enabled."
        )

    resampled_average, src_ds_mem = gdal_average(dst_ds, src_ds, src_ds_mem, thresh)
    src_dtype = src_ds_mem.GetRasterBand(1).DataType
    src_gt = src_ds_mem.GetGeoTransform()

    # required to match Legacy output
    if tmp_ds:
        _alignment(input_raster, resolution, resampled_average, src_ds_mem, src_gt, tmp_ds)

    # grab metadata from existing geotiff
    gt = dst_ds.GetGeoTransform()
    wkt = dst_ds.GetProjection()

    # TEST HERE IF EXISTING FILE HAS PYRATE METADATA. IF NOT ADD HERE
    if not ifc.DATA_TYPE in dst_ds.GetMetadata() and header is not None:
        md = shared.collate_metadata(header)
    else:
        md = dst_ds.GetMetadata()

    # update metadata for output
    for k, v in md.items():
        if k == ifc.DATA_TYPE:
            # update data type metadata
            if v == ifc.ORIG and coherence_path:
                md.update({ifc.DATA_TYPE: ifc.COHERENCE})
            elif v == ifc.ORIG and not coherence_path:
                md.update({ifc.DATA_TYPE: ifc.MULTILOOKED})
            elif v == ifc.DEM:
                md.update({ifc.DATA_TYPE: ifc.MLOOKED_DEM})
            elif v == ifc.INCIDENCE:
                md.update({ifc.DATA_TYPE: ifc.MLOOKED_INC})
            elif v == ifc.COHERENCE and coherence_path:
                pass
            elif v == ifc.MULTILOOKED and coherence_path:
                md.update({ifc.DATA_TYPE: ifc.COHERENCE})
            elif v == ifc.MULTILOOKED and not coherence_path:
                pass
            else:
                raise TypeError("Data Type metadata not recognised")

    # In-memory GDAL driver doesn't support compression so turn it off.
    creation_opts = ["compress=packbits"] if out_driver_type != "MEM" else []
    out_ds = shared.gdal_dataset(
        output_file,
        dst_ds.RasterXSize,
        dst_ds.RasterYSize,
        driver=out_driver_type,
        bands=1,
        dtype=src_dtype,
        metadata=md,
        crs=wkt,
        geotransform=gt,
        creation_opts=creation_opts,
    )

    shared.write_geotiff(resampled_average, out_ds, np.nan)


def _alignment(input_tif, new_res, resampled_average, src_ds_mem, src_gt, tmp_ds):
    """Correction step to match python multi-look/crop output to match that of
    Legacy data. Modifies the resampled_average array in place.

    Args:
      input_tif: param new_res:
      resampled_average: param src_ds_mem:
      src_gt: param tmp_ds:
      new_res: 
      src_ds_mem: 
      tmp_ds: 

    Returns:

    """
    src_ds = gdal.Open(input_tif)
    data = src_ds.GetRasterBand(1).ReadAsArray()
    xlooks = ylooks = int(new_res[0] / src_gt[1])
    xres, yres = _get_resampled_data_size(xlooks, ylooks, data)
    nrows, ncols = resampled_average.shape
    # Legacy nearest neighbor resampling for the last
    # [yres:nrows, xres:ncols] cells without nan_conversion
    # turn off nan-conversion
    src_ds_mem.GetRasterBand(1).SetNoDataValue(LOW_FLOAT32)
    # nearest neighbor resapling
    gdal.ReprojectImage(src_ds_mem, tmp_ds, "", "", gdal.GRA_NearestNeighbour)
    # only take the [yres:nrows, xres:ncols] slice
    if nrows > yres or ncols > xres:
        resampled_nearest_neighbor = tmp_ds.GetRasterBand(1).ReadAsArray()
        resampled_average[yres - nrows:, xres - ncols:] = resampled_nearest_neighbor[yres - nrows:, xres - ncols:]


def gdal_average(dst_ds, src_ds, src_ds_mem, thresh):
    """Perform subsampling of an image by averaging values

    Args:
      gdal: Dataset dst_ds: Destination gdal dataset object
      str: input_tif: Input geotif
      float: thresh: NaN fraction threshold
    :return resampled_average: resampled image data
      dst_ds: param src_ds:
      src_ds_mem: param thresh:
      src_ds: 
      thresh: 

    Returns:

    """
    src_ds_mem.GetRasterBand(2).SetNoDataValue(-100000)
    src_gt = src_ds.GetGeoTransform()
    src_ds_mem.SetGeoTransform(src_gt)
    gdal.ReprojectImage(src_ds_mem, dst_ds, "", "", gdal.GRA_Average)
    # dst_ds band2 average is our nan_fraction matrix
    nan_frac = dst_ds.GetRasterBand(2).ReadAsArray()
    resampled_average = dst_ds.GetRasterBand(1).ReadAsArray()
    resampled_average[nan_frac >= thresh] = np.nan
    return resampled_average, src_ds_mem


def _setup_source(input_tif):
    """convenience setup function for gdal_average

    Args:
      input_tif: 

    Returns:

    """
    src_ds = gdal.Open(input_tif)
    data = src_ds.GetRasterBand(1).ReadAsArray()
    src_dtype = src_ds.GetRasterBand(1).DataType
    mem_driver = gdal.GetDriverByName("MEM")
    src_ds_mem = mem_driver.Create("", src_ds.RasterXSize, src_ds.RasterYSize, 2, src_dtype)
    src_ds_mem.GetRasterBand(1).WriteArray(data)
    src_ds_mem.GetRasterBand(1).SetNoDataValue(0)
    # if data==0, then 1, else 0
    nan_matrix = np.isclose(data, 0, atol=1e-6)
    src_ds_mem.GetRasterBand(2).WriteArray(nan_matrix)
    src_ds_mem.SetGeoTransform(src_ds.GetGeoTransform())
    return src_ds, src_ds_mem


def _get_resampled_data_size(xscale, yscale, data):
    """convenience function mimicking the Legacy output size

    Args:
      xscale: param yscale:
      data: 
      yscale: 

    Returns:

    """
    xscale = int(xscale)
    yscale = int(yscale)
    ysize, xsize = data.shape
    xres, yres = int(xsize / xscale), int(ysize / yscale)
    return xres, yres
