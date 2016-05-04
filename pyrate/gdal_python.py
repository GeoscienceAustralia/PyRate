from osgeo import gdal, gdalnumeric, gdalconst
from PIL import Image, ImageDraw
import os
import numpy as np
gdal.SetCacheMax(2**15)
GDAL_WARP_MEMORY_LIMIT = 2**10


def world_to_pixel(geo_transform, x, y):
    '''
    Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate
    the pixel location of a geospatial coordinate; from:
    http://pcjericks.github.io/py-gdalogr-cookbook/raster_layers.html#clip-a-geotiff-with-shapefile
    '''
    ulX = geo_transform[0]
    ulY = geo_transform[3]
    res = geo_transform[1]
    pixel = int(np.round((x - ulX) / res))
    line = int(np.round((ulY - y) / res))
    return pixel, line


def crop(input_file, extents, gt=None, nodata=np.nan):
    '''

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
        gt              An optional GDAL GeoTransform to use instead
        nodata          The NoData value; defaults to -9999.
    '''

    def image_to_array(i):
        '''
        Converts a Python Imaging Library (PIL) array to a gdalnumeric image.
        '''
        a = gdalnumeric.fromstring(i.tobytes(), 'b')
        a.shape = i.im.size[1], i.im.size[0]
        return a

    raster = gdal.Open(input_file)
    # Can accept either a gdal.Dataset or numpy.array instance
    if not isinstance(raster, np.ndarray):
        if not gt:
            gt = raster.GetGeoTransform()
        raster = raster.ReadAsArray()
    else:
        if not gt:
            raise

    # Convert the layer extent to image pixel coordinates
    minX, minY, maxX, maxY = extents
    ulX, ulY = world_to_pixel(gt, minX, maxY)
    lrX, lrY = world_to_pixel(gt, maxX, minY)

    # Calculate the pixel size of the new image
    pxWidth = int(lrX - ulX)
    pxHeight = int(lrY - ulY)

    # If the clipping features extend out-of-bounds and ABOVE the raster...
    if gt[3] < maxY:
        # In such a case... ulY ends up being negative--can't have that!
        iY = ulY
        ulY = 0

    # Multi-band image?
    try:
        clip = raster[:, ulY:lrY, ulX:lrX]

    except IndexError:
        clip = raster[ulY:lrY, ulX:lrX]

    # Create a new geomatrix for the image
    gt2 = list(gt)
    gt2[0] = minX
    gt2[3] = maxY

    # Map points to pixels for drawing the boundary on a blank 8-bit,
    #   black and white, mask image.
    points = [(minX, minY), (maxX, minY), (maxX, maxY), (minY, maxY)]
    pixels = []

    # for p in range(pts.GetPointCount()):
    #     points.append((pts.GetX(p), pts.GetY(p)))

    for p in points:
        pixels.append(world_to_pixel(gt2, p[0], p[1]))

    raster_poly = Image.new('L', size=(pxWidth, pxHeight), color=1)
    rasterize = ImageDraw.Draw(raster_poly)
    rasterize.polygon(pixels, 0)  # Fill with zeroes

    # If the clipping features extend out-of-bounds and ABOVE the raster...
    # SB: we don't implement clipping for out-of-bounds and ABOVE the raster
    # We might need this down the line when we have looked at `maximum crop` in
    # detail
    # if gt[3] < maxY:
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
    #     gt2[3] = maxY - (maxY - gt[3])
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
    dst, resampled_proj, src, src_proj = crop_rasample_setup(extents,
                                                                input_tif,
                                                                new_res,
                                                                output_file)
    # Do the work
    gdal.ReprojectImage(src, dst, src_proj, resampled_proj,
                        gdalconst.GRA_NearestNeighbour)
    return dst.ReadAsArray()


def crop_rasample_setup(extents, input_tif, new_res, output_file,
                        dst_driver_type='GTiff', out_bands=2):
    # Source
    src_ds = gdal.Open(input_tif, gdalconst.GA_ReadOnly)
    src_proj = src_ds.GetProjection()

    # source metadata to be copied into the output
    md = src_ds.GetMetadata()

    # get the image extents
    minX, minY, maxX, maxY = extents
    gt = src_ds.GetGeoTransform()  # gt is tuple of 6 numbers

    # Create a new geotransform for the image
    gt2 = list(gt)
    gt2[0] = minX
    gt2[3] = maxY
    # We want a section of source that matches this:
    resampled_proj = src_proj
    if new_res[0]:  # if new_res is not None, it can't be zero either
        resampled_geotrans = gt2[:1] + [new_res[0]] + gt2[2:-1] + [new_res[1]]
    else:
        resampled_geotrans = gt2

    px_height, px_width = get_width_and_height(maxX, maxY, minX, minY,
                                               resampled_geotrans)

    # Output / destination
    dst = gdal.GetDriverByName(dst_driver_type).Create(
        output_file, px_width, px_height, out_bands, gdalconst.GDT_Float32)
    dst.SetGeoTransform(resampled_geotrans)
    dst.SetProjection(resampled_proj)

    for k, v in md.iteritems():
        dst.SetMetadataItem(k, v)

    return dst, resampled_proj, src_ds, src_proj


def get_width_and_height(maxX, maxY, minX, minY, resampled_geotrans):
    # modified image extents
    ulX, ulY = world_to_pixel(resampled_geotrans, minX, maxY)
    lrX, lrY = world_to_pixel(resampled_geotrans, maxX, minY)
    # Calculate the pixel size of the new image
    px_width = int(lrX - ulX)
    px_height = int(lrY - ulY)
    return px_height, px_width  # this is the same as `gdalwarp`


def crop_and_resample_average(
        input_tif, extents, new_res, output_file, thresh):
    dst_ds, resampled_proj, src_ds, src_proj = crop_rasample_setup(
        extents, input_tif, new_res, output_file,
        out_bands=2, dst_driver_type='MEM')

    src_ds = gdal.Open(input_tif)
    data = src_ds.GetRasterBand(1).ReadAsArray()

    mem_driver = gdal.GetDriverByName('MEM')
    src_ds_mem = mem_driver.Create('',
                                   src_ds.RasterXSize, src_ds.RasterYSize,
                                   2, gdalconst.GDT_Float32)

    src_ds_mem.GetRasterBand(1).WriteArray(data)
    src_ds_mem.GetRasterBand(1).SetNoDataValue(0)

    # if data==0, then 1, else 0
    nan_matrix = np.isclose(data, 0, atol=1e-6)
    src_ds_mem.GetRasterBand(2).WriteArray(nan_matrix)
    src_ds_mem.GetRasterBand(2).SetNoDataValue(-100000)
    src_ds_mem.SetGeoTransform(src_ds.GetGeoTransform())

    gdal.ReprojectImage(src_ds_mem, dst_ds, '', '', gdal.GRA_Average)

    # dst_ds band2 average is our nan_fraction matrix
    nan_frac = dst_ds.GetRasterBand(2).ReadAsArray()
    resampled_average = dst_ds.GetRasterBand(1).ReadAsArray()
    resampled_average[nan_frac >= thresh] = np.nan

    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(output_file, dst_ds.RasterXSize, dst_ds.RasterYSize,
                      1, gdalconst.GDT_Float32)

    out_ds.GetRasterBand(1).WriteArray(resampled_average)
    return resampled_average

if __name__ == '__main__':

    # import subprocess
    # subprocess.check_call(['gdalwarp', '-overwrite', '-srcnodata', 'None', '-te', '150.91', '-34.229999976', '150.949166651', '-34.17', '-q', 'out/20070709-20070813_utm.tif', '/tmp/test.tif'])
    # clipped_ref = gdal.Open('/tmp/test.tif').ReadAsArray()
    #
    # rast = gdal.Open('out/20070709-20070813_utm.tif')
    # extents = (150.91, -34.229999976, 150.949166651, -34.17)
    # clipped = crop(rast, extents)[0]
    # np.testing.assert_array_almost_equal(clipped_ref, clipped)
    #
    # # 2nd gdalwarp call
    # subprocess.check_call(['gdalwarp', '-overwrite', '-srcnodata', 'None', '-te', '150.91', '-34.229999976', '150.949166651', '-34.17', '-q', '-tr', '0.001666666', '-0.001666666', 'out/20060619-20061002_utm.tif', '/tmp/test2.tif'])
    # resampled_ds = gdal.Open('/tmp/test2.tif')
    # resampled_ref = resampled_ds.ReadAsArray()
    #
    #
    # rast = gdal.Open('out/20060619-20061002_utm.tif')
    # new_res = (0.001666666, -0.001666666)
    # # adfGeoTransform[0] /* top left x */
    # # adfGeoTransform[1] /* w-e pixel resolution */
    # # adfGeoTransform[2] /* 0 */
    # # adfGeoTransform[3] /* top left y */
    # # adfGeoTransform[4] /* 0 */
    # # adfGeoTransform[5] /* n-s pixel resolution (negative value) */
    # resampled = resample_nearest_neighbour('out/20060619-20061002_utm.tif', extents, new_res, '/tmp/resampled.tif')
    # np.testing.assert_array_almost_equal(resampled_ref, resampled)
    from pyrate.tests import common
    from pyrate import prepifg
    import tempfile
    sydney_test_ifgs = common.sydney_data_setup()
    # minX, minY, maxX, maxY = extents
    extents = [150.91, -34.229999976, 150.949166651, -34.17]
    extents_str = [str(e) for e in extents]
    x_looks = 2
    resolutions = [0.001666666, .001]
    for res in resolutions:
        print res
        for s in sydney_test_ifgs:
            print s.data_path
            # Old prepifg crop + resample + averaging
            # manual old prepifg style average with nearest neighbour
            averaged_data = prepifg.warp_old(
                s, x_looks, x_looks, extents_str, [res, -res], thresh=0.5,
                crop_out=4, verbose=False, ret_ifg=True).phase_data
            looks_path = prepifg.mlooked_path(s.data_path, x_looks,
                                              crop_out=4)
            os.remove(looks_path)
            resampled_temp_tif = tempfile.mktemp(suffix='.tif',
                                                prefix='resampled_')
            cropped_and_averaged = crop_and_resample_average(
                s.data_path, extents, [res, -res],
                resampled_temp_tif, thresh=0.5)
            print 'nan count', np.sum(np.isnan(averaged_data)), \
                np.sum(np.isnan(cropped_and_averaged))
            print os.path.exists(resampled_temp_tif), 'path exists?'
            # np.testing.assert_array_almost_equal(averaged_data,
            #                                      cropped_and_averaged)
            #os.remove(resampled_temp_tif)
