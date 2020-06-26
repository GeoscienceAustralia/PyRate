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
This Python module does post-processing steps to assemble the
stack rate and time series outputs and save as geotiff files
"""
from os.path import join, isfile
import pickle
import numpy as np
from osgeo import gdal
import subprocess
from pathlib import Path
from math import floor

from pyrate.core import shared, stack, ifgconstants as ifc, mpiops, config as cf
from pyrate.core.logger import pyratelogger as log

gdal.SetCacheMax(64)


def main(params):
    """
    PyRate merge main function. Assembles product tiles in to
    single geotiff files
    """
    # setup paths
    rows, cols = params["rows"], params["cols"]
    mpiops.run_once(_merge_stack, rows, cols, params)
    mpiops.run_once(_create_png_from_tif, params[cf.OUT_DIR])

    if params[cf.TIME_SERIES_CAL]:
        _merge_timeseries(rows, cols, params)
        # mpiops.run_once(_delete_tsincr_files, params)


def _merge_stack(rows, cols, params):
    """
    Merge stacking outputs
    """
    shape, tiles, ifgs_dict = _merge_setup(rows, cols, params)

    log.info('Merging and writing Stack Rate product geotiffs')

    # read and assemble tile outputs
    rate = assemble_tiles(shape, params[cf.TMPDIR], tiles, out_type='stack_rate')
    error = assemble_tiles(shape, params[cf.TMPDIR], tiles, out_type='stack_error')
    samples = assemble_tiles(shape, params[cf.TMPDIR], tiles, out_type='stack_samples')

    # mask pixels according to threshold
    if params[cf.LR_MAXSIG] > 0:
        rate, error = stack.mask_rate(rate, error, params[cf.LR_MAXSIG])
    else:
        log.info('Skipping stack product masking (maxsig = 0)')

    # save geotiff and numpy array files
    for out, ot in zip([rate, error, samples], ['stack_rate', 'stack_error', 'stack_samples']):
        _save_merged_files(ifgs_dict, params[cf.OUT_DIR], out, ot, savenpy=params["savenpy"])


def _merge_timeseries(rows, cols, params):
    """
    Merge time series output
    """
    log.info("Merging timeseries output")
    shape, tiles, ifgs_dict = _merge_setup(rows, cols, params)

    # load the first tsincr file to determine the number of time series tifs
    tsincr_file = join(params[cf.TMPDIR], 'tsincr_0.npy')
    tsincr = np.load(file=tsincr_file)
    # pylint: disable=no-member
    no_ts_tifs = tsincr.shape[2]

    # create 2 x no_ts_tifs as we are splitting tsincr and tscuml to all processes.
    process_tifs = mpiops.array_split(range(2 * no_ts_tifs))

    # depending on nvelpar, this will not fit in memory
    # e.g. nvelpar=100, nrows=10000, ncols=10000, 32bit floats need 40GB memory
    # 32 * 100 * 10000 * 10000 / 8 bytes = 4e10 bytes = 40 GB
    # the double for loop helps us overcome the memory limit
    log.info('Process {} writing {} timeseries tifs of '
             'total {}'.format(mpiops.rank, len(process_tifs), no_ts_tifs * 2))
    for i in process_tifs:
        if i < no_ts_tifs:
            tscum_g = assemble_tiles(shape, params[cf.TMPDIR], tiles, out_type='tscuml', index=i)
            _save_merged_files(ifgs_dict, params[cf.OUT_DIR], tscum_g, out_type='tscuml', index=i,
                               savenpy=params["savenpy"])
        else:
            i %= no_ts_tifs
            tsincr_g = assemble_tiles(shape, params[cf.TMPDIR], tiles, out_type='tsincr', index=i)
            _save_merged_files(ifgs_dict, params[cf.OUT_DIR], tsincr_g, out_type='tsincr', index=i,
                               savenpy=params["savenpy"])

    mpiops.comm.barrier()
    log.debug('Process {} finished writing {} timeseries tifs of '
             'total {}'.format(mpiops.rank, len(process_tifs), no_ts_tifs * 2))


def _create_png_from_tif(output_folder_path):
    """
    Wrapper for rate and error png/kml generation
    """
    create_png_and_kml_from_tif(output_folder_path, output_type='rate')
    create_png_and_kml_from_tif(output_folder_path, output_type='error')


def create_png_and_kml_from_tif(output_folder_path: str, output_type: str) -> None:
    """
    Function to create a preview PNG format image from a geotiff, and a KML file
    """
    log.info(f'Creating quicklook image for stack_{output_type}')
    # open raster and choose band to find min, max
    raster_path = join(output_folder_path, f"stack_{output_type}.tif")
    if not isfile(raster_path):
        raise Exception(f"stack_{output_type}.tif file not found at: " + raster_path)
    gtif = gdal.Open(raster_path)
    # find bounds of image
    west, north, east, south = "", "", "", ""
    for line in gdal.Info(gtif).split('\n'):
        if "Upper Left" in line:
            west, north = line.split(")")[0].split("(")[1].split(",")
        if "Lower Right" in line:
            east, south = line.split(")")[0].split("(")[1].split(",")
    # write KML file
    kml_file_path = join(output_folder_path, f"stack_{output_type}.kml")
    kml_file_content = f"""<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://earth.google.com/kml/2.1">
  <Document>
    <name>stack_{output_type}.kml</name>
    <GroundOverlay>
      <name>stack_{output_type}.png</name>
      <Icon>
        <href>stack_{output_type}.png</href>
      </Icon>
      <LatLonBox>
        <north> """ + north + """ </north>
        <south> """ + south + """ </south>
        <east>  """ + east + """ </east>
        <west>  """ + west + """ </west>
      </LatLonBox>
    </GroundOverlay>
  </Document>
</kml>"""
    with open(kml_file_path, "w") as f:
        f.write(kml_file_content)
    
    # Get raster statistics
    srcband = gtif.GetRasterBand(1)
    minimum, maximum, _, _ = srcband.GetStatistics(True, True)
    del gtif  # close geotiff (used to calculate statistics)
    # steps used for the colourmap, must be even (currently hard-coded to 254 resulting in 255 values)
    no_of_steps = 254
    # slightly different code required for rate map and rate error map
    if output_type == 'rate':
        # minimum value might be negative
        maximum = max(abs(minimum), abs(maximum))
        minimum = -1 * maximum
        # colours: blue -> white -> red (white==0) 
        # note that an extra value will be added for zero (i.e. white: 255 255 255)
        # generate a colourmap for odd number of values (currently hard-coded to 255)
        mid = int(no_of_steps * 0.5)
        # allocate RGB values to three numpy arrays r, g, b
        r = np.arange(0, mid) / mid
        g = r
        r = np.concatenate((r, np.ones(mid + 1)))
        g = np.concatenate((g, np.array([1]), np.flipud(g)))
        b = np.flipud(r)
        # change direction of colours (blue: positve, red: negative)
        r = np.flipud(r) * 255
        g = np.flipud(g) * 255
        b = np.flipud(b) * 255
    if output_type == 'error':
        # colours: white -> red (minimum error -> maximum error  
        # allocate RGB values to three numpy arrays r, g, b
        r = np.ones(no_of_steps+1)*255
        g = np.arange(0,no_of_steps+1)/(no_of_steps)
        g = np.flipud(g)*255
        b = g  
    # generate the colourmap file in the output folder   
    color_map_path = join(output_folder_path, f"colourmap_{output_type}.txt")
    log.info(
        'Saving colour map to file {}; min/max values: {:.2f}/{:.2f}'.format(color_map_path, minimum,
                                                                                            maximum))
    with open(color_map_path, "w") as f:
        f.write("nan 0 0 0 0\n")
        for i, value in enumerate(np.linspace(minimum, maximum, no_of_steps+1)):
            f.write("%f %f %f %f 255\n" % (value, r[i], g[i], b[i]))
    input_tif_path = join(output_folder_path, f"stack_{output_type}.tif")
    output_png_path = join(output_folder_path, f"stack_{output_type}.png")
    subprocess.check_call(["gdaldem", "color-relief", "-of", "PNG", input_tif_path, "-alpha",
                           color_map_path, output_png_path, "-nearest_color_entry"])
    log.debug(f'Finished creating quicklook image for stack_{output_type}')


def assemble_tiles(s, dir, tiles, out_type, index=None):
    """
    Function to reassemble tiles from numpy files in to a merged array

    :param tuple s: shape for merged array.
    :param str dir: path to directory containing numpy tile files.
    :param str out_type: product type string, used to construct numpy tile file name.
    :param int index: array third dimension index to extract from 3D time series array tiles.

    :return: merged_array: array assembled from all tiles.
    :rtype: ndarray
    """
    log.debug('Re-assembling tiles for {}'.format(out_type))
    # pre-allocate dest array
    merged_array = np.empty(shape=s, dtype=np.float32)

    # loop over each tile, load and slot in to correct spot
    for t in tiles:
        tile_file = Path(join(dir, out_type + '_'+str(t.index)+'.npy'))
        tile = np.load(file=tile_file)
        if index is None: #2D array
            merged_array[t.top_left_y:t.bottom_right_y, t.top_left_x:t.bottom_right_x] = tile
        else: #3D array
            merged_array[t.top_left_y:t.bottom_right_y, t.top_left_x:t.bottom_right_x] = tile[:, :, index]

    log.debug('Finished assembling tiles for {}'.format(out_type))
    return merged_array


def _save_merged_files(ifgs_dict, outdir, array, out_type, index=None, savenpy=None):
    """
    Convenience function to save PyRate geotiff and numpy array files
    """
    log.debug('Saving PyRate outputs {}'.format(out_type))
    gt, md, wkt = ifgs_dict['gt'], ifgs_dict['md'], ifgs_dict['wkt']
    epochlist = ifgs_dict['epochlist']

    if out_type in ('tsincr', 'tscuml'):
        epoch = epochlist.dates[index + 1]
        dest = join(outdir, out_type + "_" + str(epoch) + ".tif")
        npy_file = join(outdir, out_type + "_" + str(epoch) + ".npy")
        # sequence position; first time slice is #0
        md['SEQUENCE_POSITION'] = index+1
        md[ifc.EPOCH_DATE] = epoch
    else:
        dest = join(outdir, out_type + ".tif")
        npy_file = join(outdir, out_type + '.npy')
        md[ifc.EPOCH_DATE] = epochlist.dates

    if out_type == 'stack_rate':
        md[ifc.DATA_TYPE] = ifc.STACKRATE
    elif out_type == 'stack_error':
        md[ifc.DATA_TYPE] = ifc.STACKERROR
    elif out_type == 'stack_samples':
        md[ifc.DATA_TYPE] = ifc.STACKSAMP
    elif out_type == 'tsincr':
        md[ifc.DATA_TYPE] = ifc.INCR
    else: #tscuml
        md[ifc.DATA_TYPE] = ifc.CUML

    shared.write_output_geotiff(md, gt, wkt, array, dest, np.nan)
    if savenpy:
        np.save(file=npy_file, arr=array)

    log.debug('Finished saving {}'.format(out_type))


def _merge_setup(rows, cols, params):
    """
    Convenience function for Merge set up steps
    """
    # setup paths
    xlks, _, crop = cf.transform_params(params)
    base_unw_paths = []

    for p in Path(params[cf.OUT_DIR]).rglob("*rlks_*cr.tif"):
        if "dem" not in str(p):
            base_unw_paths.append(str(p))

    if "tif" in base_unw_paths[0].split(".")[1]:
        dest_tifs = base_unw_paths # cf.get_dest_paths(base_unw_paths, crop, params, xlks)
        for i, dest_tif in enumerate(dest_tifs):
            dest_tifs[i] = dest_tif.replace("_tif", "")
    else:
        dest_tifs = base_unw_paths # cf.get_dest_paths(base_unw_paths, crop, params, xlks)

    # load previously saved preread_ifgs dict
    preread_ifgs_file = join(params[cf.TMPDIR], 'preread_ifgs.pk')
    ifgs_dict = pickle.load(open(preread_ifgs_file, 'rb'))
    ifgs = [v for v in ifgs_dict.values() if isinstance(v, shared.PrereadIfg)]
    shape = ifgs[0].shape
    tiles = shared.get_tiles(dest_tifs[0], rows, cols)
    return shape, tiles, ifgs_dict


def _delete_tsincr_files(params):
    """
    Convenience function to delete tsincr files
    """
    out_dir = Path(params[cf.OUT_DIR])
    for file_path in out_dir.iterdir():
        if "tsincr" in str(file_path):
            file_path.unlink()
