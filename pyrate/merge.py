#   This Python module is part of the PyRate software package.
#
#   Copyright 2022 Geoscience Australia
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
from os.path import join, isfile, exists
import subprocess
from pathlib import Path
from typing import Optional, Tuple
import pickle
import numpy as np
from osgeo import gdal

import pyrate.constants as C
from pyrate.core import shared, stack, ifgconstants as ifc, mpiops
from pyrate.core.logger import pyratelogger as log
from pyrate.configuration import Configuration
from pyrate.core.shared import Tile

gdal.SetCacheMax(64)


def main(params: dict) -> None:
    """
    PyRate merge main function. Assembles product tiles in to
    single geotiff files
    """
    if params[C.SIGNAL_POLARITY] == 1:
        log.info("Saving output products with the same sign convention as input data")
    elif params[C.SIGNAL_POLARITY] == -1:
        log.info("Saving output products with reversed sign convention compared to the input data")
    else:
        log.warning("Check the value of signal_polarity parameter")

    out_types = []

    tsfile = join(params[C.TMPDIR], 'tscuml_0.npy')
    if exists(tsfile):
        _merge_timeseries(params, 'tscuml')
        _merge_linrate(params)
        out_types += ['linear_rate', 'linear_error', 'linear_rsquared']

        # optional save of merged tsincr products
        if params["savetsincr"] == 1:
            _merge_timeseries(params, 'tsincr')
    else:
        log.warning(f'Not merging time series products; {tsfile} does not exist')

    stfile = join(params[C.TMPDIR], 'stack_rate_0.npy')
    if exists(stfile):
        # setup paths
        mpiops.run_once(_merge_stack, params)
        out_types += ['stack_rate', 'stack_error']
    else:
        log.warning(f'Not merging stack products; {stfile} does not exist')

    if len(out_types) > 0:
        process_out_types = mpiops.array_split(out_types)
        for out_type in process_out_types:
            create_png_and_kml_from_tif(params[C.VELOCITY_DIR], output_type=out_type)
    else:
        log.warning('Exiting: no products to merge')


def _merge_stack(params: dict) -> None:
    """
    Merge stacking outputs
    """
    shape, tiles, ifgs_dict = __merge_setup(params)

    log.info('Merging and writing Stack Rate product geotiffs')

    # read and assemble tile outputs
    rate = assemble_tiles(shape, params[C.TMPDIR], tiles, out_type='stack_rate')
    error = assemble_tiles(shape, params[C.TMPDIR], tiles, out_type='stack_error')
    samples = assemble_tiles(shape, params[C.TMPDIR], tiles, out_type='stack_samples')

    # mask pixels according to threshold
    if params[C.LR_MAXSIG] > 0:
        msg = f"Masking stack_rate and stack_error maps where error is greater " \
              f"than {params[C.LR_MAXSIG]} millimetres"
        log.info(msg)
        rate, error = stack.mask_rate(rate, error, params[C.LR_MAXSIG])
    else:
        log.info('Skipping stack product masking (maxsig = 0)')

    # save geotiff and numpy array files
    for out, otype in zip([rate, error, samples], ['stack_rate', 'stack_error', 'stack_samples']):
        __save_merged_files(ifgs_dict, params, out, otype, savenpy=params["savenpy"])


def _merge_linrate(params: dict) -> None:
    """
    Merge linear rate outputs
    """
    shape, tiles, ifgs_dict = mpiops.run_once(__merge_setup, params)

    log.info('Merging and writing Linear Rate product geotiffs')

    # read and assemble tile outputs
    out_types = ['linear_' + x for x in ['rate', 'rsquared', 'error', 'intercept', 'samples']]
    process_out_types = mpiops.array_split(out_types)
    for p_out_type in process_out_types:
        out = assemble_tiles(shape, params[C.TMPDIR], tiles, out_type=p_out_type)
        __save_merged_files(ifgs_dict, params, out, p_out_type, savenpy=params["savenpy"])
    mpiops.comm.barrier()


def _merge_timeseries(params: dict, tstype: str) -> None:
    """
    Merge tiled time series outputs
    """
    log.info(f'Merging {tstype} time series outputs')
    shape, tiles, ifgs_dict = __merge_setup(params)

    # load the first time series file to determine the number of time series tifs
    ts_file = join(params[C.TMPDIR], tstype + '_0.npy')
    ts = np.load(file=ts_file)
    # pylint: disable=no-member
    no_ts_tifs = ts.shape[2]
    process_tifs = mpiops.array_split(range(no_ts_tifs))
    # depending on nvelpar, this will not fit in memory
    # e.g. nvelpar=100, nrows=10000, ncols=10000, 32bit floats need 40GB memory
    # 32 * 100 * 10000 * 10000 / 8 bytes = 4e10 bytes = 40 GB
    # the double for loop helps us overcome the memory limit
    log.info(f'Process {mpiops.rank} writing {len(process_tifs)} {tstype} ' \
             f'time series tifs of total {no_ts_tifs}')
    for i in process_tifs:
        ts_arr = assemble_tiles(shape, params[C.TMPDIR], tiles, out_type=tstype, index=i)
        __save_merged_files(
            ifgs_dict, params, ts_arr, out_type=tstype, index=i, savenpy=params["savenpy"]
        )

    mpiops.comm.barrier()
    log.debug(f'Process {mpiops.rank} finished writing {len(process_tifs)} {tstype} time series '
              f'tifs of total {no_ts_tifs}')


def create_png_and_kml_from_tif(output_folder_path: str, output_type: str) -> None:
    """
    Function to create a preview PNG format image from a geotiff, and a KML file
    """
    log.info(f'Creating quicklook image for {output_type}')
    # open raster and choose band to find min, max
    raster_path = join(output_folder_path, f"{output_type}.tif")
    if not isfile(raster_path):
        raise Exception(f"{output_type}.tif file not found at: " + raster_path)
    gtif = gdal.Open(raster_path)
    # find bounds of image
    west, north, east, south = "", "", "", ""
    for line in gdal.Info(gtif).split('\n'):
        if "Upper Left" in line:
            west, north = line.split(")")[0].split("(")[1].split(",")
        if "Lower Right" in line:
            east, south = line.split(")")[0].split("(")[1].split(",")
    # write KML file
    kml_file_path = join(output_folder_path, f"{output_type}.kml")
    kml_file_content = f"""<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://earth.google.com/kml/2.1">
  <Document>
    <name>{output_type}.kml</name>
    <GroundOverlay>
      <name>{output_type}.png</name>
      <Icon>
        <href>{output_type}.png</href>
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
    with open(kml_file_path, "w", encoding="utf-8") as f:
        f.write(kml_file_content)

    # Get raster statistics
    srcband = gtif.GetRasterBand(1)
    minimum, maximum, _, _ = srcband.GetStatistics(True, True)
    del gtif  # close geotiff (used to calculate statistics)
    # steps used for the colourmap, must be even
    # (currently hard-coded to 254 resulting in 255 values)
    no_of_steps = 254
    # slightly different code required for rate map and rate error map
    if output_type in ('stack_rate', 'linear_rate'):
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
    if output_type in ('stack_error', 'linear_error', 'linear_rsquared'):
        # colours: white -> red (minimum error -> maximum error
        # allocate RGB values to three numpy arrays r, g, b
        r = np.ones(no_of_steps + 1) * 255
        g = np.arange(0, no_of_steps + 1) / (no_of_steps)
        g = np.flipud(g) * 255
        b = g
        # generate the colourmap file in the output folder
    cmap_path = join(output_folder_path, f"colourmap_{output_type}.txt")
    log.info(f'Saving colour map to file {cmap_path}; min/max values: {minimum:.2f}/{maximum:.2f}')

    # write colourmap as text file
    with open(cmap_path, "w", encoding="utf-8") as f:
        f.write("nan 0 0 0 0\n")
        for i, value in enumerate(np.linspace(minimum, maximum, no_of_steps + 1)):
            f.write(f"{value} {r[i]} {g[i]} {b[i]} 255\n")

    # create PNG image file
    input_tif_path = join(output_folder_path, f"{output_type}.tif")
    output_png_path = join(output_folder_path, f"{output_type}.png")
    subprocess.check_call(["gdaldem", "color-relief", "-of", "PNG", input_tif_path, "-alpha",
                           cmap_path, output_png_path, "-nearest_color_entry"])
    log.debug(f'Finished creating quicklook image for {output_type}')


def assemble_tiles(
    shape: Tuple,
    tile_dir: str,
    tiles: Tile,
    out_type: str,
    index: Optional[int] = None
) -> np.ndarray:
    """
    Function to reassemble tiles from numpy files in to a merged array

    :param s: shape for merged array.
    :param tile_dir: path to directory containing numpy tile files.
    :param tiles: pyrate.core.shared.Tile Class object.
    :param out_type: product type string, used to construct numpy tile file name.
    :param index: array third dimension index to extract from 3D time series array tiles.
    :return: merged_array: array assembled from all tiles.
    """
    log.debug(f'Re-assembling tiles for {out_type}')
    # pre-allocate dest array
    merged_array = np.empty(shape=shape, dtype=np.float32)

    # loop over each tile, load and slot in to correct spot
    for tile in tiles:
        tile_file = Path(join(tile_dir, out_type + '_' + str(tile.index) + '.npy'))
        tile_data = np.load(file=tile_file)
        slice_y = slice(tile.top_left_y, tile.bottom_right_y)
        slice_x = slice(tile.top_left_x, tile.bottom_right_x)
        if index is None:  # 2D array
            merged_array[slice_y, slice_x] = tile_data
        else:  # 3D array
            merged_array[slice_y, slice_x] = tile_data[:, :, index]

    log.debug(f'Finished assembling tiles for {out_type}')
    return merged_array


out_type_md_dict = {
    'stack_rate': ifc.STACKRATE,
    'stack_error': ifc.STACKERROR,
    'stack_samples': ifc.STACKSAMP,
    'linear_rate': ifc.LINRATE,
    'linear_error': ifc.LINERROR,
    'linear_samples': ifc.LINSAMP,
    'linear_intercept': ifc.LINICPT,
    'linear_rsquared': ifc.LINRSQ,
    'tsincr': ifc.INCR,
    'tscuml': ifc.CUML
}

error_out_types = {'linear_error', 'stack_error'}
los_projection_out_types = {'tsincr', 'tscuml', 'linear_rate', 'stack_rate'}
los_projection_divisors = {
    ifc.LINE_OF_SIGHT: lambda data: 1,
    ifc.PSEUDO_VERTICAL: np.cos,
    ifc.PSEUDO_HORIZONTAL: np.sin
}


def __save_merged_files(ifgs_dict, params, array, out_type, index=None, savenpy=None):
    """
    Convenience function to save PyRate geotiff and numpy array files
    """
    ts_dir = params[C.TIMESERIES_DIR]
    vel_dir = params[C.VELOCITY_DIR]
    log.debug(f'Saving PyRate outputs {out_type}')
    gt, md, wkt = ifgs_dict['gt'], ifgs_dict['md'], ifgs_dict['wkt']
    epochlist = ifgs_dict['epochlist']

    if out_type in ('tsincr', 'tscuml'):
        epoch = epochlist.dates[index + 1]
        dest = join(ts_dir, out_type + "_" + str(epoch) + ".tif")
        npy_file = join(ts_dir, out_type + "_" + str(epoch) + ".npy")
        # sequence position; first time slice is #0
        md['SEQUENCE_POSITION'] = index + 1
        md[ifc.EPOCH_DATE] = epoch
    else:
        dest = join(vel_dir, out_type + ".tif")
        npy_file = join(vel_dir, out_type + '.npy')
        md[ifc.EPOCH_DATE] = [d.strftime('%Y-%m-%d') for d in epochlist.dates]

    md[ifc.DATA_TYPE] = out_type_md_dict[out_type]

    # apply LOS projection and sign conversion for these outputs
    los_proj = params[C.LOS_PROJECTION]

    if out_type in los_projection_out_types:
        incidence_path = Path(Configuration.geometry_files(params)['incidence_angle'])
        if incidence_path.exists():  # We can do LOS projection
            if los_proj == ifc.LINE_OF_SIGHT:
                log.info(f"Retaining Line-of-sight signal projection for file {dest}")
            elif los_proj != ifc.LINE_OF_SIGHT:
                log.info(f"Projecting Line-of-sight signal into "
                         f"{ifc.LOS_PROJECTION_OPTION[los_proj]} for file {dest}")
            incidence = shared.Geometry(incidence_path)
            incidence.open()
            array /= los_projection_divisors[los_proj](incidence.data) * params[C.SIGNAL_POLARITY]
            md[C.LOS_PROJECTION.upper()] = ifc.LOS_PROJECTION_OPTION[los_proj]
            md[C.SIGNAL_POLARITY.upper()] = params[C.SIGNAL_POLARITY]

    if out_type in error_out_types:  # add n-sigma value to error file metadata
        # TODO: move the multiplication of nsigma and error data here?
        md[C.VELERROR_NSIG.upper()] = params[C.VELERROR_NSIG]
        log.info(f"Saving file {dest} with {params[C.VELERROR_NSIG]}-sigma uncertainties")

    shared.write_output_geotiff(md, gt, wkt, array, dest, np.nan)

    # clear the extra metadata so it does not mess up other images
    if C.LOS_PROJECTION.upper() in md:
        md.pop(C.LOS_PROJECTION.upper())
    if C.SIGNAL_POLARITY.upper() in md:
        md.pop(C.SIGNAL_POLARITY.upper())
    if C.VELERROR_NSIG.upper() in md:
        md.pop(C.VELERROR_NSIG.upper())

    if savenpy:
        np.save(file=npy_file, arr=array)

    log.debug(f'Finished saving {out_type}')


def __merge_setup(params):
    """
    Convenience function for Merge set up steps
    """
    # load previously saved preread_ifgs dict
    preread_ifgs_file = Configuration.preread_ifgs(params)
    with open(preread_ifgs_file, 'rb') as file:
        ifgs_dict = pickle.load(file)
    ifgs = [v for v in ifgs_dict.values() if isinstance(v, shared.PrereadIfg)]
    shape = ifgs[0].shape
    tiles = Configuration.get_tiles(params)
    return shape, tiles, ifgs_dict
