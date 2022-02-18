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
import pickle
import numpy as np
from osgeo import gdal
import subprocess
from pathlib import Path
from typing import Optional, Tuple

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
        log.info(f"Saving output products with the same sign convention as input data")
    elif params[C.SIGNAL_POLARITY] == -1:
        log.info(f"Saving output products with reversed sign convention compared to the input data")
    else:
        log.warning(f"Check the value of signal_polarity parameter")

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
        log.warning('Not merging time series products; {} does not exist'.format(tsfile))

    stfile = join(params[C.TMPDIR], 'stack_rate_0.npy')
    if exists(stfile):
        # setup paths
        mpiops.run_once(_merge_stack, params)
        out_types += ['stack_rate', 'stack_error']
    else:
        log.warning('Not merging stack products; {} does not exist'.format(stfile))

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
        log.info(f"Masking stack_rate and stack_error maps where error is greater than {params[C.LR_MAXSIG]} millimetres")
        rate, error = stack.mask_rate(rate, error, params[C.LR_MAXSIG])
    else:
        log.info('Skipping stack product masking (maxsig = 0)')

    # save geotiff and numpy array files
    for out, ot in zip([rate, error, samples], ['stack_rate', 'stack_error', 'stack_samples']):
        __save_merged_files(ifgs_dict, params, out, ot, savenpy=params["savenpy"])


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
    log.info('Merging {} time series outputs'.format(tstype))
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
    log.info('Process {} writing {} {} time series tifs of '
             'total {}'.format(mpiops.rank, len(process_tifs), tstype, no_ts_tifs))
    for i in process_tifs:
        ts_arr = assemble_tiles(shape, params[C.TMPDIR], tiles, out_type=tstype, index=i)
        __save_merged_files(ifgs_dict, params, ts_arr, out_type=tstype, index=i, savenpy=params["savenpy"])

    mpiops.comm.barrier()
    log.debug('Process {} finished writing {} {} time series tifs of '
              'total {}'.format(mpiops.rank, len(process_tifs), tstype, no_ts_tifs))


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
    with open(kml_file_path, "w") as f:
        f.write(kml_file_content)

    # Get raster statistics
    srcband = gtif.GetRasterBand(1)
    minimum, maximum, _, _ = srcband.GetStatistics(True, True)
    del gtif  # close geotiff (used to calculate statistics)
    # steps used for the colourmap, must be even (currently hard-coded to 254 resulting in 255 values)
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
    color_map_path = join(output_folder_path, f"colourmap_{output_type}.txt")
    log.info('Saving colour map to file {}; min/max values: {:.2f}/{:.2f}'.format(
        color_map_path, minimum, maximum))

    # write colourmap as text file
    with open(color_map_path, "w") as f:
        f.write("nan 0 0 0 0\n")
        for i, value in enumerate(np.linspace(minimum, maximum, no_of_steps + 1)):
            f.write("%f %f %f %f 255\n" % (value, r[i], g[i], b[i]))

    # create PNG image file
    input_tif_path = join(output_folder_path, f"{output_type}.tif")
    output_png_path = join(output_folder_path, f"{output_type}.png")
    subprocess.check_call(["gdaldem", "color-relief", "-of", "PNG", input_tif_path, "-alpha",
                           color_map_path, output_png_path, "-nearest_color_entry"])
    log.debug(f'Finished creating quicklook image for {output_type}')


def assemble_tiles(s: Tuple, dir: str, tiles: Tile, out_type: str, index: Optional[int] = None) -> np.ndarray:
    """
    Function to reassemble tiles from numpy files in to a merged array

    :param s: shape for merged array.
    :param dir: path to directory containing numpy tile files.
    :param tiles: pyrate.core.shared.Tile Class object.
    :param out_type: product type string, used to construct numpy tile file name.
    :param index: array third dimension index to extract from 3D time series array tiles.
    :return: merged_array: array assembled from all tiles.
    """
    log.debug('Re-assembling tiles for {}'.format(out_type))
    # pre-allocate dest array
    merged_array = np.empty(shape=s, dtype=np.float32)

    # loop over each tile, load and slot in to correct spot
    for t in tiles:
        tile_file = Path(join(dir, out_type + '_' + str(t.index) + '.npy'))
        tile = np.load(file=tile_file)
        if index is None:  # 2D array
            merged_array[t.top_left_y:t.bottom_right_y, t.top_left_x:t.bottom_right_x] = tile
        else:  # 3D array
            merged_array[t.top_left_y:t.bottom_right_y, t.top_left_x:t.bottom_right_x] = tile[:, :, index]

    log.debug('Finished assembling tiles for {}'.format(out_type))
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
    ifc.PSEUDO_VERTICAL: lambda data: np.cos(data),
    ifc.PSEUDO_HORIZONTAL: lambda data: np.sin(data)
}


def __save_merged_files(ifgs_dict, params, array, out_type, index=None, savenpy=None):
    """
    Convenience function to save PyRate geotiff and numpy array files
    """
    outdir = params[C.OUT_DIR]
    ts_dir = params[C.TIMESERIES_DIR]
    vel_dir = params[C.VELOCITY_DIR]
    log.debug('Saving PyRate outputs {}'.format(out_type))
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

    if out_type in los_projection_out_types:  # apply LOS projection and sign conversion for these outputs
        incidence_path = Path(Configuration.geometry_files(params)['incidence_angle'])
        if incidence_path.exists():  # We can do LOS projection
            if params[C.LOS_PROJECTION] == ifc.LINE_OF_SIGHT:
                log.info(f"Retaining Line-of-sight signal projection for file {dest}")
            elif params[C.LOS_PROJECTION] != ifc.LINE_OF_SIGHT:
                log.info(f"Projecting Line-of-sight signal into "
                         f"{ifc.LOS_PROJECTION_OPTION[params[C.LOS_PROJECTION]]} for file {dest}")
            incidence = shared.Geometry(incidence_path)
            incidence.open()
            array /= los_projection_divisors[params[C.LOS_PROJECTION]](incidence.data) * params[C.SIGNAL_POLARITY]
            md[C.LOS_PROJECTION.upper()] = ifc.LOS_PROJECTION_OPTION[params[C.LOS_PROJECTION]]
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

    log.debug('Finished saving {}'.format(out_type))


def __merge_setup(params):
    """
    Convenience function for Merge set up steps
    """
    # load previously saved preread_ifgs dict
    preread_ifgs_file = Configuration.preread_ifgs(params)
    ifgs_dict = pickle.load(open(preread_ifgs_file, 'rb'))
    ifgs = [v for v in ifgs_dict.values() if isinstance(v, shared.PrereadIfg)]
    shape = ifgs[0].shape
    tiles = Configuration.get_tiles(params)
    return shape, tiles, ifgs_dict
