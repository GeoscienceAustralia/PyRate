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
rate and time series outputs and save as geotiff files
"""
import os
from os.path import join
import pickle as cp
import numpy as np
from osgeo import gdal
import subprocess
from pathlib import Path

from pyrate.core import shared, ifgconstants as ifc, mpiops, config as cf
from pyrate.core.shared import PrereadIfg
from pyrate.constants import REF_COLOR_MAP_PATH
from pyrate.core.config import OBS_DIR, OUT_DIR, ConfigException
from pyrate.core.logger import pyratelogger as log

gdal.SetCacheMax(64)


# Constants
MASTER_PROCESS = 0


def main(params):
    """
    PyRate merge main function. Assembles product tiles in to
    single geotiff files
    """
    # setup paths
    rows, cols = params["rows"], params["cols"]
    _merge_stack(rows, cols, params)

    if params[cf.TIME_SERIES_CAL]:
        _merge_timeseries(rows, cols, params)
        #mpiops.run_once(_delete_tsincr_files, params)

    log.info('Creating quicklook images.')
    mpiops.run_once(create_png_from_tif, params[cf.OUT_DIR])
    log.debug('Finished creating quicklook images.')


def create_png_from_tif(output_folder_path):

    # open raster and choose band to find min, max
    raster_path = os.path.join(output_folder_path, "stack_rate.tif")

    if not os.path.isfile(raster_path):
        raise Exception("stack_rate.tif file not found at: "+raster_path)
    gtif = gdal.Open(raster_path)
    srcband = gtif.GetRasterBand(1)

    west, north, east, south = "", "", "", ""
    for line in gdal.Info(gtif).split('\n'):
        if "Upper Left" in line:
            west, north = line.split(")")[0].split("(")[1].split(",")
        if "Lower Right" in line:
            east, south = line.split(")")[0].split("(")[1].split(",")

    kml_file_path = os.path.join(output_folder_path, "stack_rate.kml")
    kml_file_content = """<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://earth.google.com/kml/2.1">
  <Document>
    <name>stack_rate.kml</name>
    <GroundOverlay>
      <name>stack_rate.png</name>
      <Icon>
        <href>stack_rate.png</href>
      </Icon>
      <LatLonBox>
        <north> """+north+""" </north>
        <south> """+south+""" </south>
        <east>  """+east+""" </east>
        <west>  """+west+""" </west>
      </LatLonBox>
    </GroundOverlay>
  </Document>
</kml>"""

    with open(kml_file_path, "w") as f:
        f.write(kml_file_content)

    # Get raster statistics
    minimum, maximum, mean, stddev = srcband.GetStatistics(True, True)
    maximum = max(abs(minimum), abs(maximum))
    minimum = -1 * maximum
    step = (maximum - minimum) / 256.0

    del gtif  # manually close raster

    # read color map from utilities and write it to the output folder

    with open(REF_COLOR_MAP_PATH, "r") as f:
        color_map_list = []
        for line in f.readlines():
            color_map_list.append(line.strip().split(" "))

    no_of_data_value = len(np.arange(minimum, maximum, step))
    for i, no in enumerate(np.arange(minimum, maximum, step)):
        color_map_list[i+1][0] = str(no)

    color_map_path = os.path.join(output_folder_path, "colourmap.txt")
    with open(color_map_path, "w") as f:
        for i in range(no_of_data_value):
            f.write(' '.join(color_map_list[i]) + "\n")

    input_tif_path = os.path.join(output_folder_path, "stack_rate.tif")
    output_png_path = os.path.join(output_folder_path, "stack_rate.png")
    subprocess.check_call(["gdaldem", "color-relief", "-of", "PNG", input_tif_path, "-alpha",
                           color_map_path, output_png_path, "-nearest_color_entry"])


def _merge_stack(rows, cols, params):
    """
    Merge stacking outputs
    """
    # setup paths
    xlks, _, crop = cf.transform_params(params)
    base_unw_paths = []

    for p in Path(params[OUT_DIR]).rglob("*rlks_*cr.tif"):
        if "dem" not in str(p):
            base_unw_paths.append(str(p))

    if not base_unw_paths:
        raise ConfigException(f"No rlooked files available in {params[OBS_DIR]} or {params[OBS_DIR]}")

    if "tif" in base_unw_paths[0].split(".")[1]:
        dest_tifs = base_unw_paths # cf.get_dest_paths(base_unw_paths, crop, params, xlks)
        for i, dest_tif in enumerate(dest_tifs):
            dest_tifs[i] = dest_tif.replace("_tif","")
    else:
        dest_tifs = base_unw_paths # cf.get_dest_paths(base_unw_paths, crop, params, xlks)

    # load previously saved prepread_ifgs dict
    preread_ifgs_file = join(params[cf.TMPDIR], 'preread_ifgs.pk')
    ifgs = cp.load(open(preread_ifgs_file, 'rb'))
    tiles = shared.get_tiles(dest_tifs[0], rows, cols)

    # stacking aggregation
    if mpiops.size >= 3:
        [_save_stack(ifgs, params, tiles, out_type=t)
         for i, t in enumerate(['stack_rate', 'stack_error', 'stack_samples'])
         if i == mpiops.rank]
    else:
        if mpiops.rank == MASTER_PROCESS:
            [_save_stack(ifgs, params, tiles, out_type=t)
             for t in ['stack_rate', 'stack_error', 'stack_samples']]
    mpiops.comm.barrier()


def _save_stack(ifgs_dict, params, tiles, out_type):
    """
    Save stacking outputs
    """
    log.info('Merging PyRate outputs {}'.format(out_type))
    gt, md, wkt = ifgs_dict['gt'], ifgs_dict['md'], ifgs_dict['wkt']
    epochlist = ifgs_dict['epochlist']
    ifgs = [v for v in ifgs_dict.values() if isinstance(v, PrereadIfg)]
    dest = os.path.join(params[cf.OUT_DIR], out_type + ".tif")
    md[ifc.EPOCH_DATE] = epochlist.dates
    if out_type == 'stack_rate':
        md[ifc.DATA_TYPE] = ifc.STACKRATE
    elif out_type == 'stack_error':
        md[ifc.DATA_TYPE] = ifc.STACKERROR
    else:
        md[ifc.DATA_TYPE] = ifc.STACKSAMP

    rate = np.zeros(shape=ifgs[0].shape, dtype=np.float32)

    for t in tiles:
        rate_file = os.path.join(params[cf.TMPDIR], out_type + '_'+str(t.index)+'.npy')
        rate_file = Path(rate_file)
        rate_tile = np.load(file=rate_file)
        rate[t.top_left_y:t.bottom_right_y, t.top_left_x:t.bottom_right_x] = rate_tile
    shared.write_output_geotiff(md, gt, wkt, rate, dest, np.nan)
    npy_rate_file = os.path.join(params[cf.OUT_DIR], out_type + '.npy')
    np.save(file=npy_rate_file, arr=rate)

    log.debug('Finished PyRate merging {}'.format(out_type))


def _merge_timeseries(rows, cols, params):
    """
    Merge time series output
    """
    xlks, _, crop = cf.transform_params(params)

    base_unw_paths = []

    for p in Path(params[OUT_DIR]).rglob("*rlks_*cr.tif"):
        if "dem" not in str(p):
            base_unw_paths.append(str(p))

    if "tif" in base_unw_paths[0].split(".")[1]:
        dest_tifs = base_unw_paths  # cf.get_dest_paths(base_unw_paths, crop, params, xlks)
        for i, dest_tif in enumerate(dest_tifs):
            dest_tifs[i] = dest_tif.replace("_tif", "")
    else:
        dest_tifs = base_unw_paths  # cf.get_dest_paths(base_unw_paths, crop, params, xlks)

    output_dir = params[cf.TMPDIR]
    # load previously saved prepread_ifgs dict
    preread_ifgs_file = join(output_dir, 'preread_ifgs.pk')
    ifgs = cp.load(open(preread_ifgs_file, 'rb'))

    # metadata and projections
    gt, md, wkt = ifgs['gt'], ifgs['md'], ifgs['wkt']
    epochlist = ifgs['epochlist']
    ifgs = [v for v in ifgs.values() if isinstance(v, PrereadIfg)]

    tiles = shared.get_tiles(dest_tifs[0], rows, cols)

    # load the first tsincr file to determine the number of time series tifs
    tsincr_file = os.path.join(output_dir, 'tsincr_0.npy')

    tsincr = np.load(file=tsincr_file)

    # pylint: disable=no-member
    no_ts_tifs = tsincr.shape[2]
    # we create 2 x no_ts_tifs as we are splitting tsincr and tscuml
    # to all processes.
    process_tifs = mpiops.array_split(range(2 * no_ts_tifs))

    # depending on nvelpar, this will not fit in memory
    # e.g. nvelpar=100, nrows=10000, ncols=10000, 32bit floats need 40GB memory
    # 32 * 100 * 10000 * 10000 / 8 bytes = 4e10 bytes = 40 GB
    # the double for loop helps us overcome the memory limit
    log.info('Process {} writing {} timeseries tifs of '
             'total {}'.format(mpiops.rank, len(process_tifs), no_ts_tifs * 2))
    for i in process_tifs:
        tscum_g = np.empty(shape=ifgs[0].shape, dtype=np.float32)
        if i < no_ts_tifs:
            for n, t in enumerate(tiles):
                _assemble_tiles(i, n, t, tscum_g, output_dir, 'tscuml')
            md[ifc.EPOCH_DATE] = epochlist.dates[i + 1]
            # sequence position; first time slice is #0
            md['SEQUENCE_POSITION'] = i+1
            dest = os.path.join(params[cf.OUT_DIR],
                                'tscuml' + "_" +
                                str(epochlist.dates[i + 1]) + ".tif")
            md[ifc.DATA_TYPE] = ifc.CUML
            shared.write_output_geotiff(md, gt, wkt, tscum_g, dest, np.nan)
        else:
            tsincr_g = np.empty(shape=ifgs[0].shape, dtype=np.float32)
            i %= no_ts_tifs
            for n, t in enumerate(tiles):
                _assemble_tiles(i, n, t, tsincr_g, output_dir, 'tsincr')
            md[ifc.EPOCH_DATE] = epochlist.dates[i + 1]
            # sequence position; first time slice is #0
            md['SEQUENCE_POSITION'] = i+1
            dest = os.path.join(params[cf.OUT_DIR],
                                'tsincr' + "_" + str(
                                    epochlist.dates[i + 1]) + ".tif")
            md[ifc.DATA_TYPE] = ifc.INCR
            shared.write_output_geotiff(md, gt, wkt, tsincr_g, dest, np.nan)
    mpiops.comm.barrier()
    log.debug('Process {} finished writing {} timeseries tifs of '
             'total {}'.format(mpiops.rank, len(process_tifs), no_ts_tifs * 2))


def _assemble_tiles(i, n, tile, tsincr_g, output_dir, outtype):
    # pylint: disable=too-many-arguments
    """
    A reusable time series tile assembly function
    """
    tsincr_file = os.path.join(output_dir, '{}_{}.npy'.format(outtype, n))
    tsincr = np.load(file=tsincr_file)
    tsincr_g[tile.top_left_y:tile.bottom_right_y, tile.top_left_x:tile.bottom_right_x] = tsincr[:, :, i]


def _delete_tsincr_files(params):
    """
    Convenience function to delete tsincr files
    """
    out_dir = Path(params[cf.OUT_DIR])
    for file_path in out_dir.iterdir():
        if "tsincr" in str(file_path):
            file_path.unlink()
