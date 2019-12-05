#   This Python module is part of the PyRate software package.
#
#   Copyright 2017 Geoscience Australia
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
This Python module does post-processing steps to assemble the linear
rate and time series outputs and save as geotiff files
"""
import os
from os.path import join
import logging
import pickle as cp
import numpy as np
import gdal
import subprocess
import pathlib
import time
from constants import REF_COLOR_MAP_PATH
from core import shared, ifgconstants as ifc, mpiops, config as cf
from core.shared import PrereadIfg
gdal.SetCacheMax(64)
log = logging.getLogger(__name__)

# Constants
MASTER_PROCESS = 0

def create_png_from_tif(output_folder_path):

    # open raster and choose band to find min, max
    raster_path = os.path.join(output_folder_path, "linrate.tif")

    if not os.path.isfile(raster_path):
        raise Exception("linrate.tif file not found at: "+raster_path)
    gtif = gdal.Open(raster_path)
    srcband = gtif.GetRasterBand(1)

    west, north, east, south = "", "", "", ""
    for line in gdal.Info(gtif).split('\n'):
        if "Upper Left" in line:
            west, north = line.split(")")[0].split("(")[1].split(",")
        if "Lower Right" in line:
            east, south = line.split(")")[0].split("(")[1].split(",")

    kml_file_path = os.path.join(output_folder_path, "linrate.kml")
    kml_file_content = """<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://earth.google.com/kml/2.1">
  <Document>
    <name>linrate.kml</name>
    <GroundOverlay>
      <name>linrate.png</name>
      <Icon>
        <href>linrate.png</href>
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

    # read color map from utilises and write it to the output folder

    with open(REF_COLOR_MAP_PATH, "r") as f:
        color_map_list = []
        for line in f.readlines():
            color_map_list.append(line.strip().split(" "))

    no_of_data_value = len(np.arange(minimum, maximum, step))
    for i, no in enumerate(np.arange(minimum, maximum, step)):
        color_map_list[i+1][0] = str(no)

    color_map_path = os.path.join(output_folder_path, "colormap.txt")
    with open(color_map_path, "w") as f:
        for i in range(no_of_data_value):
            f.write(' '.join(color_map_list[i]) + "\n")

    input_tif_path = os.path.join(output_folder_path, "linrate.tif")
    output_png_path = os.path.join(output_folder_path, "linrate.png")
    subprocess.check_call(["gdaldem", "color-relief", "-of", "PNG", input_tif_path, "-alpha", color_map_path, output_png_path, "-nearest_color_entry"])


def main(params, rows, cols):
    """
    PyRate merge main function. Assembles product tiles in to
    single geotiff files
    """
    # setup paths
    _merge_linrate(rows, cols, params)
    if params[cf.TIME_SERIES_CAL]:
        _merge_timeseries(rows, cols, params)

    log.info('Start creating quicklook results.')
    output_folder_path = os.path.dirname(params["tmpdir"])
    create_png_from_tif(output_folder_path)
    log.info('Finished creating quick look results.')

from core.config import OBS_DIR
def _merge_linrate(rows, cols, params):
    """
    Merge linear rate outputs
    """
    # setup paths
    xlks, _, crop = cf.transform_params(params)
    base_unw_paths = cf.original_ifg_paths(params[cf.IFG_FILE_LIST], params[cf.OBS_DIR])

    base_unw_paths = []
    for p in pathlib.Path(params[OBS_DIR]).rglob("*rlks_*cr.tif"):
        if "dem" not in str(p):
            base_unw_paths.append(str(p))

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

    # linrate aggregation
    if mpiops.size >= 3:
        [_save_linrate(ifgs, params, tiles, out_type=t)
         for i, t in enumerate(['linrate', 'linerror', 'linsamples'])
         if i == mpiops.rank]
    else:
        if mpiops.rank == MASTER_PROCESS:
            [_save_linrate(ifgs, params, tiles, out_type=t)
             for t in ['linrate', 'linerror', 'linsamples']]


def _save_linrate(ifgs_dict, params, tiles, out_type):
    """
    Save linear rate outputs
    """
    log.info('Merging PyRate outputs {}'.format(out_type))
    gt, md, wkt = ifgs_dict['gt'], ifgs_dict['md'], ifgs_dict['wkt']
    epochlist = ifgs_dict['epochlist']
    ifgs = [v for v in ifgs_dict.values() if isinstance(v, PrereadIfg)]
    dest = os.path.join(params[cf.OUT_DIR], out_type + ".tif")
    md[ifc.EPOCH_DATE] = epochlist.dates
    if out_type == 'linrate':
        md[ifc.DATA_TYPE] = ifc.LINRATE
    elif out_type == 'linerror':
        md[ifc.DATA_TYPE] = ifc.LINERROR
    else:
        md[ifc.DATA_TYPE] = ifc.LINSAMP

    rate = np.zeros(shape=ifgs[0].shape, dtype=np.float32)

    for t in tiles:
        rate_file = os.path.join(params[cf.TMPDIR], out_type + '_'+str(t.index)+'.npy')
        rate_file = pathlib.Path(rate_file)
        rate_tile = np.load(file=rate_file)
        rate[t.top_left_y:t.bottom_right_y, t.top_left_x:t.bottom_right_x] = rate_tile
    shared.write_output_geotiff(md, gt, wkt, rate, dest, np.nan)
    npy_rate_file = os.path.join(params[cf.OUT_DIR], out_type + '.npy')
    np.save(file=npy_rate_file, arr=rate)

    log.info('Finished PyRate merging {}'.format(out_type))


def _merge_timeseries(rows, cols, params):
    """
    Merge time series output
    """
    xlks, _, crop = cf.transform_params(params)
    base_unw_paths = cf.original_ifg_paths(params[cf.IFG_FILE_LIST], params[cf.OBS_DIR])

    base_unw_paths = []
    for p in pathlib.Path(params[OBS_DIR]).rglob("*rlks_*cr.tif"):
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
    log.info('process {} will write {} ts (incr/cuml) tifs of '
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
    log.info('process {} finished writing {} ts (incr/cuml) tifs of '
             'total {}'.format(mpiops.rank, len(process_tifs), no_ts_tifs * 2))


def _assemble_tiles(i, n, tile, tsincr_g, output_dir, outtype):
    # pylint: disable=too-many-arguments
    """
    A reusable time series tile assembly function
    """
    tsincr_file = os.path.join(output_dir, '{}_{}.npy'.format(outtype, n))
    tsincr = np.load(file=tsincr_file)
    tsincr_g[tile.top_left_y:tile.bottom_right_y, tile.top_left_x:tile.bottom_right_x] = tsincr[:, :, i]
