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
This Python script applies optional multilooking and cropping to input
interferogram geotiff files.
"""
# -*- coding: utf-8 -*-
import os
from subprocess import check_call
from typing import List
from pathlib import Path
from joblib import Parallel, delayed
import numpy as np
from pyrate.core import shared, mpiops, config as cf, prepifg_helper, gamma, roipac
from pyrate.core.prepifg_helper import PreprocessError
from pyrate.core.logger import pyratelogger as log
from pyrate.core.shared import output_tiff_filename


GAMMA = 1
ROIPAC = 0
GEOTIF = 2


def main(params):
    """
    Main workflow function for preparing interferograms for PyRate.

    :param dict params: Parameters dictionary read in from the config file
    """
    # TODO: looks like ifg_paths are ordered according to ifg list
    # This probably won't be a problem because input list won't be reordered
    # and the original gamma generated list is ordered) this may not affect
    # the important pyrate stuff anyway, but might affect gen_thumbs.py.
    # Going to assume ifg_paths is ordered correcly
    # pylint: disable=too-many-branches
    shared.mpi_vs_multiprocess_logging("prepifg", params)

    ifg_paths = params[cf.INTERFEROGRAM_FILES]
    if params[cf.DEM_FILE] is not None:  # optional DEM conversion
        ifg_paths.append(params[cf.DEM_FILE_PATH])

    shared.mkdir_p(params[cf.OUT_DIR])  # create output dir

    process_base_ifgs_paths = np.array_split(ifg_paths, mpiops.size)[mpiops.rank]
    gtiff_paths = [p.converted_path for p in process_base_ifgs_paths]
    do_prepifg(gtiff_paths, params)
    mpiops.comm.barrier()
    log.info("Finished prepifg")


def do_prepifg(gtiff_paths: List[str], params: dict) -> None:
    """
    Prepare interferograms by applying multilooking/cropping operations.

    :param list gtiff_paths: List of full-res geotiffs
    :param dict params: Parameters dictionary corresponding to config file
    """
    # pylint: disable=expression-not-assigned
    parallel = params[cf.PARALLEL]

    for f in gtiff_paths:
        if not os.path.isfile(f):
            raise FileNotFoundError("Can not find geotiff: " + str(f) + ". Ensure you have converted your "
                                    "interferograms to geotiffs.")

    ifgs = [prepifg_helper.dem_or_ifg(p) for p in gtiff_paths]
    xlooks, ylooks, crop = cf.transform_params(params)
    user_exts = (params[cf.IFG_XFIRST], params[cf.IFG_YFIRST], params[cf.IFG_XLAST], params[cf.IFG_YLAST])
    exts = prepifg_helper.get_analysis_extent(crop, ifgs, xlooks, ylooks, user_exts=user_exts)
    thresh = params[cf.NO_DATA_AVERAGING_THRESHOLD]

    if params[cf.LARGE_TIFS]:
        ifg = ifgs[0]
        res_str = [xlooks * ifg.x_step, ylooks * ifg.y_step]
        res_str = ' '.join([str(e) for e in res_str])
        # if parallel:
        #     Parallel(n_jobs=params[cf.PROCESSES], verbose=50)(
        #         delayed(__prepifg_multiprocess_system)(crop, exts, gtiff_paths, res_str, params, thresh, xlooks, ylooks)
        #     )
        # else:
        for gtiff_path in gtiff_paths:
            __prepifg_multiprocess_system(crop, exts, gtiff_path, params, res_str, thresh, xlooks, ylooks)
    else:
        if parallel:
            Parallel(n_jobs=params[cf.PROCESSES], verbose=50)(
                delayed(_prepifg_multiprocessing)(p, xlooks, ylooks, exts, thresh, crop, params) for p in gtiff_paths
            )
        else:
            for gtiff_path in gtiff_paths:
                _prepifg_multiprocessing(gtiff_path, xlooks, ylooks, exts, thresh, crop, params)


def __prepare_prepifg_files(crop, exts, gtiff_paths, res, params, thresh, xlooks, ylooks):
    # calculate nan fraction matrix, crop resample average nan fraction
    for gtiff_path in gtiff_paths:
        __prepifg_multiprocess_system(crop, exts, gtiff_path, params, res, thresh, xlooks, ylooks)


def __prepifg_multiprocess_system(crop, exts, gtiff_path, params, res, thresh, xlooks, ylooks):
    p, c, l = _prepifg_multiprocessing(gtiff_path, xlooks, ylooks, exts, thresh, crop, params)
    extents = ' '.join([str(e) for e in exts])

    # change nodataval from zero
    p_unset = Path(p).with_suffix('.unset.tif')
    check_call('gdal_translate -a_nodata -99999\t'
               '{p} {out_file}'.format(p=p, out_file=p_unset), shell=True)

    # calculate nan-fraction
    out_file = Path(l).with_suffix('.nanfrac.tif')
    check_call('gdal_calc.py --overwrite -A {p}\t'
               '--calc=\"isclose(A, 0, 0.000001)\"\t'
               '--outfile={out_file}\t'
               '--NoDataValue=-99999\n'.format(p=p_unset, out_file=out_file), shell=True)
    out_file_avg = Path(l).with_suffix('.nanfrac.avg.tif')
    # crop resample/average multilooking of nan-fraction
    check_call('gdalwarp -te\t{extents}\t-tr\t{res}\t-r\taverage\t{p}\t{out_file}\n'.format(
        extents=extents, res=res, p=out_file, out_file=out_file_avg), shell=True)
    if c is not None:
        # coh masking
        # change no data value
        c_out = Path(c).with_suffix('.coh.tif')
        check_call('gdal_translate -a_nodata -99999\t'
                    '{c} {out_file}'.format(c=c, out_file=c_out), shell=True)

        coh_corrected_p = Path(p).with_suffix('.coh.corrected.tif')
        check_call('gdal_calc.py \t-A\t{c}\t-B\t{p}\t--outfile={out_file}\t'
                   '--calc=\"B*(A>={th})-99999*(A<{th})\"\t'
                   '--NoDataValue=-99999'.format(c=c_out, p=p, th=params[cf.COH_THRESH], out_file=coh_corrected_p),
                   shell=True)
        in_file = coh_corrected_p
    else:
        in_file = p
    # crop resample/average multilooking of raster\
    temp_l = Path(l).with_suffix('.temp.tif')
    check_call('gdalwarp -te\t{extents}\t-tr\t{res}\t-r\taverage\t{p}\t{l}\n'.format(
        extents=extents, res=res, p=in_file, l=temp_l), shell=True)
    # resampled_average[nan_frac >= thresh] = nodatavalue
    check_call('gdal_calc.py -A {p}\t-B {q}\t'
               '--calc=\"B*(A<{th})-99999*(A>={th})\"\t'
               '--outfile={out_file}\t'
               '--NoDataValue=-99999\n'.format(p=out_file_avg, q=temp_l, out_file=l, th=thresh), shell=True)


def _prepifg_multiprocessing(path, xlooks, ylooks, exts, thresh, crop, params):
    """
    Multiprocessing wrapper for prepifg
    """
    processor = params[cf.PROCESSOR]  # roipac, gamma or geotif
    if (processor == GAMMA) or (processor == GEOTIF):
        header = gamma.gamma_header(path, params)
    elif processor == ROIPAC:
        header = roipac.roipac_header(path, params)
    else:
        raise PreprocessError('Processor must be ROI_PAC (0) or GAMMA (1)')

    # If we're performing coherence masking, find the coherence file for this IFG.
    if params[cf.COH_MASK] and shared._is_interferogram(header):
        coherence_path = cf.coherence_paths_for(path, params, tif=True)
        coherence_thresh = params[cf.COH_THRESH]
    else:
        coherence_path = None
        coherence_thresh = None

    if params[cf.LARGE_TIFS]:
        op = output_tiff_filename(path, params[cf.OUT_DIR])
        looks_path = cf.mlooked_path(op, ylooks, crop)
        return path, coherence_path, looks_path
    else:
        prepifg_helper.prepare_ifg(path, xlooks, ylooks, exts, thresh, crop, out_path=params[cf.OUT_DIR],
                                   header=header, coherence_path=coherence_path, coherence_thresh=coherence_thresh)
