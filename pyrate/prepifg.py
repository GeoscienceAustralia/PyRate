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
from osgeo import gdal
from pyrate.core import shared, mpiops, config as cf, prepifg_helper, gamma, roipac, ifgconstants as ifc
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

    process_ifgs_paths = np.array_split(ifg_paths, mpiops.size)[mpiops.rank]

    gtiff_paths = [p.converted_path for p in process_ifgs_paths]
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
    if mpiops.size > 1:
        parallel = False

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
        log.info("Using gdal system calls to process prepifg")
        ifg = ifgs[0]
        res_str = [xlooks * ifg.x_step, ylooks * ifg.y_step]
        res_str = ' '.join([str(e) for e in res_str])
        if parallel:
            Parallel(n_jobs=params[cf.PROCESSES], verbose=50)(
                delayed(__prepifg_system)(
                    crop, exts, gtiff_path, params, res_str, thresh, xlooks, ylooks) for gtiff_path in gtiff_paths
            )
        else:
            for gtiff_path in gtiff_paths:
                __prepifg_system(crop, exts, gtiff_path, params, res_str, thresh, xlooks, ylooks)
    else:
        if parallel:
            Parallel(n_jobs=params[cf.PROCESSES], verbose=50)(
                delayed(_prepifg_multiprocessing)(p, xlooks, ylooks, exts, thresh, crop, params) for p in gtiff_paths
            )
        else:
            for gtiff_path in gtiff_paths:
                _prepifg_multiprocessing(gtiff_path, xlooks, ylooks, exts, thresh, crop, params)


COMMON_OPTIONS = "-co BLOCKXSIZE=256 -co BLOCKYSIZE=256 -co TILED=YES --config GDAL_CACHEMAX=64 -q"
COMMON_OPTIONS2 = "--co BLOCKXSIZE=256 --co BLOCKYSIZE=256 --co TILED=YES --quiet"
GDAL_CALC = 'gdal_calc_local.py'


def __prepifg_system(crop, exts, gtiff, params, res, thresh, xlooks, ylooks):
    p, c, l = _prepifg_multiprocessing(gtiff, xlooks, ylooks, exts, thresh, crop, params)
    log.info("Multilooking {p} into {l}".format(p=p, l=l))
    extents = ' '.join([str(e) for e in exts])

    if isinstance(prepifg_helper.dem_or_ifg(p), shared.DEM):
        check_call('gdalwarp {co} -te\t{extents}\t-tr\t{res}\t-r\taverage \t{p}\t{l}\n'.format(
            co=COMMON_OPTIONS, extents=extents, res=res, p=p, l=l), shell=True)
        __update_meta_data(p, c, l)
        return

    p_unset = Path(params[cf.OUT_DIR]).joinpath(Path(p).name).with_suffix('.unset.tif')
    # change nodataval from zero, also leave input geotifs unchanged if one supplies conv2tif output/geotifs
    check_call('gdal_translate {co} -a_nodata nan\t{p}\t{q}'.format(co=COMMON_OPTIONS, p=p, q=p_unset),
               shell=True)

    # calculate nan-fraction
    # TODO: use output options and datatypes to reduce size of the next two tifs
    nan_frac = Path(l).with_suffix('.nanfrac.tif')
    nan_frac_avg = Path(l).with_suffix('.nanfrac.avg.tif')
    corrected_p = Path(p_unset).with_suffix('.corrected.tif')

    if c is not None:
        # find all the nans
        log.info(f"applying coherence + nodata masking on {p}")
        check_call(f'{GDAL_CALC} {COMMON_OPTIONS2} -A {p_unset} -B {c} --outfile={nan_frac}\t'
                   f'--calc=\"logical_or((B<{params[cf.COH_THRESH]}), isclose(A,0,atol=0.000001))\"\t'
                   f'--NoDataValue=nan', shell=True)

        # coh masking
        check_call(f'{GDAL_CALC} {COMMON_OPTIONS2} --overwrite -A {p_unset} -B {c}\t'
                   f'--calc=\"A*(B>={params[cf.COH_THRESH]})'
                   f'-99999*logical_or((B<{params[cf.COH_THRESH]}),isclose(A,0,atol=0.000001))\"\t'
                   f'--outfile={corrected_p}\t'
                   f'--NoDataValue=nan', shell=True)
    else:
        log.info(f"applying nodata masking on {p}")
        check_call(f'{GDAL_CALC} {COMMON_OPTIONS2} --overwrite -A {p_unset}\t'
                   f'--calc=\"isclose(A,0,atol=0.000001)\"\t'
                   f'--outfile={nan_frac}\t'
                   f'--NoDataValue=nan', shell=True)
        check_call(f'{GDAL_CALC} {COMMON_OPTIONS2} --overwrite -A {p_unset}\t'
                   f'--calc=\"A - 99999 *isclose(A, 0, atol=0.000001)\"\t'
                   f'--outfile={corrected_p}\t'
                   f'--NoDataValue=nan', shell=True)

    # crop resample/average multilooking of nan-fraction
    check_call('gdalwarp {co} -te\t{extents}\t-tr\t{res}\t-r\taverage\t{p}\t{out_file}'.format(
        co=COMMON_OPTIONS, extents=extents, res=res, p=nan_frac, out_file=nan_frac_avg), shell=True)

    # crop resample/average multilooking of raster
    check_call('gdalwarp {co} -te\t{extents}\t-tr\t{res}\t-r\taverage \t{p}\t{l}'.format(
        co=COMMON_OPTIONS, extents=extents, res=res, p=corrected_p, l=l), shell=True)

    check_call(f'{GDAL_CALC} {COMMON_OPTIONS2} --overwrite -A {nan_frac_avg}\t-B {l}\t'
               f'--calc=\"B*(A < {thresh}) -99999*(A >= {thresh})\"\t'
               f'--outfile={l}\t'
               f'--NoDataValue=nan', shell=True)

    __update_meta_data(p_unset.as_posix(), c, l)

    # clean up
    nan_frac_avg.unlink()
    nan_frac.unlink()
    corrected_p.unlink()
    p_unset.unlink()


def __update_meta_data(p_unset, c, l):
    # update metadata
    ds = gdal.Open(p_unset)
    md = ds.GetMetadata()
    # remove data type
    v = md.pop(ifc.DATA_TYPE)
    # update data type
    if c is not None:  # it's a interferogram when COH_MASK=1
        md_str = '-mo {k}={v}'.format(k=ifc.DATA_TYPE, v=ifc.COHERENCE)
    else:
        if v == ifc.DEM:  # it's a dem
            md_str = '-mo {k}={v}'.format(k=ifc.DATA_TYPE, v=ifc.MLOOKED_DEM)
        else:  # it's an ifg
            md_str = '-mo {k}={v}'.format(k=ifc.DATA_TYPE, v=ifc.MULTILOOKED)
    for k, v in md.items():
        md_str += ' -mo {k}={v}'.format(k=k.replace(' ', '_'), v=v.replace(' ', '_'))
    check_call('gdal_edit.py -unsetmd {md} {f}'.format(md=md_str, f=l), shell=True)
    ds = None


def _prepifg_multiprocessing(path, xlooks, ylooks, exts, thresh, crop, params):
    """
    Multiprocessing wrapper for prepifg
    """
    processor = params[cf.PROCESSOR]  # roipac, gamma or geotif
    if (processor == GAMMA) or (processor == GEOTIF):
        header = gamma.gamma_header(path, params)
    elif processor == ROIPAC:
        log.info("Warning: ROI_PAC support will be deprecated in a future PyRate release")
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
