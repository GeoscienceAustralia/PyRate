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
This Python script applies optional multilooking and cropping to input
interferogram geotiff files.

There are two modes of running prepifg using pyrate:
1. A python/numpy version, which is also the default version, and
2. a `largetifs` option which can be activated using the config option largetifs: 1

The python/numpy version is recommended when the both the input interferogram and multilooked interferogram can fit
into the memory allocated to the process.

When dealing with relatively large (compared to memory available to the process) interferogram, the largetif option
can be used. This option uses a slightly modified version of the `gdal_calc.py` included with pyrate. Arbitrarily
large interferograms can be multilooked using largetifs.

The `largetifs` option uses system calls to gdal utilities and avoid loading large chunks of memory. Our tests
(tests/test_prepifg_largetifs_vs_python.py::test_prepifg_largetifs_vs_python) indicate that for two different small
datasets included in our test suite, the `largetifs` option exactly match the output of the numpy version. However,
on large practical datasets we have observed numerical differences in the multilooked output in the 3rd and 4th decimal
places (in sub mm range).

"""
# -*- coding: utf-8 -*-
import os
from subprocess import check_call
import warnings
from typing import List, Tuple, Dict
from pathlib import Path
from joblib import Parallel, delayed
import numpy as np
from osgeo import gdal

import pyrate.constants as C
from pyrate.core import shared, geometry, mpiops, prepifg_helper, gamma, roipac, ifgconstants as ifc
from pyrate.core.prepifg_helper import PreprocessError, coherence_paths_for, transform_params
from pyrate.core.logger import pyratelogger as log
from pyrate.configuration import MultiplePaths, Configuration
from pyrate.core.shared import Ifg, DEM
from pyrate.core.refpixel import convert_geographic_coordinate_to_pixel_value
from pyrate.core import ifgconstants as ifg


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

    ifg_paths = params[C.INTERFEROGRAM_FILES]
    if params[C.DEM_FILE] is not None:  # optional DEM conversion
        ifg_paths.append(params[C.DEM_FILE_PATH])

    if params[C.COH_FILE_LIST] is not None:
        ifg_paths.extend(params[C.COHERENCE_FILE_PATHS])

    if params[C.COH_FILE_LIST] is None and params[C.COH_MASK]:
        raise FileNotFoundError("Cannot apply coherence masking: no coherence file list "
                                "supplied (parameter 'cohfilelist')")

    shared.mkdir_p(params[C.OUT_DIR])  # create output dir

    user_exts = (params[C.IFG_XFIRST], params[C.IFG_YFIRST], params[
        C.IFG_XLAST], params[C.IFG_YLAST])
    xlooks, ylooks, crop = transform_params(params)
    ifgs = [prepifg_helper.dem_or_ifg(p.converted_path) for p in ifg_paths]
    exts = prepifg_helper.get_analysis_extent(crop, ifgs, xlooks, ylooks, user_exts=user_exts)

    ifg0 = ifgs[0]
    ifg0.open()

    transform = ifg0.dataset.GetGeoTransform()

    process_ifgs_paths = np.array_split(ifg_paths, mpiops.size)[mpiops.rank]
    do_prepifg(process_ifgs_paths, exts, params)

    if params[C.LT_FILE] is not None:
        log.info("Calculating and writing geometry files")
        __write_geometry_files(params, exts, transform, ifg_paths[0].sampled_path)
    else:
        log.info("Skipping geometry calculations: Lookup table not provided")

    if params[C.COH_FILE_LIST] is not None:
        log.info("Calculating and writing coherence statistics")
        mpiops.run_once(__calc_coherence_stats, params, ifg_paths[0].sampled_path)
    else:
        log.info("Skipping coherence file statistics computation.")
    log.info("Finished 'prepifg' step")
    ifg0.close()


def __calc_coherence_stats(params, ifg_path):
    coherence_files_multi_paths = params[C.COHERENCE_FILE_PATHS]
    sampled_paths = [c.sampled_path for c in coherence_files_multi_paths]
    ifgs = [Ifg(s) for s in sampled_paths]
    for i in ifgs:
        i.open()
    phase_data = np.stack([i.phase_data for i in ifgs])
    coh_stats = Configuration.coherence_stats(params)

    for stat_func, out_type in zip([np.nanmedian, np.nanmean, np.nanstd], [ifg.COH_MEDIAN, ifg.COH_MEAN, ifg.COH_STD]):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            arr = stat_func(phase_data, axis=0)
        arr[arr==0.0] = np.nan # convert exact zeros (no-data) to NaN
        dest = coh_stats[out_type]
        __save_geom_files(ifg_path, dest, arr, out_type)


def do_prepifg(multi_paths: List[MultiplePaths], exts: Tuple[float, float, float, float], params: dict) -> None:
    """
    Prepare interferograms by applying multilooking/cropping operations.

    """
    # pylint: disable=expression-not-assigned
    parallel = params[C.PARALLEL]
    if mpiops.size > 1:
        parallel = False

    for f in multi_paths:
        if not os.path.isfile(f.converted_path):
            raise FileNotFoundError("Can not find geotiff: " + str(f) + ". Ensure you have converted your "
                                    "interferograms to geotiffs.")

    if params[C.LARGE_TIFS]:
        log.info("Using gdal system calls to execute 'prepifg' step")
        ifg = prepifg_helper.dem_or_ifg(multi_paths[0].converted_path)
        ifg.open()
        xlooks, ylooks = params[C.IFG_LKSX], params[C.IFG_LKSY]
        res_str = [xlooks * ifg.x_step, ylooks * ifg.y_step]
        res_str = ' '.join([str(e) for e in res_str])
        if parallel:
            Parallel(n_jobs=params[C.PROCESSES], verbose=50)(
                delayed(__prepifg_system)(exts, gtiff_path, params, res_str) for gtiff_path in multi_paths)
        else:
            for m_path in multi_paths:
                __prepifg_system(exts, m_path, params, res_str)
    else:
        if parallel:
            Parallel(n_jobs=params[C.PROCESSES], verbose=50)(
                delayed(_prepifg_multiprocessing)(p, exts, params) for p in multi_paths
            )
        else:
            for m_path in multi_paths:
                _prepifg_multiprocessing(m_path, exts, params)
    mpiops.comm.barrier()


COMMON_OPTIONS = "-co BLOCKXSIZE=256 -co BLOCKYSIZE=256 -co TILED=YES --config GDAL_CACHEMAX=64 -q"
COMMON_OPTIONS2 = "--co BLOCKXSIZE=256 --co BLOCKYSIZE=256 --co TILED=YES --quiet"
GDAL_CALC = 'gdal_calc_local.py'


def __prepifg_system(exts, gtiff, params, res):
    thresh = params[C.NO_DATA_AVERAGING_THRESHOLD]
    p, c, l = _prepifg_multiprocessing(gtiff, exts, params)
    log.info("Multilooking {p} into {l}".format(p=p, l=l))
    extents = ' '.join([str(e) for e in exts])

    if isinstance(prepifg_helper.dem_or_ifg(p), shared.DEM):
        check_call('gdalwarp {co} -te\t{extents}\t-tr\t{res}\t-r\taverage \t{p}\t{l}\n'.format(
            co=COMMON_OPTIONS, extents=extents, res=res, p=p, l=l), shell=True)
        __update_meta_data(p, c, l, params)
        return

    p_unset = Path(params[C.OUT_DIR]).joinpath(Path(p).name).with_suffix('.unset.tif')
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
                   f'--calc=\"logical_or((B<{params[C.COH_THRESH]}), isclose(A,0,atol=0.000001))\"\t'
                   f'--NoDataValue=nan', shell=True)

        # coh masking
        check_call(f'{GDAL_CALC} {COMMON_OPTIONS2} --overwrite -A {p_unset} -B {c}\t'
                   f'--calc=\"A*(B>={params[C.COH_THRESH]}) - '
                   f'99999*logical_or((B<{params[C.COH_THRESH]}), isclose(A,0,atol=0.000001))\"\t'
                   f'--outfile={corrected_p}\t'
                   f'--NoDataValue=nan', shell=True)
    else:
        log.info(f"applying nodata masking on {p}")
        check_call(f'{GDAL_CALC} {COMMON_OPTIONS2} --overwrite -A {p_unset}\t'
                   f'--calc=\"isclose(A,0,atol=0.000001)\"\t'
                   f'--outfile={nan_frac}\t'
                   f'--NoDataValue=nan', shell=True)
        check_call(f'{GDAL_CALC} {COMMON_OPTIONS2} --overwrite -A {p_unset}\t'
                   f'--calc=\"A - 99999*isclose(A, 0, atol=0.000001)\"\t'
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

    __update_meta_data(p_unset.as_posix(), c, l, params)

    # clean up
    nan_frac_avg.unlink()
    nan_frac.unlink()
    corrected_p.unlink()
    p_unset.unlink()


def __update_meta_data(p_unset, c, l, params):
    # update metadata
    ds = gdal.Open(p_unset)
    md = ds.GetMetadata()
    # remove data type
    v = md.pop(ifc.DATA_TYPE)
    md[ifc.IFG_LKSX] = str(params[C.IFG_LKSX])
    md[ifc.IFG_LKSY] = str(params[C.IFG_LKSY])
    md[ifc.IFG_CROP] = str(params[C.IFG_CROP_OPT])
    # update data type
    if c is not None:  # it's a interferogram when COH_MASK=1
        md_str = '-mo {k}={v}'.format(k=ifc.DATA_TYPE, v=ifc.MLOOKED_COH_MASKED_IFG)
        md[ifc.COH_THRESH] = str(params[C.COH_THRESH])
    else:
        if v == ifc.DEM:  # it's a dem
            md_str = '-mo {k}={v}'.format(k=ifc.DATA_TYPE, v=ifc.MLOOKED_DEM)
        elif v == ifc.COH:
            md_str = '-mo {k}={v}'.format(k=ifc.DATA_TYPE, v=ifc.MULTILOOKED_COH)
        else:  # it's an ifg
            md_str = '-mo {k}={v}'.format(k=ifc.DATA_TYPE, v=ifc.MULTILOOKED)
    for k, v in md.items():
        md_str += ' -mo {k}={v}'.format(k=k.replace(' ', '_'), v=v.replace(' ', '_'))
    check_call('gdal_edit.py -unsetmd {md} {f}'.format(md=md_str, f=l), shell=True)
    ds = None

    # make prepifg output readonly
    Path(l).chmod(0o444)  # readonly output


def _prepifg_multiprocessing(m_path: MultiplePaths, exts: Tuple[float, float, float, float], params: dict):
    """
    Multiprocessing wrapper for prepifg
    """
    thresh = params[C.NO_DATA_AVERAGING_THRESHOLD]
    hdr = find_header(m_path, params)
    hdr[ifc.INPUT_TYPE] = m_path.input_type
    xlooks, ylooks, crop = transform_params(params)
    hdr[ifc.IFG_LKSX] = xlooks
    hdr[ifc.IFG_LKSY] = ylooks
    hdr[ifc.IFG_CROP] = crop

    # If we're performing coherence masking, find the coherence file for this IFG.
    if params[C.COH_MASK] and shared._is_interferogram(hdr):
        coherence_path = coherence_paths_for(m_path.converted_path, params, tif=True)
        coherence_thresh = params[C.COH_THRESH]
    else:
        coherence_path = None
        coherence_thresh = None

    if params[C.LARGE_TIFS]:
        return m_path.converted_path, coherence_path, m_path.sampled_path
    else:
        prepifg_helper.prepare_ifg(m_path.converted_path, xlooks, ylooks, exts, thresh, crop,
                                   out_path=m_path.sampled_path, header=hdr, coherence_path=coherence_path,
                                   coherence_thresh=coherence_thresh)
        Path(m_path.sampled_path).chmod(0o444)  # readonly output


def find_header(path: MultiplePaths, params: dict) -> Dict[str, str]:
    processor = params[C.PROCESSOR]  # roipac, gamma or geotif
    tif_path = path.converted_path
    if (processor == GAMMA) or (processor == GEOTIF):
        header = gamma.gamma_header(tif_path, params)
    elif processor == ROIPAC:
        import warnings
        warnings.warn("Warning: ROI_PAC support will be deprecated in a future PyRate release",
                      category=DeprecationWarning)
        header = roipac.roipac_header(tif_path, params)
    else:
        raise PreprocessError('Processor must be ROI_PAC (0) or GAMMA (1)')
    header[ifc.INPUT_TYPE] = path.input_type
    return header


def __write_geometry_files(params: dict, exts: Tuple[float, float, float, float],
                           transform, ifg_path: str) -> None:
    """
    Calculate geometry and save to geotiff files using the information in the
    first interferogram in the stack, i.e.:
    - rdc_azimuth.tif (azimuth radar coordinate at each pixel)
    - rdc_range.tif (range radar coordinate at each pixel)
    - azimuth_angle.tif (satellite azimuth angle at each pixel)
    - incidence_angle.tif (incidence angle at each pixel)
    - look_angle.tif (look angle at each pixel)
    """
    ifg = Ifg(ifg_path)
    ifg.open(readonly=True)

    # calculate per-pixel lon/lat
    lon, lat = geometry.get_lonlat_coords(ifg)

    # not currently implemented for ROIPAC data which breaks some tests
    # if statement can be deleted once ROIPAC is deprecated from PyRate
    if ifg.meta_data[ifc.PYRATE_INSAR_PROCESSOR] == 'ROIPAC':
        log.warning("Geometry calculations are not implemented for ROI_PAC")
        return

    # get geometry information and save radar coordinates and angles to tif files
    # using metadata of the first image in the stack
    # get pixel values of crop (needed to crop lookup table file)
    # pixel extent of cropped area (original IFG input)
    xmin, ymax = convert_geographic_coordinate_to_pixel_value(exts[0], exts[1], transform)
    xmax, ymin = convert_geographic_coordinate_to_pixel_value(exts[2], exts[3], transform)
    # xmin, xmax: columns of crop
    # ymin, ymax: rows of crop

    # calculate per-pixel radar coordinates
    az, rg = geometry.calc_radar_coords(ifg, params, xmin, xmax, ymin, ymax)

    # Read height data from DEM
    dem_file = params[C.DEM_FILE_PATH].sampled_path
    dem = DEM(dem_file)
    # calculate per-pixel look angle (also calculates and saves incidence and azimuth angles)
    lk_ang, inc_ang, az_ang, rg_dist = geometry.calc_pixel_geometry(ifg, rg, lon.data, lat.data, dem.data)

    # save radar coordinates and angles to geotiff files
    combinations = zip([az, rg, lk_ang, inc_ang, az_ang, rg_dist], C.GEOMETRY_OUTPUT_TYPES)
    shared.iterable_split(__parallelly_write, combinations, params, ifg_path)


def __parallelly_write(tup, params, ifg_path):
    array, output_type = tup
    dest = Configuration.geometry_files(params)[output_type]
    if mpiops.size > 0:
        log.debug(f"Writing {dest} using process {mpiops.rank}")
    __save_geom_files(ifg_path, dest, array, output_type)


out_type_md_dict = {
    'rdc_azimuth': ifc.RDC_AZIMUTH,
    'rdc_range': ifc.RDC_RANGE,
    'look_angle': ifc.LOOK,
    'incidence_angle': ifc.INCIDENCE,
    'azimuth_angle': ifc.AZIMUTH,
    'range_dist': ifc.RANGE_DIST,
    ifc.COH_MEAN: ifc.COH_MEAN,
    ifc.COH_MEDIAN: ifc.COH_MEDIAN,
    ifc.COH_STD: ifc.COH_STD + '_1SIGMA'
}


def __save_geom_files(ifg_path, dest, array, out_type):
    """
    Convenience function to save geometry geotiff files
    """
    log.debug('Saving PyRate outputs {}'.format(out_type))
    gt, md, wkt = shared.get_geotiff_header_info(ifg_path)
    md[ifc.DATA_TYPE] = out_type_md_dict[out_type]
    shared.remove_file_if_exists(dest)
    log.info(f"Writing geotiff: {dest}")
    shared.write_output_geotiff(md, gt, wkt, array, dest, np.nan)
