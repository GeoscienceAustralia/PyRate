from __future__ import print_function
"""
Main workflow script for PyRate
"""
import logging
import os
from os.path import join
import numpy as np
from osgeo import gdal
import pickle as cp

import pyrate.config as cf
import pyrate.linrate as linrate
import pyrate.mst as mst
import pyrate.orbital as orbital
import pyrate.prepifg as prepifg
import pyrate.refpixel as refpixel
import pyrate.timeseries as timeseries
from pyrate import algorithm
from pyrate import ifgconstants as ifc
from pyrate import matlab_mst as matlab_mst
from pyrate import reference_phase_estimation as rpe
from pyrate import vcm as vcm_module
from pyrate.shared import Ifg, write_output_geotiff, \
    pre_prepare_ifgs, prepare_ifgs_without_phase, create_tiles
from pyrate.compat import PyAPS_INSTALLED
from pyrate.nci.common_nci import save_latest_phase
from pyrate import mpiops
if PyAPS_INSTALLED:
    from pyrate import remove_aps_delay as aps

PYRATEPATH = cf.PYRATEPATH
MASTER_PROCESS = 0
# print screen output
VERBOSE = True
log = logging.getLogger(__name__)


def get_tiles(ifg_path, rows, cols):
    ifg = Ifg(ifg_path)
    ifg.open(readonly=True)
    tiles = create_tiles(ifg.shape, nrows=rows, ncols=cols)
    ifg.close()
    return tiles


class PrereadIfg:

    def __init__(self, path, nan_fraction, master, slave, time_span):
        self.path = path
        self.nan_fraction = nan_fraction
        self.master = master
        self.slave = slave
        self.time_span = time_span


def _join_dicts(dicts):
    if dicts is None:
        return
    d = {k: v for D in dicts for k, v in D.items()}
    return d


def convert_phase_to_numpy(dest_tifs, params, tiles):
    """
    1. Convert ifg phase data into numpy binary files.
    2. Save the preread_ifgs dict with information about the ifgs that are
    later used for fast loading of Ifg files in IfgPart class

    Parameters
    ----------
    dest_tifs: list
        list of destination tifs
    params: dict
        config dict
    tiles: list
        list of all Tile instances

    Returns
    -------
    preread_ifgs: str
        preread_ifgs file path
    """
    preread_ifgs_dict = {}
    process_tifs = np.array_split(dest_tifs, mpiops.size)[mpiops.rank]
    for d in process_tifs:
        ifg = save_latest_phase(d, tiles, params)
        nan_fraction = ifg.nan_fraction
        master = ifg.master
        slave = ifg.slave
        time_span = ifg.time_span
        ifg.close()

        preread_ifgs_dict[d] = {PrereadIfg(path=d,
                                           nan_fraction=nan_fraction,
                                           master=master,
                                           slave=slave,
                                           time_span=time_span)}

    preread_ifgs_dict = _join_dicts(
        mpiops.comm.gather(preread_ifgs_dict, root=0))

    preread_ifgs = join(params[cf.OUT_DIR], 'preread_ifgs.pk')

    if mpiops.rank == MASTER_PROCESS:
        cp.dump(preread_ifgs_dict, open(preread_ifgs, 'wb'))

    log.info('finish converting phase_data to numpy '
             'in process {}'.format(mpiops.rank))
    mpiops.comm.barrier()
    return preread_ifgs


def mpi_mst_calc(dest_tifs, params, tiles, preread_ifgs):
    """
    MPI function that control each process during MPI run
    :param MPI_myID:
    :param dest_tifs: paths of cropped amd resampled geotiffs
    :param mpi_log_filename:
    :param num_processors:
    :param parallel: MPI Parallel class instance
    :param params: config params dictionary
    :param mst_file: mst file (2d numpy array) save to disc
    :return:
    """

    log.info('Calculating mst')
    log.info('Calculating minimum spanning tree matrix '
             'using NetworkX method')
    process_tiles = np.array_split(tiles, mpiops.size)[mpiops.rank]

    # TODO: later change to use Tile.index instead of process indices
    process_indices = np.array_split(range(len(tiles)),
                                     mpiops.size)[mpiops.rank]

    def save_mst_tile(tile, i, preread_ifgs):
        mst_tile = mst.mst_multiprocessing(tile, dest_tifs, preread_ifgs)
        # locally save the mst_mat
        mst_file_process_n = os.path.join(
            params[cf.OUT_DIR], 'mst_mat_{}.npy'.format(i))
        np.save(file=mst_file_process_n, arr=mst_tile)

    for t, p_ind in zip(process_tiles, process_indices):
        save_mst_tile(t, p_ind, preread_ifgs)
    log.info('finished mst calculation for process {}'.format(mpiops.rank))


def process_ifgs(ifg_paths, params, rows, cols):
    """
    Top level function to perform PyRate correction steps on given ifgs
    ifgs: sequence of paths to interferrograms
    (NB: changes are saved into ifgs)
    params: dictionary of configuration parameters
    """
    tiles = get_tiles(ifg_paths[0], rows, cols)
    ifgs = pre_prepare_ifgs(ifg_paths, params)
    preread_ifgs = convert_phase_to_numpy(ifg_paths,
                                          params=params,
                                          tiles=tiles)

    if params[cf.NETWORKX_OR_MATLAB_FLAG]:   # Using matlab mst
        mst_grid = mst_calculation(ifg_paths, params)
        mpi_mst_calc(ifg_paths, params, tiles, preread_ifgs)
    else:
        ifg_instance = matlab_mst.IfgListPyRate(datafiles=ifg_paths)
        mst_grid = mst_calculation(ifg_instance, params)

    # Estimate reference pixel location
    refpx, refpy = find_reference_pixel(ifgs, params)

    # remove APS delay here, and write aps delay removed ifgs disc
    if PyAPS_INSTALLED and aps_delay_required(ifgs, params):
        ifgs = aps.remove_aps_delay(ifgs, params)
        log.info('Finished APS delay correction')

    # make sure aps correction flags are consistent
    if params[cf.APS_CORRECTION]:
        check_aps_ifgs(ifgs)

    # Estimate and remove orbit errors
    remove_orbital_error(ifgs, params)
    # open ifgs again, but without phase conversion as already converted and
    # saved to disc

    ifgs = prepare_ifgs_without_phase(ifg_paths)
    log.info('Estimating and removing phase at reference pixel')
    ref_phs, ifgs = rpe.estimate_ref_phase(ifgs, params, refpx, refpy)

    # save reference phase
    ref_phs_file = os.path.join(params[cf.OUT_DIR], 'ref_phs.npy')
    np.save(file=ref_phs_file, arr=ref_phs)

    # TODO: assign maxvar to ifg metadata (and geotiff)?
    log.info('Calculating maximum variance in interferograms')
    maxvar = [vcm_module.cvd(i, params)[0] for i in ifgs]
    maxvar_file = os.path.join(params[cf.OUT_DIR], 'maxvar.npy')
    np.save(file=maxvar_file, arr=maxvar)

    log.info('Constructing temporal variance-covariance matrix')
    vcmt = vcm_module.get_vcmt(ifgs, maxvar)

    # write vcm output to a file
    vcmt_mat_binary_file = os.path.join(
        PYRATEPATH, params[cf.OUT_DIR], 'vcmt_mat.npy')
    np.save(file=vcmt_mat_binary_file, arr=vcmt)

    if params[cf.TIME_SERIES_CAL] != 0:
        compute_time_series(ifgs, mst_grid, params, vcmt)

    # Calculate linear rate map
    rate, error, samples = calculate_linear_rate(ifgs, params, vcmt, mst_grid)

    # close all open ifgs
    for i in ifgs:
        i.close()

    log.info('PyRate workflow completed')
    return mst_grid, (refpx, refpy), maxvar, vcmt, rate, error, samples


def write_linrate_numpy_files(error, params, rate, samples):
    rate_file = os.path.join(params[cf.OUT_DIR], 'rate.npy')
    error_file = os.path.join(params[cf.OUT_DIR], 'error.npy')
    samples_file = os.path.join(params[cf.OUT_DIR], 'samples.npy')
    np.save(file=rate_file, arr=rate)
    np.save(file=error_file, arr=error)
    np.save(file=samples_file, arr=samples)


def aps_delay_required(ifgs, params):
    log.info('Removing APS delay')

    if not params[cf.APS_CORRECTION]:
        log.info('APS delay removal not required')
        return False

    # perform some general error/sanity checks
    flags = [i.dataset.GetMetadataItem(ifc.PYRATE_APS_ERROR) for i in ifgs]

    if all(flags):
        log.info('Skipped APS delay removal, ifgs are already aps corrected')
        return False
    else:
        check_aps_ifgs(ifgs)

    return True


def check_aps_ifgs(ifgs):
    flags = [i.dataset.GetMetadataItem(ifc.PYRATE_APS_ERROR) for i in ifgs]
    count = sum([f == aps.APS_STATUS for f in flags])
    if (count < len(flags)) and (count > 0):
        log.debug('Detected mix of corrected and uncorrected '
                      'APS delay in ifgs')

        for i, flag in zip(ifgs, flags):
            if flag:
                msg = '%s: prior APS delay correction detected'
            else:
                msg = '%s: no APS delay correction detected'
            logging.debug(msg % i.data_path)
        raise aps.APSException('Mixed APS removal status in ifg list')


def mst_calculation(ifg_paths_or_instance, params):
    if isinstance(ifg_paths_or_instance, list):
        ifgs = pre_prepare_ifgs(ifg_paths_or_instance, params)
        log.info('Calculating minimum spanning tree matrix using '
                 'NetworkX')

        mst_grid = mst.mst_parallel(ifgs, params)
    else:
        nan_conversion = params[cf.NAN_CONVERSION]
        assert isinstance(ifg_paths_or_instance, matlab_mst.IfgListPyRate)
        ifgs = ifg_paths_or_instance.ifgs
        for i in ifgs:
            if not i.mm_converted:
                i.nodata_value = params[cf.NO_DATA_VALUE]
                i.convert_to_mm()
        ifg_instance_updated, epoch_list = \
            matlab_mst.get_nml(ifg_paths_or_instance,
                               nodata_value=params[cf.NO_DATA_VALUE],
                               nan_conversion=nan_conversion)
        log.info('Calculating minimum spanning tree matrix '
                 'using Matlab-algorithm method')
        mst_grid = matlab_mst.matlab_mst_boolean_array(ifg_instance_updated)

    # write mst output to a file
    mst_mat_binary_file = os.path.join(
        PYRATEPATH, params[cf.OUT_DIR], 'mst_mat')
    np.save(file=mst_mat_binary_file, arr=mst_grid)

    for i in ifgs:
        i.close()
    return mst_grid


def compute_time_series(ifgs, mst_grid, params, vcmt):
    log.info('Calculating time series')
    # Calculate time series
    tsincr, tscum, tsvel = calculate_time_series(
        ifgs, params, vcmt=vcmt, mst=mst_grid)

    # tsvel_file = os.path.join(params[cf.OUT_DIR], 'tsvel.npy')
    tsincr_file = os.path.join(params[cf.OUT_DIR], 'tsincr.npy')
    tscum_file = os.path.join(params[cf.OUT_DIR], 'tscum.npy')
    np.save(file=tsincr_file, arr=tsincr)
    np.save(file=tscum_file, arr=tscum)
    # np.save(file=tsvel_file, arr=tsvel)

    # TODO: write tests for these functions
    write_timeseries_geotiff(ifgs, params, tsincr, pr_type='tsincr')
    write_timeseries_geotiff(ifgs, params, tscum, pr_type='tscuml')
    # write_timeseries_geotiff(ifgs, params, tsvel, pr_type='tsvel')
    return tsincr, tscum, tsvel


def setup_metadata(ifgs, params):
    p = os.path.join(params[cf.OUT_DIR], ifgs[0].data_path)
    ds = gdal.Open(p)
    md = ds.GetMetadata()  # get metadata for writing on output tifs
    gt = ds.GetGeoTransform()  # get geographical bounds of data
    wkt = ds.GetProjection()  # get projection of data
    epochlist = algorithm.get_epochs(ifgs)
    return epochlist, gt, md, wkt


def write_timeseries_geotiff(ifgs, params, tsincr, pr_type):
    # setup metadata for writing into result files
    epochlist, gt, md, wkt = setup_metadata(ifgs, params)
    for i in range(tsincr.shape[2]):
        md[ifc.MASTER_DATE] = epochlist.dates[i + 1]
        md['PR_SEQ_POS'] = i  # sequence position

        data = tsincr[:, :, i]
        dest = os.path.join(
            PYRATEPATH, params[cf.OUT_DIR],
            pr_type + "_" + str(epochlist.dates[i + 1]) + ".tif")
        md[ifc.PRTYPE] = pr_type
        write_output_geotiff(md, gt, wkt, data, dest, np.nan)


def insert_time_series_interpolation(ifg_instance_updated, params):

    edges = matlab_mst.get_sub_structure(ifg_instance_updated,
                                  np.zeros(len(ifg_instance_updated.id), dtype=bool))

    _, _, ntrees = matlab_mst.matlab_mst_kruskal(edges, ntrees=True)
    # if ntrees=1, no interpolation; otherwise interpolate
    if ntrees > 1:
        params[cf.TIME_SERIES_INTERP] = 1
    else:
        params[cf.TIME_SERIES_INTERP] = 0

    return params


def remove_orbital_error(ifgs, params):
    log.info('Calculating orbital error correction')

    if not params[cf.ORBITAL_FIT]:
        log.info('Orbital correction not required')
        return

    # perform some general error/sanity checks
    flags = [i.dataset.GetMetadataItem(ifc.PYRATE_ORBITAL_ERROR) for i in ifgs]

    if all(flags):
        log.info('Skipped orbital correction, ifgs already corrected')
        return
    else:
        check_orbital_ifgs(ifgs, flags)

    mlooked = None

    if (params[cf.ORBITAL_FIT_LOOKS_X] > 1 or
                params[cf.ORBITAL_FIT_LOOKS_Y] > 1):
        # resampling here to use all prior corrections to orig data
        # can't do multiprocessing without writing to disc, but can do MPI
        # due to swig pickling issue. So multiprocesing is not implemented
        mlooked_dataset = prepifg.prepare_ifgs(
            [i.data_path for i in ifgs],
            crop_opt=prepifg.ALREADY_SAME_SIZE,
            xlooks=params[cf.ORBITAL_FIT_LOOKS_X],
            ylooks=params[cf.ORBITAL_FIT_LOOKS_Y],
            thresh=params[cf.NO_DATA_AVERAGING_THRESHOLD],
            write_to_disc=False)
        mlooked = [Ifg(m[1]) for m in mlooked_dataset]

        for m in mlooked:
            m.initialize()
            m.nodata_value = params[cf.NO_DATA_VALUE]

    orbital.orbital_correction(ifgs,
                               params,
                               mlooked=mlooked)


def check_orbital_ifgs(ifgs, flags):
    count = sum([f == ifc.ORB_REMOVED for f in flags])
    if (count < len(flags)) and (count > 0):
        logging.debug('Detected mix of corrected and uncorrected '
                      'orbital error in ifgs')

        for i, flag in zip(ifgs, flags):
            if flag:
                msg = '{}: prior orbital error correction detected'.format(i)
            else:
                msg = '{}: no orbital correction detected'.format(i)
            logging.debug(msg % i.data_path)

        raise orbital.OrbitalError(msg)


def find_reference_pixel(ifgs, params):
    # unlikely, but possible the refpixel can be (0,0)
    # check if there is a pre-specified reference pixel coord
    log.info('Starting reference pixel calculation')
    refx = params[cf.REFX]
    if refx > ifgs[0].ncols - 1:
        raise ValueError("Invalid reference pixel X coordinate: %s" % refx)

    refy = params[cf.REFY]
    if refy > ifgs[0].nrows - 1:
        raise ValueError("Invalid reference pixel Y coordinate: %s" % refy)

    if refx == 0 or refy == 0:  # matlab equivalent
        refy, refx = refpixel.ref_pixel(ifgs, params)
        log.info('Found reference pixel coordinate: (%s, %s)' % (refx, refy))
    else:
        log.info('Reusing config file reference pixel (%s, %s)' % (refx, refy))

    return refx, refy


def calculate_linear_rate(ifgs, params, vcmt, mst=None):
    log.info('Calculating linear rate')

    # TODO: do these need to be checked?
    res = linrate.linear_rate(ifgs, params, vcmt, mst)
    for r in res:
        if r is None:
            raise ValueError('TODO: bad value')

    rate, error, samples = res

    write_linrate_tifs(ifgs, params, res)

    log.info('Linear rate calculated')
    return rate, error, samples


def write_linrate_tifs(ifgs, params, res):
    log.info('Writing linrate results')
    rate, error, samples = res
    epochlist, gt, md, wkt = setup_metadata(ifgs, params)
    # TODO: write tests for these functions
    dest = os.path.join(PYRATEPATH, params[cf.OUT_DIR], "linrate.tif")
    md[ifc.MASTER_DATE] = epochlist.dates
    md[ifc.PRTYPE] = 'linrate'
    write_output_geotiff(md, gt, wkt, rate, dest, np.nan)
    dest = os.path.join(PYRATEPATH, params[cf.OUT_DIR], "linerror.tif")
    md[ifc.PRTYPE] = 'linerror'
    write_output_geotiff(md, gt, wkt, error, dest, np.nan)
    write_linrate_numpy_files(error, params, rate, samples)


def calculate_time_series(ifgs, params, vcmt, mst):
    res = timeseries.time_series(ifgs, params, vcmt, mst)
    for r in res:
        if len(r.shape) != 3:
            logging.error('TODO: time series result shape is incorrect')
            raise timeseries.TimeSeriesError

    log.info('Time series calculated')
    tsincr, tscum, tsvel = res
    return tsincr, tscum, tsvel


# general function template
#
# add check for pre-existing metadata flag / skip if required
# perform calculation
# optionally save modified data to disk if required
# optionally save correction component to disk (more useful for debugging)
# set flag in dataset for correction
# write to log file


def warp_required(xlooks, ylooks, crop):
    """
    Returns True if params show rasters need to be cropped and/or resized.
    """

    if xlooks > 1 or ylooks > 1:
        return True

    if crop is None:
        return False

    return True


def working_ifg_paths(src_paths, xlooks, ylooks, cropping):
    """
    Filter. Returns paths to ifgs to process (eg. checks for mlooked tifs)
    """
    if warp_required(xlooks, ylooks, cropping):
        mlooked_unw = [cf.mlooked_path(p, xlooks, crop_out=cropping)
                       for p in src_paths]
        mlooked_paths = [os.path.splitext(m)[0]+'.tif' for m in mlooked_unw]

        if not all([os.path.exists(p) for p in mlooked_paths]):
            msg = 'Multilooked ifgs do not exist ' \
                  '(execute "run_prepifg.py" first)'
            raise IOError(msg)

        log.info('Using mlooked interferograms...')
        return mlooked_paths
    return src_paths  # multi looking not specified, work with base ifgs


def dest_ifg_paths(ifg_paths, outdir):
    """
    Returns paths to out/dest ifgs.
    """

    bases = [os.path.basename(p) for p in ifg_paths]
    return [os.path.join(outdir, p) for p in bases]


def main(config_file, rows, cols):
    base_unw_paths, dest_paths, pars = cf.get_ifg_paths(config_file)
    process_ifgs(dest_paths, pars, rows, cols)


def log_config_file(configfile, log_filename):
    output_log_file = open(log_filename, "a")
    output_log_file.write("\nConfig Settings: start\n")
    lines = open(configfile).read()
    for line in lines:
        output_log_file.write(line)
    output_log_file.write("\nConfig Settings: end\n\n")
    output_log_file.write("\n==============================================\n")
