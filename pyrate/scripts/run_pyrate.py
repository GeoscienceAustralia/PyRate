"""
Main workflow script for PyRate
"""

import datetime
import logging
import os
import sys
import numpy as np

import pyrate.config as cf
import pyrate.linrate as linrate
import pyrate.mst as mst
import pyrate.orbital as orbital
import pyrate.prepifg as prepifg
import pyrate.refpixel as refpixel
import pyrate.timeseries as timeseries
from osgeo import gdal
from pyrate import algorithm
from pyrate import ifgconstants as ifc
from pyrate import matlab_mst_kruskal as matlab_mst
from pyrate import reference_phase_estimation as rpe
from pyrate import remove_aps_delay as aps
from pyrate import vcm as vcm_module
from pyrate.shared import Ifg, write_output_geotiff, pre_prepare_ifgs

PYRATEPATH = cf.PYRATEPATH

# print screen output
VERBOSE = True


def process_ifgs(ifg_paths_or_instance, params):
    """
    Top level function to perform PyRate correction steps on given ifgs
    ifgs: sequence of paths to interferrograms (NB: changes are saved into ifgs)
    params: dictionary of configuration parameters
    """
    mst_grid = mst_calculation(ifg_paths_or_instance, params)

    # reading ifgs again, this is consistent with nci submission script
    ifgs = pre_prepare_ifgs(ifg_paths_or_instance, params)

    # Estimate reference pixel location
    refpx, refpy = find_reference_pixel(ifgs, params)

    # remove APS delay here, and write aps delay removed ifgs disc
    if aps_delay_required(ifgs, params):
        ifgs = aps.remove_aps_delay(ifgs, params)
        print 'Finished APS delay correction'

    # make sure aps correction flags are consistent
    if params[cf.APS_CORRECTION]:
        check_aps_ifgs(ifgs)

    # Estimate and remove orbit errors
    remove_orbital_error(ifgs, params)

    for i in ifgs:
        i.close()

    ifgs = pre_prepare_ifgs(ifg_paths_or_instance, params)

    write_msg('Estimating and removing phase at reference pixel')
    ref_phs, ifgs = rpe.estimate_ref_phase(ifgs, params, refpx, refpy)

    # save reference phase
    ref_phs_file = os.path.join(params[cf.OUT_DIR], 'ref_phs.npy')
    np.save(file=ref_phs_file, arr=ref_phs)

    # TODO: assign maxvar to ifg metadata (and geotiff)?
    write_msg('Calculating maximum variance in interferograms')
    maxvar = [vcm_module.cvd(i)[0] for i in ifgs]
    maxvar_file = os.path.join(params[cf.OUT_DIR], 'maxvar.npy')
    np.save(file=maxvar_file, arr=maxvar)

    write_msg('Constructing temporal variance-covariance matrix')
    vcmt = vcm_module.get_vcmt(ifgs, maxvar)

    # write vcm output to a file
    vcmt_mat_binary_file = os.path.join(
        PYRATEPATH, params[cf.OUT_DIR], 'vcmt_mat.npy')
    np.save(file=vcmt_mat_binary_file, arr=vcmt)

    p = os.path.join(params[cf.OUT_DIR], ifgs[0].data_path)
    assert os.path.exists(p) == True

    ds = gdal.Open(p)
    md = ds.GetMetadata()  # get metadata for writing on output tifs
    gt = ds.GetGeoTransform()  # get geographical bounds of data
    wkt = ds.GetProjection()  # get projection of data
    epochlist = algorithm.get_epochs(ifgs)

    if params[cf.TIME_SERIES_CAL] != 0:
        compute_time_series(epochlist, gt, ifgs, md, mst_grid, params, vcmt,
                            wkt)
    # Calculate linear rate map
    rate, error, samples = calculate_linear_rate(
                   ifgs, params, vcmt, mst=mst_grid)

    md[ifc.MASTER_DATE] = epochlist.dates
    dest = os.path.join(PYRATEPATH, params[cf.OUT_DIR], "linrate.tif")
    # remove metadata added to md in compute_time_series that doesn't
    # make sense for the following tiffs
    if 'PR_SEQ_POS' in md:
        del md['PR_SEQ_POS']
    md['PR_TYPE'] = 'linrate'
    write_output_geotiff(md, gt, wkt, rate, dest, np.nan)
    dest = os.path.join(PYRATEPATH, params[cf.OUT_DIR], "linerror.tif")

    md['PR_TYPE'] = 'linerror'
    write_output_geotiff(md, gt, wkt, error, dest, np.nan)

    write_linrate_numpy_files(error, params, rate, samples)

    # close all open ifgs
    for i in ifgs:
        i.close()

    write_msg('PyRate workflow completed')
    return mst_grid, (refpx, refpy), maxvar, vcmt, rate, error, samples


def write_linrate_numpy_files(error, params, rate, samples):
    rate_file = os.path.join(params[cf.OUT_DIR], 'rate.npy')
    error_file = os.path.join(params[cf.OUT_DIR], 'error.npy')
    samples_file = os.path.join(params[cf.OUT_DIR], 'samples.npy')
    np.save(file=rate_file, arr=rate)
    np.save(file=error_file, arr=error)
    np.save(file=samples_file, arr=samples)


def aps_delay_required(ifgs, params):
    write_msg('Removing APS delay')

    if not params[cf.APS_CORRECTION]:
        write_msg('APS delay removal not required')
        return False

    # perform some general error/sanity checks
    flags = [i.dataset.GetMetadataItem(ifc.PYRATE_APS_ERROR) for i in ifgs]

    if all(flags):
        write_msg('Skipped APS delay removal, ifgs are already aps corrected')
        return False
    else:
        check_aps_ifgs(ifgs)

    return True


def check_aps_ifgs(ifgs):
    flags = [i.dataset.GetMetadataItem(ifc.PYRATE_APS_ERROR) for i in ifgs]
    count = sum([f == aps.APS_STATUS for f in flags])
    if (count < len(flags)) and (count > 0):
        logging.debug('Detected mix of corrected and uncorrected '
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
        write_msg(
            'Calculating minimum spanning tree matrix using NetworkX method')

        mst_grid = mst.mst_parallel(ifgs, params)
    else:
        # the matlab side has not been worked for a while, may need updating
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
        write_msg(
            'Calculating minimum spanning tree matrix '
            'using Matlab-algorithm method')
        mst_grid = matlab_mst.matlab_mst_boolean_array(ifg_instance_updated)

        # Insert INTERP into the params for timeseries calculation
        # params = insert_time_series_interpolation(ifg_instance_updated, params)

    # write mst output to a file
    mst_mat_binary_file = os.path.join(
        PYRATEPATH, params[cf.OUT_DIR], 'mst_mat')
    np.save(file=mst_mat_binary_file, arr=mst_grid)

    for i in ifgs:
        i.close()
    return mst_grid


def compute_time_series(epochlist, gt, ifgs, md, mst_grid, params, vcmt, wkt):
    # Calculate time series
    tsincr, tscum, tsvel = calculate_time_series(
        ifgs, params, vcmt=vcmt, mst=mst_grid)

    # TODO: write tests for these functions
    write_timeseries_geotiff(epochlist, gt, md, params, tsincr, wkt,
                             pr_type='tsincr')
    write_timeseries_geotiff(epochlist, gt, md, params, tscum, wkt,
                             pr_type='tscuml')
    write_timeseries_geotiff(epochlist, gt, md, params, tsvel, wkt,
                             pr_type='tsvel')
    return tsincr, tscum, tsvel


def write_timeseries_geotiff(epochlist, gt, md, params, tsincr, wkt, pr_type):
    PRTYPE = 'PR_TYPE'
    for i in range(len(tsincr[0, 0, :])):
        md[ifc.MASTER_DATE] = epochlist.dates[i + 1]
        md['PR_SEQ_POS'] = i  # sequence position

        data = tsincr[:, :, i]
        dest = os.path.join(
            PYRATEPATH, params[cf.OUT_DIR],
            pr_type + "_" + str(epochlist.dates[i + 1]) + ".tif")
        md[PRTYPE] = pr_type
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
    write_msg('Calculating orbital error correction')

    if not params[cf.ORBITAL_FIT]:
        write_msg('Orbital correction not required')
        return

    # perform some general error/sanity checks
    flags = [i.dataset.GetMetadataItem(ifc.PYRATE_ORBITAL_ERROR) for i in ifgs]

    if all(flags):
        write_msg('Skipped orbital correction, ifgs already corrected')
        return
    else:
        check_orbital_ifgs(ifgs, flags)

    mlooked = None

    if (params[cf.ORBITAL_FIT_LOOKS_X] > 1 or
                params[cf.ORBITAL_FIT_LOOKS_Y] > 1):
        # resampling here to use all prior corrections to orig data
        # can't do multiprocessing without writing to disc, but can do MPI
        # due to swig pickling issue. So multiprocesing is not implemented
        mlooked_dataset = prepifg.prepare_ifgs([i.data_path for i in ifgs],
                             crop_opt=prepifg.ALREADY_SAME_SIZE,
                             xlooks=params[cf.ORBITAL_FIT_LOOKS_X],
                             ylooks=params[cf.ORBITAL_FIT_LOOKS_Y],
                             thresh=params[cf.NO_DATA_AVERAGING_THRESHOLD],
                             write_to_disc=False)[0]
        mlooked = [Ifg(m) for m in mlooked_dataset]

        for m in mlooked:
            m.initialize()
            m.nodata_value = params[cf.NO_DATA_VALUE]

    orbital.orbital_correction(ifgs,
                               params,
                               mlooked=mlooked)

    # write data to disc after orbital error correction
    for i in ifgs:
        i.dataset.SetMetadataItem(ifc.PYRATE_ORBITAL_ERROR, ifc.ORB_REMOVED)
        i.write_modified_phase()
        logging.debug('%s: orbital error removed' % i.data_path)


def check_orbital_ifgs(ifgs, flags):
    count = sum([f == ifc.ORB_REMOVED for f in flags])
    if (count < len(flags)) and (count > 0):
        logging.debug('Detected mix of corrected and uncorrected '
                      'orbital error in ifgs')

        for i, flag in zip(ifgs, flags):
            if flag:
                msg = '%s: prior orbital error correction detected'
            else:
                msg = '%s: no orbital correction detected'
            logging.debug(msg % i.data_path)

        raise orbital.OrbitalError(msg)


def find_reference_pixel(ifgs, params):
    # unlikely, but possible the refpixel can be (0,0)
    # check if there is a pre-specified reference pixel coord
    refx = params[cf.REFX]
    if refx > ifgs[0].ncols - 1:
        raise ValueError("Invalid reference pixel X coordinate: %s" % refx)

    refy = params[cf.REFY]
    if refy > ifgs[0].nrows - 1:
        raise ValueError("Invalid reference pixel Y coordinate: %s" % refy)

    if refx == 0 or refy == 0:  # matlab equivalent
        write_msg('Finding reference pixel')
        refy, refx = refpixel.ref_pixel(ifgs, params)
        write_msg('Reference pixel coordinate: (%s, %s)' % (refx, refy))
    else:
        write_msg('Reusing config file reference pixel (%s, %s)' % (refx, refy))

    return refx, refy


def calculate_linear_rate(ifgs, params, vcmt, mst=None):
    write_msg('Calculating linear rate')

    # TODO: do these need to be checked?
    res = linrate.linear_rate(ifgs, params, vcmt, mst)
    for r in res:
        if r is None:
            raise ValueError('TODO: bad value')

    rate, error, samples = res

    logging.debug('Linear rate calculated')
    return rate, error, samples


def calculate_time_series(ifgs, params, vcmt, mst):
    write_msg('Calculating time series')
    res = timeseries.time_series(ifgs, params, vcmt, mst)
    for r in res:
        if len(r.shape) != 3:
            logging.error('TODO: time series result shape is incorrect')
            raise timeseries.TimeSeriesError

    logging.debug('Time series calculated')
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


def transform_params(params):
    """
    Returns subset of config parameters for cropping and multilooking.
    """

    t_params = [cf.IFG_LKSX, cf.IFG_LKSY, cf.IFG_CROP_OPT]
    xlooks, ylooks, crop = [params[k] for k in t_params]
    return xlooks, ylooks, crop


def warp_required(xlooks, ylooks, crop):
    """
    Returns True if params show rasters need to be cropped and/or resized.
    """

    if xlooks > 1 or ylooks > 1:
        return True

    if crop is None:
        return False

    return True


def original_ifg_paths(ifglist_path):
    """
    Returns sequence of paths to files in given ifglist file.
    """

    basedir = os.path.dirname(ifglist_path)
    ifglist = cf.parse_namelist(os.path.join(PYRATEPATH, ifglist_path))
    return [os.path.join(basedir, p) for p in ifglist]


def working_ifg_paths(src_paths, xlooks, ylooks, cropping):
    """
    Filter. Returns paths to ifgs to process (eg. checks for mlooked tifs)
    """
    if warp_required(xlooks, ylooks, cropping):
        mlooked_unw = [prepifg.mlooked_path(p, xlooks, crop_out=cropping)
                       for p in src_paths]
        mlooked_paths = [os.path.splitext(m)[0]+'.tif' for m in mlooked_unw]

        if not all([os.path.exists(p) for p in mlooked_paths]):
            msg = 'Multilooked ifgs do not exist (execute "run_prepifg.py" first)'
            raise IOError(msg)

        logging.debug('Using mlooked interferograms...')
        return mlooked_paths
    return src_paths  # multi looking not specified, work with base ifgs


def dest_ifg_paths(ifg_paths, outdir):
    """
    Returns paths to out/dest ifgs.
    """

    bases = [os.path.basename(p) for p in ifg_paths]
    return [os.path.join(outdir, p) for p in bases]


def init_logging(level):
    t = datetime.datetime.now()
    path = 'pyrate_%s.log' % t.isoformat().replace(':', '_')
    fmt = '%(asctime)s %(message)s'
    datefmt = '%d/%m/%Y %I:%M:%S %p'
    logging.basicConfig(filename=path, format=fmt, datefmt=datefmt, level=level)
    logging.debug('Log started')
    return path


def get_dest_paths(base_paths, crop, params, looks):
    dest_mlooked_ifgs = [prepifg.mlooked_path(os.path.basename(q).split('.')[0]
        + '.tif', looks=looks, crop_out=crop) for q in base_paths]

    return [os.path.join(os.environ['PYRATEPATH'], params[cf.OUT_DIR], p)
            for p in dest_mlooked_ifgs]


# TODO: expand CLI documentation
# TODO: ensure clean exception handling
# TODO: add parameter error checking: induce fail fast before number crunching
def main():
    base_unw_paths, dest_paths, pars = get_ifg_paths()

    if pars[cf.NETWORKX_OR_MATLAB_FLAG]:  # Using networkx mst
        process_ifgs(dest_paths, pars)
    else:  # Using matlab mst
        ifg_instance = matlab_mst.IfgListPyRate(datafiles=dest_paths)
        process_ifgs(ifg_instance, pars)


def get_ifg_paths():
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [config-file]\nRuns PyRate workflow.')
    parser.add_option('-i', '--ifglist', type=str,
                      help='name of file containing list of interferograms')
    options, args = parser.parse_args()
    log_file_name = init_logging(logging.DEBUG)
    try:
        cfg_path = args[0] if args else 'pyrate.conf'
        global pars
        pars = cf.get_config_params(cfg_path)

    except IOError as err:
        emsg = 'Config file error: %s "%s"' % (err.strerror, err.filename)
        logging.debug(emsg)
        sys.exit(err.errno)
    ifg_file_list = options.ifglist or pars.get(cf.IFG_FILE_LIST)
    pars[cf.IFG_FILE_LIST] = ifg_file_list

    # log config file
    log_config_file(cfg_path, log_file_name)

    if ifg_file_list is None:
        emsg = 'Error {code}: Interferogram list file name not provided ' \
               'or does not exist'.format(code=2)
        logging.debug(emsg)
        raise IOError(2, emsg)
    xlks, ylks, crop = transform_params(pars)

    # base_unw_paths need to be geotiffed and multilooked by run_prepifg
    base_unw_paths = original_ifg_paths(ifg_file_list)

    # dest_paths are tifs that have been geotif converted and multilooked
    dest_paths = get_dest_paths(base_unw_paths, crop, pars, xlks)

    return base_unw_paths, dest_paths, pars


def log_config_file(configfile, log_filename):
    output_log_file = open(log_filename, "a")
    output_log_file.write("\nConfig Settings: start\n")
    lines = open(configfile).read()
    for line in lines:
        output_log_file.write(line)
    output_log_file.write("\nConfig Settings: end\n\n")
    output_log_file.write("\n===============================================\n")


def write_msg(msg):
    """
    write message to log file and screen output
    """
    logging.debug(msg)
    if VERBOSE:
        print msg


if __name__ == "__main__":
    main()
