"""
Main workflow script for PyRate

Created on 17/09/2012
@author: Ben Davies, NCI
@author: Sudipta Basak, GA
"""

import os, sys, shutil, logging, datetime
from osgeo import gdal
import numpy as np

import pyrate.mst as mst
import pyrate.prepifg as prepifg
import pyrate.refpixel as refpixel
import pyrate.orbital as orbital
import pyrate.linrate as linrate
import pyrate.timeseries as timeseries
import pyrate.config as cf
from pyrate.shared import Ifg, write_output_geotiff
from pyrate import vcm as vcm_module
from pyrate import matlab_mst_kruskal as matlab_mst
from pyrate import reference_phase_estimation as rpe
from pyrate import algorithm
from pyrate import ifgconstants as ifc


# constants for metadata flags
ORB_REMOVED = 'REMOVED'
pars = None
PYRATEPATH = cf.PYRATEPATH

# print screen output
VERBOSE = True


def process_ifgs(ifg_paths_or_instance, params):
    """
    Top level function to perform PyRate correction steps on given ifgs
    ifgs: sequence of paths to interferograms (NB: changes are saved into ifgs)
    params: dictionary of configuration parameters
    """
    ifgs, mst_grid, params = mst_calculation(ifg_paths_or_instance, params)

    # Estimate reference pixel location
    refpx, refpy = find_reference_pixel(ifgs, params)

    # remove APS delay here

    # Estimate and remove orbit errors
    if params[cf.ORBITAL_FIT] != 0:
        remove_orbital_error(ifgs, params)

    write_msg('Estimating and removing phase at reference pixel')
    _, ifgs = rpe.estimate_ref_phase(ifgs, params, refpx, refpy)

    # TODO: assign maxvar to ifg metadata (and geotiff)?
    write_msg('Calculating maximum variance in interferograms')
    maxvar = [vcm_module.cvd(i)[0] for i in ifgs]

    write_msg('Constructing temporal variance-covariance matrix')
    vcmt = vcm_module.get_vcmt(ifgs, maxvar)

    # write vcm output to a file
    vcmt_mat_binary_file = os.path.join(
        PYRATEPATH, params[cf.OUT_DIR], 'vcmt_mat')
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
    md[ifc.PYRATE_DATE] = epochlist.dates
    dest = os.path.join(PYRATEPATH, params[cf.OUT_DIR], "linrate.tif")
    write_output_geotiff(md, gt, wkt, rate, dest, np.nan)
    dest = os.path.join(PYRATEPATH, params[cf.OUT_DIR], "linerror.tif")
    write_output_geotiff(md, gt, wkt, error, dest, np.nan)

    write_msg('PyRate workflow completed')


def mst_calculation(ifg_paths_or_instance, params):
    if isinstance(ifg_paths_or_instance, list):
        ifgs = prepare_ifgs_for_networkx_mst(ifg_paths_or_instance, params)
        write_msg(
            'Calculating minimum spanning tree matrix using NetworkX method')

        mst_grid = mst.mst_parallel(ifgs, params)
        # check if mst is not a tree, then do interpolate
        if mst.mst_from_ifgs(ifgs)[1]:
            params[cf.TIME_SERIES_INTERP] = 0
        else:
            params[cf.TIME_SERIES_INTERP] = 1
    else:
        nan_conversion = params[cf.NAN_CONVERSION]
        assert isinstance(ifg_paths_or_instance, matlab_mst.IfgListPyRate)
        ifgs = ifg_paths_or_instance.ifgs
        for i in ifgs:
            if not i.mm_converted:
                i.nodata_value = params[cf.NO_DATA_VALUE]
                i.convert_to_mm()
                i.write_modified_phase()
        ifg_instance_updated, epoch_list = \
            matlab_mst.get_nml(ifg_paths_or_instance,
                               nodata_value=params[cf.NO_DATA_VALUE],
                               nan_conversion=nan_conversion)
        write_msg(
            'Calculating minimum spanning tree matrix using Matlab-algorithm method')
        mst_grid = matlab_mst.matlab_mst_boolean_array(ifg_instance_updated)

        # Insert INTERP into the params for timeseries calculation
        params = insert_time_series_interpolation(ifg_instance_updated, params)

    # make sure by now we have the time series interpolation parameter
    assert params[cf.TIME_SERIES_INTERP] is not None

    # write mst output to a file
    mst_mat_binary_file = os.path.join(
        PYRATEPATH, params[cf.OUT_DIR], 'mst_mat')
    np.save(file=mst_mat_binary_file, arr=mst_grid)

    return ifgs, mst_grid, params


def prepare_ifgs_for_networkx_mst(ifg_paths_or_instance, params):
    nan_conversion = params[cf.NAN_CONVERSION]
    ifgs = [Ifg(p) for p in ifg_paths_or_instance]
    for i in ifgs:
        if not i.is_open:
            i.open(readonly=False)
        if nan_conversion:  # nan conversion happens here in networkx mst
            i.nodata_value = params[cf.NO_DATA_VALUE]
            i.convert_to_nans()
        if not i.mm_converted:
            i.convert_to_mm()
            i.write_modified_phase()
    return ifgs


def compute_time_series(epochlist, gt, ifgs, md, mst_grid, params, vcmt, wkt):
    # Calculate time series
    tsincr, tscum, tsvel = calculate_time_series(
        ifgs, params, vcmt=vcmt, mst=mst_grid)
    for i in range(len(tsincr[0, 0, :])):
        md[ifc.PYRATE_DATE] = epochlist.dates[i + 1]
        data = tsincr[:, :, i]
        dest = os.path.join(
            PYRATEPATH, params[cf.OUT_DIR],
            "tsincr_" + str(epochlist.dates[i + 1]) + ".tif")
        write_output_geotiff(md, gt, wkt, data, dest, np.nan)

        data = tscum[:, :, i]
        dest = os.path.join(
            PYRATEPATH, params[cf.OUT_DIR],
            "tscuml_" + str(epochlist.dates[i + 1]) + ".tif")
        write_output_geotiff(md, gt, wkt, data, dest, np.nan)

        data = tsvel[:, :, i]
        dest = os.path.join(
            PYRATEPATH, params[cf.OUT_DIR],
            "tsvel_" + str(epochlist.dates[i + 1]) + ".tif")
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
    if (params[cf.ORBITAL_FIT_LOOKS_X] > 1 or params[cf.ORBITAL_FIT_LOOKS_Y] > 1)\
            and (params[cf.ORBITAL_FIT_METHOD] != 1):
        # resampling here to use all prior corrections to orig data
        mlooked_phase_data = prepifg.prepare_ifgs([i.data_path for i in ifgs],
                             crop_opt=prepifg.ALREADY_SAME_SIZE,
                             xlooks=params[cf.ORBITAL_FIT_LOOKS_X],
                             ylooks=params[cf.ORBITAL_FIT_LOOKS_Y],
                             write_to_disc=False)
        mlooked = [Ifg(m) for m in mlooked_phase_data]

        for m, i in zip(mlooked, ifgs):
            m.initialize()
            m.nodata_value = params[cf.NO_DATA_VALUE]

    orbital.orbital_correction(ifgs,
                               degree=params[cf.ORBITAL_FIT_DEGREE],
                               method=params[cf.ORBITAL_FIT_METHOD],
                               mlooked=mlooked)

    for i in ifgs:
        i.dataset.SetMetadataItem(ifc.PYRATE_ORBITAL_ERROR, ORB_REMOVED)
        logging.debug('%s: orbital error removed' % i.data_path)


def check_orbital_ifgs(ifgs, flags):
    count = sum([f == ORB_REMOVED for f in flags])
    if (count < len(flags)) and (count > 0):
        logging.debug('Detected mix of corrected and uncorrected orbital error in ifgs')

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
        refy, refx = refpixel.ref_pixel(ifgs, params[cf.REFNX],
            params[cf.REFNY], params[cf.REF_CHIP_SIZE], params[cf.REF_MIN_FRAC],
                                        params[cf.PARALLEL])
        write_msg('Reference pixel coordinate: (%s, %s)' % (refx, refy))
    else:
        write_msg('Reusing config file reference pixel (%s, %s)' % (refx, refy))

    return refx, refy


def calculate_linear_rate(ifgs, params, vcmt, mst=None):
    write_msg('Calculating linear rate')

    # MULTIPROCESSING parameters
    parallel = params[cf.PARALLEL]
    processes = params[cf.PROCESSES]

    # TODO: do these need to be checked?
    res = linrate.linear_rate(ifgs, params, vcmt, mst,
                              parallel=parallel, processes=processes)
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
    init_logging(logging.DEBUG)
    try:
        cfg_path = args[0] if args else 'pyrate.conf'
        global pars
        pars = cf.get_config_params(cfg_path)

    except IOError as err:
        emsg = 'Config file error: %s "%s"' % (err.strerror, err.filename)
        logging.debug(emsg)
        sys.exit(err.errno)
    ifg_file_list = options.ifglist or pars.get(cf.IFG_FILE_LIST)
    if ifg_file_list is None:
        emsg = 'Error {code}: Interferogram list file name not provided ' \
               'or does not exist'.format(code=2)
        logging.debug(emsg)
        raise IOError(2, emsg)
    xlks, ylks, crop = transform_params(pars)

    # base_unw_paths need to be geotiffed and multilooked by run_prepifg
    base_unw_paths = original_ifg_paths(pars[cf.IFG_FILE_LIST])

    # dest_paths are tifs that have been geotif converted and multilooked
    dest_paths = get_dest_paths(base_unw_paths, crop, pars, xlks)

    return base_unw_paths, dest_paths, pars


def write_msg(msg):
    """
    write message to log file and screen output
    """
    logging.debug(msg)
    if VERBOSE:
        print msg


if __name__ == "__main__":
    main()
