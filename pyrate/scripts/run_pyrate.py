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
from pyrate.shared import Ifg
from pyrate import vcm as vcm_module
from pyrate import matlab_mst_kruskal as matlab_mst
from pyrate import reference_phase_estimation as rpe
from pyrate import algorithm
from pyrate import ifgconstants as ifc

# constants for metadata flags
META_UNITS = 'PHASE_UNITS'
MILLIMETRES = 'MILLIMETRES'
META_ORBITAL = 'ORBITAL_ERROR'
META_REMOVED = 'REMOVED'
pars = None
PYRATEPATH = cf.PYRATEPATH


def process_ifgs(ifg_paths_or_instance, params):
    """
    Top level function to perform correction steps on given ifgs
    ifgs: sequence of paths to interferograms (NB: changes are saved into ifgs)
    params: dict of run config params
    """
    nan_conversion = int(params[cf.NAN_CONVERSION])
    # ifg_instance, nan_conversion=True

    if isinstance(ifg_paths_or_instance, list):
        ifgs = [Ifg(p) for p in ifg_paths_or_instance]

        for i in ifgs:
            if not i.is_open:
                i.open(readonly=False)
            if nan_conversion:  # nan conversion happens here in networkx mst
                i.convert_to_nans()
            if not i.mm_converted:
                i.convert_to_mm()
                i.write_modified_phase()
        mst_grid = mst.mst_matrix_ifg_indices_as_boolean_array(ifgs)
    else:
        assert isinstance(ifg_paths_or_instance, matlab_mst.IfgListPyRate)
        ifgs = ifg_paths_or_instance.ifgs
        for i in ifgs:
            if not i.mm_converted:
                i.convert_to_mm()
                i.write_modified_phase()
        ifg_instance_updated, epoch_list = \
            matlab_mst.get_nml(ifg_paths_or_instance,
                               nan_conversion=nan_conversion)
        mst_grid = matlab_mst.matlab_mst_boolean_array(ifg_instance_updated)

    # write mst output to a file
    mst_mat_binary_file = os.path.join(
        PYRATEPATH, params[cf.OUT_DIR], 'mst_mat')
    np.save(file=mst_mat_binary_file, arr=mst_grid)

    # Estimate reference pixel location
    refpx, refpy = find_reference_pixel(ifgs, params)

    # remove APS delay here

    # Estimate and remove orbit errors
    if params[cf.ORBITAL_FIT] != 0:
        remove_orbital_error(ifgs, params)

    _, ifgs = rpe.estimate_ref_phase(ifgs, params, refpx, refpy)

    # Calculate interferogram noise
    # TODO: assign maxvar to ifg metadata (and geotiff)?
    maxvar = [vcm_module.cvd(i)[0] for i in ifgs]
    
    # Calculate temporal variance-covariance matrix
    vcmt = vcm_module.get_vcmt(ifgs, maxvar)

    p = os.path.join(params[cf.OUT_DIR], ifgs[0].data_path)
    assert os.path.exists(p) == True

    ds = gdal.Open(p)
    md = ds.GetMetadata()  # get metadata for writing on output tifs
    epochlist = algorithm.get_epochs(ifgs)

    if params[cf.TIME_SERIES_CAL] != 0:

        # Calculate time series
        tsincr, tscum, tsvel = calculate_time_series(
            ifgs, params, vcmt=vcmt, mst=mst_grid)

        for i in range(len(tsincr[0, 0, :])):
            md[ifc.PYRATE_DATE] = epochlist.dates[i+1]
            data = tsincr[:, :, i]
            dest = os.path.join(
                PYRATEPATH, params[cf.OUT_DIR],
                "tsincr_" + str(epochlist.dates[i+1]) + ".tif")
            timeseries.write_geotiff_output(md, data, dest, np.nan)

            data = tscum[:, :, i]
            dest = os.path.join(
                PYRATEPATH, params[cf.OUT_DIR],
                "tscuml_" + str(epochlist.dates[i+1]) + ".tif")
            timeseries.write_geotiff_output(md, data, dest, np.nan)

            data = tsvel[:, :, i]
            dest = os.path.join(
                PYRATEPATH, params[cf.OUT_DIR],
                "tsvel_" + str(epochlist.dates[i+1]) + ".tif")
            timeseries.write_geotiff_output(md, data, dest, np.nan)

    # Calculate linear rate map
    rate, error, samples = calculate_linear_rate(
                   ifgs, params, vcmt, mst=mst_grid)
    md[ifc.PYRATE_DATE] = epochlist.dates
    dest = os.path.join(PYRATEPATH, params[cf.OUT_DIR], "linrate.tif")
    timeseries.write_geotiff_output(md, rate, dest, np.nan)
    dest = os.path.join(PYRATEPATH, params[cf.OUT_DIR], "linerror.tif")
    timeseries.write_geotiff_output(md, rate, dest, np.nan)

    # final cleanup, SB: Why do we need this?
    # while ifgs:
    #     i = ifgs.pop()
    #     i.write_phase()
    #     i.dataset.FlushCache()
    #     i = None  # force close    TODO: may need to implement close()

    logging.debug('PyRate run completed\n')


def remove_orbital_error(ifgs, params):
    if not params[cf.ORBITAL_FIT]:
        logging.debug('Orbital correction not required.')
        return
   
    logging.debug('Calculating orbital correction')

    # perform some general error/sanity checks
    flags = [i.dataset.GetMetadataItem(META_ORBITAL) for i in ifgs]

    if all(flags):
        msg = 'Skipped orbital correction, ifgs already corrected'
        logging.debug(msg)
        return
    else:
        check_orbital_ifgs(ifgs, flags)

    mlooked = None
    if params[cf.ORBITAL_FIT_LOOKS_X] > 1 or params[cf.ORBITAL_FIT_LOOKS_Y] > 1:
        # resampling here to use all prior corrections to orig data
        # TODO: avoid writing mlooked to disk by using mock ifgs/in mem arrays?
        mlooked = prepifg.prepare_ifgs(ifgs,
                                    crop_opt=prepifg.ALREADY_SAME_SIZE,
                                    xlooks=params[cf.ORBITAL_FIT_LOOKS_X],
                                    ylooks=params[cf.ORBITAL_FIT_LOOKS_Y])

    orbital.orbital_correction(ifgs,
                            degree=params[cf.ORBITAL_FIT_DEGREE],
                            method=params[cf.ORBITAL_FIT_METHOD],
                            mlooked=mlooked)

    # mlooked layers discarded as not used elsewhere
    if mlooked:
        for path in [m.data_path for m in mlooked]:
            msg = '%s: deleted (multilooked orbital correction file)'
            logging.debug(msg % path)
            os.remove(path)

    for i in ifgs:
        i.dataset.SetMetadataItem(META_ORBITAL, META_REMOVED)
        logging.debug('%s: orbital error removed' % i.data_path)


def check_orbital_ifgs(ifgs, flags):
    count = sum([f == META_REMOVED for f in flags])
    if (count < len(flags)) and (count > 0):
        msg = 'Detected mix of corrected and uncorrected orbital error in ifgs'
        logging.debug(msg)

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
        logging.debug('Calculating reference pixel')
        refy, refx = refpixel.ref_pixel(ifgs, params[cf.REFNX],
            params[cf.REFNY], params[cf.REF_CHIP_SIZE], params[cf.REF_MIN_FRAC])
    else:
        msg = 'Reusing config file reference pixel (%s, %s)'
        logging.debug(msg % (refx, refy))
        # reuse preset ref pixel

    logging.debug('Reference pixel coordinate: (%s, %s)' % (refx, refy))
    return refx, refy


def calculate_linear_rate(ifgs, params, vcmt, mst=None):
    logging.debug('Calculating linear rate')

    pthr = params[cf.LR_PTHRESH]
    nsig = params[cf.LR_NSIG]
    maxsig = params[cf.LR_MAXSIG]
    # MULTIPROCESSING parameters
    parallel = params[cf.PARALLEL]
    processes = params[cf.PROCESSES]

    # TODO: do these need to be checked?
    res = linrate.linear_rate(ifgs, vcmt, pthr, nsig, maxsig, mst,
                              parallel=parallel, processes=processes)
    for r in res:
        if r is None:
            raise ValueError('TODO: bad value')

    rate, error, samples = res

    logging.debug('Linear rate calculated')
    return rate, error, samples


def calculate_time_series(ifgs, params, vcmt, mst):
    logging.debug('Calculating time series')
    ts_pthr = params[cf.TIME_SERIES_PTHRESH]        
    res = timeseries.time_series(ifgs, ts_pthr, params, vcmt, mst)
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

    ifg_instance = matlab_mst.IfgListPyRate(datafiles=dest_paths)

    if int(pars[cf.NETWORKX_OR_MATLAB_FLAG]):
        msg = 'Running networkx mst'
        print msg
        logging.debug(msg)
        process_ifgs(dest_paths, pars)
    else:
        msg = 'Running matlab mst'
        print msg
        logging.debug(msg)
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
        print emsg
        sys.exit(err.errno)
    ifgListFile = options.ifglist or pars.get(cf.IFG_FILE_LIST)
    if ifgListFile is None:
        emsg = 'Error {code}: Interferogram list file name not provided ' \
               'or does not exist'.format(code=2)
        logging.debug(emsg)
        raise IOError(2, emsg)
    xlks, ylks, crop = transform_params(pars)
    base_unw_paths = original_ifg_paths(pars[cf.IFG_FILE_LIST])
    dest_paths = get_dest_paths(base_unw_paths, crop, pars, xlks)

    return base_unw_paths, dest_paths, pars


if __name__ == "__main__":
    main()
