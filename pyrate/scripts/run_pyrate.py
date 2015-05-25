"""
Main workflow script for PyRate

Created on 17/09/2012
@author: Ben Davies, NCI
"""

from __future__ import print_function

import os,sys,shutil,logging,datetime

import pyrate.mst as mst
import pyrate.prepifg as prepifg
import pyrate.algorithm as algorithm
import pyrate.refpixel as refpixel
import pyrate.orbital as orbital
import pyrate.linrate as linrate
import pyrate.timeseries as timeseries
import pyrate.config as cf
from pyrate.shared import Ifg

# constants for metadata flags
META_UNITS = 'PHASE_UNITS'
MILLIMETRES = 'MILLIMETRES'
META_ORBITAL = 'ORBITAL_ERROR'
META_REMOVED = 'REMOVED'


def process_ifgs(ifg_paths, params):
    """
    Top level function to perform correction steps on given ifgs
    ifgs: sequence of paths to interferograms (NB: changes are saved into ifgs)
    params: dict of run config params
    """
    ifgs = [Ifg(p) for p in ifg_paths]

    for i in ifgs:
        if not i.is_open:
            i.open(readonly=False)
        convert_wavelength(i)

    remove_orbital_error(ifgs, params)

    mst_grid = mst.mst_matrix_ifgs_only(ifgs)
    refpx, refpy = find_reference_pixel(ifgs, params)

    # TODO: VCM code. which part gets called? cvd?
    #vcm = None
    #calculate_linear_rate(ifgs, params, vcm, mst=None)

    pthresh = params[cf.TIME_SERIES_PTHRESH]
    #calculate_time_series(ifgs, pthresh, mst=None)  # TODO: check is correct MST

    # TODO: outputs?

    # final cleanup
    while ifgs:
        i = ifgs.pop()
        i.write_phase()
        i.dataset.FlushCache()
        i = None # force close    TODO: may need to implement close()

    logging.debug('PyRate run completed\n')


def convert_wavelength(ifg):
    if ifg.dataset.GetMetadataItem(META_UNITS) == MILLIMETRES:
        msg = '%s: ignored as previous wavelength conversion detected'
        logging.debug(msg % ifg.data_path)
        return

    ifg.data = algorithm.wavelength_radians_to_mm(ifg.phase_data, ifg.wavelength)
    ifg.dataset.SetMetadataItem(META_UNITS, MILLIMETRES)
    msg = '%s: converted wavelength to millimetres'
    logging.debug(msg % ifg.data_path)


def remove_orbital_error(ifgs, params):
    if not params[cf.ORBITAL_FIT]:
        logging.debug('Orbital correction not required.')
        return

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
    if count < len(flags) and count > 0:
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
    refx = params.get(params[cf.REFX])
    if refx > ifgs[0].ncols - 1:
        raise ValueError("Invalid reference pixel X coordinate: %s" % refx)

    refy = params.get(params[cf.REFY])
    if refy > ifgs[0].nrows - 1:
        raise ValueError("Invalid reference pixel Y coordinate: %s" % refy)

    if refx >= 0 and refy >= 0:
        msg = 'Reusing config file reference pixel (%s, %s)'
        logging.debug(msg % (refx, refy))
        return (refy, refx) # reuse preset ref pixel

    refy, refx = refpixel.ref_pixel(ifgs,
                                    params[cf.REFNX],
                                    params[cf.REFNY],
                                    params[cf.REF_CHIP_SIZE],
                                    params[cf.REF_MIN_FRAC])

    logging.debug('Reference pixel coordinate: (%s, %s)' % (refx, refy))
    return refx, refy


def calculate_linear_rate(ifgs, params, vcm, mst=None):
    logging.debug('Calculating linear rate')

    pthr = params[cf.LR_PTHRESH]
    nsig = params[cf.LR_NSIG]
    maxsig = params[cf.LR_MAXSIG]

    # TODO: do these need to be checked?
    res = linrate.linear_rate(ifgs, vcm, pthr, nsig, maxsig, mst)
    for r in res:
        if r is None:
            raise ValueError('TODO: bad value')

    rate, error, samples = res

    logging.debug('Linear rate calculated')
    return rate, error, samples


def calculate_time_series(ifgs, pthresh, mst):
    logging.debug('Calculating time series')
    res = timeseries.time_series(ifgs, pthresh, mst)

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

    if crop is None or crop == prepifg.ALREADY_SAME_SIZE:
        return False
    return True


def original_ifg_paths(ifglist_path):
    """
    Returns sequence of paths to files in given ifglist file.
    """

    basedir = os.path.dirname(ifglist_path)
    ifglist = cf.parse_namelist(ifglist_path)
    return [os.path.join(basedir, p) for p in ifglist]


def working_ifg_paths(src_paths, xlooks, ylooks, cropping):
    """
    Filter. Returns paths to ifgs to process (eg. checks for mlooked ifgs)
    """
    if warp_required(xlooks, ylooks, cropping):
        mlooked_paths = [prepifg.mlooked_path(p, xlooks) for p in src_paths]

        if not all([os.path.exists(p) for p in mlooked_paths]):
            msg = 'Multilooked ifgs do not exist (check config settings?)'
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


# TODO: write to alternate file if log exists
def init_logging(level):
    t = datetime.datetime.now()
    path = 'pyrate_%s_%02d_%02d.log' % (t.year, t.month, t.day)
    fmt = '%(asctime)s %(message)s'
    datefmt = '%d/%m/%Y %I:%M:%S %p'
    logging.basicConfig(filename=path, format=fmt, datefmt=datefmt, level=level)
    logging.debug('Log started')


# TODO: expand CLI documentation
# TODO: ensure clean exception handling
# TODO: add parameter error checking: induce fail fast before number crunching
def main():
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [config-file]\nRuns PyRate workflow.')
    options, args = parser.parse_args()

    init_logging(logging.DEBUG)

    try:
        cfg_path = args[0] if args else 'pyrate.conf'
        pars = cf.get_config_params(cfg_path)
    except IOError as err:
        emsg = 'Config file error: %s "%s"' % (err.strerror, err.filename)
        logging.debug(emsg)
        print(emsg)
        sys.exit(err.errno)

    # FIXME: make output ifgs here, or in process_ifgs() ?
    xlks, ylks, crop = transform_params(pars)
    base_ifg_paths = original_ifg_paths(pars[cf.IFG_FILE_LIST])
    working_paths = working_ifg_paths(base_ifg_paths, xlks, ylks, crop)
    dest_paths = dest_ifg_paths(working_paths, pars[cf.OUT_DIR])

    if not os.path.exists(pars[cf.OUT_DIR]):
        os.makedirs(pars[cf.OUT_DIR])

    # process copies of source data
    for wp, dp in zip(working_paths, dest_paths):
        shutil.copy(wp, dp)

    process_ifgs(dest_paths, pars)

if __name__ == "__main__":
    main()
