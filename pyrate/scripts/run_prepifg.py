# -*- coding: utf-8 -*-
import sys
import os
import logging
import glob
import luigi
import parmap

from pyrate.tasks.utils import pythonifyConfig
from pyrate.tasks.prepifg import PrepareInterferograms
from pyrate import prepifg
from pyrate import config as cf
from pyrate.shared import Ifg
from pyrate.scripts import run_pyrate
from pyrate import roipac
from pyrate import gamma
from pyrate.shared import write_geotiff, mkdir_p
from pyrate.tasks import gamma as gamma_task
import pyrate.ifgconstants as ifc
from pyrate.utils import tparmap
from operator import itemgetter
from osgeo import gdal
from pyrate.prepifg import mlooked_path

ROI_PAC_HEADER_FILE_EXT = 'rsc'


def main(params=None):
    """
    :param config_file: config file to use. This provides a convenient way to
     use run_prepifg from within the module.
    :return:
    """
    usage = 'Usage: python run_prepifg.py <config file>'
    if params:
        base_ifg_paths = glob.glob(os.path.join(params[cf.OBS_DIR], '*.unw'))
        LUIGI = params[cf.LUIGI]  # luigi or no luigi
        if LUIGI:
            raise cf.ConfigException('params can not be provided with luigi')
    else:  # if params not provided read from config file
        if (not params) and ((len(sys.argv) == 1)
                or (sys.argv[1] == '-h' or sys.argv[1] == '--help')):
            print usage
            return
        base_ifg_paths, _, params = run_pyrate.get_ifg_paths()
        '''
        TODO: looks like base_igf_paths are ordered according to ifg list
        Sarah says this probably won't be a problem because they only removedata from ifg list to get pyrate to work
        (i.e. things won't be reordered and the original gamma generated list is ordered)
        this may not affect the important pyrate stuff anyway, but might affect gen_thumbs.py
        going to assume base_ifg_paths is ordered correcly
        '''
        raw_config_file = sys.argv[1]

    LUIGI = params[cf.LUIGI]  # luigi or no luigi
    PROCESSOR = params[cf.PROCESSOR]  # roipac or gamma
    run_pyrate.init_logging(logging.DEBUG)

    if LUIGI:
        msg = "running luigi prepifg"
        print msg
        logging.info(msg)
        luigi.configuration.LuigiConfigParser.add_config_path(
            pythonifyConfig(raw_config_file))
        luigi.build([PrepareInterferograms()], local_scheduler=True)
    else:
        msg = "running serial prepifg"
        print msg
        logging.info(msg)
        if PROCESSOR == 0:
            roipac_prepifg(base_ifg_paths, params)
        else:
            gamma_prepifg(base_ifg_paths, params)


def roipac_prepifg(base_ifg_paths, params):
    msg = "running roipac prepifg"
    print msg
    logging.info(msg)
    xlooks, ylooks, crop = run_pyrate.transform_params(params)
    dem_file = os.path.join(params[cf.ROIPAC_RESOURCE_HEADER])
    projection = roipac.parse_header(dem_file)[ifc.PYRATE_DATUM]
    dest_base_ifgs = [os.path.join(
        params[cf.OUT_DIR], os.path.basename(q).split('.')[0] + '.tif')
                      for q in base_ifg_paths]
    for b, d in zip(base_ifg_paths, dest_base_ifgs):
        header_file = "%s.%s" % (b, ROI_PAC_HEADER_FILE_EXT)
        header = roipac.parse_header(header_file)

        # DEM already has DATUM, so get it from dem if not in header
        if ifc.PYRATE_DATUM not in header:
            header[ifc.PYRATE_DATUM] = projection
        write_geotiff(header, b, d, nodata=params[cf.NO_DATA_VALUE])
    prepifg.prepare_ifgs(
        dest_base_ifgs, crop_opt=crop, xlooks=xlooks, ylooks=ylooks)


def gamma_prepifg(base_unw_paths, params):
    msg = "running gamma prepifg"
    print msg
    logging.info(msg)
    # location of geo_tif's
    dest_base_ifgs = [os.path.join(
        params[cf.OUT_DIR], os.path.basename(q).split('.')[0] + '.tif')
                      for q in base_unw_paths]

    parallel = params[cf.PARALLEL]
    if parallel:
        print 'running gamma in parallel with {} ' \
              'processes'.format(params[cf.PROCESSES])
        import multiprocessing
        print 'found', multiprocessing.cpu_count(), 'CPUs'
        parmap.map(gamma_multiprocessing, base_unw_paths,
                   params, processes=params[cf.PROCESSES])
    else:
        for b in base_unw_paths:
            gamma_multiprocessing(b, params)
    ifgs = [Ifg(p) for p in dest_base_ifgs]
    xlooks, ylooks, crop = run_pyrate.transform_params(params)
    exts = prepifg.getAnalysisExtent(crop, ifgs, xlooks, ylooks, userExts=None)
    thresh = params[cf.NO_DATA_AVERAGING_THRESHOLD]
    if parallel:
        parmap.map(prepifg.prepare_ifg, dest_base_ifgs,
                   xlooks, ylooks, exts, thresh, crop,
                   processes=params[cf.PROCESSES])
    else:
        [prepifg.prepare_ifg(i, xlooks, ylooks, exts,
                             thresh, crop) for i in dest_base_ifgs]

    # easier just to set sequence information here
    # open all generated ifgs and set sequence metadata
    # assuming all files have finished being written to (i.e. all above tasks are complete)
    # create shallow copy just incase
    dest_base_ifgs_scp = list(dest_base_ifgs)
    mlooked_paths = [mlooked_path(dest_base_ifg, xlooks, crop) for dest_base_ifg in dest_base_ifgs_scp]
    order_fn_parts = [(os.path.split(dest_base_ifg)[1]).split('_')[0] for dest_base_ifg in dest_base_ifgs]
    zps = zip(order_fn_parts, dest_base_ifgs, mlooked_paths)
    zps_sorted = sorted(zps, key=itemgetter(0))
    for it, zp in enumerate(zps_sorted):
        _, dest_base_ifg, mled = zp
        # -----------------------------------------------
        g_ds1 = gdal.Open(dest_base_ifg)
        g_ds1.SetMetadataItem('PR_SEQ_POS', str(it))
        # -----------------------------------------------
        g_ds2 = gdal.Open(mled)
        g_ds2.SetMetadataItem('PR_SEQ_POS', str(it))

def gamma_multiprocessing(b, params):
    dem_hdr_path = params[cf.DEM_HEADER_FILE]
    DEM_HDR = gamma.parse_dem_header(dem_hdr_path)
    try:
        SLC_DIR = params[cf.SLC_DIR]
    except:
        SLC_DIR = None

    mkdir_p(params[cf.OUT_DIR])
    d = os.path.join(
        params[cf.OUT_DIR], os.path.basename(b).split('.')[0] + '.tif')

    header_paths = gamma_task.get_header_paths(b, slc_dir=SLC_DIR)
    if len(header_paths) != 2:
        raise
    hdrs = [gamma.parse_epoch_header(p) for p in header_paths]
    COMBINED = gamma.combine_headers(hdrs[0], hdrs[1], dem_hdr=DEM_HDR)
    '''
    print type(COMBINED)
    print COMBINED
    while True: pass
    '''
    COMBINED['PR_TYPE'] = 'ifg_1'   # non-cropped, non-multilooked ifg
    write_geotiff(COMBINED, b, d, nodata=params[cf.NO_DATA_VALUE])

if __name__ == '__main__':
    main()
