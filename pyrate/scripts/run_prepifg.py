# -*- coding: utf-8 -*-
import sys
import os
import logging

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
from pyrate.tasks import gamma as gamma_task
import pyrate.ifgconstants as ifc
from pyrate.utils import tparmap

ROI_PAC_HEADER_FILE_EXT = 'rsc'


def main(params=None):
    """
    :param config_file: config file to use. This provides a convenient way to
     use run_prepifg from within the module.
    :return:
    """
    if params:
        base_ifg_paths, dest_paths, _ = run_pyrate.get_ifg_paths()
    else:
        base_ifg_paths, dest_paths, params = run_pyrate.get_ifg_paths()

    LUIGI = params[cf.LUIGI]  # luigi or no luigi
    PROCESSOR = params[cf.PROCESSOR]  # roipac or gamma
    run_pyrate.init_logging(logging.DEBUG)

    usage = 'Usage: python run_prepifg.py <config file>'
    if len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print usage
        return
    raw_config_file = sys.argv[1]
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
        roipac.to_geotiff(header, b, d, nodata=params[cf.NO_DATA_VALUE])
    ifgs = [Ifg(p) for p in dest_base_ifgs]
    prepifg.prepare_ifgs(
        ifgs, crop_opt=crop, xlooks=xlooks, ylooks=ylooks)


def gamma_prepifg(base_ifg_paths, params):
    msg = "running gamma prepifg"
    print msg
    logging.info(msg)
    xlooks, ylooks, crop = run_pyrate.transform_params(params)
    dem_hdr_path = params[cf.DEM_HEADER_FILE]
    DEM_HDR = gamma.parse_dem_header(dem_hdr_path)
    try:
        SLC_DIR = params[cf.SLC_DIR]
    except:
        SLC_DIR = None

    # location of geo_tif's
    dest_base_ifgs = [os.path.join(
        params[cf.OUT_DIR], os.path.basename(q).split('.')[0] + '.tif')
                      for q in base_ifg_paths]

    parallel = params[cf.PARALLEL]
    if parallel:
        print 'running gamma in parallel with {} ' \
              'processes'.format(params[cf.PROCESSES])
        parmap.map(gamma_multiprocessing, base_ifg_paths,
                   DEM_HDR, SLC_DIR, params,
                   processes=params[cf.PROCESSES])
    else:
        for b in base_ifg_paths:
            gamma_multiprocessing(b, DEM_HDR, SLC_DIR, params)

    ifgs = [Ifg(p) for p in dest_base_ifgs]

    exts = prepifg.getAnalysisExtent(crop, ifgs, xlooks, ylooks, userExts=None)
    thresh = params[cf.NO_DATA_AVERAGING_THRESHOLD]
    verbose = False
    if parallel:  # using threadpool due to pickling issue
        tparmap.map(prepifg.prepare_ifg, ifgs,
                    xlooks, ylooks, exts, thresh, crop, verbose,
                    processes=params[cf.PROCESSES])
    else:
        [prepifg.prepare_ifg(i,
               xlooks, ylooks, exts, thresh, crop, verbose) for i in ifgs]


def gamma_multiprocessing(b, DEM_HDR, SLC_DIR, params):
    d = os.path.join(
        params[cf.OUT_DIR], os.path.basename(b).split('.')[0] + '.tif')
    header_paths = gamma_task.get_header_paths(b, slc_dir=SLC_DIR)
    if len(header_paths) != 2:
        raise
    hdrs = [gamma.parse_epoch_header(p) for p in header_paths]
    COMBINED = gamma.combine_headers(hdrs[0], hdrs[1], dem_hdr=DEM_HDR)
    gamma.to_geotiff(COMBINED, b, d,
                     nodata=params[cf.NO_DATA_VALUE])

if __name__ == '__main__':
    main()
