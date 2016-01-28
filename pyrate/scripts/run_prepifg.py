# -*- coding: utf-8 -*-
import sys
import luigi
import os
from pyrate.tasks.utils import pythonifyConfig
from pyrate.tasks.prepifg import PrepareInterferograms
from pyrate import prepifg
from pyrate import config as cf
from pyrate.shared import Ifg
from pyrate.scripts import run_pyrate
from pyrate import roipac
from pyrate import gamma
import pyrate.ifgconstants as ifc

ROI_PAC_HEADER_FILE_EXT = 'rsc'


def main():
    """
    :param config_file: config file to use. This provides a convenient way to
     use run_prepifg from within the module.
    :return:
    """
    base_ifg_paths, dest_paths, params = run_pyrate.get_ifg_paths()
    LUIGI = params[cf.LUIGI]  # luigi or no luigi
    PROCESSOR = params[cf.PROCESSOR]  # roipac or gamma

    usage = 'Usage: python run_prepifg.py <config file>'
    if len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print usage
        return
    raw_config_file = sys.argv[1]
    if LUIGI:
        print "running luigi prepifg"
        luigi.configuration.LuigiConfigParser.add_config_path(
            pythonifyConfig(raw_config_file))
        luigi.build([PrepareInterferograms()], local_scheduler=True)
    else:
        print "running serial prepifg"
        xlooks, ylooks, crop = run_pyrate.transform_params(params)

        if PROCESSOR == 0:
            print 'running roipac prepifg'
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
        else:
            print 'running gamma prepifg'
            header_files = [os.path.join(
                params[cf.OBS_DIR],
                os.path.basename(q).split('-')[0] + '_slc.par')
                            for q in base_ifg_paths]

            hdr_paths = [os.path.join(cf.PYRATEPATH, f) for f in header_files]
            hdrs = [gamma.parse_epoch_header(p) for p in hdr_paths]
            dem_hdr_path = params[cf.DEM_HEADER_FILE]
            DEM_HDR = gamma.parse_dem_header(dem_hdr_path)
            COMBINED = gamma.combine_headers(hdrs[0], hdrs[1], dem_hdr=DEM_HDR)

            dest_base_ifgs = [os.path.join(
                params[cf.OUT_DIR], os.path.basename(q).split('.')[0] + '.tif')
                          for q in base_ifg_paths]
            for b, d in zip(base_ifg_paths, dest_base_ifgs):
                gamma.to_geotiff(COMBINED, b, d,
                                 nodata=params[cf.NO_DATA_VALUE])

            ifgs = [Ifg(p) for p in dest_base_ifgs]
            prepifg.prepare_ifgs(
                ifgs, crop_opt=crop, xlooks=xlooks, ylooks=ylooks)


if __name__ == '__main__':
    main()
