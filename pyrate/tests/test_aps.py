from pyrate.tasks.roipac import _DoConvertToGeotiffRoipac

__author__ = 'sudipta'

import unittest
import tempfile
import os
import shutil
import gdal
from pyrate import remove_aps_delay as aps
from pyrate.tests import common
from pyrate import config as cf
from pyrate.scripts import run_pyrate, run_prepifg
from pyrate import ifgconstants as ifc


class MetaDataTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tif_dir = tempfile.mkdtemp()
        cls.test_conf = common.SYDNEY_TEST_CONF

        # change the required params
        cls.params = cf.get_config_params(cls.test_conf)
        cls.params[cf.OBS_DIR] = common.SYD_TEST_GAMMA
        cls.params[cf.PROCESSOR] = 1  # gamma
        file_list = cf.parse_namelist(os.path.join(common.SYD_TEST_GAMMA,
                                                   'ifms_17'))
        fd, cls.params[cf.IFG_FILE_LIST] = tempfile.mkstemp(suffix='.conf',
                                                             dir=cls.tif_dir)
        os.close(fd)
        # write a short filelist with only 3 gamma unws
        with open(cls.params[cf.IFG_FILE_LIST], 'w') as fp:
            for f in file_list[:3]:
                fp.write(os.path.join(common.SYD_TEST_GAMMA, f) + '\n')
        cls.params[cf.OUT_DIR] = cls.tif_dir
        cls.params[cf.PARALLEL] = 0
        cls.params[cf.REF_EST_METHOD] = 1
        cls.params[cf.DEM_FILE] = common.SYD_TEST_DEM_GAMMA
        # base_unw_paths need to be geotiffed and multilooked by run_prepifg
        base_unw_paths = run_pyrate.original_ifg_paths(
            cls.params[cf.IFG_FILE_LIST])
        base_unw_paths.append(common.SYD_TEST_DEM_GAMMA)

        xlks, ylks, crop = run_pyrate.transform_params(cls.params)

        # dest_paths are tifs that have been geotif converted and multilooked
        run_prepifg.gamma_prepifg(base_unw_paths, cls.params)
        base_unw_paths.pop()  # removed dem as we don't want it in ifgs list

        dest_paths = run_pyrate.get_dest_paths(
            base_unw_paths, crop, cls.params, xlks)
        cls.ifgs = common.sydney_data_setup(datafiles=dest_paths)
        aps.remove_aps_delay(cls.ifgs, cls.params)

    @classmethod
    def tearDownClass(cls):
        for i in cls.ifgs:
            i.close()
        shutil.rmtree(cls.tif_dir)

    def test_metadata_was_copied(self):
        for i in self.ifgs:
            md = i.meta_data
            self.assertIn(ifc.PYRATE_APS_ERROR, md.keys())
            self.assertIn(aps.APS_STATUS, md.values())

    def test_meta_data_was_written(self):
        for i in self.ifgs:
            i.close()  # until file is closed metadata is not written
            md = i.meta_data
            ds = gdal.Open(i.data_path)
            md_w = ds.GetMetadata()
            self.assertDictEqual(md, md_w)
            ds = None

    def test_dems_present(self):
        # geotiffed dem
        os.path.exists(os.path.join(self.params[cf.OUT_DIR],
                       os.path.splitext(common.SYD_TEST_DEM_GAMMA)[0]
                                    + '.tif'))

        # multilooked dem
        os.path.exists(os.path.join(self.params[cf.OUT_DIR],
                       os.path.splitext(common.SYD_TEST_DEM_GAMMA)[0]
                        + '_{looks}rlks_{crop}cr.tif'.format(
                           looks=self.params[cf.IFG_LKSX],
                           crop=self.params[cf.IFG_CROP_OPT])))


if __name__ == '__main__':
    unittest.main()
