from __future__ import print_function
import glob
import os
import re
import shutil
import subprocess
import tempfile
import unittest

from pyrate import config as cf
from pyrate.scripts import run_prepifg
from tests import common


class MPITests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.gamma_pypar_dir = tempfile.mkdtemp()
        cls.gamma_serial_dir = tempfile.mkdtemp()
        cls.conf_file = os.path.join(common.SYD_TEST_GAMMA,
                                     'pyrate_gamma.conf')
        cls.params = cf.get_config_params(cls.conf_file)
        cls.params[cf.OBS_DIR] = common.SYD_TEST_GAMMA
        cls.params[cf.IFG_FILE_LIST] = os.path.join(common.SYD_TEST_GAMMA,
                                                    'ifms_17')
        cls.params[cf.OUT_DIR] = cls.gamma_pypar_dir
        cls.temp_conf_file = tempfile.mktemp(suffix='.conf',
                                             dir=cls.gamma_pypar_dir)
        cf.write_config_file(cls.params, cls.temp_conf_file)

        # run mpi gamma
        cmd = "mpirun -np 4 python pyrate/nci/run_prepifg_pypar.py " + \
              cls.temp_conf_file
        subprocess.check_call(cmd, shell=True)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.gamma_pypar_dir)
        shutil.rmtree(cls.gamma_serial_dir)

    def test_mpi_vs_serial_gamma(self):

        self.params[cf.OUT_DIR] = self.gamma_serial_dir
        self.params[cf.PARALLEL] = False

        # run serial gamma
        run_prepifg.main(self.params)
        gamma_PTN = re.compile(r'\d{8}')
        mpi_tifs = []
        for i in glob.glob(os.path.join(self.gamma_pypar_dir,
                                        "*.tif")):
            if len(gamma_PTN.findall(i)) == 2:
                mpi_tifs.append(i)

        serial_tifs = []
        for i in glob.glob(os.path.join(self.gamma_serial_dir,
                                        "*.tif")):
            if len(gamma_PTN.findall(i)) == 2:
                serial_tifs.append(i)

        self.assertEqual(len(mpi_tifs), len(serial_tifs))
        for m, s in zip(mpi_tifs, serial_tifs):
            self.assertEqual(os.path.basename(m), os.path.basename(s))
        # 17 geotifs, and 17 mlooked tifs
        self.assertEqual(len(mpi_tifs), 34)

    def test_mpi_log_file_created(self):
        self.assertTrue(os.path.exists(
            os.path.join(self.gamma_pypar_dir, 'mpi.log')))


if __name__ == '__main__':
    unittest.main()
