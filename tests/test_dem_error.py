import shutil
import glob
import os
from os.path import join
from scipy.interpolate import griddata
import math
import pytest
import numpy as np
from tests import common
from pyrate.configuration import Configuration, MultiplePaths
from pyrate import prepifg, correct
import pyrate.core.config as cf
import pyrate.core.geometry as geom
from pyrate.core.dem_error import dem_error_calc_wrapper, _calculate_bperp_for_tile
from pyrate.core.shared import Ifg, Geometry, DEM, save_numpy_phase


geometry_path = common.MEXICO_TEST_DIR_GEOMETRY


@pytest.fixture(params=list(range(200)))
def point():
    x, y = np.random.randint(0, 60), np.random.randint(0, 100)
    return x, y


class TestPyRateGammaBperp:

    @classmethod
    def setup_class(cls):
        cls.params = Configuration(common.MEXICO_CONF).__dict__
        # run prepifg
        prepifg.main(cls.params)
        # copy IFGs to temp folder
        correct._copy_mlooked(cls.params)
        # read radar azimuth, range and dem tif files
        rdc_az_file = join(cls.params[cf.OUT_DIR], 'rdc_azimuth.tif')
        geom_az = Geometry(rdc_az_file)
        geom_az.open(readonly=True)
        cls.az = geom_az.geometry_data
        rdc_rg_file = join(cls.params[cf.OUT_DIR], 'rdc_range.tif')
        geom_rg = Geometry(rdc_rg_file)
        geom_rg.open(readonly=True)
        cls.rg = geom_rg.geometry_data
        dem_file = join(cls.params[cf.OUT_DIR], 'dem.tif')
        dem_data = DEM(dem_file, tile=None)
        dem_data.open(readonly=True)
        cls.dem = dem_data.height_data
        # calc bperp using pyrate funcs
        cls.pbperp = cls.pyrate_bperp()

    def gamma_bperp(self, x, y):
        """
        Calculate Bperp for specified pixel from GAMMA out files (interpolation required)
        x0, y0 is the interpolation location in azimuth and range
        """
        # round azimuth and range coordinates to closest step (500 for az, 200 for rg)
        azstep = 500
        rgstep = 200
        az = self.az[x, y]
        rg = self.rg[x, y]

        az1 = azstep * math.floor(az / azstep)
        rg1 = rgstep * math.floor(rg / rgstep)
        az2 = azstep * math.ceil(az / azstep)
        rg2 = rgstep * math.ceil(rg / rgstep)

        # four coordinates for bi-linear interpolation
        teststr1 = str(az1).rjust(6, ' ') + str(rg1).rjust(7, ' ')
        teststr2 = str(az1).rjust(6, ' ') + str(rg2).rjust(7, ' ')
        teststr3 = str(az2).rjust(6, ' ') + str(rg1).rjust(7, ' ')
        teststr4 = str(az2).rjust(6, ' ') + str(rg2).rjust(7, ' ')

        # loop through all corresponding bperp.par files in base.par list
        bperp_files = sorted(list(glob.glob(os.path.join(geometry_path, '*_bperp.par'))))
        bperp_int = np.empty(shape=(len(bperp_files)))
        for i, bperp_file in enumerate(bperp_files):
            # read Mexico city bperp file
            with open(bperp_file, 'r') as f:
                for line in f.readlines():
                    if teststr1 in line:
                        bperp1 = line.split()[7]
                    if teststr2 in line:
                        bperp2 = line.split()[7]
                    if teststr3 in line:
                        bperp3 = line.split()[7]
                    if teststr4 in line:
                        bperp4 = line.split()[7]

            # setup numpy array for bi-linear interpolation
            n = np.array([(az1, rg1, bperp1),
                          (az1, rg2, bperp2),
                          (az2, rg1, bperp3),
                          (az2, rg2, bperp4)])
            # interpolate using scipy function "griddata"
            bperp_int[i] = griddata(n[:, 0:2], n[:, 2], [(az, rg)], method='linear')

        return bperp_int

    @classmethod
    def pyrate_bperp(cls):
        """
        Calculate Bperp image for each ifg using PyRate functions
        """
        multi_paths = cls.params[cf.INTERFEROGRAM_FILES]
        tmp_paths = [ifg_path.tmp_sampled_path for ifg_path in multi_paths]
        # keep only ifg files in path list (i.e. remove coherence and dem files)
        ifg_paths = [item for item in tmp_paths if 'ifg.tif' in item]
        # read and open the first IFG in list
        ifg0_path = ifg_paths[0]
        ifg0 = Ifg(ifg0_path)
        ifg0.open(readonly=True)
        # size of ifg dataset
        # calculate per-pixel lon/lat
        lon, lat = geom.get_lonlat_coords(ifg0)
        bperp = _calculate_bperp_for_tile(ifg_paths, cls.az, cls.rg, lat, lon, cls.dem, tile=None)[0]
        return np.moveaxis(bperp, (0, 1, 2), (2, 0, 1))

    @classmethod
    def teardown_class(cls):
        shutil.rmtree(cls.params[cf.OUT_DIR], ignore_errors=True)

    def test_pyrate_bperp_matches_gamma_bperp(self, point):
        x, y = point
        az = self.az[x, y]
        rg = self.rg[x, y]

        if az < 0 or rg < 0:
            pytest.skip('skipped due to -ve az or rg')

        res = self.pbperp[x, y, :]
        exp = self.gamma_bperp(* point)
        np.testing.assert_array_almost_equal(exp, res, 2)  # max difference < 1cm


def test_calc_dem_errors():
    pass


class TestDEMErrorFilesReusedFromDisc:

    @classmethod
    def setup_class(cls):
        cls.conf = common.MEXICO_CONF
        cls.params = Configuration(cls.conf).__dict__
        prepifg.main(cls.params)
        cls.params = Configuration(cls.conf).__dict__
        multi_paths = cls.params[cf.INTERFEROGRAM_FILES]
        cls.ifg_paths = [p.tmp_sampled_path for p in multi_paths]

    @classmethod
    def teardown_class(cls):
        shutil.rmtree(cls.params[cf.OUT_DIR])

    def test_dem_error_used_from_disc_on_rerun(self):
        correct._update_params_with_tiles(self.params)
        times_written = self.__run_once()
        assert len(times_written) == len(self.ifg_paths)
        times_written_1 = self.__run_once()
        np.testing.assert_array_equal(times_written_1, times_written)

    def __run_once(self):
        dem_files = [MultiplePaths.dem_error_path(i, self.params) for i in self.ifg_paths]
        correct._copy_mlooked(self.params)
        correct._update_params_with_tiles(self.params)
        correct._create_ifg_dict(self.params)
        save_numpy_phase(self.ifg_paths, self.params)
        dem_error_calc_wrapper(self.params)
        assert all(m.exists() for m in dem_files)
        return [os.stat(o).st_mtime for o in dem_files]
