#   This Python module is part of the PyRate software package.
#
#   Copyright 2017 Geoscience Australia
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
"""
This Python module contains tests for the shared.py PyRate module.
"""

import os
import shutil
import sys
import tempfile
import unittest
from itertools import product
from numpy import isnan, where, nan
from os.path import join, basename, exists
from stat import S_IRGRP, S_IWGRP, S_IWOTH, S_IROTH, S_IRUSR, S_IWUSR

import numpy as np
from numpy.testing import assert_array_equal
from osgeo import gdal
from osgeo.gdal import Open, Dataset, UseExceptions

from tests.common import SYD_TEST_TIF, SYD_TEST_DEM_TIF, TEMPDIR
from pyrate import config as cf
from pyrate import gamma
from pyrate import ifgconstants as ifc
from pyrate import prepifg
from pyrate import shared
from pyrate.scripts import run_prepifg
from pyrate.shared import Ifg, DEM, RasterException
from pyrate.shared import cell_size, utm_zone

from tests import common

UseExceptions()

if not exists(SYD_TEST_TIF):
    sys.exit("ERROR: Missing sydney_test data for unit tests\n")


class IfgTests(unittest.TestCase):
    """Unit tests for the Ifg/interferogram class."""

    def setUp(self):
        self.ifg = Ifg(join(SYD_TEST_TIF, 'geo_060619-061002_unw.tif'))
        self.ifg.open()
        self.ifg.nodata_value = 0

    def test_headers_as_attr(self):
        for a in ['ncols', 'nrows', 'x_first', 'x_step',
                  'y_first', 'y_step', 'wavelength', 'master', 'slave']:
            self.assertTrue(getattr(self.ifg, a) is not None)

    def test_convert_to_nans(self):
        self.ifg.convert_to_nans()
        self.assertTrue(self.ifg.nan_converted)

    def test_xylast(self):
        # ensure the X|Y_LAST header element has been created
        self.assertAlmostEqual(self.ifg.x_last, 150.9491667)
        self.assertAlmostEqual(self.ifg.y_last, -34.23)

    def test_num_cells(self):
        # test cell size from header elements
        data = self.ifg.phase_band.ReadAsArray()
        ys, xs = data.shape
        exp_ncells = ys * xs
        self.assertEqual(exp_ncells, self.ifg.num_cells)

    def test_shape(self):
        self.assertEqual(self.ifg.shape, self.ifg.phase_data.shape)

    def test_nan_count(self):
        num_nan = 0
        for row in self.ifg.phase_data:
            for v in row:
                if isnan(v):
                    num_nan += 1
        if self.ifg.nan_converted:
            self.assertEqual(num_nan, self.ifg.nan_count)
        else:
            self.assertEqual(num_nan, 0)

    def test_phase_band(self):
        data = self.ifg.phase_band.ReadAsArray()
        self.assertEqual(data.shape, (72, 47) )

    def test_nan_fraction(self):
        # NB: source data lacks 0 -> NaN conversion
        data = self.ifg.phase_data
        data = where(data == 0, nan, data) # fake 0 -> nan for the count below

        # manually count # nan cells
        nans = 0
        ys, xs = data.shape
        for y, x in product(range(ys), range(xs)):
            if isnan(data[y, x]):
                nans += 1
        del data

        num_cells = float(ys * xs)
        self.assertTrue(nans > 0)
        self.assertTrue(nans <= num_cells)
        self.assertEqual(nans / num_cells, self.ifg.nan_fraction)

    def test_xy_size(self):
        self.assertFalse(self.ifg.ncols is None)
        self.assertFalse(self.ifg.nrows is None)

        # test with tolerance from base 90m cell
        # within 2% of cells over sydney?
        self.assertTrue(self.ifg.y_size > 88.0)
        self.assertTrue(self.ifg.y_size < 92.0, 'Got %s' % self.ifg.y_size)

        syd_width = 76.9 # from nearby pirate coords
        self.assertTrue(self.ifg.x_size > 0.97 * syd_width)  # ~3% tolerance
        self.assertTrue(self.ifg.x_size < 1.03 * syd_width)

    def test_centre_latlong(self):
        lat_exp = self.ifg.y_first + \
                  (int(self.ifg.nrows / 2) * self.ifg.y_step)
        long_exp = self.ifg.x_first + \
                   (int(self.ifg.ncols / 2) * self.ifg.x_step)
        self.assertEqual(lat_exp, self.ifg.lat_centre)
        self.assertEqual(long_exp, self.ifg.long_centre)

    def test_centre_cell(self):
        self.assertEqual(self.ifg.x_centre, 23)
        self.assertEqual(self.ifg.y_centre, 36)

    def test_time_span(self):
        self.assertAlmostEqual(self.ifg.time_span, 0.287474332649)

    def test_wavelength(self):
        self.assertEqual(self.ifg.wavelength, 0.0562356424)


class IfgIOTests(unittest.TestCase):

    def setUp(self):
        self.ifg = Ifg(join(SYD_TEST_TIF, 'geo_070709-070813_unw.tif'))

    def test_open(self):
        self.assertTrue(self.ifg.dataset is None)
        self.assertTrue(self.ifg.is_open is False)
        self.ifg.open(readonly=True)
        self.assertTrue(self.ifg.dataset is not None)
        self.assertTrue(self.ifg.is_open is True)
        self.assertTrue(isinstance(self.ifg.dataset, Dataset))

        # ensure open cannot be called twice
        self.failUnlessRaises(RasterException, self.ifg.open, True)

    def test_open_ifg_from_dataset(self):
        """
        Test showing open() can not be used for Ifg created with
        gdal.Dataset object as Dataset has already been read in
        """
        paths = [self.ifg.data_path]
        mlooked_phase_data = prepifg.prepare_ifgs(paths,
                             crop_opt=prepifg.ALREADY_SAME_SIZE,
                             xlooks=2,
                             ylooks=2,
                             write_to_disc=False)
        mlooked = [Ifg(m[1]) for m in mlooked_phase_data]
        self.assertRaises(RasterException, mlooked[0].open)

    def test_write(self):
        base = TEMPDIR
        src = self.ifg.data_path
        dest = join(base, basename(self.ifg.data_path))

        # shutil.copy needs to copy writeable permission from src
        os.chmod(src, S_IRGRP | S_IWGRP | S_IWOTH | S_IROTH |
                 S_IRUSR | S_IWUSR)
        shutil.copy(src, dest)
        os.chmod(src, S_IRGRP | S_IROTH | S_IRUSR)  # revert

        i = Ifg(dest)
        i.open()
        i.phase_data[0, 1:] = nan
        i.write_modified_phase()
        del i

        # reopen to ensure data/nans can be read back out
        i = Ifg(dest)
        i.open(readonly=True)
        assert_array_equal(True, isnan(i.phase_data[0, 1:]))
        i.close()
        os.remove(dest)

    def test_readonly_permission_failure(self):
        # ensure failure if opening R/O permission file as writeable/GA_Update
        self.assertRaises(IOError, self.ifg.open, False)

    def test_write_fails_on_readonly(self):
        # check readonly status is same before
        # and after open() for readonly file
        self.assertTrue(self.ifg.is_read_only)
        self.ifg.open(readonly=True)
        self.assertTrue(self.ifg.is_read_only)
        self.assertRaises(IOError, self.ifg.write_modified_phase)

    def test_phase_band_unopened_ifg(self):
        try:
            _ = self.ifg.phase_band
            self.fail("Should not be able to access band without open dataset")
        except RasterException:
            pass

    def test_nan_fraction_unopened(self):
        try:
            # NB: self.assertRaises doesn't work here (as it is a property?)
            _ = self.ifg.nan_fraction
            self.fail("Shouldn't be able to "
                      "call nan_fraction() with unopened Ifg")
        except RasterException:
            pass

    def test_phase_data_properties(self):
        # Use raw GDAL to isolate raster reading from Ifg functionality
        ds = Open(self.ifg.data_path)
        data = ds.GetRasterBand(1).ReadAsArray()
        del ds

        self.ifg.open()

        # test full array and row by row access
        assert_array_equal(data, self.ifg.phase_data)
        for y, row in enumerate(self.ifg.phase_rows):
            assert_array_equal(data[y], row)

        # test the data is cached if changed
        crd = (5,4)
        orig = self.ifg.phase_data[crd]
        self.ifg.phase_data[crd] *= 2
        nv = self.ifg.phase_data[crd] # pull new value out again
        self.assertEqual(nv, 2 * orig)


# FIXME:
# class IncidenceFileTests(unittest.TestCase):
#     'Unit tests to verify operations on GeoTIFF format Incidence rasters'
#
#     def setUp(self):
#         raise NotImplementedError
#         self.inc = Incidence(join(INCID_TEST_DIR, '128x2.tif'))
#         self.inc.open()
#
#
#     def test_incidence_data(self):
#         # check incidences rises while traversing the scene
#         data = self.inc.incidence_data
#         diff = data.ptp()
#         self.assertTrue(diff > 0.5, "Got ptp() diff of %s" % diff)
#
#         # ascending pass, values should increase from W->E across scene
#         for i in range(2):
#             d = data[i]
#             self.assertFalse((d == 0).any()) # ensure no NODATA
#             self.assertFalse((isnan(d)).any())
#
#             diff = array([d[i+1] - d[i] for i in range(len(d)-1)])
#             res = abs(diff[diff < 0])
# TODO: check if this is normal
#             self.assertTrue((res < 1e-4).all())
#
#
#     def test_azimuth_data(self):
#         # ensure azimuth is fairly constant
#
#         az = self.inc.azimuth_data
#         self.assertFalse((az == 0).all())
#         az = az[az != 0] # filter NODATA cells
#
#         # azimuth should be relatively constant
#         ptp = az.ptp()
#         self.assertTrue(ptp < 0.1, msg="min -> max diff is %s" % ptp)


class DEMTests(unittest.TestCase):
    'Unit tests to verify operations on GeoTIFF format DEMs'

    def setUp(self):
        self.ras = DEM(SYD_TEST_DEM_TIF)

    def test_create_raster(self):
        # validate header path
        self.assertTrue(os.path.exists(self.ras.data_path))

    def test_headers_as_attr(self):
        self.ras.open()
        attrs = ['ncols', 'nrows', 'x_first', 'x_step', 'y_first', 'y_step' ]

        # TODO: are 'projection' and 'datum' attrs needed?
        for a in attrs:
            self.assertTrue(getattr(self.ras, a) is not None)

    def test_is_dem(self):
        self.ras = DEM(join(SYD_TEST_TIF, 'geo_060619-061002_unw.tif'))
        self.assertFalse(hasattr(self.ras, 'datum'))

    def test_open(self):
        self.assertTrue(self.ras.dataset is None)
        self.ras.open()
        self.assertTrue(self.ras.dataset is not None)
        self.assertTrue(isinstance(self.ras.dataset, Dataset))

        # ensure open cannot be called twice
        self.failUnlessRaises(RasterException, self.ras.open)

    def test_band(self):
        # test accessing bands with open and unopened datasets
        try:
            _ = self.ras.height_band
            self.fail("Should not be able to access band without open dataset")
        except RasterException:
            pass

        self.ras.open()
        data = self.ras.height_band.ReadAsArray()
        self.assertEqual(data.shape, (72, 47))


class WriteUnwTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tif_dir = tempfile.mkdtemp()
        cls.test_conf = common.SYDNEY_TEST_CONF

        # change the required params
        cls.params = cf.get_config_params(cls.test_conf)
        cls.params[cf.OBS_DIR] = common.SYD_TEST_GAMMA
        cls.params[cf.PROCESSOR] = 1  # gamma
        file_list = list(cf.parse_namelist(os.path.join(common.SYD_TEST_GAMMA,
                                                        'ifms_17')))
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
        cls.base_unw_paths = cf.original_ifg_paths(
            cls.params[cf.IFG_FILE_LIST])
        cls.base_unw_paths.append(common.SYD_TEST_DEM_GAMMA)

        xlks, ylks, crop = cf.transform_params(cls.params)

        # dest_paths are tifs that have been geotif converted and multilooked
        run_prepifg.gamma_prepifg(cls.base_unw_paths, cls.params)
        cls.base_unw_paths.pop()  # removed dem as we don't want it in ifgs

        dest_paths = cf.get_dest_paths(
            cls.base_unw_paths, crop, cls.params, xlks)
        cls.ifgs = common.sydney_data_setup(datafiles=dest_paths)

    @classmethod
    def tearDownClass(cls):
        for i in cls.ifgs:
            i.close()
        shutil.rmtree(cls.tif_dir)

    def test_unw_contains_same_data_as_numpy_array(self):
        from datetime import time
        temp_unw = tempfile.mktemp(suffix='.unw')
        temp_tif = tempfile.mktemp(suffix='.tif')

        # setup some header files for use in write_geotif
        dem_header_file = common.SYD_TEST_DEM_HDR_GAMMA
        dem_header = gamma.parse_dem_header(dem_header_file)

        header = gamma.parse_epoch_header(
            os.path.join(common.SYD_TEST_GAMMA, '20060828_slc.par'))
        header.update(dem_header)

        # insert some dummy data so we are the dem in write_geotiff is not
        # not activated and ifg write_geotiff operation works
        header[ifc.PYRATE_TIME_SPAN] = 0
        header[ifc.SLAVE_DATE] = 0
        header[ifc.DATA_UNITS] = 'degrees'
        header[ifc.DATA_TYPE] = ifc.ORIG
        header[ifc.SLAVE_TIME] = time(10)

        # now create aritrary data
        data = np.random.rand(dem_header[ifc.PYRATE_NROWS],
                              dem_header[ifc.PYRATE_NCOLS])

        # convert numpy array to .unw
        shared.write_unw_from_data_or_geotiff(geotif_or_data=data,
                                              dest_unw=temp_unw,
                                              ifg_proc=1)
        # convert the .unw to geotif
        shared.write_geotiff(header=header, data_path=temp_unw,
                             dest=temp_tif, nodata=np.nan)

        # now compare geotiff with original numpy array
        ds = gdal.Open(temp_tif, gdal.GA_ReadOnly)
        data_lv_theta = ds.ReadAsArray()
        ds = None
        np.testing.assert_array_almost_equal(data, data_lv_theta)
        os.remove(temp_tif)
        os.remove(temp_unw)

    def test_equality_of_unw_with_geotiff(self):
        geotiffs = [os.path.join(
            self.params[cf.OUT_DIR], os.path.basename(b).split('.')[0] + '_' 
            + os.path.basename(b).split('.')[1] + '.tif')
            for b in self.base_unw_paths]

        # create .unws from geotiffs and make sure they read the same
        dest_unws = []
        for g in geotiffs:
            dest_unw = os.path.join(self.params[cf.OUT_DIR],
                         os.path.splitext(g)[0] + '.unw')
            shared.write_unw_from_data_or_geotiff(
                geotif_or_data=g, dest_unw= dest_unw, ifg_proc=1)
            dest_unws.append(dest_unw)

        new_geotiffs = [run_prepifg.gamma_multiprocessing(b, self.params)
                        for b in dest_unws]

        for g, u in zip(geotiffs, new_geotiffs):
            g_ds = gdal.Open(g)
            u_gs = gdal.Open(u)
            np.testing.assert_array_almost_equal(u_gs.ReadAsArray(),
                                                 g_ds.ReadAsArray())
            u_gs = None
            g_ds = None

    def test_unws_created_are_same_as_original(self):
        geotiffs = [os.path.join(
            self.params[cf.OUT_DIR], os.path.basename(b).split('.')[0] + '_' 
            + os.path.basename(b).split('.')[1] + '.tif')
            for b in self.base_unw_paths]

        new_base_unw_paths = []
        # create .unws from geotiffs and make sure they read the same
        for g in geotiffs:
            dest_unw = os.path.join(self.params[cf.OUT_DIR],
                                    os.path.splitext(g)[0] + '.unw')
            shared.write_unw_from_data_or_geotiff(
                geotif_or_data=g, dest_unw=dest_unw, ifg_proc=1)

            # unws created
            assert os.path.exists(dest_unw)
            new_base_unw_paths.append(dest_unw)

        # make sure we can recovert the unws to gettiffs
        new_geotiffs = [run_prepifg.gamma_multiprocessing(b, self.params)
                        for b in new_base_unw_paths]

        # assert data equal
        for g, u in zip(geotiffs, new_geotiffs):
            g_ds = gdal.Open(g)
            u_gs = gdal.Open(u)
            np.testing.assert_array_almost_equal(u_gs.ReadAsArray(),
                                                 g_ds.ReadAsArray())
            u_gs = None
            g_ds = None

    def test_roipac_raises(self):
        geotiffs = [os.path.join(
            self.params[cf.OUT_DIR], os.path.basename(b).split('.')[0] + '_' 
            + os.path.basename(b).split('.')[1] + '.tif')
            for b in self.base_unw_paths]

        for g in geotiffs[:1]:
            dest_unw = os.path.join(self.params[cf.OUT_DIR],
                                    os.path.splitext(g)[0] + '.unw')
            with self.assertRaises(NotImplementedError):
                shared.write_unw_from_data_or_geotiff(
                    geotif_or_data=g, dest_unw=dest_unw, ifg_proc=0)


class GeodesyTests(unittest.TestCase):

    def test_utm_zone(self):
        # test some different zones (collected manually)
        for lon in [174.0, 176.5, 179.999, 180.0]:
            self.assertEqual(60, utm_zone(lon))

        for lon in [144.0, 144.1, 146.3456, 149.9999]:
            self.assertEqual(55, utm_zone(lon))

        for lon in [-180.0, -179.275, -176.925]:
            self.assertEqual(1, utm_zone(lon))

        for lon in [-72.0, -66.1]:
            self.assertEqual(19, utm_zone(lon))

        for lon in [0.0, 0.275, 3.925, 5.999]:
            self.assertEqual(31, utm_zone(lon))


    def test_cell_size_polar_region(self):
        # Can't have polar area zones: see http://www.dmap.co.uk/utmworld.htm
        for lat in [-80.1, -85.0, -90.0, 84.1, 85.0, 89.9999, 90.0]:
            self.assertRaises(ValueError, cell_size, lat, 0, 0.1, 0.1)


    def test_cell_size_calc(self):
        # test conversion of X|Y_STEP to X|Y_SIZE
        x_deg = 0.000833333
        y_deg = -x_deg

        approx = 90.0 # x_deg is approx 90m
        exp_low = approx - (.15 * approx) # assumed tolerance
        exp_high = approx + (.15 * approx)

        latlons = [(10.0, 15.0), (-10.0, 15.0), (10.0, -15.0), (-10.0, -15.0),
            (178.0, 33.0), (-178.0, 33.0), (178.0, -33.0), (-178.0, -33.0) ]

        for lon, lat in latlons:
            xs, ys = cell_size(lat, lon, x_deg, y_deg)
            for s in (xs, ys):
                self.assertTrue(s > 0, msg="size=%s" % s)
                self.assertTrue(s > exp_low, msg="size=%s" % s)
                self.assertTrue(s < exp_high, msg="size=%s" % s)


if __name__ == "__main__":
    unittest.main()
