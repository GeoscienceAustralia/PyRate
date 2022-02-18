#   This Python module is part of the PyRate software package.
#
#   Copyright 2022 Geoscience Australia
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
This Python module contains tests for the Merge step of PyRate.
"""
import shutil
import os
from subprocess import check_call
import itertools
import tempfile
import pytest
from pathlib import Path

import numpy as np

import pyrate.constants as C
import pyrate.core.ifgconstants as ifc
from pyrate.merge import create_png_and_kml_from_tif, los_projection_out_types
from pyrate.merge import _merge_stack, _merge_linrate, _merge_timeseries, los_projection_divisors
from pyrate.core.shared import DEM
from pyrate.core.ifgconstants import LOS_PROJECTION_OPTION
from pyrate.configuration import Configuration, write_config_file
from tests.common import manipulate_test_conf, PY37GDAL302, MEXICO_CROPA_CONF, assert_same_files_produced


@pytest.fixture(params=list(LOS_PROJECTION_OPTION.keys()))
def los_projection(request):
    return request.param


@pytest.fixture(params=[-1, 1])
def signal_polarity(request):
    return request.param


def test_los_conversion_divisors():
    """
    Unit test to check the LOS conversions for specific incidence angles
    """
    inc = [0, 30, 45, 90]  # incidence angles in degrees

    # Test pseudo-vertical
    res = [los_projection_divisors[ifc.PSEUDO_VERTICAL](np.radians(x)) for x in inc]
    exp = [1.0, 0.8660254037844387, 0.7071067811865476, 6.123233995736766e-17]
    np.testing.assert_almost_equal(res, exp, decimal=6)

    # Test pseudo-horizontal
    res = [los_projection_divisors[ifc.PSEUDO_HORIZONTAL](np.radians(x)) for x in inc]
    exp = [0.0, 0.49999999999999994, 0.7071067811865475, 1.0]
    np.testing.assert_almost_equal(res, exp, decimal=6)

    # Test line-of-sight
    res = [los_projection_divisors[ifc.LINE_OF_SIGHT](np.radians(x)) for x in inc]
    exp = [1.0, 1.0, 1.0, 1.0]
    np.testing.assert_almost_equal(res, exp, decimal=1)


@pytest.mark.mpi
@pytest.mark.slow
@pytest.mark.skipif(not PY37GDAL302, reason="Only run in one CI env")
class TestLOSConversion:
    @classmethod
    def setup_class(cls):
        tdir = Path(tempfile.mkdtemp())
        params = manipulate_test_conf(MEXICO_CROPA_CONF, tdir)
        output_conf_file = tdir.joinpath('conf.cfg')
        output_conf = tdir.joinpath(output_conf_file)
        write_config_file(params=params, output_conf_file=output_conf)
        check_call(f"mpirun -n 3 pyrate conv2tif -f {output_conf}", shell=True)
        check_call(f"mpirun -n 3 pyrate prepifg -f {output_conf}", shell=True)
        check_call(f"mpirun -n 3 pyrate correct -f {output_conf}", shell=True)
        check_call(f"mpirun -n 3 pyrate timeseries -f {output_conf}", shell=True)
        check_call(f"mpirun -n 3 pyrate stack -f {output_conf}", shell=True)

        params = Configuration(output_conf).__dict__
        cls.params = params
        cls.tdir = tdir

    @classmethod
    def teardown_class(cls):
        shutil.rmtree(cls.tdir, ignore_errors=True)

    def test_los_conversion_comparison(self):
        """
        compare outputs in each of the los projection types
        compares sine and cosine components are larger than LOS component
        assert metadata equal except
        """
        params = self.params
        all_dirs = {}
        params[C.SIGNAL_POLARITY] = -1
        for k in LOS_PROJECTION_OPTION.keys():
            params[C.LOS_PROJECTION] = k
            k_dir = Path(params[C.OUT_DIR]).joinpath(ifc.LOS_PROJECTION_OPTION[k])
            k_dir.mkdir(exist_ok=True)
            all_dirs[k] = k_dir
            self.run_with_new_params(k_dir, params)

        signal_dir = Path(params[C.OUT_DIR]).joinpath('signal_polarity_dir')
        signal_dir.mkdir(exist_ok=True)
        all_dirs[C.SIGNAL_POLARITY] = signal_dir

        params[C.SIGNAL_POLARITY] = 1
        self.run_with_new_params(signal_dir, params)

        los_proj_dir = all_dirs[ifc.LINE_OF_SIGHT]
        pseudo_ver = all_dirs[ifc.PSEUDO_VERTICAL]
        pseudo_hor = all_dirs[ifc.PSEUDO_HORIZONTAL]
        num_files = 24 if params[C.PHASE_CLOSURE] else 26  # phase closure removes 1 interferogram
        assert len(list(los_proj_dir.glob('*.tif'))) == num_files  # 12 tsincr, 12 tscuml + 1 stack rate + 1 linear rate
        signal_polarity_reversed_pseudo_hor = all_dirs[C.SIGNAL_POLARITY]

        for tif in los_proj_dir.glob('*.tif'):
            ds = DEM(tif)
            ds_ver = DEM(pseudo_ver.joinpath(tif.name))
            ds_hor = DEM(pseudo_hor.joinpath(tif.name))
            ds_hor_sig = DEM(signal_polarity_reversed_pseudo_hor.joinpath(tif.name))
            ds.open()
            ds_ver.open()
            ds_hor.open()
            ds_hor_sig.open()

            non_nans_indices = ~np.isnan(ds.data)
            # assert division by sine and cosine always yields larger components in vertical and horizontal directions
            assert np.all(np.abs(ds.data[non_nans_indices]) <= np.abs(ds_ver.data[non_nans_indices]))
            assert np.all(np.abs(ds.data[non_nans_indices]) <= np.abs(ds_hor.data[non_nans_indices]))
            assert np.all(np.abs(ds.data[non_nans_indices]) <= np.abs(ds_hor_sig.data[non_nans_indices]))
            assert np.all(ds_hor.data[non_nans_indices] == -ds_hor_sig.data[non_nans_indices])
            ds_md = ds.dataset.GetMetadata()
            assert ds_md.pop(C.LOS_PROJECTION.upper()) == ifc.LOS_PROJECTION_OPTION[ifc.LINE_OF_SIGHT]
            ds_ver_md = ds_ver.dataset.GetMetadata()
            assert ds_ver_md.pop(C.LOS_PROJECTION.upper()) == ifc.LOS_PROJECTION_OPTION[ifc.PSEUDO_VERTICAL]
            assert ds_md == ds_ver_md
            assert ds_md.pop(C.SIGNAL_POLARITY.upper()) == '-1'
            ds_hor_md = ds_hor.dataset.GetMetadata()
            ds_hor_sig_md = ds_hor_sig.dataset.GetMetadata()
            assert ds_hor_sig_md.pop(C.SIGNAL_POLARITY.upper()) != ds_hor_md.pop(C.SIGNAL_POLARITY.upper())
            assert ds_hor_sig_md == ds_hor_md
            assert ds_hor_md.pop(C.LOS_PROJECTION.upper()) == ifc.LOS_PROJECTION_OPTION[ifc.PSEUDO_HORIZONTAL]
            assert ds_hor_sig_md.pop(C.LOS_PROJECTION.upper()) == ifc.LOS_PROJECTION_OPTION[ifc.PSEUDO_HORIZONTAL]
            assert ds_md == ds_hor_md
            assert ds_md == ds_hor_sig_md

    def run_with_new_params(self, k_dir, params):
        _merge_stack(params)
        _merge_linrate(params)
        _merge_timeseries(params, 'tscuml')
        _merge_timeseries(params, 'tsincr')

        for out_type in los_projection_out_types:
            for tif in itertools.chain(Path(params[C.VELOCITY_DIR]).glob(out_type + '*.tif'),
                                       Path(params[C.TIMESERIES_DIR]).glob(out_type + '*.tif')):
                shutil.move(tif, k_dir.joinpath(tif.name))

    def test_file_creation(self, los_projection):
        params = self.params
        params[C.LOS_PROJECTION] = los_projection
        _merge_stack(params)
        _merge_linrate(params)
        _merge_timeseries(params, 'tscuml')
        _merge_timeseries(params, 'tsincr')
        # check if color map is created
        for ot in ['stack_rate', 'stack_error', 'linear_rate', 'linear_error', 'linear_rsquared']:
            create_png_and_kml_from_tif(params[C.VELOCITY_DIR], output_type=ot)
            output_color_map_path = os.path.join(params[C.VELOCITY_DIR], f"colourmap_{ot}.txt")
            assert Path(output_color_map_path).exists(), "Output color map file not found at: " + output_color_map_path

        # check if merged files are created
        for _type, ot in itertools.product(['stack_rate', 'stack_error', 'linear_rate',
                                            'linear_error', 'linear_rsquared'], ['.tif', '.png', '.kml']):
            output_image_path = os.path.join(params[C.VELOCITY_DIR], _type + ot)
            print(f"checking {output_image_path}")
            assert Path(output_image_path).exists(), f"Output {ot} file not found at {output_image_path}"

        # check los_projection metadata
        for out_type in los_projection_out_types:
            for tif in Path(params[C.TIMESERIES_DIR]).glob(out_type + '*.tif'):
                self.__check_md(los_projection, tif.as_posix())

    @staticmethod
    def __check_md(los_projection, output_image_path):
        ifg = DEM(output_image_path)
        ifg.open()
        assert ifg.dataset.GetMetadataItem(C.LOS_PROJECTION.upper()) == LOS_PROJECTION_OPTION[los_projection]
