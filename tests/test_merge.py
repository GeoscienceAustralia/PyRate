# coding: utf-8
#   This Python module is part of the PyRate software package.
#
#   Copyright 2020 Geoscience Australia
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
from pyrate.merge import _merge_stack, _merge_linrate, _merge_timeseries
from pyrate.core.shared import DEM
from pyrate.core.ifgconstants import LOS_PROJECTION_OPTION
from pyrate.configuration import Configuration, write_config_file
from tests.common import manipulate_test_conf, PY37GDAL304, MEXICO_CROPA_CONF, assert_same_files_produced


@pytest.fixture(params=list(LOS_PROJECTION_OPTION.keys()))
def los_projection(request):
    return request.param


@pytest.fixture(scope='module')
def create_pre_merge_output():
    tdir = Path(tempfile.mkdtemp())
    params = manipulate_test_conf(MEXICO_CROPA_CONF, tdir)
    output_conf_file = tdir.joinpath('conf.cfg')
    output_conf = tdir.joinpath(output_conf_file)
    write_config_file(params=params, output_conf_file=output_conf)
    check_call(f"mpirun -n 2 pyrate conv2tif -f {output_conf}", shell=True, env=os.environ)
    check_call(f"mpirun -n 2 pyrate prepifg -f {output_conf}", shell=True, env=os.environ)
    check_call(f"mpirun -n 2 pyrate correct -f {output_conf}", shell=True, env=os.environ)
    check_call(f"mpirun -n 2 pyrate timeseries -f {output_conf}", shell=True, env=os.environ)
    check_call(f"mpirun -n 2 pyrate stack -f {output_conf}", shell=True, env=os.environ)

    params = Configuration(output_conf).__dict__
    return params


def test_los_conversion_comparison(create_pre_merge_output):
    """
    compare outputs in each of the los projection types
    compares sine and cosine components are larger than LOS component
    assert metadata equal except
    """
    params = create_pre_merge_output
    all_dirs = {}
    for k in LOS_PROJECTION_OPTION.keys():
        params[C.LOS_PROJECTION] = k
        _merge_stack(params)
        _merge_linrate(params)
        _merge_timeseries(params, 'tscuml')
        _merge_timeseries(params, 'tsincr')
        k_dir = Path(params[C.OUT_DIR]).joinpath(ifc.LOS_PROJECTION_OPTION[k])
        k_dir.mkdir(exist_ok=True)

        for out_type in los_projection_out_types:
            print(out_type)
            print(len(list(Path(params[C.OUT_DIR]).glob(out_type + '*.tif'))))
            for tif in Path(params[C.OUT_DIR]).glob(out_type + '*.tif'):
                print(tif)
                shutil.move(tif, k_dir)
        all_dirs[k] = k_dir

    los_proj_dir = all_dirs[ifc.LINE_OF_SIGHT]
    pseudo_ver = all_dirs[ifc.PSEUDO_VERTICAL]
    pseudo_hor = all_dirs[ifc.PSEUDO_HORIZONTAL]
    assert len(list(los_proj_dir.glob('*.tif'))) == 26  # 12 tsincr, 12 tscuml + 1 stack rate + 1 linear rate
    for tif in los_proj_dir.glob('*.tif'):
        ds = DEM(tif)
        ds_ver = DEM(pseudo_ver.joinpath(tif.name))
        ds_hor = DEM(pseudo_hor.joinpath(tif.name))
        ds.open()
        ds_ver.open()
        ds_hor.open()
        non_nans_indices = ~np.isnan(ds.data)
        # assert division by sine and cosine always yields larger components in vertical and horizontal directions
        assert np.all(np.abs(ds.data[non_nans_indices]) <= np.abs(ds_ver.data[non_nans_indices]))
        assert np.all(np.abs(ds.data[non_nans_indices]) <= np.abs(ds_hor.data[non_nans_indices]))
        ds_md = ds.dataset.GetMetadata()
        assert ds_md.pop(C.LOS_PROJECTION.upper()) == ifc.LOS_PROJECTION_OPTION[ifc.LINE_OF_SIGHT]
        ds_ver_md = ds_ver.dataset.GetMetadata()
        assert ds_ver_md.pop(C.LOS_PROJECTION.upper()) == ifc.LOS_PROJECTION_OPTION[ifc.PSEUDO_VERTICAL]
        assert ds_md == ds_ver_md
        ds_hor_md = ds_hor.dataset.GetMetadata()
        assert ds_hor_md.pop(C.LOS_PROJECTION.upper()) == ifc.LOS_PROJECTION_OPTION[ifc.PSEUDO_HORIZONTAL]
        assert ds_md == ds_hor_md


@pytest.mark.mpi
@pytest.mark.slow
@pytest.mark.skipif(not PY37GDAL304, reason="Only run in one CI env")
def test_file_creation(create_pre_merge_output, los_projection):
    params = create_pre_merge_output
    params[C.LOS_PROJECTION] = los_projection
    try:
        _merge_stack(params)
        _merge_linrate(params)
        _merge_timeseries(params, 'tscuml')
        _merge_timeseries(params, 'tsincr')
    except RuntimeError:
        return
    # check if color map is created
    for ot in ['stack_rate', 'stack_error', 'linear_rate', 'linear_error', 'linear_rsquared']:
        create_png_and_kml_from_tif(params[C.OUT_DIR], output_type=ot)
        output_color_map_path = os.path.join(params[C.OUT_DIR], f"colourmap_{ot}.txt")
        assert Path(output_color_map_path).exists(), "Output color map file not found at: " + output_color_map_path

    # check if merged files are created
    for _type, ot in itertools.product(['stack_rate', 'stack_error', 'linear_rate',
                                        'linear_error', 'linear_rsquared'], ['.tif', '.png', '.kml']):
        output_image_path = os.path.join(params[C.OUT_DIR], _type + ot)
        print(f"checking {output_image_path}")
        assert Path(output_image_path).exists(), f"Output {ot} file not found at {output_image_path}"

    # check los_projection metadata
    for out_type in los_projection_out_types:
        for tif in Path(params[C.OUT_DIR]).glob(out_type + '*.tif'):
            __check_md(los_projection, tif.as_posix())
    shutil.rmtree(params[C.OUT_DIR], ignore_errors=True)


def __check_md(los_projection, output_image_path):
    ifg = DEM(output_image_path)
    ifg.open()
    assert ifg.dataset.GetMetadataItem(C.LOS_PROJECTION.upper()) == LOS_PROJECTION_OPTION[los_projection]
