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
import os
from subprocess import check_call
import itertools
import tempfile
import pickle
import pytest
from pathlib import Path

import pyrate.constants as C
from pyrate.merge import create_png_and_kml_from_tif, los_projection_out_types
from pyrate.merge import _merge_stack, _merge_linrate, _merge_timeseries
from pyrate.core.shared import DEM
from pyrate.core.ifgconstants import LOS_PROJECTION_OPTION
from pyrate.configuration import Configuration, write_config_file
from tests.common import manipulate_test_conf, PY37GDAL304, MEXICO_CROPA_CONF


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
    check_call(f"mpirun -n 3 pyrate conv2tif -f {output_conf}", shell=True)
    check_call(f"mpirun -n 3 pyrate prepifg -f {output_conf}", shell=True)
    check_call(f"mpirun -n 3 pyrate correct -f {output_conf}", shell=True)
    check_call(f"mpirun -n 3 pyrate timeseries -f {output_conf}", shell=True)
    check_call(f"mpirun -n 3 pyrate stack -f {output_conf}", shell=True)

    params = Configuration(output_conf).__dict__
    return params


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


def __check_md(los_projection, output_image_path):
    ifg = DEM(output_image_path)
    ifg.open()
    assert ifg.dataset.GetMetadataItem(C.LOS_PROJECTION) == LOS_PROJECTION_OPTION[los_projection]
