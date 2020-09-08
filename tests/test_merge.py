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
import pytest
from pathlib import Path
from pyrate.merge import create_png_and_kml_from_tif
from pyrate.core import config as cf
from pyrate.merge import _merge_stack, _merge_linrate
from pyrate.configuration import Configuration, write_config_file
from tests.common import manipulate_test_conf


@pytest.fixture
def create_merge_output(tempdir, gamma_conf):
    tdir = Path(tempdir())
    params = manipulate_test_conf(gamma_conf, tdir)
    output_conf_file = tdir.joinpath('conf.cfg')
    output_conf = tdir.joinpath(output_conf_file)
    write_config_file(params=params, output_conf_file=output_conf)
    check_call(f"pyrate conv2tif -f {output_conf}", shell=True)
    check_call(f"pyrate prepifg -f {output_conf}", shell=True)
    check_call(f"pyrate correct -f {output_conf}", shell=True)
    check_call(f"pyrate timeseries -f {output_conf}", shell=True)
    check_call(f"pyrate stack -f {output_conf}", shell=True)

    params = Configuration(output_conf).__dict__
    _merge_stack(params)
    _merge_linrate(params)
    return params


@pytest.mark.slow
def test_file_creation(create_merge_output):
    params = create_merge_output

    # check if color map is created
    for ot in ['stack_rate', 'stack_error', 'linear_rate', 'linear_error', 'linear_rsquared']:
        create_png_and_kml_from_tif(params[cf.OUT_DIR], output_type=ot)
        output_color_map_path = os.path.join(params[cf.OUT_DIR], f"colourmap_{ot}.txt")
        assert Path(output_color_map_path).exists(), "Output color map file not found at: " + output_color_map_path

    # check if merged files are created
    for _type, ot in itertools.product(['stack_rate', 'stack_error', 'linear_rate',
                                        'linear_error', 'linear_rsquared'], ['.tif', '.png', '.kml']):
        output_image_path = os.path.join(params[cf.OUT_DIR], _type + ot)
        print(f"checking {output_image_path}")
        assert Path(output_image_path).exists(), f"Output {ot} file not found at {output_image_path}"
