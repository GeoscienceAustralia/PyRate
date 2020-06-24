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
from pyrate.merge import _create_png_from_tif
from pyrate.core import config as cf
from pyrate.merge import _merge_stack
from pyrate.configuration import Configuration
from tests.common import manipulate_test_conf


@pytest.fixture
def create_stack_output(tempdir, gamma_conf):
    tdir = Path(tempdir())
    params = manipulate_test_conf(gamma_conf, tdir)
    output_conf_file = tdir.joinpath('conf.cfg')
    output_conf = tdir.joinpath(output_conf_file)
    cf.write_config_file(params=params, output_conf_file=output_conf)
    check_call(f"pyrate conv2tif -f {output_conf}", shell=True)
    check_call(f"pyrate prepifg -f {output_conf}", shell=True)
    check_call(f"pyrate process -f {output_conf}", shell=True)

    params = Configuration(output_conf).__dict__
    return params, tdir


@pytest.mark.slow
def test_png_creation(create_stack_output):
    params, tdir = create_stack_output

    output_folder_path = params[cf.OUT_DIR]

    rows, cols = params["rows"], params["cols"]
    _merge_stack(rows, cols, params)
    _create_png_from_tif(output_folder_path)

    # check if color map is created
    for out_type in ['rate', 'error']:
        output_color_map_path = os.path.join(output_folder_path, f"colourmap_{out_type}.txt")
        assert Path(output_color_map_path).exists(), "Output color map file not found at: " + output_color_map_path

    # check if merged files are created
    for _type, output_type in itertools.product(["stack_rate", "stack_error"], ['.tif', '.png', '.kml']):
        output_image_path = os.path.join(output_folder_path, _type + output_type)
        print(f"checking {output_image_path}")
        assert Path(output_image_path).exists(), f"Output {output_type} file not found at {output_image_path}"
