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

import pytest
from pathlib import Path
import numpy as np

import pyrate.constants as C
from pyrate.configuration import Configuration, write_config_file
from tests.common import MEXICO_CROPA_CONF, manipulate_test_conf, PY37GDAL302, sub_process_run


@pytest.fixture
def modified_config(tempdir, get_lks=1, get_crop=1, orbfit_lks=2, orbfit_method=1, orbfit_degrees=1, ref_est_method=1):
    def modify_params(conf_file, parallel_vs_serial, output_conf_file):
        tdir = Path(tempdir())
        params = manipulate_test_conf(conf_file, tdir)

        if params[C.PROCESSOR] == 1:  # turn on coherence for gamma
            params[C.COH_MASK] = 1

        params[C.PARALLEL] = parallel_vs_serial
        params[C.PROCESSES] = 4
        params[C.APSEST] = 1
        params[C.IFG_LKSX], params[C.IFG_LKSY] = get_lks, get_lks
        params[C.REFNX], params[C.REFNY] = 2, 2

        params[C.IFG_CROP_OPT] = get_crop
        params[C.ORBITAL_FIT_LOOKS_X], params[
            C.ORBITAL_FIT_LOOKS_Y] = orbfit_lks, orbfit_lks
        params[C.ORBITAL_FIT] = 1
        params[C.ORBITAL_FIT_METHOD] = orbfit_method
        params[C.ORBITAL_FIT_DEGREE] = orbfit_degrees
        params[C.REF_EST_METHOD] = ref_est_method
        params[C.MAX_LOOP_LENGTH] = 3
        params["rows"], params["cols"] = 3, 2
        params["savenpy"] = 1
        params["notiles"] = params["rows"] * params["cols"]  # number of tiles

        # write new temp config
        output_conf = tdir.joinpath(output_conf_file)
        write_config_file(params=params, output_conf_file=output_conf)

        return output_conf, params

    return modify_params


@pytest.mark.mpi
@pytest.mark.slow
@pytest.mark.skipif((not PY37GDAL302), reason="Only run in one CI env")
def test_mpi_vs_single_process(modified_config):
    mpi_conf, m_params = modified_config(MEXICO_CROPA_CONF, 0, 'mpi_conf.conf')
    sub_process_run(f"mpirun -n 3 pyrate conv2tif -f {mpi_conf}")
    sub_process_run(f"mpirun -n 3 pyrate prepifg -f {mpi_conf}")
    sub_process_run(f"mpirun -n 3 pyrate correct -f {mpi_conf}")

    serial_conf, s_params = modified_config(MEXICO_CROPA_CONF, 0, 'single_conf.conf')
    sub_process_run(f"pyrate conv2tif -f {serial_conf}")
    sub_process_run(f"pyrate prepifg -f {serial_conf}")
    sub_process_run(f"pyrate correct -f {serial_conf}")

    parallel_conf, p_params = modified_config(MEXICO_CROPA_CONF, 1, 'parallel_conf.conf')
    sub_process_run(f"pyrate conv2tif -f {parallel_conf}")
    sub_process_run(f"pyrate prepifg -f {parallel_conf}")
    sub_process_run(f"pyrate correct -f {parallel_conf}")

    m_config = Configuration(mpi_conf)
    s_config = Configuration(serial_conf)
    p_config = Configuration(parallel_conf)
    m_closure = np.load(m_config.closure().closure)
    s_closure = np.load(s_config.closure().closure)
    p_closure = np.load(p_config.closure().closure)

    # loops
    m_loops = np.load(m_config.closure().loops, allow_pickle=True)
    s_loops = np.load(s_config.closure().loops, allow_pickle=True)
    p_loops = np.load(p_config.closure().loops, allow_pickle=True)
    m_weights = [m.weight for m in m_loops]
    s_weights = [m.weight for m in s_loops]
    p_weights = [m.weight for m in p_loops]
    np.testing.assert_array_equal(m_weights, s_weights)
    np.testing.assert_array_equal(m_weights, p_weights)

    for i, (m, s) in enumerate(zip(m_loops, s_loops)):
        assert all(m_e == s_e for m_e, s_e in zip(m.edges, s.edges))

    # closure
    np.testing.assert_array_almost_equal(np.abs(m_closure), np.abs(s_closure), decimal=4)
    np.testing.assert_array_almost_equal(np.abs(m_closure), np.abs(p_closure), decimal=4)

    # num_occurrences_each_ifg
    m_num_occurrences_each_ifg = np.load(m_config.closure().num_occurences_each_ifg, allow_pickle=True)
    s_num_occurrences_each_ifg = np.load(s_config.closure().num_occurences_each_ifg, allow_pickle=True)
    p_num_occurrences_each_ifg = np.load(p_config.closure().num_occurences_each_ifg, allow_pickle=True)
    np.testing.assert_array_equal(m_num_occurrences_each_ifg, s_num_occurrences_each_ifg)
    np.testing.assert_array_equal(m_num_occurrences_each_ifg, p_num_occurrences_each_ifg)

    # check ps
    m_ifgs_breach_count = np.load(m_config.closure().ifgs_breach_count)
    s_ifgs_breach_count = np.load(s_config.closure().ifgs_breach_count)
    p_ifgs_breach_count = np.load(p_config.closure().ifgs_breach_count)
    np.testing.assert_array_equal(m_ifgs_breach_count, s_ifgs_breach_count)
    np.testing.assert_array_equal(m_ifgs_breach_count, p_ifgs_breach_count)
