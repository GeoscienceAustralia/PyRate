import shutil
import pytest
from pathlib import Path
from subprocess import check_call, CalledProcessError, run
import numpy as np

from pyrate.configuration import Configuration, write_config_file
from pyrate.core import config as cf
from tests.common import MEXICO_CROPA_CONF, manipulate_test_conf, PYTHON3P8, sub_process_run


@pytest.fixture()
def modified_config(tempdir, get_lks=1, get_crop=1, orbfit_lks=2, orbfit_method=1, orbfit_degrees=1, ref_est_method=1):
    def modify_params(conf_file, parallel_vs_serial, output_conf_file):
        tdir = Path(tempdir())
        params = manipulate_test_conf(conf_file, tdir)

        if params[cf.PROCESSOR] == 1:  # turn on coherence for gamma
            params[cf.COH_MASK] = 1

        params[cf.PARALLEL] = parallel_vs_serial
        params[cf.PROCESSES] = 4
        params[cf.APSEST] = 1
        params[cf.IFG_LKSX], params[cf.IFG_LKSY] = get_lks, get_lks
        params[cf.REFNX], params[cf.REFNY] = 2, 2

        params[cf.IFG_CROP_OPT] = get_crop
        params[cf.ORBITAL_FIT_LOOKS_X], params[cf.ORBITAL_FIT_LOOKS_Y] = orbfit_lks, orbfit_lks
        params[cf.ORBITAL_FIT] = 1
        params[cf.ORBITAL_FIT_METHOD] = orbfit_method
        params[cf.ORBITAL_FIT_DEGREE] = orbfit_degrees
        params[cf.REF_EST_METHOD] = ref_est_method
        params["rows"], params["cols"] = 3, 2
        params["savenpy"] = 1
        params["notiles"] = params["rows"] * params["cols"]  # number of tiles

        print(params)
        # write new temp config
        output_conf = tdir.joinpath(output_conf_file)
        write_config_file(params=params, output_conf_file=output_conf)

        return output_conf, params
    return modify_params


@pytest.mark.slow
@pytest.mark.skipif(not PYTHON3P8, reason="Only run Python3.8 env")
def test_mpi_vs_single_process(modified_config):
    mpi_conf, m_params = modified_config(MEXICO_CROPA_CONF, 0, 'mpi_conf.conf')
    sub_process_run(f"mpirun -n 3 pyrate conv2tif -f {mpi_conf}")
    sub_process_run(f"mpirun -n 3 pyrate prepifg -f {mpi_conf}")
    sub_process_run(f"mpirun -n 3 pyrate correct -f {mpi_conf}")

    serial_conf, s_params = modified_config(MEXICO_CROPA_CONF, 0, 'single_conf.conf')
    sub_process_run(f"pyrate conv2tif -f {serial_conf}")
    sub_process_run(f"pyrate prepifg -f {serial_conf}")
    sub_process_run(f"pyrate correct -f {serial_conf}")

    m_config = Configuration(mpi_conf)
    s_config = Configuration(serial_conf)
    m_closure = np.load(m_config.closure().closure)
    s_closure = np.load(s_config.closure().closure)
    m_loops = np.load(m_config.closure().loops, allow_pickle=True)
    s_loops = np.load(s_config.closure().loops, allow_pickle=True)
    m_weights = [m.weight for m in m_loops]
    s_weights = [m.weight for m in s_loops]
    print(m_weights)
    print(s_weights)

    for i, (m, s) in enumerate(zip(m_loops, s_loops)):
        print(f'========================={i}=================')
        assert m.weight == s.weight
        print(m.weight, m.edges)
        print(s.weight, s.edges)
        for m_e, s_e in zip(m.edges, s.edges):
            print('-'*10)
            print(m_e)
            print(s_e)
            assert m_e == s_e

    import IPython; IPython.embed(); import sys; sys.exit()
    np.testing.assert_array_almost_equal(m_closure, s_closure)
