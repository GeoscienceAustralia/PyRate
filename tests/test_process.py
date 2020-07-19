import pytest
from pyrate import process
from pyrate.core.config import ConfigException


def test_unsupported_process_steps_raises(gamma_params):
    gamma_params['process'] = ['orbfit2', 'something_other_step']
    with pytest.raises(ConfigException):
        process.process_ifgs(gamma_params)


def test_supported_process_steps_dont_raise(gamma_params):
    supported_stpes = ['orbfit', 'refphase', 'mst', 'apscorrect', 'maxvar', 'timeseries', 'stack']
    assert all([s in gamma_params['process'] for s in supported_stpes])
    process.__validate_process_steps(params=gamma_params)


def test_process_treats_mlooked_files_as_readonly():
    return
