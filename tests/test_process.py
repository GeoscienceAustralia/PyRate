from pathlib import Path
import pytest
from pyrate import process, prepifg, conv2tif
import pyrate.core.config as cf
from pyrate.core.config import ConfigException
from tests import common


def test_unsupported_process_steps_raises(gamma_params):
    gamma_params['process'] = ['orbfit2', 'something_other_step']
    with pytest.raises(ConfigException):
        process.process_ifgs(gamma_params)


def test_supported_process_steps_dont_raise(gamma_params):
    supported_stpes = ['orbfit', 'refphase', 'mst', 'apscorrect', 'maxvar', 'timeseries', 'stack']
    assert all([s in gamma_params['process'] for s in supported_stpes])
    process.__validate_process_steps(params=gamma_params)


def test_process_treats_prepif_outputs_readonly(gamma_conf, tempdir, coh_mask):
    from pyrate.configuration import Configuration
    tdir = Path(tempdir())
    params = common.manipulate_test_conf(gamma_conf, tdir)
    params[cf.COH_MASK] = coh_mask
    params[cf.PARALLEL] = 0
    output_conf = tdir.joinpath('conf.cfg')
    cf.write_config_file(params=params, output_conf_file=output_conf)
    params = Configuration(output_conf).__dict__
    conv2tif.main(params)
    tifs = list(Path(params[cf.OUT_DIR]).glob('*_unw_ifg.tif'))
    assert len(tifs) == 17

    params = Configuration(output_conf).__dict__
    prepifg.main(params)
    cropped = list(Path(params[cf.OUT_DIR]).glob('*cr.tif'))

    if coh_mask:  # 17 + 1 dem + 17 coh files
        assert len(cropped) == 35
    else:  # 17 + 1 dem
        assert len(cropped) == 18
    # check all tifs from conv2tif are still readonly
    for t in tifs:
        assert t.stat().st_mode == 33060

    # check all prepifg outputs are readonly
    for c in cropped:
        assert c.stat().st_mode == 33060

    params = Configuration(output_conf).__dict__
    process.main(params)

    # check all after process steps multilooked files are still readonly
    for c in cropped:
        assert c.stat().st_mode == 33060
