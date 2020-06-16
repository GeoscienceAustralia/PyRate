import os
import shutil
import pytest
import glob
import copy
import itertools
from pathlib import Path

import pyrate.core.config as cf
from pyrate.core.shared import Ifg, DEM
from pyrate.core import ifgconstants as ifc
from pyrate import conv2tif, prepifg, configuration
from tests.common import manipulate_test_conf


def test_dem_and_incidence_not_converted(gamma_params):
    gp_copy = copy.deepcopy(gamma_params)
    gp_copy[cf.DEM_FILE] = None
    gp_copy[cf.APS_INCIDENCE_MAP] = None
    conv2tif.main(gp_copy)
    inc_tif = glob.glob(os.path.join(gp_copy[cf.OBS_DIR], '*inc.tif'))
    assert len(inc_tif) == 0
    dem_tif = glob.glob(os.path.join(gp_copy[cf.OBS_DIR], '*dem.tif'))
    assert len(dem_tif) == 0


def test_conv2tif_file_types(tempdir, gamma_conf):
    tdir = Path(tempdir())
    params = manipulate_test_conf(gamma_conf, tdir)
    params[cf.COH_MASK] = 1
    output_conf_file = 'conf.conf'
    output_conf = tdir.joinpath(output_conf_file)
    cf.write_config_file(params=params, output_conf_file=output_conf)
    params_s = configuration.Configuration(output_conf).__dict__
    conv2tif.main(params_s)
    ifg_files = list(Path(tdir.joinpath(params_s[cf.OUT_DIR])).glob('*_ifg.tif'))
    coh_files = list(Path(tdir.joinpath(params_s[cf.OUT_DIR])).glob('*_coh.tif'))
    dem_file = list(Path(tdir.joinpath(params_s[cf.OUT_DIR])).glob('*_dem.tif'))[0]
    # assert coherence and ifgs have correct metadata
    for i in itertools.chain(*[ifg_files, coh_files]):
        ifg = Ifg(i)
        ifg.open()
        md = ifg.meta_data
        if i.name.endswith('_ifg.tif'):
            assert md[ifc.DATA_TYPE] == ifc.ORIG
            continue
        if i.name.endswith('_coh.tif'):
            assert md[ifc.DATA_TYPE] == ifc.COH
            continue

    # assert dem has correct metadata
    dem = DEM(dem_file.as_posix())
    dem.open()
    md = dem.dataset.GetMetadata()
    assert md[ifc.DATA_TYPE] == ifc.DEM
    shutil.rmtree(tdir)


def test_tifs_placed_in_out_dir(gamma_params):
    # Test no tifs in obs dir
    tifs = glob.glob(os.path.join(gamma_params[cf.OUT_DIR], '*.tif'))
    assert len(tifs) == 0
    # Test tifs in obs dir
    conv2tif.main(gamma_params)
    tifs = glob.glob(os.path.join(gamma_params[cf.OUT_DIR], '*.tif'))
    assert len(tifs) == 18


def test_num_gamma_tifs_equals_num_unws(gamma_params):
    gtifs = conv2tif.main(gamma_params)
    # 17 unws + dem
    assert len(gtifs) == 18

    for g, _ in gtifs:   # assert all output from conv2tfi are readonly
        assert Path(g).stat().st_mode == 33060

def test_num_roipac_tifs_equals_num_unws(roipac_params):
    gtifs = conv2tif.main(roipac_params)
    # 17 unws + dem
    assert len(gtifs) == 18


def test_conversion_not_recomputed(gamma_params):

    """
    If a gtif already exists, conv2tif will return None instead
    of path to the gtif.
    """
    conv2tif.main(gamma_params)
    gtifs_and_conversion_bool = conv2tif.main(gamma_params)
    assert all([not b for gt, b in gtifs_and_conversion_bool])


def test_no_tifs_exits(gamma_params):
    with pytest.raises(Exception):
        prepifg.main(gamma_params)


def test_tifs_succeeds(gamma_params):
    conv2tif.main(gamma_params)
    prepifg.main(gamma_params)
