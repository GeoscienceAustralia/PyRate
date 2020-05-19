import os
import pytest
import glob
import copy

import pyrate.core.config as cf
from pyrate import conv2tif, prepifg


def test_dem_and_incidence_not_converted(gamma_params):
    gp_copy = copy.deepcopy(gamma_params)
    gp_copy[cf.DEM_FILE] = None
    gp_copy[cf.APS_INCIDENCE_MAP] = None
    conv2tif.main(gp_copy)
    inc_tif = glob.glob(os.path.join(gp_copy[cf.OBS_DIR], '*inc.tif'))
    assert len(inc_tif) == 0
    dem_tif = glob.glob(os.path.join(gp_copy[cf.OBS_DIR], '*dem.tif'))
    assert len(dem_tif) == 0


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
