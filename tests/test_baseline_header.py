import pytest
from pyrate.core.gamma import GammaException, combine_headers


def test_combine_header():
    hdr0 = hdr1 = dem_hdr = {}
    base_hdr = []
    with pytest.raises(GammaException):
        combine_headers(hdr0, hdr1, dem_hdr, base_hdr)
