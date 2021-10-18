import re
from pathlib import Path
import pytest
from pyrate.constants import twelve_digits_pattern, sixteen_digits_pattern
from pyrate.configuration import parse_namelist
from tests.common import IFMS16, GAMMA_SML_TEST_DIR


@pytest.fixture
def six_digits_filenames():
    return IFMS16


@pytest.fixture
def eight_digits_filenames():
    files = list(parse_namelist(Path(GAMMA_SML_TEST_DIR).joinpath('ifms_17').as_posix()))
    return files


@pytest.mark.parametrize(
    "regex_pattern,extracted_string_length",
    [
        (twelve_digits_pattern, 12+1),  # +1 due to the internal joining `-`
        (sixteen_digits_pattern, 16+1)
     ]
)
def test_file_patterns(regex_pattern, extracted_string_length, six_digits_filenames, eight_digits_filenames):

    if regex_pattern == twelve_digits_pattern:
        for f in six_digits_filenames:
            m = re.search(twelve_digits_pattern, Path(f).stem).group()
            assert len(m) == extracted_string_length
    if regex_pattern == sixteen_digits_pattern:
        for f in eight_digits_filenames:
            m = re.search(sixteen_digits_pattern, Path(f).stem).group()
            assert len(m) == extracted_string_length

