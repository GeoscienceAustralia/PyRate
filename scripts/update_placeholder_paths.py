import os
from pathlib import Path
PYRATE = Path(__file__).parent.parent


def remove_place_holder(file_path, root):
    root = os.path.dirname(root).replace("\\", "/")
    lines = []
    with open(file_path) as file_in:
        for line in file_in:
            line = line.replace("/absolute/path/to", root)
            lines.append(line)

    with open(file_path, "w") as f:
        for line in lines:
            f.write(line)


remove_place_holder(os.path.join(PYRATE, "sample_data", "input_parameters.conf"), PYRATE)
remove_place_holder(os.path.join(PYRATE, "sample_data", "input", "coherence_list.txt"), PYRATE)
remove_place_holder(os.path.join(PYRATE, "sample_data", "input", "headers_list.txt"), PYRATE)
remove_place_holder(os.path.join(PYRATE, "sample_data", "input", "interferogram_list.txt"), PYRATE)
remove_place_holder(PYRATE.joinpath('tests', "test_data", "small_test", "conf", "pyrate2.conf"), PYRATE)

# roipac
remove_place_holder(PYRATE.joinpath('tests', "test_data", "system", "roipac", "header_list.txt"), PYRATE)
remove_place_holder(PYRATE.joinpath('tests', "test_data", "small_test", "conf", "pyrate_roipac_test.conf"), PYRATE)
remove_place_holder(PYRATE.joinpath("tests", "test_data", "small_test", "roipac_obs", "headers_17"), PYRATE)
remove_place_holder(PYRATE.joinpath("tests", "test_data", "small_test", "roipac_obs", "ifms_17"), PYRATE)
remove_place_holder(os.path.join("tests", "test_data", "system", "roipac", "input_parameters.conf"), PYRATE)
remove_place_holder(os.path.join("tests", "test_data", "system", "roipac", "interferogram_list.txt"), PYRATE)


# gamma
remove_place_holder(PYRATE.joinpath('tests', "test_data", "small_test", "conf", "pyrate_gamma_test.conf"), PYRATE)
remove_place_holder(PYRATE.joinpath("tests", "test_data", "small_test", "gamma_obs", "headers"), PYRATE)
remove_place_holder(PYRATE.joinpath("tests", "test_data", "small_test", "gamma_obs", "ifms_17"), PYRATE)
remove_place_holder(os.path.join("tests", "test_data", "system", "gamma", "header_list.txt"), PYRATE)
remove_place_holder(PYRATE.joinpath("tests", "test_data", "system", "gamma", "input_parameters.conf"), PYRATE)
remove_place_holder(PYRATE.joinpath("tests", "test_data", "system", "gamma", "interferogram_list.txt"), PYRATE)


# geotiff
remove_place_holder(os.path.join("tests", "test_data", "system", "geotiff", "header_list.txt"), PYRATE)
remove_place_holder(os.path.join("tests", "test_data", "system", "geotiff", "input_parameters.conf"), PYRATE)
remove_place_holder(os.path.join("tests", "test_data", "system", "geotiff", "interferogram_list.txt"), PYRATE)
remove_place_holder(os.path.join("tests", "test_data", "small_test", "tif", "ifms_17"), PYRATE)
