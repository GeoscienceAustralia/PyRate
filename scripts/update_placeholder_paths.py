import os
from pathlib import Path

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


PYRATE = Path(__file__).parent.parent
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
input_parameters_file = os.path.join(root, "sample_data", "input_parameters.conf")
remove_place_holder(input_parameters_file, root)
coherence_list_file = os.path.join(root, "sample_data", "input", "coherence_list.txt")
remove_place_holder(coherence_list_file, root)
headers_list = os.path.join(root, "sample_data", "input", "headers_list.txt")
remove_place_holder(headers_list, root)
interferogram_list = os.path.join(root, "sample_data", "input", "interferogram_list.txt")
remove_place_holder(interferogram_list, root)
pyrate2_conf = PYRATE.joinpath('tests', "test_data", "small_test", "conf", "pyrate2.conf")
remove_place_holder(pyrate2_conf, PYRATE)

# roipac
roipac_headers = PYRATE.joinpath('tests', "test_data", "system", "roipac", "header_list.txt")
remove_place_holder(roipac_headers, PYRATE)
pyrate_roipac_test_conf = PYRATE.joinpath('tests', "test_data", "small_test", "conf", "pyrate_roipac_test.conf")
remove_place_holder(pyrate_roipac_test_conf, PYRATE)
roipac_headers = PYRATE.joinpath("tests", "test_data", "small_test", "roipac_obs", "headers_17")
remove_place_holder(roipac_headers, PYRATE)
roipac_ifms = PYRATE.joinpath("tests", "test_data", "small_test", "roipac_obs", "ifms_17")
remove_place_holder(roipac_ifms, PYRATE)


# gamma
pyrate_gamma_test_conf = PYRATE.joinpath('tests', "test_data", "small_test", "conf", "pyrate_gamma_test.conf")
remove_place_holder(pyrate_gamma_test_conf, PYRATE)
gamma_headers = PYRATE.joinpath("tests", "test_data", "small_test", "gamma_obs", "headers")
remove_place_holder(gamma_headers, PYRATE)
gamma_obs = PYRATE.joinpath("tests", "test_data", "small_test", "gamma_obs", "ifms_17")
remove_place_holder(gamma_obs, PYRATE)
