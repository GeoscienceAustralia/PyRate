from pathlib import Path
import os

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
pyrate2_conf = PYRATE.joinpath('tests', "test_data", "small_test", "conf", "pyrate2.conf")
remove_place_holder(pyrate2_conf, PYRATE)
roipac_headers = PYRATE.joinpath('tests', "test_data", "system", "roipac", "header_list.txt")
remove_place_holder(roipac_headers, PYRATE)
