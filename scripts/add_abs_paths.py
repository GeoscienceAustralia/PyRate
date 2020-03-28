import os


def remove_place_holder(file_path):
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    file_path = os.path.join(root, file_path)
    root = os.path.dirname(root).replace("\\", "/")
    lines = []
    with open(file_path) as file_in:
        for line in file_in:
            line = line.replace("/absolute/path/to", root)
            lines.append(line)

    with open(file_path, "w") as f:
        for line in lines:
            f.write(line)


remove_place_holder(os.path.join("sample_data", "input_parameters.conf"))
remove_place_holder(os.path.join("sample_data", "input", "coherence_list.txt"))
remove_place_holder(os.path.join("sample_data", "input", "headers_list.txt"))
remove_place_holder(os.path.join("sample_data", "input", "interferogram_list.txt"))
remove_place_holder(os.path.join("tests", "test_data", "small_test", "conf", "pyrate1.conf"))
remove_place_holder(os.path.join("tests", "test_data", "small_test", "conf", "pyrate2.conf"))
remove_place_holder(os.path.join("tests", "test_data", "small_test", "conf", "pyrate_gamma_test.conf"))
remove_place_holder(os.path.join("tests", "test_data", "small_test", "conf", "pyrate_roipac_test.conf"))
remove_place_holder(os.path.join("tests", "test_data", "small_test", "gamma_obs", "headers"))
remove_place_holder(os.path.join("tests", "test_data", "small_test", "gamma_obs", "ifms_17"))
remove_place_holder(os.path.join("tests", "test_data", "small_test", "roipac_obs", "ifms_17"))
remove_place_holder(os.path.join("tests", "test_data", "small_test", "roipac_obs", "headers_17"))
remove_place_holder(os.path.join("tests", "test_data", "small_test", "tif", "ifms_17"))
remove_place_holder(os.path.join("tests", "test_data", "system", "gamma", "header_list.txt"))
remove_place_holder(os.path.join("tests", "test_data", "system", "gamma", "input_parameters.conf"))
remove_place_holder(os.path.join("tests", "test_data", "system", "gamma", "interferogram_list.txt"))
remove_place_holder(os.path.join("tests", "test_data", "system", "geotiff", "header_list.txt"))
remove_place_holder(os.path.join("tests", "test_data", "system", "geotiff", "input_parameters.conf"))
remove_place_holder(os.path.join("tests", "test_data", "system", "geotiff", "interferogram_list.txt"))
remove_place_holder(os.path.join("tests", "test_data", "system", "roipac", "header_list.txt"))
remove_place_holder(os.path.join("tests", "test_data", "system", "roipac", "input_parameters.conf"))
remove_place_holder(os.path.join("tests", "test_data", "system", "roipac", "interferogram_list.txt"))
