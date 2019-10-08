import re
import configparser


def pprint_config_parser(config_parser):
    for section_name in config_parser.sections():
        # print('Section:', section_name)
        # print('Parameters:', config_parser.options(section_name))
        for name, value in config_parser.items(section_name):
            print('  {} : {}'.format(name, value))
        print()


def read_raw_header_file(header_path):
    file_content = "[root]\n"
    with open(header_path) as f:
        for line in f.readlines():
            if ":" in line:
                file_content = file_content + line

    header = configparser.RawConfigParser()
    header.read_string(file_content)

    return header


def find_all_substring(pattern, string):
    match = re.search(pattern, str(string))
    if match:
        return match[0]
    return None


def find_substring(pattern, string):
    match = re.search(pattern, str(string))
    if match:
        return match[0]
    return None


def sort_headers(header_paths, interferograms_path):

    headers = {}
    interferogram_epochs = find_all_substring(r'\d{8}', interferograms_path)

    if len(interferogram_epochs) < 2:
        Exception("At least two epoches are required in the  interferogram file name:"+str(interferograms_path))

    if len(header_paths) < 2:
        Exception("At least two header files are required.")

    for header_path in header_paths:
        epoch = int(find_substring(r'\d{8}', header_path))
        if epoch and str(epoch) in interferogram_epochs:
            headers[epoch] = header_path

    epochs = headers.keys()

    return headers[max(epochs)], headers[min(epochs)]


def format_header_value(raw_value, pattern, date_type, default=None):

    match = re.match(pattern, raw_value)
    if match:
        return date_type(match[0])
    return default
