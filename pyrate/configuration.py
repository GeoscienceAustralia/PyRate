#   This Python module is part of the PyRate software package.
#
#   Copyright 2020 Geoscience Australia
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
"""
This Python module contains utilities to validate user input parameters
parsed in a PyRate configuration file.
"""
from configparser import ConfigParser
from pathlib import Path, PurePath
from pyrate.constants import NO_OF_PARALLEL_PROCESSES
from pyrate.default_parameters import PYRATE_DEFAULT_CONFIGURATION
from pyrate.core.algorithm import factorise_integer
from pyrate.core.shared import extract_epochs_from_filename, InputTypes
from pyrate.core.config import parse_namelist, ConfigException


def set_parameter_value(data_type, input_value, default_value, required, input_name):
    if len(input_value) < 1:
        input_value = None
        if required:  # pragma: no cover
            raise ValueError("A required parameter is missing value in input configuration file: " + str(input_name))

    if input_value is not None:
        if str(data_type) in "path":
            return Path(input_value)
        return data_type(input_value)
    return default_value


def validate_parameter_value(input_name, input_value, min_value=None, max_value=None, possible_values=None):
    if isinstance(input_value, PurePath):
        if not Path.exists(input_value):  # pragma: no cover
            raise ValueError("Given path: " + str(input_value) + " does not exist.")
    if input_value is not None:
        if min_value is not None:
            if input_value < min_value:  # pragma: no cover
                raise ValueError(
                    "Invalid value for " + str(input_name) + " supplied: " + str(input_value) + ". Provide a value greater than or equal to " + str(min_value) + ".")
    if input_value is not None:
        if max_value is not None:
            if input_value > max_value:  # pragma: no cover
                raise ValueError(
                    "Invalid value for " + str(input_name) + " supplied: " + str(input_value) + ". Provide a value less than or equal to " + str(max_value) + ".")

    if possible_values is not None:
        if input_value not in possible_values:  # pragma: no cover
            raise ValueError(
                "Invalid value for " + str(input_name) + " supplied: " + str(input_value) + ". Provide a value from: " + str(possible_values) + ".")
    return True


def validate_file_list_values(file_list, no_of_epochs):
    if file_list is None:  # pragma: no cover
        raise ValueError("No value supplied for input file list: " + str(file_list))

    files = parse_namelist(file_list)

    for f in files:
        if not Path(f).exists():  # pragma: no cover
            raise ConfigException(f"{f} does not exist")
        else:
            matches = extract_epochs_from_filename(filename_with_epochs=f)
            if len(matches) < no_of_epochs:  # pragma: no cover
                raise ConfigException(f"the number of epochs in {f} names are less the required number: {no_of_epochs}")


class MultiplePaths:
    def __init__(self, out_dir: str, file_name: str, ifglksx: int = 1, ifgcropopt: int = 1,
                 input_type: InputTypes = InputTypes.IFG):
        self.input_type = input_type
        b = Path(file_name)
        if b.suffix == ".tif":
            self.unwrapped_path = None
            converted_path = b  # original file
            self.sampled_path = Path(out_dir).joinpath(
                b.stem + '_' + str(ifglksx) + "rlks_" + str(ifgcropopt) + "cr.tif").as_posix()
        else:
            self.unwrapped_path = b.as_posix()
            converted_path = Path(out_dir).joinpath(
                b.stem.split('.')[0] + '_' + b.suffix[1:] + input_type.value).with_suffix('.tif')
            self.sampled_path = converted_path.with_name(
                converted_path.stem + '_' + str(ifglksx) + "rlks_" + str(ifgcropopt) + "cr.tif").as_posix()
        self.converted_path = converted_path.as_posix()

    def __str__(self):  # pragma: no cover
        return"""
            unwrapped_path = """ + self.unwrapped_path+""" 
            converted_path = """ + self.converted_path+""" 
            sampled_path = """ + self.sampled_path+"""    
            """


class Configuration:

    def __init__(self, config_file_path):

        parser = ConfigParser()
        parser.optionxform = str
        # mimic header to fulfil the requirement for configparser
        with open(config_file_path) as stream:
            parser.read_string("[root]\n" + stream.read())

        for key, value in parser._sections["root"].items():
            self.__dict__[key] = value

        # make output path, if not provided will error
        Path(self.outdir).mkdir(exist_ok=True, parents=True)

        # Validate required parameters exist.

        required = {k for k, v in PYRATE_DEFAULT_CONFIGURATION.items() if v['Required']}

        if not required.issubset(self.__dict__):  # pragma: no cover
            raise ValueError("Required configuration parameters: " + str(
                required.difference(self.__dict__)) + " are missing from input config file.")

        # handle control parameters
        for parameter_name in PYRATE_DEFAULT_CONFIGURATION:
            param_value = self.__dict__[parameter_name] if parameter_name in required or \
                                                           parameter_name in self.__dict__ else ''

            self.__dict__[parameter_name] = set_parameter_value(PYRATE_DEFAULT_CONFIGURATION[parameter_name]["DataType"],
                                                                param_value,
                                                                PYRATE_DEFAULT_CONFIGURATION[parameter_name]["DefaultValue"],
                                                                PYRATE_DEFAULT_CONFIGURATION[parameter_name]["Required"],
                                                                parameter_name)
            validate_parameter_value(parameter_name, self.__dict__[parameter_name],
                                     PYRATE_DEFAULT_CONFIGURATION[parameter_name]["MinValue"],
                                     PYRATE_DEFAULT_CONFIGURATION[parameter_name]["MaxValue"],
                                     PYRATE_DEFAULT_CONFIGURATION[parameter_name]["PossibleValues"])

        # bespoke parameter validation
        if self.refchipsize % 2 != 1:  # pragma: no cover
            if self.refchipsize - 1 > 1:
                # Configuration parameters refchipsize must be odd
                # values too large (>101) will slow down the process without significant gains in results.
                self.refchipsize = self.refchipsize - 1

        # calculate rows and cols if not supplied
        if hasattr(self, 'rows') and hasattr(self, 'cols'):
            self.rows, self.cols = int(self.rows), int(self.cols)
        else:
            self.rows, self.cols = [int(num) for num in factorise_integer(NO_OF_PARALLEL_PROCESSES)]

        # create a temporary directory if not supplied
        if not hasattr(self, 'tmpdir'):
            self.tmpdir = Path(self.outdir).joinpath("tmpdir")
        else:
            self.tmpdir = Path(self.tmpdir)

        self.tmpdir.mkdir(parents=True, exist_ok=True)

        # var no longer used
        self.APS_ELEVATION_EXT = None
        self.APS_INCIDENCE_EXT = None
        self.apscorrect = 0
        self.apsmethod = 0
        self.elevationmap = None
        self.incidencemap = None

        # define parallel processes that will run
        self.NUMEXPR_MAX_THREADS = str(NO_OF_PARALLEL_PROCESSES)

        # Validate file names supplied in list exist and contain correct epochs in file names
        if self.cohfilelist is not None:
            # if self.processor != 0:  # not roipac
            validate_file_list_values(self.cohfilelist, 1)
            self.coherence_file_paths = self.__get_files_from_attr('cohfilelist', input_type=InputTypes.COH)

        self.header_file_paths = self.__get_files_from_attr('hdrfilelist', input_type=InputTypes.HEADER)

        self.interferogram_files = self.__get_files_from_attr('ifgfilelist')

        self.dem_file = MultiplePaths(self.outdir, self.demfile, self.ifglksx, self.ifgcropopt,
                                      input_type=InputTypes.DEM)

        # backward compatibility for string paths
        for key in self.__dict__:
            if isinstance(self.__dict__[key], PurePath):
                self.__dict__[key] = str(self.__dict__[key])

    def __get_files_from_attr(self, attr, input_type=InputTypes.IFG):
        val = self.__getattribute__(attr)
        files = parse_namelist(val)
        return [MultiplePaths(self.outdir, p, self.ifglksx, self.ifgcropopt, input_type=input_type) for p in files]
