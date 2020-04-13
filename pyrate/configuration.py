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
from configparser import ConfigParser
from pathlib import Path, PurePath
import re
from pyrate.constants import NO_OF_PARALLEL_PROCESSES
from pyrate.default_parameters import PYRATE_DEFAULT_CONFIGRATION
from pyrate.core.user_experience import break_number_into_factors


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
        if input_name in "outdir":
            input_value = input_value.parents[0]
        if not Path.exists(input_value):  # pragma: no cover
            raise ValueError("Given path: " + str(input_value) + " dose not exist.")
    if input_value is not None:
        if min_value is not None:
            if input_value < min_value:  # pragma: no cover
                raise ValueError(
                    "Invalid value for " + str(input_name) + " supplied: " + str(input_value) + ". Please provided a valid value greater than " + str(min_value) + ".")
    if input_value is not None:
        if max_value is not None:
            if input_value > max_value:  # pragma: no cover
                raise ValueError(
                    "Invalid value for " + str(input_name) + " supplied: " + str(input_value) + ". Please provided a valid value less than " + str(max_value) + ".")

    if possible_values is not None:
        if input_value not in possible_values:  # pragma: no cover
            raise ValueError(
                "Invalid value for " + str(input_name) + " supplied: " + str(input_value) + ". Please provided a valid value from with in: " + str(possible_values) + ".")
    return True


def validate_file_list_values(dir, fileList, no_of_epochs):
    if dir is None:  # pragma: no cover
        raise ValueError("No value supplied for input directory: " + str(dir))
    if fileList is None:  # pragma: no cover
        raise ValueError("No value supplied for input file list: " + str(fileList))

    for path_str in fileList.read_text().split('\n'):
        # ignore empty lines in file
        if len(path_str) > 1:
            if not Path(dir).joinpath(path_str).exists():  # pragma: no cover
                raise ValueError("Give file name: " + str(path_str) + " dose not exist at: " + str(dir))
            else:
                matches = re.findall(r"(\d{8})", str(dir / path_str))
                if len(matches) < no_of_epochs:  # pragma: no cover
                    raise ValueError("For the given file name: " + str(path_str) + " in: " + str(
                        dir) + " the number of epochs in file names are less the required number:" + str(
                        no_of_epochs))


class MultiplePaths:
    def __init__(self, baseDir, baseName, ifglksx=0, ifgcropopt=0):
        if ".tif" in baseName:
            self.unwrapped_path = None
            self.converted_path = Path(baseDir).joinpath(baseName).as_posix()
            self.sampled_path = self.converted_path.split(".tif")[0] + "_" + str(ifglksx) + "rlks_" + str(ifgcropopt) + "cr.tif"
        else:
            self.unwrapped_path = Path(baseDir).joinpath(baseName).as_posix()
            baseName = baseName.split(".")[0] + "_" + baseName.split(".")[1] + ".tif"

            self.converted_path = Path(baseDir).joinpath(baseName).as_posix()
            self.sampled_path = self.converted_path.split(".tif")[0] + "_" + str(ifglksx) + "rlks_" + str(ifgcropopt) + "cr.tif"

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

        # Validate required parameters exist.
        if not set(PYRATE_DEFAULT_CONFIGRATION).issubset(self.__dict__):  # pragma: no cover
            raise ValueError("Required configuration parameters: " + str(
                set(PYRATE_DEFAULT_CONFIGRATION).difference(self.__dict__)) + " are missing from input config file.")

        # handle control parameters
        for parameter_name in PYRATE_DEFAULT_CONFIGRATION:
            self.__dict__[parameter_name] = set_parameter_value(PYRATE_DEFAULT_CONFIGRATION[parameter_name]["DataType"],
                                                                self.__dict__[parameter_name],
                                                                PYRATE_DEFAULT_CONFIGRATION[parameter_name]["DefaultValue"],
                                                                PYRATE_DEFAULT_CONFIGRATION[parameter_name]["Required"],
                                                                parameter_name)
            validate_parameter_value(parameter_name, self.__dict__[parameter_name],
                                     PYRATE_DEFAULT_CONFIGRATION[parameter_name]["MinValue"],
                                     PYRATE_DEFAULT_CONFIGRATION[parameter_name]["MaxValue"],
                                     PYRATE_DEFAULT_CONFIGRATION[parameter_name]["PossibleValues"])

        # bespoke parameter validation

        if self.refchipsize % 2 != 1:  # pragma: no cover
            raise ValueError("Configuration parameters refchipsize must be odd: " + str(self.refchipsize))

        # calculate rows and cols if not supplied
        if any([~hasattr(self, 'rows'), ~hasattr(self, 'cols')]):
            self.rows, self.cols = [int(no) for no in break_number_into_factors(NO_OF_PARALLEL_PROCESSES)]

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

        if self.processor != 0:   # not roipac
            validate_file_list_values(self.obsdir, self.ifgfilelist, 2)
            validate_file_list_values(self.slcFileDir, self.slcfilelist, 1)

        self.coherence_file_paths = None
        if self.cohfiledir is not None and self.cohfilelist is not None:
            # FIXME: tests don't reach here due to configs not having cohfilelist
            if self.processor != 0:  # not roipac
                validate_file_list_values(self.cohfiledir, self.cohfilelist, 1)
            self.coherence_file_paths = []
            for path_str in self.cohfilelist.read_text().split('\n'):
                # ignore empty lines in file
                if len(path_str) > 1:
                    self.coherence_file_paths.append(MultiplePaths(self.cohfiledir, path_str, self.ifglksx, self.ifgcropopt))

        self.header_file_paths = []
        for path_str in self.slcfilelist.read_text().split('\n'):
            # ignore empty lines in file
            if len(path_str) > 1:
                self.header_file_paths.append(MultiplePaths(self.slcFileDir, path_str, self.ifglksx, self.ifgcropopt))

        self.__get_interferrograms_files()

        self.dem_file = MultiplePaths(self.demfile.parents[0], self.demfile.name, self.ifglksx, self.ifgcropopt)

        # backward compatibility for string paths
        for key in self.__dict__:
            if isinstance(self.__dict__[key], PurePath):
                self.__dict__[key] = str(self.__dict__[key])

    def __get_interferrograms_files(self):
        self.interferogram_files = []
        for path_str in self.ifgfilelist.read_text().split('\n'):
            # ignore empty lines in file
            if len(path_str) > 1:
                self.interferogram_files.append(MultiplePaths(self.outdir, path_str, self.ifglksx, self.ifgcropopt))
