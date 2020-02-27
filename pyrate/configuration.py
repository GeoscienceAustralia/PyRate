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
import pathlib
import re
from configparser import ConfigParser

from constants import NO_OF_PARALLEL_PROCESSES
from core.user_experience import break_number_into_factors
from default_parameters import PYRATE_DEFAULT_CONFIGRATION


def set_parameter_value(data_type, input_value, default_value, required, input_name):
    """

    Args:
      data_type: param input_value:
      default_value: param required:
      input_name: 
      input_value: 
      required: 

    Returns:

    """
    if len(input_value) < 1:
        input_value = None
        if required:
            raise ValueError("A required parameter is missing value in input configuration file: " + str(input_name))
    try:
        if input_value is not None:
            if str(data_type) in "path":
                return pathlib.Path(input_value)
            return data_type(input_value)
    except ValueError:
        raise ValueError("A given parameter: "+str(input_name)+" has invalid data type: " +str(input_value)+" in input configuration file.")
    return default_value


def validate_parameter_value(input_name, input_value, min_value=None, max_value=None, possible_values=None):
    """
    Args:
      input_name: param input_value:
      min_value: Default value = None)
      max_value: Default value = None)
      possible_values: Default value = None)
      input_value: 

    Returns:
    """
    if isinstance(input_value, pathlib.PurePath):
        if input_name in "outdir":
            input_value = input_value.parents[0]
        if not pathlib.Path.exists(input_value):
            raise ValueError("Given path: " + str(input_value) + " dose not exist.")
    if input_value is not None:
        if min_value is not None:
            if input_value < min_value:
                raise ValueError(
                    "Invalid value for "
                    + str(input_name)
                    + " supplied: "
                    + str(input_value)
                    + ". Please provided a valid value greater than "
                    + str(min_value)
                    + "."
                )
    if input_value is not None:
        if max_value is not None:
            if input_value > max_value:
                raise ValueError(
                    "Invalid value for "
                    + str(input_name)
                    + " supplied: "
                    + str(input_value)
                    + ". Please provided a valid value less than "
                    + str(max_value)
                    + "."
                )

    if possible_values is not None:
        if input_value not in possible_values:
            raise ValueError(
                "Invalid value for "
                + str(input_name)
                + " supplied: "
                + str(input_value)
                + ". Please provided a valid value from with in: "
                + str(possible_values)
                + "."
            )
    return True


def validate_file_list_values(fileList, no_of_epochs):
    """
    Args:
      fileList:
      no_of_epochs:
    Returns:
    """
    if fileList is None:
        raise ValueError("No value supplied for input file list: " + str(fileList))

    for file_path_str in fileList.read_text().split("\n"):
        # ignore empty lines in file
        if len(file_path_str) > 1:

            if not pathlib.Path.exists(pathlib.Path(file_path_str)):
                raise ValueError("Give file name: " + str(file_path_str)+" does not exist.")
            else:
                matches = re.findall("(\d{8})", file_path_str)
                if len(matches) < no_of_epochs:
                    matches = re.findall("(\d{6})", file_path_str)
                    if len(matches) < no_of_epochs:
                        raise ValueError(
                            "For the given file name: "
                            + str(file_path_str)
                            + " the number of epochs in file names are less the required number:"
                            + str(no_of_epochs)
                        )


class MultiplePaths:
    def __init__(self, outDir, baseFilePath, ifglksx=0, ifgcropopt=0):
        """

        Args:
          baseDir: param baseName:
          ifglksx: param ifgcropopt:

        Returns:

        """
        if ".tif" in str(baseFilePath):
            self.unwrapped_path = None
            self.converted_path = str(baseFilePath)
            self.sampled_path = self.converted_path.split(".tif")[0] + "_" + str(ifglksx) + "rlks_" + str(
                ifgcropopt) + "cr.tif"
        else:
            self.unwrapped_path = str(baseFilePath)
            baseName = baseFilePath.name.split(".")[0] + "_" + baseFilePath.name.split(".")[1] + ".tif"

            self.converted_path = str(outDir / baseName)
            self.sampled_path = self.converted_path.split(".tif")[0] + "_" + str(ifglksx) + "rlks_" + str(
                ifgcropopt) + "cr.tif"

    def __str__(self):
        return (
                """
    unwrapped_path = """
                + self.unwrapped_path
                + """ 
converted_path = """
                + self.converted_path
                + """ 
sampled_path = """
                + self.sampled_path
                + """    
"""
        )


class Configuration:
    def __init__(self, config_file_path):
        """
        Args:
          config_file_path:

        Returns:

        """
        parser = ConfigParser()
        parser.optionxform = str
        # mimic header to fulfil the requirement for configparser
        with open(config_file_path) as stream:
            parser.read_string("[root]\n" + stream.read())

        for key, value in parser._sections["root"].items():
            self.__dict__[key] = value

        # Validate required parameters exist.
        if not set(PYRATE_DEFAULT_CONFIGRATION).issubset(self.__dict__):
            raise ValueError(
                "Required configuration parameters: "
                + str(set(PYRATE_DEFAULT_CONFIGRATION).difference(self.__dict__))
                + " are missing from input config file."
            )

        # handle control parameters
        for parameter_name in PYRATE_DEFAULT_CONFIGRATION:
            self.__dict__[parameter_name] = set_parameter_value(
                PYRATE_DEFAULT_CONFIGRATION[parameter_name]["DataType"],
                self.__dict__[parameter_name],
                PYRATE_DEFAULT_CONFIGRATION[parameter_name]["DefaultValue"],
                PYRATE_DEFAULT_CONFIGRATION[parameter_name]["Required"],
                parameter_name,
            )
            validate_parameter_value(
                parameter_name,
                self.__dict__[parameter_name],
                PYRATE_DEFAULT_CONFIGRATION[parameter_name]["MinValue"],
                PYRATE_DEFAULT_CONFIGRATION[parameter_name]["MaxValue"],
                PYRATE_DEFAULT_CONFIGRATION[parameter_name]["PossibleValues"],
            )

        # bespoke parameter validation
        if self.refchipsize % 2 != 1:
            if self.refchipsize - 1 > 1:
                # Configuration parameters refchipsize must be odd
                # values too large (>101) will slow down the process without significant gains in results.
                self.refchipsize = self.refchipsize - 1

        self.rows, self.cols = [int(no) for no in break_number_into_factors(NO_OF_PARALLEL_PROCESSES)]

        # create a temporary directory
        self.tmpdir = self.outdir / "tmpdir"
        self.tmpdir.mkdir(parents=True, exist_ok=True)
        self.tmpdir = str(self.tmpdir)

        # var no longer used
        self.APS_ELEVATION_EXT = None
        self.APS_INCIDENCE_EXT = None
        self.apscorrect = 0
        self.apsmethod = 0
        self.elevationmap = None
        self.incidencemap = None

        # define parallel processes that will run
        self.NUMEXPR_MAX_THREADS = str(NO_OF_PARALLEL_PROCESSES)

        self.coherence_file_paths = None
        if self.cohfilelist is not None:
            validate_file_list_values(self.cohfilelist, 1)
            self.coherence_file_paths = []
            for coherence_file_path in self.cohfilelist.read_text().split("\n"):
                # ignore empty lines in file
                if len(coherence_file_path) > 1:
                    coherence_file_path = pathlib.Path(coherence_file_path)
                    self.coherence_file_paths.append(
                        MultiplePaths(self.outdir, coherence_file_path, self.ifglksx, self.ifgcropopt)
                    )

        # Validate file names supplied in list exist and contain correct epochs in file names
        validate_file_list_values(self.slcfilelist, 1)

        self.header_file_paths = []
        for header_file_path in self.slcfilelist.read_text().split("\n"):
            # ignore empty lines in file
            if len(header_file_path) > 1:
                header_file_path = pathlib.Path(header_file_path)
                self.header_file_paths.append(MultiplePaths(self.outdir, header_file_path, self.ifglksx, self.ifgcropopt))

        validate_file_list_values(self.ifgfilelist, 2)

        self.interferogram_files = []
        for interferogram_file_path in self.ifgfilelist.read_text().split("\n"):
            # ignore empty lines in file
            if len(interferogram_file_path) > 1:
                interferogram_file_path = pathlib.Path(interferogram_file_path)
                self.interferogram_files.append(MultiplePaths(self.outdir, interferogram_file_path, self.ifglksx, self.ifgcropopt))

        self.dem_file = MultiplePaths(self.outdir, self.demfile, self.ifglksx, self.ifgcropopt)

        # backward compatibility for string paths
        for key in self.__dict__:
            if isinstance(self.__dict__[key], pathlib.PurePath):
                self.__dict__[key] = str(self.__dict__[key])


if __name__ == "__main__":
    config_file_path = "C:\\Users\\sheec\\Desktop\\Projects\\PyRate\\sample_data\\input_parameters.conf"
    config = Configuration(config_file_path)
    for interPath in config.interferogram_files:
        print(interPath)
