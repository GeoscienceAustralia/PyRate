#   This Python module is part of the PyRate software package.
#
#   Copyright 2017 Geoscience Australia
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
This Python module is a Luigi wrapper for converting input data to geotiff.
"""
import luigi
import pyrate.config as config
from pyrate.tasks.utils import InputParam
from pyrate.tasks.roipac import ConvertToGeotiff as ConvertToGeotiffRoipac
from pyrate.tasks.gamma import ConvertToGeotiff as ConvertToGeotiffGamma

ROIPAC_PROCESSOR = 0
GAMMA_PROCESSOR = 1


class ConvertToGeotiff(luigi.WrapperTask):
    """
    Luigi wrapper class for both roipac and gamma geotiff conversion
    """
    processor = luigi.IntParameter(config_path=InputParam(config.PROCESSOR))

    def requires(self):
        if self.processor == int(ROIPAC_PROCESSOR):
            return [ConvertToGeotiffRoipac()]
        elif self.processor == int(GAMMA_PROCESSOR):
            return [ConvertToGeotiffGamma()]
        else:
            raise luigi.parameter.ParameterException(
                'invalid value for parameter {}'.format(config.PROCESSOR))
