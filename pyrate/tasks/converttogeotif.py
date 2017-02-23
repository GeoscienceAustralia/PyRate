"""
Convert to geotif wrapper using luigi
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
