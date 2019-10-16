import os
import unittest
import pytest
import glob
import copy

import pyrate.core.config as cf
from pyrate import converttogeotif, prepifg
from pyrate.core import gdal_python


class CoherenceMaskingTest(unittest.TestCase):


    def test_coherence_files_not_converted(self):

        # create a sample gdal dataset
        # create a coherence mask dataset
        # create a artificial masked dataset

        # use the gdal_python.coherence_masking to find the actual mask dataset

        # compare the artificial masked and actual masked datasets
        self.assertTrue(True)



if __name__ == '__main__':
    unittest.main()
