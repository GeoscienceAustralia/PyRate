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

        # create a gdal data set
        # create a coherence mask dataset
        # create a artificial resultant data set

        # use the gdal_python.coherence_masking to find the resultant mask data set

        # compare the artificial and actual resultant data set
        self.assertTrue(True)



if __name__ == '__main__':
    unittest.main()
