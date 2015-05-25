'''
Tests for the PyRate's Variance/Covariance matrix functionality.

Created on 14/03/2013

.. codeauthor:: Ben Davies
'''

import unittest
from numpy import array
from numpy.testing import assert_array_almost_equal

from pyrate.vcm import cvd, get_vcmt
from pyrate.tests.common import sydney5_mock_ifgs, sydney5_ifgs



class TestCVD(unittest.TestCase):
    def test_cvd(self):
        ifgs = sydney5_ifgs()

        for i in ifgs:
            i.open()

        if bool((i.phase_data == 0).all()) is True:
            raise Exception("All zero")

        maxvar, alpha = cvd(i)
        assert maxvar is not None
        assert alpha is not None
        print "maxvar: %s, alpha: %s" % (maxvar, alpha)
        print "\n"


    def test_get_vcmt(self):
        ifgs = sydney5_mock_ifgs(5,9)
        maxvar = [8.486, 12.925, 6.313, 0.788, 0.649 ]

        # from Octave with MG's edits to Hua's bug
        exp = array([[8.486, 5.23645, 0.0, 0.0, 0.0],
            [5.23645, 12.925,  4.51651,  1.59569,  0.0],
            [0.0, 4.51651, 6.313, 1.1152, 0.0],
            [0.0, 1.59569, 1.1152, 0.788, -0.35757],
            [0.0, 0.0, 0.0, -0.35757, 0.649]])

        act = get_vcmt(ifgs, maxvar)
        assert_array_almost_equal(act, exp, decimal=3)



if __name__ == "__main__":
    unittest.main()
