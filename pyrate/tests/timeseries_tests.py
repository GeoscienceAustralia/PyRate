'''
Collection of tests for validating PyRate's time series analysis code.
Author: Vanessa Newey
'''

import unittest
from numpy import nan, mean, std, isnan, asarray, diff, where
from numpy.testing import assert_array_almost_equal
from datetime import date, timedelta

from common import sydney_data_setup, MockIfg
from pyrate.timeseries import time_series 
from pyrate.config import ConfigException

from pyrate.config import TS_PTHRESH


def default_params():
    return { TS_PTHRESH : 10 }


class TimeSeriesTests(unittest.TestCase):
    '''Verifies error checking capabilities of the time_series function'''


    def setUp(self):
        self.ifgs = sydney_data_setup()



    def test_time_series(self):
        '''
        Checks that the code works the same as the pirate MatLab code
        '''
        params=default_params()
        
        tsincr,tscum,tsvel = time_series(self.ifgs, pthresh=params[TS_PTHRESH])
        expected= asarray([-10.7925865, -2.24628509, -11.07861038, -7.68240018,
        -8.66006356, -4.17442079, -10.34031155, -2.37418019,
        -5.84745442, -7.25372609, -8.10866434, -9.39923284])
        
        assert_array_almost_equal(tscum[10,10,:],expected)

        
    def test_time_series_unit(self):

        imaster = asarray([1,1,2,2,3,3,4,5])
        islave = asarray([2,4,3,4,5,6,6,6])
        timeseries = asarray([0.0,0.1,0.6,0.8,1.1,1.3])
        phase = asarray([0.5,4,2.5,3.5,2.5,3.5,2.5,1])
        
        now = date.today()
        
        dates =[now + timedelta(days=(t*365.25))  for t in timeseries]
        dates.sort()
        master =[dates[m_num -1] for m_num in imaster]
        slave =[dates[s_num -1] for s_num in islave]
        
        self.ifgs = [SinglePixelIfg(m,s,p) for m,s,p in zip(master,slave,phase)]
        params = { TS_PTHRESH : 0 }
        tsincr,tscum,tsvel = time_series(self.ifgs, pthresh=params[TS_PTHRESH])
        expected= asarray([[[ 0.50,  3.0,  4.0,  5.5,  6.5]]])
        assert_array_almost_equal(tscum,expected)


class SinglePixelIfg(object):

    def __init__(self,master,slave,phase):
        self.phase_data = asarray([phase])  
        self.master =master
        self.slave =slave
        self.nrows =1
        self.ncols =1      
        
    def convert_to_nans(self, val=0):
        '''
        Converts given values in phase data to NaNs
        val - value to convert, default is 0
        '''
        self.phase_data = where(self.phase_data == val, nan, self.phase_data)
        self.nan_converted = True 
        
            
if __name__ == "__main__":
    unittest.main()