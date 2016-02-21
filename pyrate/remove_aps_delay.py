__author__ = 'Sudipta Basak'
__date_created__ = '4/02/16'

import PyAPS as pa
import numpy as np
import sys
import os
from pyrate.scripts import run_pyrate
from pyrate.tests.common import sydney_data_setup
from pyrate.tests.common import SYD_TEST_DEM_UNW
from pyrate.scripts import run_prepifg
import re
PTN = re.compile(r'\d{6}')
PYRATEPATH = os.environ['PYRATEPATH']
ECMWF_DIR = os.path.join(PYRATEPATH, 'ECMWF')
ECMWF_PRE = 'ERA-Int_'
ECMWF_EXT = '_12.grib'


if __name__ == "__main__":
    sys.argv = ['run_pyrate.py', 'pyrate.conf']
    run_prepifg.main()

    sys.argv = ['run_pyrate.py', 'pyrate.conf']
    base_ifg_paths, dest_paths, params = run_pyrate.get_ifg_paths()
    ifgs = sydney_data_setup(datafiles=dest_paths)

    list_of_dates_for_grb_download = []

    for i in ifgs[:1]:  # demo for only one ifg
        # adding 20 to dates here, so dates before 2000 won't work
        # TODO: fix pre 2000 dates
        add_these = ['20'+i for i in PTN.findall(os.path.basename(i.data_path))]
        list_of_dates_for_grb_download += add_these
        first_grb = os.path.join(ECMWF_DIR, ECMWF_PRE + add_these[0] + ECMWF_EXT)
        second_grb = os.path.join(ECMWF_DIR, ECMWF_PRE + add_these[1] + ECMWF_EXT)

        # download .grb file if does not exist
        if not (os.path.exists(first_grb) and os.path.exists(second_grb)):
            # download weather files at 12:00 UTC (other options 00:00, 06:00, 18:00)
            pa.ecmwf_download(add_these, '12', 'ECMWF')

        """ using rdr coordinates to remove APS """
        aps1 = pa.PyAPS_rdr(
            os.path.join(ECMWF_DIR, ECMWF_PRE + add_these[0] + ECMWF_EXT),
            SYD_TEST_DEM_UNW, grib='ECMWF', verb=True,
            demfmt='HGT', demtype=np.int16)
        aps2 = pa.PyAPS_rdr(
            os.path.join(ECMWF_DIR, ECMWF_PRE + add_these[1] + ECMWF_EXT),
            SYD_TEST_DEM_UNW, grib='ECMWF', verb=True,
            demfmt='HGT', demtype=np.int16)

        phs1 = np.zeros((aps1.ny, aps1.nx))
        phs2 = np.zeros((aps2.ny, aps2.nx))

        print 'Without Lat Lon files'
        # using random incidence angle
        aps1.getdelay(phs1, inc=23.0)  # MG `to describe how to find incidence
        aps2.getdelay(phs2, inc=23.0)
        aps_delay = phs2-phs1  # delay in meters as we don't provide wavelength
        print aps_delay

        # TODO: lat lon correction when lat and lon files are available
        # aps1.getgeodelay(phs1, inc=23.0, wvl=0.056,
        #   lat=os.path.join(PYAPS_EXAMPLES, 'lat.flt'),
        #   lon=os.path.join(PYAPS_EXAMPLES, 'lon.flt'))
        # aps2.getgeodelay(phs2, inc=23.0, wvl=0.056,
        #   lat=os.path.join(PYAPS_EXAMPLES, 'lat.flt'),
        #   lon=os.path.join(PYAPS_EXAMPLES, 'lon.flt'))
        # LLphs = phs2-phs1

        """ using geo coordinates to remove APS """
        aps1 = pa.PyAPS_geo(
            os.path.join(ECMWF_DIR, ECMWF_PRE + add_these[0] + ECMWF_EXT),
            SYD_TEST_DEM_UNW, grib='ECMWF', verb=True,
            demfmt='HGT', demtype=np.int16)
        aps2 = pa.PyAPS_geo(
            os.path.join(ECMWF_DIR, ECMWF_PRE + add_these[1] + ECMWF_EXT),
            SYD_TEST_DEM_UNW, grib='ECMWF', verb=True,
            demfmt='HGT', demtype=np.int16)

        phs1 = np.zeros((aps1.ny, aps1.nx))
        phs2 = np.zeros((aps2.ny, aps2.nx))

        print 'Without Lat Lon files'
        # using random incidence angle
        aps1.getdelay(phs1, inc=23.0)  # MG `to describe how to find incidence
        aps2.getdelay(phs2, inc=23.0)
        aps_delay = phs2-phs1  # delay in meters as we don't provide wavelength
        print aps_delay
