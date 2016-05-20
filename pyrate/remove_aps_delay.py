"""
This is a demo of how to remove pyaps delay using PyAPS
Usage: first run pyrate/scripts/run_prepifg.py using your config file
Then use the command:

"""
__author__ = 'Sudipta Basak'
__date_created__ = '4/02/16'
import PyAPS as pa
import numpy as np
import sys
import os
import re
import glob2
# from pyrate.tests.common import sydney_data_setup
from pyrate.tests.common import SYD_TEST_DEM_ROIPAC
from pyrate import config as cf
from pyrate import ifgconstants as ifc
from pyrate import prepifg
from pyrate import gamma


PYRATEPATH = os.environ['PYRATEPATH']
ECMWF_DIR = os.path.join(PYRATEPATH, 'ECMWF')
ECMWF_PRE = 'ERA-Int_'
ECMWF_EXT = '_12.grib'
APS_STATUS = 'REMOVED'
GEOTIFF = 'GEOTIFF'

def remove_aps_delay(ifgs, params):
    list_of_dates_for_grb_download = []

    incidence_angle = None

    for ifg in ifgs:  # demo for only one ifg
        if params[cf.PROCESSOR] == 1:  # gamma
            PTN = re.compile(r'\d{8}')
            date_pair = [i for i in PTN.findall(os.path.basename(ifg.data_path))]
        elif params[cf.PROCESSOR] == 1:  # roipac
            # adding 20 to dates here, so dates before 2000 won't work
            # TODO: fix pre 2000 dates
            PTN = re.compile(r'\d{8}')
            date_pair = ['20' + i for i in
                         PTN.findall(os.path.basename(ifg.data_path))]
        else:
            raise AttributeError('processor needs to be gamma(1) or roipac(0)')

        list_of_dates_for_grb_download += date_pair
        first_grb = os.path.join(ECMWF_DIR,
                                 ECMWF_PRE + date_pair[0] + ECMWF_EXT)
        second_grb = os.path.join(ECMWF_DIR,
                                  ECMWF_PRE + date_pair[1] + ECMWF_EXT)

        # download .grb file if does not exist
        if not (os.path.exists(first_grb) and os.path.exists(second_grb)):
            # download weather files at 12:00 UTC (other options 00:00, 06:00, 18:00)
            pa.ecmwf_download(date_pair, '12', 'ECMWF')

        # rdr_correction(date_pair)

        # TODO: lat lon correction when lat and lon files are available
        # aps1.getgeodelay(phs1, inc=23.0, wvl=0.056,
        #   lat=os.path.join(PYAPS_EXAMPLES, 'lat.flt'),
        #   lon=os.path.join(PYAPS_EXAMPLES, 'lon.flt'))
        # aps2.getgeodelay(phs2, inc=23.0, wvl=0.056,
        #   lat=os.path.join(PYAPS_EXAMPLES, 'lat.flt'),
        #   lon=os.path.join(PYAPS_EXAMPLES, 'lon.flt'))
        # LLphs = phs2-phs1

        # no need to calculate incidence angle for all ifgs, they are the same
        if incidence_angle is None:
            incidence_angle = get_incidence_angle(date_pair, params)

        aps_delay = geo_correction(date_pair, params, incidence_angle)
        ifg.phase_data -= aps_delay  # remove delay
        ifg.dataset.SetMetadataItem(ifc.PYRATE_APS_ERROR, APS_STATUS)
        ifg.write_modified_phase()


def rdr_correction(date_pair):
    """ using rdr coordinates to remove APS """
    raise NotImplementedError("This has not been implemented yet for PyRate")
    # aps1 = pa.PyAPS_rdr(
    #     os.path.join(ECMWF_DIR, ECMWF_PRE + date_pair[0] + ECMWF_EXT),
    #     SYD_TEST_DEM_ROIPAC, grib='ECMWF', verb=True,
    #     demfmt='HGT', demtype=np.int16)
    # aps2 = pa.PyAPS_rdr(
    #     os.path.join(ECMWF_DIR, ECMWF_PRE + date_pair[1] + ECMWF_EXT),
    #     SYD_TEST_DEM_ROIPAC, grib='ECMWF', verb=True,
    #     demfmt='HGT', demtype=np.int16)
    # phs1 = np.zeros((aps1.ny, aps1.nx))
    # phs2 = np.zeros((aps2.ny, aps2.nx))
    # print 'Without Lat Lon files'
    # # using random incidence angle
    # aps1.getdelay(phs1, inc=23.0)  # TODO: MG to describe how to find incidence
    # aps2.getdelay(phs2, inc=23.0)
    # aps_delay = phs2 - phs1  # delay in meters as we don't provide wavelength
    # return aps_delay


def geo_correction(date_pair, params, incidence_angle):

    dem_file = params[cf.DEM_FILE]
    geotif_dem = os.path.join(
        params[cf.OUT_DIR], os.path.basename(dem_file).split('.')[0] + '.tif')

    mlooked_dem = prepifg.mlooked_path(geotif_dem,
                                       looks=params[cf.IFG_LKSX],
                                       crop_out=params[cf.IFG_CROP_OPT])

    # make sure mlooked dem exist
    if not os.path.exists(mlooked_dem):
        raise prepifg.PreprocessError('mlooked dem was not found.'
                                      'Please run prepifg.')

    dem_header = gamma.parse_dem_header(params[cf.DEM_HEADER_FILE])
    lat, lon, nx, ny = return_pyaps_lat_lon(dem_header)


    """ using geo coordinates to remove APS """

    aps1 = pa.PyAPS_geo(
        os.path.join(ECMWF_DIR, ECMWF_PRE + date_pair[0] + ECMWF_EXT),
        mlooked_dem, grib='ECMWF', verb=True,
        demfmt=GEOTIFF, demtype=np.float32, dem_header=(lon, lat, nx, ny))
    aps2 = pa.PyAPS_geo(
        os.path.join(ECMWF_DIR, ECMWF_PRE + date_pair[1] + ECMWF_EXT),
        mlooked_dem, grib='ECMWF', verb=True,
        demfmt=GEOTIFF, demtype=np.float32, dem_header=(lon, lat, nx, ny))
    phs1 = np.zeros((aps1.ny, aps1.nx))
    phs2 = np.zeros((aps2.ny, aps2.nx))
    print 'Without Lat Lon files'
    # using random incidence angle
    if params[cf.APS_METHOD] == 1:
        aps1.getdelay(phs1, inc=incidence_angle)
        aps2.getdelay(phs2, inc=incidence_angle)
    elif params[cf.APS_METHOD] == 2:
        pass
    else:
        raise
    aps_delay = phs2 - phs1  # delay in meters as we don't provide wavelength
    return aps_delay


def get_incidence_angle(date_pair, params):
    # incidence angle exists in unw header files, not in dem
    SLC_DIR = params[cf.SLC_DIR] if params[cf.SLC_DIR] else \
        params[cf.OBS_DIR]
    header_path = glob2.glob(os.path.join(
        SLC_DIR, '**/*%s*slc.par' % date_pair[0]))[0]
    header = gamma.parse_epoch_header(header_path)
    incidence_angle = header[ifc.INCIDENCE_ANGLE]
    return incidence_angle


def return_pyaps_lat_lon(dem_header):
    nx, ny = dem_header[ifc.PYRATE_NCOLS], dem_header[ifc.PYRATE_NROWS]
    lat = np.zeros((2, 1))
    lon = np.zeros((2, 1))
    lat[1] = dem_header[ifc.PYRATE_LAT]
    lon[0] = dem_header[ifc.PYRATE_LONG]
    if lon[0] < 0:
        lon[0] += 360.0
    dx = np.float(dem_header['X_STEP'])
    dy = np.float(dem_header['Y_STEP'])
    lat[0] = lat[1] + dy * ny
    lon[1] = lon[0] + dx * nx
    return lat, lon, nx, ny

