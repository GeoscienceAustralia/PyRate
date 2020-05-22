import os
import stat
import tempfile
import numpy as np
from osgeo import osr
from osgeo import gdal
from pathlib import Path
from copy import copy

import pyrate.core.shared
from pyrate.core.shared import Ifg
from pyrate.core import gdal_python
from pyrate.core import config as cf
from pyrate.core import prepifg_helper
from pyrate.core import ifgconstants as ifc
from pyrate.configuration import MultiplePaths
from pyrate import conv2tif


def test_small_data_coherence(gamma_params):

    gamma_params[cf.COH_MASK] = 1

    ifg_multilist = copy(gamma_params[cf.INTERFEROGRAM_FILES])
    conv2tif.main(gamma_params)

    for i in ifg_multilist:
        p = Path(i.converted_path)
        p.chmod(0o664)  # assign write permission as conv2tif output is readonly
        ifg = pyrate.core.shared.dem_or_ifg(data_path=p.as_posix())
        if not isinstance(ifg, Ifg):
            continue
        ifg.open()
        # now do coherence masking and compare
        ifg = pyrate.core.shared.dem_or_ifg(data_path=p.as_posix())
        ifg.open()
        converted_coh_file_path = cf.coherence_paths_for(p, gamma_params, tif=True)
        gdal_python.coherence_masking(ifg.dataset,
                                      coherence_file_path=converted_coh_file_path,
                                      coherence_thresh=gamma_params[cf.COH_THRESH]
                                      )
        nans = np.isnan(ifg.phase_data)
        coherence_path = cf.coherence_paths_for(p, gamma_params, tif=True)
        cifg = Ifg(coherence_path)
        cifg.open()
        cifg_below_thrhold = cifg.phase_data < gamma_params[cf.COH_THRESH]
        np.testing.assert_array_equal(nans, cifg_below_thrhold)


def test_coherence_files_not_converted():
    # define constants
    NO_DATA_VALUE = 0
    driver = gdal.GetDriverByName('GTiff')

    # create a sample gdal dataset

    # sample gdal dataset
    sample_gdal_filename = "dataset_01122000.tif"
    options = ['PROFILE=GeoTIFF']
    sample_gdal_dataset = driver.Create(sample_gdal_filename, 5, 5, 1, gdal.GDT_Float32, options=options)
    srs = osr.SpatialReference()
    wkt_projection = srs.ExportToWkt()
    sample_gdal_dataset.SetProjection(wkt_projection)

    sample_gdal_band = sample_gdal_dataset.GetRasterBand(1)
    sample_gdal_band.SetNoDataValue(NO_DATA_VALUE)
    sample_gdal_band.WriteArray(np.arange(25).reshape(5, 5))
    sample_gdal_dataset.SetMetadataItem(ifc.MASTER_DATE, '2019-10-20')
    sample_gdal_dataset.SetMetadataItem(ifc.SLAVE_DATE, '2019-11-01')
    sample_gdal_dataset.SetMetadataItem(ifc.PYRATE_WAVELENGTH_METRES, '10.05656')
    sample_gdal_dataset.FlushCache()
    sample_gdal_dataset = None
    ifg = Ifg(sample_gdal_filename)
    ifg.open()

    # create a coherence mask dataset
    tmpdir = tempfile.mkdtemp()
    out_dir = Path(tmpdir) # we won't be creating any output coherence mask files as there are already GeoTIFFs
    coherence_mask_filename = MultiplePaths(out_dir, Path("mask_dataset_01122000.tif").as_posix())
    coherence_mask_dataset = driver.Create(coherence_mask_filename.converted_path, 5, 5, 1, gdal.GDT_Float32)
    srs = osr.SpatialReference()
    wkt_projection = srs.ExportToWkt()
    coherence_mask_dataset.SetProjection(wkt_projection)
    coherence_mask_band = coherence_mask_dataset.GetRasterBand(1)
    coherence_mask_band.SetNoDataValue(NO_DATA_VALUE)
    arr = np.arange(0, 75, 3).reshape(5, 5) / 100.0
    arr[3, 4] = 0.25  # insert some random lower than threshold number
    arr[4, 2] = 0.20  # insert some random lower than threshold number

    coherence_mask_band.WriteArray(arr)
    # del the tmp handler datasets created
    del coherence_mask_dataset
    # create an artificial masked dataset
    expected_result_array = np.nan_to_num(
        np.array(
            [
                [np.nan, np.nan, np.nan, np.nan, np.nan],
                [np.nan, np.nan, np.nan, np.nan, np.nan],
                [10.0, 11.0, 12.0, 13.0, 14.0],
                [15.0, 16.0, 17.0, 18.0, np.nan],
                [20.0, 21.0, np.nan, 23.0, 24.0],
            ]
        )
    )

    # use the gdal_python.coherence_masking to find the actual mask dataset
    coherence_thresh = 0.3

    gdal_python.coherence_masking(ifg.dataset, coherence_mask_filename.converted_path, coherence_thresh)

    sample_gdal_array = np.nan_to_num(ifg.phase_data)

    # compare the artificial masked and actual masked datasets
    np.testing.assert_array_equal(sample_gdal_array, expected_result_array)

    # del the tmp datasets created
    os.remove(coherence_mask_filename.converted_path)

    ifg.close()
    os.remove(sample_gdal_filename)
