import os
import unittest
import pytest
import glob
import copy

import pyrate.core.config as cf
from pyrate import converttogeotif, prepifg
from pyrate.core import gdal_python

import os
import gdal
import numpy as np
import osr

NO_DATA_VALUE = 0
driver = gdal.GetDriverByName( 'GTiff' )

# sample gdal dataset
sample_gdal_filename = 'sample_gdal_dataset.tif'
sample_gdal_dataset = driver.Create(sample_gdal_filename, 5, 5, 1, gdal.GDT_Float32)
srs = osr.SpatialReference()
wkt_projection = srs.ExportToWkt()
sample_gdal_dataset.SetProjection(wkt_projection)

sample_gdal_band = sample_gdal_dataset.GetRasterBand(1)
sample_gdal_band.SetNoDataValue(NO_DATA_VALUE)
sample_gdal_band.WriteArray(np.arange(25).reshape(5, 5))


print("sample_gdal_band: \n", sample_gdal_band.ReadAsArray())



# coherence mask dataset
coherence_mask_filename = 'coherence_mask_dataset.tif'
coherence_mask_dataset = driver.Create(coherence_mask_filename, 5, 5, 1, gdal.GDT_Float32)
srs = osr.SpatialReference()
wkt_projection = srs.ExportToWkt()
sample_gdal_dataset.SetProjection(wkt_projection)

coherence_mask_band = coherence_mask_dataset.GetRasterBand(1)
coherence_mask_band.SetNoDataValue(NO_DATA_VALUE)
coherence_mask_band.WriteArray(np.arange(0,75,3).reshape(5, 5)/100.0)

print("coherence_mask_band: \n", coherence_mask_band.ReadAsArray())

threshold = 0.3
print("threshold: ", threshold)

gdal_python.coherence_masking(sample_gdal_dataset, coherence_mask_dataset, threshold)


print("Result: \n", sample_gdal_dataset.GetRasterBand(1).ReadAsArray())

del coherence_mask_dataset
os.remove(coherence_mask_filename)

del sample_gdal_dataset
os.remove(sample_gdal_filename)