#  This Python module is part of the PyRate software package.
#
#  Copyright 2020 Geoscience Australia
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#  http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
import os
from pyrate.core.gdal_python import coherence_masking
from pyrate.configuration import Configuration
from osgeo import gdal, osr, gdalconst
import numpy as np

# fix seed so same test dataset is being created each time this script is run
np.random.seed(0)

config_file = "C:/Users/sheec/Desktop/Projects/PyRate/sample_data/input_parameters.conf"
config_file = os.path.abspath(config_file)
config = Configuration(config_file)

for interferogram_file in config.interferogram_files:
    input_file = interferogram_file.converted_path
    output_mask_file = interferogram_file.converted_path.replace(".tif","_mask.tif")
    output_mask_file = interferogram_file.converted_path.replace("input\\interferograms", "input\\coherence")

    gtif = gdal.Open(input_file)
    # gtifMetadata = gtif.GetMetadata()

    band = 1
    gtifBand = gtif.GetRasterBand(band)

    if gtifBand is None:
        print("No band found")

    stats = gtifBand.GetStatistics(True, True)
    transform = gtif.GetGeoTransform()

    raster_info = {
        "min": stats[0],
        "max": stats[1],
        "mean": stats[2],
        "stdDev": stats[3],
        "noDataValue": gtifBand.GetNoDataValue(),
        "scale": gtifBand.GetScale(),
        "unitType": gtifBand.GetUnitType(),
        "xOrigin": transform[0],
        "yOrigin": transform[3],
        "pixelWidth": transform[1],
        "pixelHeight": transform[5],
        "cols":gtif.RasterXSize,
        "rows":gtif.RasterYSize
    }

    syntactic_data_array = np.random.random((raster_info["rows"], raster_info["cols"]))

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(output_mask_file, raster_info["cols"], raster_info["rows"], 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform((raster_info["xOrigin"], raster_info["pixelWidth"], 0, raster_info["yOrigin"], 0, raster_info["pixelHeight"]))

    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(syntactic_data_array)

    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromWkt(gtif.GetProjectionRef())
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()

    # close dataset
    gtif = None

# # masking
# source_dataset_path = "C:/Users/sheec/Desktop/Projects/PyRate/sample_data/input/interferograms/20060828-20061211_utm_unw_1rlks_4cr.tif"
# source_dataset = gdal.Open(source_dataset_path)
# coherence_thresh = 0.16919
# coherence_file_paths = config.coherence_file_paths
# coherence_masking(source_dataset, coherence_file_paths, coherence_thresh)
