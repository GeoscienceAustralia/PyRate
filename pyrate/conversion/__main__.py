from osgeo import osr, gdal
import numpy as np
import struct
from datetime import datetime, timedelta
import statistics
from utilities import *
from constants import *
from pyrate.configuration import Configuration
import time
import multiprocessing as mp
import pathlib
import gamma
import roipac

if __name__ == "__main__":

    start_time = time.time()
    gdal.SetCacheMax(GDAL_CACHE_MAX)

    # Input parameters
    config_file_path = "C:/Users/sheec/Desktop/Projects/PyRate/input_parameters.conf"
    config = Configuration(config_file_path)

    if config.processor:

        # Init multiprocessing.Pool()
        pool = mp.Pool(mp.cpu_count())

        # Running pools
        destination_paths = pool.map(gamma.convert_gamma_interferogram, [(config, interferogram_path) for interferogram_path in config.interferogram_paths])

        # Closing pools
        pool.close()

        destination_paths_string = ""
        for destination_path in destination_paths:
            destination_paths_string += str(destination_path) + "\n"

        parameters = (config, config.dem_path)
        destination_path = gamma.convert_dem_interferogram(parameters)
        destination_paths_string += str(destination_path) + "\n"

        config.output_tiff_list.write_text(destination_paths_string)

    else:
        # Init multiprocessing.Pool()
        pool = mp.Pool(mp.cpu_count())

        # Running pools
        destination_paths = pool.map(roipac.convert_roipac_interferogram, [(config, interferogram_path) for interferogram_path in config.interferogram_paths])

        # Closing pools
        pool.close()

        destination_paths_string = ""
        for destination_path in destination_paths:
            destination_paths_string += str(destination_path) + "\n"

        parameters = (config, config.dem_header_path)
        destination_path = roipac.convert_dem_interferogram(parameters)
        destination_paths_string += str(destination_path) + "\n"

        config.output_tiff_list.write_text(destination_paths_string)

    print("--- %s seconds ---" % (time.time() - start_time))
