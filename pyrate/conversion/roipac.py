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


def convert_dem_interferogram(parameters):

    config, interferogram_path = parameters
    destination_path = pathlib.Path(interferogram_path.parent.as_posix() + '/' + interferogram_path.stem + '_dem.tif')

    # read dem headers
    dem_header = read_raw_header_file(config.dem_header_path)

    # get extent for the dataset
    x_size = format_header_value(dem_header["root"]["width"], r'[\+\-]?\d+', int)
    y_size = format_header_value(dem_header["root"]["nlines"], r'[\+\-]?\d+', int)

    longitude = format_header_value(dem_header["root"]["corner_lon"], r'[\+\-]?\d+[.]\d+', float)
    x_step = format_header_value(dem_header["root"]["post_lon"], r'[\+\-]?\d+[.]\d+[e][\+\-]\d+', float)

    latitude = format_header_value(dem_header["root"]["corner_lat"], r'[\+\-]?\d+[.]\d+', float)
    y_step = format_header_value(dem_header["root"]["post_lat"], r'[\+\-]?\d+[.]\d+[e][\-]\d+', float)

    # create an empty dataset
    driver = gdal.GetDriverByName("GTiff")
    output_dataset = driver.Create(str(destination_path), x_size, y_size, DATA_BANDS, gdal.GDT_Float32, options=['compress=packbits'])

    # add geo-spatial info
    geo_transform = [longitude, x_step, 0, latitude, 0, y_step]
    output_dataset.SetGeoTransform(geo_transform)
    srs = osr.SpatialReference()
    wkt_projection = srs.ExportToWkt()
    output_dataset.SetProjection(wkt_projection)
    band = output_dataset.GetRasterBand(1)
    band.SetNoDataValue(NO_DATA_VALUE)

    # create output dataset
    fmtstr = '<' + ('h' * x_size)  # roipac DEM is little endian signed int16
    bytes_per_col = 2
    row_bytes = x_size * bytes_per_col

    with open(interferogram_path, 'rb') as f:
        # Read the input array byte by byte and write it to the new dataset
        for y in range(y_size):
            data = struct.unpack(fmtstr, f.read(row_bytes))
            # write data to geo-tiff
            band.WriteArray(np.array(data).reshape(1, x_size), yoff=y)

    # manual close dataset
    band.FlushCache()  # Write to disk
    output_dataset = None
    del output_dataset
    print("Finish processing interferogram", interferogram_path.stem)
    return destination_path


def convert_roipac_interferogram(parameters):

    config, interferogram_path = parameters
    print("Start  processing interferogram", interferogram_path.stem)
    # for the given interferogram find the master and slave header files
    header_master_path, header_slave_path = sort_headers(config.header_paths, interferogram_path)
    destination_path = pathlib.Path(interferogram_path.parent.as_posix() + '/' + interferogram_path.stem + '_unw.tif')

    # read headers
    master_header = read_raw_header_file(header_master_path)
    slave_header = read_raw_header_file(header_slave_path)

    # read dem headers
    dem_header = read_raw_header_file(config.dem_header_path)

    # get extent for the dataset
    x_size = format_header_value(dem_header["root"]["width"], r'[\+\-]?\d+', int)
    y_size = format_header_value(dem_header["root"]["nlines"], r'[\+\-]?\d+', int)

    longitude = format_header_value(dem_header["root"]["corner_lon"], r'[\+\-]?\d+[.]\d+', float)
    x_step = format_header_value(dem_header["root"]["post_lon"], r'[\+\-]?\d+[.]\d+[e][\+\-]\d+', float)

    latitude = format_header_value(dem_header["root"]["corner_lat"], r'[\+\-]?\d+[.]\d+', float)
    y_step = format_header_value(dem_header["root"]["post_lat"], r'[\+\-]?\d+[.]\d+[e][\-]\d+', float)

    # create an empty dataset
    driver = gdal.GetDriverByName("GTiff")
    output_dataset = driver.Create(str(destination_path), x_size, y_size, DATA_BANDS, gdal.GDT_Float32, options=['compress=packbits'])

    # add geo-spatial info
    geo_transform = [longitude, x_step, 0, latitude, 0, y_step]
    output_dataset.SetGeoTransform(geo_transform)
    srs = osr.SpatialReference()
    wkt_projection = srs.ExportToWkt()
    output_dataset.SetProjection(wkt_projection)
    band = output_dataset.GetRasterBand(1)
    band.SetNoDataValue(NO_DATA_VALUE)

    # create metadata
    master_frequency = format_header_value(master_header["root"]["radar_frequency"], r'[\+\-]?\d+[.]\d+[e][\+\-]\d+', float, 0.00001)
    slave_frequency = format_header_value(master_header["root"]["radar_frequency"], r'[\+\-]?\d+[.]\d+[e][\+\-]\d+', float, 0.00001)

    mean_frequency = statistics.mean([master_frequency, slave_frequency])
    round_spaces = len(str(mean_frequency).split(".")[1]) - 1
    mean_frequency = round(mean_frequency, round_spaces)

    wavelength_metres = SPEED_OF_LIGHT_METRES_PER_SECOND / mean_frequency

    master_date = datetime.strptime(master_header["root"]["date"], "%Y %d %M")
    slave_date = datetime.strptime(slave_header["root"]["date"], "%Y %d %M")
    time_span_year = (master_date.day - slave_date.day) / DAYS_PER_YEAR

    master_time = str(timedelta(seconds=float(re.match(r'\d+[\.]?\d+', master_header["root"]["center_time"])[0]))).split('.')[0]
    slave_time = str(timedelta(seconds=float(re.match(r'\d+[\.]?\d+', slave_header["root"]["center_time"])[0]))).split('.')[0]

    incidence_angle_master = float(re.match(r'\d+[.]?\d+', master_header["root"]["incidence_angle"])[0])
    incidence_angle_slave = float(re.match(r'\d+[.]?\d+', slave_header["root"]["incidence_angle"])[0])
    incidence_angle_mean = statistics.mean([incidence_angle_master, incidence_angle_slave])
    round_spaces = len(str(incidence_angle_mean).split(".")[1]) - 1
    incidence_angle_mean = round(incidence_angle_mean, round_spaces)

    print(master_time, slave_time)

    metadata = {
        'WAVELENGTH_METRES': wavelength_metres,
        'TIME_SPAN_YEAR': time_span_year,
        'INSAR_PROCESSOR': 'GAMMA',
        'MASTER_DATE': master_date.strftime("%Y-%d-%M"),
        'SLAVE_DATE': slave_date.strftime("%Y-%d-%M"),
        'DATA_UNITS': 'RADIANS',
        'DATA_TYPE': 'ORIGINAL_IFG',
        'MASTER_TIME': master_time,
        'SLAVE_TIME': slave_time,
        'INCIDENCE_DEGREES': incidence_angle_mean
    }

    if metadata is not None:
        for k, v in metadata.items():
            output_dataset.SetMetadataItem(k, str(v))

    # create output dataset
    fmtstr = '<' + ('f' * x_size)  # roipac ifgs are little endian float32s
    bytes_per_col = 4
    row_bytes = x_size * bytes_per_col

    with open(interferogram_path, 'rb') as f:
        # Read the input array byte by byte and write it to the new dataset
        for y in range(y_size):
            f.seek(row_bytes, 1)  # skip interleaved band 1
            data = struct.unpack(fmtstr, f.read(row_bytes))
            # write data to geo-tiff
            band.WriteArray(np.array(data).reshape(1, x_size), yoff=y)

    # manual close dataset
    band.FlushCache()  # Write to disk
    output_dataset = None
    del output_dataset
    print("Finish processing interferogram", interferogram_path.stem)
    return destination_path


if __name__ == "__main__":

    start_time = time.time()
    gdal.SetCacheMax(GDAL_CACHE_MAX)

    # Input parameters
    config_file_path = "C:/Users/sheec/Desktop/Projects/PyRate/input_parameters_test.conf"
    config = Configuration(config_file_path)

    # Init multiprocessing.Pool()
    pool = mp.Pool(mp.cpu_count())

    # Running pools
    destination_paths = pool.map(convert_roipac_interferogram, [(config, interferogram_path) for interferogram_path in config.interferogram_paths])

    # Closing pools
    pool.close()

    destination_paths_string = ""
    for destination_path in destination_paths:
        destination_paths_string += str(destination_path) + "\n"


    parameters = (config, config.dem_header_path)
    destination_path = convert_dem_interferogram(parameters)
    destination_paths_string += str(destination_path) + "\n"

    config.output_tiff_list.write_text(destination_paths_string)
    print("--- %s seconds ---" % (time.time() - start_time))