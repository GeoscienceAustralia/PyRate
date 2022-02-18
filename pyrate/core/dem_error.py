#   This Python module is part of the PyRate software package.
#
#   Copyright 2022 Geoscience Australia
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
"""
This Python module implements the calculation of correction for residual topographic effects (a.k.a. as DEM errors).
"""
# pylint: disable=invalid-name, too-many-locals, too-many-arguments
import os
import numpy as np
from typing import Tuple, Optional
from os.path import join
from pathlib import Path

import pyrate.constants as C
from pyrate.core import geometry, shared, mpiops, ifgconstants as ifc
from pyrate.core.logger import pyratelogger as log
from pyrate.core.shared import Ifg, Geometry, DEM, Tile, tiles_split
from pyrate.core.timeseries import TimeSeriesError
from pyrate.configuration import MultiplePaths, Configuration
from pyrate.merge import assemble_tiles


def dem_error_calc_wrapper(params: dict) -> None:
    """
    MPI wrapper for DEM error correction
    :param params: Dictionary of PyRate configuration parameters.
    """
    if not params[C.DEMERROR]:
        log.info("DEM error correction not required")
        return

    # geometry information needed to calculate Bperp for each pixel using first IFG in list
    ifg_paths = [ifg_path.tmp_sampled_path for ifg_path in params[C.INTERFEROGRAM_FILES]]

    # check if DEM error correction is already available
    if mpiops.run_once(__check_and_apply_demerrors_found_on_disc, ifg_paths, params):
        log.warning("Reusing DEM error correction from previous run!!!")
    else:
        log.info("Calculating DEM error correction")
        # read and open the first IFG in list
        ifg0_path = ifg_paths[0]
        ifg0 = Ifg(ifg0_path)
        ifg0.open(readonly=True)

        # not currently implemented for ROIPAC data which breaks some tests
        # if statement can be deleted once ROIPAC is deprecated from PyRate
        if ifg0.meta_data[ifc.PYRATE_INSAR_PROCESSOR] == 'ROIPAC':
            log.warning("Geometry calculations are not implemented for ROI_PAC")
            return

        if params[C.BASE_FILE_LIST] is None:
            log.warning("No baseline files supplied: DEM error has not been computed")
            return

        # todo: subtract other corrections (e.g. orbital) from displacement phase before estimating the DEM error

        # the following code is a quick way to do the bperp calculation, but is not identical to the GAMMA output
        # where the near range of the first SLC is used for each pair.
        # calculate look angle for interferograms (using the Near Range of the first SLC)
        # look_angle = geometry.calc_local_geometry(ifg0, None, rg, lon, lat, params)
        # bperp = geometry.calc_local_baseline(ifg0, az, look_angle)

        # split full arrays in to tiles for parallel processing
        tiles_split(_process_dem_error_per_tile, params)

        preread_ifgs = params[C.PREREAD_IFGS]
        # write dem error and correction values to file
        mpiops.run_once(_write_dem_errors, ifg_paths, params, preread_ifgs)
        shared.save_numpy_phase(ifg_paths, params)

        log.debug('Finished DEM error correction step')


def _process_dem_error_per_tile(tile: Tile, params: dict) -> None:
    """
    Convenience function for processing DEM error in tiles
    :param tile: pyrate.core.shared.Tile Class object.
    :param params: Dictionary of PyRate configuration parameters.
    """
    ifg_paths = [ifg_path.tmp_sampled_path for ifg_path in params[C.INTERFEROGRAM_FILES]]
    ifg0_path = ifg_paths[0]
    ifg0 = Ifg(ifg0_path)
    ifg0.open(readonly=True)
    # read lon and lat values of multi-looked ifg (first ifg only)
    lon, lat = geometry.get_lonlat_coords(ifg0)
    # read azimuth and range coords and DEM from tif files generated in prepifg
    geom_files = Configuration.geometry_files(params)
    rdc_az_file = geom_files['rdc_azimuth']
    geom_az = Geometry(rdc_az_file)
    rdc_rg_file = geom_files['rdc_range']
    geom_rg = Geometry(rdc_rg_file)
    dem_file = params[C.DEM_FILE_PATH].sampled_path
    dem = DEM(dem_file)
    preread_ifgs = params[C.PREREAD_IFGS]
    threshold = params[C.DE_PTHR]
    ifg_parts = [shared.IfgPart(p, tile, preread_ifgs, params) for p in ifg_paths]
    lon_parts = lon(tile)
    lat_parts = lat(tile)
    az_parts = geom_az(tile)
    rg_parts = geom_rg(tile)
    dem_parts = dem(tile)
    log.debug(f"Calculating per-pixel baseline for tile {tile.index} during DEM error correction")
    bperp, look_angle, range_dist = _calculate_bperp_wrapper(ifg_paths, az_parts, rg_parts,
                                                             lat_parts, lon_parts, dem_parts)
    log.debug(f"Calculating DEM error for tile {tile.index} during DEM error correction")

    # mst_tile = np.load(Configuration.mst_path(params, tile.index))
    # calculate the DEM error estimate and the correction values for each IFG
    # current implementation uses the look angle and range distance matrix of the primary SLC in the last IFG
    # todo: check the impact of using the same information from another SLC
    dem_error, dem_error_correction, _ = calc_dem_errors(ifg_parts, bperp, look_angle, range_dist, threshold)
    # dem_error contains the estimated DEM error for each pixel (i.e. the topographic change relative to the DEM)
    # size [row, col]
    # dem_error_correction contains the correction value for each interferogram
    # size [num_ifg, row, col]
    # save tiled data in tmpdir
    np.save(file=os.path.join(params[C.TMPDIR], 'dem_error_{}.npy'.format(tile.index)), arr=dem_error)
    # swap the axes of 3D array to fit the style used in function assemble_tiles
    tmp_array = np.moveaxis(dem_error_correction, 0, -1)
    # new dimension is [row, col, num_ifg]
    # save tiled data into tmpdir
    np.save(file=os.path.join(params[C.TMPDIR], 'dem_error_correction_{}.npy'.format(tile.index)), arr=tmp_array)

    # Calculate and save the average perpendicular baseline for the tile
    bperp_avg =  np.nanmean(bperp, axis=(1, 2), dtype=np.float64)
    np.save(file=os.path.join(params[C.TMPDIR], 'bperp_avg_{}.npy'.format(tile.index)), arr=bperp_avg)


def _calculate_bperp_wrapper(ifg_paths: list, az_parts: np.ndarray, rg_parts: np.ndarray,
                             lat_parts: np.ndarray, lon_parts: np.ndarray, dem_parts: np.ndarray,
                             ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Wrapper function to calculate the perpendicular baseline for each pixel in each interferogram.
    :param ifg_paths: list of pyrate.core.shared.Ifg Class objects.
    :param az_parts: azimuth coordinate (i.e. line) for each pixel.
    :param rg_parts: range coordinate (i.e. column) for each pixel.
    :param lat_parts: latitude for each pixel.
    :param lon_parts: longitude for each pixel.
    :param dem_parts: DEM height for each pixel.
    :return: bperp: perpendicular baseline for each pixel and interferogram.
    :return: look_angle: look angle for each pixel.
    :return: range_dist: range distance measurement for each pixel.
    """
    nifgs = len(ifg_paths)
    bperp = np.empty((nifgs, lon_parts.shape[0], lon_parts.shape[1])) * np.nan
    # calculate per-pixel perpendicular baseline for each IFG
    # loop could be avoided by approximating the look angle for the first Ifg
    for ifg_num, ifg_path in enumerate(ifg_paths):
        ifg = Ifg(ifg_path)
        ifg.open(readonly=True)
        # calculate look angle for interferograms (using the Near Range of the primary SLC)
        look_angle, _, _, range_dist = geometry.calc_pixel_geometry(ifg, rg_parts, lon_parts, lat_parts, dem_parts)
        bperp[ifg_num, :, :] = geometry.calc_local_baseline(ifg, az_parts, look_angle)
    return bperp, look_angle, range_dist


def calc_dem_errors(ifgs: list, bperp: np.ndarray, look_angle: np.ndarray, range_dist: np.ndarray,
                    threshold: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Function to calculate the DEM error for each pixel using least-squares adjustment of phase data and
    perpendicular baseline. The least-squares adjustment co-estimates the velocities.
    - *nrows* is the number of rows in the ifgs, and
    - *ncols* is the  number of columns in the ifgs.

    :param ifgs: list of interferogram class objects.
    :param bperp: Per-pixel perpendicular baseline for each interferogram
    :param look_angle: Per-pixel look angle
    :param range_dist: Per-pixel range distance measurement
    :param threshold: minimum number of redundant phase values at pixel (config parameter de_pthr)
    :return: dem_error: estimated per-pixel dem error (nrows x ncols)
    :return: dem_error_correction: DEM error correction for each pixel and interferogram (nifgs x nrows x ncols)
    :return: vel: velocity estimate for each pixel (nrows x ncols)
    """
    ifg_data, mst, ncols, nrows, ifg_time_span = _per_tile_setup(ifgs)
    if threshold < 4:
        msg = f"pixel threshold too low (i.e. <4) resulting in non-redundant DEM error estimation"
        raise DEMError(msg)
    # pixel-by-pixel calculation
    # preallocate empty arrays for results
    dem_error = np.empty([nrows, ncols], dtype=np.float32) * np.nan
    vel = np.empty([nrows, ncols], dtype=np.float32) * np.nan
    # nested loops to loop over the 2 image dimensions
    for row in range(nrows):
        for col in range(ncols):
            # calc DEM error for each pixel with valid Bperp and ifg phase data
            # check pixel for non-redundant ifgs
            sel = np.nonzero(mst[:, row, col])[0]  # trues in mst are chosen
            if len(sel) >= threshold:  # given threshold for number of valid pixels in time series
                # phase observations (in mm)
                y = ifg_data[sel, row, col]
                bperp_pix = bperp[sel, row, col]
                # If NaN: skip pixel
                if np.isnan(y).any() or np.isnan(bperp_pix).any():
                    continue
                # using the actual geometry of a particular IFG would be possible but is likely not signif. different
                # geom = bperp_pix / (range_dist[row, col] * np.sin(look_angle[row, col]))
                time_span = ifg_time_span[sel]
                # new covariance matrix using actual number of observations
                m = len(sel)
                # Design matrix of least-squares system, velocities are co-estimated as done in StaMPS or MintPy
                A = np.column_stack((np.ones(m), bperp_pix, time_span))
                # solve ordinary least-squares system using numpy.linalg.lstsq function
                xhat, res, rnk, s = np.linalg.lstsq(A, y, rcond=None)
                # dem error estimate for the pixel is the second parameter (cf. design matrix)
                dem_error[row][col] = xhat[1]
                # velocity estimate for the pixel is the third parameter (cf. design matrix)
                vel[row][col] = xhat[2]

    # calculate correction value for each IFG by multiplying the least-squares estimate with the Bperp value
    dem_error_correction = np.multiply(dem_error, bperp)
    # calculate metric difference to DEM by multiplying the estimate with the per-pixel geometry
    # (i.e. range distance and look angle, see Eq. (2.4.12) in Hanssen (2001))
    # also scale by -0.001 since the phase observations are in mm with positive values away from the sensor
    dem_error = np.multiply(dem_error, np.multiply(range_dist, np.sin(look_angle))) * (-0.001)

    return dem_error, dem_error_correction, vel


def _per_tile_setup(ifgs: list) -> Tuple[np.ndarray, np.ndarray, int, int, np.ndarray]:
    """
    Convenience function for setting up DEM error computation parameters
    :param ifgs: list of interferogram class objects.
    :return: ifg_data: phase observations for each pixel and interferogram
    :return: mst: an array of booleans representing valid ifg connections (i.e. the minimum spanning tree)
    :return: ncols: number of columns
    :return: nrows: number of rows
    :return: ifg_time_span: date difference for each interferogram
    """
    if len(ifgs) < 1:
        msg = 'Time series requires 2+ interferograms'
        raise TimeSeriesError(msg)

    nrows = ifgs[0].nrows
    ncols = ifgs[0].ncols
    nifgs = len(ifgs)
    ifg_data = np.zeros((nifgs, nrows, ncols), dtype=np.float32)
    for ifg_num in range(nifgs):
        ifg_data[ifg_num] = ifgs[ifg_num].phase_data # phase data is read from numpy files
    mst = ~np.isnan(ifg_data)
    ifg_time_span = np.zeros((nifgs))
    for ifg_num in range(nifgs):
        ifg_time_span[ifg_num] = ifgs[ifg_num].time_span

    return ifg_data, mst, ncols, nrows, ifg_time_span


def _write_dem_errors(ifg_paths: list, params: dict, preread_ifgs: dict) -> None:
    """
    Convenience function for writing DEM error (one file) and DEM error correction for each IFG to disc
    :param ifg_paths: List of interferogram class objects.
    :param params: Dictionary of PyRate configuration parameters.
    :param preread_ifgs: Dictionary of interferogram metadata.
    """
    tiles = params[C.TILES]

    # re-assemble tiles and save into dem_error dir
    shape = preread_ifgs[ifg_paths[0]].shape

    # save dem error as geotiff file in out directory
    gt, md, wkt = shared.get_geotiff_header_info(ifg_paths[0])
    md[ifc.DATA_TYPE] = ifc.DEM_ERROR
    dem_error = assemble_tiles(shape, params[C.TMPDIR], tiles, out_type='dem_error')
    dem_error_file = os.path.join(params[C.DEM_ERROR_DIR], 'dem_error.tif')
    shared.remove_file_if_exists(dem_error_file)
    shared.write_output_geotiff(md, gt, wkt, dem_error, dem_error_file, np.nan)

    # read the average bperp vals for each ifg and each tile
    bperp = np.empty(shape=(len(tiles), len(ifg_paths)), dtype=np.float64)
    for t in tiles:
        bperp_file = Path(join(params[C.TMPDIR], 'bperp_avg_' + str(t.index) + '.npy'))
        arr = np.load(file=bperp_file)
        bperp[t.index, :] = arr

    # loop over all ifgs
    idx = 0
    for ifg_path in ifg_paths:
        ifg = Ifg(ifg_path)
        ifg.open()
        shared.nan_and_mm_convert(ifg, params) # ensure we have phase data in mm
        # read dem error correction file from tmpdir
        dem_error_correction_ifg = assemble_tiles(shape, params[C.TMPDIR], tiles, out_type='dem_error_correction',
                                                  index=idx)
        # calculate average bperp value across all tiles for the ifg
        bperp_val = np.nanmean(bperp[:, idx])
        dem_error_correction_on_disc = MultiplePaths.dem_error_path(ifg.data_path, params)
        np.save(file=dem_error_correction_on_disc, arr=dem_error_correction_ifg)
        idx += 1

        # subtract DEM error from the ifg
        ifg.phase_data -= dem_error_correction_ifg
        _save_dem_error_corrected_phase(ifg, bperp_val)


def __check_and_apply_demerrors_found_on_disc(ifg_paths: list, params: dict) -> bool:
    """
    Convenience function to check if DEM error correction files have already been produced in a previous run
    :param ifg_paths: List of interferogram class objects.
    :param params: Dictionary of PyRate configuration parameters.
    :return: bool value: True if dem error files found on disc, otherwise False.
    """
    saved_dem_err_paths = [MultiplePaths.dem_error_path(ifg_path, params) for ifg_path in ifg_paths]
    for d, i in zip(saved_dem_err_paths, ifg_paths):
        if d.exists():
            dem_corr = np.load(d)
            if isinstance(i, str):
                # are paths
                ifg = Ifg(i)
                ifg.open()
            else:
                ifg = i
            ifg.phase_data -= dem_corr
            # set geotiff meta tag and save phase to file
            # TODO: calculate avg bperp and add to metadata even for reused DEM error correction
            _save_dem_error_corrected_phase(ifg)

    return all(d.exists() for d in saved_dem_err_paths)


def _save_dem_error_corrected_phase(ifg: Ifg, bperp: Optional[np.float64] = None) -> None:
    """
    Convenience function to update metadata and save latest phase after DEM error correction
    :param ifg: pyrate.core.shared.Ifg Class object
    :param bperp: perpendicular baseline value for adding to geotiff metadata.
    """
    # update geotiff tags after DEM error correction
    ifg.dataset.SetMetadataItem(ifc.PYRATE_DEM_ERROR, ifc.DEM_ERROR_REMOVED)
    if bperp is not None:
        ifg.dataset.SetMetadataItem(ifc.PYRATE_BPERP, str(bperp))
    ifg.write_modified_phase()
    ifg.close()


class DEMError(Exception):
    """
    Generic exception for DEM errors.
    """
