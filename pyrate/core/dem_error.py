#   This Python module is part of the PyRate software package.
#
#   Copyright 2020 Geoscience Australia
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
from os.path import join
import pickle as cp
from numpy.linalg import inv, LinAlgError, lstsq

from pyrate.core import geometry, shared, mpiops, config as cf, ifgconstants as ifc
from pyrate.core.logger import pyratelogger as log
from pyrate.core.shared import Ifg, Geometry
from pyrate.core.timeseries import TimeSeriesError
from pyrate.configuration import Configuration, MultiplePaths
from pyrate.merge import assemble_tiles


def dem_error_calc_wrapper(params: dict) -> None:
    """
    MPI wrapper for DEM error correction
    """
    if not params[cf.DEMERROR]:
        log.info("DEM error correction not required")
        return

    # geometry information needed to calculate Bperp for each pixel using first IFG in list
    ifg_paths = [ifg_path.tmp_sampled_path for ifg_path in params[cf.INTERFEROGRAM_FILES]]

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

        if params[cf.BASE_FILE_LIST] is None:
            log.warning("No baseline files supplied: DEM error has not been computed")
            return

        # read azimuth and range coords from tif files generated in prepifg
        rdc_az_file = join(params[cf.OUT_DIR], 'rdc_azimuth.tif')
        geom_az = Geometry(rdc_az_file)
        geom_az.open(readonly=True)
        az = geom_az.geometry_data
        rdc_rg_file = join(params[cf.OUT_DIR], 'rdc_range.tif')
        geom_rg = Geometry(rdc_rg_file)
        geom_rg.open(readonly=True)
        rg = geom_rg.geometry_data

        log.info("Calculating per-pixel baseline")

        # split into tiles to calculate DEM error correction
        tiles = params[cf.TILES]
        preread_ifgs = params[cf.PREREAD_IFGS]
        # todo: subtract other corrections (e.g. orbital) from displacement phase before estimating the DEM error
        threshold = params[cf.DE_PTHR]

        # read lon and lat values of multi-looked ifg (first ifg only)
        lon, lat = geometry.get_lonlat_coords(ifg0)

        # the following code is a quick way to do the bperp calculation, but is not identical to the GAMMA output
        # where the near range of the first SLC is used for each pair.
        # calculate look angle for interferograms (using the Near Range of the first SLC)
        # look_angle = geometry.calc_local_geometry(ifg0, None, rg, lon, lat, params)
        # bperp = geometry.calc_local_baseline(ifg0, az, look_angle)

        # process in tiles; cut arrays to size
        process_tiles = mpiops.array_split(tiles)
        for t in process_tiles:
            ifg_parts = [shared.IfgPart(p, t, preread_ifgs, params) for p in ifg_paths]
            lon_parts = lon[t.top_left_y:t.bottom_right_y, t.top_left_x:t.bottom_right_x]
            lat_parts = lat[t.top_left_y:t.bottom_right_y, t.top_left_x:t.bottom_right_x]
            az_parts = az[t.top_left_y:t.bottom_right_y, t.top_left_x:t.bottom_right_x]
            rg_parts = rg[t.top_left_y:t.bottom_right_y, t.top_left_x:t.bottom_right_x]

            nifgs = len(ifg_paths)
            bperp = np.empty((nifgs, lon_parts.shape[0], lon_parts.shape[1])) * np.nan
            # calculate per-pixel perpendicular baseline for each IFG
            # loop could be avoided by approximating the look angle for the first Ifg
            for ifg_num, ifg_path in enumerate(ifg_paths):
                ifg = Ifg(ifg_path)
                ifg.open(readonly=True)
                # calculate look angle for interferograms (using the Near Range of the primary SLC)
                look_angle, _, _, range_dist = geometry.calc_pixel_geometry(ifg, rg_parts, lon_parts,
                                                                            lat_parts, params, tile=t)
                bperp[ifg_num, :, :] = geometry.calc_local_baseline(ifg, az_parts, look_angle)

            log.debug('Calculating DEM error for tile {} during DEM error correction'.format(t.index))
            # mst_tile = np.load(Configuration.mst_path(params, t.index))
            # calculate the DEM error estimate and the correction values for each IFG
            # current implementation uses the look angle and range distance matrix of the primary SLC in the last IFG
            # todo: check the impact of using the same information from another SLC
            dem_error, dem_error_correction, _ = calc_dem_errors(ifg_parts, bperp, look_angle, range_dist, threshold)
            # dem_error contains the estimated DEM error for each pixel (i.e. the topographic change relative to the DEM)
            # size [row, col]
            # dem_error_correction contains the correction value for each interferogram
            # size [num_ifg, row, col]

            # save tiled data in tmpdir
            np.save(file=os.path.join(params[cf.TMPDIR], 'dem_error_{}.npy'.format(t.index)), arr=dem_error)
            # swap the axes of 3D array to fit the style used in function assemble_tiles
            tmp_array = np.moveaxis(dem_error_correction, 0, -1)
            # new dimension is [row, col, num_ifg]
            # save tiled data into tmpdir
            np.save(file=os.path.join(params[cf.TMPDIR], 'dem_error_correction_{}.npy'.format(t.index)),
                    arr=tmp_array)

        # wait for all processes to finish
        mpiops.comm.barrier()

        # write dem error and correction values to file
        mpiops.run_once(_write_dem_errors, ifg_paths, params, preread_ifgs, tiles)
        shared.save_numpy_phase(ifg_paths, params)

        log.debug('Finished DEM error correction step')


def calc_dem_errors(ifgs, bperp, look_angle, range_dist, threshold):
    """
    Calculates the per-pixel DEM error using least-squares adjustment of phase data and
    perpendicular baseline. The least-squares adjustment co-estimates the velocities.

        - *nrows* is the number of rows in the ifgs,
        - *ncols* is the  number of columns in the ifgs, and

    :param list ifgs: list of interferogram class objects.
    :param ndarray bperp: Per-pixel perpendicular baseline for each interferogram
    :param ndarray look_angle: Per-pixel look angle
    :param ndarray range_dist: Per-pixel range distance measurement
    :param int threshold: minimum number of redundant phase values at pixel (config parameter de_pthr)

    :return: ndarray dem_error: estimated per-pixel dem error (nrows x ncols)
    :return: ndarray dem_error_correction: DEM error correction for each pixel and interferogram (nifgs x nrows x ncols)
    :return: ndarray vel: velocity estimate for each pixel (nrows x ncols)
    """
    ifg_data, mst, ncols, nrows, bperp_data, ifg_time_span = _per_tile_setup(ifgs, bperp)
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
                bperp_pix = bperp_data[sel, row, col]
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
    dem_error_correction = np.multiply(dem_error, bperp_data)
    # calculate metric difference to DEM by multiplying the estimate with the per-pixel geometry
    # (i.e. range distance and look angle, see Eq. (2.4.12) in Hanssen (2001))
    # also scale by -0.001 since the phase observations are in mm with positive values away from the sensor
    dem_error = np.multiply(dem_error, np.multiply(range_dist, np.sin(look_angle))) * (-0.001)

    return dem_error, dem_error_correction, vel


def _per_tile_setup(ifgs, bperp):
    """
    Convenience function for setting up DEM error computation parameters

    :param list ifgs: list of interferogram class objects.
    :param ndarray bperp: Per-pixel perpendicular baseline for each interferogram
    """
    if len(ifgs) < 1:
        msg = 'Time series requires 2+ interferograms'
        raise TimeSeriesError(msg)

    nrows = ifgs[0].nrows
    ncols = ifgs[0].ncols
    nifgs = len(ifgs)
    ifg_data = np.zeros((nifgs, nrows, ncols), dtype=np.float32)
    bperp_data = np.zeros((nifgs, nrows, ncols), dtype=np.float32)
    for ifg_num in range(nifgs):
        ifg_data[ifg_num] = ifgs[ifg_num].phase_data
        bperp_data[ifg_num] = bperp[ifg_num]
    mst = ~np.isnan(ifg_data)

    ifg_time_span = np.zeros((nifgs))
    for ifg_num in range(nifgs):
        ifg_time_span[ifg_num] = ifgs[ifg_num].time_span

    return ifg_data, mst, ncols, nrows, bperp_data, ifg_time_span


def _write_dem_errors(ifg_paths, params, preread_ifgs, tiles):
    """
    Convenience function for writing DEM error (one file) and DEM error correction for each IFG to disc
    """
    # re-assemble tiles and save into dem_error dir
    shape = preread_ifgs[ifg_paths[0]].shape

    # save dem error as geotiff file in out directory
    gt, md, wkt = shared.get_geotiff_header_info(ifg_paths[0])
    md[ifc.DATA_TYPE] = ifc.DEM_ERROR
    dem_error = assemble_tiles(shape, params[cf.TMPDIR], tiles, out_type='dem_error')
    dem_error_file = os.path.join(params[cf.OUT_DIR], 'dem_error.tif')
    shared.remove_file_if_exists(dem_error_file)
    shared.write_output_geotiff(md, gt, wkt, dem_error, dem_error_file, np.nan)

    # loop over all ifgs
    idx = 0
    for ifg_path in ifg_paths:
        ifg = Ifg(ifg_path)
        ifg.open()
        # read dem error correction file from tmpdir
        dem_error_correction_ifg = assemble_tiles(shape, params[cf.TMPDIR], tiles, out_type='dem_error_correction',
                                                  index=idx)
        idx += 1
        dem_error_correction_on_disc = MultiplePaths.dem_error_path(ifg.data_path, params)
        np.save(file=dem_error_correction_on_disc, arr=dem_error_correction_ifg)

        # subtract DEM error from the ifg
        ifg.phase_data -= dem_error_correction_ifg
        _save_dem_error_corrected_phase(ifg)


def __check_and_apply_demerrors_found_on_disc(ifg_paths, params):
    """
    Convenience function to check if DEM error correction files have already been produced in a previous run
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
            _save_dem_error_corrected_phase(ifg)
    return all(d.exists() for d in saved_dem_err_paths)


def _save_dem_error_corrected_phase(ifg):
    """
    Convenience function to update metadata and save latest phase after DEM error correction
    """
    # update geotiff tags after DEM error correction
    ifg.dataset.SetMetadataItem(ifc.PYRATE_DEM_ERROR, ifc.DEM_ERROR_REMOVED)
    ifg.write_modified_phase()
    ifg.close()


class DEMError(Exception):
    """
    Generic exception for DEM errors.
    """
