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
from numpy.linalg import inv, LinAlgError

from pyrate.core import geometry, shared, mpiops, config as cf, ifgconstants as ifc
from pyrate.core.logger import pyratelogger as log
from pyrate.core.shared import Ifg, Geometry
from pyrate.configuration import Configuration
from pyrate.configuration import MultiplePaths
from pyrate.merge import assemble_tiles


def dem_error_calc_wrapper(params: dict) -> None:
    """
    MPI wrapper for DEM error correction
    """
    if params[cf.DEMERROR]:
        log.info('Doing DEM error correction')
    else:
        log.info('DEM error correction not required')
        return
    # geometry information needed to calculate Bperp for each pixel using first IFG in list
    ifg_paths = [ifg_path.tmp_sampled_path for ifg_path in params[cf.INTERFEROGRAM_FILES]]
    ifg0_path = ifg_paths[0]
    ifg0 = Ifg(ifg0_path)
    ifg0.open(readonly=True)

    # not currently implemented for ROIPAC data which breaks some tests
    # if statement can be deleted once ROIPAC is deprecated from PyRate
    if not ifg0.meta_data[ifc.PYRATE_INSAR_PROCESSOR] == 'ROIPAC':

        if params[cf.BASE_FILE_LIST] is None:
            msg = f"No baseline files supplied: DEM error correction not computed"
            raise DEMError(msg)

        if params[cf.LT_FILE] is None:
            msg = f"No lookup table file supplied: DEM error correction not computed"
            raise DEMError(msg)

        log.info('Calculating per-pixel baseline')

        # read radar azimuth and range tif files
        rdc_az_file = join(params[cf.OUT_DIR], 'rdc_azimuth.tif')
        geom_az = Geometry(rdc_az_file)
        geom_az.open(readonly=True)
        az = geom_az.geometry_data
        rdc_rg_file = join(params[cf.OUT_DIR], 'rdc_range.tif')
        geom_rg = Geometry(rdc_rg_file)
        geom_rg.open(readonly=True)
        rg = geom_rg.geometry_data

        # split into tiles to calculate DEM error correction
        if not Configuration.vcmt_path(params).exists():
            raise FileNotFoundError("VCMT is not found on disc. Have you run the 'correct' step?")
        params[cf.PREREAD_IFGS] = cp.load(open(Configuration.preread_ifgs(params), 'rb'))
        params[cf.VCMT] = np.load(Configuration.vcmt_path(params))
        params[cf.TILES] = Configuration.get_tiles(params)
        tiles = params[cf.TILES]
        preread_ifgs = params[cf.PREREAD_IFGS]
        vcmt = params[cf.VCMT]
        threshold = params[cf.DE_PTHR]

        # read lon and lat values of multi-looked ifg (first ifg only)
        lon, lat = geometry.get_lonlat_coords(ifg0)
        # cut rg and az to tile size

    # the following code is a quick way to do the bperp calculation, but is not identical to the GAMMA output
    # where the near range of the first SLC is used for each pair.
    # calculate look angle for interferograms (using the Near Range of the first SLC)
    #look_angle = geometry.calc_local_geometry(ifg0, None, rg, lon, lat, params)
    #bperp = geometry.calc_local_baseline(ifg0, az, look_angle, params)
    #print(bperp.shape)

        # process in tiles
        process_tiles = mpiops.array_split(tiles)
        for t in process_tiles:
            ifg_parts = [shared.IfgPart(p, t, preread_ifgs, params) for p in ifg_paths]
            lon_parts = lon[t.top_left_y:t.bottom_right_y, t.top_left_x:t.bottom_right_x]
            lat_parts = lat[t.top_left_y:t.bottom_right_y, t.top_left_x:t.bottom_right_x]
            az_parts = az[t.top_left_y:t.bottom_right_y, t.top_left_x:t.bottom_right_x]
            rg_parts = rg[t.top_left_y:t.bottom_right_y, t.top_left_x:t.bottom_right_x]

            nifgs = len(ifg_paths)
            bperp = np.empty((nifgs, lon_parts.shape[0], lon_parts.shape[1])) * np.nan
            ifg_num = 0
            # calculate per-pixel perpendicular baseline for each IFG
            for ifg_path in ifg_paths: # loop could be avoided by approximating the look angle for the first Ifg
                ifg = Ifg(ifg_path)
                ifg.open(readonly=True)
                # calculate look angle for interferograms (using the Near Range of the primary SLC)
                look_angle = geometry.calc_local_geometry(ifg, None, rg_parts, lon_parts, lat_parts, params)
                bperp[ifg_num, :, :] = geometry.calc_local_baseline(ifg, az_parts, look_angle, params)
                ifg_num += 1


            log.debug('Calculating DEM error for tile {} during DEM error correction'.format(t.index))
            #mst_tile = np.load(Configuration.mst_path(params, t.index))
            dem_error, dem_error_correction = calc_dem_errors(ifg_parts, bperp, threshold, vcmt)
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
            np.save(file=os.path.join(params[cf.TMPDIR], 'dem_error_correction_{}.npy'.format(t.index)), \
                    arr=tmp_array)

        # wait for all processes to finish
        mpiops.comm.barrier()

        # write dem error and correction values to file
        mpiops.run_once(_write_dem_errors, ifg_paths, params, preread_ifgs, tiles)


    log.info('Finished DEM error correction')


def calc_dem_errors(ifgs, bperp, threshold, vcmt):
    """
    Calculates the per-pixel DEM error using least-squares adjustment of phase data and
    perpendicular baseline. The least-squares adjustment co-estimates the velocities.

        - *nrows* is the number of rows in the ifgs,
        - *ncols* is the  number of columns in the ifgs, and
        - *nepochs* is the number of unique epochs (dates)

    :param list ifgs: list of interferogram class objects.
    :param dict params: Dictionary of configuration parameters
    :param ndarray bperp: Per-pixel perpendicular baseline for each interferogram
    :param ndarray vcmt: Positive definite temporal variance covariance matrix
    :param ndarray mst: [optional] Minimum spanning tree array

    :return: ndarray dem_error: estimated per-pixel dem error (nrows x ncols)
    :return: ndarray dem_error_correction: DEM error correction for each pixel and interferogram (nifgs x nrows x ncols)
    """
    ifg_data, mst, ncols, nrows, bperp_data, ifg_time_span = _perpixel_setup(ifgs, bperp)
    nifgs = ifg_data.shape[0]
    if threshold < 4:
        msg = f"pixel threshold too low (i.e. <4) resulting in singularities in DEM error estimation"
        raise DEMError(msg)
    # pixel-by-pixel calculation
    # preallocate empty arrays for results
    dem_error = np.empty([nrows, ncols], dtype=np.float32) * np.nan
    # nested loops to loop over the 2 image dimensions
    for row in range(nrows):
        for col in range(ncols):
            # calc DEM error for each pixel with valid Bperp and ifg phase data
            # todo use a threshold for the minimum number of ifgs per pixel
            # check pixel for non-redundant ifgs
            sel = np.nonzero(mst[:, row, col])[0]  # trues in mst are chosen
            if len(sel) >= threshold:
                ifgv = ifg_data[sel, row, col]
                bperp_pix = bperp_data[sel, row, col]
                time_span = ifg_time_span[sel]
                # new covariance matrix using actual number of observations
                m = len(sel)
                Qyy = np.eye(m)
                #Qyy[:m, :m] = vcmt[sel, np.vstack(sel)]
                # Design matrix of least-squares system
                A = np.column_stack((np.ones(m), bperp_pix , time_span))
                y = np.vstack(ifgv)
                # solve weighted least-squares system (not available in scipy!)
                try:
                    inv(Qyy) # Var-cov matrix of observations could be singular
                    N = np.matmul(np.matmul(np.transpose(A), inv(Qyy)), A)
                    Qxx = inv(N)
                    xhat = np.matmul(np.matmul(np.matmul(Qxx, np.transpose(A)), inv(Qyy)), y)
                    dem_error[row][col] = xhat[1]
                except LinAlgError as err: # nan value for DEM error in case of singular Qyy matrix
                    if 'Singular matrix' in str(err):
                        dem_error[row][col] = np.nan
    dem_error_correction = np.multiply(dem_error, bperp_data)

    return dem_error, dem_error_correction


def _perpixel_setup(ifgs, bperp):
    """
    Convenience function for setting up time series computation parameters

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
    bperp_data =  np.zeros((nifgs, nrows, ncols), dtype=np.float32)
    for ifg_num in range(nifgs):
         ifg_data[ifg_num] = ifgs[ifg_num].phase_data
         bperp_data[ifg_num] = bperp[ifg_num]
    mst = ~np.isnan(ifg_data)

    ifg_time_span = np.zeros((nifgs))
    for ifg_num in range(nifgs):
        ifg_time_span[ifg_num] = ifgs[ifg_num].time_span

    return ifg_data, mst, ncols, nrows, bperp_data, ifg_time_span


def _write_dem_errors(ifg_paths, params, preread_ifgs, tiles):

    # re-assemble tiles and save into dem_error dir
    shape = preread_ifgs[ifg_paths[0]].shape

    # save dem error as geotiff file in out directory
    gt, md, wkt = shared.get_geotiff_header_info(ifg_paths[0])
    md[ifc.EPOCH_DATE] = None  # needs to have a value in write_output_geotiff
    md[ifc.DATA_TYPE] = ifc.DEM_ERROR
    dem_error = assemble_tiles(shape, params[cf.TMPDIR], tiles, out_type='dem_error', index=None)
    dem_error_file = os.path.join(params[cf.OUT_DIR], 'dem_error.tif')
    geometry.remove_file_if_exists(dem_error_file)
    shared.write_output_geotiff(md, gt, wkt, dem_error, dem_error_file, np.nan)

    # loop over all ifgs
    idx = 0
    for ifg_path in ifg_paths:
        ifg = Ifg(ifg_path)
        ifg.open(readonly=True)
        # read dem error correction file from tmpdir (size
        dem_error_correction_ifg = assemble_tiles(shape, params[cf.TMPDIR], tiles, out_type='dem_error_correction', \
                                                  index=idx)
        idx += 1
        dem_error_correction_on_disc = MultiplePaths.dem_error_path(ifg.data_path, params)
        np.save(file=dem_error_correction_on_disc, arr=dem_error_correction_ifg)


def _remove_file_if_exists(file):
    """
        Function to remove a geometry file if it already exists.
    """
    try:
        os.remove(file)
    except OSError:
        pass


class DEMError(Exception):
    """
    Generic exception for DEM errors.
    """