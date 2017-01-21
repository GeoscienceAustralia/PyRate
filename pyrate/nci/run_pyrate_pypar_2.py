from __future__ import print_function
import os
import sys
from operator import itemgetter
import numpy as np
from collections import namedtuple
from osgeo import gdal

from pyrate import config as cf
from pyrate import ifgconstants as ifc
from pyrate import reference_phase_estimation as rpe
from pyrate import shared
from pyrate import vcm as vcm_module
from pyrate.nci.parallel import Parallel
from pyrate.scripts import run_pyrate
from pyrate.scripts.run_pyrate import write_msg
from pyrate import orbital
gdal.SetCacheMax(64)


# Constants
MASTER_PROCESS = 0
data_path = 'DATAPATH'

PrereadIfg = namedtuple('PrereadIfg',
                        'path nan_fraction master slave time_span')


def main(params, config_file=sys.argv[1]):

    # setup paths
    xlks, ylks, crop = run_pyrate.transform_params(params)
    base_unw_paths = run_pyrate.original_ifg_paths(params[cf.IFG_FILE_LIST])
    dest_tifs = run_pyrate.get_dest_paths(base_unw_paths, crop, params, xlks)

    # Setting up parallelisation
    parallel = Parallel(True)
    rank = parallel.rank
    ref_pixel_file = os.path.join(params[cf.OUT_DIR], 'ref_pixel.npy')
    refpx, refpy = np.load(ref_pixel_file)
    print('Found reference pixel', refpx, refpy)

    orb_fit_calc_mpi(dest_tifs, parallel, params)

    parallel.barrier()
    output_dir = params[cf.OUT_DIR]

    # save phase data and phase_sum used in the reference phase estimation
    if rank == MASTER_PROCESS:
        phase_sum = 0
        for d in dest_tifs:
            ifg = shared.Ifg(d)
            ifg.open()
            ifg.nodata_value = params[cf.NO_DATA_VALUE]
            phase_sum += ifg.phase_data
            ifg.save_numpy_phase(numpy_file=os.path.join(
                output_dir, os.path.basename(d).split('.')[0] + '.npy'))
            ifg.close()
        comp = np.isnan(phase_sum)  # this is the same as in Matlab
        comp = np.ravel(comp, order='F')  # this is the same as in Matlab
        np.save(file=os.path.join(output_dir, 'comp.npy'), arr=comp)

    parallel.barrier()
    # estimate and remove reference phase
    ref_phase_estimation_mpi(rank, dest_tifs, parallel, params, refpx, refpy)

    parallel.barrier()
    maxvar_vcm_mpi(rank, dest_tifs, parallel, params)

    parallel.finalize()


def ref_phase_estimation_mpi(MPI_myID, ifg_paths, parallel, params,
                             refpx, refpy):
    # TODO: may benefit from tiling and using a median of median algorithms
    write_msg('Finding and removing reference phase')
    num_processors = parallel.size
    no_ifgs = len(ifg_paths)
    process_indices = parallel.calc_indices(no_ifgs)
    process_ifgs = [itemgetter(p)(ifg_paths) for p in process_indices]
    process_ref_phs = np.zeros(len(process_ifgs))
    output_dir = params[cf.OUT_DIR]
    if params[cf.REF_EST_METHOD] == 1:
        for n, p in enumerate(process_ifgs):
            process_ref_phs[n] = ref_phase_method1_dummy(p, output_dir)
            print('finished processing {} of process total {}, ' \
                  'of overall {}'.format(n, len(process_ifgs), no_ifgs))

    elif params[cf.REF_EST_METHOD] == 2:
        for n, p in enumerate(process_ifgs):
            process_ref_phs[n] = ref_phase_method2_dummy(params, p,
                                                         refpx, refpy)
            print('finished processing {} of process total {}, ' \
                  'of overall {}'.format(n, len(process_ifgs), no_ifgs))
    else:
        raise cf.ConfigException('Ref phase estimation method must be 1 or 2')

    ref_phs_file = os.path.join(params[cf.OUT_DIR], 'ref_phs.npy')

    if MPI_myID == MASTER_PROCESS:
        all_indices = parallel.calc_all_indices(no_ifgs)
        ref_phs = np.zeros(no_ifgs)
        ref_phs[process_indices] = process_ref_phs

        for i in range(1, num_processors):
            process_ref_phs = \
                parallel.receive(source=i, tag=-1,
                                 return_status=False)
            ref_phs[all_indices[i]] = process_ref_phs
        np.save(file=ref_phs_file, arr=ref_phs)
    else:
        # send reference phase data to master process
        parallel.send(process_ref_phs, destination=MASTER_PROCESS,
                      tag=MPI_myID)


def ref_phase_method2_dummy(params, ifg_path, refpx, refpy):
    half_chip_size = int(np.floor(params[cf.REF_CHIP_SIZE] / 2.0))
    chipsize = 2 * half_chip_size + 1
    thresh = chipsize * chipsize * params[cf.REF_MIN_FRAC]
    output_dir = params[cf.OUT_DIR]
    numpy_file = os.path.join(
        output_dir, os.path.basename(ifg_path).split('.')[0] + '.npy')
    phase_data = np.load(numpy_file)
    ref_phs = rpe.est_ref_phase_method2_multi(phase_data,
                                              half_chip_size,
                                              refpx, refpy, thresh)
    phase_data -= ref_phs
    ifg = shared.Ifg(ifg_path)
    ifg.open()
    md = ifg.meta_data
    md[ifc.REF_PHASE] = ifc.REF_PHASE_REMOVED
    ifg.write_modified_phase(data=phase_data)
    ifg.close()
    return ref_phs


def ref_phase_method1_dummy(ifg_path, output_dir):
    comp_file = os.path.join(output_dir, 'comp.npy')
    comp = np.load(comp_file)
    numpy_file = os.path.join(
        output_dir, os.path.basename(ifg_path).split('.')[0] + '.npy')
    phase_data = np.load(numpy_file)
    ref_phs = rpe.est_ref_phase_method1_multi(phase_data, comp)
    phase_data -= ref_phs
    ifg = shared.Ifg(ifg_path)
    ifg.open()
    md = ifg.meta_data
    md[ifc.REF_PHASE] = ifc.REF_PHASE_REMOVED
    ifg.write_modified_phase(data=phase_data)
    ifg.close()
    return ref_phs


def maxvar_vcm_mpi(rank, ifg_paths, parallel, params):
    print('calculating maxvar and vcm')
    num_processors = parallel.size
    no_ifgs = len(ifg_paths)
    process_indices = parallel.calc_indices(no_ifgs)
    process_ifgs = [itemgetter(p)(ifg_paths) for p in process_indices]
    process_maxvar = []
    for n, i in enumerate(process_ifgs):
        print('calculating maxvar for {} of process ifgs {} of ' \
              'total {}'.format(n, len(process_ifgs), no_ifgs))
        # TODO: cvd calculation is still pretty slow - revisit
        process_maxvar.append(vcm_module.cvd(i, params)[0])
    maxvar_file = os.path.join(params[cf.OUT_DIR], 'maxvar.npy')
    vcmt_file = os.path.join(params[cf.OUT_DIR], 'vcmt.npy')
    if rank == MASTER_PROCESS:
        all_indices = parallel.calc_all_indices(no_ifgs)
        maxvar = np.empty(no_ifgs)
        maxvar[all_indices[MASTER_PROCESS]] = process_maxvar

        for i in range(1, num_processors):
            process_maxvar = parallel.receive(source=i, tag=-1,
                                              return_status=False)
            maxvar[all_indices[i]] = process_maxvar
        np.save(file=maxvar_file, arr=maxvar)
    else:
        parallel.send(process_maxvar, destination=MASTER_PROCESS, tag=rank)
    parallel.barrier()
    maxvar = np.load(maxvar_file)
    ifgs = shared.prepare_ifgs_without_phase(ifg_paths, params)
    vcmt = vcm_module.get_vcmt(ifgs, maxvar)
    if rank == MASTER_PROCESS:
        np.save(file=vcmt_file, arr=vcmt)

    return maxvar, vcmt


def orb_fit_calc_mpi(ifg_paths, parallel, params):
    print('calculating orbfit correction')
    if params[cf.ORBITAL_FIT_METHOD] != 1:
        raise cf.ConfigException('For now orbfit method must be 1')

    # ifgs = shared.prepare_ifgs_without_phase(ifg_paths, params)
    no_ifgs = len(ifg_paths)
    process_indices = parallel.calc_indices(no_ifgs)
    process_ifgs = [itemgetter(p)(ifg_paths) for p in process_indices]

    mlooked = None
    # TODO: MPI orbfit method 2
    orbital.orbital_correction(process_ifgs, params, mlooked=mlooked)
    print('finished orbfit calculation for process {}'.format(parallel.rank))


if __name__ == '__main__':
    # read in the config file, and params for the simulation
    params = run_pyrate.get_ifg_paths()[2]
    main(params)
