import os, luigi, pickle
import pyrate.config as config
from pyrate.prepifg import (
    Ifg,
    getAnalysisExtent,
    mlooked_path,
    prepare_ifg,
    PreprocessError)
from pyrate.tasks.converttogeotif import ConvertToGeotiff
from pyrate.tasks.utils import (
    IfgListMixin,
    InputParam,
    RasterParam)
from pyrate.scripts.run_pyrate import warp_required


class GetAnalysisExtents(IfgListMixin, luigi.Task):
    crop_opt = luigi.IntParameter(config_path=InputParam(config.IFG_CROP_OPT))
    ifgx_first = luigi.FloatParameter(default=None,
                                     config_path=InputParam(config.IFG_XFIRST))
    ifgy_first = luigi.FloatParameter(default=None,
                                     config_path=InputParam(config.IFG_YFIRST))
    ifgx_last = luigi.FloatParameter(default=None,
                                    config_path=InputParam(config.IFG_XLAST))
    ifgy_last = luigi.FloatParameter(default=None,
                                    config_path=InputParam(config.IFG_YLAST))
    xlooks = luigi.IntParameter(config_path=InputParam(config.IFG_LKSX))
    ylooks = luigi.IntParameter(config_path=InputParam(config.IFG_LKSY))

    def requires(self):
        return [ConvertToGeotiff()]

    def run(self):
        userExts = (self.ifgx_first, self.ifgy_first, self.ifgx_last, self.ifgy_last)

        if not all(userExts):
            if self.crop_opt == 3:
                raise PreprocessError('No custom cropping extents specified')
            userExts = None

        ifgs = [Ifg(path) for path in self.ifgTiffList()]

        extents = getAnalysisExtent(
            self.crop_opt,
            ifgs,
            self.xlooks,
            self.ylooks,
            userExts)

        with open(self.extentsFileName, 'wb') as extFile:
            pickle.dump(extents, extFile)

    def output(self):
        return luigi.file.LocalTarget(self.extentsFileName)


class PrepareInterferogram(IfgListMixin, luigi.WrapperTask):
    """
    Produces multilooked/resampled data files for PyRate analysis.

    :param ifgs: sequence of Ifg objs (DEM obj may be included for processing)
    :param crop_opt: integer cropping type option (see config)
    :param xlooks: multilooking factor for the X axis
    :param ylooks: Y axis multilooking factor
    :param float thresh: (0.0, 1.0). Controls NaN handling when resampling to coarser grids.
         Value is the proportion above which the number of NaNs in an area is
         considered invalid. thresh=0 resamples to NaN if 1 or more contributing
         cells are NaNs. At 0.25, it resamples to NaN if 1/4 or more contributing
         cells are NaNs. At 1.0, areas are resampled to NaN only if all
         contributing cells are NaNs.
    :param user_exts: CustomExts tuple with user sepcified lat long corners
    :param verbose: Controls level of gdalwarp output
    """

    ifg = RasterParam()
    thresh = luigi.FloatParameter(config_path=InputParam(
        config.NO_DATA_AVERAGING_THRESHOLD))
    crop_opt = luigi.IntParameter(config_path=InputParam(config.IFG_CROP_OPT))
    xlooks = luigi.IntParameter(config_path=InputParam(config.IFG_LKSX))
    ylooks = luigi.IntParameter(config_path=InputParam(config.IFG_LKSY))
    # verbose = luigi.BooleanParameter(default=True, significant=False)

    def requires(self):
        return [GetAnalysisExtents()]

    def run(self):
        with open(self.extentsFileName, 'rb') as extFile:
            extents = pickle.load(extFile)
        prepare_ifg(
            self.ifg.data_path,
            self.xlooks,
            self.ylooks,
            extents,
            self.thresh,
            self.crop_opt
            )
        self.ifg.close()

    def output(self):
        if warp_required(self.xlooks, self.ylooks, self.crop_opt):
            return luigi.file.LocalTarget(
                mlooked_path(self.ifg.data_path, self.ylooks, self.crop_opt))
        else:
            return []

    def complete(self):
        if self.output():
            return super(luigi.WrapperTask, self).complete()
        else:
            # then this didn't do anything, so check that the requres are complete
            # ... which is exactly what luigi.WrapperTask does.
            # TODO: This requires knowledge of prepare_ifg, that is opaque. Address with refactoring.
            return super(PrepareInterferogram, self).complete()


class PrepareInterferograms(IfgListMixin, luigi.WrapperTask):
    def __init__(self, *args, **kwargs):
        super(PrepareInterferograms, self).__init__(*args, **kwargs)
        self.extentsRemoved = False

    def requires(self):
        return [PrepareInterferogram(ifg=Ifg(path)) for path in self.ifgTiffList()]

    def run(self):
        try:
            if os.path.exists(self.extentsFileName):
                os.remove(self.extentsFileName)
        except:
            raise PrepifgException(
                'Extents file was not found in the desired '
                'location: {}'.format(self.extentsFileName),
                'Make sure your paths are setup correctly in config file')

        self.extentsRemoved = True

    def complete(self):
        return self.extentsRemoved and super(PrepareInterferograms, self).complete()


class PrepifgException(Exception):

    pass

