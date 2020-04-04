import os
import pathlib
from pyrate.core import config as cf
from pyrate import CONV2TIF, PREPIFG, PROCESS, MERGE  # Step names


def delete_tsincr_files(params):

    outDirPath = pathlib.Path(params["outdir"])
    for filePath in outDirPath.iterdir():
        if "tsincr" in str(filePath):
            filePath.unlink()


if __name__ == "__main__":
    config_file = "input_parameters.conf"
    config_file = os.path.abspath(config_file)
    params = cf.get_config_params(config_file, step=MERGE)
    delete_tsincr_files(params)
