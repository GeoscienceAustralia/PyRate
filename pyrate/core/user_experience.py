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


def break_number_into_factors(n, memo={}, left=2):
    if (n, left) in memo:
        return memo[(n, left)]
    if left == 1:
        return (n, [n])
    i = 2
    best = n
    bestTuple = [n]
    while i * i <= n:
        if n % i == 0:
            rem = break_number_into_factors(n / i, memo, left - 1)
            if rem[0] + i < best:
                best = rem[0] + i
                bestTuple = [i] + rem[1]
        i += 1
    return bestTuple