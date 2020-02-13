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
import pathlib


def delete_tsincr_files(params):
    """Args:
    
    Args:

    Args:
      params: 

    Returns:
      

    """
    outDirPath = pathlib.Path(params["outdir"])
    for filePath in outDirPath.iterdir():
        if "tsincr" in str(filePath):
            filePath.unlink()


def break_number_into_factors(n, memo={}, left=2):
    """

    Args:
      n: param memo:  (Default value = {})
      left: Default value = 2)
      memo: (Default value = {})

    Returns:

    """
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

    # handle edge case when only one processor is available
    if bestTuple == [4]:
        return 2, 2

    if len(bestTuple) == 1:
        bestTuple.append(1)
        return bestTuple

    return bestTuple


def add_ref_pixel_location_to_interferogram_metadata():
    """ """
    None
