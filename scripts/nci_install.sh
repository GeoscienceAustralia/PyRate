#!/bin/bash
#   This script is part of the PyRate software package.
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

rm -rf  ~/PyRateVenv
rm -rf  ~/PyRate
git clone https://github.com/GeoscienceAustralia/PyRate.git --branch test --single-branch ~/PyRate
module purge
module load python3/3.7.4
module load gdal/3.0.2
module load openmpi/2.1.6
export PYTHONPATH=/apps/gdal/3.0.2/lib64/python3.7/site-packages:$PYTHONPATH
python3 -m venv ~/PyRateVenv
source ~/PyRateVenv/bin/activate
cd ~/PyRate
python setup.py install
