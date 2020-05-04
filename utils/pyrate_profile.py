# This Python module is part of the PyRate software package
#
# Copyright 2020 Geoscience Australia
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
This script creates the pycallgraph file output_file='pyrate_with_roipac.png'.
This script can be run from the 'PyRate' directory.
Change the name of the config file as required.
"""

import sys
from pycallgraph import PyCallGraph
from pycallgraph import Config
from pycallgraph import GlobbingFilter
from pycallgraph.output import GraphvizOutput
from pyrate import process

config = Config()
config.trace_filter = GlobbingFilter(exclude=[
    'pycallgraph.*',
    '*.secret_function',
])

graphviz = GraphvizOutput(output_file='pyrate_with_roipac.png')
config = Config(max_depth=6, groups=False, threaded=True)

# sys.argv[0]: name of this script
# sys.argv[1]: name of the config file
sys.argv = ['pyrate_profile.py', 'pyrate.conf']

with PyCallGraph(output=graphviz, config=config):
    process.main()
