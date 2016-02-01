#!/usr/bin/env python
"""
This script creates the pycallgraph file output_file='pyrate_with_roipac.png'.
This script can be run from the 'PyRate' directory.
Change the name of the config file as required.
"""
__author__ = 'Sudipta Basak'
__date_created__ = '1/02/16'

import sys
from pycallgraph import PyCallGraph
from pycallgraph import Config
from pycallgraph import GlobbingFilter
from pycallgraph.output import GraphvizOutput
from pyrate.scripts import run_pyrate

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
    run_pyrate.main()
