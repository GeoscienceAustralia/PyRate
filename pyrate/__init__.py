# -*- coding: utf-8 -*-
import os

"""
PyRate
------------

Main package for PyRate.
"""


# Turn off MPI warning
os.environ['OMPI_MCA_btl_base_warn_component_unused'] = '0'

# Step name constants
CONV2TIF = 'conv2tif'
PREPIFG = 'prepifg'
PROCESS = 'process'
MERGE = 'merge'
