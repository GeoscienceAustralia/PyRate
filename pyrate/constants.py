import os
__version__ = "0.4.0"
CLI_DESCRIPTION = """
PyRate workflow: 

    Step 1: conv2tif
    Step 2: prepifg
    Step 3: process
    Step 4: merge 

Refer to https://geoscienceaustralia.github.io/PyRate/usage.html for 
more details.
"""
CONV2TIF = 'conv2tif'
PREPIFG = 'prepifg'
PROCESS = 'process'
MERGE = 'merge'

REF_COLOR_MAP_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "utils", "colormap.txt")
