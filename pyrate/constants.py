from pathlib import Path

PYRATEPATH = Path(__file__).parent.parent


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

REF_COLOR_MAP_PATH = PYRATEPATH.joinpath("utils", "colormap.txt")
