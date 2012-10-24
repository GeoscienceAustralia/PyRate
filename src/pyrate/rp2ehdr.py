#!/usr/bin/env python
# Quick and dirty script to generate EHdr format header files from ROIPAC data

import os, sys
from pyrate.roipac import to_ehdr_header

usage = "Usage: rp2ehdr.py [geo****.unw.rsc] file" 
if len(sys.argv) != 2:
	sys.stderr.write(usage)

to_ehdr_header(sys.argv[1])
