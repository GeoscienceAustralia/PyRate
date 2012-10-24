#!/usr/bin/env python
# Quick and dirty script to generate EHdr format header files from ROIPAC data
# Ben Davies, ANUSF/NCI  24/10/12

import sys
from pyrate.roipac import to_ehdr_header


usage = "Usage: rp2ehdr.py [ROIPAC unw.rsc files ...]\n" 
if len(sys.argv) < 2:
	sys.stderr.write(usage)
	sys.exit()

for path in sys.argv[1:]:
	to_ehdr_header(path)
