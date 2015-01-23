#!/usr/bin/env python
"""combinefiles.py <output_file> <input_file> ...

Append files listed into one file for uploading to UCSC genome browser.
"""

import sys, os
from fnmatch import fnmatch

if len(sys.argv) < 3:
	print __doc__
	sys.exit(0)

outlines = []
outfile = open(sys.argv[1], 'w')

for file in sys.argv[2:]:
	outlines += open(file).readlines()
	
for line in outlines:
	outfile.write(line)