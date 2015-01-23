#!/usr/bin/env python
"""genenamesgenbank.py <genbank file>
Script to output genenames from a genbank file"""

import sys

if len(sys.argv) != 2:
	print 'Invalid arguments'
	print __doc__
	sys.exit(0)
	
file = open(sys.argv[1])
for line in file:
	if line:
		words = line.split()
		if len(words) > 1:
			if words[0] == 'LOCUS':
				print words[1]