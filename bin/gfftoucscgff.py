#!/usr/bin/env python
"""gfftoucscgff.py <gff_file> [position]

Script to output gff file with correct name and format for loading into UCSC genome browser as 
a custom track.

position is an optional string to set the position of the genome browser.
Defaults to chr2L:13500000-14000000

Returns a file with name gff_file.out
"""
import sys
from feature import *

if len(sys.argv) != 3 or len(sys.argv) != 2:
	print __doc__
	sys.exit(0)
	
if len(sys.argv) == 3:
	header = 'browser position ' + sys.argv[2] + '\nbrowser hide all\ntrack name='
else:
	header = 'browser position chr2L:13500000-14000000\nbrowser hide all\ntrack name='

infile = [(line.split(None, 1)[0], line.split(None, 1)[1]) for line in open(sys.argv[1])]
outfile = open(sys.argv[1].split('/')[-1] + '.out', 'w')
exons, nc = '',''

for (ref, rest) in infile:
	if rest.split()[1] == 'exon':
		exons += ('chr2L\t' + rest[0:-1] + ref + '\n')
	else:
		nc += ('chr2L\t' + rest[0:-1] + ref + '\n')

if exons:
	outfile.write(header + sys.argv[1] + 'exon description="' + sys.argv[1] + 'exon" visibility=2\n' + exons)
if nc:
	outfile.write(header + sys.argv[1] + 'nc description="' + sys.argv[1] + 'nc" visibility=2\n' + nc)
