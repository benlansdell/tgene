#!/usr/bin/env python
"""inttowig.py <intensity file>

Script to output intensities to .wig format for use in UCSC browser
"""

import sys

if len(sys.argv) != 2:
	print __doc__
	sys.exit(0)

infile = open(sys.argv[1]).readlines()
outfile = open(sys.argv[1].split('/')[-1] + '.wig', 'w')
start = 13506243
end = 16425263
#	200 base wide points graph at every 300 bases, 50 pixel high graph
#	autoScale off and viewing range set to [0:1000]
#	priority = 30 positions this as the third graph
#	Note, one-relative coordinate system in use for this format
header = 'track type=wiggle_0 name="SMC data" description="SMC data from tiling arrays" visibility=full autoScale=off viewLimits=-50:100 color=0,200,100 maxHeightPixels=100:50:20 graphType=points priority=30\nfixedStep chrom=chr2L start=13506243 step=10 span=5\n'
outfile.write(header)
pos = 0
for line in infile:
	if line:
		if line[0] != '>':
			if not pos%10:
				outfile.write(line)
			pos += 1
