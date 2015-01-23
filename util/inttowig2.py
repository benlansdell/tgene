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
start = 10506243
end = 13506243
#	200 base wide points graph at every 300 bases, 50 pixel high graph
#	autoScale off and viewing range set to [0:1000]
#	priority = 30 positions this as the third graph
#	Note, one-relative coordinate system in use for this format
header = 'track type=wiggle_0 name="SMC data" description="SMC data from tiling arrays" visibility=full autoScale=off viewLimits=-50:100 color=0,200,100 maxHeightPixels=100:50:20 graphType=points priority=30\nvariableStep chrom=chr2L span=5\n'
outfile.write(header)
for line in infile:
	if line:
		pos = int(line.split()[0])
		val = line.split()[1]
		if pos > start and pos < end:
			outfile.write(str(pos) + '\t' + val + '\n')