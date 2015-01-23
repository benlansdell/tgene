#!/usr/bin/env python
"""parseintsoutput.py <ints file> <start> <end>
Script to parse intensities from aroma to only include a specific region and to remove repeated probes.
"""
import sys

if len(sys.argv) != 4:
    print 'Invalid arguments'
    print __doc__
    sys.exit(0)

words = [line.split() for line in open(sys.argv[1])]
output = open(sys.argv[1] + '.adh', 'w')
start, end = int(sys.argv[2]), int(sys.argv[3])
currentpos = 0

for (str_position, str_intensity) in words:
    position = int(str_position)+13506243
    if position != currentpos and position > start and position < end:
        output.write(str(position-start) + '\t' + str_intensity + '\n')
    currentpos = position