#!/usr/bin/env python
"""averageprobes.py <intensities> ...
Script to average a set of probe intensities across replicates. Assumes there are 3 replicates -> number of input files must be a multiple of three
"""
import sys

length = len(sys.argv) - 1
if length%3:
    print 'Invalid arguments'
    print __doc__
    sys.exit(0)

lines = [open(file).readlines() for file in sys.argv[1:]]

for file_idx in range(0, length, 3):
    outfile = open(sys.argv[file_idx+1] + str(file_idx), 'w')
    for pos_idx in range(len(lines[file_idx])):
        position = lines[file_idx][pos_idx].split()[0]
        intensity = sum([float(lines[a][pos_idx].split()[1]) for a in range(file_idx, file_idx+3)])/3
        outfile.write(position + '\t' + str(intensity) + '\n')
    outfile.close()