#!/usr/bin/env python
"""smoothtiling.py <tiling> <header> <offset>
Script to convert discrete outputs for each probe into a smoothed per nucleotide score for use in GHMM.
Tiling should be a textfile list with two columns, position and intensity. Header is the fasta-style header output.
Assumes tiles are not-overlapping.
"""
import sys
from numpy import inf

def interpolate(length):
	retval = ([inf]*29 + [0])*(length/30+2)
	return retval[0:length]
	

if len(sys.argv) != 4:
    print 'Invalid arguments'
    print __doc__
    sys.exit(0)
    
fileout = open(sys.argv[1]+'.out', 'w')
fileout.write('>' + sys.argv[2] + '\n')
offset = int(sys.argv[3])
start = int(sys.argv[2].split(':')[1])
end = int(sys.argv[2].split(':')[2])

probelength = 25
probeintensities = [(int(line.split()[0])+offset, float(line.split()[1])) for line in open(sys.argv[1])]

minpos = min(probeintensities)[0]
maxpos = max(probeintensities)[0]
smoothedintensities = [inf]*(maxpos+probelength-minpos)
oldpos = offset

for (pos, intensity) in probeintensities:
    smoothedintensities[pos-minpos+int(probelength/2):pos-minpos+int(probelength/2)+1] = [intensity]
    if pos - oldpos > 100:
        smoothedintensities[oldpos-minpos:pos-minpos] = interpolate(pos-oldpos)
    oldpos = pos

for a in range(start, end+1):
    try:
        fileout.write(str(smoothedintensities[a-minpos]) + '\n')
    except: pass