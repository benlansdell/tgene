#!/usr/bin/env python
"""smoothtiling.py <tiling> <header>
Script to convert discrete outputs for each probe into a smoothed per nucleotide score for use in GHMM.
Tiling should be a textfile list with two columns, position and intensity. Header is the fasta-style header output.
Assumes tiles are not-overlapping.
"""
import sys

def smoothed(inta, intb, length):
    #Linear interpolant
    #return [inta + (intb-inta)*(float(a)/length) for a in range(length)]
    #Contstant interpolant
    return [float(intb+inta)/2 for a in range(length)]

if len(sys.argv) != 3:
    print 'Invalid arguments'
    print __doc__
    sys.exit(0)
    
fileout = open(sys.argv[1]+'.out', 'w')
fileout.write('>' + sys.argv[2] + '\n')
start = int(sys.argv[2].split(':')[1])
end = int(sys.argv[2].split(':')[2])

probelength = 25
probeintensities = [(int(line.split()[0]), float(line.split()[1])) for line in open(sys.argv[1])]

minpos = min(probeintensities)[0]
maxpos = max(probeintensities)[0]
minint = min([a[1] for a in probeintensities])
#Set to a missing emission score of 0. 
#if minint > 0: minint = 0.0
#print minint
smoothedintensities = [minint]*(maxpos+probelength-minpos)
oldpos = oldint = 0

for (pos, intensity) in probeintensities:
    smoothedintensities[pos-minpos:pos-minpos+probelength] = [intensity]*probelength
#    print oldpos+minpos, pos
#    if oldpos+minpos < pos and oldpos+minpos >= start and pos <= end: print oldpos+minpos, pos
    if oldpos and oldint and pos-minpos - oldpos < 100:
        smoothedintensities[oldpos:pos-minpos] = smoothed(oldint, intensity, pos-minpos-oldpos)
    oldpos = pos-minpos+probelength
    oldint = intensity

for a in range(start, end+1):
    try:
        fileout.write(str(smoothedintensities[a-minpos]) + '\n')
    except: pass