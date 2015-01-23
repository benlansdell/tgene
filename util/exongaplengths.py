#!/usr/bin/env python
"""exongaplength.py <tiling data> <gap list> <annotation>
Script to compute a list of tiling gap lengths that occur within exons. gap list is a file with a list of two columns --
the begin and end of each gap respectively.
"""
import sys
from common.feature import Features
from common.expressionstruct import Expression

if len(sys.argv) != 4:
    print 'Invalid arguments'
    print __doc__
    sys.exit(0)
    
features = Features(sys.argv[3])
expression = Expression(sys.argv[1])

header = expression.expression.keys()[0]
(ref, start, end) = header.split(':')[0:3]
startref, endref = int(start), int(end)

exonstring = features.exonString(ref, startref, endref)
gapstring = [0]*(endref - startref + 1)

for line in open(sys.argv[2]):
    (gstart, gend) = line.split()[0:2]
    gstart, gend = int(gstart), int(gend)
    gapstring[gstart-startref:gend-startref] = [1]*(gend-gstart)

for a in range(len(gapstring)):
    if not (gapstring[a] and exonstring[a] == 'F'):
        gapstring[a] = 0
        
ingap = True
gaplength = 0

for a in gapstring:
    if ingap and a == 0:
        print gaplength
        gaplength = 0
        ingap = False
    if ingap:
        gaplength += 1
    if a == 1:
        ingap = True
        
        