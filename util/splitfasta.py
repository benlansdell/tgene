#!/usr/bin/env python
"""splitfasta.py <fasta file> <intensity file> <gff file>
Script to split (single sequence)fasta file into smaller components, such that doesn't cut gene in half."""

import sys, fasta
from common.feature import Features
from common.expressionstruct import Expression

if len(sys.argv) != 4:
    print "Invalid arguments"
    print __doc__
    sys.exit(0)

ncomponents = 20
features = Features(sys.argv[3])
head,seq = fasta.load(sys.argv[1])
startref = int(head.split(':')[1])
expression = Expression(sys.argv[2])
keys = expression.expression.keys()
intensities = expression[keys[0]]

splitpoints = []

outfile = open(sys.argv[2][:-3] + '.fly.combined', 'w')


for ref in features:
    currentpos= 0
    #Find long enough split points
    for generef in features[ref]:
        if currentpos and features[ref][generef].min - currentpos > 2000:
            splitpoints.append((features[ref][generef].min + currentpos)/2)
        currentpos = features[ref][generef].max
    #Split sequence and write to files
    numsplits = len(splitpoints)
    points = [splitpoints[int(n*1.0*numsplits/ncomponents)] for n in range(1,ncomponents)]
    points = [startref] + points + [startref + len(seq)]
    print points
    for n in range(ncomponents):
        seqn = seq[points[n]-startref:points[n+1]-startref]
        intn = intensities[points[n]-startref:points[n+1]-startref]
        headn = head.split(':')[0] + ':' + str(points[n]) + ':' + str(points[n]+len(seqn)-1)
        fasta.save(sys.argv[1][:-3] + '.' + str(n) + '.' + str(ncomponents) + '.fa', headn, seqn)
        outfile.write('>' + headn + '\n')
        [outfile.write(str(a) + '\n') for a in intn]