#!/usr/bin/env python
"""codingscore2.py <fasta file> <ghmm cfg> <gff file>
In windows of step evaluates log(P(S|q = coding) / P(S|q = noncoding))
and outputs to file
"""

from common.GHMM import *
from common.feature import *
import fasta, sys

if len(sys.argv) != 4:
    print 'Invalid arguments'
    print __doc__
    sys.exit(0)

head,seq = fasta.load(sys.argv[1])
ghmm = GHMM(sys.argv[2])
features = Features(sys.argv[3])
sequence = sequenceDict(seq, 5)
outfile = open(sys.argv[1] + '.scr', 'w')

posscorecnc = []
posscore = []

for ref in features:
    pos = 0
    for generef in features[ref]:
        gene = features[ref][generef]
        genecoords = gene.coords
        for (start, end) in genecoords:
            coding1 = max(ghmm.content['eintn'].probEmit(sequence, start, end+1, 0, '+'), ghmm.content['eintn'].probEmit(sequence, start, end+1, 1, '+'), ghmm.content['eintn'].probEmit(sequence, start, end+1, 2, '+'))
            coding2 = max(ghmm.content['esing'].probEmit(sequence, start, end+1, 0, '+'), ghmm.content['esing'].probEmit(sequence, start, end+1, 1, '+'), ghmm.content['esing'].probEmit(sequence, start, end+1, 2, '+'))
            coding3 = max(ghmm.content['eterm'].probEmit(sequence, start, end+1, 0, '+'), ghmm.content['eterm'].probEmit(sequence, start, end+1, 1, '+'), ghmm.content['eterm'].probEmit(sequence, start, end+1, 2, '+'))
            coding4 = max(ghmm.content['einit'].probEmit(sequence, start, end+1, 0, '+'), ghmm.content['einit'].probEmit(sequence, start, end+1, 1, '+'), ghmm.content['einit'].probEmit(sequence, start, end+1, 2, '+'))
            noncoding1 = ghmm.content['inter'].probEmit(sequence, start, end+1, 0, '+')
            noncoding2 = max(ghmm.content['intrn'].probEmit(sequence, start, end+1, 0, '+'), ghmm.content['intrn'].probEmit(sequence, start, end+1, 1, '+'), ghmm.content['intrn'].probEmit(sequence, start, end+1, 2, '+'))
            coding = max(coding1, coding2, coding3, coding4)
            noncoding = max(noncoding1, noncoding2)
            score2 = coding-noncoding
            coding1 = max(ghmm.content['eintn'].probEmit(sequence, pos, start, 0, '+'), ghmm.content['eintn'].probEmit(sequence, pos, start, 1, '+'), ghmm.content['eintn'].probEmit(sequence, pos, start, 2, '+'))
            coding2 = max(ghmm.content['esing'].probEmit(sequence, pos, start, 0, '+'), ghmm.content['esing'].probEmit(sequence, pos, start, 1, '+'), ghmm.content['esing'].probEmit(sequence, pos, start, 2, '+'))
            coding3 = max(ghmm.content['eterm'].probEmit(sequence, pos, start, 0, '+'), ghmm.content['eterm'].probEmit(sequence, pos, start, 1, '+'), ghmm.content['eterm'].probEmit(sequence, pos, start, 2, '+'))
            coding4 = max(ghmm.content['einit'].probEmit(sequence, pos, start, 0, '+'), ghmm.content['einit'].probEmit(sequence, pos, start, 1, '+'), ghmm.content['einit'].probEmit(sequence, pos, start, 2, '+'))
            noncoding1 = ghmm.content['inter'].probEmit(sequence, pos, start, 0, '+')
            noncoding2 = max(ghmm.content['intrn'].probEmit(sequence, pos, start, 0, '+'), ghmm.content['intrn'].probEmit(sequence, pos, start, 1, '+'), ghmm.content['intrn'].probEmit(sequence, pos, start, 2, '+'))
            coding = max(coding1, coding2, coding3, coding4)
            noncoding = max(noncoding1, noncoding2)
            score1 = coding-noncoding
            posscore.append((pos, score1/(start-pos)))
            posscore.append((start, score1/(start-pos)))
            posscore.append((start, score2/(end-start)))
            posscore.append((end+1, score2/(end-start)))
            pos = end+1
start = len(sequence)
coding1 = max(ghmm.content['eintn'].probEmit(sequence, pos, start, 0, '+'), ghmm.content['eintn'].probEmit(sequence, pos, start, 1, '+'), ghmm.content['eintn'].probEmit(sequence, pos, start, 2, '+'))
coding2 = max(ghmm.content['esing'].probEmit(sequence, pos, start, 0, '+'), ghmm.content['esing'].probEmit(sequence, pos, start, 1, '+'), ghmm.content['esing'].probEmit(sequence, pos, start, 2, '+'))
coding3 = max(ghmm.content['eterm'].probEmit(sequence, pos, start, 0, '+'), ghmm.content['eterm'].probEmit(sequence, pos, start, 1, '+'), ghmm.content['eterm'].probEmit(sequence, pos, start, 2, '+'))
coding4 = max(ghmm.content['einit'].probEmit(sequence, pos, start, 0, '+'), ghmm.content['einit'].probEmit(sequence, pos, start, 1, '+'), ghmm.content['einit'].probEmit(sequence, pos, start, 2, '+'))
noncoding1 = ghmm.content['inter'].probEmit(sequence, pos, start, 0, '+')
noncoding2 = max(ghmm.content['intrn'].probEmit(sequence, pos, start, 0, '+'), ghmm.content['intrn'].probEmit(sequence, pos, start, 1, '+'), ghmm.content['intrn'].probEmit(sequence, pos, start, 2, '+'))
coding = max(coding1, coding2, coding3, coding4)
noncoding = max(noncoding1, noncoding2)
score1 = coding-noncoding
posscore.append((pos, score1/(start-pos)))
posscore.append((start, score1/(start-pos)))

for (a, score) in posscore:
    outfile.write(str(a) + '\t' + str(score) + '\n')