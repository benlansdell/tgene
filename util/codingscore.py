#!/usr/bin/env python
"""codingscore.py <fasta file> <ghmm cfg> <gff file>
In windows of step evaluates log(P(S|q = coding) / P(S|q = noncoding))
and outputs to file
"""

from common.GHMM import *
from common.feature import *
import fasta, sys

def smoothPoints(score):
    newscore = []
    for a in range(smoothingwindow, len(score)-smoothingwindow):
        scrsum = sum([score[a-smoothingwindow:a+smoothingwindow+1][i][1] for i in range(2*smoothingwindow+1)])
        newscr = scrsum / (2*smoothingwindow+1)
        newscore.append((score[a][0], newscr))
    return newscore

if len(sys.argv) != 4:
    print 'Invalid arguments'
    print __doc__
    sys.exit(0)

step = 50
smoothingwindow = 8
head,seq = fasta.load(sys.argv[1])
ghmm = GHMM(sys.argv[2])
features = Features(sys.argv[3])
sequence = SequenceDict(seq, 5)
outfile = open(sys.argv[1] + '.scr', 'w')

posscorecnc = []
posscore = []

for a in range(0, len(seq)-step, step):
    coding1 = max(ghmm.content['eintn'].probEmit(sequence, a, a+step, 0, '+'), ghmm.content['eintn'].probEmit(sequence, a, a+step, 1, '+'), ghmm.content['eintn'].probEmit(sequence, a, a+step, 2, '+'))
    coding2 = max(ghmm.content['esing'].probEmit(sequence, a, a+step, 0, '+'), ghmm.content['esing'].probEmit(sequence, a, a+step, 1, '+'), ghmm.content['esing'].probEmit(sequence, a, a+step, 2, '+'))
    coding3 = max(ghmm.content['eterm'].probEmit(sequence, a, a+step, 0, '+'), ghmm.content['eterm'].probEmit(sequence, a, a+step, 1, '+'), ghmm.content['eterm'].probEmit(sequence, a, a+step, 2, '+'))
    coding4 = max(ghmm.content['einit'].probEmit(sequence, a, a+step, 0, '+'), ghmm.content['einit'].probEmit(sequence, a, a+step, 1, '+'), ghmm.content['einit'].probEmit(sequence, a, a+step, 2, '+'))
    noncoding1 = ghmm.content['inter'].probEmit(sequence, a, a+step, 0, '+')
    noncoding2 = max(ghmm.content['intrn'].probEmit(sequence, a, a+step, 0, '+'), ghmm.content['intrn'].probEmit(sequence, a, a+step, 1, '+'), ghmm.content['intrn'].probEmit(sequence, a, a+step, 2, '+'))
    coding = max(coding1, coding2, coding3, coding4)
    noncoding = max(noncoding1, noncoding2)
    score = coding-noncoding
    posscorecnc.append((a, coding, noncoding))
    posscore.append((a, score))

posscore = smoothPoints(posscore)

for ref in features:
    for gene in features[ref]:
        if features[ref][gene].strand == '+':
            for (start,end) in features[ref][gene].coords:
                posscore.append((start, 10))
                posscore.append((end, 10))
posscore.sort()

for (a, score) in posscore:
    outfile.write(str(a) + '\t' + str(score) + '\n')