#!/usr/bin/env python
"""overlaptransfrag.py <transfrag file>
Takes transfrags from different time points and overlaps them to produce gff"""

import sys
from common.feature import Features, Feature

if len(sys.argv) != 2:
    print "Invalid arguments"
    print __doc__
    sys.exit(0)

startref = 10506243
endref = 13506243
#nucleotide level only...

lines = [line.split() for line in open(sys.argv[1])]
features = Features()

count = 0
for words in lines:
    if words:
        if words[0] != 'Time':
            count += 1
            start = int(words[0])
            end = int(words[1])
            if start > startref and start < endref:
                features.addFeature(Feature('Adh', 'manak', ['transfrag'], [[start,end]], None, '+', None, 'Gene ' + str(count)))

features.writeGff(sys.argv[1] + '.out')