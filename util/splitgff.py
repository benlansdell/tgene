#!/usr/bin/env python
"""splitgff.py <list file> <gff file>
Script to split (single reference)gff file into smaller components, coordinates given as a two column list in 'list file'."""

import sys
from common.feature import Features

if len(sys.argv) != 3:
    print "Invalid arguments"
    print __doc__
    sys.exit(0)

features = Features(sys.argv[2])
coords = [(int(line.split()[0]), int(line.split()[1])) for line in open(sys.argv[1])]

count = 0
for (start, end) in coords:
    newfeatures = Features()
    for ref in features:
        for gene in features[ref]:
            if features[ref][gene].min >= start and features[ref][gene].max <= end:
                feature = features[ref][gene]
                feature.ref = ref + str(count)# + ':' + str(start) + ':' + str(end) 
                newfeatures.addFeature(feature)
#    newfeatures.writeGff(sys.argv[2] + '.' + str(count), 'a')
    newfeatures.writeGff(sys.argv[2] + '.combined', None, 'a')
    count += 1
