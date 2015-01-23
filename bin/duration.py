#!/usr/bin/env python
"""duration.py <gff_file> <output_file>

For state length distribution estimation, generates list of lengths of features in gff file for further processing.
"""

from feature import *
import sys

if len(sys.argv) != 3:
    print 'Invalid arguments', sys.argv
    print __doc__
    sys.exit(0)
    
features = Features(sys.argv[1])
outfile = open(sys.argv[2], 'w')

for (start,end) in features.exonList():
    outfile.write(str(end - start + 1) + '\n')