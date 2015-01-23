#!/usr/bin/env python
"""
genepred.py <fasta_file> <expression_file> <ghmm_cfg_file> <output_file>

Gene prediction using given sequence, expression, ghmm config file.

fasta_file must be in format:
	>ref:start:end
	ATGCTACG...
	
expression_file must be in format:
	>ref:start:end
	inf
	inf
	23.9
	...
	
ghmm_cgf_file must be in format described in the template.cfg file

Returns a gff file containing predicted coding and non-coding features.
"""

import sys
import fasta
from decoding import *

if len(sys.argv) != 5:
    print 'Invalid arguments'
    print __doc__
    sys.exit(0)
    
maxrefs = 100

fastaseqs = fasta.loadMfa(sys.argv[1])
expression = Expression(sys.argv[2])
ghmm = GHMM(sys.argv[3])
predictions = Features()
numrefs = min(maxrefs, len(fastaseqs))
a = 1

for (header, seq) in fastaseqs:
    print str(a) + '/' + str(numrefs) + ' - Predicting genes in ' + sys.argv[1]+', '+ header
    try:
        intensities = expression[header]
    except KeyError:
        print 'Cannot find ', header, ' in ', sys.argv[2]; sys.exit(0)
    if len(intensities) != len(seq):
        raise ValueError, 'Intensities and sequence are different lengths.'+ str(len(intensities))+','+ str(len(seq))
    parse = decoder(ghmm, SequenceDict(seq, 5), expression[header])
    parsetoFeatures(ghmm, parse, predictions, header)
    if a >= numrefs: break
    else: a += 1

predictions.writeGff(sys.argv[4])
