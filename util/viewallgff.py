#!/usr/bin/env python
"""viewallgff.py <fasta_file> <gff_file>

Script to output sequences, combined into one gene, from a list of gff features
"""

import fasta, sys
from common.feature import *
from common.sequence import *

if len(sys.argv) != 3:
	print __doc__
	sys.exit(0)

genes = Features(sys.argv[2])
fastaseqs = fasta.loadMfa(sys.argv[1])
flank = 0

for (head,seq) in fastaseqs:
	header = head.split(':')
	startref = int(header[1])
	if head in genes: ref = head; print ref
	else: ref = header[0]
	if ref in genes:
		for name in genes[ref]:
			gene = genes[ref][name]
			print '>' + name
			codingsequence = ''
			coords= gene.coords[:]
			if gene.strand == '-':
				coords.reverse()
			for coord in coords:
				if gene.strand == '-':
					codingsequence = seq[coord[0]-startref-flank:coord[1]+1-startref+flank] + codingsequence
				if gene.strand == '+':
					codingsequence += seq[coord[0]-startref-flank:coord[1]+1-startref+flank]
			if gene.strand == '-': codingsequence = reverseComplement(codingsequence)
			print codingsequence
