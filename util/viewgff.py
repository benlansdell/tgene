#!/usr/bin/env python
"""viewgff.py <fasta_file> <gff_file>

Outputs sequences from a list of gff features. ref from fasta file must match gff features.
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
	print startref
	if head in genes: ref = head; print ref
	else: ref = header[0]
	if ref in genes:
		for name in genes[ref]:
			gene = genes[ref][name]
			print '**', name, '**\n'
			codingsequence = ''
			coords= gene.coords[:]
			if gene.strand == '-':
				coords.reverse()
			for coord in coords:
				print gene.ref + ',', gene.name, gene.strand, coord, ':'
				if gene.strand == '-':
					codingsequence = seq[coord[0]-startref-flank:coord[1]+1-startref+flank] + codingsequence
					print reverseComplement(seq[coord[0]-startref-flank:coord[1]+1-startref+flank])
				if gene.strand == '+':
					codingsequence += seq[coord[0]-startref-flank:coord[1]+1-startref+flank]
					print seq[coord[0]-startref-flank:coord[1]+1-startref+flank]
			if gene.strand == '-': codingsequence = reverseComplement(codingsequence)
			try: print '\nEncodes: ', translate(codingsequence), '\n'
			except ValueError, msg: print msg

