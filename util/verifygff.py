#!/usr/bin/env python
"""viewgff.py <fasta_file> <gff_file>

Script to check genes contain no inframe stop codons and encode a protein.
Removes ones that aren't valid, output to fasta_file.out
"""

import fasta, sys
from common.feature import *
from common.sequence import *

if len(sys.argv) != 3:
	print __doc__
	sys.exit(0)

genes = Features(sys.argv[2], 'CDS')
fastaseqs = fasta.loadMfa(sys.argv[1])
flank = 0
toremove = []

for (head,seq) in fastaseqs:
	head = head.split(':')
	ref = head[0]
	startref = int(head[1])
	if ref in genes:
		toremove = []
		for name in genes[ref]:
			gene = genes[ref][name]
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
			try: translation = translate(codingsequence)
			except:
				translation = 'M'
				print "Removing ", name, " because length is no a multiple of 3"
				if not name in toremove: toremove.append(name)
			if translation[0] != 'M':
				print "Removing ", name, " because doesn't start with M aa"
				if not name in toremove: toremove.append(name)
			elif '*' in translation[:-1]:
				print "Removing ", name, " because contains an inframe stop codon"
				if not name in toremove: toremove.append(name)
			elif translation[-1] != '*':
				print "Removing ", name, " because doesn't end with a stop codon"
				if not name in toremove: toremove.append(name)
		for name in toremove:
			genes[ref].pop(name)
genes.writeGff(sys.argv[2] + '.out')