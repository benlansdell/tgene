#!/usr/bin/env python
"""gffformat.py <gff_file>

Script to format gff training files

Filters for CDS predictions in Adh region

Outputs file to gff_file.out
"""
from common.feature import *
import sys

exons = Features(sys.argv[1], ['CDS'])
newfeatures = Features()
for ref in exons:
	name = ref.split(':')[0]
	#start = int(ref.split(":")[1])
	for generef in exons[ref]:
		gene = exons[ref][generef]
		gene.ref = name
		newfeatures.addFeature(gene)
		#if gene.strand == '+':
		#	gene.coords[-1][1] += 3
		#if gene.strand == '-':
		#	gene.coords[0][0] -= 3
		#for a in range(len(gene.coords)):
		#	gene.coords[a][0] += start-1
		#	gene.coords[a][1] += start-1
newfeatures.writeGff(sys.argv[1] + '.out')
		
