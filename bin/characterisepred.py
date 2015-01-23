#!/usr/bin/env python
"""characterisepred.py <predicted_gff> <startref> <endref>

Program to rate performance of gene predictions given set of actual coordinates
and predicted coordinates.

NOT WORKING"""

import sys
from feature import *

def genestring(fullgenes, exongenes, expressednames, startref, endref):
	string = (endref-startref+1)*'N'
	expressedutrcount = 0
	lengths = []
	for ref in fullgenes:
		for generef in fullgenes[ref]:
			gene = fullgenes[ref][generef]
			start = gene.min - startref
			end = gene.max - startref
			exonlength = end - start + 1
			exonstring = (exonlength)*'U'
			if generef in expressednames and not exonlength in lengths:
				expressedutrcount += exonlength
				lengths.append(exonlength)
				#print generef, '+', exonlength
			string = replaceString(string, exonstring, start)	
	lengths = []
	for ref in exongenes:
		for generef in exongenes[ref]:
			gene = exongenes[ref][generef]
			start = gene.min - startref
			end = gene.max - startref
			exonlength = end - start + 1
			exonstring = (exonlength)*'I'
			if generef in expressednames and not exonlength in lengths:
				expressedutrcount -= exonlength
				lengths.append(exonlength)
				#print generef, '-', exonlength
			string = replaceString(string, exonstring, start)			
			for [exonmin, exonmax] in gene.coords:
				start = exonmin - startref
				end = exonmax - startref
				exonlength = end - start + 1
				exonstring = (exonlength)*'E'
				string = replaceString(string, exonstring, start)
	return (string, expressedutrcount)


def exonstring(genes, startref, endref):
	string = (endref-startref+1)*'N'
	for ref in genes:
		for generef in genes[ref]:
			gene = genes[ref][generef]
			for exon in gene.coords:
				start = exon[0] - startref
				end = exon[1] - startref
				exonlength = end - start + 1
				exonstring = (exonlength)*'F'
				string = replaceString(string, exonstring, start)			
	return string

if len(sys.argv) != 4:
	print 'Invalid arguments'
	print __doc__
	sys.exit(0)

annotationfile = '/home/users/lab0605/lansdell/TileGene/data/annotation/adh.flybase.2004.gff'
expressedfile = '/home/users/lab0605/lansdell/TileGene/data/annotation/adh.flybase.2004.combined.gff.expressed'

fullgenes = Features(annotationfile)
exongenes = Features(annotationfile, 'CDS')
expressed = Features(expressedfile)
predictions = Features(sys.argv[1])

startref = int(sys.argv[2])
endref = int(sys.argv[3])
expressednames = []
for ref in expressed:
	for gene in expressed[ref]:
		expressednames.append(gene)
(annotatedstring, expressedutrcount) = genestring(fullgenes, exongenes, expressednames, startref, endref)
predictedstring = exonstring(predictions, startref, endref)

nexon = [annotatedstring[a] == 'E' and predictedstring[a] == 'F' for a in range(len(predictedstring))].count(True)
ninter = [annotatedstring[a] == 'N' and predictedstring[a] == 'F' for a in range(len(predictedstring))].count(True)
nintron = [annotatedstring[a] == 'I' and predictedstring[a] == 'F' for a in range(len(predictedstring))].count(True)
nutr = [annotatedstring[a] == 'U' and predictedstring[a] == 'F' for a in range(len(predictedstring))].count(True)
ntotal = [predictedstring[a] == 'F' for a in range(len(predictedstring))].count(True)

print 'Using:', annotationfile
print '% in exons:', float(nexon)*100/ntotal
print '% in intergenic:', float(ninter)*100/ntotal
print '% in introns:', float(nintron)*100/ntotal
print '% in UTR:', float(nutr)*100/ntotal
print '% of expressed UTR correctly predicted (sensitivity):', 100*float(nutr)/expressedutrcount
print '% of full UTR correctly predicted (sensitivity):', 100*float(nutr)/annotatedstring.count('U')
print 'Total:', ntotal
print 'Expressed UTR', expressedutrcount
print 'Total UTR', annotatedstring.count('U')