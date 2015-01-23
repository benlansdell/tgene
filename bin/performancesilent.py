#!/usr/bin/env python
"""
performancesilent.py <pred_prefix> <actual_sequence> <actual_gff> <output_file>

Program to rate performance of several gene predictions given set of actual coordinates
and predicted coordinates. Only compute nucleotide scores if fasta_file is provided.
All features in gff are loaded.

Fasta file header must be in format:
>name:start:end

Returns a list of genes correctly predicted, as well as sensivitiy, specificity
at nuc, exon and gene level and splice site stats.
"""

import sys
import getopt
import os
import fasta
import string
from fnmatch import fnmatch
from numpy import zeros
from sequence import *
from feature import *

def exonString(genes, startref, endref):
	"""Returns string of length length, with F for feature, N for not 
	"""
	string = (endref-startref+1)*'N'
	for gene in genes:
		for exon in genes[gene].coords:
			start = exon[0] - startref
			end = exon[1] - startref
			exonlength = end - start + 1
			exonstring = (exonlength)*'F'
			string = replaceString(string, exonstring, start)
			
	return string

def verifyFasta(head,seq,pred):
	"""Verifies if all splice sites, start/stops match consensus sequence. 
	"""
	return True

def computeNucStats(trues, preds):
	"""Returns vector with nuc stats. Namely TP, TN, FP, FN count.
	Overlaps two strings returned by ExonString to compute stats.
	Ignores strand.	Slow.
	"""
	#TP, TN, FP, FN
	count = zeros(4)
	for a in range(len(trues)):
		if (trues[a] == 'F') & (preds[a] == 'F'):
			count[0] += 1
		elif (trues[a] == 'N') & (preds[a] == 'N'):
			count[1] += 1
		elif (trues[a] == 'N') & (preds[a] == 'F'):
			count[2] += 1
		elif (trues[a] == 'F') & (preds[a] == 'N'):
			count[3] += 1
	return count

def classifyExonStats(actual, predicted, ref, classifiedexons):
	"""Classifies exons into types ME, WE, PE, CE"""
	actualexons = actual.exonList(ref)
	predictedexons = predicted.exonList(ref)
	exon = [{'Wrong':0, 'Partial':0, 'Correct':0}, {'Wrong':0, 'Partial':0, 'Correct':0}]
	
	for (astart, aend) in actualexons:
		category = 'Wrong'
		for (pstart, pend) in predictedexons:
			if (pstart, pend) == (astart, aend):
				category = 'Correct'
				break
			elif (pstart <= astart and pend >= astart) or (astart <= pstart and aend >= pend):
				category ='Partial'
				break
		exon[0][category] += 1
	for (pstart, pend) in predictedexons:
		category = 'Wrong'
		for (astart, aend) in actualexons:
			if (pstart, pend) == (astart, aend):
				category = 'Correct'
				break
			elif (astart <= pstart and aend >= pend) or (pstart <= astart and pend >= astart):
				category ='Partial'
				break
		exon[1][category] += 1
	classifiedexons['ME'] += exon[0]['Wrong']
	classifiedexons['WE'] += exon[1]['Wrong']
	classifiedexons['PE'] += (exon[1]['Partial'] + exon[0]['Partial'])/2
	classifiedexons['CE'] += exon[1]['Correct']
		

def computeSpliceStats(actual, predicted, ref, splicestats):
	actualcounts = {'ATG':[], 'GT':[], 'AG':[], 'TGA':[]}
	predictedcounts = {'ATG':[], 'GT':[], 'AG':[], 'TGA':[]}
	for gene in actual[ref]:
		for (start, end) in actual[ref][gene].coords:
			if start == actual[ref][gene].min:
				if actual[ref][gene].strand == '+':
					actualcounts['ATG'].append(start)
				else:
					actualcounts['TGA'].append(start)
			else:
				if actual[ref][gene].strand == '+':
					actualcounts['AG'].append(start)
				else:
					actualcounts['GT'].append(start)
			if end == actual[ref][gene].max:
				if actual[ref][gene].strand == '+':
					actualcounts['TGA'].append(end)
				else:
					actualcounts['ATG'].append(end)
			else:
				if actual[ref][gene].strand == '+':
					actualcounts['GT'].append(end)
				else:
					actualcounts['AG'].append(end)
	for gene in predicted[ref]:
		for (start, end) in predicted[ref][gene].coords:
			if start == predicted[ref][gene].min:
				if predicted[ref][gene].strand == '+':
					predictedcounts['ATG'].append(start)
				else:
					predictedcounts['TGA'].append(start)
			else:
				if predicted[ref][gene].strand == '+':
					predictedcounts['AG'].append(start)
				else:
					predictedcounts['GT'].append(start)
			if end == predicted[ref][gene].max:
				if predicted[ref][gene].strand == '+':
					predictedcounts['TGA'].append(end)
				else:
					predictedcounts['ATG'].append(end)
			else:
				if predicted[ref][gene].strand == '+':
					predictedcounts['GT'].append(end)
				else:
					predictedcounts['AG'].append(end)			

	for type in actualcounts:
		count = [0,0,0,0]
		for coord in actualcounts[type]:
			if coord in predictedcounts[type]:
				count[0] += 1
		count[3] += len(actualcounts[type]) - count[0]
		count[2] += len(predictedcounts[type]) - count[0]
		if not type in splicestats: splicestats[type] = count
		else: splicestats[type] = [splicestats[type][a] + count[a] for a in range(4)]

def removeDuplicates(list):
	retlist = []
	for item in list:
		if not item in retlist:
			retlist.append(item)
	return retlist
	
def computeExonStats(actual, predicted, ref):
	"""
	Returns vector with exon stats. Namely TP, TN, FP, FN count.
	"""
	#TP, TN, FP, FN
	count = zeros(4)
	#Note that duplicate exons are removed from this count, but not from the gene count...
	actualexon = removeDuplicates(actual.exonList(ref))
	predictedexon = predicted.exonList(ref)
	#Count true positives, the rest can be inferred from there
	for exon in actualexon:
		if predictedexon.count(exon) > 0:
			count[0] += 1
	count[3] = len(actualexon) - count[0]
	count[2] = len(predictedexon) - count[0]
	return count

def computeGeneStats(actual, predicted, ref):
	count = zeros(4)
	for gene in actual[ref]:
		if actual[ref][gene] in predicted[ref].values():
			#print 'Match:', gene
			count[0] += 1
	count[3] = len(actual[ref]) - count[0]
	count[2] = len(predicted[ref]) - count[0]
	return count

def printExonItem(dict,key):
	ME = int(dict[key]['ME'])
	WE = int(dict[key]['WE'])
	CE = int(dict[key]['CE'])
	PE = int(dict[key]['PE'])
	totalPredicted = float(WE + CE + PE)
	totalActual = float(ME + CE + PE)
	try:
		strng = key + '\t' + str(100*ME/totalActual) + '\t' + str(100*WE/totalPredicted) + '\t' + str(PE) + '\t' + str(CE) + '\n'
	except:
		strng = key + '\t' + '--' + '\t' + '--' + '\t' + str(PE) + '\t' + str(CE) + '\n'
	return strng

def printItem(dict, key):
	strng = key + '\t'
	l = list(dict[key])
	for item in l:
		strng += str(item) + '\t'
	strng += str(100*float(l[0])/(l[0] + l[3])) + '\t' + str(100*float(l[0])/(l[0] + l[2]))
	strng += '\n'
	return strng

def main(fastafile, predfile, actualfile):
	readfasta = True
	#Load each gff
	nucstats, exonstats, genestats = zeros(4), zeros(4), zeros(4)
	exonclassified = {'ME':0, 'WE':0, 'PE':0, 'CE':0}
	splicestats = {}
	predicted = Features(predfile, 'exon')
	actual = Features(actualfile)
	syncReferences(actual, predicted)
	#Nucleotide stats
	if readfasta == True:
		try:
			fastaseqs = fasta.loadMfa(fastafile)
		except StandardError:
			print "Failed to load fasta sequence."
			sys.exit(2)
		for (head,seq) in fastaseqs:
			header = head.split(':')
			#print header, head, head in actual
			if header[0] in actual: title = header[0]
			if head in actual: title = head
			startref = string.atoi(header[1])
			endref = startref + len(seq) - 1
			#verifyFasta(head,seq,predlist)
			if title in actual:
				actualnuc = exonString(actual[title], startref, endref)
				predictednuc = exonString(predicted[title], startref, endref)
				stats = computeNucStats(actualnuc, predictednuc)
				nucstats += stats
			else:
				print 'Sequence from ', title, ' has no actual gene annotation. Skipping.'
	
	#Exon, gene stats
	for ref in actual:
		computeSpliceStats(actual, predicted, ref, splicestats)
		classifyExonStats(actual, predicted, ref, exonclassified)
		exonstats += computeExonStats(actual, predicted, ref)
		genestats += computeGeneStats(actual, predicted, ref)
	return (nucstats, exonstats, genestats, exonclassified)

if __name__ == "__main__":
	#Read arguments
	if len(sys.argv) != 5:
		print __doc__
		sys.exit(0)
	prefix = sys.argv[1]
	fastafile = sys.argv[2]
	actualfile = sys.argv[3]
	outputfile = open(sys.argv[4], 'w')
	nucstats = {}
	exonstats = {}
	genestats = {}
	exonclassified = {}
	
	files = filter(lambda x: fnmatch(x, '*.gff') and x[0:len(prefix)] == prefix, os.listdir('.'))
	
	for file in files:
		(nuc, exon, gene, exonclass) = main(fastafile, file, actualfile)
		nucstats[file] = nuc
		exonstats[file] = exon
		genestats[file] = gene
		exonclassified[file] = exonclass
	outputfile.write("Nucleotide stats: TP TN FP FN SN SP\n")
	for file in nucstats:
		outputfile.write(printItem(nucstats, file))
	outputfile.write("\nExon stats: TP TN FP FN SN SP\n")
	for file in nucstats:
		outputfile.write(printItem(exonstats, file))
	outputfile.write('\nExon classification stats: ME% WE% PE CE\n')
	for file in nucstats:
		outputfile.write(printExonItem(exonclassified, file))
	outputfile.write("\nGene stats: TP TN FP FN SN SP\n")
	for file in nucstats:
		outputfile.write(printItem(genestats, file))
	
