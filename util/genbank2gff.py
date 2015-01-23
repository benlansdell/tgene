#!/usr/bin/env python
"""genbank2gff.py <genbank file>
Script to output CDS from a genbank file"""

import sys, fasta, string
from common.feature import *

def parseCoordsLine(string):
	#join(2613..3702,3779..6055,6119..6816)
	#join(2613..3702,3779..6055,6119..6816,
	#19077..19284,24304..24446,25225..25374,29379..29542,
	#29661..29824,32733..32845)
	#29661..29824
	#join(8413..8598,13951..14262,18008..18074,22016..22264,2531
	#print string
	if 'join' in string:
		string = string[5:]
	if string[-1] in [')']:
		string = string[:-1]
	return string

def splitCoords(line):
#	print line
	words = line.split(',')
	words = [word.split('..') for word in words]
#	print words
	coords = [[int(word[0]), int(word[1])] for word in words]
	return coords

def parseSeqLine(line):
	#1 acatccacgg caacatcgac ggagtcttcg agtggatctc ccccgagggt gtccatgtgc
	seqline = ''
	for a in line:
		if a in string.lowercase:
			seqline += a
	return string.upper(seqline)

if len(sys.argv) != 2:
	print 'Invalid arguments'
	print __doc__
	sys.exit(0)
	
file = open(sys.argv[1])
features = Features()
fastawriter = fasta.MfaWriter(sys.argv[1] + '.fa')
reference = "augustus_train"
genenumber = 0
mode = None

for line in file:
	if line:
		#print line
		words = line.split()
		if words:
			if words[0][0] == '/':
				if mode == 'FASTA':
					head = reference + str(genenumber) + ':1:' + str(len(seq)-1)
					fastawriter.write(head, seq)
				if mode == 'CDS':
					intcoords = splitCoords(coords)
					#print intcoords
					features.addFeature(Feature(reference + str(genenumber), reference, ['CDS']*len(intcoords), intcoords, [None]*len(intcoords), '+', [None]*len(intcoords), str(genenumber)))
				mode = None
			#print line
			if mode == 'CDS':
				coords += parseCoordsLine(words[0])
			if mode == 'FASTA':
				seq += parseSeqLine(line)
			if words[0] == 'CDS' and len(words) == 2:
				mode = 'CDS'
				genenumber += 1
				coords = parseCoordsLine(words[1])
			if words[0] == 'ORIGIN':
				mode = 'FASTA'
				seq = ''
				
features.writeGff(sys.argv[1] + '.gff')