#!/usr/bin/env python
"""extractsensor.py <fasta_file> <gff_file> <output_prefix>

Extracts start/stop codon and splice site sequences from provided fasta and gff files. Fasta file must be in format:
>ref:start:end
ATGCTACG...

Features in gff file must have same ref as gff file.

It is recommended that util/verifygff.py is used to check that fasta file and gff file match and give valid gene structures.

Outputs a number of files, named with .AGs, .GTs, etc. appeneded to 'output_prefix'
"""

import sys, fasta
from feature import *
from sequence import reverseComplement

def writeSeq(ref, seq, pos, type, strand):
    preoffset, postoffset = offsets[type]
    ntype = type
    if strand == '-':
        if ntype == 'ATG': type = 'TGA'
        elif ntype == 'TGA': type = 'ATG'
        elif ntype == 'GT': type = 'AG'
        else: type = 'GT'
        postoffset += 2
    else:
        preoffset += 2
    retseq = seq[pos-preoffset:pos+postoffset]
    if strand == '-': retseq = reverseComplement(retseq)
    head = ref + ':' + str(pos-preoffset) + ':' + str(pos+postoffset-1) + ':' + strand
    outfiles[type].write(head, retseq) 

if len(sys.argv) != 4:
    print 'Invalid arguments'
    print __doc__
    sys.exit(0)
    
fastaseqs = fasta.loadMfa(sys.argv[1])
features = Features(sys.argv[2])

outprefix = sys.argv[3]
outfiles = {'AG': fasta.MfaWriter(outprefix + '.AGs'), 'GT': fasta.MfaWriter(outprefix + '.GTs'),
               'ATG': fasta.MfaWriter(outprefix + '.ATGs'), 'TGA': fasta.MfaWriter(outprefix + '.TGAs')}
plusATG = 33
minusATG = 3
plusTGA = 12
minusTGA = 3
plusGT = 12
minusGT = 3
plusAG = 12
minusAG = 3
offsets = {'ATG': [plusATG,3+minusATG], 'GT': [minusGT-1,plusGT+3], 'AG': [plusAG+2,minusAG], 'TGA': [minusTGA+2,plusTGA+1]}    
#offsets = {'ATG': [6,9], 'GT': [5,9], 'AG': [8,6], 'TGA': [8,7]}    
#offsets = {'ATG': [9,12], 'GT': [8,12], 'AG': [11,9], 'TGA': [11,10]}    
   
for (head, seq) in fastaseqs:
    head = head.split(':')
    ref = head[0]
    startref = int(head[1])
    if ref in features:
        for generef in features[ref]:
            gene = features[ref][generef]
            for (start, end) in gene.coords:
                if (start == gene.min) and (end == gene.max):
                    writeSeq(ref, seq, start-startref, 'ATG', gene.strand)
                    writeSeq(ref, seq, end-startref, 'TGA', gene.strand)
                elif (start == gene.min):
                    writeSeq(ref, seq, start-startref, 'ATG', gene.strand)
                    writeSeq(ref, seq, end-startref, 'GT', gene.strand)
                elif (end == gene.max):
                    writeSeq(ref, seq, start-startref, 'AG', gene.strand)
                    writeSeq(ref, seq, end-startref, 'TGA', gene.strand)
                else:
                    writeSeq(ref, seq, start-startref, 'AG', gene.strand)
                    writeSeq(ref, seq, end-startref, 'GT', gene.strand)