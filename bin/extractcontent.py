#!/usr/bin/env python
"""extractcontent.py <fasta_file> <gff_file> <output_prefix>

Extracts exon, intron and intergenic sequence (and positions) from provided fasta and gff files. Fasta file must be in format:
>ref:start:end
ATGCTACG...

Features in gff file must have same ref as gff file.

It is recommended that util/verifygff.py is used to check that fasta file and gff file match and give valid gene structures.

Outputs a number of files, named with .einit, .eterm, etc. appeneded to 'output_prefix'
"""

import sys, fasta
from feature import *
from sequence import reverseComplement

def swap(frame):
    if frame == 0: return 2
    elif frame == 1: return 1
    return 0

if __name__ == "__main__":
    
    if len(sys.argv) != 4:
        print 'Invalid arguments'
        print __doc__
        sys.exit(0)
    
    genes = Features(sys.argv[2])
    fastaseqs = fasta.loadMfa(sys.argv[1])
    outfile = sys.argv[3]
    outgff = outfile + '.gff'
    outfa = outfile + '.fa'
        
    outfiles = {}
    outfiles['einit'] = open(outgff+'.einit', 'w'), fasta.MfaWriter(outfa+'.einit')
    outfiles['eintn'] = open(outgff+'.eintn', 'w'), fasta.MfaWriter(outfa+'.eintn')
    outfiles['esing'] = open(outgff+'.esing', 'w'), fasta.MfaWriter(outfa+'.esing')
    outfiles['eterm'] = open(outgff+'.eterm', 'w'), fasta.MfaWriter(outfa+'.eterm')
    outfiles['intrn'] = open(outgff+'.intrn', 'w'), fasta.MfaWriter(outfa+'.intrn')
    outfiles['inter'] = open(outgff+'.inter', 'w'), fasta.MfaWriter(outfa+'.inter')
    
    
    for (head,seq) in fastaseqs:
        head = head.split(":")
        ref = head[0]
        startref = int(head[1])
        endref = startref + len(seq) - 1
        split = Features(None, {})
        if ref in genes:
            prev = startref-1
            for generef in genes[ref]:
                length = 0
                gene = genes[ref][generef]
                if gene.strand == '+': forward, backward = 'einit','eterm'
                else: forward, backward = 'eterm','einit'
                coords = gene.coords
                for coord in coords:
                    frame = str(length%3)
                    length += coord[1] - coord[0] + 1
                    if gene.strand == '-': frame = str(swap((length-1)%3))
                    if (coord[0] == gene.min) and (coord[1] == gene.max):
                        split.addFeature(Feature(gene.ref+'_inter', gene.source, ['inter'], [[prev+1, coord[0]-1]], gene.score, '+', '.', gene.name))
                        split.addFeature(Feature(gene.ref+'_esing', gene.source, ['esing'], [coord], gene.score, gene.strand, frame, gene.name))
                    elif coord[0] == gene.min:
                        split.addFeature(Feature(gene.ref+'_inter', gene.source, ['inter'], [[prev+1, coord[0]-1]], gene.score, '+', '.', gene.name))
                        split.addFeature(Feature(gene.ref+'_'+forward, gene.source, [forward], [coord], gene.score, gene.strand, frame, gene.name))
                    elif coord[1] == gene.max:
                        split.addFeature(Feature(gene.ref+'_intrn', gene.source, ['intrn'], [[prev+1, coord[0]-1]], gene.score, gene.strand, '.', gene.name))
                        split.addFeature(Feature(gene.ref+'_'+backward, gene.source, [backward], [coord], gene.score, gene.strand, frame, gene.name))
                    else:
                        split.addFeature(Feature(gene.ref+'_intrn', gene.source, ['intrn'], [[prev+1, coord[0]-1]], gene.score, gene.strand, '.', gene.name))
                        split.addFeature(Feature(gene.ref+'_eintn', gene.source, ['eintn'], [coord], gene.score, gene.strand, frame, gene.name))
                    prev = coord[1]
            split.addFeature(Feature(gene.ref+'_inter', gene.source, ['inter'], [[prev+1, endref]], gene.score, '+', '.', gene.name))
            for typeref in split:
                type = typeref.split('_')[-1]
                if outfiles.has_key(type):
                    split.writeGff(outfiles[type][0], typeref)
                    for generef in split[typeref]:
                        gene = split[typeref][generef]
                        for a in range(len(gene.coords)):
                            start = gene.coords[a][0]
                            end = gene.coords[a][1]
                            ref = gene.ref
                            strand = gene.strand
                            outhead = gene.ref + ':' + str(start) + ':' + str(end) + ':' + strand + ':' + gene.frame[a]
                            outseq = seq[start-startref:end-startref+1]
                            if strand == '-': outseq = reverseComplement(outseq)
                            outfiles[gene.type[0]][1].write(outhead, outseq)