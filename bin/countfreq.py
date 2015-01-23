#!/usr/bin/env python
"""countfreq.py <fasta_file> <output_file> <order>

For MLE training parameters - counts frequency of specified order emissions in provided sequence.
Fasta file must be in format:
>ref:start:end
ATCGATCAG...
"""

import fasta, sys, string
from sequence import *

def addCount(seq, counter):
    if len(seq) == 0:
        return
    if not seq[-1] in counter:
        counter[seq[-1]] = {}
    counter = counter[seq[-1]]
    if -1 in counter:
        counter[-1] += 1
    else:
        counter[-1] = 1    
    addCount(seq[:-1], counter)
    
def listCounts(counter, builtstring, builtdict):
    for key in counter:
        if key == -1:
            builtdict[builtstring] = counter[key]
        else:
            listCounts(counter[key], key+builtstring, builtdict)

def orderListRecurse(builtlist, templist, appendlist, order):
    builtlist += templist
    if order == 0:
        return builtlist
    else:
        newlist = [a+b for a in templist for b in appendlist]
        return orderListRecurse(builtlist, newlist, appendlist, order -1)
    
def orderList(builtlist, order):
    return orderListRecurse([], builtlist[:], builtlist[:], order)
    
def writeCounts(counter, filename):
    outfile = open(filename, 'w')
    keylist = orderList(['A', 'C', 'G', 'T'], order)
    #outfile.write('Frequency counts for ' + sys.argv[1] + ' order ' + str(order) + '\n\tA\tC\tG\tT\n')
    for key in keylist:
        if key in counter:
            val = str(counter[key])
        else:
            val = '0'
        first, last = key[:-1], key[-1]
        if last == 'T': outline = val + '\n' 
        elif last == 'A': outline = first + '\t' + val + '\t'
        else: outline = val + '\t'
        outfile.write(outline)
        
if len(sys.argv) != 4:
    print 'Not valid arguments'
    print __doc__
    sys.exit(0)
    
fastaseqs = fasta.loadMfa(sys.argv[1])
order = string.atoi(sys.argv[3])
countdict = {}
compactcount = {}

print 'Counting frequencies in', sys.argv[1]

for (head,seq) in fastaseqs:    
    for a in range(len(seq)):
        addCount(seq[max(0,a-order):a+1], countdict)
listCounts(countdict, '', compactcount)    
writeCounts(compactcount, sys.argv[2])