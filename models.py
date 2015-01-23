"""models.py
Probabilistic models used in GHMM"""

from io import *
from sequence import *

def parseProbList(problist):
    if len(problist) == 5: prefix = problist[0]
    else: prefix = ''
    freqlist = [int(a) for a in problist[-4:]]
    return (prefix, freqlist)

def readProbTable(file, pseudocountused):
    if pseudocountused: pseudocount = 1
    else: pseudocount = 0
    probtable = {0:{}, 1:{}, 2:{}, 3:{}, 4:{}, 5:{}}
    for line in open(file):
        (prefix, freqlist) = parseProbList(line.split())
        order = len(prefix)
        freqsum = max(sum(freqlist) + 4*pseudocount, 1)
        for a in range(len(freqlist)):
            seq = higherOrder(seqToInt(prefix)+[a])
            probtable[order][seq] = log(float(freqlist[a] + pseudocount)/freqsum)
    return probtable

class HMC:
    """Homogeneous markov chain"""
    def __init__(self, files, pseudo, order):
        self.pseudocount = pseudo
        self.order = order
        filelist = files.split()[0].split(',')
        self.nfiles = len(filelist)
        if self.nfiles != 1:
            raise ValueError, 'HMC5 requires 1 file: ' + files
        self.prob = readProbTable(filelist[0], self.pseudocount)
    
    def probEmit(self, sequence, start, end, phase, strand):
        prob = 0
        if strand == '-':
            start, end = len(sequence)-end, len(sequence)-start
        start = max(0, start)
        end = min(end, len(sequence))
        if end-start == 1:
            order = min(self.order, start)
            prob = self.prob[order][sequence[strand][order][start]]
        else:
            for a in range(start, end):
                order = min(self.order, a-start)
                prob += self.prob[order][sequence[strand][order][a]]
        return prob
    
class P3IMC:
    """3 period inhomogeneous markov chain"""
    def __init__(self, files, pseudo, order):
        self.pseudocount = pseudo
        self.order = order
        self.prob = {}
        filelist = files.split()[0].split(',')
        self.nfiles = len(filelist)
        if self.nfiles != 3:
            raise ValueError, 'P3IMC5 requires 3 files: ' + files
        for a in range(3):
            self.prob[a] = readProbTable(filelist[a], self.pseudocount)
    
    def probEmit(self, sequence, start, end, phase, strand):
        prob = 0
        if strand == '-':
            start, end = len(sequence)-end, len(sequence)-start
        start = max(start,0)
        end = min(end, len(sequence))
        if end-start == 1:
            order = min(self.order, start)
            prob = self.prob[phase%3][order][sequence[strand][order][start]]
        else:
            for a in range(start, end):
                order = min(self.order, a-start)
                prob += self.prob[phase%3][order][sequence[strand][order][a]]
                phase += 1
        return prob
    
class IMC:
    """Inhomogeneous markov chain"""
    def __init__(self, files, pseudo, order):
        self.pseudocount = pseudo
        self.order = order
        self.prob = {}
        filelist = files.split()[0].split(',')
        self.nfiles = len(filelist)
        for a in range(self.nfiles):
            self.prob[a] = readProbTable(filelist[a], self.pseudocount)
    
    def probEmit(self, sequence, start, end, phase, strand):
        prob = 0
        if strand == '-':
            start, end = len(sequence)-end, len(sequence)-start
        start = max(start, 0)
        end = min(end, len(sequence))
        if end-start > self.nfiles:
            raise ValueError, 'Sequence longer than nfiles, ' + str(nfiles)
        for a in range(start, end):
#            order = min(self.order, a-start)
            order = min(self.order, start)
            prob += self.prob[a-start][order][sequence[strand][order][a]]
        return prob

class Null:
    def __init__(self, files, pseudo, order):
        pass
    def probEmit(self, sequence, start, end, phase, strand):
        return 0.0
    
modeldict = {'IMC':IMC, 'P3IMC':P3IMC, 'HMC':HMC, 'NULL':Null}