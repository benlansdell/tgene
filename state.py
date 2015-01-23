"""state.py
Holds class and functions for states in GHMM framework"""

from lengthdist import *
from models import *
from tiling import *

class State:
    def __init__(self, id, config, ghmm):
        self.id = id
        #The weight assigned to tiling probability scores...
        self.weight = ghmm.weight
        self.newweight = ghmm.newweight
        self.symbol = config['symbol']
        self.coding = config['coding']
        self.phase = config['phase']
        self.strand = config['strand']        
        self.content = ghmm.content[config['content']]
        self.tiling = ghmm.tilings[config['tiling']]
        if config['signals']:
            signalrefs = config['signals'].split(',')
            self.frontsignal = ghmm.signals[signalrefs[0]]
            self.backsignal = ghmm.signals[signalrefs[1]]
            self.hassignal = True
        else:
            self.frontsignal = None
            self.backsignal = None
            self.hassignal = False
        distribpara = config['lengthdist'].split(',',1)
        self.duration = lengthdict[distribpara[0]](distribpara[1])
        
    def signalOffset(self, start, end, swap):
        if self.hassignal:
            #if self.strand == '-' and swap: front = self.backsignal; back = self.frontsignal
            front = self.frontsignal; back = self.backsignal
            if front.coninclude: startoffset = start + front.position
            else: startoffset = start + front.position + front.consensuslength
            if back.coninclude: endoffset = end - back.length + back.position + back.consensuslength
            else: endoffset = end - back.length + back.position
            if self.strand == '-' and swap: startoffset, endoffset = start + end-endoffset, start + end-startoffset    
            return (startoffset, endoffset)
        else:
            return (start, end)
            
    def probDuration(self, l):
        return self.duration.probLength(l)
    
    def probCumDuration(self):
        return self.duration.probCumLength()
        
    def probEmit(self, sequence, start, end, phase, expression):
        if phase == None: inphase = self.phase
        else: inphase = phase
        tileprob = self.tiling.probEmit(expression, start, end)
        if self.hassignal:
            if self.strand == '+':
                front = self.frontsignal; back = self.backsignal
            else:
                front = self.backsignal; back = self.frontsignal
            contentstart = start + front.length
            contentend = end - back.length
            if not self.newweight: return (front.probEmit(sequence, start, self.strand) + self.content.probEmit(sequence, contentstart, contentend, inphase, self.strand) + back.probEmit(sequence, contentend, self.strand))*(1-self.weight) + self.weight*tileprob
            else: return (front.probEmit(sequence, start, self.strand) + self.content.probEmit(sequence, contentstart, contentend, inphase, self.strand) + back.probEmit(sequence, contentend, self.strand)) + self.weight*tileprob
        else:
            if not self.newweight: return (1-self.weight)*self.content.probEmit(sequence, start, end, inphase, self.strand) + self.weight*tileprob
            else: return self.content.probEmit(sequence, start, end, inphase, self.strand) + self.weight*tileprob
            
class Signal:
    def __init__(self, id, config):
        self.id = id
        self.position = config['position']
        self.coninclude = config['coninclude']
        self.length = config['length']
        self.consensus = [seqToInt(a) for a in config['consensus'].split(',')]
        self.consensuslength = len(self.consensus[0])
        self.consensusposition = self.position + self.consensuslength
        self.model = modeldict[config['model']](config['infiles'], config['pseudo'], config['order'])
        
    def checkSignal(self, sequence, strand, pos):
        if sequence[strand][0][pos+self.position:pos+self.position+self.consensuslength] in self.consensus:
            if pos + self.length <= len(sequence):
                return True
        return False
    
    def probEmit(self, sequence, pos, strand):
        return self.model.probEmit(sequence, pos, pos+self.length, 0, strand)
      
class Content:
    def __init__(self, id, config):
        self.id = id
        self.model = modeldict[config['model']](config['infiles'], config['pseudo'], config['order'])
    
    def probEmit(self, sequence, start, end, phase, strand):
        return self.model.probEmit(sequence, start, end, phase, strand)
    
class Tiling:
    def __init__(self, id, config):
        self.id = id
        self.model = tilingmodeldict[config['model']](config['infiles'])
    
    def probEmit(self, expression, start, end):
        return self.model.probEmit(expression, start, end)
        