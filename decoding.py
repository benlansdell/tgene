"""decode.py
Decoding algorithms and related functions"""

from GHMM import *
from feature import *
import sys

def stdOutput(output):
    sys.stdout.write(output)
    sys.stdout.flush()

class ViterbiStructure:
    def __init__(self, ghmm, sequence, expression):
        self.length = len(sequence)
        self.ghmm = ghmm
        self.sequence = sequence
        self.expression = expression
        self.bestprevprob = -inf*ones((self.length+1, ghmm.nstates))
        self.bestprevncstate = zeros((self.length+1, ghmm.nstates))
        self.bestprevcstate = zeros((self.length+1, ghmm.nstates))
        self.bestprevpos = zeros((self.length+1, ghmm.nstates))
        for ncidx in ghmm.ncstates:
            ncstate = ghmm.states[ncidx]
            self.bestprevprob[0,ncidx] = ghmm.initial[ncidx] + ncstate.probEmit(sequence,0,1, None, expression) + ncstate.probCumDuration()
        
    def updateStructure(self, codingstates, backpos, postidx):
        changestateprob = [-inf]
        poststate = self.ghmm.states[postidx]
        for (preidx,idx,frontpos) in codingstates[backpos][postidx]:
            if preidx != idx:
                prestate = self.ghmm.states[preidx]
                state = self.ghmm.states[idx]
                if state.strand == '+': phase = prestate.phase
                else: phase = poststate.phase
                yemit = state.probEmit(self.sequence, frontpos, backpos, phase, self.expression)
                ydur = state.probDuration(backpos-frontpos)
                changestateprob.append([self.bestprevprob[frontpos-1, preidx] + prestate.probDuration(0) + self.ghmm.transition[preidx,idx] + ydur + yemit + self.ghmm.transition[idx,postidx] + poststate.probCumDuration() + poststate.probEmit(self.sequence, backpos, backpos+1, None, self.expression), preidx, idx, frontpos-1])                
                #if yemit != -inf and self.bestprevprob[frontpos-1,preidx] != -inf: print 'Start:', frontpos, 'End:', backpos, '\t(x,y,z)', prestate.symbol, state.symbol, poststate.symbol, '\tDelta:', self.bestprevprob[frontpos-1, preidx], '\t\tY emit:', yemit
            else:
                prestate = self.ghmm.states[preidx]
                changestateprob.append([self.bestprevprob[frontpos-1, preidx] + prestate.probDuration(0) + self.ghmm.transition[preidx,postidx] + poststate.probCumDuration() + poststate.probEmit(self.sequence, backpos, backpos+1, None, self.expression), preidx, preidx, frontpos-1])
        keepstateprob = [self.bestprevprob[backpos-1, postidx] + poststate.probCumDuration() + poststate.probEmit(self.sequence, backpos,backpos+1, None, self.expression), postidx, postidx, backpos-1]
        self.bestprevprob[backpos,postidx],self.bestprevncstate[backpos,postidx],self.bestprevcstate[backpos,postidx],self.bestprevpos[backpos,postidx] = max(keepstateprob, max(changestateprob))
        
    def terminateStructure(self, codingstates):
        for postidx in self.ghmm.ncstates:
            poststate = self.ghmm.states[postidx]
            changestateprob = [-inf]
            for (preidx, idx, frontpos) in codingstates[self.length][postidx]:
                if preidx != idx:
                    prestate = self.ghmm.states[preidx]
                    state = self.ghmm.states[idx]
                    if state.strand == '+': phase = prestate.phase
                    else: phase = poststate.phase
                    yemit = state.probEmit(self.sequence, frontpos, self.length, phase, self.expression)
                    ydur = state.probDuration(self.length-frontpos)
                    changestateprob.append([self.bestprevprob[frontpos-1, preidx] + prestate.probDuration(0) + self.ghmm.transition[preidx,idx] + ydur + yemit + self.ghmm.transition[idx,postidx], postidx, idx, frontpos-1])                
                else:
                    prestate = self.ghmm.states[preidx]
                    changestateprob.append([self.bestprevprob[frontpos-1, preidx] + prestate.probDuration(0) + self.ghmm.transition[preidx,postidx], preidx, preidx, frontpos-1])            
            keepstateprob = [self.bestprevprob[self.length-1, postidx], postidx, postidx, self.length-1]
            self.bestprevprob[self.length,postidx],self.bestprevncstate[self.length,postidx],self.bestprevcstate[self.length,postidx],self.bestprevpos[self.length,postidx] = max(keepstateprob, max(changestateprob))
            
    def generateParse(self):
        parse = [-1]*self.length
        oldpos = self.length
        pos = max([[self.bestprevprob[self.length,ncidx], self.bestprevpos[self.length,ncidx]] for ncidx in self.ghmm.ncstates])[1]
        ncstate = max([[self.bestprevprob[self.length,ncidx], self.bestprevncstate[self.length,ncidx]] for ncidx in self.ghmm.ncstates])[1]
        cstate = max([[self.bestprevprob[self.length,ncidx], self.bestprevcstate[self.length,ncidx]] for ncidx in self.ghmm.ncstates])[1]
        while pos != 0:
            parse[int(pos)+1:int(oldpos)] = [int(cstate)]*int(oldpos-1-pos)
            parse[int(pos)] = int(ncstate)
            oldpos = pos
            pos = self.bestprevpos[oldpos, ncstate]
            cstate = self.bestprevcstate[oldpos, ncstate]
            ncstate = self.bestprevncstate[oldpos, ncstate]
        parse[1:int(oldpos)] = [cstate]*int(oldpos-1)
        parse[0] = int(ncstate)
        return parse
    
def generateSignalList(ghmm, sequence):
    signals = {'+':[], '-':[]}
    stopsignalframe = {'+': {0:[], 1:[], 2:[]}, '-':{0:[], 1:[], 2:[]}}
    for strand in sequence:
        for pos in range(len(sequence)):
            for signal in ghmm.signals:
                if ghmm.signals[signal].checkSignal(sequence, strand, pos):
                    if signal == ghmm.stopcodon.id:
                        signalpos = pos + ghmm.stopcodon.position
                        stopsignalframe[strand][signalpos%3].append(signalpos)
                    signals[strand].append([pos, signal])
    #print signals
    return (signals, stopsignalframe)

def updateStopList(ghmm, sequence, strand, stopsignalframe, pos, signal):
    #Remove all stop signals behind current signal position
    maximumstop = 0
    for frame in stopsignalframe[strand]:
        for stoppos in stopsignalframe[strand][frame][:]:
            if stoppos < pos + ghmm.signals[signal].consensusposition:
                stopsignalframe[strand][frame].pop(0)
            else:
                break
        if not stopsignalframe[strand][frame]:
            maximumstop = len(sequence)
        else:
            maximumstop = max(maximumstop, stopsignalframe[strand][frame][0])
    return maximumstop

def updateBackSignals(ghmm, sequence, backsignals, frontpos):
    #Remove all back signals behind current signal position
    for [backpos, backsignal] in backsignals[:]:
        if backpos <= frontpos:
            backsignals.pop(0)
        else:
            break
        
def generateCNCTuples(ghmm, sequence, strand, frontbacksignals):
    ctuples = []
    cnctuples = []
    for [frontpos, frontsignal, backpos, backsignal] in frontbacksignals:
        for state in ghmm.cstates:
            if (ghmm.states[state].strand == strand) and (ghmm.states[state].frontsignal.id == frontsignal) and (ghmm.states[state].backsignal.id == backsignal):
                ctuples.append([state, frontpos, backpos+ghmm.signals[backsignal].length])
    for [state, frontpos, backpos] in ctuples:
        for prestate in ghmm.ncstates:
            for poststate in ghmm.ncstates:
                if ghmm.sparsity[prestate][state] and ghmm.sparsity[state][poststate]:
                    if strand == '-': beforestate, afterstate = poststate, prestate
                    else: beforestate, afterstate = prestate, poststate
                    cnctuples.append([beforestate, state, afterstate, frontpos, backpos])
    return cnctuples

def eclipsed(ghmm, strand, stopsignalframe, frame, state, pos):
    #print pos, frame, stopsignalframe[strand]
    if stopsignalframe[strand][frame]:
        if ghmm.states[state].backsignal == ghmm.stopcodon and pos == stopsignalframe[strand][frame][0] + ghmm.stopcodon.consensuslength:
            return False
        else:
            return pos >= stopsignalframe[strand][frame][0] + ghmm.stopcodon.consensuslength
    else:
        return False

def generateCodingStateList(ghmm, sequence):
    codingstates = {}
    for pos in range(len(sequence)+1):
        codingstates[pos] = {}
        for state in ghmm.ncstates:
            codingstates[pos][state] = []
    (frontsignals, stopsignalframe) = generateSignalList(ghmm, sequence)
    for strand in sequence:
        backsignals = frontsignals[strand][:]
        for [frontpos, frontsignal] in frontsignals[strand]:
            frontbacksignals = []
            maximumstop = updateStopList(ghmm, sequence, strand, stopsignalframe, frontpos, frontsignal)
            updateBackSignals(ghmm, sequence, backsignals, frontpos)
            for [backpos, backsignal] in backsignals:
                if backpos <= maximumstop+ghmm.stopcodon.consensuslength:
                    if frontpos + ghmm.signals[frontsignal].length <= backpos:
                        frontbacksignals.append([frontpos, frontsignal, backpos, backsignal])
                else:
                    break
            cnctuples = generateCNCTuples(ghmm, sequence, strand, frontbacksignals)
            for [prestate, state, poststate, frontpos, backpos] in cnctuples:
                codingstart, codingend = ghmm.states[state].signalOffset(frontpos, backpos, False)
                changeinphase = (ghmm.states[poststate].phase - ghmm.states[prestate].phase)%3
                frame = (codingstart - ghmm.states[prestate].phase)%3
                #print strand, codingstart, codingend, changeinphase, frame, 'pre/cur/post state:', prestate, state, poststate, 'maxstop', maximumstop
                if changeinphase == (codingend-codingstart)%3 and not eclipsed(ghmm, strand, stopsignalframe, frame, state, codingend):
                    if strand == '-':
                        prestate, poststate = poststate, prestate
                        frontpos, backpos = len(sequence) - backpos, len(sequence) - frontpos
                    codingstates[backpos][poststate].append([prestate, state, frontpos])
    return codingstates

def addNonCodingStates(ghmm, codingstates):
    validtransitions = {}
    for poststate in ghmm.ncstates:
        for prestate in ghmm.ncstates:
            if ghmm.sparsity[prestate][poststate]:
                if poststate in validtransitions:
                    validtransitions[poststate].append(prestate)
                else:
                    validtransitions[poststate] = []
                    validtransitions[poststate].append(prestate)
    for pos in range(len(codingstates)):
        for poststate in validtransitions:
            for prestate in validtransitions[poststate]:
                codingstates[pos][poststate].append([prestate, prestate, pos])

def decoder(ghmm, sequence, expression):
    print 'Preprocessing sequence'
    codingstates = generateCodingStateList(ghmm, sequence)
    addNonCodingStates(ghmm, codingstates)
    stdOutput('Decoding')
    viterbiStructure = ViterbiStructure(ghmm, sequence, expression)
    length = len(sequence)
    pos = 0
    for backpos in range(1, len(sequence)):
        if int(10*backpos/length) != pos:
            stdOutput('.')
            pos = int(10*backpos/length)
        for postidx in ghmm.ncstates:
            viterbiStructure.updateStructure(codingstates, backpos, postidx)
    viterbiStructure.terminateStructure(codingstates)
    print 'done'
    return viterbiStructure.generateParse()

def parsetoFeatures(ghmm, parse, features, header):
    statelist = []
    start = 0
    currentstate = None
    splithead = header.split(':')
    ref, startref, endref = splithead[0:3]
    #print parse
    for a in range(len(parse)):
        if currentstate != None:
            if currentstate != parse[a]:
                statelist.append((start, a-1, currentstate))
                currentstate = parse[a]
                start = a
        else:
            currentstate = parse[a]
    ngene = 0      
    for (start, end, state) in statelist:
        if state == 0: ngene += 1
        if ghmm.states[state].coding:
            (startoffset, endoffset) = ghmm.states[state].signalOffset(start, end, True)
            features.addFeature(Feature(ref, 'GHMM', ['exon'], [[startoffset+int(startref), endoffset+int(startref)]], None, ghmm.states[state].strand, None, 'Gene ' + str(ngene)))
        elif ghmm.states[state].tiling.id == 'Expressed':
            features.addFeature(Feature(ref, 'GHMM', ['ncoding'], [[int(startref) + start, int(startref) + end]], None, ghmm.states[state].strand, None, 'Gene ' + str(ngene)))

def test():
    a = GHMM('./ghmm7.cfg')
    testseq = SequenceDict(('AAAAAAAAAAAAAAGATGCTGATAGCTTATAAAGATAGAAAAAAAAAAAAAAAAAATAAGCGGTCGCTTCAAAAAAAAAAAAA'), 5)
    signals, stopframes = generateSignalList(a, testseq)
    codingstates = generateCodingStateList(a, testseq)
    addNonCodingStates(a, codingstates)
    for t in range(len(testseq)+1):
        for ncstate in a.ncstates:
            for [prevstate, state, front] in codingstates[t][ncstate]:
                seq = intToSeq(testseq['+'][0][front:t])
                print 'Prev state:', a.states[prevstate].symbol, 'State:', a.states[state].symbol, 'Post state:', a.states[ncstate].symbol, 'Front/back:', front, t, seq
