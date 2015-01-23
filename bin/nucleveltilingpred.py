#!/usr/bin/env python
"""nucleveltilingpred.py <expression_file> <ghmm_cfg_file> <output_file>

Performs transcript mapping at nucleotide level using given ghmm parameters.
Uses a HMM rather than GHMM equivalent.

expression_file must be in format:
	>ref:start:end
	inf
	inf
	23.9
	...
	
ghmm_cgf_file must be in format described in the template.cfg file

Returns a gff file containing predicted transcribed and non-transcribed features
"""

from GHMM import *
from expressionstruct import *
from feature import *

def decoder(ghmm, express):
    T = len(express)
    transit = log(array([[1-9.4154599970058844e-03, 9.4154599970058844e-03],
        [0.0062224423584052134,1-0.0062224423584052134]]))
    initial = ghmm.initial
    delta = -inf*ones((T, ghmm.nstates))
    psi = -inf*ones((T, ghmm.nstates))
    qstar = -1*ones(T)
    #Initialise, t=0 position
    for idx in range(ghmm.nstates):
        state = ghmm.states[idx]
        delta[0, idx] = initial[idx]+state.tiling.probEmit(express, 0, 1)#/30
        psi[0,idx] = 0
    #Recursive, for 1<t<T-1. most of time spent here
    for t in range(1,T):
        for i in range(ghmm.nstates):
            tempdelta = []
            for j in range(ghmm.nstates):
                state = ghmm.states[i]
                tempdelta.append([delta[t-1,j]+transit[j,i]+state.tiling.probEmit(express, t, t+1), j])#remove /30 from probEmit...not sure if it's needed
                print 'Pos:', t, 'Curr State:', i, 'Old state:', j, 'Score:', delta[t-1,j]+transit[j,i]+state.tiling.probEmit(express, t, t+1)/30, 'This score:', transit[j,i]+state.tiling.probEmit(express, t, t+1)/30
            delta[t,i],psi[t,i] = max(tempdelta)
    #Terminate
    pstar = max(delta[T-1,:])
    qstar[T-1] = argmax(delta[T-1,:])
    #Backtrack
    for t in range(T-2,-1,-1):
        qstar[t] = psi[t+1,qstar[t+1]]
    return qstar

def parsetoFeatures(ghmm, parse, predictions, header):
    infeature = False
    print parse
    startref = int(header.split(':')[1])
    ref = header.split(':')[0]
    print 'Parsing', ref
    start = 0
    end = 0
    number = 0
    for idx in range(len(parse)):
        if parse[idx]:
            if infeature:
                end = idx + startref
            else:
                start = idx + startref
                end = idx + startref + 25
            infeature = True
        else:
            if infeature:
                number += 1
                predictions.addFeature(Feature(ref, 'nuclevelHMM', ['noncoding'], [[start, end]], None, '+', None, 'Gene ' + str(number)))
            infeature = False
    if infeature:
        number += 1
        predictions.addFeature(Feature(ref, 'nuclevelHMM', ['noncoding'], [[start, end]], None, '+', None, 'Gene ' + str(number)))

if len(sys.argv) != 4:
    print 'Invalid arguments'
    print __doc__
    sys.exit(0)
    
expression = Expression(sys.argv[1])
ghmm = GHMM(sys.argv[2])
predictions = Features()

for header in expression.expression:
    parse = decoder(ghmm, expression[header])
    parsetoFeatures(ghmm, parse, predictions, header)

predictions.writeGff(sys.argv[3])