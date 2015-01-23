#!/usr/bin/env python
"""probeleveltilingpred.py <expression_file> <ghmm_cfg_file> <output_file> <offset> <factor>

Performs transcript mapping at probe level using given ghmm parameters. 
Expression file must be two columns, position and intensity.

Performs transcript mapping at probe level using given ghmm parameters.
Uses a HMM rather than GHMM equivalent.

expression_file must be in format:
3	90.4
34	93.2
61	87.3
...
	
ghmm_cgf_file must be in format described in the template.cfg file

Some expression files are not properly formatted so:

offset is added to predicted features

Intensity scores are multiplied by factor

These should otherwise default to 0 and 1 respectively.

Returns a gff file containing predicted transcribed and non-transcribed features"""

from GHMM import *
from expressionstruct import *
from feature import *

def decoder(ghmm, express):
    T = len(express)
    transit = log(array([[0.95,0.05],
        [0.20 ,0.80]]))
    initial = ghmm.initial
    delta = -inf*ones((T, ghmm.nstates))
    psi = -inf*ones((T, ghmm.nstates))
    qstar = -1*ones(T)
    #Initialise, t=0 position
    for idx in range(ghmm.nstates):
        state = ghmm.states[idx]
        minval = state.tiling.model.minval
        delta[0, idx] = initial[idx]+state.tiling.model.prob[int(factor*express[0])-minval]
        psi[0,idx] = 0
    #Recursive, for 1<t<T-1. most of time spent here
    for t in range(1,T):
        for i in range(ghmm.nstates):
            tempdelta = []
            for j in range(ghmm.nstates):
                state = ghmm.states[i]
                minval = state.tiling.model.minval
                tempdelta.append([delta[t-1,j]+transit[j,i]+state.tiling.model.prob[int(factor*express[t])-minval], j])
                #print 'Pos:', t, 'Curr State:', i, 'Old state:', j, 'Total score:', delta[t-1,j]+transit[j,i]+state.tiling.model.prob[int(express[t])], 'This score:', transit[j,i]+state.tiling.model.prob[int(express[t])]
            delta[t,i],psi[t,i] = max(tempdelta)
    #Terminate
    pstar = max(delta[T-1,:])
    qstar[T-1] = argmax(delta[T-1,:])
    #Backtrack
    for t in range(T-2,-1,-1):
        qstar[t] = psi[t+1,qstar[t+1]]
    return qstar

def parsetoFeatures(parse, expression, predictions):
    infeature = False
    startref = int(sys.argv[4])
    start = 0
    end = 0
    number = 0
    for idx in range(len(parse)):
        if parse[idx]:
            if infeature:
                end = expression.positions[idx] + startref
            else:
                start = expression.positions[idx] + startref
                end = expression.positions[idx] + 25 + startref
            infeature = True
        else:
            if infeature:
                number += 1
                predictions.addFeature(Feature('Adh', 'probelevelHMM', ['noncoding'], [[start, end]], None, '+', None, 'Gene ' + str(number)))
            infeature = False
    if infeature:
        number += 1
        predictions.addFeature(Feature('Adh', 'probelevelHMM', ['noncoding'], [[start, end]], None, '+', None, 'Gene ' + str(number)))

if len(sys.argv) != 6:
    print 'Invalid arguments'
    print __doc__
    sys.exit(0)

factor = int(sys.argv[5])
expression = ExpressionData([sys.argv[1]])
ghmm = GHMM(sys.argv[2])
predictions = Features()
print 'Predicting genes in ' + sys.argv[1] + '\n'
#rows = (array(expression.positions) > 13506243) * (array(expression.positions) < 16425263)
parse = decoder(ghmm, expression.expression)
parsetoFeatures(parse, expression, predictions)
predictions.writeGff(sys.argv[3])