from feature import Features
from io import *
import sys

from maths import *

class ExpressionData:
    def __init__(self, exprsfiles, annotation=None):
        """exprsfiles a list of files into a numpy array object. Each file must consist of two columns, one with position, one with intensity.
        arrayindex is a list of indices such that replicates will have the same value.
        annotation is a gff file containing annotation to label probes."""
        exprs = []
        self.positions = []
        loadedpositions = False
        self.names = []
        for filename in exprsfiles:
            exprs.append([])
            self.names.append(filename)
            for line in open(filename): 
                words = line.split()
                if 1:#int(words[0]) > 13506243 and int(words[0]) < 16425263:
                    if not loadedpositions: self.positions.append(int(words[0]))
                    intensity = float(words[1])
                    exprs[-1].append(intensity)
            loadedpositions = True
        self.expression = array(exprs).transpose()
        self.annotateExpression(annotation)
        
    def __len__(self):
        return len(self.positions)
    
    def annotateExpression(self, annotation=None, probelength=25):
        """Computes a boolean vector indicating which probes are (entirely) within annotated genes."""
        withingene = []
        if type(annotation) == str:         
            exonlist = Features(annotation).exonList()
        elif type(annotation) == list:
            exonlist = annotation
        else:
            self.annotation = [False]*len(self)
            return

        for probep in self.positions:
            probe = probep
            value = False
            if exonlist:
                while exonlist[0][1] <= probe:
                    exonlist.pop(0)
                    if not exonlist: break
                if exonlist:
                    if exonlist[0][0] <= probe and exonlist[0][1] >= probe+probelength:
                        value = True
            withingene.append(value)
        self.annotation = withingene
        return withingene

    def listWithinGenes(self, withingene = True, output = sys.stdout, idx = None):
        output = smartopen(output, 'w')
        if idx == None: idx = range(len(self.names))
        elif type(idx) == int: idx = [idx]
        for j in idx:
                for i in range(len(self)):
                    if self.annotation[i] == withingene:
                        print>>output, self.expression[i,j]

class Expression:
    def __init__(self, infile):
        expression = smartopen(infile)
        self.expression = {}
        currentkey = ''
        for line in expression:
            if line:
                if line[0] == '>':
                    currentkey = line[1:-1]
                elif currentkey in self.expression:
                    self.expression[currentkey].append(float(line))
                else:
                    self.expression[currentkey] = []
                    self.expression[currentkey].append(float(line))
    def __getitem__(self, key):
        return self.expression[key]