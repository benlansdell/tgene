"""lengthdist.py
State length distributions"""

from maths import *
from io import *

class LengthDist:
    """Base length distribution class. GHMM module will provide constructor with string
    'params' taken from ghmm.cfg file."""
    def probLength(self, length):
        """Log of probability of emitting a state of length"""        
        raise NotImplementedError
    def probCumLength(self):
        """Log of probability of continuing in current state for 1bp"""    
        raise NotImplementedError

class Geometric(LengthDist):
    """Geometric distribution. Parameterised by probability of success(change of state)"""
    def __init__(self, params):
        """params = prob, where prob is probability of success"""
        probsuccess = float(params)
        if probsuccess < 0 or probsuccess > 1:
            raise ValueError, 'Input probability ' + str(probsuccess) + ' needs to be 0<p<1'
        self.prob = log(probsuccess)
        self.probcomplement = log(1-probsuccess)
    def probLength(self, length):
        if length < 0:
            raise ValueError, 'Length ' + str(length) + ' needs to be >= 0'
        return self.probcomplement*length + self.prob
    def probCumLength(self):
        return self.probcomplement

class Fixed(LengthDist):
    """params = length, where length is fixed length"""
    def __init__(self, params):
        length = int(params)
        if length < 1:
            raise ValueError, 'Input length ' + str(params) + ' needs to be > 0'
        self.length = params
    def probLength(self, length):
        if length == self.length: return 0
        else: return -inf
        
class EmpiricalKDE(LengthDist):
    """Empirical kernel density estimates. params = filename, file to read in data - as a list of lengths"""
    #The question is: what to do about lengths past the maximum? Need geometric decay...add later
    def __init__(self, params):
        print 'Loading length distributions from', params
        self.lengths = [int(line) for line in open(params)]
        (self.density, self.minval) = kernelDensity(self.lengths, 1.5, 0.05, 250)
        if self.minval < 0: raise ValueError, 'Negative length encountered in length distribution estimates.'
    def probLength(self, length):
        if length >= len(self.density): return log(0)
        else: return log(self.density[length])
        
lengthdict = {'geometric':Geometric, 'fixed': Fixed, 'empiricalkde':EmpiricalKDE}