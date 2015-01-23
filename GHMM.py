"""GHMM.py
contains class for GHMM framework"""

from state import *
from string import lowercase, uppercase

class Config:
    """Manages and loads data from ghmm.cfg file"""
    def __init__(self, filename):
        self.ghmm = {}
        self.states = {}
        self.signals = {}
        self.content = {}
        self.tilings = {}
        self._readConfig(filename)
        
    def _parseLine(self, line):
        parse = line.split("#", 1)
        if len(parse) == 0:
            return '', ''
        elif len(parse) == 1:
            return parse[0], ''
        else:
            return parse
    
    def _convertValue(self, value):
        try: value = int(value)
        except: pass
        if value == 'None': value = None
        if value == 'True': value = True
        if value == 'False': value = False
        return value
    
    def _keyCase(self, key):
        lowercount = sum([i in lowercase for i in key])
        uppercount = sum([i in uppercase for i in key])
        if lowercount + uppercount == len(key):
            if lowercount == 0:
                return 'UPPER'
            elif uppercount == 0:
                return 'lower'
            else:
                return 'Camel'
        return 'other'

    def _readConfig(self, filename):
        """Config file contains all details of GHMM, including state submodels and signal sensors"""
        infile = smartopen(filename)
        current = self.ghmm
        for line in infile:
            code, comment = self._parseLine(line)
            keyvalues = code.split()
            if keyvalues:
                if keyvalues[0]:
                    if keyvalues[0][0] == '[':
                        currentkey = keyvalues[0].split(']')[0][1:]
                        try:
                            self.states[int(currentkey)] = {}; current = self.states[int(currentkey)]
                        except:
                            keycase = self._keyCase(currentkey)
                            if keycase == 'UPPER':
                                self.signals[currentkey] = {}; current = self.signals[currentkey]
                            elif keycase == 'lower':
                                self.content[currentkey] = {}; current = self.content[currentkey]
                            elif keycase == 'Camel':
                                self.tilings[currentkey] = {}; current = self.tilings[currentkey]
                            else:
                                raise ValueError, 'Invalid attribute key in .cfg file. Must be numeric, uppercase, lowercase or camelcase: ' + currentkey
            if len(keyvalues) == 3:
                if keyvalues[1] == '=':
                    key = keyvalues[0]
                    value = self._convertValue(keyvalues[2])
                    current[key] = value

class GHMM:
    """GHMM class contains model topology in transition array and 
    states, signals and cotent to be loaded"""
    def __init__(self, filename):
        self.states = {}
        self.signals = {}
        self.content = {}
        self.tilings = {}
        config = Config(filename)
        if 'weight' in config.ghmm: self.weight = float(config.ghmm['weight'])
        else: self.weight = 0.0
        if 'newweight' in config.ghmm: self.newweight = config.ghmm['newweight']
        else: self.newweight = False
        print 'Setting newweight to', self.newweight
        for key in config.tilings:
            self.tilings[key] = Tiling(key, config.tilings[key])
        for key in config.signals:
            self.signals[key] = Signal(key, config.signals[key])
        for key in config.content:
            self.content[key] = Content(key, config.content[key])
        for key in config.states:
            self.states[key] = State(key, config.states[key], self)
        self.stopcodon = self.signals[config.ghmm['stopcodon']]
        self.transition = self._loadMatrix(config.ghmm['transition'])
        self.initial = self._loadMatrix(config.ghmm['initial'])
        self.nstates = config.ghmm['nstates']
        self._loadSparse()
        self._loadCNC()
    
    def _loadMatrix(self, file):
        transition = readMatrix(file)
        if arrayDimension(transition) == 1:
            rowsum = sum(transition)
            if rowsum == 0: rowsum = 1
            transition /= rowsum
            return log(transition)
        for row in range(len(transition)):
            rowsum = sum(transition[row])
            if rowsum == 0: rowsum = 1
            transition[row] = transition[row] / rowsum
        return log(transition)
            
    def _loadCNC(self):
        """Seperate coding/noncoding states"""
        self.ncstates = []
        self.cstates = []
        for state in self.states:
            if not self.states[state].coding:
                self.ncstates.append(state)
            else:
                self.cstates.append(state)
                  
    def _loadSparse(self):
        self.sparsity = zeros((self.nstates, self.nstates))
        for i in range(self.nstates):
            for j in range(self.nstates):
                if self.transition[i][j] != -inf:
                    self.sparsity[i][j] = 1
