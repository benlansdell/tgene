"""For dealing with GFF files/objects"""

from io import smartopen
from sequence import replaceString

class Feature:
	"""A GFF feature class, each feature is a gene. Coords stores exon coordinates as list."""
	def __init__(self, ref, source, etype, coords, score, strand, frame, name):
		self.ref = ref
		self.source = source
		self.strand = strand
		self.name = name
		self.coords = coords
		if type(etype) == list: self.type = etype
		else: self.type = [etype]
		if type(score) == list: self.score = score
		else: self.score = [score]
		if type(frame) == list: self.frame = frame
		else: self.frame = [frame]
		self.min = min([c[0] for c in coords])
		self.max = max([c[1] for c in coords])
		self.coords.sort()
		
	def __eq__(self, other):
		if (self.coords == other.coords) and (self.ref == other.ref):
			return True
		else:	
			return False
		
	def __cmp__(self, other):
		if self.ref < other.ref:
			return -1
		if self.ref > other.ref:
			return 1
		if self.min < other.min:
			return -1
		if self.min > other.min:
			return 1
		return 0

	def __add__(self, other):
		#Adds 2 genes with same reference and name
		return Feature(self.ref, self.source, self.type + other.type, self.coords + other.coords, self.score + other.score, self.strand, self.frame + other.frame, self.name)
	
	def __repr__(self):
		stringrep = ''
		for a in range(len(self.coords)):
			if self.frame[a] == None: frame = '.'
			else: frame = str(self.frame[a])
			if self.score[a] == None: score = '.'
			else: score = str(self.score[a])
			stringrep += self.ref + '\t' + self.source + '\t' + self.type[a] + '\t' + str(self.coords[a][0]) + '\t' + str(self.coords[a][1]) + '\t' + score + '\t' + self.strand + '\t' + frame + '\t' + self.name + '\n'
		return stringrep
	
class FeatureDict(dict):
	"""Stores features in each reference. Overloads iterator so that features are returned in order of coordinate"""
	def __iterator(self):
		keys = [(self[generef].min, self[generef].name) for generef in self.keys()]
		keys.sort()
		keys = [a[1] for a in keys]
		for key in keys:
			yield key
		
	def __iter__(self):
		return self.__iterator()
		
class Features:
	"""Class to hold a dictionary of features read from a GFF file. Usage: features[ref][gene_name]"""
	def __init__(self, filename=None, types = '', features=None, ):
		if not features: self.features = {}
		else: self.features = features
		if filename:
			self.readGff(filename, types)
			
	def __repr__(self):
		stringrep = ''
		for ref in self.features:
			for gene in self.features[ref]:
				stringrep += str(self.features[ref][gene])
		return stringrep[:-1]
	
	def __iterfeature(self):
		for a in self.features:
			yield a

	def __iter__(self):
		return self.__iterfeature()
	
	def __getitem__(self, key):
		return self.features[key]
	
	def __setitem__(self, key, value):
		self.features[key] = value
	
	def addFeature(self, feature):
		ref = feature.ref
		name = feature.name
		if ref in self.features:
			if name in self.features[ref]:
				self.features[ref][name] += feature
			else:
				self.features[ref][name] = feature
		else:
			self.features[ref] = FeatureDict()
			self.features[ref][name] = feature
	
	def __splitLine(self, line):
		words = line.split('\t', 8)
		if len(words[8]):
			words[8] = words[8][:-1]
		return words
	
	def readGff(self, filename, types):
		"""Reads features from file. Will read in GFF, not GTF, so last column (gene name etc)
		must be identical between exons."""
		file = smartopen(filename)
		for line in file:
			if len(line) > 1:
				if line[0] != '#':
					words = self.__splitLine(line)
					try:
						ref = words[0]
						source = words[1]
						type = words[2]
						a = int(words[3])
						b = int(words[4])
						if words[5] == '.': score = None
						else: score = float(words[5])
						strand = words[6]
						frame = words[7]
						if words[7] == '.': frame = None
						else: frame = int(words[7])
						name = words[8]
					except StandardError, msg:
						raise StandardError, msg, 'Invalid GFF file format'
					start,end = min(a,b), max(a,b)
					if type in types or not types:
						self.addFeature(Feature(ref, source, type, [[start,end]], score, strand, frame, name))
						
	def writeGff(self, filename, reflist=None, mode = 'w'):
		"""Outputs list of gene features to gff file"""
		output = smartopen(filename, mode)
		if not reflist: reflist = self.features.keys()
		if type(reflist) == str: reflist = [reflist]
		for ref in reflist:
			if ref in self.features:
				for gene in self.features[ref]:
					output.write(str(self.features[ref][gene]))
		return 
		
	def exonList(self, reflist=None):
		"""Converts feature[ref]'s gene features to a list of exon coord pairs"""
		list = []
		if not reflist: reflist = self.features.keys()
		if type(reflist) == str: reflist = [reflist]
		for ref in reflist:
			for gene in self.features[ref]:
				for exon in self.features[ref][gene].coords:
					list.append([exon[0], exon[1]])
		list.sort()
		return list
	
	def exonString(self, ref, startref, endref):
		"""Returns string of length length, with F for feature, N for not 
		"""
		string = (endref-startref+1)*'N'
		for gene in self.features[ref]:
			for exon in self.features[ref][gene].coords:
				start = exon[0] - startref
				end = exon[1] - startref
				exonlength = end - start + 1
				exonstring = (exonlength)*'F'
				string = replaceString(string, exonstring, start)
		return string

#Needed?
def syncReferences(first, second):
	"""Ensures each feature object has the other's references, even if empty."""
	for ref in first:
		if ref not in second:
			second[ref] = {}
	for ref in second:
		if ref not in first:
			first[ref] = {}