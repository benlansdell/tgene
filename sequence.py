"""sequence.py
Common sequence functions"""

seqdict = { 'A':0, 'C':1, 'G':2, 'T':3,
	     0:'A', 1:'C', 2:'G', 3:'T' }
complement = { 'A':'T', 'C':'G', 'G':'C', 'T':'A',
	      0:3, 1:2, 2:1, 3:0}

aminoacids = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }

nuclist = ['A', 'C', 'G', 'T']

def seqToInt(seq):
	"""Converts seq to list of integers.
	
	@param seq: String of uppercase ACGT
	@returns: List of numbers 0123
	"""
	seqlist = []
	for nt in seq:
 		if not nt in 'ACGT':
			raise ValueError, 'Invalid sequence input; upper case ACGT only.'
		seqlist.append(seqdict[nt])
	return seqlist
	
def intToSeq(seqlist):
	"""Converts list of integers to DNA sequence
	
	@param seqlist: List of numbers 0123
	@returns: String of uppercase ACGT
	"""
	seq = ''
	for int in seqlist:
		if not int in [0,1,2,3]:
			raise ValueError, 'Invalid sequence input; [0,1,2,3] only.'
		seq = seq + seqdict[int]
	return seq

def translate(seq):
	"""Translate DNA sequence to amino-acid sequence
	
	@param seq: list of integers 0123 or string of ACGT
	@returns: string of amino-acid residues
	"""
	if type(seq) == list: sequence = intToSeq(seq)
	else: sequence = seq
	translation = ''
	if len(seq)%3 != 0:
		raise ValueError, 'Sequence length must be multiple of 3:' + str(len(seq))
	for nt in sequence:
		if not nt in 'ACGT':
			raise ValueError, 'Invalid seuqence input; ACGT or [0,1,2,3] only.'
	for idx in range(len(seq)/3):
		translation += aminoacids[sequence[3*idx:3*idx+3]]
	return translation

def reverse(seq):
	"""Reverse sequence
	
	@param seq: List of integers 0123 or string of ACGT
	@returns: Reversed sequence
	"""
	string = False
	if type(seq) == str:
		string = True
		seq = list(seq)
	seq.reverse()
	if string: seq = ''.join(seq)
	return seq

def reverseComplement(seq):
	"""Reverse and complement sequence

	@param seq: List of integers 0123 or string of ACGT
	@returns: Reversed complement sequence
	"""
	string = False
	if type(seq) == str: string = True
	try:
		seq = [complement[a] for a in seq]
	except KeyError: raise ValueError, 'Invalid sequence input; ACGT or [0,1,2,3] only.'
	if string: seq = ''.join(seq)
	return reverse(seq)

def replaceString(astring, repstring, pos):
	"""Replace part of string with another string
	
	String slicing with a little more error checking.
	
	@param astring: original string
	@param repstring: replacement string
	@param pos: position to insert repstring
	@returns: astring with repstring inserted
	"""
	length = len(astring)
	if (pos >= 0) and (pos < len(astring)): 
		astring = astring[0:pos] + repstring + astring[pos+len(repstring):]
	return astring[0:length]

def _decimal(seq):
	"""Converts a base-four number to a decimal"""
	value = 0
	count = 1
	for a in seq:
		value += count*a
		count *= 4
		if not a in [0,1,2,3]:
			raise ValueError, 'Invalid sequence input; [0,1,2,3] only.'

	return value

def _basefour(value, order):
	"""Converts a decimal to an base-four integer sequence. DEPRECATED"""
	if value == None:
		return None
	seq = []
	denom = pow(4,order)
	for a in range(order):
		seq = [value/denom] + seq
		value = value - denom*(value/denom)
		denom /= 4
	seq = [value] + seq
	return seq

def _genHigherOrder(seq, order):
	"""Generates list of symbols corresponding to higher order emissions in sequence"""
	sequencen = []
	if len(seq) < order+1: raise ValueError, 'Sequence too short for order.'
	for a in range(len(seq)):
		if a-order < 0:
			sequencen.append(None)
		else:
			sequencen.append(_decimal(seq[a-order:a+1]))
	return sequencen

def _genLowerOrder(seq, order):
	"""From a list of higher order emissions, generate lower (0th) order emissions. DEPRECATED"""
	overlapseq = []
	retseq = [0]*len(seq)
	for a in seq:
		overlapseq.append(_basefour(a, order))
	for a in range(len(overlapseq)):
		if not overlapseq[a]: continue
		retseq = replaceString(retseq, overlapseq[a], a-order)
	return retseq

def _higherOrderDict(seq, maxorder):
	"""From a string or list of 0th order emissions create a dictionary of form
	dict[order] = list of emissions"""
	dict = {}
	for order in range(maxorder + 1):
		dict[order] = _genHigherOrder(seq, order)
	return dict

class SequenceDict:
	"""Store sequence as a dictionary of higher-order emissions.
	
	Stored as SequenceDict[strand][order]
	"""
	def __init__(self, seq, maxorder):
		"""@param seq: List of 0123 or string of ACGT
		@param maxorder: Highest order to store emissions 
		"""
		self.dictionary = {}
		if type(seq) == str: sequence = seqToInt(seq)
		else: sequence = seq
		self.dictionary['+'] = _higherOrderDict(sequence, maxorder)
		self.dictionary['-'] = _higherOrderDict(reverseComplement(sequence), maxorder)
		
	def __len__(self):
		return len(self.dictionary['+'][0])
	
	def __repr__(self):
		return intToSeq(self.dictionary['+'][0])
	
	def __iter__(self):
		return self._iter()
	
	def _iter(self):
		for key in self.dictionary:
			yield key
	
	def __getitem__(self, key):
		return self.dictionary[key]
	
	def __setitem__(self, key, value):
		self.dictionary[key] = value
	
