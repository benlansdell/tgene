"""Maths functions"""

from numpy import array, zeros, ones, eye, random, ndarray, hstack, fromstring, inf, log, random, sum, sqrt, argmax, float128, exp, pi

normalcdf = [float128(float(line)) for line in open('/home/users/lab0605/lansdell/TGene/stdnormal.data')]

def arrayDimension(a):
	"""Returns dimension of a numpy array"""
	if a.shape[0] == 0:
		return 0
	else:
		return len(a.shape)
	
def chooseRand(dist):
	"""Given a non-negative vector sum with 1(a discrete probability distribution) return randomly chosen element"""
	if sum(dist) != 1: raise ValueError, 'Distribution does not sum to 1'
	if reduce(lambda x,y:x*y, dist) < 0: raise ValueError, 'Non-negative probabilty'
	randstate = random.multinomial(1, dist)
	pos = map(lambda c,d:c * d, randstate, range(len(dist)))
	return reduce (lambda c,d: c + d, pos)

def normalCdf(x, mu, sigma):
	"""Normal cumulative distribution function"""
	#return 0.5 + 0.5*erf((x-mu)/sigma/sqrt(2))
	z = float(x - mu)/sigma
	if -10 <= z <= 0:
		return normalcdf[int(z*1000)+30000]
	elif 0 < z <= 10:
		return 1-normalCdf(-z, 0, 1)
	elif z < -10:
		return 0.0
	else:
		return 1.0

def discreteNormal(mu, sigma, xmin, xmax):
	"""Returns a vector (Pr(N = x),...)_{0<=x<xmax} where N ~ N(mu, sigma) is a discrete, normal distribution, scaled to be between xmin, xmax"""
	density = []
	scalefactor = normalCdf(xmax, mu, sigma) - normalCdf(xmin, mu, sigma)
	for x in range(xmin, xmax):
		density.append(normalCdf(x+1, mu, sigma) - normalCdf(x, mu, sigma))
	return array(density)/scalefactor

def normalpdf(x, mu, sigma):
	"""Gaussian density function"""
	return exp(-1*(x-mu)*(x-mu)/2/sigma)/sqrt(2*pi*sigma)

def ctsNormal(mu, sigma, minlist, maxlist):
	"""Returns a vector (f(x),...)_{0<=x<xmax} where x is an integer, N ~ N(mu, sigma) is a continutuous, normal distribution, scaled to be between xmin, xmax"""
	density = []
	for x in range(minlist, maxlist):
		density.append(normalpdf(x+1, mu, sigma)/(normalCdf(maxlist, mu, sigma) - normalCdf(minlist, mu, sigma)))
	return array(density)
	
def continuousKernelDensity(list, sigma):
	"""Returns a vector of length max(list), the density f(x)_{0<x<max}, similar to kernelDensity, but sums densities, not probabilities"""
	if min(list) < 0: minlist = int(min(list))
	else: minlist = 0
	maxlist = int(max(list))+1
	density = zeros(maxlist-minlist, float128)
	for val in list:
		density += ctsNormal(val, sigma, minlist, maxlist)
	return (density/len(list), minlist)

def kernelDensity(list, scale_bounds=1, c=0.05, d=None):
	"""Returns a discrete density vector of length max based on guassian kernal density estimates of data in list"""
	#Hard to take log(prob) using this sum method, means that really small probabilities are set to zero, is this a problem?
	if min(list) < 0: minlist = int(min(list)*scale_bounds)
	else: minlist = 0
	maxlist = int(max(list)*scale_bounds) + 1
	density = zeros(maxlist-minlist, float128)
	s = list[:]; s.sort()
	if d == None:
		N = len(list)
		d = (s[int(N*0.75)] - s[int(N*0.25)])
		print "Setting d to", d
	for val in list:
		sigma = abs(val)*float(c) + float(d)
		density += discreteNormal(val, sigma, minlist, maxlist)
	return (density/sum(density), minlist)