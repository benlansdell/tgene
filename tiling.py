"""tiling.py
Contains models for tiling data."""

from numpy import correlate, corrcoef, column_stack, array
from numpy.random import rand, random
from scipy.stats.stats import percentileofscore

def computeIntStats(data, maxprobes):
	return map(sum, data.expression[0:maxprobes])

#def percentile(scores):
#	testscores = array(range(-100,100))/float(1)
#	testpctiles = map(lambda x: percentileofscore(scores, x, 100), list(testscores))
#	percentiles = []
#	for score in scores:
#		print int((score + 100))
#		percentiles.append(testpctiles[int((score+100)*1)])
#	return percentiles
	#return map(lambda x: percentileofscore(scores, x, 100), scores)

def percentile(scores):
	rankscores = [(scores[a], a) for a in range(len(scores))]
	rankscores.sort()
	percentilescores = [(rankscores[a][1], a) for a in range(len(scores))]
	percentilescores.sort()
	return [a[1]/float(len(scores)) for a in percentilescores]
	
def sumints(data, position1):
	return sum(data.expression[position1%len(data)])

def unbiasedcorr(data, position1, position2):
    #Takes 
    x = data.expression[position1%len(data)]
    y = data.expression[position2%len(data)]
    randx = rand(12,3)/10000
    randy = rand(12,3)/10000
    x.shape = (12,3)
    y.shape = (12,3)
    x = x + randx
    y = y + randy
    xi = sum(x, axis = 1)/3
    xi = xi.transpose()
    xi = column_stack([xi, xi, xi])
    yi = sum(y, axis = 1)/3
    yi = yi.transpose()
    yi = column_stack([yi, yi, yi])
    xij = sum(x)/36
    yij = sum(y)/36
    ssbx = sum([(xi[i,j] - xij)*(xi[i,j] - xij) for i in range(12) for j in range(3)])/11
    sswx = sum([(x[i,j] - xi[i,j])*(x[i,j] - xi[i,j]) for i in range(12) for j in range(3)])/24
    ssby = sum([(yi[i,j] - yij)*(yi[i,j] - yij) for i in range(12) for j in range(3)])/11
    sswy = sum([(y[i,j] - yi[i,j])*(y[i,j] - yi[i,j]) for i in range(12) for j in range(3)])/24
    spb = sum([(xi[i,j] - xij)*(yi[i,j] - yij) for i in range(12) for j in range(3)])/11
    spw = sum([(x[i,j] - xi[i,j])*(y[i,j] - yi[i,j]) for i in range(12) for j in range(3)])/24
    varxi = (ssbx - sswx)/3
    vareta = (ssby - sswy)/3
    cov = (spb - spw)/3
    if varxi <= 0 or vareta <= 0:
        corr = 0
    else:
        corr = cov/sqrt(varxi * vareta)
    return corr
    
def fstats(x):
    #Takes 
    randx = rand(3,12)/10000
    x = x + randx    
    xi = sum(x, axis = 1)/12
    xj = sum(x, axis = 0)/3
#    print xi, xj
    xij = sum(x)/36
    ssbx = sum([(xj[j] - xij)*(xj[j] - xij) for i in range(3) for j in range(12)])
    sswx = sum([(x[i,j] - xi[i])*(x[i,j] - xi[i]) for i in range(3) for j in range(12)])
    fscore = ssbx/sswx
    return fscore    

def spf(data, position1, w):
    
    x1 = [data.expression[(position1+i)%len(data)] for i in range(w)]
    x2 = [data.expression[(position1-i)%len(data)] for i in range(w)]
    return max(fstats(array(x1)), fstats(array(x2)))
    
def correlation(data, position1, position2):
    #What about negative values? Negative correlation in this case is bad...
    #What to do when one std deviation is zero?? Ignore that point...
    #it tells you nothing about correlation. Or...just set it to zero. Assume uncorrelated.
    x = data.expression[position1%len(data)]
    y = data.expression[position2%len(data)]
    randx = (2*rand(len(x))-1)/10000
    randy = (2*rand(len(x))-1)/10000
    x = x + randx
    y = y + randy
    cceff = corrcoef(x, y)[0,1]
    if cceff <0 or cceff >= 0: pass
    else: cceff = 0#2*random()-1
    
    stringout = str(cceff) + '\t'
    for val in list(x):
        stringout += str(val) + '\t'
    stringout += '\n' + str(cceff) + '\t'
    for val in list(y):
        stringout += str(val) + '\t'
    if cceff < 0.8 and cceff > 0.3: pass#print stringout
    return cceff

def computeSMC(data, position, w):
    smcleft = [correlation(data, position + a, position) for a in range(-w, 0)]
    smcright = [correlation(data, position + a, position) for a in range(1,w+1)]
    #smcleft = [unbiasedcorr(data, position + a, position) for a in range(-w, 0)]
    #smcright = [unbiasedcorr(data, position + a, position) for a in range(1,w+1)]
    #score = spf(data, position, w)
    #return score
    #return sumints(data, position)
    return 100*max(sum(smcleft)/w, sum(smcright)/w)
    
def computeMISMC(SMCpctile, Intpctile):
    mismc = []
    for a in range(len(SMCpctile)):
        mismc.append(100*max(SMCpctile[a], Intpctile[a]))
    return mismc

def computeSMCstats(data, maxprobes, w):
    return [computeSMC(data, pos, w) for pos in range(maxprobes)]

def analyseData(data, filename, maxprobe = None):
        """Output statistics to a file"""
        if not maxprobe: maxprobe = len(data)
        maxprobes = min(maxprobe, len(data))
        outfile = open(filename, 'w')
        #outfile.write('Position\tSMC\tSMC%\tInt\tInt%\tWithin-gene\n')
        SMCstats = computeSMCstats(data, maxprobes, 3)
        Intstats = computeIntStats(data, maxprobes)
        print 'Computed stats'
        #SMCpctile = percentile(SMCstats)
        #Intpctile = percentile(Intstats)
        #print 'Computed percentiles'
        SMCpctile = percentile(SMCstats)
        #print 'Computed SMC %'
        Intpctile = percentile(Intstats)
        #print 'Computed percentiles'
        MISMC = computeMISMC(SMCpctile, Intpctile)
        print 'Writing to file'
        for a in range(maxprobes):
                #outfile.write(str(data.positions[a]) + '\t' + str(SMCstats[a]) + '\t' + str(SMCpctile[a]) + '\t' + str(Intstats[a]) + '\t' + str(Intpctile[a]) + '\t' + str(MISMC[a]) + '\t' + str(data.annotation[a]) + '\n')
                outfile.write(str(data.positions[a]) + '\t' + str(MISMC[a]) + '\t' + str(Intstats[a]) + '\t' + str(data.annotation[a]) + '\n')

class EXPRSMC:
    """Expression Markov chain"""
    def __init__(self, files):
        filelist = files.split()[0].split(',')
        self.nfiles = len(filelist)
        if self.nfiles != 1:
            raise ValueError, 'EXPRSMC requires 1 file: ' + files
        #(self.prob, self.minval) = kernelDensity([int(float(line.split()[0])) for line in open(filelist[0])], 6, 0, 20)
        #(self.prob, self.minval) = kernelDensity([int(float(line.split()[0])) for line in open(filelist[0])], 1, 0, 7)
        (self.prob, self.minval) = continuousKernelDensity([int(float(line.split()[0])) for line in open(filelist[0])], 7)
        self.prob = log(self.prob)
        self.emissionprob = log(float(1)/30)
        self.noemissionprob = log(float(29)/30)
    def probEmit(self, expression, start, end):
        #Need to add geometric cut-off tail...
        score = 0
        for val in expression[start:end]:
            if val != inf:
                if int(val) >= len(self.prob)-self.minval or int(val) < self.minval: return -inf
                else: 
                    try: score += self.prob[int(val)-self.minval] + self.emissionprob
                    except: print val, self.minval, len(self.prob)
            else: 
                score += self.noemissionprob
        return score
    
class EXPRSNONE:
    """Expression probabilities for when tiling data isn't modeled"""
    def __init__(self, files):
        pass
    def probEmit(self, expression, start, end):
        return 0.0

tilingmodeldict = {'EXPRSMC':EXPRSMC, 'EXPRSNONE':EXPRSNONE}
