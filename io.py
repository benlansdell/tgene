from maths import *

def smartopen(filehandle, mode='r'):
	"""Opens a file, unless passed an already open filehandle
	
	@param filehandle: String or filehandle
	@param mode: Filemode
	@returns: Open filehandle
	"""
	if type(filehandle) == str:
		return open(filehandle, mode)
	else:
		return filehandle

def writeMatrix(mat, filehandle):
	"""Writes numpy array to a file
	
	@param mat: numpy vector or matrix to write
	@param filehandle: Filehandle
	@returns: True if successful
	"""
	if arrayDimension(mat) == 1:
		mat = array([mat])
	if arrayDimension(mat) in (0,1,2):
		outputmatrix = hstack((mat.shape, mat.ravel()))
	else:
		raise ValueError, 'Array dimension = ' + str(arrayDimension(mat)) + '). Empty, 1d and 2d matrices only.'
	output = smartopen(filehandle, 'w')	
	columns = len(outputmatrix)
	if mat.shape == (0,):
		output.write('0.0\n')
		output.close()
		return True
	for column in range(columns):
		output.write(str(outputmatrix[column]) + '\t')
	output.write('\n')
	output.close()
	return True
		
def readMatrix(filehandle):
	"""Read a matrix from file
	
	@param filehandle: File to read matrix from
	@returns: Numpy vector or matrix"""
	infile = smartopen(filehandle)
	line = infile.readline()
	if line == '0.0\n': return array([])
	inputrow = fromstring(line, float, -1, '\t')
	rows,columns = inputrow[0],inputrow[1]
	matrix = inputrow[2:2+rows*columns]	
	matrix.shape = (rows, columns)
	if (rows == 1):
		matrix = matrix[0]
	return matrix