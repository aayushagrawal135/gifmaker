import numpy as np
import cv2

def findMed(extA, updn, forwardMed, iter):
	newA = []
	#if high limit
	if updn == 1:
		for x in range(len(extA)):
			if extA[x] > forwardMed:
				newA.append(extA[x])
	
	#if low limit
	else:
		for x in range(len(extA)):
			if(extA[x] < forwardMed):
				newA.append(extA[x])

	median = np.median(newA)
	if iter == 0:
		return median
	else:
		return findMed(newA, updn, median, iter-1)

def zStretch(A, iter):
	#matrix A is 2D, make it 1D by appending columns one after the other to make a single row
	rows, cols= np.shape(A)
	extA = np.reshape(A, (rows*cols))
	
	firstMed = np.median(extA)

	lowLimit = findMed(extA, -1, firstMed, iter)
	highLimit = findMed(extA, 1, firstMed, iter)

	for x in range(rows):
		for y in range(cols):
			if A[x][y] > highLimit:
				A[x][y] = highLimit
			elif A[x][y] < lowLimit:
				A[x][y] = lowLimit

	#noralize
	A = A - np.amin(A)
	A = np.divide(A, np.amax(A))
	A = np.multiply(A, 255)

	A = A.astype(np.uint8)
	return A

# converts 2D array of rows*cols into 1D array of #elements = rows*columns 
def vecize(inMat):
	rows, cols = np.shape(inMat)
	inMat = np.transpose(inMat)
	outvec = np.reshape(inMat, rows*cols)
	return outvec

def toGray(filename):
	gray_img = cv2.cvtColor(filename, cv2.COLOR_RGB2GRAY)
	return gray_img