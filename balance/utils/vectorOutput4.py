# ---------------------------------------------------------------------------- #
def generateOutputMatrix(vectorList, delimeter='\t'):
	
	outputMatrix = []
	outputString = ''
	i = 0
	
	
	
	while i < len(vectorList[0]):
		
		outputString = ''
		j = 0
		
		while j < len(vectorList):
			outputString += str(vectorList[j][i])
			if j < len(vectorList) - 1:
				outputString += delimeter
			j += 1
		
		outputString += '\n'
		outputMatrix.append(outputString)
		i += 1
	
	return outputMatrix
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
def generateOutputMatrixWithHeaders(vectorList, headers, delimeter='\t'):
	
	if len(headers) != len(vectorList):
		print("Header array length != vector list length")
		return
	
	outputMatrix = []
	
	k = 0
	headerString = ''
	while k < len(headers):
		headerString += str(headers[k]) + delimeter
		k += 1
	
	headerString += '\n'
	outputMatrix.append(headerString)
	
	outputString = ''
	i = 0
	
	while i < len(vectorList[0]):
		
		outputString = ''
		j = 0
		
		while j < len(vectorList):
			outputString += str(vectorList[j][i])
			if j < len(vectorList) - 1:
				outputString += delimeter
			
			j += 1
		
		outputString += '\n'
		outputMatrix.append(outputString)
		i += 1
	
	return outputMatrix
# ---------------------------------------------------------------------------- #




# ---------------------------------------------------------------------------- #
def writeOutputMatrix(fileName, outputMatrix):
	outputFile = open(fileName, 'w')
	i = 0
	while i < len(outputMatrix):
		outputFile.write(outputMatrix[i])
		i+=1
	return	
# ---------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
def Write_SMatrix(fileName, sMatrixT, ioStatus, compounds, reactions):

	outputStr = ''
	
	# Output the first line, which contains the compound names. Remember it starts with an empty
	# cell for the column where the reactions will go. 
	outputStr += ','
	i = 0
	while i < len(compounds):
		outputStr += compounds[i]
		if i < len(compounds) - 1:
			outputStr += ','
		i += 1
	outputStr += '\n'
	
	# Output the iostatus data. Remember it starts with an empty cell for the column where the 	
	# reactions will go. 
	outputStr += ','
	i = 0
	while i < len(ioStatus):
		outputStr += ioStatus[i]
		if i < len(ioStatus) - 1:
			outputStr += ','
		i += 1
	outputStr += '\n'
	
	i = 0
	while i < len(reactions):
		outputStr += reactions[i] + ','
		j = 0
		while j < len(compounds):
			outputStr += str(sMatrixT[i][j])
			if j < len(compounds) - 1:
				outputStr += ','
			j += 1
		outputStr += '\n'
		i += 1
	
	oHandle = open(fileName, 'w')
	oHandle.write(outputStr)
	oHandle.close()
	
	return
# ----------------------------------------------------------------------------------------------- #

