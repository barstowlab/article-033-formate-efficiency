import re
import os.path
import pdb

from numpy import unique, int8, array

from utils.balanceUtils import SolveFluxBalanceEquation, ConvertIndexedSMatrix, ImportIOStatus, \
ImportReactionFile, GenerateIndexedSMatrixT, PrintStoichiometry, GenerateMergedIOStatusList, \
Get_IO_Status_for_Compound_List, ExportUniqueCompoundsWithIOStatus, GenerateIOStatusList, ParseReactionList

from utils.vectorOutput4 import generateOutputMatrixWithHeaders, writeOutputMatrix, Write_SMatrix
from utils.specutils12 import ensure_dir

# P2O7 = Diphosphate
# TrxOx = Thioredoxin disulfide
# PO43- = Phosphate

reactionList = \
[\
'THF + HCO2- + ATP -> ADP + PO43- + 10-Formyl-THF', \
\
'10-Formyl-THF + Hplus -> 5-10-Methenyl-THF + H2O', \
\
'5-10-Methenyl-THF + NADPH -> 5-10-Methylene-THF + NADPplus', \
\
'5-10-Methylene-THF + NH3 + CO2 + NADH + Hplus -> Glycine + THF + NADplus', \
\
'5-10-Methylene-THF + Glycine + H2O -> THF + L-Serine', \
\
'L-Serine -> Pyruvate + NH3', \
\
'2*Fd1ox + Pyruvate + CoA -> 2*Fd1red + Acetyl-CoA + CO2 + 2*Hplus', \
\
'2*Acetyl-CoA -> CoA + Acetoacetyl-CoA', \
\
'Acetoacetyl-CoA + NADPH + Hplus -> Hydroxybutanoyl-CoA + NADPplus', \
\
'Hydroxybutanoyl-CoA -> Crotonoyl-CoA + H2O', \
\
'Crotonoyl-CoA + NADH + Hplus -> Butanoyl-CoA + NADplus', \
\
'Butanoyl-CoA + NADH + Hplus -> Butanal + CoA + NADplus', \
\
'Butanal + NADH + Hplus -> 1-Butanol + NADplus', \
\
'NADH + NADPplus -> NADPH + NADplus', \
\
'2*Fd1red + Hplus + NADplus -> 2*Fd1ox + NADH' \
]









reactantsToGet = ['ATP', 'NADH', 'Fdred', 'Sulfate', 'N2', 'CO2', 'HCO3-', 'HCO2-', 'CO','Fd1red'] 

uniqueCompounds, reactions, sMatrixTKeyIndexed = ParseReactionList(reactionList, reactionArrow='->')

uniqueCompoundsIOStatusProposed = GenerateIOStatusList(uniqueCompounds)

 
uniqueCompoundsIOStatusEdit = \
[['1-Butanol', 'Target'],
['10-Formyl-THF', 'Intermediate'],
['5-10-Methenyl-THF', 'Intermediate'],
['5-10-Methylene-THF', 'Intermediate'],
['ADP', 'Input/Output'],
['ATP', 'Input'],
['Acetoacetyl-CoA', 'Intermediate'],
['Acetyl-CoA', 'Intermediate'],
['Butanal', 'Intermediate'],
['Butanoyl-CoA', 'Intermediate'],
['CO2', 'Input/Output'],
['CoA', 'Input/Output'],
['Crotonoyl-CoA', 'Intermediate'],
['Fd1ox', 'Intermediate'],
['Fd1red', 'Intermediate'],
['Glycine', 'Intermediate'],
['H2O', 'Input/Output'],
['HCO2-', 'Input'],
['Hplus', 'Input'],
['Hydroxybutanoyl-CoA', 'Intermediate'],
['L-Serine', 'Intermediate'],
['NADH', 'Input'],
['NADPH', 'Intermediate'],
['NADPplus', 'Intermediate'],
['NADplus', 'Input/Output'],
['NH3', 'Intermediate'],
['PO43-', 'Input/Output'],
['Pyruvate', 'Intermediate'],
['THF', 'Input/Output']]

 


i = 0
ioStatusEdit = []
while i < len(uniqueCompoundsIOStatusEdit):
    ioStatusEdit.append(uniqueCompoundsIOStatusEdit[i][1])
    i += 1
    
sMatrixT = ConvertIndexedSMatrix(sMatrixTKeyIndexed, uniqueCompounds)
sMatrix = sMatrixT.transpose()


[fVectorOpt, cDotVectorOpt, cDotVectorOptNorm, result] = \
SolveFluxBalanceEquation(sMatrix, reactions, uniqueCompounds, ioStatusEdit)

i = 0
printEverything = True
while i < len(uniqueCompounds):
	ioStatus = ioStatusEdit[i]
	compound = uniqueCompounds[i]
	if printEverything == False and (compound in reactantsToGet) or ioStatus == 'Target':
		print(uniqueCompounds[i] + ':\t' + ioStatus + '\t' + "%.1f" % cDotVectorOptNorm[i])
	elif printEverything == True:
		print(uniqueCompounds[i] + ':\t' + ioStatus + '\t' + "%.1f" % cDotVectorOptNorm[i])
	i += 1

print()
print('Important Compounds')
i = 0
printEverything = False
while i < len(uniqueCompounds):
	ioStatus = ioStatusEdit[i]
	compound = uniqueCompounds[i]
	if printEverything == False and (compound in reactantsToGet) or ioStatus == 'Target':
		print(uniqueCompounds[i] + ':\t' + ioStatus + '\t' + "%.1f" % cDotVectorOptNorm[i])
	elif printEverything == True:
		print(uniqueCompounds[i] + ':\t' + ioStatus + '\t' + "%.1f" % cDotVectorOptNorm[i])
	i += 1
