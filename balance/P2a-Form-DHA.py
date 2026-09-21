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
'CoA + HCO2- + ATP -> Formyl-CoA + AMP + P2O7', \
\
'Formyl-CoA + NADH + Hplus -> Formaldehyde + NADplus + CoA', \
\
'3*Formaldehyde -> DHA', \
\
'DHA + ATP -> DHAP + ADP', \
\
'DHAP -> GA3P', \
\
'GA3P + NADPplus + H2O -> Glycerate3P + NADPH + Hplus', \
\
'Glycerate3P -> Glycerate2P', \
\
'Glycerate2P -> PPyruvate + H2O', \
\
'ADP + PPyruvate -> ATP + Pyruvate', \
\
'2*Fd1ox + Pyruvate + CoA -> 2*Fd1red + Acetyl-CoA + CO2 + 2*Hplus', \
\
'2*Fd1ox + NADH -> 2*Fd1red + Hplus + NADplus', \
\
'2*Acetyl-CoA -> CoA + Acetoacetyl-CoA', \
\
'Acetoacetyl-CoA + NADH + Hplus -> Hydroxybutanoyl-CoA + NADplus', \
\
'Hydroxybutanoyl-CoA -> Crotonoyl-CoA + H2O', \
\
'Crotonoyl-CoA + NADH + Hplus -> Butanoyl-CoA + NADplus', \
\
'Butanoyl-CoA + NADH + Hplus -> Butanal + CoA + NADplus', \
\
'Butanal + NADH + Hplus -> 1-Butanol + NADplus', \
\
'ATP + AMP -> 2*ADP', \
\
'NADH + NADPplus -> NADPH + NADplus' \
]









reactantsToGet = ['ATP', 'NADH', 'Fdred', 'Sulfate', 'N2', 'CO2', 'HCO3-', 'HCO2-', 'CO','Fd1red'] 

uniqueCompounds, reactions, sMatrixTKeyIndexed = ParseReactionList(reactionList, reactionArrow='->')

uniqueCompoundsIOStatusProposed = GenerateIOStatusList(uniqueCompounds)

# Changing IO status to EXACTLY match other code
uniqueCompoundsIOStatusEdit = \
[['1-Butanol', 'Target'],
 ['ADP', 'Input/Output'],
 ['AMP', 'Intermediate'],
 ['ATP', 'Input/Output'],
 ['Acetoacetyl-CoA', 'Intermediate'],
 ['Acetyl-CoA', 'Intermediate'],
 ['Butanal', 'Intermediate'],
 ['Butanoyl-CoA', 'Intermediate'],
 ['CO2', 'Input/Output'],
 ['CoA', 'Intermediate'],
 ['Crotonoyl-CoA', 'Intermediate'],
 ['DHA', 'Intermediate'],
 ['DHAP', 'Intermediate'],
 ['Fd1ox', 'Intermediate'],
 ['Fd1red', 'Intermediate'],
 ['Formaldehyde', 'Intermediate'],
 ['Formyl-CoA', 'Intermediate'],
 ['GA3P', 'Intermediate'],
 ['Glycerate2P', 'Intermediate'],
 ['Glycerate3P', 'Intermediate'],
 ['H2O', 'Input/Output'],
 ['HCO2-', 'Input/Output'],
 ['Hplus', 'Input/Output'],
 ['Hydroxybutanoyl-CoA', 'Intermediate'],
 ['NADH', 'Input/Output'],
 ['NADPH', 'Intermediate'],
 ['NADPplus', 'Intermediate'],
 ['NADplus', 'Input/Output'],
 ['P2O7', 'Input/Output'],
 ['PPyruvate', 'Intermediate'],
 ['Pyruvate', 'Intermediate']]

 


 


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
