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

# In this version, I've added back the formate dehydrogenase reaction. 

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
'Glycine + PO43- + TrxRed -> Acetyl-phosphate + NH3 + TrxOx + H2O', \
\
'ADP + Acetyl-phosphate -> ATP + Acetate', \
\
'ATP + Acetate + CoA -> AMP + P2O7 + Acetyl-CoA', \
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
'ATP + AMP -> 2*ADP', \
\
'HCO2- + NADplus -> CO2 + NADH + Hplus', \
\
'NADH + NADPplus -> NADPH + NADplus', \
\
'TrxOx + NADPH -> TrxRed + NADPplus + Hplus' \
]






reactantsToGet = ['ATP', 'NADH', 'Fdred', 'Sulfate', 'N2', 'CO2', 'HCO3-', 'HCO2-', 'CO','Fd1red'] 

uniqueCompounds, reactions, sMatrixTKeyIndexed = ParseReactionList(reactionList, reactionArrow='->')

uniqueCompoundsIOStatusProposed = GenerateIOStatusList(uniqueCompounds)

 
uniqueCompoundsIOStatusEdit = \
[['1-Butanol', 'Target'],
['10-Formyl-THF', 'Intermediate'],
['5-10-Methenyl-THF', 'Intermediate'],
['5-10-Methylene-THF', 'Intermediate'],
['ADP', 'Output'],
['AMP', 'Intermediate'],
['ATP', 'Input'],
['Acetate', 'Intermediate'],
['Acetoacetyl-CoA', 'Intermediate'],
['Acetyl-CoA', 'Intermediate'],
['Acetyl-phosphate', 'Intermediate'],
['Butanal', 'Intermediate'],
['Butanoyl-CoA', 'Intermediate'],
['CO2', 'Intermediate'],
['CoA', 'Output'],
['Crotonoyl-CoA', 'Intermediate'],
['Glycine', 'Intermediate'],
['H2O', 'Input/Output'],
['HCO2-', 'Input/Output'],
['Hplus', 'Input/Output'],
['Hydroxybutanoyl-CoA', 'Intermediate'],
['NADH', 'Input'],
['NADPH', 'Intermediate'],
['NADPplus', 'Intermediate'],
['NADplus', 'Input/Output'],
['NH3', 'Intermediate'],
['P2O7', 'Output'],
['PO43-', 'Output'],
['THF', 'Output'],
['TrxOx', 'Intermediate'],
['TrxRed', 'Intermediate']
]

 


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
