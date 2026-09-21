import re
import os.path
import pdb

from numpy import unique, int8, array

from utils.balanceUtils import SolveFluxBalanceEquation, ConvertIndexedSMatrix, ImportIOStatus, \
ImportReactionFile, GenerateIndexedSMatrixT, PrintStoichiometry, GenerateMergedIOStatusList, \
Get_IO_Status_for_Compound_List, ExportUniqueCompoundsWithIOStatus, GenerateIOStatusList, ParseReactionList

from utils.vectorOutput4 import generateOutputMatrixWithHeaders, writeOutputMatrix, Write_SMatrix
from utils.specutils12 import ensure_dir

reactionList = \
[\
'CO2 + Ribulose15BP + H2O -> 2*3PG', \
'3PG + ATP -> 13BPG + ADP', \
'13BPG + NADH + Hplus -> G3P + P + NADplus', \
'G3P + GlyceronePhosph -> F16BP', \
'F16BP + H2O -> F6P + P', \
'F6P + G3P -> E4P + X5P', \
'G3P -> GlyceronePhosph', \
'E4P + GlyceronePhosph -> SH17BP', \
'X5P -> Ribulose5P', \
'SH17BP + H2O -> SH7P + P', \
'SH7P + G3P -> Ribose5P + X5P', \
'Ribose5P -> Ribulose5P', \
'ATP + Ribulose5P -> ADP + Ribulose15BP', \
'3PG -> 2PG', \
'2PG -> PEP + H2O', \
'ADP + PEP -> Pyruvate + ATP', \
'2*Fd1ox + Pyruvate + CoA -> 2*Fd1red + Acetyl-CoA + CO2 + 2*Hplus', \
'NADH + NADPplus -> NADPH + NADplus', \
'2*Fd1red + Hplus + NADplus -> 2*Fd1ox + NADH', \
'2*Acetyl-CoA -> CoA + Acetoacetyl-CoA', \
'Acetoacetyl-CoA + NADPH + Hplus -> Hydroxybutanoyl-CoA + NADPplus', \
'Hydroxybutanoyl-CoA -> Crotonoyl-CoA + H2O', \
'Crotonoyl-CoA + NADH + Hplus -> Butanoyl-CoA + NADplus', \
'Butanoyl-CoA + NADH + Hplus -> Butanal + CoA + NADplus', \
'Butanal + NADH + Hplus -> 1-Butanol + NADplus' \
]

reactantsToGet = ['ATP', 'NADH', 'Fdred', 'Sulfate', 'N2', 'CO2', 'HCO3-', 'HCO2-', \
'CO','Fd1red'] 

uniqueCompounds, reactions, sMatrixTKeyIndexed = ParseReactionList(reactionList, reactionArrow='->')

uniqueCompoundsIOStatusProposed = GenerateIOStatusList(uniqueCompounds)
uniqueCompoundsIOStatusProposed

uniqueCompoundsIOStatusEdit = \
[['1-Butanol', 'Target'],
 ['13BPG', 'Intermediate'],
 ['2PG', 'Intermediate'],
 ['3PG', 'Intermediate'],
 ['ADP', 'Input/Output'],
 ['ATP', 'Input/Output'],
 ['Acetoacetyl-CoA', 'Intermediate'],
 ['Acetyl-CoA', 'Intermediate'],
 ['Butanal', 'Intermediate'],
 ['Butanoyl-CoA', 'Intermediate'],
 ['CO2', 'Input/Output'],
 ['CoA', 'Intermediate'],
 ['Crotonoyl-CoA', 'Intermediate'],
 ['E4P', 'Intermediate'],
 ['F16BP', 'Intermediate'],
 ['F6P', 'Intermediate'],
 ['Fd1ox', 'Intermediate'],
 ['Fd1red', 'Intermediate'],
 ['G3P', 'Intermediate'],
 ['GlyceronePhosph', 'Intermediate'],
 ['H2O', 'Input/Output'],
 ['Hplus', 'Input/Output'],
 ['Hydroxybutanoyl-CoA', 'Intermediate'],
 ['NADH', 'Input/Output'],
 ['NADPH', 'Intermediate'],
 ['NADPplus', 'Intermediate'],
 ['NADplus', 'Input/Output'],
 ['P', 'Input/Output'],
 ['PEP', 'Intermediate'],
 ['Pyruvate', 'Intermediate'],
 ['Ribose5P', 'Intermediate'],
 ['Ribulose15BP', 'Intermediate'],
 ['Ribulose5P', 'Intermediate'],
 ['SH17BP', 'Intermediate'],
 ['SH7P', 'Intermediate'],
 ['X5P', 'Intermediate']]


i = 0
ioStatusEdit = []
while i < len(uniqueCompoundsIOStatusEdit):
    ioStatusEdit.append(uniqueCompoundsIOStatusEdit[i][1])
    i += 1

sMatrixT = ConvertIndexedSMatrix(sMatrixTKeyIndexed, uniqueCompounds)
sMatrix = sMatrixT.transpose()


[fVectorOpt, cDotVectorOpt, cDotVectorOptNorm, result] = \
SolveFluxBalanceEquation(sMatrix, reactions, uniqueCompounds, ioStatusEdit)

[fVectorOpt, cDotVectorOpt, cDotVectorOptNorm, result]

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
