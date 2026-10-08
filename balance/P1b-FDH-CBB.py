import re
import os.path
import pdb

from numpy import unique, int8, array

from utils.balanceUtils import SolveFluxBalanceEquation, ConvertIndexedSMatrix, ImportIOStatus, \
ImportReactionFile, GenerateIndexedSMatrixT, PrintStoichiometry, GenerateMergedIOStatusList, \
Get_IO_Status_for_Compound_List, ExportUniqueCompoundsWithIOStatus, GenerateIOStatusList, ParseReactionList

from utils.vectorOutput4 import generateOutputMatrixWithHeaders, writeOutputMatrix, Write_SMatrix
from utils.specutils12 import ensure_dir

# Calvin-Benson-Bassham cycle, version b with formate dehydrogenase (Dataset S1 sheet A, reactions 1-19).
# Compound names match Dataset S1 exactly. The reaction arrow is written '->'
# and stoichiometric coefficients '2*X' because that is what ParseReactionList
# expects; the published tables use '→' and '2×X' for the same reactions.
# Reactions after the pathway are the acetyl-CoA to butanol module
# (Dataset S1 sheet G, reactions 1-6).

reactionList = \
[\
'CO₂ + D-Ribulose-1,5-bisphosphate + H₂O -> 2*D-Glycerate-3P + 2*H⁺', \
'D-Glycerate-3P + ATP -> D-Glycerate-1,3-bisphosphate + ADP', \
'D-Glycerate-1,3-bisphosphate + NADPH + H⁺ -> D-Glyceraldehyde-3P + HPO₄²⁻ + NADP⁺', \
'D-Glyceraldehyde-3P + Dihydroxyacetone-P -> D-Fructose-1,6-bisphosphate', \
'D-Fructose-1,6-bisphosphate + H₂O -> D-Fructose-6P + HPO₄²⁻', \
'D-Fructose-6P + D-Glyceraldehyde-3P -> D-Erythrose-4P + D-Xylulose-5P', \
'D-Glyceraldehyde-3P -> Dihydroxyacetone-P', \
'D-Erythrose-4P + Dihydroxyacetone-P -> Sedoheptulose-1,7-bisphosphate', \
'D-Xylulose-5P -> D-Ribulose-5P', \
'Sedoheptulose-1,7-bisphosphate + H₂O -> Sedoheptulose-7P + HPO₄²⁻', \
'Sedoheptulose-7P + D-Glyceraldehyde-3P -> D-Ribose-5P + D-Xylulose-5P', \
'D-Ribose-5P -> D-Ribulose-5P', \
'ATP + D-Ribulose-5P -> ADP + D-Ribulose-1,5-bisphosphate + H⁺', \
'D-Glycerate-3P -> D-Glycerate-2P', \
'D-Glycerate-2P -> Phosphoenolpyruvate + H₂O', \
'ADP + Phosphoenolpyruvate + H⁺ -> ATP + Pyruvate', \
'Pyruvate + CoA + NAD⁺ -> Acetyl-CoA + CO₂ + NADH', \
'NADH + NADP⁺ -> NADPH + NAD⁺', \
'HCO₂⁻ + NAD⁺ -> CO₂ + NADH', \
'2*Acetyl-CoA -> CoA + Acetoacetyl-CoA', \
'Acetoacetyl-CoA + NADH + H⁺ -> (S)-3-Hydroxybutanoyl-CoA + NAD⁺', \
'(S)-3-Hydroxybutanoyl-CoA -> Crotonoyl-CoA + H₂O', \
'Crotonoyl-CoA + NADPH + H⁺ -> Butanoyl-CoA + NADP⁺', \
'Butanoyl-CoA + NADH + H⁺ -> Butanal + CoA + NAD⁺', \
'Butanal + NADH + H⁺ -> 1-Butanol + NAD⁺' \
]

reactantsToGet = ['ATP', 'NADH', 'Fdred', 'Sulfate', 'N2', 'CO₂', 'HCO₃⁻', 'HCO₂⁻', \
'CO', 'Fd1red']

uniqueCompounds, reactions, sMatrixTKeyIndexed = ParseReactionList(reactionList, reactionArrow='->')

uniqueCompoundsIOStatusProposed = GenerateIOStatusList(uniqueCompounds)

# Changing IO status to EXACTLY match other code
uniqueCompoundsIOStatusEdit = \
[['(S)-3-Hydroxybutanoyl-CoA', 'Intermediate'],
 ['1-Butanol', 'Target'],
 ['ADP', 'Input/Output'],
 ['ATP', 'Input/Output'],
 ['Acetoacetyl-CoA', 'Intermediate'],
 ['Acetyl-CoA', 'Intermediate'],
 ['Butanal', 'Intermediate'],
 ['Butanoyl-CoA', 'Intermediate'],
 ['CO₂', 'Intermediate'],
 ['CoA', 'Intermediate'],
 ['Crotonoyl-CoA', 'Intermediate'],
 ['D-Erythrose-4P', 'Intermediate'],
 ['D-Fructose-1,6-bisphosphate', 'Intermediate'],
 ['D-Fructose-6P', 'Intermediate'],
 ['D-Glyceraldehyde-3P', 'Intermediate'],
 ['D-Glycerate-1,3-bisphosphate', 'Intermediate'],
 ['D-Glycerate-2P', 'Intermediate'],
 ['D-Glycerate-3P', 'Intermediate'],
 ['D-Ribose-5P', 'Intermediate'],
 ['D-Ribulose-1,5-bisphosphate', 'Intermediate'],
 ['D-Ribulose-5P', 'Intermediate'],
 ['D-Xylulose-5P', 'Intermediate'],
 ['Dihydroxyacetone-P', 'Intermediate'],
 ['HCO₂⁻', 'Input/Output'],
 ['HPO₄²⁻', 'Input/Output'],
 ['H⁺', 'Input/Output'],
 ['H₂O', 'Input/Output'],
 ['NADH', 'Input/Output'],
 ['NADPH', 'Intermediate'],
 ['NADP⁺', 'Intermediate'],
 ['NAD⁺', 'Input/Output'],
 ['Phosphoenolpyruvate', 'Intermediate'],
 ['Pyruvate', 'Intermediate'],
 ['Sedoheptulose-1,7-bisphosphate', 'Intermediate'],
 ['Sedoheptulose-7P', 'Intermediate']]


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
