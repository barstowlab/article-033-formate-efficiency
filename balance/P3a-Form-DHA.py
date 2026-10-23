import re
import os.path
import pdb

from numpy import unique, int8, array

from utils.balanceUtils import SolveFluxBalanceEquation, ConvertIndexedSMatrix, ImportIOStatus, \
ImportReactionFile, GenerateIndexedSMatrixT, PrintStoichiometry, GenerateMergedIOStatusList, \
Get_IO_Status_for_Compound_List, ExportUniqueCompoundsWithIOStatus, GenerateIOStatusList, ParseReactionList

from utils.vectorOutput4 import generateOutputMatrixWithHeaders, writeOutputMatrix, Write_SMatrix
from utils.specutils12 import ensure_dir

# Formolase pathway, dihydroxyacetone variant, version a (Dataset S1 sheet C, reactions 2-13).
# Compound names match Dataset S1 exactly. The reaction arrow is written '->'
# and stoichiometric coefficients '2*X' because that is what ParseReactionList
# expects; the published tables use '→' and '2×X' for the same reactions.
# Reactions after the pathway are the acetyl-CoA to butanol module
# (Dataset S1 sheet G, reactions 1-6).

reactionList = \
[\
'CoA + HCO₂⁻ + ATP -> Formyl-CoA + AMP + HP₂O₇³⁻', \
'Formyl-CoA + NADH + H⁺ -> Formaldehyde + NAD⁺ + CoA', \
'3*Formaldehyde -> Dihydroxyacetone', \
'Dihydroxyacetone + ATP -> Dihydroxyacetone-P + ADP + H⁺', \
'Dihydroxyacetone-P -> D-Glyceraldehyde-3P', \
'D-Glyceraldehyde-3P + NADP⁺ + H₂O -> D-Glycerate-3P + NADPH + 2*H⁺', \
'D-Glycerate-3P -> D-Glycerate-2P', \
'D-Glycerate-2P -> Phosphoenolpyruvate + H₂O', \
'ADP + Phosphoenolpyruvate + H⁺ -> ATP + Pyruvate', \
'Pyruvate + CoA + NAD⁺ -> Acetyl-CoA + CO₂ + NADH', \
'ATP + AMP -> 2*ADP', \
'NADH + NADP⁺ -> NADPH + NAD⁺', \
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
 ['AMP', 'Intermediate'],
 ['ATP', 'Input/Output'],
 ['Acetoacetyl-CoA', 'Intermediate'],
 ['Acetyl-CoA', 'Intermediate'],
 ['Butanal', 'Intermediate'],
 ['Butanoyl-CoA', 'Intermediate'],
 ['CO₂', 'Input/Output'],
 ['CoA', 'Intermediate'],
 ['Crotonoyl-CoA', 'Intermediate'],
 ['D-Glyceraldehyde-3P', 'Intermediate'],
 ['D-Glycerate-2P', 'Intermediate'],
 ['D-Glycerate-3P', 'Intermediate'],
 ['Dihydroxyacetone', 'Intermediate'],
 ['Dihydroxyacetone-P', 'Intermediate'],
 ['Formaldehyde', 'Intermediate'],
 ['Formyl-CoA', 'Intermediate'],
 ['HCO₂⁻', 'Input/Output'],
 ['HP₂O₇³⁻', 'Input/Output'],
 ['H⁺', 'Input/Output'],
 ['H₂O', 'Input/Output'],
 ['NADH', 'Input/Output'],
 ['NADPH', 'Intermediate'],
 ['NADP⁺', 'Intermediate'],
 ['NAD⁺', 'Input/Output'],
 ['Phosphoenolpyruvate', 'Intermediate'],
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
