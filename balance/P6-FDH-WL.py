import re
import os.path
import pdb

from numpy import unique, int8, array

from utils.balanceUtils import SolveFluxBalanceEquation, ConvertIndexedSMatrix, ImportIOStatus, \
ImportReactionFile, GenerateIndexedSMatrixT, PrintStoichiometry, GenerateMergedIOStatusList, \
Get_IO_Status_for_Compound_List, ExportUniqueCompoundsWithIOStatus, GenerateIOStatusList, ParseReactionList

from utils.vectorOutput4 import generateOutputMatrixWithHeaders, writeOutputMatrix, Write_SMatrix
from utils.specutils12 import ensure_dir

# Wood-Ljungdahl pathway with formate dehydrogenase (Dataset S1 sheet F, reactions 1-9).
# Compound names match Dataset S1 exactly. The reaction arrow is written '->'
# and stoichiometric coefficients '2*X' because that is what ParseReactionList
# expects; the published tables use '→' and '2×X' for the same reactions.
# Reactions after the pathway are the acetyl-CoA to butanol module
# (Dataset S1 sheet G, reactions 1-6).

reactionList = \
[\
'HCO₂⁻ + NAD⁺ -> CO₂ + NADH', \
'THF + HCO₂⁻ + ATP -> ADP + HPO₄²⁻ + 10-Formyl-THF', \
'10-Formyl-THF + H⁺ -> 5,10-Methenyl-THF + H₂O', \
'5,10-Methenyl-THF + NADPH -> 5,10-Methylene-THF + NADP⁺', \
'5,10-Methylene-THF + NADH + H⁺ -> 5-Methyl-THF + NAD⁺', \
'CO₂ + 2*Fd1red + 2*H⁺ -> CO + H₂O + 2*Fd1ox', \
'5-Methyl-THF + CoA + CO -> Acetyl-CoA + THF', \
'NADH + NADP⁺ -> NADPH + NAD⁺', \
'2*Fd1ox + NADH -> 2*Fd1red + H⁺ + NAD⁺', \
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
 ['10-Formyl-THF', 'Intermediate'],
 ['5,10-Methenyl-THF', 'Intermediate'],
 ['5,10-Methylene-THF', 'Intermediate'],
 ['5-Methyl-THF', 'Intermediate'],
 ['ADP', 'Input/Output'],
 ['ATP', 'Input/Output'],
 ['Acetoacetyl-CoA', 'Intermediate'],
 ['Acetyl-CoA', 'Intermediate'],
 ['Butanal', 'Intermediate'],
 ['Butanoyl-CoA', 'Intermediate'],
 ['CO', 'Intermediate'],
 ['CO₂', 'Intermediate'],
 ['CoA', 'Intermediate'],
 ['Crotonoyl-CoA', 'Intermediate'],
 ['Fd1ox', 'Intermediate'],
 ['Fd1red', 'Intermediate'],
 ['HCO₂⁻', 'Input/Output'],
 ['HPO₄²⁻', 'Input/Output'],
 ['H⁺', 'Input/Output'],
 ['H₂O', 'Input/Output'],
 ['NADH', 'Input/Output'],
 ['NADPH', 'Intermediate'],
 ['NADP⁺', 'Intermediate'],
 ['NAD⁺', 'Input/Output'],
 ['THF', 'Intermediate']]


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
