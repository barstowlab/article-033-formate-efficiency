#!/opt/local/bin/python3.7

# ------------------------------------------------------------------------------------------------ #
# fig-formate-to-butanol.py
# Calculates efficiency of electricity and CO2 to formate and then to butanol. 
# 
# Buz Barstow
# Last updated by Buz Barstow on 2026-01-28
# ------------------------------------------------------------------------------------------------ #

from rewiredcarbon.scenario import ImportScenarioTable, CalculateScenarioEfficiencies, \
Plot_Efficiency_Bargraph, Generate_EfficienciesDict_Keys_Sorted_by_Efficiency, \
Export_Efficiency_Bargraph

from rewiredcarbon.utils import ensure_dir

from os.path import join

# ------------------------------------------------------------------------------------------------ #
# Get input and outputs all set up

scenarioTableFileName = 'input/fig-efficiency-sota.csv'
outputFilenameEff = 'output/fig-efficiency-sota.csv'
outputFilenameFuelMassEff = 'output/fig-efficiency-sota-mass-eff.csv'

ensure_dir(outputFilenameEff)
ensure_dir(outputFilenameFuelMassEff)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Do the calculation 

# Import the input data
scenarioDict = ImportScenarioTable(scenarioTableFileName)

# Calculate efficiencies
efficienciesDict = CalculateScenarioEfficiencies(scenarioDict)
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Plot the results

keysArray = list(efficienciesDict.keys())

# Plot the energy consumption per unit mass
# Note, in the paper, we plot out energy consumption per mole. We do this conversion in a 
# Excel spreadsheet, which is included in this repository, We will get this code to output in terms
# of energy per mole in a later version. 
Plot_Efficiency_Bargraph(efficienciesDict, 'effTotalElectricalFuelMassEfficiency', \
'effTotalElectricalFuelMassEfficiency_lowerError', \
'effTotalElectricalFuelMassEfficiency_upperError', keysToPlot=keysArray)


# Output fuel mass efficiency in (grams of product per joule of input energy). 
Export_Efficiency_Bargraph(outputFilenameFuelMassEff, efficienciesDict, scenarioDict, \
'effTotalElectricalFuelMassEfficiency', 'effTotalElectricalFuelMassEfficiency_lowerError', \
'effTotalElectricalFuelMassEfficiency_upperError', keysToPlot=keysArray)


# Plot the energy conversion efficiency
Plot_Efficiency_Bargraph(efficienciesDict, 'effTotalElectricalToFuel', \
'effTotalElectricalToFuel_lowerError', 'effTotalElectricalToFuel_upperError', keysToPlot=keysArray)

Export_Efficiency_Bargraph(outputFilenameEff, efficienciesDict, scenarioDict, \
'effTotalElectricalToFuel', 'effTotalElectricalToFuel_lowerError', \
'effTotalElectricalToFuel_upperError', keysToPlot=keysArray)

# ------------------------------------------------------------------------------------------------ #