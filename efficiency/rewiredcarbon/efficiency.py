# ------------------------------------------------------------------------------------------------ #
# efficiency.py
# Buz Barstow
# Created: 2017-10-26
# Last Update: 2018-06-27 by Farshid Salimijazi
# Last Update: 2018-07-22 by Buz Barstow
# Last Update: 2018-07-26 by Farshid Salimijazi
# Update: 2019-06-01 by Buz Barstow
# Update: 2019-06-15 by Buz Barstow
# Update: 2020-04-17 by Buz Barstow
# Update: 2020-06-02 by Farshid Salimijazi
# Update: 2025-01-25 by Buz Barstow
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# On 2025-01-02 I (Buz) seriously simplified this code. I feel like it was pretty bloated, and did
# a lot of things that I wasn't using (relating to scale up by stirring), which made upgrading
# the electrochemical parts of the code difficult. \
# 
# Earlier in 2024 Buz made a minor correction to the electrochemical CO2 calculations as well. 
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# For symbol definitions see the file SymbolTable.pdf, included with this repository.
# Note that LaTeX macro names, are the same as symbols used here, with some minor mismatches. 
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Variables used in this file

# Input parameters
# vRedoxCellOne - Potential difference across first electrochemical cell in Volts
# vBiasCellOne - Additional external bias voltage across cell one in Volts
# vRedoxCellTwo - Potential difference across second electrochemical cell (bio-cell) in Volts
# vBiasCellTwo - Additional external bias voltage across cell 2 in Volts

# vCellTwo - Potential across second electrochemical cell in Volts

# efficiencyCurrentToFirstCell - Efficiency of conversion of current to first cell
# efficiencyCurrentToSecondCell - Efficiency of conversion of current to second cell
# efficiencyCarbonToSecondCell - Efficiency of conversion of fixed carbon from first to second cell
# Assumed to be 1.0.

# carbonsPerPrimaryFix - Number of carbon atoms per primary fixation


# electronsPerPrimaryFix - Number of electrons needed per carbon to reduce CO2 to 
# primary fixation product in first cell

# electronsNADHforFuel - Number of electrons needed to make NADH to reduce CO2 to fuel

# energyPerFuelMolecule - Energy extracted when fuel molecule is combusted in standard conditions
# in Joules

# carbonsPerFuel - Number of carbon atoms per fuel

# efficiencyElectronsToHyd - Efficiency of conversion of current to H2. Assumed to be 1.0. 

# electricalPower - Scaling power number for calculation in Watts. Efficiency should not change
# if this is changed unless you get close to a very small number (electronsPerFuelInH) number of 
# electrons. This is around 24 to 40 electrons. 

# vRedox - Redox potential difference of electrochemical cell in Volts

# vBias - Additional external bias voltage needed to make H2 in Volts

# vMembrane - Potential across inner membrane in microbe in Volts

# vAcceptor - Potential of terminal electron acceptor in Volts vs. Standard Hydrogen Electrode

# vMtr - Potential of Mtr EEU complex in Volts vs. Standard Hydrogen Electrode

# vQuinone - Potential of quinone molecule in Volts vs. Standard Hydrogen Electrode

# NADHforFuel - Number of NADH to reduce CO2 to fuel

# FdForFuel - Number of Ferredoxin to reduce CO2 to fuel

# ATPforFuel - Number of ATP to make fuel from CO2
# energyPerFuelMolecule - Energy extracted when fuel molecule is combusted in standard conditions
# in Joules

# efficiencyElectronsToMtr - Efficiency of conversion of current to reduced Mtr. Assumed to be 1.0.

# efficiencyElectronsMtrToQuinone - Efficiency of conversion of current to reduced Mtr. 
# Assumed to be 1.0.

# numberOfProtonsPumpedInForATP - Number of protons needed to be pumped through ATP synthase to 
# catalyze ADP + Pi -> ATP. When set to -1, we calculate based on the membrane potential and ...

# deltaGADPATP - The free energy difference between ATP and ADP

# cellDensity - Density of planktonic cells in electrochemical cell (usually used for H2-mediated
# cases. Unit is cells per m^3. 
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def Energy_Per_Fuel_Molecule(specificEnergyFuelIn, molecularWeightFuelIn):
# Calculates the combustion energy stored in each fuel molecule
# specificEnergyFuelIn - Energy density of fuel in Megajoules per Kg
# molecularWeightFuelIn - Molecular weight of compound in grams per mole


	from rewiredcarbon.physicalconstants import avogadroNumber as N_a
	
	energyPerMol = (specificEnergyFuelIn * 1000) * molecularWeightFuelIn
	energyPerFuelMolecule = energyPerMol/N_a

	return energyPerFuelMolecule
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def ProtonPumping_and_Electron_Requirements_for_EEU_to_Fuel(vMembrane, vAcceptor, vMtr, vQuinone, \
deltaGADPATP, NADHforFuel, FdForFuel, ATPforFuel, numberOfProtonsPumpedInForATP=-1, \
maxAllowedProtonsPumpedOutPerElectronDown=-1):

	from math import floor, ceil
	from rewiredcarbon.physicalconstants import elementaryCharge as e
	from rewiredcarbon.electrochemicalconstants import vNADH, vFd
	import pdb
	
	# -------------------------------------------------------------------------------------------- #
	# ATP
	# Calculate number of protons needed to be pumped in to generate 1 ATP, and number of electrons
	# needed to be sent downhill to generate proton gradient
	# See Farshid #1 pages 29 and 30

	if numberOfProtonsPumpedInForATP == -1:
		# If the number of protons pumped into make a single ATP is set to -1, we calculate it
		# from the membrane potential
		numberOfProtonsPumpedInForATP = ceil(abs(deltaGADPATP/(e*vMembrane)))
		# print("Number of protons pumped in for ATP: " + str(numberOfProtonsPumpedInForATP))
# 		print("Transmembrane voltage: " + str(vMembrane))
		# Otherwise, you can override this by setting numberOfProtonsPumpedInForATP to a positive
		# value
	
	# Calculate the maximum possible number of protons that could be pumped out by sending an
	# electron downhill to terminal electron acceptor. 
	maxPossibleProtonsPumpedOutPerElectronDown = floor(abs((vQuinone - vAcceptor) / vMembrane))

	if maxAllowedProtonsPumpedOutPerElectronDown == -1:
		numberOfProtonsPumpedOutPerElectronDownInEEU = maxPossibleProtonsPumpedOutPerElectronDown
		# Default behavior
		# Unlimited protons pumped out per electron. 
	
	elif maxPossibleProtonsPumpedOutPerElectronDown > maxAllowedProtonsPumpedOutPerElectronDown:
		numberOfProtonsPumpedOutPerElectronDownInEEU = maxAllowedProtonsPumpedOutPerElectronDown
	
	elif maxPossibleProtonsPumpedOutPerElectronDown <= maxAllowedProtonsPumpedOutPerElectronDown:
		numberOfProtonsPumpedOutPerElectronDownInEEU = maxPossibleProtonsPumpedOutPerElectronDown
		
		
	try:
		electronsDownATPforFuelInEEU = (ATPforFuel * numberOfProtonsPumpedInForATP)\
		/ numberOfProtonsPumpedOutPerElectronDownInEEU
	except:
		pdb.set_trace()
	# -------------------------------------------------------------------------------------------- #

	# -------------------------------------------------------------------------------------------- #
	# Reductants: NADH and Ferredoxin
	# See Farshid #1 pages 29 and 30
	# Calculate proton pumping needed to send electrons uphill to NADH and Ferredoxin, and electrons
	# that need to be sent downhill to generate the necessary proton gradient. 
	
	numberOfProtonsPumpedInPerElectronForNADH = ceil(abs((vNADH - vQuinone)/vMembrane))
	electronsNADHforFuel = 2*NADHforFuel
	
	numberOfProtonsPumpedInPerElectronForFd = ceil(abs((vFd - vQuinone)/vMembrane))
	numberOfProtonsPumpedInPerFd = 2*numberOfProtonsPumpedInPerElectronForFd
	electronsFdforFuel = 2*FdForFuel

	electronsDownNADHforFuel = (electronsNADHforFuel * numberOfProtonsPumpedInPerElectronForNADH) \
	/numberOfProtonsPumpedOutPerElectronDownInEEU
	
	electronsDownFdForFuel = (electronsFdforFuel * numberOfProtonsPumpedInPerElectronForFd) \
	/numberOfProtonsPumpedOutPerElectronDownInEEU
	# -------------------------------------------------------------------------------------------- #

	# -------------------------------------------------------------------------------------------- #
	# Calculate the total number of electrons we need to make one fuel molecule through EEU.
	# Defined in Farshid #1 page 30.
	electronsPerFuelInEEU = electronsNADHforFuel + electronsFdforFuel \
	+ electronsDownATPforFuelInEEU + electronsDownNADHforFuel + electronsDownFdForFuel
	
	# Calculate the number of protons that get pumped out and back in to make the ATP and reductants
	numberOfProtonsPumpedOut = \
	(electronsDownATPforFuelInEEU * numberOfProtonsPumpedOutPerElectronDownInEEU) \
	+ (electronsDownNADHforFuel * numberOfProtonsPumpedOutPerElectronDownInEEU) \
	+ (electronsDownFdForFuel * numberOfProtonsPumpedOutPerElectronDownInEEU)
	# -------------------------------------------------------------------------------------------- #


	return numberOfProtonsPumpedOut, electronsPerFuelInEEU, numberOfProtonsPumpedInForATP
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def ProtonPumping_and_Electron_Requirements_for_Hydrogen_to_Fuel(vMembrane, vAcceptor, \
deltaGADPATP, NADHforFuel, FdForFuel, ATPforFuel, numberOfProtonsPumpedInForATP=-1, \
maxAllowedProtonsPumpedOutPerElectronDown=-1):
	
	from rewiredcarbon.electrochemicalconstants import vH2
	from rewiredcarbon.physicalconstants import elementaryCharge as e
	from math import floor, ceil
	
	# Calculate the number of protons that you get pumped outside the membrane to store energy by 
	# sending one electron downhill.
	# See Farshid #32
	
	# Calculate the number of electrons the bug has to send downhill to the acceptor
	# in order to make all of the ATP it needs for one fuel molecule.

	# Number of ATPs per fuel molecule * number of protons that you pump inside the membrane to 
	# make one ATP / number of protons pumped out by sending one electron downhill to acceptor
	# See Farshid #1 page 30.
	
	if numberOfProtonsPumpedInForATP == -1:
		# If the number of protons pumped into make a single ATP is set to -1, we calculate it
		# from the membrane potential
		numberOfProtonsPumpedInForATP = ceil(abs(deltaGADPATP/(e * vMembrane)))
		# Otherwise, you can override this by setting numberOfProtonsPumpedInForATP to a positive
		# value
	
	# Calculate the maximum possible number of protons that could be pumped out by sending an
	# electron downhill to terminal electron acceptor. 
	maxPossibleProtonsPumpedOutPerElectronDown = floor(abs((vH2 - vAcceptor) / vMembrane))

	if maxAllowedProtonsPumpedOutPerElectronDown == -1:
		numberOfProtonsPumpedOutPerElectronDownInH = maxPossibleProtonsPumpedOutPerElectronDown
		# Default behavior
		# Unlimited protons pumped out per electron. 
	
	elif maxPossibleProtonsPumpedOutPerElectronDown > maxAllowedProtonsPumpedOutPerElectronDown:
		numberOfProtonsPumpedOutPerElectronDownInH = maxAllowedProtonsPumpedOutPerElectronDown
	
	elif maxPossibleProtonsPumpedOutPerElectronDown <= maxAllowedProtonsPumpedOutPerElectronDown:
		numberOfProtonsPumpedOutPerElectronDownInH = maxPossibleProtonsPumpedOutPerElectronDown
		
		
	try:
		electronsDownATPforFuelInH = (ATPforFuel * numberOfProtonsPumpedInForATP)\
		/ numberOfProtonsPumpedOutPerElectronDownInH
	except:
		pdb.set_trace()
	
	electronsPerFuelInH = 2*NADHforFuel + 2*FdForFuel \
	+ ATPforFuel*(numberOfProtonsPumpedInForATP/numberOfProtonsPumpedOutPerElectronDownInH)
	
	electronsDownATPforFuelInH = (ATPforFuel * numberOfProtonsPumpedInForATP)\
	/ numberOfProtonsPumpedOutPerElectronDownInH
	
	numberOfProtonsPumpedOut = \
	electronsDownATPforFuelInH * numberOfProtonsPumpedOutPerElectronDownInH

	
	return numberOfProtonsPumpedOut, electronsPerFuelInH, numberOfProtonsPumpedInForATP
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def ProtonPumping_and_Electron_Requirements_for_LoP_to_Fuel(vMediator, vMembrane, \
vAcceptor, deltaGADPATP, NADHforFuel, FdForFuel, ATPforFuel, numberOfProtonsPumpedInForATP=-1, \
maxAllowedProtonsPumpedOutPerElectronDown=-1):

	from rewiredcarbon.physicalconstants import elementaryCharge as e
	from math import floor, ceil
	
	# Calculate the number of protons that you get pumped outside the membrane to store energy by 
	# sending one electron downhill from a low potential mediator, and the number of electrons the 
	# bug has to send downhill to the acceptor in order to make all of the ATP it needs for one
	# fuel molecule.

	# This is a generalized version of ProtonPumping_and_Electron_Requirements_for_Hydrogen_to_Fuel
	# that applies to low potential (i.e., with a potential lower than NADH) mediators. 

	# Number of ATPs per fuel molecule * number of protons that you pump inside the membrane to 
	# make one ATP / number of protons pumped out by sending one electron downhill to acceptor
	# See Farshid #1 page 30.
	
	# vMediator - The redox potential of the mediator in Volts vs. SHE.

	print('291')
	
	if numberOfProtonsPumpedInForATP == -1:
		# If the number of protons pumped into make a single ATP is set to -1, we calculate it
		# from the membrane potential
		numberOfProtonsPumpedInForATP = ceil(abs(deltaGADPATP/(e * vMembrane)))
		# Otherwise, you can override this by setting numberOfProtonsPumpedInForATP to a positive
		# value
	
	# Calculate the maximum possible number of protons that could be pumped out by sending an
	# electron downhill to terminal electron acceptor. 
	maxPossibleProtonsPumpedOutPerElectronDown = floor(abs((vMediator - vAcceptor) / vMembrane))

	if maxAllowedProtonsPumpedOutPerElectronDown == -1:
		numberOfProtonsPumpedOutPerElectronDownInLoP = maxPossibleProtonsPumpedOutPerElectronDown
		# Default behavior
		# Unlimited protons pumped out per electron. 
	
	elif maxPossibleProtonsPumpedOutPerElectronDown > maxAllowedProtonsPumpedOutPerElectronDown:
		numberOfProtonsPumpedOutPerElectronDownInLoP = maxAllowedProtonsPumpedOutPerElectronDown
	
	elif maxPossibleProtonsPumpedOutPerElectronDown <= maxAllowedProtonsPumpedOutPerElectronDown:
		numberOfProtonsPumpedOutPerElectronDownInLoP = maxPossibleProtonsPumpedOutPerElectronDown
		
		
	try:
		electronsDownATPforFuelInLoP = (ATPforFuel * numberOfProtonsPumpedInForATP)\
		/ numberOfProtonsPumpedOutPerElectronDownInLoP
	except:
		pdb.set_trace()
	
	print('322')
	
	electronsPerFuelInLoP = 2*NADHforFuel + 2*FdForFuel \
	+ ATPforFuel*(numberOfProtonsPumpedInForATP/numberOfProtonsPumpedOutPerElectronDownInLoP)
	
	electronsDownATPforFuelInH = (ATPforFuel * numberOfProtonsPumpedInForATP)\
	/ numberOfProtonsPumpedOutPerElectronDownInLoP
	
	numberOfProtonsPumpedOut = \
	electronsDownATPforFuelInLoP * numberOfProtonsPumpedOutPerElectronDownInLoP

	print('333')

	
	return numberOfProtonsPumpedOut, electronsPerFuelInLoP, numberOfProtonsPumpedInForATP
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def Calculate_Fuel_Mass_Efficiency(molecularWeightFuelMolecule, totalFuelProductionRate, \
availableElectricalPower, totalElectricalPower, totalInputPower):
	
	from rewiredcarbon.physicalconstants import avogadroNumber

	# Calculate the mass production rate of fuel (grams per second)
	massPerFuelMolecule = molecularWeightFuelMolecule/avogadroNumber
	totalFuelMassProductionRate = totalFuelProductionRate*massPerFuelMolecule
	
	available_Electrical_to_Fuel_Mass_Efficiency = \
	availableElectricalPower / totalFuelMassProductionRate
	
	total_Electrical_to_to_Fuel_Mass_Efficiency = \
	totalElectricalPower / totalFuelMassProductionRate

	total_Input_Power_to_to_Fuel_Mass_Efficiency = totalInputPower / totalFuelMassProductionRate

	return available_Electrical_to_Fuel_Mass_Efficiency, \
	total_Electrical_to_to_Fuel_Mass_Efficiency, total_Input_Power_to_to_Fuel_Mass_Efficiency

# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Update_Output_Dict_with_Fuel_Mass_Efficiency(outputDict, molecularWeightFuelMolecule, \
totalFuelProductionRate, availableElectricalPower, totalElectricalPower, totalInputPower):

	import pdb
	
	# Update the output dict with the fuel mass efficiency (grams of product per joule of input
	# power). 
		
	available_Electrical_to_Fuel_Mass_Efficiency, total_Electrical_to_to_Fuel_Mass_Efficiency, \
	total_Input_Power_to_to_Fuel_Mass_Efficiency = \
	Calculate_Fuel_Mass_Efficiency(molecularWeightFuelMolecule, totalFuelProductionRate, \
	availableElectricalPower, totalElectricalPower, totalInputPower)

	outputDict['effAvailElectricalFuelMassEfficiency'] = \
	available_Electrical_to_Fuel_Mass_Efficiency
	outputDict['effTotalElectricalFuelMassEfficiency'] = \
	total_Electrical_to_to_Fuel_Mass_Efficiency
	outputDict['effTotalPowerFuelMassEfficiency'] = \
	total_Input_Power_to_to_Fuel_Mass_Efficiency
	
	return	
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Efficiencies_Current_TotalEnergyContentOfFuel_Hydrogen_BioCO2(availableElectricalPower, \
totalElectricalPower, totalInputPower, voltageCellTwo, efficiencyElectronsToHyd, \
energyPerFuelMolecule, electronsPerFuelInH, molecularWeightFuelMolecule=-1):

	from rewiredcarbon.physicalconstants import elementaryCharge as e
	import pdb
	
	# Calculate the currents in the electrochemical cell. Defined in Buz Medium #1 page 22 and 
	# Farshid #1 page 32.
	# Updated so that cell voltage is just one variable, voltageCellTwo. 
	currentCell = availableElectricalPower / voltageCellTwo
	hydrogenCurrent = efficiencyElectronsToHyd * currentCell

	
	# Calculate the total rate of fuel production (fuel molecules per second)
	# See Farshid #1 page 32 and Buz Medium #1 pages 22-23. 
	totalFuelProductionRateInHyd = hydrogenCurrent / (e * electronsPerFuelInH)
	totalEnergyContentOfFuel = totalFuelProductionRateInHyd * energyPerFuelMolecule
		
	
	# Calculate the H2 to fuel efficiencies. Defined in Farshid #1 page 37.
	
	available_Electrical_to_Hyd_to_Fuel_Efficiency = \
	totalEnergyContentOfFuel / availableElectricalPower

	total_Electrical_to_Hyd_to_Fuel_Efficiency = \
	totalEnergyContentOfFuel / totalElectricalPower

	total_Input_Power_to_Hyd_to_Fuel_Efficiency = totalEnergyContentOfFuel / totalInputPower

	outputDict = {}
	outputDict['effAvailElectricalToFuel'] = available_Electrical_to_Hyd_to_Fuel_Efficiency
	outputDict['effTotalElectricalToFuel'] = total_Electrical_to_Hyd_to_Fuel_Efficiency
	outputDict['effTotalPowerToFuel'] = total_Input_Power_to_Hyd_to_Fuel_Efficiency
	outputDict['totalEnergyContentOfFuel'] = totalEnergyContentOfFuel
	outputDict['hydrogenCurrent'] = hydrogenCurrent
	outputDict['currentCell'] = currentCell
	
	
	# Calculate the mass production rate of fuel (grams per second) if information is available
	if molecularWeightFuelMolecule > -1:
		Update_Output_Dict_with_Fuel_Mass_Efficiency(outputDict, molecularWeightFuelMolecule, \
		totalFuelProductionRateInHyd, availableElectricalPower, totalElectricalPower, \
		totalInputPower)

	
	return outputDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Efficiencies_Current_TotalEnergyContentOfFuel_LoP_BioCO2(availableElectricalPower, \
totalElectricalPower, totalInputPower, voltageCellTwo, efficiencyCurrentToSecondCell, \
energyPerFuelMolecule, electronsPerFuelInLoP, molecularWeightFuelMolecule=-1):

	from rewiredcarbon.physicalconstants import elementaryCharge as e
	import pdb
	
	# Calculate the currents in the electrochemical cell. Defined in Buz Medium #1 page 22 and 
	# Farshid #1 page 32.
	# Updated so that cell voltage is just one variable, voltageCellTwo. 
	currentCell = availableElectricalPower / voltageCellTwo
	loPCurrent = efficiencyCurrentToSecondCell * currentCell

	
	# Calculate the total rate of fuel production (fuel molecules per second)
	# See Farshid #1 page 32 and Buz Medium #1 pages 22-23. 
	totalFuelProductionRateInLoP = loPCurrent / (e * electronsPerFuelInLoP)
	totalEnergyContentOfFuel = totalFuelProductionRateInLoP * energyPerFuelMolecule
		
	
	# Calculate the low potential mediator to fuel efficiencies. Defined in Farshid #1 page 37.
	
	available_Electrical_to_LoP_to_Fuel_Efficiency = \
	totalEnergyContentOfFuel / availableElectricalPower

	total_Electrical_to_LoP_to_Fuel_Efficiency = \
	totalEnergyContentOfFuel / totalElectricalPower

	total_Input_Power_to_LoP_to_Fuel_Efficiency = totalEnergyContentOfFuel / totalInputPower

	outputDict = {}
	outputDict['effAvailElectricalToFuel'] = available_Electrical_to_LoP_to_Fuel_Efficiency
	outputDict['effTotalElectricalToFuel'] = total_Electrical_to_LoP_to_Fuel_Efficiency
	outputDict['effTotalPowerToFuel'] = total_Input_Power_to_LoP_to_Fuel_Efficiency
	outputDict['totalEnergyContentOfFuel'] = totalEnergyContentOfFuel
	outputDict['loPCurrent'] = loPCurrent
	outputDict['currentCell'] = currentCell
	
	
	# Calculate the mass production rate of fuel (grams per second) if information is available
	if molecularWeightFuelMolecule > -1:
		Update_Output_Dict_with_Fuel_Mass_Efficiency(outputDict, molecularWeightFuelMolecule, \
		totalFuelProductionRateInLoP, availableElectricalPower, totalElectricalPower, \
		totalInputPower)

	
	return outputDict
# ------------------------------------------------------------------------------------------------ #






# ------------------------------------------------------------------------------------------------ #
def Efficiencies_Current_TotalEnergyContentOfFuel_EEU_BioCO2(availableElectricalPower, \
totalElectricalPower, totalInputPower, voltageCellTwo, efficiencyElectronsToMtr, \
efficiencyElectronsMtrToQuinone, energyPerFuelMolecule, electronsPerFuelInEEU, \
molecularWeightFuelMolecule):

	from rewiredcarbon.physicalconstants import elementaryCharge as e
	
	# Calculate the current in the electrochemical cell.
	currentCell = availableElectricalPower / voltageCellTwo

	# This the electrical current (think electrons) carried from the electrochemical cell to the
	# Mtr EEU complex, and then to quinone pool. Measured in Amps. 
	# Defined in Farshid #1 page 28. 
	mtrCurrent = efficiencyElectronsToMtr * currentCell
	quinoneCurrent = efficiencyElectronsMtrToQuinone * mtrCurrent

	# Calculate the total rate of fuel production (fuel molecules per second)
	# See Farshid #1 page 30. 
	totalFuelProductionRateInEEU = quinoneCurrent / (e * electronsPerFuelInEEU)
	
	# Calculate the energy content of the fuel made in one second (in Joules)
	totalEnergyContentOfFuel = totalFuelProductionRateInEEU * energyPerFuelMolecule


	# Calculate the EEU to fuel efficiencies 
	# Defined in Farshid #1 page 37.

	available_Electrical_to_EEU_to_Fuel_Efficiency = \
	totalEnergyContentOfFuel / availableElectricalPower

	total_Electrical_to_EEU_to_Fuel_Efficiency = \
	totalEnergyContentOfFuel / totalElectricalPower
	
	total_Input_Power_to_EEU_to_Fuel_Efficiency = totalEnergyContentOfFuel / totalInputPower


	outputDict = {}
	outputDict['effAvailElectricalToFuel'] = available_Electrical_to_EEU_to_Fuel_Efficiency
	outputDict['effTotalElectricalToFuel'] = total_Electrical_to_EEU_to_Fuel_Efficiency
	outputDict['effTotalPowerToFuel'] = total_Input_Power_to_EEU_to_Fuel_Efficiency
	outputDict['totalEnergyContentOfFuel'] = totalEnergyContentOfFuel
	outputDict['mtrCurrent'] = mtrCurrent
	outputDict['quinoneCurrent'] = quinoneCurrent
	outputDict['currentCell'] = currentCell
	
	# Calculate the mass production rate of fuel (grams per second) if information is available
	Update_Output_Dict_with_Fuel_Mass_Efficiency(outputDict, molecularWeightFuelMolecule, \
	totalFuelProductionRateInEEU, availableElectricalPower, totalElectricalPower, \
	totalInputPower)


	return outputDict
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def Efficiencies_Current_TotalEnergyContentOfFuel_Hydrogen_ElectrochemicalCO2(\
availableElectricalPower, totalElectricalPower, totalInputPower, vCellOne, vCellTwo, \
energyPerFuelMolecule, primaryFixPerFuel, electronsToAddInHyd, electronsPerPrimaryFix, \
carbonsPerPrimaryFix, efficiencyCurrentToFirstCell, efficiencyCurrentToSecondCell, \
efficiencyCarbonToSecondCell, molecularWeightFuelMolecule):

	from rewiredcarbon.physicalconstants import elementaryCharge as e
	import pdb
	
	# Calculate the currents through first and second electrochemical cells.
    # Defined in Farshid #1 pages 36 and 37.
    
    # This was later corrected on 2024-12-19 by Buz Barstow. I found that this was incorrectly 
    # predicting (over-estimating) efficiency for two ostensibly similar models with and without
    # carbon recycling. carbonsPerFuel should really be primaryFixPerFuel.
 	
	if carbonsPerPrimaryFix > 0:
		ICellOneDivICellTwo = \
		(primaryFixPerFuel * electronsPerPrimaryFix * efficiencyCurrentToSecondCell)\
		/(electronsToAddInHyd * efficiencyCurrentToFirstCell * carbonsPerPrimaryFix \
		* efficiencyCarbonToSecondCell)
	elif carbonsPerPrimaryFix == 0:
		ICellOneDivICellTwo = 0
	
	ICellTwo = abs(totalElectricalPower)/((ICellOneDivICellTwo * vCellOne) + vCellTwo)
	ICellOne = ICellTwo * ICellOneDivICellTwo
	
	hydrogenCurrent = ICellTwo * efficiencyCurrentToSecondCell
	
	# Calculate the total rate of fuel production (fuel molecules per second)
	# See Farshid #1 page 30. 
	totalFuelProductionRateInHyd = hydrogenCurrent / (e * electronsToAddInHyd)
	
	# Calculate the energy content of the fuel made in one second (in Joules)
	totalEnergyContentOfFuel = totalFuelProductionRateInHyd * energyPerFuelMolecule
	
	
	# Calculate the EEU to fuel efficiencies. Defined in Farshid #1 page 37.
	available_Electrical_to_Hyd_to_Fuel_Efficiency = \
	totalEnergyContentOfFuel / availableElectricalPower
	
	total_Electrical_to_Hyd_to_Fuel_Efficiency = \
	totalEnergyContentOfFuel / totalElectricalPower
	
	total_Input_Power_to_Hyd_to_Fuel_Efficiency = \
	totalEnergyContentOfFuel / totalInputPower
	
	
	outputDict = {}
	outputDict['effAvailElectricalToFuel'] = available_Electrical_to_Hyd_to_Fuel_Efficiency
	outputDict['effTotalElectricalToFuel'] = total_Electrical_to_Hyd_to_Fuel_Efficiency
	outputDict['effTotalPowerToFuel'] = total_Input_Power_to_Hyd_to_Fuel_Efficiency
	outputDict['totalEnergyContentOfFuel'] = totalEnergyContentOfFuel
	outputDict['ICellOne'] = ICellOne
	outputDict['ICellTwo'] = ICellTwo
	outputDict['hydrogenCurrent'] = hydrogenCurrent
	outputDict['currentCell'] = ICellTwo
	
	# Calculate the mass production rate of fuel (grams per second) if information is available
	Update_Output_Dict_with_Fuel_Mass_Efficiency(outputDict, molecularWeightFuelMolecule, \
	totalFuelProductionRateInHyd, availableElectricalPower, totalElectricalPower, \
	totalInputPower)


	return outputDict
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def Efficiencies_Current_TotalEnergyContentOfFuel_LoP_ElectrochemicalCO2(\
availableElectricalPower, totalElectricalPower, totalInputPower, vCellOne, vCellTwo, \
energyPerFuelMolecule, primaryFixPerFuel, electronsToAddInLoP, electronsPerPrimaryFix, \
carbonsPerPrimaryFix, efficiencyCurrentToFirstCell, efficiencyCurrentToSecondCell, \
efficiencyCarbonToSecondCell, molecularWeightFuelMolecule):

	from rewiredcarbon.physicalconstants import elementaryCharge as e
	import pdb
	
	# Calculate the currents through first and second electrochemical cells.
    # Defined in Farshid #1 pages 36 and 37.
    
    # This was later corrected on 2024-12-19 by Buz Barstow. I found that this was incorrectly 
    # predicting (over-estimating) efficiency for two ostensibly similar models with and without
    # carbon recycling. carbonsPerFuel should really be primaryFixPerFuel.
 	
	if carbonsPerPrimaryFix > 0:
		ICellOneDivICellTwo = \
		(primaryFixPerFuel * electronsPerPrimaryFix * efficiencyCurrentToSecondCell)\
		/(electronsToAddInLoP * efficiencyCurrentToFirstCell * carbonsPerPrimaryFix \
		* efficiencyCarbonToSecondCell)
	elif carbonsPerPrimaryFix == 0:
		ICellOneDivICellTwo = 0
	
	ICellTwo = abs(totalElectricalPower)/((ICellOneDivICellTwo * vCellOne) + vCellTwo)
	ICellOne = ICellTwo * ICellOneDivICellTwo
	
	loPCurrent = ICellTwo * efficiencyCurrentToSecondCell
	
	# Calculate the total rate of fuel production (fuel molecules per second)
	# See Farshid #1 page 30. 
	totalFuelProductionRateInLoP = loPCurrent / (e * electronsToAddInLoP)
	
	# Calculate the energy content of the fuel made in one second (in Joules)
	totalEnergyContentOfFuel = totalFuelProductionRateInLoP * energyPerFuelMolecule
	
	
	# Calculate the EEU to fuel efficiencies. Defined in Farshid #1 page 37.
	available_Electrical_to_LoP_to_Fuel_Efficiency = \
	totalEnergyContentOfFuel / availableElectricalPower
	
	total_Electrical_to_LoP_to_Fuel_Efficiency = \
	totalEnergyContentOfFuel / totalElectricalPower
	
	total_Input_Power_to_LoP_to_Fuel_Efficiency = \
	totalEnergyContentOfFuel / totalInputPower
	
	
	outputDict = {}
	outputDict['effAvailElectricalToFuel'] = available_Electrical_to_LoP_to_Fuel_Efficiency
	outputDict['effTotalElectricalToFuel'] = total_Electrical_to_LoP_to_Fuel_Efficiency
	outputDict['effTotalPowerToFuel'] = total_Input_Power_to_LoP_to_Fuel_Efficiency
	outputDict['totalEnergyContentOfFuel'] = totalEnergyContentOfFuel
	outputDict['ICellOne'] = ICellOne
	outputDict['ICellTwo'] = ICellTwo
	outputDict['hydrogenCurrent'] = loPCurrent
	outputDict['currentCell'] = ICellTwo
	
	# Calculate the mass production rate of fuel (grams per second) if information is available
	Update_Output_Dict_with_Fuel_Mass_Efficiency(outputDict, molecularWeightFuelMolecule, \
	totalFuelProductionRateInLoP, availableElectricalPower, totalElectricalPower, \
	totalInputPower)


	return outputDict
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def Efficiencies_Current_TotalEnergyContentOfFuel_EEU_ElectrochemicalCO2(totalElectricalPower, \
totalInputPower, vCellOne, vCellTwo, energyPerFuelMolecule, primaryFixPerFuel, electronsToAddInEEU,\
electronsPerPrimaryFix, carbonsPerPrimaryFix, efficiencyCurrentToFirstCell, \
efficiencyCurrentToSecondCell, efficiencyCarbonToSecondCell, \
molecularWeightFuelMolecule, efficiencyElectronsMtrToQuinone=1.0):

    # This was later corrected on 2024-12-19 by Buz Barstow. I found that this was incorrectly 
    # predicting (over-estimating) efficiency for two ostensibly similar models with and without
    # carbon recycling. carbonsPerFuel should really be primaryFixPerFuel.

	from rewiredcarbon.physicalconstants import elementaryCharge as e
	import pdb
	
	print('513')
	
	if 	carbonsPerPrimaryFix > 0:
		ICellOneDivICellTwo = \
		(primaryFixPerFuel * electronsPerPrimaryFix * efficiencyCurrentToSecondCell)\
		/(electronsToAddInEEU * efficiencyCurrentToFirstCell * carbonsPerPrimaryFix \
		* efficiencyCarbonToSecondCell)
	elif carbonsPerPrimaryFix == 0:
		ICellOneDivICellTwo = 0
	
	ICellTwo = efficiencyCurrentToSecondCell*\
	abs(totalElectricalPower)/((ICellOneDivICellTwo * vCellOne) + vCellTwo)
	ICellOne = ICellTwo * ICellOneDivICellTwo
	
	print('526')
	
	quinoneCurrent = ICellTwo * efficiencyElectronsMtrToQuinone

	# Calculate the total rate of fuel production (fuel molecules per second)
	# See Farshid #1 page 30. 
	totalFuelProductionRateInEEU = quinoneCurrent / (e * electronsToAddInEEU)
	
	# Calculate the energy content of the fuel made in one second (in Joules)
	totalEnergyContentOfFuel = totalFuelProductionRateInEEU * energyPerFuelMolecule
	
	print('537')
	
	
	# Calculate the EEU to fuel efficiencies. Defined in Farshid #1 page 37.
	available_Electrical_to_EEU_to_Fuel_Efficiency = \
	totalEnergyContentOfFuel / totalElectricalPower
	
	total_Electrical_to_EEU_to_Fuel_Efficiency = \
	totalEnergyContentOfFuel / totalElectricalPower
	
	total_Input_Power_to_EEU_to_Fuel_Efficiency = \
	totalEnergyContentOfFuel / totalInputPower
	
	print('548')
	
	outputDict = {}
	outputDict['effAvailElectricalToFuel'] = available_Electrical_to_EEU_to_Fuel_Efficiency
	outputDict['effTotalElectricalToFuel'] = total_Electrical_to_EEU_to_Fuel_Efficiency
	outputDict['effTotalPowerToFuel'] = total_Input_Power_to_EEU_to_Fuel_Efficiency
	outputDict['totalEnergyContentOfFuel'] = totalEnergyContentOfFuel
	outputDict['ICellOne'] = ICellOne
	outputDict['ICellTwo'] = ICellTwo
	outputDict['quinoneCurrent'] = quinoneCurrent
	outputDict['currentCell'] = ICellTwo
	
	# Calculate the mass production rate of fuel (grams per second) if information is available
	Update_Output_Dict_with_Fuel_Mass_Efficiency(outputDict, molecularWeightFuelMolecule, \
	totalFuelProductionRateInEEU, totalElectricalPower, totalElectricalPower, \
	totalInputPower)

	print('565')

	return outputDict
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def Efficiency_Hydrogen_BioCO2_to_Fuel_No_ScaleUp(voltageCellTwo, vMembrane, vAcceptor, \
NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, efficiencyElectronsToHyd=1.0, 
numberOfProtonsPumpedInForATP=-1, stirPower=0.0, totalElectricalPower=330.0, \
totalInputPower=1000.0, maxAllowedProtonsPumpedOutPerElectronDown=-1, \
molecularWeightFuelMolecule=-1):
# Calculate electrical to hydrogen to fuel efficiency 
# This is a really simple version where we don't do any scale up calculations beyond considering
# what the stir power does to efficiency

	from rewiredcarbon.electrochemicalconstants import deltaGADPATP
	from rewiredcarbon.physicalconstants import elementaryCharge as e
	import sys
	import pdb
	from copy import deepcopy

	# First, calculate proton pumping and electron per fuel molecule requirements. This is 
	# independent of any scale up considerations. 

	numberOfProtonsPumpedOut, electronsPerFuelInH, numberOfProtonsPumpedInForATP = \
	\
	ProtonPumping_and_Electron_Requirements_for_Hydrogen_to_Fuel(vMembrane, vAcceptor, \
	deltaGADPATP, NADHforFuel, FdForFuel, ATPforFuel, \
	numberOfProtonsPumpedInForATP=numberOfProtonsPumpedInForATP, \
	maxAllowedProtonsPumpedOutPerElectronDown=maxAllowedProtonsPumpedOutPerElectronDown)

	availableElectricalPower = totalElectricalPower - stirPower
	
	
	outputDict = \
	Efficiencies_Current_TotalEnergyContentOfFuel_Hydrogen_BioCO2(availableElectricalPower, \
	totalElectricalPower, totalInputPower, voltageCellTwo, efficiencyElectronsToHyd, \
	energyPerFuelMolecule, electronsPerFuelInH, \
	molecularWeightFuelMolecule=molecularWeightFuelMolecule)
	
	efficiencyDict = deepcopy(outputDict)
	
	efficiencyDict['protonsPumpedOut'] = numberOfProtonsPumpedOut
	efficiencyDict['electronsPerFuel'] = electronsPerFuelInH
	efficiencyDict['numberOfProtonsPumpedInForATP'] = numberOfProtonsPumpedInForATP
	
	return efficiencyDict
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def Efficiency_LoP_BioCO2_to_Fuel_No_ScaleUp(voltageCellTwo, vMediator, vMembrane, vAcceptor, \
NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, efficiencyCurrentToSecondCell=1.0, 
numberOfProtonsPumpedInForATP=-1, stirPower=0.0, totalElectricalPower=330.0, \
totalInputPower=1000.0, maxAllowedProtonsPumpedOutPerElectronDown=-1, \
molecularWeightFuelMolecule=-1):
# Calculate electrical to low potential mediator to fuel efficiency 

	from rewiredcarbon.electrochemicalconstants import deltaGADPATP
	from rewiredcarbon.physicalconstants import elementaryCharge as e
	import sys
	import pdb
	from copy import deepcopy

	# First, calculate proton pumping and electron per fuel molecule requirements. This is 
	# independent of any scale up considerations. 
	
	print('831')

	numberOfProtonsPumpedOut, electronsPerFuelInLoP, numberOfProtonsPumpedInForATP = \
	\
	ProtonPumping_and_Electron_Requirements_for_LoP_to_Fuel(vMediator, vMembrane, vAcceptor, \
	deltaGADPATP, NADHforFuel, FdForFuel, ATPforFuel, \
	numberOfProtonsPumpedInForATP=numberOfProtonsPumpedInForATP, \
	maxAllowedProtonsPumpedOutPerElectronDown=maxAllowedProtonsPumpedOutPerElectronDown)

	availableElectricalPower = totalElectricalPower - stirPower
	
	
	outputDict = \
	Efficiencies_Current_TotalEnergyContentOfFuel_LoP_BioCO2(availableElectricalPower, \
	totalElectricalPower, totalInputPower, voltageCellTwo, efficiencyCurrentToSecondCell, \
	energyPerFuelMolecule, electronsPerFuelInLoP, \
	molecularWeightFuelMolecule=molecularWeightFuelMolecule)
	
	print('849')
	
	efficiencyDict = deepcopy(outputDict)
	
	efficiencyDict['protonsPumpedOut'] = numberOfProtonsPumpedOut
	efficiencyDict['electronsPerFuel'] = electronsPerFuelInLoP
	efficiencyDict['numberOfProtonsPumpedInForATP'] = numberOfProtonsPumpedInForATP
	
	return efficiencyDict
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def Efficiency_EEU_BioCO2_to_Fuel_No_ScaleUp(voltageCellTwo, vMembrane, vAcceptor, vMtr, vQuinone, \
NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, molecularWeightFuelMolecule, \
efficiencyElectronsToMtr=1.0, efficiencyElectronsMtrToQuinone=1.0, \
numberOfProtonsPumpedInForATP=-1, totalElectricalPower=330.0, totalInputPower=1000.0, \
maxAllowedProtonsPumpedOutPerElectronDown=-1):

	from rewiredcarbon.electrochemicalconstants import deltaGADPATP
	import pdb
	from copy import deepcopy
	
	if vMtr > vQuinone:
		print('Cannot transfer electron down due to cytochrome Mtr being more positive than' \
		+ ' menaquinone')
		sys.exit()

	availableElectricalPower = totalElectricalPower
	
	numberOfProtonsPumpedOut, electronsPerFuelInEEU, numberOfProtonsPumpedInForATP = \
	ProtonPumping_and_Electron_Requirements_for_EEU_to_Fuel(vMembrane, vAcceptor, vMtr, \
	vQuinone, deltaGADPATP, NADHforFuel, FdForFuel, ATPforFuel, \
	numberOfProtonsPumpedInForATP=numberOfProtonsPumpedInForATP, \
	maxAllowedProtonsPumpedOutPerElectronDown=maxAllowedProtonsPumpedOutPerElectronDown)

	print('643')
	
	outputDict = \
	Efficiencies_Current_TotalEnergyContentOfFuel_EEU_BioCO2(availableElectricalPower, \
	totalElectricalPower, totalInputPower, voltageCellTwo, efficiencyElectronsToMtr, \
	efficiencyElectronsMtrToQuinone, energyPerFuelMolecule, electronsPerFuelInEEU, \
	molecularWeightFuelMolecule)
	
	efficiencyDict = deepcopy(outputDict)

	efficiencyDict['protonsPumpedOut'] = numberOfProtonsPumpedOut
	efficiencyDict['electronsPerFuel'] = electronsPerFuelInEEU
	efficiencyDict['numberOfProtonsPumpedInForATP'] = numberOfProtonsPumpedInForATP

	
	return efficiencyDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Efficiency_Hydrogen_ElectrochemCO2_to_Fuel_No_ScaleUp(voltageCellOne, voltageCellTwo, \
vMembrane, vAcceptor, NADHforFuel, FdForFuel, ATPforFuel, \
energyPerFuelMolecule, primaryFixPerFuel, electronsPerPrimaryFix, carbonsPerPrimaryFix, \
molecularWeightFuelMolecule, \
efficiencyCurrentToFirstCell=1, efficiencyCurrentToSecondCell=1, efficiencyCarbonToSecondCell=1, \
numberOfProtonsPumpedInForATP=-1, availableElectricalPower=330.0, totalInputPower=1000.0):
	
	print('670')
	
	from rewiredcarbon.electrochemicalconstants import deltaGADPATP
	from copy import deepcopy
	import pdb
	 
	if efficiencyCurrentToFirstCell == 0 or efficiencyCurrentToSecondCell == 0 or \
	efficiencyCarbonToSecondCell == 0:
		return 0
	
	
	numberOfProtonsPumpedOut, electronsToAddInHyd, numberOfProtonsPumpedInForATP = \
	ProtonPumping_and_Electron_Requirements_for_Hydrogen_to_Fuel(vMembrane, vAcceptor, \
	deltaGADPATP, NADHforFuel, FdForFuel, ATPforFuel, \
	numberOfProtonsPumpedInForATP=numberOfProtonsPumpedInForATP)
	
	print('686')
	
	outputDict = \
	Efficiencies_Current_TotalEnergyContentOfFuel_Hydrogen_ElectrochemicalCO2(\
	availableElectricalPower, availableElectricalPower, totalInputPower, voltageCellOne, \
	voltageCellTwo, energyPerFuelMolecule, primaryFixPerFuel, electronsToAddInHyd, \
	electronsPerPrimaryFix, carbonsPerPrimaryFix, efficiencyCurrentToFirstCell, \
	efficiencyCurrentToSecondCell, efficiencyCarbonToSecondCell, molecularWeightFuelMolecule)
	
	print('696')
	
	efficiencyDict = deepcopy(outputDict)
	efficiencyDict['protonsPumpedOut'] = numberOfProtonsPumpedOut
	efficiencyDict['electronsPerFuel'] = electronsToAddInHyd
	 
	 
	return efficiencyDict
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def Efficiency_LoP_ElectrochemCO2_to_Fuel_No_ScaleUp(voltageCellOne, voltageCellTwo, vMediator, \
vMembrane, vAcceptor, NADHforFuel, FdForFuel, ATPforFuel, \
energyPerFuelMolecule, primaryFixPerFuel, electronsPerPrimaryFix, carbonsPerPrimaryFix, \
molecularWeightFuelMolecule, \
efficiencyCurrentToFirstCell=1, efficiencyCurrentToSecondCell=1, efficiencyCarbonToSecondCell=1, \
numberOfProtonsPumpedInForATP=-1, availableElectricalPower=330.0, totalInputPower=1000.0):
	
	print('In Efficiency_LoP_ElectrochemCO2_to_Fuel_No_ScaleUp')
	
	from rewiredcarbon.electrochemicalconstants import deltaGADPATP
	from copy import deepcopy
	import pdb
	 
	if efficiencyCurrentToFirstCell == 0 or efficiencyCurrentToSecondCell == 0 or \
	efficiencyCarbonToSecondCell == 0:
		return 0
	
	
	numberOfProtonsPumpedOut, electronsToAddInLoP, numberOfProtonsPumpedInForATP = \
	ProtonPumping_and_Electron_Requirements_for_LoP_to_Fuel(vMediator, vMembrane, vAcceptor, \
	deltaGADPATP, NADHforFuel, FdForFuel, ATPforFuel, \
	numberOfProtonsPumpedInForATP=numberOfProtonsPumpedInForATP)
	
	print('974')
	
	outputDict = \
	Efficiencies_Current_TotalEnergyContentOfFuel_LoP_ElectrochemicalCO2(\
	availableElectricalPower, availableElectricalPower, totalInputPower, voltageCellOne, \
	voltageCellTwo, energyPerFuelMolecule, primaryFixPerFuel, electronsToAddInLoP, \
	electronsPerPrimaryFix, carbonsPerPrimaryFix, efficiencyCurrentToFirstCell, \
	efficiencyCurrentToSecondCell, efficiencyCarbonToSecondCell, molecularWeightFuelMolecule)
	
	print('982')
	
	efficiencyDict = deepcopy(outputDict)
	efficiencyDict['protonsPumpedOut'] = numberOfProtonsPumpedOut
	efficiencyDict['electronsPerFuel'] = electronsToAddInLoP
	 
	 
	return efficiencyDict
# ------------------------------------------------------------------------------------------------ #






# ------------------------------------------------------------------------------------------------ #
def Efficiency_EEU_ElectrochemCO2_to_Fuel_No_ScaleUp(voltageCellOne, voltageCellTwo, vMembrane, \
vAcceptor, vMtr, vQuinone, NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, \
primaryFixPerFuel, electronsPerPrimaryFix, carbonsPerPrimaryFix, molecularWeightFuelMolecule, \
efficiencyCurrentToFirstCell=1, efficiencyCurrentToSecondCell=1, efficiencyCarbonToSecondCell=1, \
numberOfProtonsPumpedInForATP=-1, totalElectricalPower=330.0, totalInputPower=1000.0, \
efficiencyElectronsMtrToQuinone=1.0):
	
	from rewiredcarbon.electrochemicalconstants import deltaGADPATP
	import pdb
	from copy import deepcopy
	
	print('716')
	
	if vMtr > vQuinone:
		print('Cannot transfer electron down due to cytochrome Mtr being more positive than' \
		+ ' menaquinone')
		sys.exit()
	
	if efficiencyCurrentToFirstCell == 0 or efficiencyCurrentToSecondCell == 0 or \
	efficiencyCarbonToSecondCell == 0:
		return 0
	
	
	numberOfProtonsPumpedOut, electronsToAddInEEU, numberOfProtonsPumpedInForATP = \
	ProtonPumping_and_Electron_Requirements_for_EEU_to_Fuel(vMembrane, vAcceptor, vMtr, vQuinone, \
	deltaGADPATP, NADHforFuel, FdForFuel, ATPforFuel, \
	numberOfProtonsPumpedInForATP=numberOfProtonsPumpedInForATP)
	
	print('733')
	
	outputDict = \
	Efficiencies_Current_TotalEnergyContentOfFuel_EEU_ElectrochemicalCO2(totalElectricalPower, \
	totalInputPower, voltageCellOne, voltageCellTwo, energyPerFuelMolecule, primaryFixPerFuel, \
	electronsToAddInEEU,\
	electronsPerPrimaryFix, carbonsPerPrimaryFix, efficiencyCurrentToFirstCell, \
	efficiencyCurrentToSecondCell, efficiencyCarbonToSecondCell, \
	molecularWeightFuelMolecule, efficiencyElectronsMtrToQuinone=1.0)
	
	efficiencyDict = deepcopy(outputDict)	
	efficiencyDict['protonsPumpedOut'] = numberOfProtonsPumpedOut
	efficiencyDict['electronsPerFuel'] = electronsToAddInEEU
	
	return efficiencyDict
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def Import_and_Calculate_vCellTwo(scenarioData, tol=0.001):
# I'm adding this in January 2025 as keeping track of whole cell voltages and biases becomes
# increasingly difficult as we consider formate as a mediator and carbon source, and think more
# about real world electrochemical cells. 	

	import sys
	
	try:
		voltageCellTwoCathode = float(scenarioData['voltageCellTwoCathode'])
		voltageCellTwoAnode = float(scenarioData['voltageCellTwoAnode'])
		voltageCellTwoCathodeBias = float(scenarioData['voltageCellTwoCathodeBias'])
		voltageCellTwoAnodeBias = float(scenarioData['voltageCellTwoAnodeBias'])
		voltageCellTwoOhmic = float(scenarioData['voltageCellTwoOhmic'])
	
		vRedoxCalc = abs(float(voltageCellTwoCathode) - float(voltageCellTwoAnode))
		vBiasCalc = abs(float(voltageCellTwoCathodeBias) + float(voltageCellTwoAnodeBias))
		voltageCellTwoCalc = vRedoxCalc + vBiasCalc + voltageCellTwoOhmic
		cellBiasVoltagesDefined = True
	except:
		cellBiasVoltagesDefined = False
	
	
	try:
		voltageCellTwo = float(scenarioData['voltageCellTwo'])
		wholeCellVoltageDefined = True
	except:
		wholeCellVoltageDefined = False
		
	
	if cellBiasVoltagesDefined and wholeCellVoltageDefined:
		if abs(voltageCellTwoCalc - voltageCellTwo) < tol:
			returnVoltage = voltageCellTwo
		else:
			print("Inconsistency between calculated and entered cell two voltages.")
			sys.exit()
	
	elif cellBiasVoltagesDefined == False and wholeCellVoltageDefined == False:
		print("Whole cell two voltage data is missing.")
		sys.exit()
	
	elif cellBiasVoltagesDefined == False and wholeCellVoltageDefined == True:
		returnVoltage = voltageCellTwo
	
	elif cellBiasVoltagesDefined == True and wholeCellVoltageDefined == False:
		returnVoltage = voltageCellTwoCalc
			
	
	return returnVoltage
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Import_and_Calculate_vCellOne(scenarioData, tol=0.001):
# I'm adding this in January 2025 as keeping track of whole cell voltages and biases becomes
# increasingly difficult as we consider formate as a mediator and carbon source, and think more
# about real world electrochemical cells. 	

	import sys
	
	try:
		voltageCellOneCathode = float(scenarioData['voltageCellOneCathode'])
		voltageCellOneAnode = float(scenarioData['voltageCellOneAnode'])
		voltageCellOneCathodeBias = float(scenarioData['voltageCellOneCathodeBias'])
		voltageCellOneAnodeBias = float(scenarioData['voltageCellOneAnodeBias'])
		voltageCellOneOhmic = float(scenarioData['voltageCellOneOhmic'])
	
		vRedoxCalc = abs(float(voltageCellOneCathode) - float(voltageCellOneAnode))
		vBiasCalc = abs(float(voltageCellOneCathodeBias) + float(voltageCellOneAnodeBias))
		voltageCellOneCalc = vRedoxCalc + vBiasCalc + voltageCellOneOhmic
		cellBiasVoltagesDefined = True
	except:
		cellBiasVoltagesDefined = False
	
	
	try:
		voltageCellOne = float(scenarioData['voltageCellOne'])
		wholeCellVoltageDefined = True
	except:
		wholeCellVoltageDefined = False
		
	
	if cellBiasVoltagesDefined and wholeCellVoltageDefined:
		if abs(voltageCellOneCalc - voltageCellOne) < tol:
			returnVoltage = voltageCellOne
		else:
			print("Inconsistency between calculated and entered cell one voltages.")
			sys.exit()
	
	elif cellBiasVoltagesDefined == False and wholeCellVoltageDefined == False:
		print("Whole cell one voltage data is missing.")
		sys.exit()
	
	elif cellBiasVoltagesDefined == False and wholeCellVoltageDefined == True:
		returnVoltage = voltageCellOne
	
	elif cellBiasVoltagesDefined == True and wholeCellVoltageDefined == False:
		returnVoltage = voltageCellOneCalc
			
	
	return returnVoltage
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def Process_Hydrogen_with_BioCO2_Scenario(scenarioData):
	
	import pdb
	from rewiredcarbon.scaleup import ImportAndExtraplotePowerNumberCurve
	import sys
	
	print('In H2')
	
	# These are variables that every scale up scenario needs
	NADHforFuel = float(scenarioData['NADHforFuel'])
	FdForFuel = float(scenarioData['FdForFuel'])
	ATPforFuel = float(scenarioData['ATPforFuel'])	
	vAcceptor = float(scenarioData['vAcceptor'])
	energyPerFuelMolecule = float(scenarioData['energyPerFuelMolecule'])
	efficiencyCurrentToSecondCell = float(scenarioData['efficiencyCurrentToSecondCell'])
	vMembrane = float(scenarioData['voltageMembrane'])/1000.0
	totalElectricalPower = float(scenarioData['totalElectricalPower'])
	totalInputPower = float(scenarioData['totalInputPower'])
	
	voltageCellTwo = Import_and_Calculate_vCellTwo(scenarioData)	
	
	try:
		maxAllowedProtonsPumpedOutPerElectronDown \
		= int(scenarioData['maxAllowedProtonsPumpedOutPerElectronDown'])
	except:
		maxAllowedProtonsPumpedOutPerElectronDown =-1
		
		
	try:
		molecularWeightFuelMolecule = float(scenarioData['molecularWeightFuelMolecule'])
	except:
		molecularWeightFuelMolecule = -1
			
	
	efficiencyDict = \
	Efficiency_Hydrogen_BioCO2_to_Fuel_No_ScaleUp(voltageCellTwo, vMembrane, \
	vAcceptor, NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, \
	efficiencyElectronsToHyd=efficiencyCurrentToSecondCell, numberOfProtonsPumpedInForATP=-1, 
	stirPower=0, totalElectricalPower=totalElectricalPower, \
	totalInputPower=totalInputPower, maxAllowedProtonsPumpedOutPerElectronDown=\
	maxAllowedProtonsPumpedOutPerElectronDown, \
	molecularWeightFuelMolecule=molecularWeightFuelMolecule)
	
		
	return efficiencyDict 
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Process_LoP_with_BioCO2_Scenario(scenarioData):
	
	import pdb
	from rewiredcarbon.scaleup import ImportAndExtraplotePowerNumberCurve
	import sys
	
	print('Process_LoP_with_BioCO2_Scenario')
	
	# These are variables that every scale up scenario needs
	NADHforFuel = float(scenarioData['NADHforFuel'])
	FdForFuel = float(scenarioData['FdForFuel'])
	ATPforFuel = float(scenarioData['ATPforFuel'])	
	vAcceptor = float(scenarioData['vAcceptor'])
	vMediator = float(scenarioData['vMediator'])
	energyPerFuelMolecule = float(scenarioData['energyPerFuelMolecule'])
	efficiencyCurrentToSecondCell = float(scenarioData['efficiencyCurrentToSecondCell'])
	vMembrane = float(scenarioData['voltageMembrane'])/1000.0
	totalElectricalPower = float(scenarioData['totalElectricalPower'])
	totalInputPower = float(scenarioData['totalInputPower'])
	
	voltageCellTwo = Import_and_Calculate_vCellTwo(scenarioData)	
	
	try:
		maxAllowedProtonsPumpedOutPerElectronDown \
		= int(scenarioData['maxAllowedProtonsPumpedOutPerElectronDown'])
	except:
		maxAllowedProtonsPumpedOutPerElectronDown =-1
		
		
	try:
		molecularWeightFuelMolecule = float(scenarioData['molecularWeightFuelMolecule'])
	except:
		molecularWeightFuelMolecule = -1
			
	print('1237')
	
	efficiencyDict = \
	Efficiency_LoP_BioCO2_to_Fuel_No_ScaleUp(voltageCellTwo, vMediator, vMembrane, \
	vAcceptor, NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, \
	efficiencyCurrentToSecondCell=efficiencyCurrentToSecondCell, numberOfProtonsPumpedInForATP=-1, 
	stirPower=0, totalElectricalPower=totalElectricalPower, \
	totalInputPower=totalInputPower, maxAllowedProtonsPumpedOutPerElectronDown=\
	maxAllowedProtonsPumpedOutPerElectronDown, \
	molecularWeightFuelMolecule=molecularWeightFuelMolecule)
	
	print('1248')
		
	return efficiencyDict 
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def Process_EEU_with_BioCO2_Scenario(scenarioData):
	
	import pdb
	
	print('858')
	
	
	NADHforFuel = float(scenarioData['NADHforFuel'])
	FdForFuel = float(scenarioData['FdForFuel'])
	ATPforFuel = float(scenarioData['ATPforFuel'])	
	vAcceptor = float(scenarioData['vAcceptor'])
	energyPerFuelMolecule = float(scenarioData['energyPerFuelMolecule'])
	vMtr = float(scenarioData['vMtr'])
	vQuinone = float(scenarioData['vQuinone'])
	vMembrane = float(scenarioData['voltageMembrane'])/1000.0
	efficiencyCurrentToSecondCell = float(scenarioData['efficiencyCurrentToSecondCell'])
	totalElectricalPower = float(scenarioData['totalElectricalPower'])
	totalInputPower = float(scenarioData['totalInputPower'])
	
	voltageCellTwo = Import_and_Calculate_vCellTwo(scenarioData)	
	
	try:
		maxAllowedProtonsPumpedOutPerElectronDown \
		= int(scenarioData['maxAllowedProtonsPumpedOutPerElectronDown'])
	except:
		maxAllowedProtonsPumpedOutPerElectronDown =-1
		
	try:
		molecularWeightFuelMolecule = float(scenarioData['molecularWeightFuelMolecule'])
	except:
		molecularWeightFuelMolecule = -1
		
	print('886')
	
	efficiencyDict = \
	Efficiency_EEU_BioCO2_to_Fuel_No_ScaleUp(voltageCellTwo, vMembrane, \
	vAcceptor, vMtr, vQuinone, NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, \
	molecularWeightFuelMolecule, \
	efficiencyElectronsToMtr=efficiencyCurrentToSecondCell, efficiencyElectronsMtrToQuinone=1.0, \
	numberOfProtonsPumpedInForATP=-1, totalElectricalPower=totalElectricalPower, \
	totalInputPower=totalInputPower, maxAllowedProtonsPumpedOutPerElectronDown=\
	maxAllowedProtonsPumpedOutPerElectronDown)
		
	print('897')
	
	return efficiencyDict
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def Process_Hydrogen_with_ElectrochemicalCO2_Scenario(scenarioData):
	
	print('960')
	print('In Electrochem H2')
	
	NADHforFuel = float(scenarioData['NADHforFuel'])
	FdForFuel = float(scenarioData['FdForFuel'])
	ATPforFuel = float(scenarioData['ATPforFuel'])	
	vAcceptor = float(scenarioData['vAcceptor'])
	energyPerFuelMolecule = float(scenarioData['energyPerFuelMolecule'])
	carbonsPerFuel = float(scenarioData['carbonsPerFuel'])
	
	electronsPerPrimaryFix = float(scenarioData['electronsPerPrimaryFix'])
	carbonsPerPrimaryFix = float(scenarioData['carbonsPerPrimaryFix'])
	primaryFixPerFuel = float(scenarioData['primaryFixPerFuel'])
	
	print('973')
	
	voltageCellTwo = Import_and_Calculate_vCellTwo(scenarioData)	
	voltageCellOne = Import_and_Calculate_vCellOne(scenarioData)	
	
	
	print('979')
	
	vMembrane = float(scenarioData['voltageMembrane'])/1000.
	
	totalElectricalPower = float(scenarioData['totalElectricalPower'])
	totalInputPower = float(scenarioData['totalInputPower'])

	efficiencyCurrentToFirstCell = float(scenarioData['efficiencyCurrentToFirstCell'])
	efficiencyCurrentToSecondCell = float(scenarioData['efficiencyCurrentToSecondCell'])
	
	
	try:
		molecularWeightFuelMolecule = float(scenarioData['molecularWeightFuelMolecule'])
	except:
		molecularWeightFuelMolecule = -1

	
	print('996')
	
	efficiencyDict = \
	Efficiency_Hydrogen_ElectrochemCO2_to_Fuel_No_ScaleUp(voltageCellOne, voltageCellTwo, \
	vMembrane, vAcceptor, NADHforFuel, FdForFuel, ATPforFuel, \
	energyPerFuelMolecule, primaryFixPerFuel, electronsPerPrimaryFix, carbonsPerPrimaryFix, \
	molecularWeightFuelMolecule, \
	efficiencyCurrentToFirstCell=efficiencyCurrentToFirstCell, \
	efficiencyCurrentToSecondCell=efficiencyCurrentToSecondCell, efficiencyCarbonToSecondCell=1, \
	numberOfProtonsPumpedInForATP=-1, availableElectricalPower=totalElectricalPower, \
	totalInputPower=totalInputPower)

	print('1008')
	
	return efficiencyDict
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def Process_LoP_with_ElectrochemicalCO2_Scenario(scenarioData):
# For low potential mediators
	
	print('Process_LoP_with_ElectrochemicalCO2_Scenario')
	
	NADHforFuel = float(scenarioData['NADHforFuel'])
	FdForFuel = float(scenarioData['FdForFuel'])
	ATPforFuel = float(scenarioData['ATPforFuel'])	
	vAcceptor = float(scenarioData['vAcceptor'])
	vMediator = float(scenarioData['vMediator'])
	energyPerFuelMolecule = float(scenarioData['energyPerFuelMolecule'])
	carbonsPerFuel = float(scenarioData['carbonsPerFuel'])
	
	electronsPerPrimaryFix = float(scenarioData['electronsPerPrimaryFix'])
	carbonsPerPrimaryFix = float(scenarioData['carbonsPerPrimaryFix'])
	primaryFixPerFuel = float(scenarioData['primaryFixPerFuel'])
	
	
	voltageCellTwo = Import_and_Calculate_vCellTwo(scenarioData)	
	voltageCellOne = Import_and_Calculate_vCellOne(scenarioData)	
	
		
	vMembrane = float(scenarioData['voltageMembrane'])/1000.
	
	totalElectricalPower = float(scenarioData['totalElectricalPower'])
	totalInputPower = float(scenarioData['totalInputPower'])

	efficiencyCurrentToFirstCell = float(scenarioData['efficiencyCurrentToFirstCell'])
	efficiencyCurrentToSecondCell = float(scenarioData['efficiencyCurrentToSecondCell'])
	
	
	try:
		molecularWeightFuelMolecule = float(scenarioData['molecularWeightFuelMolecule'])
	except:
		molecularWeightFuelMolecule = -1

	
	print('1399')
	
	efficiencyDict = \
	Efficiency_LoP_ElectrochemCO2_to_Fuel_No_ScaleUp(voltageCellOne, voltageCellTwo, vMediator, \
	vMembrane, vAcceptor, NADHforFuel, FdForFuel, ATPforFuel, \
	energyPerFuelMolecule, primaryFixPerFuel, electronsPerPrimaryFix, carbonsPerPrimaryFix, \
	molecularWeightFuelMolecule, \
	efficiencyCurrentToFirstCell=efficiencyCurrentToFirstCell, \
	efficiencyCurrentToSecondCell=efficiencyCurrentToSecondCell, efficiencyCarbonToSecondCell=1, \
	numberOfProtonsPumpedInForATP=-1, availableElectricalPower=totalElectricalPower, \
	totalInputPower=totalInputPower)

	print('1411')
	
	return efficiencyDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Process_EEU_with_ElectrochemicalCO2_Scenario(scenarioData):
	
	import pdb
	import sys
	
	print('957')
	print('In Electrochem EEU')
	
	NADHforFuel = float(scenarioData['NADHforFuel'])
	FdForFuel = float(scenarioData['FdForFuel'])
	ATPforFuel = float(scenarioData['ATPforFuel'])	
	
	vAcceptor = float(scenarioData['vAcceptor'])
	energyPerFuelMolecule = float(scenarioData['energyPerFuelMolecule'])
	carbonsPerFuel = float(scenarioData['carbonsPerFuel'])
	vMtr = float(scenarioData['vMtr'])
	vQuinone = float(scenarioData['vQuinone'])
	
	electronsPerPrimaryFix = float(scenarioData['electronsPerPrimaryFix'])
	carbonsPerPrimaryFix = float(scenarioData['carbonsPerPrimaryFix'])
	
	print('1044')
	voltageCellTwo = Import_and_Calculate_vCellTwo(scenarioData)	
	print('1046')
	voltageCellOne = Import_and_Calculate_vCellOne(scenarioData)	
	print('1048')
	
	vMembrane = float(scenarioData['voltageMembrane'])/1000.
	
	totalElectricalPower = float(scenarioData['totalElectricalPower'])
	totalInputPower = float(scenarioData['totalInputPower'])
	
	efficiencyCurrentToFirstCell = float(scenarioData['efficiencyCurrentToFirstCell'])
	efficiencyCurrentToSecondCell = float(scenarioData['efficiencyCurrentToSecondCell'])
	
	
	try:
		molecularWeightFuelMolecule = float(scenarioData['molecularWeightFuelMolecule'])
	except:
		molecularWeightFuelMolecule = -1

		
	primaryFixPerFuel = float(scenarioData['primaryFixPerFuel'])
	
	print('1059')
	
	efficiencyDict = \
	Efficiency_EEU_ElectrochemCO2_to_Fuel_No_ScaleUp(voltageCellOne, voltageCellTwo, vMembrane, \
	vAcceptor, vMtr, vQuinone, NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, \
	primaryFixPerFuel, electronsPerPrimaryFix, carbonsPerPrimaryFix, molecularWeightFuelMolecule, \
	efficiencyCurrentToFirstCell=efficiencyCurrentToFirstCell, \
	efficiencyCurrentToSecondCell=efficiencyCurrentToSecondCell, 
	efficiencyCarbonToSecondCell=1, totalElectricalPower=totalElectricalPower, \
	totalInputPower=totalInputPower)
	
	
	return efficiencyDict
# ------------------------------------------------------------------------------------------------ #

